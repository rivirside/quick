"""
Pipeline orchestration for eTrial validation workflow.

Coordinates execution of all validation modules, handles dependencies,
and produces final decisions.
"""

import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Type, Union
from loguru import logger

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    PipelineResult,
    Decision,
    PipelineError,
)
from etrial.core.config import Config, ThresholdManager
from etrial.core.audit import AuditTrail


class ValidationPipeline:
    """
    Orchestrates the complete validation workflow.

    Manages:
    - Module registration and dependency ordering
    - Configuration and threshold management
    - Audit trail and reproducibility
    - Decision gating (KILL stops pipeline for critical modules)
    - Result aggregation and reporting
    """

    def __init__(
        self,
        config: Optional[Config] = None,
        audit_trail: Optional[AuditTrail] = None,
        enable_gating: bool = True,
    ):
        """
        Initialize validation pipeline.

        Args:
            config: Configuration object (uses default if None)
            audit_trail: Audit trail instance (creates new if None)
            enable_gating: Whether to stop pipeline on KILL decisions
        """
        self.config = config or Config.default()
        self.audit_trail = audit_trail or AuditTrail(
            output_dir=self.config.output_dir / "audit"
        )
        self.enable_gating = enable_gating

        # Module registry
        self.modules: Dict[str, ValidationModule] = {}
        self.module_order: List[str] = []

        # Threshold management
        self.threshold_manager: Optional[ThresholdManager] = None

        # Track config hash for reproducibility
        self.audit_trail.set_config_hash(self.config.get_hash())

        logger.info(f"Pipeline initialized with config hash: {self.config.get_hash()}")

    @classmethod
    def from_config(cls, config_path: Union[str, Path]) -> "ValidationPipeline":
        """
        Create pipeline from configuration file.

        Args:
            config_path: Path to YAML configuration file

        Returns:
            Configured ValidationPipeline
        """
        config = Config.from_yaml(config_path)

        # Load thresholds if they exist
        thresholds_path = Path(config_path).parent / "thresholds.yaml"
        if thresholds_path.exists():
            config.load_thresholds(thresholds_path)

        pipeline = cls(config=config)

        # Set up threshold manager
        if config._thresholds:
            pipeline.threshold_manager = ThresholdManager(config._thresholds)

        return pipeline

    def register_module(
        self,
        module: Union[ValidationModule, Type[ValidationModule]],
        name: Optional[str] = None,
        position: Optional[int] = None,
    ) -> None:
        """
        Register a validation module.

        Args:
            module: ValidationModule instance or class
            name: Optional custom name (uses module class name if None)
            position: Optional position in execution order (appends if None)
        """
        # Instantiate if class was passed
        if isinstance(module, type):
            module_config = self.config.get_module_config(
                name or module.__name__.lower().replace('module', '')
            )
            module = module(config=module_config)

        module_name = name or module.name

        # Check if module is enabled
        if not self.config.is_module_enabled(module_name):
            logger.info(f"Module {module_name} is disabled, skipping registration")
            return

        self.modules[module_name] = module

        # Handle execution order
        if position is not None:
            self.module_order.insert(position, module_name)
        else:
            self.module_order.append(module_name)

        logger.info(f"Module registered: {module_name} (position {len(self.module_order)})")

        # Register module version in audit trail
        self.audit_trail.register_model_version(module_name, module.version)

    def validate(
        self,
        candidate: TherapeuticCandidate,
        modules: Optional[List[str]] = None,
        stop_on_kill: Optional[bool] = None,
        **kwargs
    ) -> PipelineResult:
        """
        Run full validation pipeline on a candidate.

        Args:
            candidate: TherapeuticCandidate to validate
            modules: Optional list of specific modules to run (runs all if None)
            stop_on_kill: Override enable_gating for this run
            **kwargs: Additional arguments passed to modules

        Returns:
            PipelineResult with all module results and final decision
        """
        start_time = time.time()
        stop_on_kill = stop_on_kill if stop_on_kill is not None else self.enable_gating

        # Log validation start
        self.audit_trail.log_event(
            "validation_started",
            candidate_hash=candidate.get_hash(),
            data={
                "candidate_name": candidate.name,
                "modality": candidate.modality.value,
                "target": candidate.target,
            }
        )

        logger.info(f"Starting validation for {candidate.name} ({candidate.modality.value})")

        # Determine which modules to run
        modules_to_run = modules if modules else self.module_order

        # Filter out unregistered modules
        modules_to_run = [m for m in modules_to_run if m in self.modules]

        if not modules_to_run:
            raise PipelineError("No modules available to run")

        # Execute modules in order
        module_results: Dict[str, ValidationResult] = {}
        execution_order: List[str] = []

        for module_name in modules_to_run:
            module = self.modules[module_name]

            logger.info(f"Running module: {module_name}")

            # Log module start
            self.audit_trail.log_event(
                "module_started",
                module_name=module_name,
                candidate_hash=candidate.get_hash(),
            )

            try:
                # Run validation
                module_start = time.time()
                result = module.validate(candidate, **kwargs)
                module_runtime = time.time() - module_start

                # Update runtime if not already set
                if result.runtime_seconds is None:
                    result.runtime_seconds = module_runtime

                module_results[module_name] = result
                execution_order.append(module_name)

                # Log module completion
                self.audit_trail.log_event(
                    "module_completed",
                    module_name=module_name,
                    candidate_hash=candidate.get_hash(),
                    data={
                        "decision": result.decision.value,
                        "runtime_seconds": module_runtime,
                        "n_metrics": len(result.metrics),
                    }
                )

                logger.info(
                    f"Module {module_name} completed: {result.decision.value} "
                    f"({module_runtime:.2f}s)"
                )

                # Check for gating
                if stop_on_kill and result.decision == Decision.KILL:
                    if self._is_critical_module(module_name):
                        logger.warning(
                            f"Critical module {module_name} returned KILL, stopping pipeline"
                        )
                        break

            except Exception as e:
                logger.error(f"Module {module_name} failed: {e}")

                # Log module failure
                self.audit_trail.log_event(
                    "module_failed",
                    module_name=module_name,
                    candidate_hash=candidate.get_hash(),
                    data={"error": str(e)}
                )

                raise PipelineError(f"Module {module_name} failed: {e}") from e

        # Determine overall decision
        overall_decision = self._determine_overall_decision(module_results)

        # Calculate total runtime
        total_runtime = time.time() - start_time

        # Create pipeline result
        pipeline_result = PipelineResult(
            candidate=candidate,
            overall_decision=overall_decision,
            module_results=module_results,
            execution_order=execution_order,
            total_runtime_seconds=total_runtime,
            config_hash=self.config.get_hash(),
            audit_trail_id=self.audit_trail.trail_id,
        )

        # Log validation completion
        self.audit_trail.log_event(
            "validation_completed",
            candidate_hash=candidate.get_hash(),
            data={
                "overall_decision": overall_decision.value,
                "total_runtime_seconds": total_runtime,
                "modules_run": len(module_results),
            }
        )

        logger.info(
            f"Validation completed: {overall_decision.value} "
            f"({total_runtime:.2f}s, {len(module_results)} modules)"
        )

        # Save audit trail
        self.audit_trail.save()

        return pipeline_result

    def _is_critical_module(self, module_name: str) -> bool:
        """
        Check if a module is critical (KILL should stop pipeline).

        Args:
            module_name: Name of the module

        Returns:
            True if critical, False otherwise
        """
        if not self.threshold_manager:
            # Default critical modules
            return module_name in ['structure', 'specificity', 'safety']

        # Check gates configuration
        gates = self.threshold_manager.thresholds.get('gates', {})
        module_gate = gates.get(module_name, {})

        return module_gate.get('fail_stops_pipeline', False)

    def _determine_overall_decision(
        self,
        module_results: Dict[str, ValidationResult]
    ) -> Decision:
        """
        Determine overall decision based on module results.

        Logic:
        1. If any critical module is KILL → overall KILL
        2. If too many modules are REVISE → overall KILL
        3. If any module is REVISE → overall REVISE
        4. Otherwise → overall PASS

        Args:
            module_results: Dict of module results

        Returns:
            Overall Decision
        """
        if not module_results:
            return Decision.INFORMATIVE

        # Extract module decisions
        module_decisions = {
            name: result.decision.value
            for name, result in module_results.items()
        }

        # Use threshold manager if available
        if self.threshold_manager:
            decision_str = self.threshold_manager.get_overall_decision(module_decisions)
            return Decision[decision_str]

        # Fallback logic
        decisions = list(module_decisions.values())

        # Check for any KILL
        if 'KILL' in decisions:
            return Decision.KILL

        # Check for too many REVISE
        if decisions.count('REVISE') >= 4:
            return Decision.KILL

        # Check for any REVISE
        if 'REVISE' in decisions:
            return Decision.REVISE

        # All PASS
        return Decision.PASS

    def generate_dossier(
        self,
        pipeline_result: PipelineResult,
        output_path: Optional[Path] = None,
        format: str = "pdf",
    ) -> Path:
        """
        Generate validation dossier (report).

        Args:
            pipeline_result: Pipeline result to document
            output_path: Optional custom output path
            format: Output format ('pdf', 'html', 'json')

        Returns:
            Path to generated dossier
        """
        from etrial.core.reporting import ReportGenerator

        if output_path is None:
            timestamp = pipeline_result.timestamp.strftime("%Y%m%d_%H%M%S")
            filename = f"{pipeline_result.candidate.name}_dossier_{timestamp}.{format}"
            output_path = self.config.output_dir / "reports" / filename

        reporter = ReportGenerator(config=self.config.reporting)
        dossier_path = reporter.generate_report(
            pipeline_result,
            output_path=output_path,
            format=format,
        )

        logger.info(f"Dossier generated: {dossier_path}")

        return dossier_path

    def get_module(self, name: str) -> Optional[ValidationModule]:
        """Get a registered module by name."""
        return self.modules.get(name)

    def list_modules(self) -> List[str]:
        """List all registered modules in execution order."""
        return self.module_order.copy()

    def __repr__(self) -> str:
        return f"<ValidationPipeline: {len(self.modules)} modules, config={self.config.get_hash()}>"
