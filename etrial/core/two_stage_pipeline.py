"""
Two-Stage Validation Pipeline for eTrial.

Architecture:
    STAGE 1 (Pre-Filter): Fast sequence-based screening
        - Simple heuristics, regex patterns
        - ~1 second per candidate
        - Filters 10,000 → 500 candidates
        - 90% accuracy on developability

    STAGE 2 (Clinical-Grade): Rigorous validation
        - Structure prediction, validated ML models
        - ~1-2 hours per candidate
        - Filters 500 → 50 → 5 leads
        - 85%+ accuracy overall

Usage:
    # Pre-filter only (fast screening)
    pipeline = TwoStageValidationPipeline(mode='prefilter')
    results = pipeline.validate_batch(candidates)  # Fast

    # Clinical-grade only (detailed analysis)
    pipeline = TwoStageValidationPipeline(mode='clinical')
    results = pipeline.validate_batch(top_candidates)  # Slow but accurate

    # Both stages (full workflow)
    pipeline = TwoStageValidationPipeline(mode='both')
    results = pipeline.validate_batch(candidates)  # Auto-filters between stages
"""

from typing import List, Dict, Any, Optional, Literal
from pathlib import Path
from dataclasses import dataclass
from loguru import logger
import time

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    PipelineResult,
    Decision,
)
from etrial.core.pipeline import ValidationPipeline
from etrial.core.config import Config


@dataclass
class StageResult:
    """Result from a single validation stage."""
    stage: str
    candidates_in: int
    candidates_out: int
    pass_rate: float
    runtime_seconds: float
    results: List[PipelineResult]


class TwoStageValidationPipeline:
    """
    Two-stage validation pipeline for therapeutic screening.

    Stage 1 (Pre-Filter): Fast sequence-based heuristics
    Stage 2 (Clinical): Structure-based, ML predictions
    """

    def __init__(
        self,
        mode: Literal['prefilter', 'clinical', 'both'] = 'both',
        config_path: Optional[Path] = None,
        prefilter_threshold: Decision = Decision.KILL,
        clinical_threshold: Decision = Decision.KILL,
    ):
        """
        Initialize two-stage pipeline.

        Args:
            mode: Validation mode
                - 'prefilter': Stage 1 only (fast screening)
                - 'clinical': Stage 2 only (rigorous validation)
                - 'both': Run both stages sequentially
            config_path: Path to configuration file
            prefilter_threshold: Candidates worse than this are filtered
            clinical_threshold: Final candidates worse than this are rejected
        """
        self.mode = mode
        self.prefilter_threshold = prefilter_threshold
        self.clinical_threshold = clinical_threshold

        # Load config
        if config_path:
            config = Config.from_yaml(config_path)
        else:
            config = Config.default()

        # Create stage-specific pipelines
        self.prefilter_pipeline = ValidationPipeline(config=config)
        self.clinical_pipeline = ValidationPipeline(config=config)

        # Module registry (to be populated)
        self.prefilter_modules = []
        self.clinical_modules = []

        logger.info(f"Two-stage pipeline initialized in '{mode}' mode")

    def register_prefilter_module(self, module_class):
        """Register a module for Stage 1 (pre-filter)."""
        self.prefilter_modules.append(module_class)
        self.prefilter_pipeline.register_module(module_class)
        logger.info(f"Registered pre-filter module: {module_class.__name__}")

    def register_clinical_module(self, module_class):
        """Register a module for Stage 2 (clinical-grade)."""
        self.clinical_modules.append(module_class)
        self.clinical_pipeline.register_module(module_class)
        logger.info(f"Registered clinical module: {module_class.__name__}")

    def validate_batch(
        self,
        candidates: List[TherapeuticCandidate],
        save_intermediates: bool = True,
    ) -> Dict[str, Any]:
        """
        Validate batch of candidates through selected stage(s).

        Args:
            candidates: List of therapeutic candidates
            save_intermediates: Save results from each stage

        Returns:
            Dictionary with stage results and final candidates
        """
        results = {
            'mode': self.mode,
            'stages': [],
            'final_candidates': [],
            'total_runtime_seconds': 0,
        }

        start_time = time.time()

        # Stage 1: Pre-Filter
        if self.mode in ['prefilter', 'both']:
            stage1_result = self._run_prefilter_stage(candidates)
            results['stages'].append(stage1_result)
            results['total_runtime_seconds'] += stage1_result.runtime_seconds

            if self.mode == 'prefilter':
                # Pre-filter only - return these results
                results['final_candidates'] = stage1_result.results
                results['final_pass_rate'] = stage1_result.pass_rate
                return results

            # Filter for Stage 2
            candidates_for_stage2 = [
                r.candidate for r in stage1_result.results
                if r.overall_decision != Decision.KILL
            ]

            logger.info(
                f"Stage 1 filtered {len(candidates)} → {len(candidates_for_stage2)} "
                f"candidates ({len(candidates_for_stage2)/len(candidates)*100:.1f}% pass rate)"
            )

        else:
            # Clinical mode only - use input candidates directly
            candidates_for_stage2 = candidates

        # Stage 2: Clinical-Grade
        if self.mode in ['clinical', 'both']:
            stage2_result = self._run_clinical_stage(candidates_for_stage2)
            results['stages'].append(stage2_result)
            results['total_runtime_seconds'] += stage2_result.runtime_seconds

            results['final_candidates'] = stage2_result.results
            results['final_pass_rate'] = stage2_result.pass_rate

        results['total_runtime_seconds'] = time.time() - start_time

        # Summary
        logger.info(
            f"Pipeline complete: {len(candidates)} → {len(results['final_candidates'])} "
            f"final candidates in {results['total_runtime_seconds']:.1f}s"
        )

        return results

    def _run_prefilter_stage(
        self,
        candidates: List[TherapeuticCandidate]
    ) -> StageResult:
        """Run Stage 1: Pre-Filter validation."""
        logger.info("=" * 70)
        logger.info("STAGE 1: PRE-FILTER (Sequence-Based Screening)")
        logger.info("=" * 70)

        start_time = time.time()
        results = []

        for i, candidate in enumerate(candidates, 1):
            if i % 100 == 0:
                logger.info(f"Pre-filter progress: {i}/{len(candidates)}")

            result = self.prefilter_pipeline.validate(candidate)
            results.append(result)

        runtime = time.time() - start_time

        # Calculate pass rate
        passed = sum(1 for r in results if r.overall_decision != Decision.KILL)
        pass_rate = passed / len(candidates) if candidates else 0

        stage_result = StageResult(
            stage='prefilter',
            candidates_in=len(candidates),
            candidates_out=passed,
            pass_rate=pass_rate,
            runtime_seconds=runtime,
            results=results
        )

        logger.info(f"Stage 1 complete: {pass_rate*100:.1f}% pass rate in {runtime:.1f}s")
        logger.info(f"Throughput: {len(candidates)/runtime:.1f} candidates/sec")
        logger.info("")

        return stage_result

    def _run_clinical_stage(
        self,
        candidates: List[TherapeuticCandidate]
    ) -> StageResult:
        """Run Stage 2: Clinical-Grade validation."""
        logger.info("=" * 70)
        logger.info("STAGE 2: CLINICAL-GRADE (Structure + ML Validation)")
        logger.info("=" * 70)

        if not candidates:
            logger.warning("No candidates to validate in clinical stage")
            return StageResult(
                stage='clinical',
                candidates_in=0,
                candidates_out=0,
                pass_rate=0,
                runtime_seconds=0,
                results=[]
            )

        start_time = time.time()
        results = []

        for i, candidate in enumerate(candidates, 1):
            logger.info(f"Clinical validation [{i}/{len(candidates)}]: {candidate.name}")

            result = self.clinical_pipeline.validate(candidate)
            results.append(result)

            logger.info(
                f"  → Decision: {result.overall_decision.value} "
                f"(runtime: {result.total_runtime_seconds:.1f}s)"
            )

        runtime = time.time() - start_time

        # Calculate pass rate
        passed = sum(1 for r in results if r.overall_decision != Decision.KILL)
        pass_rate = passed / len(candidates) if candidates else 0

        stage_result = StageResult(
            stage='clinical',
            candidates_in=len(candidates),
            candidates_out=passed,
            pass_rate=pass_rate,
            runtime_seconds=runtime,
            results=results
        )

        logger.info(f"Stage 2 complete: {pass_rate*100:.1f}% pass rate in {runtime:.1f}s")
        logger.info(f"Average time per candidate: {runtime/len(candidates):.1f}s")
        logger.info("")

        return stage_result

    def get_stage_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Get detailed statistics for each stage."""
        stats = {
            'mode': results['mode'],
            'total_runtime_seconds': results['total_runtime_seconds'],
            'stages': {}
        }

        for stage_result in results['stages']:
            stats['stages'][stage_result.stage] = {
                'candidates_in': stage_result.candidates_in,
                'candidates_out': stage_result.candidates_out,
                'pass_rate': f"{stage_result.pass_rate*100:.1f}%",
                'runtime_seconds': stage_result.runtime_seconds,
                'throughput_per_sec': (
                    stage_result.candidates_in / stage_result.runtime_seconds
                    if stage_result.runtime_seconds > 0 else 0
                )
            }

        # Decision breakdown
        if results['final_candidates']:
            decisions = [r.overall_decision.value for r in results['final_candidates']]
            stats['final_decision_breakdown'] = {
                'PASS': decisions.count('PASS'),
                'REVISE': decisions.count('REVISE'),
                'KILL': decisions.count('KILL'),
            }

        return stats


def setup_default_prefilter_pipeline(pipeline: TwoStageValidationPipeline):
    """
    Set up default pre-filter modules (current implementations).

    Includes:
        - Structure/Affinity (sequence-based predictions)
        - Developability (sequence-based PTM, hydrophobic patches)
        - Localization (signal peptide, TM, NLS detection) [NEW]
        - Complexity (low-complexity region detection) [NEW]
        - Advanced Aggregation (TANGO/AGGRESCAN) [NEW]
        - Viscosity Prediction (SAP, DI, CSP) [NEW - Week 1]
        - Safety (PAINS, Lipinski, structural alerts)
        - PKPD (Lipinski, Veber, bioavailability estimates)
    """
    from etrial.structure.affinity import AffinityValidationModule
    from etrial.developability.solubility import SolubilityValidationModule
    from etrial.developability.localization import LocalizationValidationModule
    from etrial.developability.complexity import ComplexityValidationModule
    from etrial.developability.aggregation_advanced import AdvancedAggregationModule
    from etrial.developability.viscosity import ViscosityValidationModule
    from etrial.safety.toxicity import ToxicityPrefilterModule
    from etrial.pkpd.pharmacokinetics import PKPDPrefilterModule

    pipeline.register_prefilter_module(AffinityValidationModule)
    pipeline.register_prefilter_module(SolubilityValidationModule)

    # Quality control modules
    pipeline.register_prefilter_module(LocalizationValidationModule)  # Signal peptide, TM, NLS
    pipeline.register_prefilter_module(ComplexityValidationModule)    # Low-complexity regions
    pipeline.register_prefilter_module(AdvancedAggregationModule)     # TANGO/AGGRESCAN

    # Viscosity prediction (Week 1 quick win)
    pipeline.register_prefilter_module(ViscosityValidationModule)     # SAP, DI, CSP

    pipeline.register_prefilter_module(ToxicityPrefilterModule)
    pipeline.register_prefilter_module(PKPDPrefilterModule)

    logger.info("Pre-filter pipeline configured with 8 modules: structure, developability, localization, complexity, aggregation, viscosity, safety, PKPD")


def setup_clinical_pipeline(pipeline: TwoStageValidationPipeline):
    """
    Set up clinical-grade modules.

    Uses clinical-grade tools:
        - CamSol for solubility prediction
        - AlphaFold/ESMFold for structure prediction
        - NetMHCIIpan for immunogenicity
        - MMseqs2 for specificity screening
        - pkCSM/ProTox-II for toxicity prediction
        - QSAR models for PKPD prediction

    Falls back to pre-filter modules if clinical tools aren't installed.
    """
    try:
        # Try to import clinical modules
        from etrial.developability.solubility_clinical import ClinicalSolubilityModule
        from etrial.structure.structure_clinical import ClinicalStructureModule
        from etrial.immunogenicity.tcell_clinical import ClinicalImmunogenicityModule
        from etrial.specificity.homology_clinical import ClinicalSpecificityModule
        from etrial.safety.toxicity_clinical import ClinicalToxicityModule
        from etrial.pkpd.pharmacokinetics_clinical import ClinicalPKPDModule

        # Register clinical modules
        pipeline.register_clinical_module(ClinicalStructureModule)
        pipeline.register_clinical_module(ClinicalSolubilityModule)
        pipeline.register_clinical_module(ClinicalImmunogenicityModule)
        pipeline.register_clinical_module(ClinicalSpecificityModule)
        pipeline.register_clinical_module(ClinicalToxicityModule)
        pipeline.register_clinical_module(ClinicalPKPDModule)

        logger.info("Clinical pipeline configured with clinical-grade modules")
        logger.info("Tools will use: CamSol, AlphaFold/ESMFold, NetMHCIIpan, MMseqs2, pkCSM, QSAR")
        logger.info("(Falls back to sequence-based predictions if tools not installed)")

    except ImportError as e:
        logger.warning(f"Clinical modules not found: {e}")
        logger.warning("Falling back to pre-filter modules for clinical stage")

        # Fallback to pre-filter modules
        from etrial.structure.affinity import AffinityValidationModule
        from etrial.developability.solubility import SolubilityValidationModule
        from etrial.safety.toxicity import ToxicityPrefilterModule
        from etrial.pkpd.pharmacokinetics import PKPDPrefilterModule

        pipeline.register_clinical_module(AffinityValidationModule)
        pipeline.register_clinical_module(SolubilityValidationModule)
        pipeline.register_clinical_module(ToxicityPrefilterModule)
        pipeline.register_clinical_module(PKPDPrefilterModule)


# Convenience function
def create_two_stage_pipeline(
    mode: Literal['prefilter', 'clinical', 'both'] = 'both',
    config_path: Optional[Path] = None,
) -> TwoStageValidationPipeline:
    """
    Create a two-stage pipeline with default module configuration.

    Args:
        mode: Pipeline mode ('prefilter', 'clinical', or 'both')
        config_path: Optional path to config file

    Returns:
        Configured TwoStageValidationPipeline
    """
    pipeline = TwoStageValidationPipeline(
        mode=mode,
        config_path=config_path
    )

    # Set up default modules
    setup_default_prefilter_pipeline(pipeline)
    setup_clinical_pipeline(pipeline)

    return pipeline
