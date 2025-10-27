"""
Configuration management for eTrial.

Handles loading, validating, and accessing configuration from YAML files.
"""

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Optional, Union
import yaml
from pydantic import BaseModel, Field, field_validator


class Config(BaseModel):
    """
    Main configuration class for eTrial pipeline.

    Loads from YAML files and provides validated access to settings.
    """

    # Global settings
    version: str = Field(default="0.1.0")
    random_seed: int = Field(default=42)
    n_jobs: int = Field(default=-1)
    gpu_enabled: bool = Field(default=True)
    output_dir: Path = Field(default=Path("outputs"))
    cache_dir: Path = Field(default=Path(".cache"))
    log_level: str = Field(default="INFO")

    # Module enable/disable
    modules: Dict[str, bool] = Field(default_factory=dict)

    # Module-specific configurations
    structure: Dict[str, Any] = Field(default_factory=dict)
    specificity: Dict[str, Any] = Field(default_factory=dict)
    developability: Dict[str, Any] = Field(default_factory=dict)
    immunogenicity: Dict[str, Any] = Field(default_factory=dict)
    safety: Dict[str, Any] = Field(default_factory=dict)
    pkpd: Dict[str, Any] = Field(default_factory=dict)
    benchmarking: Dict[str, Any] = Field(default_factory=dict)

    # Reporting configuration
    reporting: Dict[str, Any] = Field(default_factory=dict)

    # Audit trail configuration
    audit: Dict[str, Any] = Field(default_factory=dict)

    # Thresholds (loaded separately)
    _thresholds: Optional[Dict[str, Any]] = None

    class Config:
        arbitrary_types_allowed = True

    @field_validator('output_dir', 'cache_dir')
    @classmethod
    def create_directory(cls, v):
        """Create directory if it doesn't exist."""
        if v is not None:
            v = Path(v)
            v.mkdir(parents=True, exist_ok=True)
        return v

    @classmethod
    def from_yaml(cls, config_path: Union[str, Path]) -> "Config":
        """
        Load configuration from YAML file.

        Args:
            config_path: Path to YAML configuration file

        Returns:
            Config instance

        Raises:
            FileNotFoundError: If config file doesn't exist
            yaml.YAMLError: If YAML is malformed
        """
        config_path = Path(config_path)

        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, 'r') as f:
            data = yaml.safe_load(f)

        # Handle nested structure
        global_settings = data.pop('global', {})

        # Merge global settings with other data
        config_data = {**global_settings, **data}

        return cls(**config_data)

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> "Config":
        """Create configuration from dictionary."""
        return cls(**config_dict)

    @classmethod
    def default(cls) -> "Config":
        """Create default configuration."""
        return cls()

    def load_thresholds(self, thresholds_path: Union[str, Path]) -> None:
        """
        Load decision thresholds from YAML file.

        Args:
            thresholds_path: Path to thresholds YAML file
        """
        thresholds_path = Path(thresholds_path)

        if not thresholds_path.exists():
            raise FileNotFoundError(f"Thresholds file not found: {thresholds_path}")

        with open(thresholds_path, 'r') as f:
            self._thresholds = yaml.safe_load(f)

    def get_thresholds(self, module_name: str) -> Dict[str, Any]:
        """
        Get decision thresholds for a specific module.

        Args:
            module_name: Name of the validation module

        Returns:
            Dictionary of thresholds for the module
        """
        if self._thresholds is None:
            raise ValueError("Thresholds not loaded. Call load_thresholds() first.")

        return self._thresholds.get(module_name, {})

    def get_module_config(self, module_name: str) -> Dict[str, Any]:
        """
        Get configuration for a specific module.

        Args:
            module_name: Name of the module (e.g., 'structure', 'specificity')

        Returns:
            Module configuration dictionary
        """
        return getattr(self, module_name, {})

    def is_module_enabled(self, module_name: str) -> bool:
        """
        Check if a module is enabled.

        Args:
            module_name: Name of the module

        Returns:
            True if enabled, False otherwise
        """
        return self.modules.get(module_name, True)

    def get_hash(self) -> str:
        """
        Generate hash of configuration for reproducibility.

        Returns:
            SHA256 hash (first 16 chars)
        """
        # Create simplified dict for hashing (exclude paths and mutable settings)
        hash_data = {
            'version': self.version,
            'random_seed': self.random_seed,
            'modules': self.modules,
            'structure': self.structure,
            'specificity': self.specificity,
            'developability': self.developability,
            'immunogenicity': self.immunogenicity,
            'safety': self.safety,
            'pkpd': self.pkpd,
        }

        json_str = json.dumps(hash_data, sort_keys=True)
        return hashlib.sha256(json_str.encode()).hexdigest()[:16]

    def to_yaml(self, output_path: Union[str, Path]) -> None:
        """
        Save configuration to YAML file.

        Args:
            output_path: Path to save YAML file
        """
        output_path = Path(output_path)

        # Convert to dict, handling Path objects
        data = self.model_dump()
        for key in ['output_dir', 'cache_dir']:
            if key in data and data[key] is not None:
                data[key] = str(data[key])

        with open(output_path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return self.model_dump()

    def __repr__(self) -> str:
        enabled_modules = [k for k, v in self.modules.items() if v]
        return f"<Config: {len(enabled_modules)} modules enabled, hash={self.get_hash()}>"


class ThresholdManager:
    """
    Helper class for managing and applying decision thresholds.
    """

    def __init__(self, thresholds: Dict[str, Any]):
        """
        Initialize threshold manager.

        Args:
            thresholds: Dictionary of thresholds loaded from YAML
        """
        self.thresholds = thresholds

    def get_module_thresholds(self, module_name: str) -> Dict[str, Any]:
        """Get thresholds for a specific module."""
        # Handle underscore vs hyphen in module names
        module_name_normalized = module_name.replace('-', '_').replace(' ', '_').lower()

        for key in self.thresholds:
            if key.replace('-', '_').replace(' ', '_').lower() == module_name_normalized:
                return self.thresholds[key]

        return {}

    def get_metric_thresholds(
        self,
        module_name: str,
        metric_name: str
    ) -> Dict[str, list]:
        """
        Get thresholds for a specific metric.

        Returns dict with 'PASS', 'REVISE', 'KILL' keys containing [min, max] ranges.
        """
        module_thresholds = self.get_module_thresholds(module_name)

        # Handle nested structure
        for section_name, section_data in module_thresholds.items():
            if isinstance(section_data, dict) and metric_name in section_data:
                return section_data[metric_name]

        return {}

    def apply_threshold(
        self,
        module_name: str,
        metric_name: str,
        value: Union[float, int, bool]
    ) -> str:
        """
        Apply threshold logic to determine decision.

        Args:
            module_name: Name of the module
            metric_name: Name of the metric
            value: Metric value to evaluate

        Returns:
            Decision string: 'PASS', 'REVISE', or 'KILL'
        """
        thresholds = self.get_metric_thresholds(module_name, metric_name)

        if not thresholds:
            return "INFORMATIVE"

        # Check each decision tier
        for decision in ['KILL', 'REVISE', 'PASS']:
            if decision in thresholds:
                range_def = thresholds[decision]

                # Handle boolean thresholds
                if isinstance(range_def, bool):
                    if value == range_def:
                        return decision
                    continue

                # Handle range thresholds [min, max]
                if isinstance(range_def, list) and len(range_def) == 2:
                    min_val, max_val = range_def

                    # Check if within range
                    if min_val is not None and value < min_val:
                        continue
                    if max_val is not None and value > max_val:
                        continue

                    return decision

        return "INFORMATIVE"

    def get_overall_decision(
        self,
        module_decisions: Dict[str, str]
    ) -> str:
        """
        Determine overall decision based on module decisions.

        Args:
            module_decisions: Dict of module_name -> decision

        Returns:
            Overall decision: 'PASS', 'REVISE', or 'KILL'
        """
        decision_logic = self.thresholds.get('decision', {})

        # Check for KILL conditions
        kill_conditions = decision_logic.get('KILL', {}).get('conditions', [])

        # If any critical module is KILL
        if 'KILL' in module_decisions.values():
            return 'KILL'

        # If too many REVISE
        revise_count = list(module_decisions.values()).count('REVISE')
        max_revise = decision_logic.get('REVISE', {}).get('max_kill_modules', 3)

        if revise_count > max_revise:
            return 'KILL'

        # If any REVISE
        if 'REVISE' in module_decisions.values():
            return 'REVISE'

        # All PASS
        return 'PASS'
