"""
Base classes and data models for eTrial validation framework.

This module defines the core abstractions used throughout the platform:
- TherapeuticCandidate: Representation of a drug candidate
- ValidationModule: Base class for all validation modules
- ValidationResult: Standardized result format
- Decision: Pass/Revise/Kill decision enumeration
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from pydantic import BaseModel, Field, field_validator
import hashlib
import json


# ============================================================================
# Enumerations
# ============================================================================


class Modality(str, Enum):
    """Therapeutic modality types."""
    ANTIBODY = "antibody"
    NANOBODY = "nanobody"
    PEPTIDE = "peptide"
    SMALL_MOLECULE = "small_molecule"
    PROTEIN = "protein"
    BISPECIFIC = "bispecific"
    ADC = "antibody_drug_conjugate"
    OTHER = "other"


class Decision(str, Enum):
    """Validation decision outcomes."""
    PASS = "PASS"           # Ready for wet-lab validation
    REVISE = "REVISE"       # Addressable issues, recommend optimization
    KILL = "KILL"           # Fundamental flaws, recommend termination
    INFORMATIVE = "INFORMATIVE"  # For non-gating modules


class ValidationStatus(str, Enum):
    """Status of a validation run."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


# ============================================================================
# Data Models
# ============================================================================


class TherapeuticCandidate(BaseModel):
    """
    Representation of a therapeutic candidate to be validated.

    Attributes:
        name: Candidate identifier
        sequence: Amino acid sequence (for biologics/peptides) or SMILES (for small molecules)
        modality: Type of therapeutic
        target: Target protein/pathway name
        structure_file: Optional PDB file path
        metadata: Additional candidate information
    """
    name: str = Field(..., description="Unique identifier for the candidate")
    sequence: Optional[str] = Field(None, description="AA sequence or SMILES string")
    modality: Modality = Field(..., description="Therapeutic modality")
    target: str = Field(..., description="Target name or identifier")
    target_sequence: Optional[str] = Field(None, description="Target protein sequence")
    structure_file: Optional[Path] = Field(None, description="Path to structure file (PDB, etc.)")
    heavy_chain: Optional[str] = Field(None, description="Heavy chain sequence (for antibodies)")
    light_chain: Optional[str] = Field(None, description="Light chain sequence (for antibodies)")
    smiles: Optional[str] = Field(None, description="SMILES string (for small molecules)")
    fasta_file: Optional[Path] = Field(None, description="Path to FASTA file")
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Additional metadata")

    @field_validator('structure_file')
    @classmethod
    def validate_structure_file(cls, v):
        """Validate structure file exists if provided."""
        if v is not None:
            v = Path(v)
            if not v.exists():
                raise ValueError(f"Structure file not found: {v}")
        return v

    def get_hash(self) -> str:
        """Generate unique hash for this candidate."""
        data = {
            "name": self.name,
            "sequence": self.sequence,
            "modality": self.modality.value,
            "target": self.target,
        }
        return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()[:16]

    def is_biologic(self) -> bool:
        """Check if this is a biologic therapeutic."""
        return self.modality in [
            Modality.ANTIBODY,
            Modality.NANOBODY,
            Modality.PROTEIN,
            Modality.BISPECIFIC,
            Modality.ADC,
        ]

    def is_peptide(self) -> bool:
        """Check if this is a peptide therapeutic."""
        return self.modality == Modality.PEPTIDE

    def is_small_molecule(self) -> bool:
        """Check if this is a small molecule therapeutic."""
        return self.modality == Modality.SMALL_MOLECULE

    def get_primary_sequence(self) -> Optional[str]:
        """
        Get the primary sequence for validation.

        Returns the most appropriate sequence based on modality:
        - Antibodies: heavy chain (or concatenated if both present)
        - Small molecules: SMILES string
        - Others: sequence field
        """
        if self.is_small_molecule():
            return self.smiles or self.sequence
        elif self.modality == Modality.ANTIBODY:
            if self.heavy_chain and self.light_chain:
                return self.heavy_chain  # Use heavy chain as primary
            elif self.heavy_chain:
                return self.heavy_chain
            elif self.light_chain:
                return self.light_chain
            else:
                return self.sequence
        else:
            return self.sequence

    @classmethod
    def from_fasta(
        cls,
        fasta_file: Path,
        name: Optional[str] = None,
        modality: Modality = Modality.PROTEIN,
        target: str = "unknown",
        **kwargs
    ) -> "TherapeuticCandidate":
        """
        Create candidate from FASTA file.

        Args:
            fasta_file: Path to FASTA file
            name: Optional name (uses FASTA header if not provided)
            modality: Therapeutic modality
            target: Target name
            **kwargs: Additional arguments

        Returns:
            TherapeuticCandidate instance
        """
        from Bio import SeqIO

        fasta_file = Path(fasta_file)

        if not fasta_file.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

        # Parse FASTA
        records = list(SeqIO.parse(fasta_file, "fasta"))

        if not records:
            raise ValueError(f"No sequences found in {fasta_file}")

        # Use first sequence
        record = records[0]
        sequence = str(record.seq)
        fasta_name = name or record.id

        return cls(
            name=fasta_name,
            sequence=sequence,
            modality=modality,
            target=target,
            fasta_file=fasta_file,
            **kwargs
        )

    @classmethod
    def from_pdb(
        cls,
        pdb_file: Path,
        name: Optional[str] = None,
        modality: Modality = Modality.PROTEIN,
        target: str = "unknown",
        **kwargs
    ) -> "TherapeuticCandidate":
        """
        Create candidate from PDB structure file.

        Args:
            pdb_file: Path to PDB file
            name: Optional name (uses PDB ID if not provided)
            modality: Therapeutic modality
            target: Target name
            **kwargs: Additional arguments

        Returns:
            TherapeuticCandidate instance
        """
        pdb_file = Path(pdb_file)

        if not pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        pdb_name = name or pdb_file.stem

        # Try to extract sequence from PDB
        try:
            from Bio.PDB import PDBParser
            from Bio.PDB.Polypeptide import PPBuilder

            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_name, pdb_file)

            # Build polypeptide
            ppb = PPBuilder()
            peptides = ppb.build_peptides(structure)

            # Get longest chain sequence
            sequences = [str(pp.get_sequence()) for pp in peptides]
            sequence = max(sequences, key=len) if sequences else None

        except ImportError:
            # BioPython not available
            sequence = None

        return cls(
            name=pdb_name,
            sequence=sequence,
            modality=modality,
            target=target,
            structure_file=pdb_file,
            **kwargs
        )

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        name: str,
        target: str = "unknown",
        **kwargs
    ) -> "TherapeuticCandidate":
        """
        Create small molecule candidate from SMILES string.

        Args:
            smiles: SMILES string
            name: Compound name
            target: Target name
            **kwargs: Additional arguments

        Returns:
            TherapeuticCandidate instance
        """
        return cls(
            name=name,
            smiles=smiles,
            sequence=smiles,  # Store in sequence field too for compatibility
            modality=Modality.SMALL_MOLECULE,
            target=target,
            **kwargs
        )


class MetricResult(BaseModel):
    """Individual metric result with value, threshold, and decision."""
    name: str = Field(..., description="Metric name")
    value: Union[float, int, bool, str, None] = Field(..., description="Metric value")
    unit: Optional[str] = Field(None, description="Unit of measurement")
    threshold_pass: Optional[Union[float, int, bool]] = Field(None, description="Pass threshold")
    threshold_revise: Optional[Union[float, int, bool]] = Field(None, description="Revise threshold")
    decision: Optional[Decision] = Field(None, description="Decision for this metric")
    confidence: Optional[float] = Field(None, ge=0.0, le=1.0, description="Prediction confidence")
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Additional metric data")


class ValidationResult(BaseModel):
    """
    Standardized result from a validation module.

    All validation modules return results in this format for consistency.
    """
    module_name: str = Field(..., description="Name of the validation module")
    candidate_hash: str = Field(..., description="Hash of the candidate")
    decision: Decision = Field(..., description="Overall module decision")
    metrics: List[MetricResult] = Field(default_factory=list, description="Individual metric results")
    summary: str = Field(..., description="Human-readable summary")
    recommendations: List[str] = Field(default_factory=list, description="Actionable recommendations")
    risks: List[str] = Field(default_factory=list, description="Identified risks")
    warnings: List[str] = Field(default_factory=list, description="Warnings")
    artifacts: Dict[str, Path] = Field(default_factory=dict, description="Generated artifacts (plots, PDBs, etc.)")
    runtime_seconds: Optional[float] = Field(None, description="Module runtime")
    timestamp: datetime = Field(default_factory=datetime.now, description="Completion timestamp")
    version: str = Field(..., description="Module version")
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Additional result data")

    def get_metric(self, name: str) -> Optional[MetricResult]:
        """Retrieve a specific metric by name."""
        for metric in self.metrics:
            if metric.name == name:
                return metric
        return None

    def get_metric_value(self, name: str) -> Any:
        """Get value of a specific metric."""
        metric = self.get_metric(name)
        return metric.value if metric else None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return self.model_dump()

    def to_json(self) -> str:
        """Convert to JSON string."""
        return self.model_dump_json(indent=2)


class PipelineResult(BaseModel):
    """Complete validation pipeline results."""
    candidate: TherapeuticCandidate = Field(..., description="Validated candidate")
    overall_decision: Decision = Field(..., description="Final decision")
    module_results: Dict[str, ValidationResult] = Field(default_factory=dict, description="Results by module")
    execution_order: List[str] = Field(default_factory=list, description="Module execution order")
    total_runtime_seconds: float = Field(..., description="Total pipeline runtime")
    timestamp: datetime = Field(default_factory=datetime.now, description="Pipeline completion time")
    config_hash: str = Field(..., description="Hash of configuration used")
    audit_trail_id: Optional[str] = Field(None, description="Audit trail identifier")

    def get_module_result(self, module_name: str) -> Optional[ValidationResult]:
        """Get results for a specific module."""
        return self.module_results.get(module_name)

    def get_failed_modules(self) -> List[str]:
        """Get list of modules that resulted in KILL decision."""
        return [
            name for name, result in self.module_results.items()
            if result.decision == Decision.KILL
        ]

    def get_revise_modules(self) -> List[str]:
        """Get list of modules that resulted in REVISE decision."""
        return [
            name for name, result in self.module_results.items()
            if result.decision == Decision.REVISE
        ]

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return self.model_dump()

    def to_json(self, filepath: Optional[Path] = None) -> str:
        """Convert to JSON, optionally save to file."""
        json_str = self.model_dump_json(indent=2)
        if filepath:
            Path(filepath).write_text(json_str)
        return json_str


# ============================================================================
# Base Classes
# ============================================================================


class ValidationModule(ABC):
    """
    Abstract base class for all validation modules.

    All validation modules (structure, specificity, developability, etc.)
    inherit from this class and implement the validate() method.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize validation module.

        Args:
            config: Module-specific configuration dictionary
        """
        self.config = config or {}
        self.name = self.__class__.__name__
        self.version = "0.1.0"
        self._setup()

    def _setup(self) -> None:
        """Module-specific setup (override in subclasses)."""
        pass

    @abstractmethod
    def validate(self, candidate: TherapeuticCandidate, **kwargs) -> ValidationResult:
        """
        Run validation on a therapeutic candidate.

        Args:
            candidate: The therapeutic candidate to validate
            **kwargs: Additional module-specific parameters

        Returns:
            ValidationResult with decision, metrics, and recommendations
        """
        pass

    def _create_result(
        self,
        candidate: TherapeuticCandidate,
        decision: Decision,
        metrics: List[MetricResult],
        summary: str,
        recommendations: Optional[List[str]] = None,
        risks: Optional[List[str]] = None,
        warnings: Optional[List[str]] = None,
        artifacts: Optional[Dict[str, Path]] = None,
        runtime_seconds: Optional[float] = None,
        **kwargs
    ) -> ValidationResult:
        """Helper to create a standardized ValidationResult."""
        return ValidationResult(
            module_name=self.name,
            candidate_hash=candidate.get_hash(),
            decision=decision,
            metrics=metrics,
            summary=summary,
            recommendations=recommendations or [],
            risks=risks or [],
            warnings=warnings or [],
            artifacts=artifacts or {},
            runtime_seconds=runtime_seconds,
            version=self.version,
            metadata=kwargs
        )

    def _apply_thresholds(
        self,
        metric_value: Union[float, int],
        thresholds: Dict[str, List[Optional[Union[float, int]]]]
    ) -> Decision:
        """
        Apply threshold logic to determine decision.

        Args:
            metric_value: The computed metric value
            thresholds: Dict with 'PASS', 'REVISE', 'KILL' keys, each containing [min, max]

        Returns:
            Decision based on which threshold range the value falls into
        """
        for decision_name in ['KILL', 'REVISE', 'PASS']:
            if decision_name in thresholds:
                min_val, max_val = thresholds[decision_name]

                # Check if value is within range
                if min_val is not None and metric_value < min_val:
                    continue
                if max_val is not None and metric_value > max_val:
                    continue

                return Decision[decision_name]

        # Default to INFORMATIVE if no threshold matched
        return Decision.INFORMATIVE

    def __repr__(self) -> str:
        return f"<{self.name} v{self.version}>"


# ============================================================================
# Exceptions
# ============================================================================


class ValidationError(Exception):
    """Base exception for validation errors."""
    pass


class ConfigurationError(ValidationError):
    """Raised when configuration is invalid."""
    pass


class PipelineError(ValidationError):
    """Raised when pipeline execution fails."""
    pass


class ModuleError(ValidationError):
    """Raised when a validation module fails."""
    pass
