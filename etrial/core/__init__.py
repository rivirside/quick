"""Core framework for eTrial validation pipeline."""

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    Decision,
    Modality,
)
from etrial.core.config import Config
from etrial.core.pipeline import ValidationPipeline
from etrial.core.audit import AuditTrail
from etrial.core.reporting import ReportGenerator

__all__ = [
    "TherapeuticCandidate",
    "ValidationModule",
    "ValidationResult",
    "Decision",
    "Modality",
    "Config",
    "ValidationPipeline",
    "AuditTrail",
    "ReportGenerator",
]
