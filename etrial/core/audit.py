"""
Audit trail and reproducibility tracking for eTrial.

Ensures all validation runs are fully reproducible by tracking:
- Input hashes
- Configuration versions
- Model versions
- Random seeds
- Runtime environment
- Execution timestamps
"""

import hashlib
import json
import platform
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import uuid

from loguru import logger


@dataclass
class AuditEntry:
    """Single entry in the audit trail."""

    entry_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    timestamp: datetime = field(default_factory=datetime.now)
    event_type: str = ""  # e.g., "validation_started", "module_completed"
    module_name: Optional[str] = None
    candidate_hash: Optional[str] = None
    data: Dict[str, Any] = field(default_factory=dict)
    input_hashes: Dict[str, str] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "entry_id": self.entry_id,
            "timestamp": self.timestamp.isoformat(),
            "event_type": self.event_type,
            "module_name": self.module_name,
            "candidate_hash": self.candidate_hash,
            "data": self.data,
            "input_hashes": self.input_hashes,
            "metadata": self.metadata,
        }


@dataclass
class EnvironmentSnapshot:
    """Snapshot of the execution environment."""

    python_version: str = field(default_factory=lambda: sys.version)
    platform_name: str = field(default_factory=lambda: platform.platform())
    platform_version: str = field(default_factory=lambda: platform.version())
    processor: str = field(default_factory=lambda: platform.processor())
    hostname: str = field(default_factory=lambda: platform.node())
    etrial_version: str = "0.1.0"
    installed_packages: Dict[str, str] = field(default_factory=dict)
    gpu_info: Optional[Dict[str, Any]] = None
    timestamp: datetime = field(default_factory=datetime.now)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "python_version": self.python_version,
            "platform": self.platform_name,
            "platform_version": self.platform_version,
            "processor": self.processor,
            "hostname": self.hostname,
            "etrial_version": self.etrial_version,
            "installed_packages": self.installed_packages,
            "gpu_info": self.gpu_info,
            "timestamp": self.timestamp.isoformat(),
        }

    @classmethod
    def capture(cls, include_packages: bool = True) -> "EnvironmentSnapshot":
        """
        Capture current environment snapshot.

        Args:
            include_packages: Whether to include installed package versions

        Returns:
            EnvironmentSnapshot
        """
        snapshot = cls()

        if include_packages:
            try:
                import pkg_resources
                snapshot.installed_packages = {
                    pkg.key: pkg.version
                    for pkg in pkg_resources.working_set
                }
            except ImportError:
                logger.warning("pkg_resources not available, skipping package inventory")

        # Try to get GPU info
        try:
            import torch
            if torch.cuda.is_available():
                snapshot.gpu_info = {
                    "available": True,
                    "count": torch.cuda.device_count(),
                    "devices": [
                        {
                            "name": torch.cuda.get_device_name(i),
                            "capability": torch.cuda.get_device_capability(i),
                        }
                        for i in range(torch.cuda.device_count())
                    ],
                }
        except ImportError:
            snapshot.gpu_info = {"available": False, "reason": "torch not installed"}

        return snapshot


class AuditTrail:
    """
    Manages audit trail for validation runs.

    Ensures full reproducibility by tracking all inputs, configurations,
    model versions, and execution details.
    """

    def __init__(
        self,
        output_dir: Optional[Path] = None,
        auto_save: bool = True,
        include_environment: bool = True,
    ):
        """
        Initialize audit trail.

        Args:
            output_dir: Directory to save audit logs
            auto_save: Automatically save entries as they're added
            include_environment: Capture environment snapshot
        """
        self.trail_id = str(uuid.uuid4())
        self.created_at = datetime.now()
        self.output_dir = Path(output_dir) if output_dir else Path("outputs/audit")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.auto_save = auto_save

        self.entries: List[AuditEntry] = []
        self.environment: Optional[EnvironmentSnapshot] = None
        self.config_hash: Optional[str] = None
        self.model_versions: Dict[str, str] = {}

        if include_environment:
            self.environment = EnvironmentSnapshot.capture()

        logger.info(f"Audit trail initialized: {self.trail_id}")

    def log_event(
        self,
        event_type: str,
        module_name: Optional[str] = None,
        candidate_hash: Optional[str] = None,
        data: Optional[Dict[str, Any]] = None,
        input_hashes: Optional[Dict[str, str]] = None,
        **kwargs
    ) -> AuditEntry:
        """
        Log an event to the audit trail.

        Args:
            event_type: Type of event (e.g., 'validation_started', 'module_completed')
            module_name: Name of the module (if applicable)
            candidate_hash: Hash of the candidate being validated
            data: Event-specific data
            input_hashes: Hashes of input data for reproducibility
            **kwargs: Additional metadata

        Returns:
            Created AuditEntry
        """
        entry = AuditEntry(
            event_type=event_type,
            module_name=module_name,
            candidate_hash=candidate_hash,
            data=data or {},
            input_hashes=input_hashes or {},
            metadata=kwargs,
        )

        self.entries.append(entry)

        if self.auto_save:
            self._save_entry(entry)

        logger.debug(f"Audit event logged: {event_type} [{module_name}]")
        return entry

    def set_config_hash(self, config_hash: str) -> None:
        """Set the configuration hash for this run."""
        self.config_hash = config_hash
        self.log_event("config_set", data={"config_hash": config_hash})

    def register_model_version(self, model_name: str, version: str) -> None:
        """
        Register a model version being used.

        Args:
            model_name: Name of the model (e.g., 'alphafold', 'netmhciipan')
            version: Version string
        """
        self.model_versions[model_name] = version
        self.log_event(
            "model_registered",
            data={"model_name": model_name, "version": version}
        )

    def hash_file(self, filepath: Union[str, Path]) -> str:
        """
        Generate SHA256 hash of a file.

        Args:
            filepath: Path to file

        Returns:
            Hex digest of file hash
        """
        filepath = Path(filepath)
        hasher = hashlib.sha256()

        with open(filepath, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                hasher.update(chunk)

        file_hash = hasher.hexdigest()

        self.log_event(
            "file_hashed",
            data={"filepath": str(filepath), "hash": file_hash}
        )

        return file_hash

    def hash_data(self, data: Any) -> str:
        """
        Generate hash of arbitrary data.

        Args:
            data: Data to hash (must be JSON-serializable)

        Returns:
            Hex digest of data hash
        """
        json_str = json.dumps(data, sort_keys=True)
        return hashlib.sha256(json_str.encode()).hexdigest()

    def _save_entry(self, entry: AuditEntry) -> None:
        """Save a single entry to disk."""
        entry_file = self.output_dir / f"entry_{entry.entry_id}.json"

        with open(entry_file, 'w') as f:
            json.dump(entry.to_dict(), f, indent=2)

    def save(self, filepath: Optional[Path] = None) -> Path:
        """
        Save complete audit trail to file.

        Args:
            filepath: Optional custom filepath (defaults to timestamped file)

        Returns:
            Path to saved file
        """
        if filepath is None:
            timestamp = self.created_at.strftime("%Y%m%d_%H%M%S")
            filepath = self.output_dir / f"audit_trail_{timestamp}_{self.trail_id[:8]}.json"

        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        audit_data = {
            "trail_id": self.trail_id,
            "created_at": self.created_at.isoformat(),
            "config_hash": self.config_hash,
            "model_versions": self.model_versions,
            "environment": self.environment.to_dict() if self.environment else None,
            "entries": [entry.to_dict() for entry in self.entries],
        }

        with open(filepath, 'w') as f:
            json.dump(audit_data, f, indent=2)

        logger.info(f"Audit trail saved: {filepath}")
        return filepath

    @classmethod
    def load(cls, filepath: Union[str, Path]) -> "AuditTrail":
        """
        Load audit trail from file.

        Args:
            filepath: Path to saved audit trail

        Returns:
            Loaded AuditTrail instance
        """
        filepath = Path(filepath)

        with open(filepath, 'r') as f:
            data = json.load(f)

        trail = cls(auto_save=False, include_environment=False)
        trail.trail_id = data["trail_id"]
        trail.created_at = datetime.fromisoformat(data["created_at"])
        trail.config_hash = data.get("config_hash")
        trail.model_versions = data.get("model_versions", {})

        # Reconstruct entries
        for entry_data in data.get("entries", []):
            entry = AuditEntry(
                entry_id=entry_data["entry_id"],
                timestamp=datetime.fromisoformat(entry_data["timestamp"]),
                event_type=entry_data["event_type"],
                module_name=entry_data.get("module_name"),
                candidate_hash=entry_data.get("candidate_hash"),
                data=entry_data.get("data", {}),
                input_hashes=entry_data.get("input_hashes", {}),
                metadata=entry_data.get("metadata", {}),
            )
            trail.entries.append(entry)

        logger.info(f"Audit trail loaded: {filepath}")
        return trail

    def get_events_by_type(self, event_type: str) -> List[AuditEntry]:
        """Get all entries of a specific event type."""
        return [e for e in self.entries if e.event_type == event_type]

    def get_events_by_module(self, module_name: str) -> List[AuditEntry]:
        """Get all entries for a specific module."""
        return [e for e in self.entries if e.module_name == module_name]

    def get_timeline(self) -> List[Dict[str, Any]]:
        """
        Get chronological timeline of events.

        Returns:
            List of events sorted by timestamp
        """
        timeline = []
        for entry in sorted(self.entries, key=lambda e: e.timestamp):
            timeline.append({
                "timestamp": entry.timestamp.isoformat(),
                "event": entry.event_type,
                "module": entry.module_name,
                "data": entry.data,
            })
        return timeline

    def verify_reproducibility(self, other: "AuditTrail") -> Dict[str, bool]:
        """
        Verify if another run is reproducible from this one.

        Args:
            other: Another AuditTrail to compare against

        Returns:
            Dict with reproducibility checks
        """
        return {
            "config_match": self.config_hash == other.config_hash,
            "model_versions_match": self.model_versions == other.model_versions,
            "environment_match": (
                self.environment.to_dict() if self.environment else None
            ) == (
                other.environment.to_dict() if other.environment else None
            ),
        }

    def __repr__(self) -> str:
        return f"<AuditTrail {self.trail_id[:8]}: {len(self.entries)} entries>"
