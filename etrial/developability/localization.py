"""
Localization Signal Detection Module.

Detects unwanted targeting signals that could affect expression/localization:
    - Signal peptides (secretion)
    - Transmembrane segments
    - Nuclear localization signals (NLS)
    - Mitochondrial targeting sequences

These signals can cause:
    - Expression problems (protein sent to wrong compartment)
    - Low yields (degradation in wrong location)
    - Activity loss (wrong post-translational modifications)

Speed: ~1ms per sequence
Tools: SignalP 6.0 (optional), TMHMM (optional), regex patterns (fallback)
"""

import time
import subprocess
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
from loguru import logger

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    MetricResult,
    Decision,
)


class LocalizationValidationModule(ValidationModule):
    """
    Detect unwanted localization signals in therapeutic sequences.

    Checks for:
    - Signal peptides (N-terminal, typically 15-30 aa)
    - Transmembrane segments (hydrophobic stretches)
    - Nuclear localization signals (NLS)
    - Other targeting sequences
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self.name = "localization_signals"
        self.version = "1.0.0"

    def _setup(self) -> None:
        """Check for SignalP and TMHMM tools."""
        logger.info(f"Setting up {self.name} module")

        # Check for SignalP
        try:
            result = subprocess.run(
                ["signalp", "-h"],
                capture_output=True,
                timeout=5
            )
            self.signalp_available = result.returncode == 0
            if self.signalp_available:
                logger.info("SignalP available for signal peptide prediction")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            self.signalp_available = False
            logger.warning("SignalP not available - using pattern-based detection")

        # Check for TMHMM
        try:
            result = subprocess.run(
                ["tmhmm", "-h"],
                capture_output=True,
                timeout=5
            )
            self.tmhmm_available = result.returncode == 0
            if self.tmhmm_available:
                logger.info("TMHMM available for TM prediction")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            self.tmhmm_available = False
            logger.warning("TMHMM not available - using hydrophobicity-based detection")

    def validate(
        self,
        candidate: TherapeuticCandidate,
        **kwargs
    ) -> ValidationResult:
        """
        Detect localization signals.

        Args:
            candidate: Therapeutic candidate to validate

        Returns:
            ValidationResult with detected signals
        """
        start_time = time.time()
        logger.info(f"Running {self.name} validation for {candidate.name}")

        metrics: List[MetricResult] = []
        warnings: List[str] = []
        recommendations: List[str] = []
        risks: List[str] = []

        # Get sequence
        sequence = candidate.get_primary_sequence()

        if not sequence:
            return self._create_result(
                candidate=candidate,
                decision=Decision.INFORMATIVE,
                metrics=[],
                summary="No sequence available for localization analysis",
                warnings=["Sequence required"],
                runtime_seconds=time.time() - start_time,
            )

        # Skip small molecules
        if candidate.is_small_molecule():
            return self._create_result(
                candidate=candidate,
                decision=Decision.INFORMATIVE,
                metrics=[],
                summary="Small molecule - localization signals not applicable",
                runtime_seconds=time.time() - start_time,
            )

        # ==================================================================
        # STEP 1: Signal Peptide Detection
        # ==================================================================

        if self.signalp_available:
            signal_result = self._predict_with_signalp(sequence)
        else:
            signal_result = self._detect_signal_peptide_pattern(sequence)

        has_signal_peptide = signal_result['has_signal']
        signal_prob = signal_result.get('probability', 0.0)
        signal_position = signal_result.get('cleavage_position', 0)

        metrics.append(MetricResult(
            name="signal_peptide_probability",
            value=signal_prob,
            unit="probability",
            threshold_pass=0.30,  # Flag if > 30% probability
            threshold_revise=0.50,  # Warning if > 50%
            decision=Decision.PASS if signal_prob < 0.30 else (
                Decision.REVISE if signal_prob < 0.50 else Decision.KILL
            ),
            metadata={'method': signal_result['method'], 'position': signal_position}
        ))

        if has_signal_peptide:
            warnings.append(
                f"Signal peptide detected at position 1-{signal_position} "
                f"(probability: {signal_prob:.2f})"
            )
            recommendations.append(
                "Remove signal peptide if cytosolic expression intended, "
                "or verify it matches expression system (bacterial vs mammalian)"
            )
            risks.append(
                "Signal peptide may cause: "
                "(1) secretion into periplasm/ER, (2) expression problems, "
                "(3) incorrect folding if cleaved"
            )

        # ==================================================================
        # STEP 2: Transmembrane Segment Detection
        # ==================================================================

        if self.tmhmm_available:
            tm_result = self._predict_with_tmhmm(sequence)
        else:
            tm_result = self._detect_tm_segments_hydrophobicity(sequence)

        num_tm_segments = tm_result['num_segments']
        tm_segments = tm_result['segments']

        metrics.append(MetricResult(
            name="transmembrane_segments",
            value=num_tm_segments,
            unit="count",
            threshold_pass=0,
            threshold_revise=1,
            decision=Decision.PASS if num_tm_segments == 0 else (
                Decision.REVISE if num_tm_segments <= 1 else Decision.KILL
            ),
            metadata={'method': tm_result['method'], 'segments': tm_segments}
        ))

        if num_tm_segments > 0:
            warnings.append(
                f"Found {num_tm_segments} transmembrane segment(s): {tm_segments}"
            )
            recommendations.append(
                "Remove TM segments if soluble protein intended, "
                "or verify membrane protein expression system compatibility"
            )
            risks.append(
                "TM segments may cause: "
                "(1) membrane insertion, (2) aggregation in cytosol, "
                "(3) inclusion bodies, (4) low expression"
            )

        # ==================================================================
        # STEP 3: Nuclear Localization Signal (NLS)
        # ==================================================================

        nls_result = self._detect_nls_patterns(sequence)
        num_nls = nls_result['num_nls']
        nls_positions = nls_result['positions']

        metrics.append(MetricResult(
            name="nuclear_localization_signals",
            value=num_nls,
            unit="count",
            threshold_pass=0,
            threshold_revise=1,
            decision=Decision.PASS if num_nls == 0 else Decision.REVISE,
            metadata={'positions': nls_positions}
        ))

        if num_nls > 0:
            warnings.append(
                f"Found {num_nls} nuclear localization signal(s) at: {nls_positions}"
            )
            recommendations.append(
                "Remove NLS if cytosolic activity intended"
            )

        # ==================================================================
        # STEP 4: Mitochondrial Targeting Sequence
        # ==================================================================

        mito_result = self._detect_mitochondrial_targeting(sequence)
        has_mito = mito_result['has_targeting']
        mito_prob = mito_result['probability']

        metrics.append(MetricResult(
            name="mitochondrial_targeting_probability",
            value=mito_prob,
            unit="probability",
            threshold_pass=0.30,
            threshold_revise=0.50,
            decision=Decision.PASS if mito_prob < 0.30 else Decision.REVISE,
        ))

        if has_mito:
            warnings.append(
                f"Potential mitochondrial targeting sequence (probability: {mito_prob:.2f})"
            )

        # ==================================================================
        # Overall Decision
        # ==================================================================

        decisions = [m.decision for m in metrics if m.decision]

        if Decision.KILL in decisions:
            overall_decision = Decision.KILL
            summary = (
                f"Multiple localization signals detected. "
                f"Likely unsuitable for intended expression system."
            )
        elif Decision.REVISE in decisions:
            overall_decision = Decision.REVISE
            summary = (
                f"Localization signal(s) detected. "
                f"Verify compatibility with expression system or remove."
            )
        else:
            overall_decision = Decision.PASS
            summary = "No problematic localization signals detected."

        runtime = time.time() - start_time

        result = self._create_result(
            candidate=candidate,
            decision=overall_decision,
            metrics=metrics,
            summary=summary,
            recommendations=recommendations,
            risks=risks,
            warnings=warnings,
            runtime_seconds=runtime,
        )

        logger.info(f"{self.name} validation completed: {overall_decision.value} ({runtime:.3f}s)")

        return result

    def _predict_with_signalp(self, sequence: str) -> Dict[str, Any]:
        """Use SignalP 6.0 for signal peptide prediction."""
        # TODO: Implement SignalP integration when tool is available
        # For now, fall back to pattern-based
        return self._detect_signal_peptide_pattern(sequence)

    def _detect_signal_peptide_pattern(self, sequence: str) -> Dict[str, Any]:
        """
        Detect signal peptide using pattern-based approach.

        Signal peptides typically:
        - 15-30 amino acids at N-terminus
        - Positively charged N-region (1-5 aa)
        - Hydrophobic H-region (7-15 aa)
        - Neutral C-region with cleavage site
        """
        if len(sequence) < 20:
            return {'has_signal': False, 'probability': 0.0, 'method': 'pattern'}

        # Check N-terminal 30 residues
        n_term = sequence[:30]

        # Check for positive charge at start
        n_region = n_term[:5]
        positive_charge = sum(1 for aa in n_region if aa in 'KR')

        # Check for hydrophobic stretch
        hydrophobic = set('LVIFMWA')
        max_hydrophobic_stretch = 0
        current_stretch = 0

        for aa in n_term[5:25]:
            if aa in hydrophobic:
                current_stretch += 1
                max_hydrophobic_stretch = max(max_hydrophobic_stretch, current_stretch)
            else:
                current_stretch = 0

        # Scoring
        probability = 0.0

        if positive_charge >= 1:
            probability += 0.20

        if max_hydrophobic_stretch >= 7:
            probability += 0.40
        elif max_hydrophobic_stretch >= 5:
            probability += 0.20

        # Check for AXA cleavage motif
        import re
        if re.search(r'[AVGS]-[AVGS]-[AVGS]', n_term[15:25]):
            probability += 0.20

        has_signal = probability >= 0.40
        cleavage_position = 15 + max_hydrophobic_stretch if has_signal else 0

        return {
            'has_signal': has_signal,
            'probability': probability,
            'cleavage_position': cleavage_position,
            'method': 'pattern_based'
        }

    def _predict_with_tmhmm(self, sequence: str) -> Dict[str, Any]:
        """Use TMHMM for TM prediction."""
        # TODO: Implement TMHMM integration when tool is available
        return self._detect_tm_segments_hydrophobicity(sequence)

    def _detect_tm_segments_hydrophobicity(self, sequence: str) -> Dict[str, Any]:
        """
        Detect transmembrane segments using Kyte-Doolittle hydrophobicity.

        TM segments are typically 20-25 hydrophobic residues.
        """
        # Kyte-Doolittle scale
        kd_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }

        window_size = 19
        tm_threshold = 1.6  # Hydrophobicity threshold

        segments = []

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            hydrophobicity = sum(kd_scale.get(aa, 0) for aa in window) / window_size

            if hydrophobicity >= tm_threshold:
                # Check if this extends an existing segment
                if segments and i <= segments[-1][1] + 5:
                    segments[-1] = (segments[-1][0], i + window_size)
                else:
                    segments.append((i, i + window_size))

        # Filter segments < 15 aa (likely not real TM)
        segments = [(start, end) for start, end in segments if end - start >= 15]

        segment_strings = [f"{start+1}-{end}" for start, end in segments]

        return {
            'num_segments': len(segments),
            'segments': segment_strings,
            'method': 'kyte_doolittle'
        }

    def _detect_nls_patterns(self, sequence: str) -> Dict[str, Any]:
        """
        Detect nuclear localization signals.

        Classic NLS patterns:
        - Monopartite: KKKRK (4-5 basic residues)
        - Bipartite: KRxxxxxxxxxxKKKK (2 basic + 10aa spacer + 5 basic)
        """
        import re

        positions = []

        # Monopartite NLS: 4+ basic residues in 5 aa window
        for match in re.finditer(r'[KR]{4,5}', sequence):
            positions.append(f"{match.start()+1}-{match.end()}")

        # Bipartite NLS: KR...10-12aa...KKKK
        for match in re.finditer(r'[KR]{2}.{10,12}[KR]{3,4}', sequence):
            positions.append(f"{match.start()+1}-{match.end()}")

        return {
            'num_nls': len(positions),
            'positions': positions
        }

    def _detect_mitochondrial_targeting(self, sequence: str) -> Dict[str, Any]:
        """
        Detect mitochondrial targeting sequence.

        Mitochondrial presequences are typically:
        - 15-50 aa at N-terminus
        - Enriched in Arg, Leu, Ser, Ala
        - Positive charge
        - Form amphipathic helix
        """
        if len(sequence) < 25:
            return {'has_targeting': False, 'probability': 0.0}

        n_term = sequence[:40]

        # Count enriched residues
        enriched = sum(1 for aa in n_term if aa in 'RLSA')
        enrichment = enriched / len(n_term)

        # Positive charge
        positive = sum(1 for aa in n_term if aa in 'KR')
        negative = sum(1 for aa in n_term if aa in 'DE')
        net_charge = positive - negative

        probability = 0.0

        if enrichment >= 0.50:
            probability += 0.30
        elif enrichment >= 0.40:
            probability += 0.15

        if net_charge >= 4:
            probability += 0.30
        elif net_charge >= 2:
            probability += 0.15

        # Lack of acidic residues
        if negative == 0:
            probability += 0.10

        has_targeting = probability >= 0.40

        return {
            'has_targeting': has_targeting,
            'probability': probability
        }
