"""
Low-Complexity Region Detection Module.

Detects problematic sequence regions:
    - Low-complexity regions (homopolymers, simple repeats)
    - Compositional bias (unusual AA distributions)
    - Tandem repeats
    - Proline-rich regions

These regions can cause:
    - Expression problems (ribosome stalling)
    - Disorder/lack of structure
    - Aggregation
    - False positives in screening assays

Speed: ~5ms per sequence
Tools: SEG, CAST (optional), regex patterns (fallback)
"""

import time
from typing import Any, Dict, List, Optional, Tuple
from collections import Counter
from loguru import logger

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    MetricResult,
    Decision,
)


class ComplexityValidationModule(ValidationModule):
    """
    Detect low-complexity and problematic sequence regions.

    Checks for:
    - Homopolymers (AAAAA, LLLLL, etc.)
    - Simple repeats (KAKAKA, TGTGTG)
    - Compositional bias (>50% of one amino acid)
    - Tandem repeats
    - Proline-rich regions
    - Low Shannon entropy
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self.name = "sequence_complexity"
        self.version = "1.0.0"

    def _setup(self) -> None:
        """Initialize complexity analysis tools."""
        logger.info(f"Setting up {self.name} module")

    def validate(
        self,
        candidate: TherapeuticCandidate,
        **kwargs
    ) -> ValidationResult:
        """
        Detect low-complexity regions.

        Args:
            candidate: Therapeutic candidate to validate

        Returns:
            ValidationResult with complexity metrics
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
                summary="No sequence available for complexity analysis",
                warnings=["Sequence required"],
                runtime_seconds=time.time() - start_time,
            )

        # Skip small molecules
        if candidate.is_small_molecule():
            return self._create_result(
                candidate=candidate,
                decision=Decision.INFORMATIVE,
                metrics=[],
                summary="Small molecule - sequence complexity not applicable",
                runtime_seconds=time.time() - start_time,
            )

        # ==================================================================
        # STEP 1: Homopolymer Detection
        # ==================================================================

        homopolymers = self._detect_homopolymers(sequence)
        num_homopolymers = len(homopolymers)

        metrics.append(MetricResult(
            name="homopolymer_count",
            value=num_homopolymers,
            unit="count",
            threshold_pass=0,
            threshold_revise=2,
            decision=Decision.PASS if num_homopolymers == 0 else (
                Decision.REVISE if num_homopolymers <= 2 else Decision.KILL
            ),
            metadata={'homopolymers': homopolymers}
        ))

        if num_homopolymers > 0:
            warnings.append(
                f"Found {num_homopolymers} homopolymer(s): {homopolymers}"
            )
            recommendations.append(
                "Break up homopolymers by introducing synonymous substitutions "
                "(maintain similar properties)"
            )
            risks.append(
                "Homopolymers can cause: "
                "(1) expression problems, (2) disorder, (3) aggregation"
            )

        # ==================================================================
        # STEP 2: Simple Repeat Detection
        # ==================================================================

        repeats = self._detect_simple_repeats(sequence)
        num_repeats = len(repeats)

        metrics.append(MetricResult(
            name="simple_repeat_count",
            value=num_repeats,
            unit="count",
            threshold_pass=1,
            threshold_revise=3,
            decision=Decision.PASS if num_repeats <= 1 else (
                Decision.REVISE if num_repeats <= 3 else Decision.KILL
            ),
            metadata={'repeats': repeats}
        ))

        if num_repeats > 1:
            warnings.append(
                f"Found {num_repeats} simple repeat(s): {repeats[:5]}"  # Show first 5
            )

        # ==================================================================
        # STEP 3: Compositional Bias
        # ==================================================================

        bias_result = self._detect_compositional_bias(sequence)
        max_bias = bias_result['max_fraction']
        biased_aa = bias_result['biased_aa']

        metrics.append(MetricResult(
            name="max_aa_fraction",
            value=max_bias,
            unit="fraction",
            threshold_pass=0.30,
            threshold_revise=0.40,
            decision=Decision.PASS if max_bias < 0.30 else (
                Decision.REVISE if max_bias < 0.40 else Decision.KILL
            ),
            metadata={'amino_acid': biased_aa}
        ))

        if max_bias >= 0.30:
            warnings.append(
                f"Compositional bias: {biased_aa} = {max_bias:.1%} of sequence"
            )
            risks.append(
                f"High {biased_aa} content may indicate disorder or low complexity"
            )

        # ==================================================================
        # STEP 4: Shannon Entropy (Sequence Complexity)
        # ==================================================================

        entropy = self._calculate_shannon_entropy(sequence)

        # Normalize entropy to 0-1 (max entropy = log2(20) ≈ 4.32 for 20 amino acids)
        normalized_entropy = entropy / 4.32

        metrics.append(MetricResult(
            name="sequence_entropy",
            value=normalized_entropy,
            unit="score",
            threshold_pass=0.70,  # Low entropy = low complexity
            threshold_revise=0.50,
            decision=Decision.PASS if normalized_entropy >= 0.70 else (
                Decision.REVISE if normalized_entropy >= 0.50 else Decision.KILL
            ),
        ))

        if normalized_entropy < 0.70:
            warnings.append(
                f"Low sequence complexity (entropy: {normalized_entropy:.2f})"
            )
            recommendations.append(
                "Increase sequence diversity by introducing substitutions"
            )

        # ==================================================================
        # STEP 5: Proline-Rich Regions
        # ==================================================================

        proline_regions = self._detect_proline_rich_regions(sequence)
        num_proline_regions = len(proline_regions)

        metrics.append(MetricResult(
            name="proline_rich_regions",
            value=num_proline_regions,
            unit="count",
            threshold_pass=0,
            threshold_revise=1,
            decision=Decision.PASS if num_proline_regions == 0 else Decision.REVISE,
            metadata={'regions': proline_regions}
        ))

        if num_proline_regions > 0:
            warnings.append(
                f"Found {num_proline_regions} proline-rich region(s): {proline_regions}"
            )
            risks.append(
                "Proline-rich regions can cause rigidity and disorder"
            )

        # ==================================================================
        # STEP 6: Tandem Repeats
        # ==================================================================

        tandem_repeats = self._detect_tandem_repeats(sequence)
        num_tandem = len(tandem_repeats)

        metrics.append(MetricResult(
            name="tandem_repeats",
            value=num_tandem,
            unit="count",
            threshold_pass=0,
            threshold_revise=1,
            decision=Decision.PASS if num_tandem == 0 else Decision.REVISE,
            metadata={'repeats': tandem_repeats}
        ))

        if num_tandem > 0:
            warnings.append(
                f"Found {num_tandem} tandem repeat(s): {tandem_repeats}"
            )

        # ==================================================================
        # STEP 7: Poly-Basic/Poly-Acidic Regions
        # ==================================================================

        charged_regions = self._detect_poly_charged_regions(sequence)
        num_charged = len(charged_regions)

        metrics.append(MetricResult(
            name="poly_charged_regions",
            value=num_charged,
            unit="count",
            threshold_pass=1,
            threshold_revise=3,
            decision=Decision.PASS if num_charged <= 1 else (
                Decision.REVISE if num_charged <= 3 else Decision.KILL
            ),
            metadata={'regions': charged_regions}
        ))

        if num_charged > 1:
            warnings.append(
                f"Found {num_charged} poly-charged region(s): {charged_regions}"
            )
            risks.append(
                "Poly-charged regions can cause: "
                "(1) cytotoxicity (poly-basic), (2) aggregation, (3) disorder"
            )

        # ==================================================================
        # Overall Decision
        # ==================================================================

        decisions = [m.decision for m in metrics if m.decision]

        if Decision.KILL in decisions:
            overall_decision = Decision.KILL
            summary = (
                "Multiple low-complexity or problematic regions detected. "
                "Sequence redesign strongly recommended."
            )
        elif Decision.REVISE in decisions:
            overall_decision = Decision.REVISE
            summary = (
                "Low-complexity or problematic regions detected. "
                "Review and consider sequence modifications."
            )
        else:
            overall_decision = Decision.PASS
            summary = "Sequence complexity is acceptable."

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

    def _detect_homopolymers(self, sequence: str, min_length: int = 5) -> List[str]:
        """
        Detect homopolymers (runs of same amino acid).

        Args:
            sequence: Protein sequence
            min_length: Minimum homopolymer length to report

        Returns:
            List of homopolymer strings with positions
        """
        import re

        homopolymers = []

        for aa in set(sequence):
            pattern = f'{aa}{{{min_length},}}'
            for match in re.finditer(pattern, sequence):
                homopolymers.append(
                    f"{aa}×{match.end()-match.start()} at {match.start()+1}-{match.end()}"
                )

        return homopolymers

    def _detect_simple_repeats(self, sequence: str, min_repeats: int = 3) -> List[str]:
        """
        Detect simple repeats (XY)n, (XYZ)n, etc.

        Args:
            sequence: Protein sequence
            min_repeats: Minimum number of repeat units

        Returns:
            List of repeat strings with positions
        """
        import re

        repeats = []

        # Check 2-mers
        for i in range(len(sequence) - 5):
            motif = sequence[i:i+2]
            pattern = f'({re.escape(motif)}){{{min_repeats},}}'
            match = re.match(pattern, sequence[i:])
            if match:
                length = match.end()
                repeats.append(
                    f"({motif})×{length//2} at {i+1}-{i+length}"
                )

        # Check 3-mers
        for i in range(len(sequence) - 8):
            motif = sequence[i:i+3]
            pattern = f'({re.escape(motif)}){{{min_repeats},}}'
            match = re.match(pattern, sequence[i:])
            if match:
                length = match.end()
                repeats.append(
                    f"({motif})×{length//3} at {i+1}-{i+length}"
                )

        # Deduplicate overlapping repeats
        return list(set(repeats))[:10]  # Return max 10

    def _detect_compositional_bias(self, sequence: str) -> Dict[str, Any]:
        """
        Detect compositional bias (overrepresented amino acids).

        Returns:
            dict with max_fraction and biased_aa
        """
        counts = Counter(sequence)
        total = len(sequence)

        max_count = max(counts.values())
        max_fraction = max_count / total

        biased_aa = [aa for aa, count in counts.items() if count == max_count][0]

        return {
            'max_fraction': max_fraction,
            'biased_aa': biased_aa,
            'composition': dict(counts)
        }

    def _calculate_shannon_entropy(self, sequence: str) -> float:
        """
        Calculate Shannon entropy of sequence.

        Higher entropy = more complex/diverse sequence
        """
        import math

        counts = Counter(sequence)
        total = len(sequence)

        entropy = 0.0
        for count in counts.values():
            p = count / total
            entropy -= p * math.log2(p)

        return entropy

    def _detect_proline_rich_regions(
        self,
        sequence: str,
        window: int = 20,
        threshold: float = 0.30
    ) -> List[str]:
        """
        Detect proline-rich regions (>30% Pro in 20aa window).

        Returns:
            List of proline-rich region positions
        """
        regions = []

        for i in range(len(sequence) - window + 1):
            seg = sequence[i:i+window]
            pro_fraction = seg.count('P') / window

            if pro_fraction >= threshold:
                regions.append(f"{i+1}-{i+window} ({pro_fraction:.0%} Pro)")

        # Merge overlapping regions
        if not regions:
            return []

        merged = [regions[0]]
        for region in regions[1:]:
            last_end = int(merged[-1].split('-')[1].split()[0])
            curr_start = int(region.split('-')[0])

            if curr_start <= last_end + 5:
                # Merge
                merged[-1] = f"{merged[-1].split('-')[0]}-{region.split('-')[1]}"
            else:
                merged.append(region)

        return merged

    def _detect_tandem_repeats(
        self,
        sequence: str,
        min_unit_size: int = 4,
        min_repeats: int = 2
    ) -> List[str]:
        """
        Detect tandem repeats of longer motifs (4+ aa).

        Returns:
            List of tandem repeat strings
        """
        import re

        repeats = []

        # Check 4-10 aa motifs
        for unit_size in range(min_unit_size, 11):
            for i in range(len(sequence) - unit_size * min_repeats):
                motif = sequence[i:i+unit_size]
                pattern = f'({re.escape(motif)}){{{min_repeats},}}'
                match = re.match(pattern, sequence[i:])

                if match:
                    length = match.end()
                    num_repeats = length // unit_size
                    repeats.append(
                        f"({motif})×{num_repeats} at {i+1}-{i+length}"
                    )

        return list(set(repeats))[:5]  # Max 5

    def _detect_poly_charged_regions(
        self,
        sequence: str,
        window: int = 10,
        threshold: float = 0.60
    ) -> List[str]:
        """
        Detect regions with >60% charged residues.

        Returns:
            List of poly-charged region positions
        """
        regions = []

        for i in range(len(sequence) - window + 1):
            seg = sequence[i:i+window]
            charged = sum(1 for aa in seg if aa in 'DEKR')
            charged_fraction = charged / window

            if charged_fraction >= threshold:
                # Classify as poly-basic or poly-acidic
                basic = sum(1 for aa in seg if aa in 'KR')
                acidic = sum(1 for aa in seg if aa in 'DE')

                if basic > acidic:
                    charge_type = "poly-basic"
                else:
                    charge_type = "poly-acidic"

                regions.append(
                    f"{charge_type} {i+1}-{i+window} ({charged_fraction:.0%} charged)"
                )

        # Merge overlapping
        if not regions:
            return []

        merged = []
        for region in regions:
            if not merged:
                merged.append(region)
            else:
                # Simple deduplication (not perfect merging)
                if region not in merged:
                    merged.append(region)

        return merged[:10]  # Max 10
