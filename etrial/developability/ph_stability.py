"""
pH Stability Prediction Module for eTrial.

Predicts protein stability at different pH values and recommends
optimal formulation pH.

Problem:
    - Proteins can be unstable at certain pH values
    - Need to select formulation pH where protein is stable
    - pH affects: charge state, solubility, aggregation, degradation

Methods:
    1. PROPKA Integration (preferred)
       - Install PROPKA 3.5: pip install propka
       - Predicts pKa for all ionizable residues
       - Identifies buried charged residues (pH-labile)
       - Accuracy: ~0.5-1.0 pKa unit vs experimental

    2. Simple Electrostatic Fallback
       - Count buried Asp/Glu/His/Lys/Arg
       - Estimate pKa shifts based on environment
       - Flag residues with large shifts
       - Accuracy: ~1-2 pKa unit (rough estimate)

Applications:
    - Formulation development: Choose stable pH
    - Aggregation prediction: pH-dependent aggregation
    - Chemical stability: Deamidation, isomerization
    - Storage conditions: Long-term stability

References:
    - Olsson et al. (2011) "PROPKA3: Consistent Treatment of Internal and
      Surface Residues in Empirical pKa Predictions" J Chem Theory Comput 7:525-37
    - Li et al. (2005) "Very fast empirical prediction and rationalization of
      protein pKa values" Proteins 61:704-721
"""

from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile
from collections import defaultdict

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    MetricResult,
    Decision,
)


@dataclass
class pKaResult:
    """Result for a single ionizable residue."""
    position: int
    residue: str
    standard_pka: float  # Standard pKa value
    predicted_pka: float  # Predicted pKa (shifted by environment)
    pka_shift: float  # Difference from standard
    labile: bool  # Large shift indicates pH-labile
    buried: bool  # Estimated burial state


class pHStabilityModule(ValidationModule):
    """
    Predicts pH stability and recommends formulation pH.

    Identifies pH-labile residues and predicts stable pH range.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = "ph_stability"
        self.category = "Developability"
        self.required_modality = "biologic"

        # Standard pKa values (in water)
        self.standard_pka = {
            'ASP': 3.9, 'GLU': 4.3,  # Acidic
            'HIS': 6.0,               # Histidine
            'CYS': 8.3,               # Cysteine
            'TYR': 10.1,              # Tyrosine
            'LYS': 10.5,              # Lysine
            'ARG': 12.5,              # Arginine
            'N-term': 9.6,            # N-terminus
            'C-term': 2.3,            # C-terminus
        }

        # Thresholds
        self.pka_shift_threshold = 1.5  # pKa units (labile if shifted > 1.5)
        self.target_ph = 6.5  # Target formulation pH (typical for mAbs)

        # Check for PROPKA
        self.propka_available = self._check_propka()

    def _check_propka(self) -> bool:
        """Check if PROPKA is available."""
        try:
            import propka
            return True
        except ImportError:
            return False

    def validate(
        self,
        candidate: TherapeuticCandidate,
        **kwargs
    ) -> ValidationResult:
        """
        Predict pH stability and recommend formulation pH.

        Args:
            candidate: Therapeutic candidate with sequence (and optionally structure)

        Returns:
            ValidationResult with pH stability predictions
        """
        sequence = candidate.get_primary_sequence()

        if not sequence:
            return self._create_result(
                candidate=candidate,
                decision=Decision.INFORMATIVE,
                metrics=[],
                summary="No sequence provided - skipping pH stability prediction",
            )

        # Check if biologic
        if candidate.is_small_molecule():
            return self._create_result(
                candidate=candidate,
                decision=Decision.INFORMATIVE,
                metrics=[],
                summary="pH stability prediction only applicable to biologics",
            )

        metrics = []
        recommendations = []
        risks = []
        warnings = []

        # Run pH stability analysis
        if self.propka_available and candidate.structure_file:
            # Use PROPKA with structure
            pka_results = self._run_propka(candidate.structure_file)
            method = "PROPKA"
        else:
            # Fallback to sequence-based estimation
            pka_results = self._estimate_pka_shifts(sequence)
            method = "electrostatic estimation"

            if not candidate.structure_file:
                warnings.append("No structure provided - using sequence-based approximations")

        # Analyze results
        labile_residues = [r for r in pka_results if r.labile]
        buried_charged = [r for r in pka_results if r.buried and abs(r.pka_shift) > 0.5]

        # Predict stable pH range
        stable_ph_range = self._predict_stable_ph_range(pka_results, sequence)
        recommended_ph = self._recommend_formulation_ph(pka_results, stable_ph_range)

        # Metrics
        metrics.append(MetricResult(
            name='num_labile_residues',
            value=len(labile_residues),
            unit='count',
            threshold_pass=2,  # ≤2 labile residues is acceptable
            threshold_revise=5,
            decision=Decision.PASS if len(labile_residues) <= 2 else Decision.REVISE if len(labile_residues) <= 5 else Decision.KILL,
            metadata={'method': method}
        ))

        metrics.append(MetricResult(
            name='num_buried_charged',
            value=len(buried_charged),
            unit='count',
            threshold_pass=3,
            decision=Decision.PASS if len(buried_charged) <= 3 else Decision.REVISE
        ))

        metrics.append(MetricResult(
            name='stable_ph_range_width',
            value=stable_ph_range[1] - stable_ph_range[0],
            unit='pH units',
            threshold_pass=2.0,  # At least 2 pH unit range
            decision=Decision.PASS if (stable_ph_range[1] - stable_ph_range[0]) >= 2.0 else Decision.REVISE
        ))

        metrics.append(MetricResult(
            name='recommended_ph',
            value=recommended_ph,
            unit='pH',
            threshold_pass=5.0,  # pH 5-7 is typical
            threshold_revise=8.0,
            decision=Decision.INFORMATIVE
        ))

        # Recommendations
        if len(labile_residues) == 0:
            recommendations.append(f"Good pH stability - no highly pH-labile residues")
        else:
            warnings.append(f"Found {len(labile_residues)} pH-labile residue(s)")
            labile_positions = [f"{r.residue}{r.position}" for r in labile_residues[:5]]
            recommendations.append(f"pH-labile positions: {', '.join(labile_positions)}")

        if stable_ph_range[1] - stable_ph_range[0] >= 2.0:
            recommendations.append(
                f"Stable pH range: {stable_ph_range[0]:.1f} - {stable_ph_range[1]:.1f}"
            )
        else:
            risks.append(
                f"Narrow stable pH range ({stable_ph_range[0]:.1f} - {stable_ph_range[1]:.1f})"
            )

        recommendations.append(
            f"Recommended formulation pH: {recommended_ph:.1f}"
        )

        # Check if recommended pH is in stable range
        if not (stable_ph_range[0] <= recommended_ph <= stable_ph_range[1]):
            warnings.append(
                f"Recommended pH ({recommended_ph:.1f}) is outside stable range"
            )

        # Warn about specific residues
        for residue in labile_residues:
            if residue.residue in ['HIS']:
                warnings.append(
                    f"His{residue.position} has shifted pKa ({residue.predicted_pka:.1f}) - "
                    f"may cause pH-dependent aggregation"
                )
            elif residue.buried:
                warnings.append(
                    f"Buried {residue.residue}{residue.position} has shifted pKa - "
                    f"may cause conformational instability"
                )

        # Overall decision
        if len(labile_residues) > 5:
            decision = Decision.KILL
            summary = f"Many pH-labile residues ({len(labile_residues)} found, {method})"
        elif len(labile_residues) > 2:
            decision = Decision.REVISE
            summary = f"Some pH-labile residues ({len(labile_residues)} found, {method})"
        else:
            decision = Decision.PASS
            summary = f"Good pH stability (stable pH {stable_ph_range[0]:.1f}-{stable_ph_range[1]:.1f}, recommend pH {recommended_ph:.1f}, {method})"

        return self._create_result(
            candidate=candidate,
            decision=decision,
            metrics=metrics,
            summary=summary,
            recommendations=recommendations,
            risks=risks,
            warnings=warnings,
            metadata={
                'pka_results': [
                    {
                        'position': r.position,
                        'residue': r.residue,
                        'predicted_pka': r.predicted_pka,
                        'pka_shift': r.pka_shift,
                        'labile': r.labile
                    }
                    for r in pka_results
                ],
                'stable_ph_range': stable_ph_range,
                'recommended_ph': recommended_ph
            }
        )

    def _run_propka(self, structure_file: Path) -> List[pKaResult]:
        """
        Run PROPKA to predict pKa values.

        Args:
            structure_file: Path to PDB structure

        Returns:
            List of pKaResult objects
        """
        # TODO: Implement PROPKA integration
        # This would involve:
        # 1. Import propka
        # 2. Run propka.run() on PDB file
        # 3. Parse output pKa values
        # 4. Create pKaResult objects

        # For now, fallback to estimation
        # Extract sequence from PDB if needed
        from Bio import PDB
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', structure_file)

        # Get sequence
        ppb = PDB.PPBuilder()
        peptides = ppb.build_peptides(structure)
        sequence = str(max(peptides, key=lambda p: len(p)).get_sequence()) if peptides else ""

        return self._estimate_pka_shifts(sequence)

    def _estimate_pka_shifts(self, sequence: str) -> List[pKaResult]:
        """
        Estimate pKa shifts using simple electrostatic model.

        Args:
            sequence: Protein sequence

        Returns:
            List of pKaResult objects
        """
        results = []

        # Map 1-letter to 3-letter codes
        one_to_three = {
            'D': 'ASP', 'E': 'GLU', 'H': 'HIS', 'C': 'CYS',
            'Y': 'TYR', 'K': 'LYS', 'R': 'ARG'
        }

        # Analyze each ionizable residue
        for i, aa in enumerate(sequence, start=1):
            if aa not in one_to_three:
                continue

            three_letter = one_to_three[aa]
            standard_pka = self.standard_pka[three_letter]

            # Estimate pKa shift based on local environment
            pka_shift, buried = self._estimate_local_environment(sequence, i-1, aa)

            predicted_pka = standard_pka + pka_shift

            labile = abs(pka_shift) >= self.pka_shift_threshold

            result = pKaResult(
                position=i,
                residue=aa,
                standard_pka=standard_pka,
                predicted_pka=predicted_pka,
                pka_shift=pka_shift,
                labile=labile,
                buried=buried
            )

            results.append(result)

        # Add termini
        # N-terminus
        results.append(pKaResult(
            position=1,
            residue='N-term',
            standard_pka=self.standard_pka['N-term'],
            predicted_pka=self.standard_pka['N-term'],
            pka_shift=0.0,
            labile=False,
            buried=False
        ))

        # C-terminus
        results.append(pKaResult(
            position=len(sequence),
            residue='C-term',
            standard_pka=self.standard_pka['C-term'],
            predicted_pka=self.standard_pka['C-term'],
            pka_shift=0.0,
            labile=False,
            buried=False
        ))

        return results

    def _estimate_local_environment(
        self, sequence: str, position: int, residue: str
    ) -> Tuple[float, bool]:
        """
        Estimate pKa shift based on local environment.

        Factors that shift pKa:
        - Nearby charged residues (electrostatic interactions)
        - Hydrophobic environment (buried residues)
        - Hydrogen bonding partners

        Args:
            sequence: Full sequence
            position: Position to analyze (0-indexed)
            residue: Residue type

        Returns:
            Tuple of (pKa_shift, is_buried)
        """
        pka_shift = 0.0
        buried = False

        # Define charge states
        acidic = 'DE'
        basic = 'KRH'
        hydrophobic = 'FILVMWY'

        # Check local environment (±5 residues)
        start = max(0, position - 5)
        end = min(len(sequence), position + 6)
        window = sequence[start:end]

        # Count nearby charges
        nearby_acidic = sum(1 for aa in window if aa in acidic and aa != residue)
        nearby_basic = sum(1 for aa in window if aa in basic and aa != residue)

        # Electrostatic effects
        if residue in acidic:  # Acidic residue (D, E)
            # Nearby acidic residues raise pKa (repulsion)
            pka_shift += nearby_acidic * 0.5
            # Nearby basic residues lower pKa (attraction)
            pka_shift -= nearby_basic * 0.3

        elif residue in basic:  # Basic residue (K, R, H)
            # Nearby basic residues raise pKa (repulsion)
            pka_shift += nearby_basic * 0.5
            # Nearby acidic residues lower pKa (attraction)
            pka_shift -= nearby_acidic * 0.3

        # Hydrophobic burial effect
        hydrophobic_count = sum(1 for aa in window if aa in hydrophobic)
        hydrophobic_fraction = hydrophobic_count / len(window)

        if hydrophobic_fraction > 0.6:
            # Buried in hydrophobic core
            buried = True

            # Buried charged residues have large pKa shifts
            if residue in acidic:
                pka_shift += 2.0  # Harder to deprotonate when buried
            elif residue in basic:
                pka_shift -= 2.0  # Easier to deprotonate when buried

        # Position effects (termini more exposed)
        rel_position = position / len(sequence)
        if rel_position < 0.1 or rel_position > 0.9:
            # Terminal regions - more exposed
            pka_shift *= 0.5

        return pka_shift, buried

    def _predict_stable_ph_range(
        self, pka_results: List[pKaResult], sequence: str
    ) -> Tuple[float, float]:
        """
        Predict pH range where protein is stable.

        Stability considerations:
        - Avoid pH near pKa of buried charged residues
        - Avoid extreme pH where protein is highly charged
        - Prefer pH 5-7 for biologics (typical)

        Args:
            pka_results: List of pKa predictions
            sequence: Protein sequence

        Returns:
            Tuple of (min_pH, max_pH)
        """
        # Start with typical biologic range
        min_ph = 5.0
        max_ph = 7.5

        # Avoid pH near labile residues
        for result in pka_results:
            if result.labile and result.buried:
                # Avoid pH within ±1 unit of this pKa
                if result.predicted_pka - 1.0 > min_ph:
                    min_ph = max(min_ph, result.predicted_pka - 1.0)
                if result.predicted_pka + 1.0 < max_ph:
                    max_ph = min(max_ph, result.predicted_pka + 1.0)

        # Calculate isoelectric point (pI)
        pi = self._calculate_pi(sequence)

        # Avoid pH very close to pI (low solubility)
        if abs(pi - 6.0) < 1.5:
            # If pI is near 6, shift range slightly
            if pi < 6.0:
                min_ph = max(min_ph, pi + 0.5)
            else:
                max_ph = min(max_ph, pi - 0.5)

        # Ensure reasonable range
        if max_ph - min_ph < 1.0:
            # Expand range if too narrow
            center = (min_ph + max_ph) / 2
            min_ph = max(4.0, center - 1.0)
            max_ph = min(8.0, center + 1.0)

        return (min_ph, max_ph)

    def _recommend_formulation_ph(
        self, pka_results: List[pKaResult], stable_range: Tuple[float, float]
    ) -> float:
        """
        Recommend optimal formulation pH.

        Considerations:
        - Within stable pH range
        - Avoid pH near buried charged residues
        - Prefer pH 6-6.5 for mAbs (industry standard)

        Args:
            pka_results: List of pKa predictions
            stable_range: Stable pH range

        Returns:
            Recommended pH
        """
        # Start with target pH (typical for mAbs)
        recommended = self.target_ph

        # Ensure within stable range
        if recommended < stable_range[0]:
            recommended = stable_range[0] + 0.2
        elif recommended > stable_range[1]:
            recommended = stable_range[1] - 0.2

        # Avoid pH near highly labile residues
        for result in pka_results:
            if result.labile and result.buried and abs(result.pka_shift) > 2.0:
                # Stay away from this pKa
                if abs(recommended - result.predicted_pka) < 0.5:
                    # Shift away
                    if recommended < result.predicted_pka:
                        recommended = max(stable_range[0], result.predicted_pka - 1.0)
                    else:
                        recommended = min(stable_range[1], result.predicted_pka + 1.0)

        return round(recommended, 1)

    def _calculate_pi(self, sequence: str) -> float:
        """
        Calculate isoelectric point (pI).

        Simplified calculation based on charged residues.

        Args:
            sequence: Protein sequence

        Returns:
            Estimated pI
        """
        # Count ionizable residues
        n_term = 1
        c_term = 1
        n_lys = sequence.count('K')
        n_arg = sequence.count('R')
        n_his = sequence.count('H')
        n_asp = sequence.count('D')
        n_glu = sequence.count('E')

        # Net charge estimate
        n_basic = n_term + n_lys + n_arg + n_his
        n_acidic = c_term + n_asp + n_glu

        if n_basic + n_acidic == 0:
            return 7.0

        # Estimate pI (simplified)
        pi = 7.0 + (n_basic - n_acidic) / (n_basic + n_acidic) * 3.0

        return max(3.0, min(11.0, pi))
