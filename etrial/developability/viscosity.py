"""
Viscosity Prediction Module for eTrial.

Predicts high-concentration viscosity for monoclonal antibodies (mAbs).

Problem:
    - mAbs at >100 mg/mL can have high viscosity (>20 cP)
    - Prevents subcutaneous injection (need <20 cP)
    - 30-40% of mAbs have viscosity issues

Algorithms:
    1. SAP (Spatial Aggregation Propensity) - Sharma et al. 2014
       - Surface patch analysis
       - Score = net charge × hydrophobic surface area
       - SAP < 20 → low viscosity
       - SAP > 40 → high viscosity

    2. Developability Index - Raybould et al. 2019
       - Combines multiple factors:
         - pI (isoelectric point)
         - PSR (patch surface ratio)
         - Fv charge symmetry
         - Hydrophobic patches
       - DI > 0 → good developability
       - DI < -0.5 → poor developability

    3. Charge Symmetry Parameter - Tomar et al. 2016
       - Asymmetric charge distribution → higher viscosity
       - CSP < 0.3 → low viscosity
       - CSP > 0.5 → high viscosity

References:
    - Sharma VK et al. (2014) "In silico selection of therapeutic antibodies for
      development: Viscosity, clearance, and chemical stability." PNAS 111(52):18601-6.
    - Raybould MIJ et al. (2019) "Five computational developability guidelines for
      therapeutic antibody profiling." PNAS 116(10):4025-4030.
    - Tomar DS et al. (2016) "In-silico prediction of concentration-dependent viscosity
      curves for monoclonal antibody solutions." mAbs 8(1):137-44.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass
import math
from collections import Counter

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    MetricResult,
    Decision,
)


@dataclass
class ViscosityPrediction:
    """Results from viscosity prediction."""
    sap_score: float
    developability_index: float
    charge_symmetry_parameter: float
    predicted_viscosity_cp: float
    risk_level: str  # 'low', 'medium', 'high'


class ViscosityValidationModule(ValidationModule):
    """
    Predicts high-concentration viscosity for biologics.

    Focuses on mAbs at >100 mg/mL concentration.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = "viscosity_prediction"
        self.category = "Developability"
        self.required_modality = "biologic"

        # Thresholds (from literature, adjusted for sequence-based predictions)
        self.sap_threshold_pass = 20.0
        self.sap_threshold_revise = 40.0
        self.di_threshold_pass = 0.0
        self.di_threshold_revise = -0.5
        # CSP thresholds adjusted - some asymmetry is normal in functional proteins
        self.csp_threshold_pass = 0.7   # Most antibodies have some asymmetry
        self.csp_threshold_revise = 0.9  # Very high asymmetry is problematic

    def validate(
        self,
        candidate: TherapeuticCandidate,
        **kwargs
    ) -> ValidationResult:
        """
        Predict viscosity for high-concentration formulation.

        Args:
            candidate: Therapeutic candidate to validate

        Returns:
            ValidationResult with viscosity predictions
        """
        sequence = candidate.get_primary_sequence()

        if not sequence:
            return self._create_result(
                candidate=candidate,
                decision=Decision.PASS,
                metrics=[],
                summary="No sequence provided - skipping viscosity prediction",
            )

        # Check if this is a biologic
        if candidate.is_small_molecule():
            return self._create_result(
                candidate=candidate,
                decision=Decision.PASS,
                metrics=[],
                summary="Viscosity prediction only applicable to biologics",
            )

        metrics = []
        risks = []
        warnings = []

        # Calculate viscosity predictions
        prediction = self._predict_viscosity(sequence)

        # SAP Score
        sap_decision = (Decision.KILL if prediction.sap_score >= self.sap_threshold_revise
                       else Decision.REVISE if prediction.sap_score >= self.sap_threshold_pass
                       else Decision.PASS)

        metrics.append(MetricResult(
            name='sap_score',
            value=prediction.sap_score,
            unit='score',
            threshold_pass=self.sap_threshold_pass,
            threshold_revise=self.sap_threshold_revise,
            decision=sap_decision,
        ))

        if prediction.sap_score >= self.sap_threshold_revise:
            risks.append(f"High SAP score ({prediction.sap_score:.1f}) - likely high viscosity")
        elif prediction.sap_score >= self.sap_threshold_pass:
            warnings.append(f"Moderate SAP score ({prediction.sap_score:.1f}) - may have viscosity issues")

        # Developability Index
        di_decision = (Decision.KILL if prediction.developability_index < self.di_threshold_revise
                      else Decision.REVISE if prediction.developability_index < self.di_threshold_pass
                      else Decision.PASS)

        metrics.append(MetricResult(
            name='developability_index',
            value=prediction.developability_index,
            unit='index',
            threshold_pass=self.di_threshold_pass,
            threshold_revise=self.di_threshold_revise,
            decision=di_decision,
        ))

        if prediction.developability_index < self.di_threshold_revise:
            risks.append(f"Poor developability index ({prediction.developability_index:.2f})")
        elif prediction.developability_index < self.di_threshold_pass:
            warnings.append(f"Suboptimal developability index ({prediction.developability_index:.2f})")

        # Charge Symmetry Parameter
        csp_decision = (Decision.KILL if prediction.charge_symmetry_parameter >= self.csp_threshold_revise
                       else Decision.REVISE if prediction.charge_symmetry_parameter >= self.csp_threshold_pass
                       else Decision.PASS)

        metrics.append(MetricResult(
            name='charge_symmetry_parameter',
            value=prediction.charge_symmetry_parameter,
            unit='csp',
            threshold_pass=self.csp_threshold_pass,
            threshold_revise=self.csp_threshold_revise,
            decision=csp_decision,
        ))

        if prediction.charge_symmetry_parameter >= self.csp_threshold_revise:
            risks.append(f"High charge asymmetry ({prediction.charge_symmetry_parameter:.2f}) - viscosity risk")
        elif prediction.charge_symmetry_parameter >= self.csp_threshold_pass:
            warnings.append(f"Moderate charge asymmetry ({prediction.charge_symmetry_parameter:.2f})")

        # Predicted viscosity
        visc_decision = (Decision.KILL if prediction.predicted_viscosity_cp > 30.0
                        else Decision.REVISE if prediction.predicted_viscosity_cp > 20.0
                        else Decision.PASS)

        metrics.append(MetricResult(
            name='predicted_viscosity',
            value=prediction.predicted_viscosity_cp,
            unit='cP',
            threshold_pass=20.0,
            threshold_revise=30.0,
            decision=visc_decision,
        ))

        if prediction.predicted_viscosity_cp > 30.0:
            risks.append(f"Very high predicted viscosity ({prediction.predicted_viscosity_cp:.1f} cP at 150 mg/mL)")
        elif prediction.predicted_viscosity_cp > 20.0:
            warnings.append(f"High predicted viscosity ({prediction.predicted_viscosity_cp:.1f} cP at 150 mg/mL)")

        # Overall decision
        decision = self._make_decision(prediction)

        summary = (
            f"Viscosity prediction: {prediction.risk_level} risk "
            f"(SAP: {prediction.sap_score:.1f}, DI: {prediction.developability_index:.2f}, "
            f"CSP: {prediction.charge_symmetry_parameter:.2f}, "
            f"η: {prediction.predicted_viscosity_cp:.1f} cP)"
        )

        return self._create_result(
            candidate=candidate,
            decision=decision,
            metrics=metrics,
            summary=summary,
            risks=risks,
            warnings=warnings,
        )

    def _predict_viscosity(self, sequence: str) -> ViscosityPrediction:
        """
        Predict viscosity using SAP, DI, and CSP.

        Args:
            sequence: Protein sequence

        Returns:
            ViscosityPrediction with all parameters
        """
        # Calculate individual components
        sap_score = self._calculate_sap_score(sequence)
        developability_index = self._calculate_developability_index(sequence)
        charge_symmetry = self._calculate_charge_symmetry(sequence)

        # Predict viscosity (empirical correlation)
        # Baseline viscosity at 150 mg/mL
        viscosity = 5.0  # cP baseline

        # SAP contribution (major factor)
        if sap_score > 20:
            viscosity += (sap_score - 20) * 0.5

        # DI contribution (moderate factor)
        if developability_index < 0:
            viscosity += abs(developability_index) * 10

        # CSP contribution (minor factor - only very high asymmetry matters)
        if charge_symmetry > 0.9:  # Only extreme asymmetry
            viscosity += (charge_symmetry - 0.9) * 50
        elif charge_symmetry > 0.7:  # Moderate asymmetry
            viscosity += (charge_symmetry - 0.7) * 10

        # Risk level
        if viscosity < 15 and sap_score < 20 and developability_index > 0:
            risk_level = 'low'
        elif viscosity > 25 or sap_score > 40 or developability_index < -0.5:
            risk_level = 'high'
        else:
            risk_level = 'medium'

        return ViscosityPrediction(
            sap_score=sap_score,
            developability_index=developability_index,
            charge_symmetry_parameter=charge_symmetry,
            predicted_viscosity_cp=viscosity,
            risk_level=risk_level
        )

    def _calculate_sap_score(self, sequence: str) -> float:
        """
        Calculate Spatial Aggregation Propensity (SAP) score.

        SAP = net_charge × hydrophobic_surface_area

        Approximation (sequence-based):
        - Estimate surface exposure based on hydrophobicity windows
        - Calculate net charge
        - Combine into SAP score

        Args:
            sequence: Protein sequence

        Returns:
            SAP score (typically 0-100)
        """
        # Hydrophobicity scale (Kyte-Doolittle)
        hydrophobicity = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
            'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
            'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
            'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
        }

        # Calculate surface hydrophobic patches (sliding window)
        window_size = 9
        max_hydrophobic_patch = 0.0
        hydrophobic_patch_count = 0

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]
            avg_hydrophobicity = sum(hydrophobicity.get(aa, 0) for aa in window) / window_size

            if avg_hydrophobicity > 1.5:  # Threshold for hydrophobic patch
                max_hydrophobic_patch = max(max_hydrophobic_patch, avg_hydrophobicity)
                hydrophobic_patch_count += 1

        # Estimate surface hydrophobic area (arbitrary units)
        surface_hydrophobic_area = max_hydrophobic_patch * hydrophobic_patch_count

        # Calculate net charge at pH 7.0
        positive = sequence.count('K') + sequence.count('R') + sequence.count('H') * 0.5  # His ~50% protonated
        negative = sequence.count('D') + sequence.count('E')
        net_charge = abs(positive - negative)

        # SAP score (normalized)
        sap_score = (net_charge * surface_hydrophobic_area) / len(sequence) * 10

        return sap_score

    def _calculate_developability_index(self, sequence: str) -> float:
        """
        Calculate Developability Index (DI).

        DI combines multiple factors:
        - pI (isoelectric point)
        - PSR (patch surface ratio)
        - Fv charge symmetry
        - Hydrophobic patches

        Higher DI = better developability
        DI > 0 = good
        DI < -0.5 = poor

        Args:
            sequence: Protein sequence

        Returns:
            Developability index (typically -2 to +2)
        """
        di = 0.0

        # Calculate pI (isoelectric point)
        pi = self._calculate_pi(sequence)

        # Prefer pI in range 7-9 for mAbs
        if 7.0 <= pi <= 9.0:
            di += 0.5
        elif pi < 6.0 or pi > 10.0:
            di -= 0.5

        # Calculate patch surface ratio (PSR)
        # Ratio of hydrophobic patches to total surface
        psr = self._calculate_psr(sequence)

        # Prefer PSR < 0.3 (not too many hydrophobic patches)
        if psr < 0.3:
            di += 0.5
        elif psr > 0.5:
            di -= 0.5

        # Check for charge clusters (bad for developability)
        charge_clusters = self._count_charge_clusters(sequence)
        if charge_clusters > 3:
            di -= 0.3

        # Check for high positive charge density (common viscosity issue)
        positive_ratio = (sequence.count('K') + sequence.count('R')) / len(sequence)
        if positive_ratio > 0.15:  # >15% positive residues
            di -= 0.4

        # Check for good distribution of charged residues
        charge_distribution = self._calculate_charge_distribution_uniformity(sequence)
        if charge_distribution > 0.7:  # Well distributed
            di += 0.3

        return di

    def _calculate_charge_symmetry(self, sequence: str) -> float:
        """
        Calculate Charge Symmetry Parameter (CSP).

        Measures asymmetry in charge distribution.
        High asymmetry → higher viscosity.

        CSP < 0.7 = good (some asymmetry is normal)
        CSP > 0.9 = poor (extreme asymmetry)

        Args:
            sequence: Protein sequence

        Returns:
            Charge symmetry parameter (0-1)
        """
        # Split sequence into halves
        mid = len(sequence) // 2
        first_half = sequence[:mid]
        second_half = sequence[mid:]

        # Calculate net charge in each half
        def net_charge(seq):
            positive = seq.count('K') + seq.count('R') + seq.count('H') * 0.5
            negative = seq.count('D') + seq.count('E')
            return positive - negative

        charge_first = net_charge(first_half)
        charge_second = net_charge(second_half)

        # Count total number of charged residues (not net charge)
        def count_charged_residues(seq):
            return seq.count('K') + seq.count('R') + seq.count('D') + seq.count('E') + seq.count('H')

        total_charged = count_charged_residues(sequence)

        # If very few charged residues, asymmetry is not a concern
        if total_charged < 5:
            return 0.0

        # CSP = |charge_first - charge_second| / (|charge_first| + |charge_second|)
        # But handle the case where net charges cancel out
        total_magnitude = abs(charge_first) + abs(charge_second)

        if total_magnitude < 1.0:
            # Charges cancel almost perfectly - check if there's local asymmetry
            # Use absolute charges instead
            positive_first = first_half.count('K') + first_half.count('R') + first_half.count('H')
            positive_second = second_half.count('K') + second_half.count('R') + second_half.count('H')
            total_positive = positive_first + positive_second

            if total_positive < 2:
                return 0.0  # Very few charged residues overall

            # Calculate asymmetry based on positive charge distribution
            csp = abs(positive_first - positive_second) / total_positive
        else:
            csp = abs(charge_first - charge_second) / total_magnitude

        return csp

    def _calculate_pi(self, sequence: str) -> float:
        """
        Estimate isoelectric point (pI).

        Simplified calculation based on charged residues.

        Args:
            sequence: Protein sequence

        Returns:
            Estimated pI
        """
        # Count ionizable residues
        n_term = 1  # N-terminus
        c_term = 1  # C-terminus
        n_lys = sequence.count('K')
        n_arg = sequence.count('R')
        n_his = sequence.count('H')
        n_asp = sequence.count('D')
        n_glu = sequence.count('E')
        n_cys = sequence.count('C')
        n_tyr = sequence.count('Y')

        # pKa values
        pka = {
            'n_term': 9.6,
            'c_term': 2.3,
            'lys': 10.5,
            'arg': 12.5,
            'his': 6.0,
            'asp': 3.9,
            'glu': 4.3,
            'cys': 8.3,
            'tyr': 10.1
        }

        # Simplified pI calculation
        # pI ≈ (sum of basic pKa) - (sum of acidic pKa) / 2 + 7
        basic_pka_sum = (pka['n_term'] + n_lys * pka['lys'] +
                         n_arg * pka['arg'] + n_his * pka['his'])
        acidic_pka_sum = (pka['c_term'] + n_asp * pka['asp'] + n_glu * pka['glu'])

        # Net charge estimate
        n_basic = n_term + n_lys + n_arg + n_his
        n_acidic = c_term + n_asp + n_glu

        if n_basic + n_acidic == 0:
            return 7.0

        # Estimate pI (simplified)
        pi = 7.0 + (n_basic - n_acidic) / (n_basic + n_acidic) * 3.0

        return max(3.0, min(11.0, pi))  # Clamp to reasonable range

    def _calculate_psr(self, sequence: str) -> float:
        """
        Calculate Patch Surface Ratio (PSR).

        Ratio of hydrophobic patches to total surface.

        Args:
            sequence: Protein sequence

        Returns:
            PSR (0-1)
        """
        # Hydrophobicity scale
        hydrophobicity = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
            'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
            'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
            'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
        }

        # Count hydrophobic patches
        window_size = 7
        hydrophobic_patches = 0
        total_windows = 0

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]
            avg_hydrophobicity = sum(hydrophobicity.get(aa, 0) for aa in window) / window_size

            total_windows += 1
            if avg_hydrophobicity > 1.0:  # Threshold for hydrophobic patch
                hydrophobic_patches += 1

        if total_windows == 0:
            return 0.0

        psr = hydrophobic_patches / total_windows
        return psr

    def _count_charge_clusters(self, sequence: str) -> int:
        """
        Count charge clusters (consecutive charged residues).

        Args:
            sequence: Protein sequence

        Returns:
            Number of charge clusters
        """
        clusters = 0
        in_cluster = False
        cluster_length = 0

        for aa in sequence:
            if aa in 'DEKR':  # Charged residues
                cluster_length += 1
                if cluster_length >= 3:  # Cluster threshold
                    if not in_cluster:
                        clusters += 1
                        in_cluster = True
            else:
                cluster_length = 0
                in_cluster = False

        return clusters

    def _calculate_charge_distribution_uniformity(self, sequence: str) -> float:
        """
        Calculate uniformity of charge distribution.

        Higher value = more uniform distribution (better).

        Args:
            sequence: Protein sequence

        Returns:
            Uniformity score (0-1)
        """
        # Split into 10 segments
        n_segments = min(10, len(sequence) // 10)
        if n_segments < 2:
            return 1.0  # Too short to assess

        segment_length = len(sequence) // n_segments
        segment_charges = []

        for i in range(n_segments):
            start = i * segment_length
            end = start + segment_length
            segment = sequence[start:end]

            # Net charge in segment
            positive = segment.count('K') + segment.count('R')
            negative = segment.count('D') + segment.count('E')
            net_charge = abs(positive - negative)

            segment_charges.append(net_charge)

        # Calculate coefficient of variation
        if not segment_charges or sum(segment_charges) == 0:
            return 1.0

        mean_charge = sum(segment_charges) / len(segment_charges)
        if mean_charge == 0:
            return 1.0

        variance = sum((c - mean_charge) ** 2 for c in segment_charges) / len(segment_charges)
        std_dev = math.sqrt(variance)
        cv = std_dev / mean_charge if mean_charge > 0 else 0

        # Convert to uniformity score (1 = uniform, 0 = very non-uniform)
        uniformity = 1.0 / (1.0 + cv)

        return uniformity

    def _make_decision(self, prediction: ViscosityPrediction) -> Decision:
        """
        Make overall decision based on viscosity prediction.

        Args:
            prediction: ViscosityPrediction results

        Returns:
            Decision (PASS/REVISE/KILL)
        """
        # KILL criteria
        if (prediction.sap_score >= self.sap_threshold_revise or
            prediction.developability_index < self.di_threshold_revise or
            prediction.charge_symmetry_parameter >= self.csp_threshold_revise or
            prediction.predicted_viscosity_cp > 30.0):
            return Decision.KILL

        # REVISE criteria
        if (prediction.sap_score >= self.sap_threshold_pass or
            prediction.developability_index < self.di_threshold_pass or
            prediction.charge_symmetry_parameter >= self.csp_threshold_pass or
            prediction.predicted_viscosity_cp > 20.0):
            return Decision.REVISE

        return Decision.PASS
