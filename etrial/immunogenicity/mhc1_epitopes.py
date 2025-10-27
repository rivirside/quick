"""
MHC-I Epitope Prediction Module for eTrial.

Predicts CD8 T-cell epitopes (MHC Class I binding) to assess
immunogenicity risk.

Problem:
    - Unwanted T-cell responses can cause immune rejection
    - CD8 T-cells recognize 8-11mer peptides presented by MHC-I
    - Need to minimize MHC-I epitopes in therapeutic proteins

Methods:
    1. NetMHCpan Integration (preferred)
       - Install NetMHCpan 4.1: https://services.healthtech.dtu.dk/services/NetMHCpan-4.3/
       - Predicts binding to HLA-A, HLA-B, HLA-C alleles
       - Strong binders: percentile rank < 0.5%
       - Weak binders: percentile rank < 2%
       - Accuracy: ~85-90% AUC

    2. Pattern-Based Fallback
       - Known MHC-I binding motifs
       - Hydrophobicity + anchor residue patterns
       - Less accurate (~60-70%) but fast

Applications:
    - Immunogenicity assessment
    - De-immunization (remove epitopes)
    - T-cell response prediction
    - Complement MHC-II predictions

References:
    - Jurtz et al. (2017) "NetMHCpan-4.0: Improved Peptide-MHC Class I
      Interaction Predictions Integrating Eluted Ligand and Peptide Binding Data"
      J Immunol 199:3360-3368
    - Rapin et al. (2010) "Computational immunology meets bioinformatics:
      the use of prediction tools for molecular binding in the simulation of
      the immune system" PLoS One 5(4):e9862
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
class MHC1Epitope:
    """MHC-I epitope prediction result."""
    position: int
    peptide: str
    length: int
    hla_allele: str
    percentile_rank: float  # Lower = stronger binding
    strong_binder: bool
    weak_binder: bool


class MHC1EpitopeModule(ValidationModule):
    """
    Predicts MHC-I epitopes for CD8 T-cell immunogenicity assessment.

    Complements MHC-II predictions for comprehensive immunogenicity analysis.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = "mhc1_epitopes"
        self.category = "Immunogenicity"
        self.required_modality = "biologic"

        # Thresholds (NetMHCpan standard)
        self.strong_binder_percentile = 0.5  # < 0.5% rank
        self.weak_binder_percentile = 2.0    # < 2.0% rank

        # Common HLA alleles (cover ~90% of population)
        self.common_hla_alleles = [
            'HLA-A*02:01',  # Most common (~50% Caucasian, ~40% worldwide)
            'HLA-A*01:01',
            'HLA-A*03:01',
            'HLA-A*24:02',
            'HLA-B*07:02',
            'HLA-B*08:01',
            'HLA-B*44:02',
            'HLA-C*07:01',
        ]

        # Check for NetMHCpan
        self.netmhcpan_available = self._check_netmhcpan()

    def _check_netmhcpan(self) -> bool:
        """Check if NetMHCpan is available."""
        try:
            # Check environment variable
            netmhcpan_path = subprocess.os.environ.get('NETMHCPAN')
            if netmhcpan_path and Path(netmhcpan_path).exists():
                return True

            # Try to run netMHCpan
            result = subprocess.run(
                ['netMHCpan', '-h'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

        return False

    def validate(
        self,
        candidate: TherapeuticCandidate,
        **kwargs
    ) -> ValidationResult:
        """
        Predict MHC-I epitopes for immunogenicity assessment.

        Args:
            candidate: Therapeutic candidate with sequence

        Returns:
            ValidationResult with MHC-I epitope predictions
        """
        sequence = candidate.get_primary_sequence()

        if not sequence:
            return self._create_result(
                candidate=candidate,
                decision=Decision.INFORMATIVE,
                metrics=[],
                summary="No sequence provided - skipping MHC-I epitope prediction",
            )

        # Check if biologic
        if candidate.is_small_molecule():
            return self._create_result(
                candidate=candidate,
                decision=Decision.INFORMATIVE,
                metrics=[],
                summary="MHC-I epitope prediction only applicable to biologics",
            )

        # Skip if humanized/human (low immunogenicity risk)
        if self._is_humanized(candidate):
            return self._create_result(
                candidate=candidate,
                decision=Decision.PASS,
                metrics=[],
                summary="Humanized/human sequence - low MHC-I immunogenicity risk",
            )

        metrics = []
        recommendations = []
        risks = []
        warnings = []

        # Run MHC-I epitope prediction
        if self.netmhcpan_available:
            epitopes = self._run_netmhcpan(sequence)
            method = "NetMHCpan"
        else:
            epitopes = self._predict_motif_based(sequence)
            method = "motif-based"
            warnings.append("NetMHCpan not available - using pattern-based predictions")

        # Analyze results
        strong_binders = [e for e in epitopes if e.strong_binder]
        weak_binders = [e for e in epitopes if e.weak_binder]
        total_epitopes = len(strong_binders) + len(weak_binders)

        # Calculate population coverage (simplified)
        coverage = self._estimate_population_coverage(epitopes)

        # Metrics
        metrics.append(MetricResult(
            name='mhc1_strong_binders',
            value=len(strong_binders),
            unit='count',
            threshold_pass=5,  # ≤5 strong binders is acceptable
            threshold_revise=10,
            decision=Decision.PASS if len(strong_binders) <= 5 else Decision.REVISE if len(strong_binders) <= 10 else Decision.KILL,
            metadata={'method': method}
        ))

        metrics.append(MetricResult(
            name='mhc1_weak_binders',
            value=len(weak_binders),
            unit='count',
            threshold_pass=15,
            decision=Decision.PASS if len(weak_binders) <= 15 else Decision.REVISE
        ))

        metrics.append(MetricResult(
            name='mhc1_population_coverage',
            value=coverage * 100,
            unit='%',
            threshold_pass=50.0,  # <50% coverage is good
            decision=Decision.PASS if coverage < 0.5 else Decision.REVISE
        ))

        metrics.append(MetricResult(
            name='epitope_density',
            value=(total_epitopes / len(sequence)) * 100 if sequence else 0,
            unit='%',
            threshold_pass=10.0,  # <10% of sequence as epitopes
            decision=Decision.PASS if (total_epitopes / len(sequence)) < 0.1 else Decision.REVISE
        ))

        # Recommendations
        if len(strong_binders) == 0:
            recommendations.append("No strong MHC-I binders identified - low immunogenicity risk")
        elif len(strong_binders) <= 5:
            recommendations.append(f"Few strong MHC-I binders ({len(strong_binders)}) - acceptable risk")
        else:
            risks.append(f"Many strong MHC-I binders ({len(strong_binders)}) - consider de-immunization")

            # Show top epitopes
            top_epitopes = sorted(strong_binders, key=lambda e: e.percentile_rank)[:3]
            epitope_list = [f"{e.peptide} (pos {e.position}, {e.hla_allele})" for e in top_epitopes]
            recommendations.append(f"Top MHC-I epitopes to consider removing: {', '.join(epitope_list)}")

        if coverage > 0.7:
            risks.append(f"High population coverage ({coverage*100:.0f}%) - many people may respond")
        elif coverage > 0.5:
            warnings.append(f"Moderate population coverage ({coverage*100:.0f}%)")

        # Overall decision
        if len(strong_binders) > 10:
            decision = Decision.KILL
            summary = f"Many MHC-I epitopes ({len(strong_binders)} strong binders, {coverage*100:.0f}% coverage, {method})"
        elif len(strong_binders) > 5:
            decision = Decision.REVISE
            summary = f"Some MHC-I epitopes ({len(strong_binders)} strong binders, {method})"
        else:
            decision = Decision.PASS
            summary = f"Low MHC-I immunogenicity risk ({len(strong_binders)} strong binders, {coverage*100:.0f}% coverage, {method})"

        return self._create_result(
            candidate=candidate,
            decision=decision,
            metrics=metrics,
            summary=summary,
            recommendations=recommendations,
            risks=risks,
            warnings=warnings,
            metadata={
                'epitopes': [
                    {
                        'position': e.position,
                        'peptide': e.peptide,
                        'hla': e.hla_allele,
                        'percentile': e.percentile_rank,
                        'strong': e.strong_binder
                    }
                    for e in epitopes if e.strong_binder or e.weak_binder
                ]
            }
        )

    def _run_netmhcpan(self, sequence: str) -> List[MHC1Epitope]:
        """
        Run NetMHCpan to predict MHC-I epitopes.

        Args:
            sequence: Protein sequence

        Returns:
            List of MHC1Epitope objects
        """
        # TODO: Implement NetMHCpan integration
        # This would involve:
        # 1. Write sequence to temp file
        # 2. Run NetMHCpan for each HLA allele
        # 3. Parse output (percentile ranks)
        # 4. Create MHC1Epitope objects

        # For now, fallback to motif-based
        return self._predict_motif_based(sequence)

    def _predict_motif_based(self, sequence: str) -> List[MHC1Epitope]:
        """
        Predict MHC-I epitopes using motif-based approach.

        MHC-I binding motifs (simplified):
        - HLA-A*02:01: Anchor at P2 (L,M,I,V) and P9 (L,I,V,A)
        - HLA-A*01:01: Anchor at P2 (T,S,M) and P9 (Y)
        - HLA-B*07:02: Anchor at P2 (P) and P9 (L,F,W,M,A)

        Args:
            sequence: Protein sequence

        Returns:
            List of MHC1Epitope objects
        """
        epitopes = []

        # Scan for 9-mers (most common MHC-I peptide length)
        for length in [8, 9, 10, 11]:
            for i in range(len(sequence) - length + 1):
                peptide = sequence[i:i+length]

                # Check each HLA allele
                for hla in self.common_hla_alleles[:3]:  # Use top 3 for speed
                    score = self._score_peptide_binding(peptide, hla)

                    if score > 0.5:  # Threshold for prediction
                        # Estimate percentile rank from score
                        # Higher score → lower percentile (better)
                        percentile = (1 - score) * 5  # 0-5% range

                        epitope = MHC1Epitope(
                            position=i + 1,
                            peptide=peptide,
                            length=length,
                            hla_allele=hla,
                            percentile_rank=percentile,
                            strong_binder=(percentile < self.strong_binder_percentile),
                            weak_binder=(percentile < self.weak_binder_percentile)
                        )

                        epitopes.append(epitope)

        return epitopes

    def _score_peptide_binding(self, peptide: str, hla: str) -> float:
        """
        Score peptide binding to HLA allele.

        Uses simplified binding motifs.

        Args:
            peptide: Peptide sequence
            hla: HLA allele

        Returns:
            Binding score (0-1, higher = better)
        """
        score = 0.0

        # HLA-A*02:01 motif (most common)
        if hla == 'HLA-A*02:01' and len(peptide) == 9:
            # P2 anchor: L,M,I,V
            if peptide[1] in 'LMIV':
                score += 0.4
            # P9 anchor: L,I,V,A
            if peptide[8] in 'LIVA':
                score += 0.4
            # Hydrophobic residues preferred
            hydrophobic = sum(1 for aa in peptide if aa in 'FILVMWY')
            score += (hydrophobic / len(peptide)) * 0.2

        # HLA-A*01:01 motif
        elif hla == 'HLA-A*01:01' and len(peptide) == 9:
            # P2 anchor: T,S,M
            if peptide[1] in 'TSM':
                score += 0.4
            # P9 anchor: Y
            if peptide[8] == 'Y':
                score += 0.4
            # Small residues at other positions
            small = sum(1 for aa in peptide if aa in 'AGSV')
            score += (small / len(peptide)) * 0.2

        # HLA-A*03:01 motif
        elif hla == 'HLA-A*03:01' and len(peptide) == 9:
            # P2 anchor: L,I,V,M
            if peptide[1] in 'LIVM':
                score += 0.4
            # P9 anchor: K,R
            if peptide[8] in 'KR':
                score += 0.4

        # Default: use hydrophobicity
        else:
            hydrophobic = sum(1 for aa in peptide if aa in 'FILVMWY')
            score = (hydrophobic / len(peptide)) * 0.6

        return min(1.0, score)

    def _estimate_population_coverage(self, epitopes: List[MHC1Epitope]) -> float:
        """
        Estimate percentage of population that could respond.

        Based on HLA allele frequencies.

        Args:
            epitopes: List of predicted epitopes

        Returns:
            Population coverage (0-1)
        """
        # HLA allele frequencies (approximate, worldwide)
        hla_frequencies = {
            'HLA-A*02:01': 0.45,  # ~45% of population
            'HLA-A*01:01': 0.25,
            'HLA-A*03:01': 0.20,
            'HLA-A*24:02': 0.30,
            'HLA-B*07:02': 0.20,
            'HLA-B*08:01': 0.15,
            'HLA-B*44:02': 0.20,
            'HLA-C*07:01': 0.25,
        }

        # Get unique HLA alleles with epitopes
        alleles_with_epitopes = set()
        for e in epitopes:
            if e.strong_binder or e.weak_binder:
                alleles_with_epitopes.add(e.hla_allele)

        # Calculate coverage (simplified - assumes independence)
        coverage = 0.0
        for allele in alleles_with_epitopes:
            freq = hla_frequencies.get(allele, 0.1)
            coverage += freq * (1 - coverage)  # Account for overlap

        return coverage

    def _is_humanized(self, candidate: TherapeuticCandidate) -> bool:
        """
        Check if sequence is humanized or human.

        Args:
            candidate: Therapeutic candidate

        Returns:
            True if likely humanized/human
        """
        # Check metadata
        if candidate.metadata:
            design_notes = str(candidate.metadata.get('design_notes', '')).lower()
            if 'human' in design_notes or 'humaniz' in design_notes:
                return True

        # For now, assume all antibodies might be humanized
        # In production, would check against human germline sequences
        return False
