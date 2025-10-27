"""
Developability: Solubility and Aggregation Assessment

This module validates:
- Solubility predictions (CamSol, OptSol)
- Aggregation propensity (AggreScan3D)
- PTM liabilities (deamidation, oxidation, isomerization)
- Expression predictions

TODO: Integrate CamSol, AggreScan3D, expression predictors
"""

import time
from typing import Any, Dict, List, Optional
from loguru import logger

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    MetricResult,
    Decision,
)


class SolubilityValidationModule(ValidationModule):
    """
    Validates solubility, aggregation, and manufacturability.

    This is a STUB implementation showing the pattern.
    In production, this would run:
    - CamSol for solubility prediction
    - AggreScan3D for aggregation hotspots
    - PTM liability scanners
    - Expression predictors
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self.name = "developability"
        self.version = "0.1.0"

    def _setup(self) -> None:
        """Initialize developability tools."""
        logger.info(f"Setting up {self.name} module")

        # Configuration
        self.check_liabilities = self.config.get('liabilities', {})
        self.target_concentration = self.config.get('expression', {}).get(
            'target_concentration_mg_ml', 100
        )

    def validate(
        self,
        candidate: TherapeuticCandidate,
        **kwargs
    ) -> ValidationResult:
        """
        Run developability validation.

        Args:
            candidate: Therapeutic candidate to validate

        Returns:
            ValidationResult with developability metrics
        """
        start_time = time.time()
        logger.info(f"Running {self.name} validation for {candidate.name}")

        # Skip if small molecule
        if candidate.is_small_molecule():
            logger.info("Small molecule detected - using alternative ADMET pipeline")
            # TODO: Route to small molecule ADMET module
            return self._small_molecule_stub(candidate)

        metrics: List[MetricResult] = []

        # ==================================================================
        # STEP 1: Solubility Prediction
        # ==================================================================
        logger.info("Predicting solubility...")

        # TODO: Run CamSol
        camsol_score = 1.2  # >1.0 = soluble (dummy value)

        if camsol_score >= 1.0:
            sol_decision = Decision.PASS
        elif camsol_score >= 0.5:
            sol_decision = Decision.REVISE
        else:
            sol_decision = Decision.KILL

        metrics.append(MetricResult(
            name="camsol_score",
            value=camsol_score,
            unit="score",
            threshold_pass=1.0,
            threshold_revise=0.5,
            decision=sol_decision,
        ))

        # ==================================================================
        # STEP 2: Aggregation Propensity
        # ==================================================================
        logger.info("Analyzing aggregation propensity...")

        # Get primary sequence
        sequence = candidate.get_primary_sequence() or ""

        # Real hydrophobic patch detection (sequence-based)
        hydrophobic_analysis = self._detect_hydrophobic_patches(sequence)

        # Calculate aggregation decision based on patch analysis
        # Note: Some hydrophobic clustering is NORMAL in functional proteins
        # Only flag severe cases
        max_patch_size = hydrophobic_analysis['max_patch_size']
        total_hydrophobic = hydrophobic_analysis['hydrophobic_residues']
        seq_length = len(sequence) if sequence else 1

        if max_patch_size >= 30:
            agg_decision = Decision.KILL  # Extreme continuous hydrophobicity
        elif max_patch_size >= 12 and total_hydrophobic / seq_length > 0.5:
            agg_decision = Decision.KILL  # Large patches + high overall hydrophobicity
        elif max_patch_size >= 8:
            agg_decision = Decision.REVISE  # Moderate-large clustered hydrophobicity
        elif max_patch_size >= 5 and hydrophobic_analysis['n_hotspots'] >= 3:
            agg_decision = Decision.REVISE  # Multiple moderate patches
        else:
            agg_decision = Decision.PASS  # Normal hydrophobicity distribution

        # Estimated AggreScan-like score (negative = good, positive = bad)
        # Based on hydrophobic content and clustering
        hydrophobic_fraction = hydrophobic_analysis['hydrophobic_residues'] / len(sequence) if sequence else 0
        estimated_aggrescan = (hydrophobic_fraction - 0.3) * 200  # Rough approximation

        metrics.append(MetricResult(
            name="hydrophobic_patches",
            value=hydrophobic_analysis['n_hotspots'],
            unit="count",
            threshold_pass=0,
            decision=agg_decision,
            metadata={
                "max_patch_size": max_patch_size,
                "positions": hydrophobic_analysis['patch_positions'],
                "severity": hydrophobic_analysis['severity']
            },
        ))

        metrics.append(MetricResult(
            name="estimated_aggregation_score",
            value=round(estimated_aggrescan, 1),
            unit="score",
            threshold_pass=10,
            threshold_revise=50,
            decision=Decision.PASS if estimated_aggrescan < 10 else (
                Decision.REVISE if estimated_aggrescan < 50 else Decision.KILL
            ),
        ))

        # ==================================================================
        # STEP 3: PTM Liability Scanning
        # ==================================================================
        logger.info("Scanning for PTM liabilities...")

        # Get primary sequence (works for all modalities)
        sequence = candidate.get_primary_sequence() or ""

        if not sequence:
            logger.warning("No sequence available for PTM scanning")
            sequence = ""

        # Real PTM scanning using regex patterns
        ptm_results = self._scan_ptm_liabilities(sequence)

        metrics.extend([
            MetricResult(
                name="n_glycosylation_sites",
                value=ptm_results['n_glycosylation'],
                unit="count",
                threshold_pass=2,
                decision=Decision.PASS if ptm_results['n_glycosylation'] <= 2 else (
                    Decision.KILL if ptm_results['n_glycosylation'] >= 5 else Decision.REVISE
                ),
                metadata={"positions": ptm_results['glyc_positions']},
            ),
            MetricResult(
                name="n_deamidation_sites",
                value=ptm_results['n_deamidation'],
                unit="count",
                threshold_pass=1,
                decision=Decision.PASS if ptm_results['n_deamidation'] <= 1 else (
                    Decision.KILL if ptm_results['n_deamidation'] >= 7 else Decision.REVISE
                ),
                metadata={"positions": ptm_results['deamid_positions']},
            ),
            MetricResult(
                name="n_oxidation_sites",
                value=ptm_results['n_methionine'],
                unit="count",
                threshold_pass=4,
                decision=Decision.PASS if ptm_results['n_methionine'] <= 4 else Decision.REVISE,
                metadata={"positions": ptm_results['met_positions']},
            ),
            MetricResult(
                name="n_free_cysteines",
                value=ptm_results['n_free_cys'],
                unit="count",
                threshold_pass=0,
                decision=Decision.KILL if ptm_results['n_free_cys'] > 0 else Decision.PASS,
                metadata={"total_cys": ptm_results['total_cys'], "positions": ptm_results['cys_positions']},
            ),
            MetricResult(
                name="n_isomerization_sites",
                value=ptm_results['n_isomerization'],
                unit="count",
                threshold_pass=2,
                decision=Decision.PASS if ptm_results['n_isomerization'] <= 2 else Decision.REVISE,
                metadata={"positions": ptm_results['isom_positions']},
            ),
        ])

        # ==================================================================
        # STEP 4: Expression Prediction
        # ==================================================================
        logger.info("Predicting expression and formulation properties...")

        # TODO: Run expression predictors
        predicted_expression = "high"  # (dummy value)
        predicted_pI = 7.5  # (dummy value)
        predicted_viscosity = 12  # cP at 100 mg/mL (dummy value)

        metrics.extend([
            MetricResult(
                name="predicted_expression_level",
                value=predicted_expression,
                decision=Decision.PASS if predicted_expression == "high" else Decision.REVISE,
            ),
            MetricResult(
                name="predicted_pI",
                value=predicted_pI,
                unit="pH",
                decision=Decision.PASS if 6.0 <= predicted_pI <= 8.5 else Decision.REVISE,
                metadata={"optimal_range": "6.0-8.5 pH"},
            ),
            MetricResult(
                name="predicted_viscosity",
                value=predicted_viscosity,
                unit="cP",
                threshold_pass=20,
                decision=Decision.PASS if predicted_viscosity <= 20 else Decision.REVISE,
            ),
        ])

        # ==================================================================
        # STEP 5: Overall Decision
        # ==================================================================

        decisions = [m.decision for m in metrics if m.decision]

        if Decision.KILL in decisions:
            overall_decision = Decision.KILL
            summary = "Critical developability issues detected (low solubility or free cysteines)."

        elif Decision.REVISE in decisions:
            overall_decision = Decision.REVISE
            n_revise = decisions.count(Decision.REVISE)
            summary = f"Developability issues identified ({n_revise} metrics need improvement)."

        else:
            overall_decision = Decision.PASS
            summary = "Good developability profile with predicted high solubility and expression."

        # Recommendations
        recommendations = []
        risks = []

        if camsol_score < 1.0:
            recommendations.append(
                "Consider surface engineering to improve solubility (target CamSol >1.0)"
            )

        if hydrophobic_analysis['n_hotspots'] > 0:
            severity = hydrophobic_analysis['severity']
            max_patch = hydrophobic_analysis['max_patch_size']
            recommendations.append(
                f"Reduce {severity} hydrophobic clustering: {hydrophobic_analysis['n_hotspots']} patch(es) detected (max size: {max_patch} residues)"
            )
            if hydrophobic_analysis['patch_positions']:
                top_patch = hydrophobic_analysis['patch_positions'][0]
                recommendations.append(
                    f"  → Largest patch at position {top_patch[0]}-{top_patch[1]}: {top_patch[2]}"
                )

        if ptm_results['n_glycosylation'] > 2:
            recommendations.append(
                f"Consider removing {ptm_results['n_glycosylation'] - 2} N-glycosylation sites at positions: {ptm_results['glyc_positions'][:3]}"
            )

        if ptm_results['n_deamidation'] > 3:
            recommendations.append(
                f"Reduce deamidation risk by mutating {ptm_results['n_deamidation']} NG/NN/NS sites at positions: {ptm_results['deamid_positions'][:3]}"
            )

        if ptm_results['n_free_cys'] > 0:
            risks.append(
                f"Unpaired cysteine detected ({ptm_results['total_cys']} total Cys) - CRITICAL: risk of aggregation and disulfide scrambling"
            )

        if ptm_results['n_methionine'] > 4:
            recommendations.append(
                f"Consider reducing methionine content ({ptm_results['n_methionine']} Met residues) to minimize oxidation risk"
            )

        if ptm_results['n_isomerization'] > 2:
            recommendations.append(
                f"Reduce isomerization risk by mutating {ptm_results['n_isomerization']} DG motifs"
            )

        if predicted_viscosity > 20:
            recommendations.append(
                "High viscosity predicted - consider charge optimization for formulation"
            )

        # Calculate runtime
        runtime = time.time() - start_time

        # Create result
        result = self._create_result(
            candidate=candidate,
            decision=overall_decision,
            metrics=metrics,
            summary=summary,
            recommendations=recommendations,
            risks=risks,
            runtime_seconds=runtime,
        )

        logger.info(f"{self.name} validation completed: {overall_decision.value} ({runtime:.2f}s)")

        return result

    def _scan_ptm_liabilities(self, sequence: str) -> Dict[str, Any]:
        """
        Scan sequence for post-translational modification liabilities.

        Args:
            sequence: Amino acid sequence

        Returns:
            Dictionary with PTM counts and positions
        """
        import re

        results = {
            'n_glycosylation': 0,
            'glyc_positions': [],
            'n_deamidation': 0,
            'deamid_positions': [],
            'n_methionine': 0,
            'met_positions': [],
            'n_free_cys': 0,
            'total_cys': 0,
            'cys_positions': [],
            'n_isomerization': 0,
            'isom_positions': [],
        }

        if not sequence:
            return results

        # N-glycosylation: NXS or NXT motifs (where X != P)
        # Pattern: N[^P][ST]
        glyc_pattern = r'N[^P][ST]'
        for match in re.finditer(glyc_pattern, sequence):
            results['n_glycosylation'] += 1
            results['glyc_positions'].append((match.start(), match.group()))

        # Deamidation: NG, NN, NS motifs
        deamid_patterns = ['NG', 'NN', 'NS']
        for pattern in deamid_patterns:
            for match in re.finditer(pattern, sequence):
                results['n_deamidation'] += 1
                results['deamid_positions'].append((match.start(), match.group()))

        # Methionine oxidation: Count all M
        for match in re.finditer('M', sequence):
            results['n_methionine'] += 1
            results['met_positions'].append(match.start())

        # Cysteine analysis
        results['total_cys'] = sequence.count('C')
        for match in re.finditer('C', sequence):
            results['cys_positions'].append(match.start())

        # Check for unpaired cysteines
        # Only flag as unpaired if it's clearly problematic:
        # - Single Cys in short peptide (<35 residues) → likely unpaired
        # - Multiple odd Cys (3, 5, 7...) in medium peptide → suspicious
        # - For longer sequences or single Cys in antibody-length sequences, ignore
        seq_length = len(sequence)
        if results['total_cys'] == 1 and seq_length < 35:
            # Single Cys in short peptide - almost certainly unpaired → KILL
            results['n_free_cys'] = 1
        elif results['total_cys'] >= 3 and results['total_cys'] % 2 == 1 and seq_length < 100:
            # Multiple odd Cys (3, 5, 7...) - likely unpaired → flag for review
            results['n_free_cys'] = results['total_cys'] % 2  # Will be 1 for odd
        else:
            # Even number, or single Cys in longer sequence (may be inter-chain)
            results['n_free_cys'] = 0

        # Isomerization: DG motifs (aspartic acid isomerization)
        for match in re.finditer('DG', sequence):
            results['n_isomerization'] += 1
            results['isom_positions'].append((match.start(), match.group()))

        return results

    def _detect_hydrophobic_patches(self, sequence: str) -> Dict[str, Any]:
        """
        Detect hydrophobic patches in sequence using sliding window.

        Hydrophobic residues: V, I, L, F, M, W (highly hydrophobic)

        Args:
            sequence: Amino acid sequence

        Returns:
            Dictionary with patch analysis
        """
        results = {
            'n_hotspots': 0,
            'max_patch_size': 0,
            'patch_positions': [],
            'severity': 'none',
            'hydrophobic_residues': 0,
        }

        if not sequence or len(sequence) < 3:
            return results

        # Define hydrophobic residues
        hydrophobic = set('VILFMW')

        # Count total hydrophobic residues
        results['hydrophobic_residues'] = sum(1 for aa in sequence if aa in hydrophobic)

        # Detect long stretches of consecutive/nearly-consecutive hydrophobic residues
        # Only flag truly problematic aggregation-prone regions
        patches = []
        i = 0

        while i < len(sequence):
            if sequence[i] in hydrophobic:
                # Start of potential patch
                patch_start = i
                hydrophobic_count = 0
                j = i

                # Extend while we see hydrophobic residues (allowing 1-2 non-hydrophobic gaps)
                consecutive_non_hydrophobic = 0

                while j < len(sequence):
                    if sequence[j] in hydrophobic:
                        hydrophobic_count += 1
                        consecutive_non_hydrophobic = 0
                        j += 1
                    elif consecutive_non_hydrophobic < 2:  # Allow small gaps
                        consecutive_non_hydrophobic += 1
                        j += 1
                    else:
                        break  # End of patch

                patch_length = j - patch_start
                patch_end = j - 1

                # Only record if it's a significant aggregation-prone region
                # Require >=5 hydrophobic out of <=8 residues (high density)
                if hydrophobic_count >= 5 and patch_length <= hydrophobic_count + 2:
                    patches.append({
                        'start': patch_start,
                        'end': patch_end,
                        'size': hydrophobic_count,
                        'sequence': sequence[patch_start:patch_end+1]
                    })

                i = j
            else:
                i += 1

        # Analyze patches
        if patches:
            results['n_hotspots'] = len(patches)
            results['max_patch_size'] = max(p['size'] for p in patches)
            results['patch_positions'] = [(p['start'], p['end'], p['sequence']) for p in patches]

            # Determine severity
            max_size = results['max_patch_size']
            if max_size >= 20:
                results['severity'] = 'severe'
            elif max_size >= 6:
                results['severity'] = 'moderate'
            elif max_size >= 3:
                results['severity'] = 'mild'
            else:
                results['severity'] = 'low'
        else:
            results['severity'] = 'none'

        return results

    def _small_molecule_stub(self, candidate: TherapeuticCandidate) -> ValidationResult:
        """Placeholder for small molecule ADMET."""
        return self._create_result(
            candidate=candidate,
            decision=Decision.INFORMATIVE,
            metrics=[],
            summary="Small molecule ADMET module not yet implemented.",
            recommendations=["Implement small molecule ADMET pipeline (Lipinski, PAINS, etc.)"],
        )
