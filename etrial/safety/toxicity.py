"""
Pre-Filter Safety Module - Fast Toxicity Screening.

Uses RDKit for fast structural alerts and basic toxicity filters:
    - PAINS (Pan Assay Interference Compounds)
    - Lipinski's Rule of Five
    - Brenk structural alerts
    - Reactive functional groups

Speed: ~1ms per molecule
Accuracy: ~70% (filters obvious toxic structures)

For clinical-grade predictions, use toxicity_clinical.py
"""

from typing import List, Dict, Any, Optional
from loguru import logger

from etrial.core.base import (
    TherapeuticCandidate,
    ValidationModule,
    ValidationResult,
    Decision,
    MetricResult,
)


class ToxicityPrefilterModule(ValidationModule):
    """
    Pre-filter toxicity screening using structural alerts.

    Fast screening for:
        - PAINS (promiscuous compounds)
        - Lipinski violations (drug-likeness)
        - Structural alerts (Brenk filters)
        - Reactive groups
    """

    def _setup(self) -> None:
        """Check for RDKit availability."""
        logger.info("Setting up ToxicityPrefilterModule")

        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski
            self.rdkit_available = True
            logger.info("RDKit available for structural analysis")
        except ImportError:
            self.rdkit_available = False
            logger.warning(
                "RDKit not available. Install with: pip install rdkit-pypi\n"
                "Using basic sequence-based heuristics only."
            )

    def validate(
        self,
        candidate: TherapeuticCandidate
    ) -> ValidationResult:
        """
        Run pre-filter toxicity screening.

        Args:
            candidate: Therapeutic candidate

        Returns:
            ValidationResult with toxicity flags
        """
        logger.info(f"Running pre-filter toxicity screening for {candidate.name}")

        metrics: List[MetricResult] = []
        risks: List[str] = []
        warnings: List[str] = []
        recommendations: List[str] = []

        # ================================================================
        # Small Molecule Toxicity Screening
        # ================================================================

        if candidate.is_small_molecule():
            if self.rdkit_available:
                tox_result = self._screen_small_molecule(candidate)
            else:
                tox_result = self._basic_small_molecule_check(candidate)
                warnings.append("RDKit not available - using basic checks only")

            # PAINS flags
            if tox_result.get('pains_alerts', 0) > 0:
                metrics.append(MetricResult(
                    name="pains_alerts",
                    value=tox_result['pains_alerts'],
                    unit="count",
                    threshold_pass=0,
                    threshold_revise=1,
                    decision=Decision.KILL if tox_result['pains_alerts'] > 1 else Decision.REVISE,
                ))
                risks.append(f"PAINS alerts: {tox_result['pains_alerts']} promiscuous structures detected")

            # Lipinski violations
            lipinski_violations = tox_result.get('lipinski_violations', 0)
            metrics.append(MetricResult(
                name="lipinski_violations",
                value=lipinski_violations,
                unit="count",
                threshold_pass=1,
                threshold_revise=2,
                decision=Decision.PASS if lipinski_violations <= 1 else (
                    Decision.REVISE if lipinski_violations <= 2 else Decision.KILL
                ),
            ))

            if lipinski_violations > 1:
                warnings.append(f"Lipinski violations: {lipinski_violations} (drug-likeness issues)")

            # Structural alerts (Brenk)
            if tox_result.get('structural_alerts', 0) > 0:
                metrics.append(MetricResult(
                    name="structural_alerts",
                    value=tox_result['structural_alerts'],
                    unit="count",
                    threshold_pass=0,
                    threshold_revise=2,
                    decision=Decision.REVISE if tox_result['structural_alerts'] < 3 else Decision.KILL,
                ))
                risks.append(f"Structural alerts: {tox_result['structural_alerts']} toxic substructures")

            # Reactive groups
            if tox_result.get('reactive_groups', 0) > 0:
                metrics.append(MetricResult(
                    name="reactive_groups",
                    value=tox_result['reactive_groups'],
                    unit="count",
                    threshold_pass=0,
                    threshold_revise=1,
                    decision=Decision.REVISE if tox_result['reactive_groups'] < 2 else Decision.KILL,
                ))
                warnings.append(f"Reactive groups: {tox_result['reactive_groups']} (instability risk)")

            # Add recommendations
            if tox_result['pains_alerts'] > 0:
                recommendations.append(
                    "PAINS detected - likely assay interference. Consider scaffold redesign."
                )

            if lipinski_violations > 2:
                recommendations.append(
                    f"Multiple Lipinski violations ({lipinski_violations}). "
                    "Consider: (1) reducing molecular weight, (2) improving lipophilicity"
                )

        # ================================================================
        # Biologics Toxicity Screening
        # ================================================================

        else:
            # For biologics, check sequence-based toxicity flags
            sequence = candidate.get_primary_sequence()

            if sequence:
                biologic_result = self._screen_biologic(sequence)

                # Poly-basic regions (cell penetration, toxicity)
                if biologic_result.get('poly_basic_regions', 0) > 0:
                    metrics.append(MetricResult(
                        name="poly_basic_regions",
                        value=biologic_result['poly_basic_regions'],
                        unit="count",
                        threshold_pass=1,
                        threshold_revise=3,
                        decision=Decision.PASS if biologic_result['poly_basic_regions'] <= 1 else (
                            Decision.REVISE if biologic_result['poly_basic_regions'] <= 3 else Decision.KILL
                        ),
                    ))

                    if biologic_result['poly_basic_regions'] > 1:
                        warnings.append(
                            f"Poly-basic regions: {biologic_result['poly_basic_regions']} "
                            "(potential cytotoxicity)"
                        )

                # Hydrophobic patches (membrane disruption)
                if biologic_result.get('large_hydrophobic_patches', 0) > 0:
                    metrics.append(MetricResult(
                        name="hydrophobic_patches",
                        value=biologic_result['large_hydrophobic_patches'],
                        unit="count",
                        threshold_pass=1,
                        threshold_revise=3,
                        decision=Decision.PASS if biologic_result['large_hydrophobic_patches'] <= 1 else Decision.REVISE,
                    ))

                # Add recommendations for biologics
                if biologic_result['poly_basic_regions'] > 2:
                    recommendations.append(
                        "High poly-basic content - risk of cell membrane disruption. "
                        "Consider charge balancing."
                    )
            else:
                return self._create_result(
                    candidate=candidate,
                    decision=Decision.INFORMATIVE,
                    metrics=[],
                    summary="No sequence/structure available for toxicity screening",
                    warnings=["No sequence or structure available"],
                )

        # ================================================================
        # Overall Decision
        # ================================================================

        if not metrics:
            # No toxicity flags detected
            return self._create_result(
                candidate=candidate,
                decision=Decision.PASS,
                metrics=[],
                summary="No toxicity flags detected in pre-filter screening",
            )

        # Aggregate decisions
        decisions = [m.decision for m in metrics]

        if Decision.KILL in decisions:
            overall_decision = Decision.KILL
        elif Decision.REVISE in decisions:
            overall_decision = Decision.REVISE
        else:
            overall_decision = Decision.PASS

        # Create summary
        n_flags = len([d for d in decisions if d in [Decision.KILL, Decision.REVISE]])
        summary = f"Pre-filter toxicity screening: {n_flags} flags detected"

        return self._create_result(
            candidate=candidate,
            decision=overall_decision,
            metrics=metrics,
            summary=summary,
            risks=risks,
            warnings=warnings,
            recommendations=recommendations,
        )

    def _screen_small_molecule(self, candidate: TherapeuticCandidate) -> Dict[str, Any]:
        """
        Screen small molecule for toxicity using RDKit.

        Returns:
            dict with pains_alerts, lipinski_violations, structural_alerts, reactive_groups
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, FilterCatalog

            # Parse SMILES
            smiles = candidate.smiles or candidate.sequence
            if not smiles:
                return {'pains_alerts': 0, 'lipinski_violations': 0, 'structural_alerts': 0, 'reactive_groups': 0}

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"Could not parse SMILES: {smiles}")
                return {'pains_alerts': 0, 'lipinski_violations': 0, 'structural_alerts': 0, 'reactive_groups': 0}

            result = {
                'pains_alerts': 0,
                'lipinski_violations': 0,
                'structural_alerts': 0,
                'reactive_groups': 0,
            }

            # ============================================================
            # 1. PAINS Filters
            # ============================================================
            try:
                params = FilterCatalog.FilterCatalogParams()
                params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
                catalog = FilterCatalog.FilterCatalog(params)

                matches = catalog.GetMatches(mol)
                result['pains_alerts'] = len(matches)

            except Exception as e:
                logger.warning(f"PAINS filtering failed: {e}")

            # ============================================================
            # 2. Lipinski's Rule of Five
            # ============================================================
            violations = 0

            mw = Descriptors.MolWt(mol)
            if mw > 500:
                violations += 1

            logp = Descriptors.MolLogP(mol)
            if logp > 5:
                violations += 1

            hbd = Lipinski.NumHDonors(mol)
            if hbd > 5:
                violations += 1

            hba = Lipinski.NumHAcceptors(mol)
            if hba > 10:
                violations += 1

            result['lipinski_violations'] = violations

            # ============================================================
            # 3. Structural Alerts (Brenk)
            # ============================================================
            # Simplified Brenk filters - common toxic substructures
            toxic_smarts = [
                '[N;!R]=[N;!R]',  # Azo groups
                '[S,s][S,s]',  # Disulfides
                '[N;!R][N;!R][N;!R]',  # Triazenes
                'C(=O)O[N;!R]',  # Hydroxamic acids
                '[Cl,Br,I][C;!R]=[C;!R]',  # Halogenated alkenes
                'N=C=S',  # Isothiocyanates
                'C(=O)Cl',  # Acid chlorides
                '[N;R0][N;R0]C(=O)',  # Hydrazides
            ]

            alerts = 0
            for smarts in toxic_smarts:
                patt = Chem.MolFromSmarts(smarts)
                if patt and mol.HasSubstructMatch(patt):
                    alerts += 1

            result['structural_alerts'] = alerts

            # ============================================================
            # 4. Reactive Groups
            # ============================================================
            reactive_smarts = [
                '[C;R0](=O)[C;R0](=O)',  # 1,2-Dicarbonyls
                '[N;R0]=[C;R0]=[N;R0]',  # Carbodiimides
                'C#N',  # Nitriles (many instances)
                '[N+](=O)[O-]',  # Nitro groups
                'C(=O)C(=O)',  # Alpha-dicarbonyl
            ]

            reactive = 0
            for smarts in reactive_smarts:
                patt = Chem.MolFromSmarts(smarts)
                if patt:
                    matches = mol.GetSubstructMatches(patt)
                    if len(matches) > 2:  # More than 2 instances = reactive
                        reactive += 1

            result['reactive_groups'] = reactive

            return result

        except Exception as e:
            logger.error(f"Small molecule screening failed: {e}")
            return {'pains_alerts': 0, 'lipinski_violations': 0, 'structural_alerts': 0, 'reactive_groups': 0}

    def _basic_small_molecule_check(self, candidate: TherapeuticCandidate) -> Dict[str, Any]:
        """Fallback: Basic checks without RDKit."""
        # Can't do much without RDKit
        return {
            'pains_alerts': 0,
            'lipinski_violations': 0,
            'structural_alerts': 0,
            'reactive_groups': 0,
        }

    def _screen_biologic(self, sequence: str) -> Dict[str, Any]:
        """
        Screen biologic (protein/peptide) for toxicity.

        Returns:
            dict with poly_basic_regions, large_hydrophobic_patches
        """
        import re

        result = {
            'poly_basic_regions': 0,
            'large_hydrophobic_patches': 0,
        }

        # Poly-basic regions (cytotoxicity risk)
        poly_basic_patterns = [
            r'K{5,}',  # Poly-lysine
            r'R{5,}',  # Poly-arginine
            r'[KR]{7,}',  # Mixed poly-basic
        ]

        for pattern in poly_basic_patterns:
            matches = re.findall(pattern, sequence)
            result['poly_basic_regions'] += len(matches)

        # Large hydrophobic patches (membrane disruption)
        hydrophobic = set('VILFMWY')
        window_size = 10
        hydrophobic_threshold = 8  # 80% hydrophobic

        patches = 0
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i+window_size]
            hydrophobic_count = sum(1 for aa in window if aa in hydrophobic)
            if hydrophobic_count >= hydrophobic_threshold:
                patches += 1

        result['large_hydrophobic_patches'] = patches

        return result


def check_rdkit_available() -> bool:
    """Check if RDKit is installed."""
    try:
        import rdkit
        return True
    except ImportError:
        return False
