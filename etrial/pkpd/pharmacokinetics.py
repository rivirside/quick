"""
PKPD Pre-Filter Module.

Fast rule-based predictions for initial screening:
    - Oral bioavailability (Veber rules, Lipinski)
    - BBB permeability (CNS drugs)
    - Plasma protein binding
    - Basic clearance/half-life estimates

Speed: ~1ms per molecule
Accuracy: ~60-70% (filters obvious poor PK candidates)

For clinical-grade predictions, use pharmacokinetics_clinical.py
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


class PKPDPrefilterModule(ValidationModule):
    """
    Pre-filter validation of pharmacokinetics using rule-based predictions.

    Uses fast heuristics:
    - Lipinski's Rule of Five
    - Veber rules (rotatable bonds, PSA)
    - CNS MPO score (for brain penetration)
    - Basic bioavailability estimates
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self.name = "pkpd_prefilter"
        self.version = "1.0.0"

    def _setup(self) -> None:
        """Initialize RDKit tools."""
        logger.info(f"Setting up {self.name} module")

        # Check for RDKit
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski
            self.rdkit_available = True
            logger.info("RDKit available for PKPD analysis")
        except ImportError:
            self.rdkit_available = False
            logger.warning("RDKit not available - using basic sequence analysis")

    def validate(
        self,
        candidate: TherapeuticCandidate,
        **kwargs
    ) -> ValidationResult:
        """
        Run fast pre-filter PKPD validation.

        Args:
            candidate: Therapeutic candidate to validate

        Returns:
            ValidationResult with predicted PK metrics
        """
        start_time = time.time()
        logger.info(f"Running {self.name} validation for {candidate.name}")

        metrics: List[MetricResult] = []
        warnings: List[str] = []

        # ==================================================================
        # Small Molecules: Use Lipinski, Veber, etc.
        # ==================================================================

        if candidate.is_small_molecule():
            if not self.rdkit_available:
                return self._create_result(
                    candidate=candidate,
                    decision=Decision.INFORMATIVE,
                    metrics=[],
                    summary="RDKit not available for PKPD pre-filter",
                    warnings=["Install RDKit: pip install rdkit"],
                    runtime_seconds=time.time() - start_time,
                )

            # Analyze SMILES
            smiles = candidate.smiles
            if not smiles:
                return self._create_result(
                    candidate=candidate,
                    decision=Decision.INFORMATIVE,
                    metrics=[],
                    summary="No SMILES available for PKPD analysis",
                    warnings=["Small molecule requires SMILES"],
                    runtime_seconds=time.time() - start_time,
                )

            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, Crippen

            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return self._create_result(
                    candidate=candidate,
                    decision=Decision.KILL,
                    metrics=[],
                    summary="Invalid SMILES structure",
                    warnings=["Could not parse SMILES"],
                    runtime_seconds=time.time() - start_time,
                )

            # ===============================================================
            # STEP 1: Lipinski's Rule of Five (Drug-likeness)
            # ===============================================================

            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)

            lipinski_violations = 0
            if mw > 500:
                lipinski_violations += 1
            if logp > 5:
                lipinski_violations += 1
            if hbd > 5:
                lipinski_violations += 1
            if hba > 10:
                lipinski_violations += 1

            metrics.append(MetricResult(
                name="lipinski_violations",
                value=lipinski_violations,
                unit="count",
                threshold_pass=0,
                threshold_revise=2,
                decision=Decision.PASS if lipinski_violations == 0 else (
                    Decision.REVISE if lipinski_violations <= 2 else Decision.KILL
                ),
                metadata={'mw': mw, 'logp': logp, 'hbd': hbd, 'hba': hba}
            ))

            if lipinski_violations > 0:
                warnings.append(
                    f"Lipinski violations: {lipinski_violations} "
                    f"(MW={mw:.0f}, LogP={logp:.1f}, HBD={hbd}, HBA={hba})"
                )

            # ===============================================================
            # STEP 2: Veber Rules (Oral Bioavailability)
            # ===============================================================

            rotatable_bonds = Lipinski.NumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)

            veber_pass = (rotatable_bonds <= 10 and tpsa <= 140)

            metrics.append(MetricResult(
                name="veber_compliance",
                value=1.0 if veber_pass else 0.0,
                unit="pass/fail",
                threshold_pass=1.0,
                threshold_revise=0.0,
                decision=Decision.PASS if veber_pass else Decision.REVISE,
                metadata={'rotatable_bonds': rotatable_bonds, 'tpsa': tpsa}
            ))

            if not veber_pass:
                warnings.append(
                    f"Veber rule violations: RotBonds={rotatable_bonds}, TPSA={tpsa:.0f} Ų"
                )

            # ===============================================================
            # STEP 3: Oral Bioavailability Estimate
            # ===============================================================

            # Simple heuristic: Lipinski + Veber compliance
            bioavailability = 0.70  # Baseline

            if lipinski_violations > 0:
                bioavailability -= 0.15 * lipinski_violations

            if not veber_pass:
                bioavailability -= 0.20

            # High lipophilicity = poor absorption
            if logp > 5:
                bioavailability -= 0.15
            elif logp < 0:
                bioavailability -= 0.10  # Too hydrophilic

            bioavailability = max(0.0, min(1.0, bioavailability))

            metrics.append(MetricResult(
                name="predicted_oral_bioavailability",
                value=bioavailability,
                unit="fraction",
                threshold_pass=0.30,
                threshold_revise=0.15,
                decision=Decision.PASS if bioavailability >= 0.30 else (
                    Decision.REVISE if bioavailability >= 0.15 else Decision.KILL
                ),
                metadata={'method': 'rule_based'}
            ))

            # ===============================================================
            # STEP 4: BBB Permeability (for CNS drugs)
            # ===============================================================

            # BBB penetration requires:
            # - MW < 450
            # - LogP 1-3
            # - TPSA < 90
            # - Low charge

            bbb_score = 1.0

            if mw > 450:
                bbb_score -= 0.30
            if logp < 1 or logp > 3:
                bbb_score -= 0.25
            if tpsa > 90:
                bbb_score -= 0.30
            if hbd + hba > 8:
                bbb_score -= 0.15

            bbb_score = max(0.0, min(1.0, bbb_score))

            metrics.append(MetricResult(
                name="bbb_permeability_score",
                value=bbb_score,
                unit="score",
                threshold_pass=0.60,
                threshold_revise=0.30,
                decision=Decision.PASS if bbb_score >= 0.60 else (
                    Decision.REVISE if bbb_score >= 0.30 else Decision.INFORMATIVE
                ),
                metadata={'note': 'Only relevant for CNS drugs'}
            ))

            # ===============================================================
            # STEP 5: Plasma Protein Binding Estimate
            # ===============================================================

            # High LogP and aromatic rings = high protein binding
            ppb = 0.50  # Baseline 50%

            if logp > 3:
                ppb += 0.15
            if logp > 5:
                ppb += 0.15

            aromatic_rings = Lipinski.NumAromaticRings(mol)
            if aromatic_rings >= 3:
                ppb += 0.10

            ppb = max(0.0, min(1.0, ppb))

            metrics.append(MetricResult(
                name="plasma_protein_binding",
                value=ppb,
                unit="fraction",
                threshold_pass=0.95,  # High PPB can be problematic
                threshold_revise=0.99,
                decision=Decision.PASS if ppb < 0.95 else (
                    Decision.REVISE if ppb < 0.99 else Decision.KILL
                ),
                metadata={'method': 'logp_based'}
            ))

            if ppb > 0.90:
                warnings.append(
                    f"High plasma protein binding predicted ({ppb:.0%}) - "
                    "may affect free drug concentration"
                )

        # ==================================================================
        # Biologics: Use sequence-based heuristics
        # ==================================================================

        else:
            sequence = candidate.get_primary_sequence()

            if not sequence:
                return self._create_result(
                    candidate=candidate,
                    decision=Decision.INFORMATIVE,
                    metrics=[],
                    summary="No sequence available for biologic PKPD",
                    warnings=["Sequence required for biologics"],
                    runtime_seconds=time.time() - start_time,
                )

            # Biologics have very different PK than small molecules
            # Key factors:
            # 1. Size (MW) - larger = slower clearance
            # 2. Charge - affects tissue distribution
            # 3. Glycosylation sites (affects half-life via FcRn)
            # 4. Aggregation propensity (affects clearance)

            mw_kda = len(sequence) * 0.11  # ~110 Da per amino acid

            # Half-life estimate (very rough)
            # Antibodies: ~21 days (FcRn recycling)
            # Peptides: minutes to hours
            # Proteins: hours to days

            if candidate.modality.value == "ANTIBODY":
                half_life_hours = 21 * 24  # ~21 days
                bioavailability = 0.70  # Good for mAbs (IV or SC)
            elif candidate.modality.value == "PEPTIDE":
                if len(sequence) < 20:
                    half_life_hours = 0.5  # Very short
                    bioavailability = 0.05  # Poor oral
                else:
                    half_life_hours = 2.0
                    bioavailability = 0.10
            else:
                # Generic protein
                half_life_hours = 12.0
                bioavailability = 0.30

            metrics.append(MetricResult(
                name="predicted_half_life",
                value=half_life_hours,
                unit="hours",
                threshold_pass=1.0,
                threshold_revise=0.5,
                decision=Decision.PASS if half_life_hours >= 1.0 else (
                    Decision.REVISE if half_life_hours >= 0.5 else Decision.KILL
                ),
                metadata={'modality': candidate.modality.value, 'method': 'heuristic'}
            ))

            metrics.append(MetricResult(
                name="predicted_bioavailability",
                value=bioavailability,
                unit="fraction",
                threshold_pass=0.50,
                threshold_revise=0.10,
                decision=Decision.PASS if bioavailability >= 0.50 else (
                    Decision.REVISE if bioavailability >= 0.10 else Decision.INFORMATIVE
                ),
                metadata={'note': 'Biologics typically require injection'}
            ))

            warnings.append(
                f"Biologic PKPD is highly complex. Pre-filter estimates: "
                f"half-life ~{half_life_hours:.1f}h, bioavailability ~{bioavailability:.0%}"
            )

        # ==================================================================
        # Overall Decision
        # ==================================================================

        decisions = [m.decision for m in metrics if m.decision]

        if Decision.KILL in decisions:
            overall_decision = Decision.KILL
            summary = "Pre-filter: Poor predicted PK properties. Likely unsuitable for development."
        elif Decision.REVISE in decisions:
            overall_decision = Decision.REVISE
            summary = "Pre-filter: Suboptimal PK properties. Clinical PKPD validation recommended."
        else:
            overall_decision = Decision.PASS
            summary = "Pre-filter: Acceptable PK properties predicted."

        # Recommendations
        recommendations = []
        risks = []

        if candidate.is_small_molecule():
            if lipinski_violations > 2:
                recommendations.append(
                    "Multiple Lipinski violations - consider structural modifications "
                    "to improve drug-likeness"
                )
                risks.append("Poor oral bioavailability expected")

            if not veber_pass:
                recommendations.append(
                    "Veber rule violations - reduce rotatable bonds or TPSA "
                    "to improve oral absorption"
                )

        warnings.append(
            "Pre-filter uses rule-based estimates only. "
            "Clinical stage will use QSAR models for accurate PK predictions."
        )

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

        logger.info(f"{self.name} validation completed: {overall_decision.value} ({runtime:.2f}s)")

        return result
