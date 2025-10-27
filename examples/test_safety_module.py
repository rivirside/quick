"""
Test Safety Module - Toxicity Screening Demo.

Tests both pre-filter and clinical-grade toxicity modules.
"""

from etrial import TherapeuticCandidate, Modality
from etrial.safety.toxicity import ToxicityPrefilterModule
from etrial.safety.toxicity_clinical import ClinicalToxicityModule


def main():
    print("=" * 70)
    print("eTrial Safety Module Test")
    print("=" * 70)
    print()

    # ========================================
    # Test Candidates
    # ========================================

    candidates = []

    # 1. Good drug (Aspirin)
    candidates.append(TherapeuticCandidate(
        name="Aspirin",
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        modality=Modality.SMALL_MOLECULE,
        target="COX-1/COX-2",
    ))

    # 2. Known toxic (contains azo group)
    candidates.append(TherapeuticCandidate(
        name="Toxic_Azo_Compound",
        smiles="c1ccc(N=Nc2ccccc2)cc1",  # Azobenzene (azo dye, toxic)
        modality=Modality.SMALL_MOLECULE,
        target="Test",
    ))

    # 3. Lipinski violations (too large)
    candidates.append(TherapeuticCandidate(
        name="Too_Large_Molecule",
        smiles="CC(C)CC1=CC=C(C=C1)C(C)C(=O)OC2=CC=CC=C2C(=O)O" + "C" * 20,  # Extended
        modality=Modality.SMALL_MOLECULE,
        target="Test",
    ))

    # 4. Biologic with poly-basic region (cytotoxic peptide)
    candidates.append(TherapeuticCandidate(
        name="Poly_Basic_Peptide",
        sequence="MKKKKKKKRRRRRRLVFMWY",  # Poly-lysine/arginine region
        modality=Modality.PEPTIDE,
        target="Cell membrane",
    ))

    # 5. Safe biologic (normal antibody fragment)
    candidates.append(TherapeuticCandidate(
        name="Safe_Antibody",
        sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS",
        modality=Modality.ANTIBODY,
        target="Antigen",
    ))

    # ========================================
    # Test Pre-Filter Module
    # ========================================

    print("=" * 70)
    print("PART 1: Pre-Filter Toxicity Screening (Fast)")
    print("=" * 70)
    print()

    prefilter = ToxicityPrefilterModule()

    for candidate in candidates:
        result = prefilter.validate(candidate)

        status = "✓" if result.decision.value == "PASS" else ("⚠" if result.decision.value == "REVISE" else "✗")

        print(f"{status} {candidate.name:25s} → {result.decision.value:7s}")

        if result.risks:
            for risk in result.risks:
                print(f"     ⚠ {risk}")

        if result.warnings:
            for warning in result.warnings:
                print(f"     ⚡ {warning}")

        # Show key metrics
        for metric in result.metrics:
            if metric.value > 0:
                print(f"     📊 {metric.name}: {metric.value} {metric.unit}")

        print()

    # ========================================
    # Test Clinical Module
    # ========================================

    print("=" * 70)
    print("PART 2: Clinical-Grade Toxicity Prediction")
    print("=" * 70)
    print()

    clinical = ClinicalToxicityModule()

    print("Note: Clinical APIs not yet integrated - using RDKit-based estimates")
    print()

    for candidate in candidates:
        result = clinical.validate(candidate)

        status = "✓" if result.decision.value == "PASS" else ("⚠" if result.decision.value == "REVISE" else "✗")

        print(f"{status} {candidate.name:25s} → {result.decision.value:7s}")

        if result.risks:
            for risk in result.risks:
                print(f"     ⚠ {risk}")

        # Show key metrics
        for metric in result.metrics:
            print(f"     📊 {metric.name}: {metric.value:.2f} {metric.unit}")

        if result.recommendations:
            for rec in result.recommendations[:1]:  # Show first recommendation only
                print(f"     💡 {rec[:80]}...")

        print()

    # ========================================
    # Summary
    # ========================================

    print("=" * 70)
    print("✓ Safety Module Test Complete!")
    print("=" * 70)
    print()
    print("Summary:")
    print("  • Pre-filter uses RDKit for fast structural alerts")
    print("  • Clinical module uses web APIs (pkCSM, ProTox-II) when available")
    print("  • Both modules gracefully fallback to simpler methods")
    print()
    print("Next steps:")
    print("  1. Integrate pkCSM API for clinical predictions")
    print("  2. Test with larger compound libraries")
    print("  3. Validate against experimental toxicity data")
    print()


if __name__ == "__main__":
    main()
