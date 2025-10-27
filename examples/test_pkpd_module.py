"""
Test PKPD Module - Verify Pharmacokinetics Predictions.

Tests both pre-filter and clinical PKPD modules with various test compounds.
"""

from etrial import TherapeuticCandidate, Modality
from etrial.pkpd.pharmacokinetics import PKPDPrefilterModule
from etrial.pkpd.pharmacokinetics_clinical import ClinicalPKPDModule


def main():
    print("=" * 70)
    print("eTrial PKPD Module Test")
    print("=" * 70)
    print()

    # Create test candidates
    candidates = [
        # Drug-like small molecule (good oral bioavailability)
        TherapeuticCandidate(
            name="Aspirin",
            smiles="CC(=O)Oc1ccccc1C(=O)O",
            modality=Modality.SMALL_MOLECULE,
            target="COX-2",
        ),

        # Large lipophilic molecule (poor oral bioavailability)
        TherapeuticCandidate(
            name="Large_Lipophilic_Drug",
            smiles="CC(C)Cc1ccc(C(C)C(=O)O)cc1C1=C(c2ccc(C(C)C)cc2)c2ccccc2C(c2ccccc2)C1c1ccccc1",
            modality=Modality.SMALL_MOLECULE,
            target="Test",
        ),

        # Small hydrophilic molecule (poor membrane permeability)
        TherapeuticCandidate(
            name="Hydrophilic_Peptide_Mimic",
            smiles="NC(=O)C(N)CCCCNC(=N)N",  # Arginine analog
            modality=Modality.SMALL_MOLECULE,
            target="Test",
        ),

        # Antibody (very long half-life)
        TherapeuticCandidate(
            name="Therapeutic_Antibody",
            sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS",
            modality=Modality.ANTIBODY,
            target="PD-1",
        ),

        # Short peptide (very short half-life)
        TherapeuticCandidate(
            name="Short_Peptide",
            sequence="MKFLKFSL",
            modality=Modality.PEPTIDE,
            target="Test",
        ),
    ]

    # ======================================================================
    # PART 1: Pre-Filter PKPD (Fast Screening)
    # ======================================================================

    print("=" * 70)
    print("PART 1: Pre-Filter PKPD Screening (Fast)")
    print("=" * 70)
    print()

    prefilter = PKPDPrefilterModule()

    for candidate in candidates:
        result = prefilter.validate(candidate)

        # Print result
        decision_symbol = "✓" if result.decision.value == "PASS" else (
            "⚠" if result.decision.value == "REVISE" else "✗"
        )

        print(f"{decision_symbol} {candidate.name:30s} → {result.decision.value:6s}")

        # Show warnings (first 2 only)
        for warning in result.warnings[:2]:
            print(f"     ⚠ {warning}")

        # Show key metrics
        for metric in result.metrics[:4]:  # Show top 4 metrics
            print(f"     📊 {metric.name}: {metric.value:.2f} {metric.unit}")

        print()

    # ======================================================================
    # PART 2: Clinical PKPD (Accurate Predictions)
    # ======================================================================

    print("=" * 70)
    print("PART 2: Clinical-Grade PKPD Prediction")
    print("=" * 70)
    print()

    clinical = ClinicalPKPDModule()

    for candidate in candidates:
        result = clinical.validate(candidate)

        # Print result
        decision_symbol = "✓" if result.decision.value == "PASS" else (
            "⚠" if result.decision.value == "REVISE" else "✗"
        )

        print(f"{decision_symbol} {candidate.name:30s} → {result.decision.value:6s}")

        # Show warnings
        for warning in result.warnings[:2]:
            print(f"     ⚠ {warning}")

        # Show ALL metrics for clinical
        for metric in result.metrics:
            print(f"     📊 {metric.name}: {metric.value:.2f} {metric.unit}")

        # Show recommendations
        for rec in result.recommendations[:1]:
            print(f"     💡 {rec[:80]}...")

        print()

    # ======================================================================
    # Summary
    # ======================================================================

    print("=" * 70)
    print("✓ PKPD Module Test Complete!")
    print("=" * 70)
    print()
    print("Summary:")
    print("  • Pre-filter uses Lipinski, Veber rules for fast screening")
    print("  • Clinical uses QSAR models (Lombardo 2018, Obach 2008)")
    print("  • Biologics use empirical PK models (FcRn, proteolysis)")
    print()
    print("Next steps:")
    print("  1. Integrate into two-stage pipeline")
    print("  2. Validate against experimental PK data")
    print("  3. Tune QSAR coefficients if needed")
    print()


if __name__ == "__main__":
    main()
