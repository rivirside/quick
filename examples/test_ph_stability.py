"""
Test pH stability module.

Tests:
1. Normal antibody (good pH stability)
2. Histidine-rich sequence (pH-dependent)
3. Buried charged residues (pH-labile)
4. Acidic sequence (low pI)
5. Basic sequence (high pI)
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from etrial.core.base import TherapeuticCandidate, Modality, Decision
from etrial.developability.ph_stability import pHStabilityModule


def print_result(name: str, result):
    """Pretty print validation result."""
    print(f"\n{'=' * 70}")
    print(f"Test: {name}")
    print(f"{'=' * 70}")
    print(f"Decision: {result.decision.value}")
    print(f"Summary: {result.summary}")

    if result.metrics:
        print(f"\nMetrics:")
        for metric in result.metrics:
            status = "✅" if metric.decision == Decision.PASS else "ℹ️" if metric.decision == Decision.INFORMATIVE else "⚠️" if metric.decision == Decision.REVISE else "❌"
            print(f"  {status} {metric.name}: {metric.value:.1f} {metric.unit}")

    if result.recommendations:
        print(f"\nRecommendations:")
        for rec in result.recommendations:
            print(f"  → {rec}")

    if result.risks:
        print(f"\nRisks:")
        for risk in result.risks:
            print(f"  ❌ {risk}")

    if result.warnings:
        print(f"\nWarnings:")
        for warning in result.warnings:
            print(f"  ⚠️  {warning}")

    # Show pH-labile residues
    if 'pka_results' in result.metadata:
        labile = [r for r in result.metadata['pka_results'] if r['labile']]
        if labile:
            print(f"\n⚠️  pH-Labile Residues:")
            for r in labile[:5]:
                print(f"  {r['residue']}{r['position']}: pKa = {r['predicted_pka']:.1f} (shift: {r['pka_shift']:+.1f})")

    print()


def main():
    """Run pH stability tests."""

    module = pHStabilityModule()

    # Test 1: Adalimumab VH (normal antibody)
    adalimumab_vh = TherapeuticCandidate(
        name="Adalimumab_VH",
        sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        modality=Modality.ANTIBODY,
        target="TNF-alpha"
    )
    result1 = module.validate(adalimumab_vh)
    print_result("Adalimumab VH (Normal)", result1)

    # Test 2: Histidine-rich sequence (pH-dependent)
    his_rich = TherapeuticCandidate(
        name="Histidine_Rich",
        sequence="HHHHHAAAHHHHHAAAHHHHHAAAHHHHHGGG",
        modality=Modality.PEPTIDE,
        target="Test"
    )
    result2 = module.validate(his_rich)
    print_result("Histidine-Rich Peptide", result2)

    # Test 3: Buried charged residues (simulated with hydrophobic context)
    buried_charged = TherapeuticCandidate(
        name="Buried_Charged",
        sequence="FFFFDFFFFKFFFFFEFFFFFRF",
        modality=Modality.PEPTIDE,
        target="Test"
    )
    result3 = module.validate(buried_charged)
    print_result("Buried Charged Residues", result3)

    # Test 4: Acidic sequence (many D/E, low pI)
    acidic = TherapeuticCandidate(
        name="Acidic_Protein",
        sequence="DDDEEEDDDEEEDDDEEEDDDEEE",
        modality=Modality.PROTEIN,
        target="Test"
    )
    result4 = module.validate(acidic)
    print_result("Acidic Protein (Low pI)", result4)

    # Test 5: Basic sequence (many K/R, high pI)
    basic = TherapeuticCandidate(
        name="Basic_Protein",
        sequence="KKKRRRKKKKRRRKKKKRRRKKK",
        modality=Modality.PROTEIN,
        target="Test"
    )
    result5 = module.validate(basic)
    print_result("Basic Protein (High pI)", result5)

    # Test 6: Balanced sequence
    balanced = TherapeuticCandidate(
        name="Balanced_Protein",
        sequence="AVDEFGHIKLMNPQRSTVWYAVDEFGHIKLMNPQRSTVWY",
        modality=Modality.PROTEIN,
        target="Test"
    )
    result6 = module.validate(balanced)
    print_result("Balanced Protein", result6)

    # Test 7: Small molecule (should skip)
    small_molecule = TherapeuticCandidate(
        name="Aspirin",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        modality=Modality.SMALL_MOLECULE,
        target="COX"
    )
    result7 = module.validate(small_molecule)
    print_result("Small Molecule (Should Skip)", result7)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    tests = [
        ("Adalimumab VH", result1, Decision.PASS),
        ("Histidine-Rich", result2, Decision.PASS),  # His clusters OK if not buried
        ("Buried Charged", result3, Decision.REVISE),  # Some labile residues
        ("Acidic Protein", result4, Decision.KILL),  # Clustered charges = labile
        ("Basic Protein", result5, Decision.KILL),  # Clustered charges = labile
        ("Balanced Protein", result6, Decision.PASS),
        ("Small Molecule (skip)", result7, Decision.INFORMATIVE),
    ]

    passed = 0
    for name, result, expected in tests:
        match = "✅" if result.decision == expected else "❌"
        print(f"{match} {name}: {result.decision.value} (expected {expected.value})")
        if result.decision == expected:
            passed += 1

    print(f"\nPassed: {passed}/{len(tests)}")

    if passed == len(tests):
        print("\n🎉 All tests passed!")
    else:
        print(f"\n⚠️  {len(tests) - passed} test(s) failed")


if __name__ == "__main__":
    main()
