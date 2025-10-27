"""
Test alanine scanning module.

Tests:
1. Known antibody with good binding interface
2. Peptide with few hotspots
3. Short sequence (should still work)
4. Sequence with many aromatics/charged residues (many hotspots)
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from etrial.core.base import TherapeuticCandidate, Modality, Decision
from etrial.structure.alanine_scanning import AlanineScanningModule


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
            status = "✅" if metric.decision == Decision.PASS else "ℹ️" if metric.decision == Decision.INFORMATIVE else "❌"
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

    # Show top hotspots
    if 'hotspots' in result.metadata:
        hotspots = [h for h in result.metadata['hotspots'] if h['hotspot']]
        if hotspots:
            print(f"\n🔥 Top Hotspots:")
            for h in sorted(hotspots, key=lambda x: x['ddg'], reverse=True)[:5]:
                print(f"  {h['aa']}{h['position']}: ΔΔG = {h['ddg']:.1f} kcal/mol")

    print()


def main():
    """Run alanine scanning tests."""

    module = AlanineScanningModule()

    # Test 1: Adalimumab VH (good antibody with proper binding interface)
    adalimumab_vh = TherapeuticCandidate(
        name="Adalimumab_VH",
        sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        modality=Modality.ANTIBODY,
        target="TNF-alpha"
    )
    result1 = module.validate(adalimumab_vh)
    print_result("Adalimumab VH (FDA-approved)", result1)

    # Test 2: Peptide with few hot spots (mostly Ala/Gly)
    weak_peptide = TherapeuticCandidate(
        name="Weak_Peptide",
        sequence="AAAGGGAAAGGGAAAGGGAAAGGGA",
        modality=Modality.PEPTIDE,
        target="Test"
    )
    result2 = module.validate(weak_peptide)
    print_result("Weak Peptide (mostly Ala/Gly)", result2)

    # Test 3: Peptide rich in hotspot residues (Arg, Trp, Tyr, Phe)
    hotspot_rich = TherapeuticCandidate(
        name="Hotspot_Rich_Peptide",
        sequence="RWYFRHWYRFHWRYFHWRYF",
        modality=Modality.PEPTIDE,
        target="Test"
    )
    result3 = module.validate(hotspot_rich)
    print_result("Hotspot-Rich Peptide (Arg/Trp/Tyr/Phe)", result3)

    # Test 4: Natural CDR sequence with mixed composition
    cdr_sequence = TherapeuticCandidate(
        name="CDR3_Sequence",
        sequence="ARGNYYDSSGYFDYWGQGTLVTVSS",
        modality=Modality.PEPTIDE,
        target="Test"
    )
    result4 = module.validate(cdr_sequence)
    print_result("CDR3 Sequence (Mixed)", result4)

    # Test 5: Small molecule (should skip)
    small_molecule = TherapeuticCandidate(
        name="Aspirin",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        modality=Modality.SMALL_MOLECULE,
        target="COX"
    )
    result5 = module.validate(small_molecule)
    print_result("Small Molecule (Should Skip)", result5)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    tests = [
        ("Adalimumab VH", result1, Decision.PASS),
        ("Weak Peptide", result2, Decision.KILL),
        ("Hotspot-Rich Peptide", result3, Decision.PASS),
        ("CDR3 Sequence", result4, Decision.PASS),
        ("Small Molecule (skip)", result5, Decision.INFORMATIVE),
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
