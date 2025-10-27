"""
Test MHC-I epitope prediction module.

Tests:
1. Normal antibody (Adalimumab VH - should have some epitopes but acceptable)
2. Epitope-rich sequence (peptides with strong HLA-A*02:01 motifs)
3. Humanized antibody (should skip or pass with low risk)
4. Weak binding sequence (no strong motifs)
5. Small molecule (should skip)
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from etrial.core.base import TherapeuticCandidate, Modality, Decision
from etrial.immunogenicity.mhc1_epitopes import MHC1EpitopeModule


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

    # Show top epitopes
    if 'epitopes' in result.metadata and result.metadata['epitopes']:
        epitopes = result.metadata['epitopes']
        strong_epitopes = [e for e in epitopes if e['strong']]
        if strong_epitopes:
            print(f"\n⚠️  Top MHC-I Epitopes (Strong Binders):")
            for e in strong_epitopes[:5]:
                print(f"  {e['peptide']} (pos {e['position']}, {e['hla']}, rank {e['percentile']:.2f}%)")

    print()


def main():
    """Run MHC-I epitope prediction tests."""

    module = MHC1EpitopeModule()

    # Test 1: Adalimumab VH (normal antibody)
    adalimumab_vh = TherapeuticCandidate(
        name="Adalimumab_VH",
        sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        modality=Modality.ANTIBODY,
        target="TNF-alpha"
    )
    result1 = module.validate(adalimumab_vh)
    print_result("Adalimumab VH (Normal Antibody)", result1)

    # Test 2: Epitope-rich sequence (strong HLA-A*02:01 motifs)
    # HLA-A*02:01: P2 = L,M,I,V and P9 = L,I,V,A
    # Create peptides with strong motifs
    epitope_rich = TherapeuticCandidate(
        name="Epitope_Rich",
        sequence="MLVVVVVVLMLVVVVVVLMLVVVVVVLMLVVVVVVL",  # Many L,M,V at P2 and P9 positions
        modality=Modality.PEPTIDE,
        target="Test"
    )
    result2 = module.validate(epitope_rich)
    print_result("Epitope-Rich Peptide", result2)

    # Test 3: Humanized antibody (should pass with low risk)
    humanized = TherapeuticCandidate(
        name="Humanized_Antibody",
        sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        modality=Modality.ANTIBODY,
        target="CD20",
        metadata={'design_notes': 'Humanized from mouse anti-CD20'}
    )
    result3 = module.validate(humanized)
    print_result("Humanized Antibody", result3)

    # Test 4: Weak binding sequence (no strong motifs)
    # Avoid L,M,I,V at P2 and P9 positions
    weak_binding = TherapeuticCandidate(
        name="Weak_Binding",
        sequence="GGGAAAGGGAAAGGGAAAGGGAAAGGGAAA",  # G,A are weak binders
        modality=Modality.PEPTIDE,
        target="Test"
    )
    result4 = module.validate(weak_binding)
    print_result("Weak Binding Peptide", result4)

    # Test 5: Moderate epitopes (some binding)
    moderate = TherapeuticCandidate(
        name="Moderate_Epitopes",
        sequence="KLFVVVILYKLSLILRYKLHLFLKY",  # Mix of motifs
        modality=Modality.PEPTIDE,
        target="Test"
    )
    result5 = module.validate(moderate)
    print_result("Moderate Epitopes", result5)

    # Test 6: Small molecule (should skip)
    small_molecule = TherapeuticCandidate(
        name="Aspirin",
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        modality=Modality.SMALL_MOLECULE,
        target="COX"
    )
    result6 = module.validate(small_molecule)
    print_result("Small Molecule (Should Skip)", result6)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    tests = [
        ("Adalimumab VH", result1, Decision.PASS),  # Some epitopes but acceptable
        ("Epitope-Rich", result2, Decision.KILL),  # Too many epitopes (25 strong binders > 10 threshold)
        ("Humanized Antibody", result3, Decision.PASS),  # Humanized = low risk
        ("Weak Binding", result4, Decision.PASS),  # Few/no epitopes
        ("Moderate Epitopes", result5, Decision.PASS),  # Some epitopes but acceptable
        ("Small Molecule", result6, Decision.INFORMATIVE),  # Should skip
    ]

    passed = 0
    for name, result, expected in tests:
        # Allow some flexibility for PASS vs REVISE (depends on exact epitope count)
        if name in ["Adalimumab VH", "Moderate Epitopes"]:
            # Could be PASS or REVISE depending on epitope count
            match_ok = result.decision in [Decision.PASS, Decision.REVISE]
            match = "✅" if match_ok else "❌"
            print(f"{match} {name}: {result.decision.value} (expected {expected.value} or REVISE)")
            if match_ok:
                passed += 1
        else:
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
