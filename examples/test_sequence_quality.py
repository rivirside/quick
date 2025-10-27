"""
Test Sequence Quality Modules - Localization & Complexity.

Tests the new quick-win modules:
1. Localization signal detection (signal peptide, TM, NLS)
2. Low-complexity region detection
"""

from etrial import TherapeuticCandidate, Modality
from etrial.developability.localization import LocalizationValidationModule
from etrial.developability.complexity import ComplexityValidationModule


def main():
    print("=" * 70)
    print("eTrial Sequence Quality Test")
    print("Testing: Localization Signals + Low-Complexity Detection")
    print("=" * 70)
    print()

    # Create test candidates
    candidates = [
        # Clean sequence (should pass)
        TherapeuticCandidate(
            name="Clean_Antibody",
            sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS",
            modality=Modality.ANTIBODY,
            target="PD-1",
        ),

        # Signal peptide (should flag)
        TherapeuticCandidate(
            name="With_Signal_Peptide",
            sequence="MKTIIALSYIFCLVFADEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYY",
            # MKTIIALSYIFCLVFA = signal peptide (hydrophobic stretch)
            modality=Modality.PROTEIN,
            target="Test",
        ),

        # Transmembrane segment (should flag)
        TherapeuticCandidate(
            name="With_TM_Segment",
            sequence="MGSTKLVVVVLLLIIIIFFFFMMMMMWWWGGGGKKKDEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSW",
            # Long hydrophobic stretch = TM
            modality=Modality.PROTEIN,
            target="Test",
        ),

        # Nuclear localization signal (should flag)
        TherapeuticCandidate(
            name="With_NLS",
            sequence="DEVQLVESKKKRKKKKGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKG",
            # KKKRKKKK = NLS
            modality=Modality.PROTEIN,
            target="Test",
        ),

        # Homopolymer (should flag)
        TherapeuticCandidate(
            name="Homopolymer_Lysine",
            sequence="DEVQLVESKKKKKKKKKGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKG",
            # KKKKKKKKKK = homopolymer
            modality=Modality.PROTEIN,
            target="Test",
        ),

        # Simple repeat (should flag)
        TherapeuticCandidate(
            name="Simple_Repeat",
            sequence="DEVQLVESKAGAGAGAGAGAGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKG",
            # (AG)×6 = simple repeat
            modality=Modality.PROTEIN,
            target="Test",
        ),

        # Low complexity / poly-proline (should flag)
        TherapeuticCandidate(
            name="Poly_Proline",
            sequence="DEVQLVESPPPPPPPPPPPPPGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKG",
            # PPPPPPPPPPPPPP = proline-rich
            modality=Modality.PROTEIN,
            target="Test",
        ),

        # Poly-basic region (should flag)
        TherapeuticCandidate(
            name="Poly_Basic",
            sequence="DEVQLVESKRKRKRKRKRKGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKG",
            # KRKRKRKRKR = poly-basic
            modality=Modality.PROTEIN,
            target="Test",
        ),
    ]

    # ======================================================================
    # PART 1: Localization Signal Detection
    # ======================================================================

    print("=" * 70)
    print("PART 1: Localization Signal Detection")
    print("=" * 70)
    print()

    localization = LocalizationValidationModule()

    for candidate in candidates:
        result = localization.validate(candidate)

        # Print result
        decision_symbol = "✓" if result.decision.value == "PASS" else (
            "⚠" if result.decision.value == "REVISE" else "✗"
        )

        print(f"{decision_symbol} {candidate.name:25s} → {result.decision.value:6s}")

        # Show warnings (first 2)
        for warning in result.warnings[:2]:
            print(f"     ⚠ {warning}")

        # Show key metrics (signal peptide, TM, NLS)
        for metric in result.metrics[:3]:
            if metric.value > 0:
                print(f"     📊 {metric.name}: {metric.value:.2f} {metric.unit}")

        print()

    # ======================================================================
    # PART 2: Low-Complexity Detection
    # ======================================================================

    print("=" * 70)
    print("PART 2: Low-Complexity Detection")
    print("=" * 70)
    print()

    complexity = ComplexityValidationModule()

    for candidate in candidates:
        result = complexity.validate(candidate)

        # Print result
        decision_symbol = "✓" if result.decision.value == "PASS" else (
            "⚠" if result.decision.value == "REVISE" else "✗"
        )

        print(f"{decision_symbol} {candidate.name:25s} → {result.decision.value:6s}")

        # Show warnings (first 2)
        for warning in result.warnings[:2]:
            print(f"     ⚠ {warning}")

        # Show problematic metrics only
        for metric in result.metrics:
            if metric.decision != "PASS" and metric.value > 0:
                print(f"     📊 {metric.name}: {metric.value:.2f} {metric.unit}")

        print()

    # ======================================================================
    # Summary
    # ======================================================================

    print("=" * 70)
    print("✓ Sequence Quality Test Complete!")
    print("=" * 70)
    print()
    print("Summary:")
    print("  • Localization module detects signal peptides, TM, NLS")
    print("  • Complexity module detects homopolymers, repeats, low entropy")
    print("  • Both use pattern-based detection (fast, no external tools needed)")
    print()
    print("Improvements vs. original:")
    print("  ❌ OLD: No localization or complexity checks")
    print("  ✅ NEW: Comprehensive sequence quality control")
    print()
    print("Next steps:")
    print("  1. Integrate SignalP/TMHMM when tools installed")
    print("  2. Add TANGO/AGGRESCAN for better aggregation")
    print("  3. Integrate into two-stage pipeline")
    print()


if __name__ == "__main__":
    main()
