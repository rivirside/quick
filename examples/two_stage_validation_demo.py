#!/usr/bin/env python3
"""
Demo: Two-Stage Validation Pipeline

Shows the difference between:
1. Pre-filter stage (fast, sequence-based)
2. Clinical-grade stage (slow, validated ML)
"""

from pathlib import Path
from etrial import TherapeuticCandidate, Modality
from etrial.core.two_stage_pipeline import create_two_stage_pipeline

def create_test_candidates():
    """Create a diverse set of test candidates."""
    return [
        # Good candidate
        TherapeuticCandidate(
            name="Good_Antibody",
            modality=Modality.ANTIBODY,
            target="HER2",
            heavy_chain="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS",
            metadata={"expected": "PASS"}
        ),

        # Aggregation risk
        TherapeuticCandidate(
            name="Aggregation_Risk",
            modality=Modality.PEPTIDE,
            target="test",
            sequence="VVVVVVIIIIILLLLLLFFFFVVVVVIIIIILLLLLFFFFFF",
            metadata={"expected": "KILL"}
        ),

        # Free cysteine
        TherapeuticCandidate(
            name="Free_Cysteine",
            modality=Modality.PEPTIDE,
            target="test",
            sequence="MKTAYIAKCRQISSQRKRK",
            metadata={"expected": "KILL"}
        ),

        # Minor PTM issues
        TherapeuticCandidate(
            name="PTM_Issues",
            modality=Modality.PEPTIDE,
            target="test",
            sequence="MKTAYIAKQRNGSQRKNNKR",
            metadata={"expected": "REVISE"}
        ),

        # Trastuzumab (approved drug)
        TherapeuticCandidate(
            name="Trastuzumab",
            modality=Modality.ANTIBODY,
            target="HER2",
            heavy_chain="EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
            light_chain="DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
            metadata={"expected": "PASS", "approved": True}
        ),
    ]


def demo_prefilter_only():
    """Demo: Pre-filter stage only (fast screening)."""
    print("\n" + "=" * 70)
    print("DEMO 1: PRE-FILTER ONLY (Fast Sequence Screening)")
    print("=" * 70)

    # Create pipeline in pre-filter mode
    pipeline = create_two_stage_pipeline(
        mode='prefilter',
        config_path=Path("config/default_config.yaml")
    )

    # Create test candidates
    candidates = create_test_candidates()

    print(f"\nScreening {len(candidates)} candidates...")

    # Run validation
    results = pipeline.validate_batch(candidates)

    # Print results
    print("\nResults:")
    for result in results['final_candidates']:
        expected = result.candidate.metadata.get('expected', 'Unknown')
        match = "✓" if result.overall_decision.value in [expected] or expected == "Unknown" else "✗"
        print(f"  {match} {result.candidate.name:20s} → {result.overall_decision.value:8s} (expected: {expected})")

    # Statistics
    stats = pipeline.get_stage_statistics(results)
    print(f"\nStatistics:")
    print(f"  Total runtime: {stats['total_runtime_seconds']:.2f}s")
    print(f"  Throughput: {stats['stages']['prefilter']['throughput_per_sec']:.1f} candidates/sec")
    print(f"  Pass rate: {stats['stages']['prefilter']['pass_rate']}")


def demo_both_stages():
    """Demo: Both stages (full workflow)."""
    print("\n" + "=" * 70)
    print("DEMO 2: BOTH STAGES (Pre-Filter → Clinical-Grade)")
    print("=" * 70)

    # Create pipeline in 'both' mode
    pipeline = create_two_stage_pipeline(
        mode='both',
        config_path=Path("config/default_config.yaml")
    )

    # Create test candidates
    candidates = create_test_candidates()

    print(f"\nValidating {len(candidates)} candidates through both stages...")

    # Run validation
    results = pipeline.validate_batch(candidates)

    # Print stage-by-stage results
    print("\nStage-by-Stage Flow:")
    for stage_result in results['stages']:
        print(f"\n  {stage_result.stage.upper()}:")
        print(f"    Input: {stage_result.candidates_in} candidates")
        print(f"    Output: {stage_result.candidates_out} candidates")
        print(f"    Pass rate: {stage_result.pass_rate*100:.1f}%")
        print(f"    Runtime: {stage_result.runtime_seconds:.2f}s")

    # Final results
    print("\nFinal Candidates:")
    for result in results['final_candidates']:
        if result.overall_decision != 'KILL':
            print(f"  • {result.candidate.name:20s} → {result.overall_decision.value}")

    # Overall statistics
    stats = pipeline.get_stage_statistics(results)
    print(f"\nOverall Statistics:")
    print(f"  Total runtime: {stats['total_runtime_seconds']:.2f}s")
    print(f"  Final decision breakdown: {stats.get('final_decision_breakdown', {})}")


def demo_clinical_only():
    """Demo: Clinical-grade only (detailed analysis of pre-filtered candidates)."""
    print("\n" + "=" * 70)
    print("DEMO 3: CLINICAL-GRADE ONLY (Detailed Analysis)")
    print("=" * 70)

    # Create pipeline in clinical mode
    pipeline = create_two_stage_pipeline(
        mode='clinical',
        config_path=Path("config/default_config.yaml")
    )

    # Use only the good candidates (simulating pre-filter output)
    candidates = [
        TherapeuticCandidate(
            name="Trastuzumab",
            modality=Modality.ANTIBODY,
            target="HER2",
            heavy_chain="EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
            light_chain="DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
        ),
    ]

    print(f"\nRunning clinical-grade validation on {len(candidates)} pre-filtered candidates...")
    print("(Note: Clinical modules not yet implemented - this will be slow when added)")

    # Run validation
    results = pipeline.validate_batch(candidates)

    print("\nResults:")
    for result in results['final_candidates']:
        print(f"  {result.candidate.name}: {result.overall_decision.value}")


def main():
    """Run all demos."""
    print("\n" + "🧬" * 35)
    print("eTrial Two-Stage Validation Pipeline Demo")
    print("🧬" * 35)

    try:
        # Demo 1: Pre-filter only (fast)
        demo_prefilter_only()

        # Demo 2: Both stages (full workflow)
        demo_both_stages()

        # Demo 3: Clinical only (detailed)
        demo_clinical_only()

        print("\n" + "=" * 70)
        print("✓ All demos complete!")
        print("=" * 70)

        print("\nNext Steps:")
        print("  1. Add clinical-grade modules (CamSol, AlphaFold, etc.)")
        print("  2. Integrate with your own data")
        print("  3. Benchmark against experimental results")

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
