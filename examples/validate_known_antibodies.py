#!/usr/bin/env python3
"""
Validate Known FDA-Approved Antibodies

This script demonstrates the pipeline workflow using real antibodies:
- Adalimumab (Humira)
- Trastuzumab (Herceptin)
- Pembrolizumab (Keytruda)
- Nivolumab (Opdivo)
- Bevacizumab (Avastin)
- Rituximab (Rituxan)

Expected: These should mostly PASS (they're FDA-approved!)
"""

import sys
from pathlib import Path
from Bio import SeqIO

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from etrial.core.base import TherapeuticCandidate, Modality, Decision
from etrial.core.two_stage_pipeline import create_two_stage_pipeline


def main():
    print("=" * 70)
    print("eTrial: Validating Known FDA-Approved Antibodies")
    print("=" * 70)
    print()

    # Load known antibodies from FASTA
    fasta_file = Path(__file__).parent.parent / "tests" / "test_data" / "known_antibodies.fasta"

    if not fasta_file.exists():
        print(f"ERROR: FASTA file not found: {fasta_file}")
        return 1

    print(f"Loading antibodies from: {fasta_file}")
    print()

    # Parse FASTA file
    records = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"Found {len(records)} antibody chains")
    print()

    # Create candidates (only heavy chains for this demo)
    candidates = []
    for record in records:
        # Parse header: >Name|TradeName|Target|Status
        parts = record.description.split('|')
        name = parts[0]

        # Only process heavy chains (for simplicity)
        if 'Heavy' not in name:
            continue

        trade_name = parts[1] if len(parts) > 1 else "Unknown"
        target = parts[2] if len(parts) > 2 else "Unknown"

        candidate = TherapeuticCandidate(
            name=name,
            sequence=str(record.seq),
            modality=Modality.ANTIBODY,
            target=target,
            metadata={
                "trade_name": trade_name,
                "status": "FDA-approved"
            }
        )
        candidates.append(candidate)

    print(f"Validating {len(candidates)} antibody heavy chains:")
    for c in candidates:
        print(f"  - {c.name} ({c.metadata['trade_name']})")
    print()

    # Create pre-filter pipeline
    print("Creating validation pipeline (pre-filter mode)...")
    pipeline = create_two_stage_pipeline(mode='prefilter')
    print()

    # Run validation
    print("Running validation...")
    print("-" * 70)
    results = pipeline.validate_batch(candidates)
    print("-" * 70)
    print()

    # Summary statistics
    stage_result = results['stages'][0]
    print(f"Validation Complete!")
    print(f"  Input: {stage_result.candidates_in} candidates")
    print(f"  Output: {stage_result.candidates_out} passed")
    print(f"  Pass Rate: {stage_result.pass_rate*100:.1f}%")
    print(f"  Runtime: {stage_result.runtime_seconds:.2f}s")
    print(f"  Throughput: {stage_result.candidates_in/stage_result.runtime_seconds:.1f} candidates/sec")
    print()

    # Detailed results
    print("=" * 70)
    print("Detailed Results")
    print("=" * 70)
    print()

    for result in results['final_candidates']:
        trade_name = result.candidate.metadata.get('trade_name', 'Unknown')
        decision_icon = "✅" if result.overall_decision == Decision.PASS else "⚠️" if result.overall_decision == Decision.REVISE else "❌"

        print(f"{decision_icon} {result.candidate.name}")
        print(f"   Trade Name: {trade_name}")
        print(f"   Target: {result.candidate.target}")
        print(f"   Decision: {result.overall_decision.value}")
        print(f"   Runtime: {result.total_runtime_seconds:.3f}s")

        # Show any failures or warnings
        failed_modules = []
        warning_modules = []

        for module_name, module_result in result.module_results.items():
            if module_result.decision == Decision.KILL:
                failed_modules.append(module_name)
            elif module_result.decision == Decision.REVISE:
                warning_modules.append(module_name)

        if failed_modules:
            print(f"   ❌ Failed: {', '.join(failed_modules)}")
        if warning_modules:
            print(f"   ⚠️  Warnings: {', '.join(warning_modules)}")

        # Show key metrics
        metrics_to_show = [
            ('viscosity_prediction', 'predicted_viscosity', 'Viscosity', 'cP'),
            ('viscosity_prediction', 'sap_score', 'SAP Score', ''),
            ('structure_affinity', 'predicted_kd', 'KD', 'nM'),
            ('advanced_aggregation', 'tango_aggregation_score', 'TANGO', ''),
        ]

        print(f"   📊 Key Metrics:")
        for module, metric, label, unit in metrics_to_show:
            mod_result = result.module_results.get(module)
            if mod_result:
                value = mod_result.get_metric_value(metric)
                if value is not None:
                    unit_str = f" {unit}" if unit else ""
                    print(f"      {label}: {value:.1f}{unit_str}")

        print()

    # Summary by decision
    print("=" * 70)
    print("Summary by Decision")
    print("=" * 70)

    decisions = {}
    for result in results['final_candidates']:
        dec = result.overall_decision.value
        decisions[dec] = decisions.get(dec, 0) + 1

    for decision, count in decisions.items():
        pct = count / len(results['final_candidates']) * 100
        print(f"  {decision}: {count}/{len(results['final_candidates'])} ({pct:.1f}%)")

    print()

    # Analysis
    print("=" * 70)
    print("Analysis")
    print("=" * 70)
    print()

    total = len(results['final_candidates'])
    passed = decisions.get('PASS', 0)
    revise = decisions.get('REVISE', 0)
    killed = decisions.get('KILL', 0)

    print(f"✅ PASS: {passed}/{total} ({passed/total*100:.0f}%)")
    print(f"   These antibodies have no major red flags and would proceed to wet-lab")
    print()

    print(f"⚠️  REVISE: {revise}/{total} ({revise/total*100:.0f}%)")
    print(f"   These have addressable issues (e.g., high charge, moderate aggregation)")
    print(f"   Could be improved with sequence optimization")
    print()

    print(f"❌ KILL: {killed}/{total} ({killed/total*100:.0f}%)")
    print(f"   These have fundamental flaws that would require redesign")
    print()

    if passed + revise >= total * 0.8:
        print("🎉 RESULT: Most FDA-approved antibodies pass validation!")
        print("   This suggests the pipeline criteria are appropriately calibrated.")
    else:
        print("⚠️  RESULT: Many FDA-approved antibodies flagged")
        print("   This may indicate overly strict thresholds in some modules.")
        print("   Consider threshold tuning with experimental data.")

    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
