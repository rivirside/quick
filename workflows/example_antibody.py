#!/usr/bin/env python3
"""
Example workflow: Antibody validation

This demonstrates how to programmatically run eTrial validation
on an antibody therapeutic.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from etrial import TherapeuticCandidate, ValidationPipeline, Config, Modality
from etrial.structure.affinity import AffinityValidationModule
from etrial.developability.solubility import SolubilityValidationModule


def main():
    print("=" * 70)
    print("eTrial Example: Antibody Validation")
    print("=" * 70)
    print()

    # ==================================================================
    # Step 1: Define your therapeutic candidate
    # ==================================================================

    print("Step 1: Defining therapeutic candidate...")

    candidate = TherapeuticCandidate(
        name="Anti-PD-L1-001",
        modality=Modality.ANTIBODY,
        target="PD-L1",
        sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        metadata={
            "indication": "NSCLC",
            "design_notes": "Humanized IgG1, targeting PD-L1 checkpoint",
        }
    )

    print(f"  Name: {candidate.name}")
    print(f"  Target: {candidate.target}")
    print(f"  Modality: {candidate.modality.value}")
    print()

    # ==================================================================
    # Step 2: Load configuration
    # ==================================================================

    print("Step 2: Loading configuration...")

    config_path = Path("config/default_config.yaml")
    thresholds_path = Path("config/thresholds.yaml")

    if config_path.exists():
        config = Config.from_yaml(config_path)
        if thresholds_path.exists():
            config.load_thresholds(thresholds_path)
    else:
        print("  Using default configuration (config files not found)")
        config = Config.default()

    print(f"  Config hash: {config.get_hash()}")
    print()

    # ==================================================================
    # Step 3: Initialize pipeline and register modules
    # ==================================================================

    print("Step 3: Initializing validation pipeline...")

    pipeline = ValidationPipeline(config=config)

    # Register validation modules
    pipeline.register_module(AffinityValidationModule(config.structure))
    pipeline.register_module(SolubilityValidationModule(config.developability))

    # In production, you would register all modules:
    # pipeline.register_module(SpecificityModule(config.specificity))
    # pipeline.register_module(ImmunogenicityModule(config.immunogenicity))
    # pipeline.register_module(SafetyModule(config.safety))
    # pipeline.register_module(PKPDModule(config.pkpd))

    print(f"  Registered {len(pipeline.modules)} modules: {pipeline.list_modules()}")
    print()

    # ==================================================================
    # Step 4: Run validation
    # ==================================================================

    print("Step 4: Running validation...")
    print("-" * 70)

    result = pipeline.validate(candidate)

    print("-" * 70)
    print()

    # ==================================================================
    # Step 5: Review results
    # ==================================================================

    print("Step 5: Validation Results")
    print()

    print(f"Overall Decision: {result.overall_decision.value}")
    print(f"Total Runtime: {result.total_runtime_seconds:.2f}s")
    print()

    print("Module Results:")
    for module_name, module_result in result.module_results.items():
        print(f"  - {module_name}: {module_result.decision.value}")
        print(f"    {module_result.summary}")

        if module_result.recommendations:
            print(f"    Recommendations:")
            for rec in module_result.recommendations:
                print(f"      • {rec}")

        if module_result.risks:
            print(f"    Risks:")
            for risk in module_result.risks:
                print(f"      ⚠ {risk}")
        print()

    # ==================================================================
    # Step 6: Generate report
    # ==================================================================

    print("Step 6: Generating validation dossier...")

    report_path = pipeline.generate_dossier(result, format="html")

    print(f"  Report saved: {report_path}")
    print()

    # ==================================================================
    # Step 7: Access individual metrics
    # ==================================================================

    print("Step 7: Accessing specific metrics...")
    print()

    # Get affinity metrics
    affinity_result = result.get_module_result("structure_affinity")
    if affinity_result:
        kd = affinity_result.get_metric_value("predicted_kd")
        print(f"  Predicted KD: {kd} nM")

        buried_sasa = affinity_result.get_metric_value("buried_sasa")
        print(f"  Buried SASA: {buried_sasa} Ų")

    # Get developability metrics
    dev_result = result.get_module_result("developability")
    if dev_result:
        camsol = dev_result.get_metric_value("camsol_score")
        print(f"  CamSol Score: {camsol}")

        n_glyc = dev_result.get_metric_value("n_glycosylation_sites")
        print(f"  N-glycosylation sites: {n_glyc}")

    print()

    # ==================================================================
    # Step 8: Export results
    # ==================================================================

    print("Step 8: Exporting results...")

    # Save JSON
    json_path = config.output_dir / "results" / f"{candidate.name}_results.json"
    json_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_json(json_path)

    print(f"  JSON results: {json_path}")
    print()

    print("=" * 70)
    print("Validation Complete!")
    print("=" * 70)

    # Return exit code based on decision
    if result.overall_decision == "PASS":
        return 0
    elif result.overall_decision == "REVISE":
        return 1
    else:
        return 2


if __name__ == "__main__":
    sys.exit(main())
