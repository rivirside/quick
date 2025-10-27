#!/usr/bin/env python3
"""
Command-line interface for running eTrial validations.

Usage:
    python run_validation.py --name AB-001 --sequence EVQLV... --target PD-L1 --modality antibody
    python run_validation.py --config my_candidate.yaml
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from etrial import TherapeuticCandidate, ValidationPipeline, Config, Modality
from etrial.structure.affinity import AffinityValidationModule
from etrial.developability.solubility import SolubilityValidationModule
from loguru import logger


def main():
    parser = argparse.ArgumentParser(
        description="Run eTrial validation on a therapeutic candidate",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Candidate definition
    parser.add_argument(
        "--name",
        type=str,
        help="Candidate name/identifier",
    )
    parser.add_argument(
        "--sequence",
        type=str,
        help="Amino acid sequence (for biologics/peptides) or SMILES (for small molecules)",
    )
    parser.add_argument(
        "--target",
        type=str,
        help="Target protein name",
    )
    parser.add_argument(
        "--modality",
        type=str,
        choices=[m.value for m in Modality],
        help="Therapeutic modality",
    )
    parser.add_argument(
        "--structure",
        type=Path,
        help="Path to structure file (PDB)",
    )

    # Configuration
    parser.add_argument(
        "--config",
        type=Path,
        default="config/default_config.yaml",
        help="Path to configuration YAML file (default: config/default_config.yaml)",
    )
    parser.add_argument(
        "--thresholds",
        type=Path,
        default="config/thresholds.yaml",
        help="Path to thresholds YAML file (default: config/thresholds.yaml)",
    )

    # Output
    parser.add_argument(
        "--output",
        type=Path,
        help="Output directory for results",
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["html", "pdf", "json"],
        default="html",
        help="Report format (default: html)",
    )

    # Execution options
    parser.add_argument(
        "--modules",
        type=str,
        nargs="+",
        help="Specific modules to run (runs all enabled if not specified)",
    )
    parser.add_argument(
        "--no-gating",
        action="store_true",
        help="Disable decision gating (run all modules even if KILL)",
    )

    args = parser.parse_args()

    # ==================================================================
    # Load Configuration
    # ==================================================================

    logger.info("=" * 60)
    logger.info("eTrial - In-Silico Clinical Trial Platform")
    logger.info("=" * 60)

    if args.config.exists():
        logger.info(f"Loading configuration from {args.config}")
        config = Config.from_yaml(args.config)
    else:
        logger.warning(f"Config file not found: {args.config}. Using defaults.")
        config = Config.default()

    # Load thresholds
    if args.thresholds.exists():
        logger.info(f"Loading thresholds from {args.thresholds}")
        config.load_thresholds(args.thresholds)
    else:
        logger.warning(f"Thresholds file not found: {args.thresholds}")

    # Set output directory
    if args.output:
        config.output_dir = args.output

    # ==================================================================
    # Create Candidate
    # ==================================================================

    logger.info("Creating therapeutic candidate...")

    candidate = TherapeuticCandidate(
        name=args.name,
        sequence=args.sequence,
        modality=Modality(args.modality),
        target=args.target,
        structure_file=args.structure,
    )

    logger.info(f"Candidate: {candidate.name}")
    logger.info(f"  Modality: {candidate.modality.value}")
    logger.info(f"  Target: {candidate.target}")
    logger.info(f"  Hash: {candidate.get_hash()}")

    # ==================================================================
    # Initialize Pipeline
    # ==================================================================

    logger.info("Initializing validation pipeline...")

    pipeline = ValidationPipeline(
        config=config,
        enable_gating=not args.no_gating,
    )

    # Register modules
    # In production, this would auto-discover or be configured
    logger.info("Registering validation modules...")

    pipeline.register_module(AffinityValidationModule(config.structure))
    pipeline.register_module(SolubilityValidationModule(config.developability))

    # TODO: Register other modules as they're implemented
    # pipeline.register_module(SpecificityModule(config.specificity))
    # pipeline.register_module(ImmunogenicityModule(config.immunogenicity))
    # pipeline.register_module(SafetyModule(config.safety))
    # pipeline.register_module(PKPDModule(config.pkpd))

    logger.info(f"Registered {len(pipeline.modules)} modules")

    # ==================================================================
    # Run Validation
    # ==================================================================

    logger.info("")
    logger.info("Starting validation...")
    logger.info("-" * 60)

    try:
        result = pipeline.validate(
            candidate=candidate,
            modules=args.modules,
        )

        logger.info("-" * 60)
        logger.info(f"Validation completed!")
        logger.info(f"Overall Decision: {result.overall_decision.value}")
        logger.info(f"Runtime: {result.total_runtime_seconds:.2f}s")
        logger.info(f"Modules run: {len(result.module_results)}")

        # ==================================================================
        # Generate Report
        # ==================================================================

        logger.info("")
        logger.info("Generating validation dossier...")

        report_path = pipeline.generate_dossier(
            result,
            format=args.format,
        )

        logger.info(f"Report saved: {report_path}")

        # ==================================================================
        # Print Summary
        # ==================================================================

        logger.info("")
        logger.info("=" * 60)
        logger.info("VALIDATION SUMMARY")
        logger.info("=" * 60)

        from etrial.core.reporting import ReportGenerator
        reporter = ReportGenerator()
        summary = reporter.generate_executive_summary(result)
        print(summary)

        logger.info("")
        logger.info(f"Full report: {report_path}")

        # Exit code based on decision
        if result.overall_decision == "PASS":
            sys.exit(0)
        elif result.overall_decision == "REVISE":
            sys.exit(1)
        else:  # KILL
            sys.exit(2)

    except Exception as e:
        logger.error(f"Validation failed: {e}")
        raise


if __name__ == "__main__":
    main()
