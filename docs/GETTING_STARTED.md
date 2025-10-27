# Getting Started with eTrial

eTrial is a comprehensive in-silico validation platform for therapeutic candidates. This guide will help you get up and running quickly.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Configuration](#configuration)
4. [Running Validations](#running-validations)
5. [Understanding Results](#understanding-results)
6. [Next Steps](#next-steps)

## Installation

### Basic Installation

```bash
# Clone or download the eTrial repository
cd eTrial

# Install core dependencies
pip install -e .

# Or install with specific extras
pip install -e ".[all]"
```

### Optional Dependencies

For specific modules, you may need additional dependencies:

```bash
# Structure modeling (AlphaFold, MD simulations)
pip install -e ".[structure,md]"

# Machine learning models (ESM, NetMHC)
pip install -e ".[ml]"

# PK/PD modeling
pip install -e ".[pkpd]"

# Report generation (PDF support)
pip install -e ".[reporting]"
```

### Cloud GPU Setup (DigitalOcean)

Since you're using DigitalOcean GPUs, ensure you have:

```bash
# CUDA toolkit (for GPU acceleration)
# Check CUDA is available
nvidia-smi

# Install GPU-accelerated PyTorch
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# Optional: Install cupy for GPU-accelerated numpy
pip install cupy-cuda11x
```

## Quick Start

### 1. Define Your Therapeutic Candidate

```python
from etrial import TherapeuticCandidate, Modality

candidate = TherapeuticCandidate(
    name="AB-001",
    modality=Modality.ANTIBODY,
    target="PD-L1",
    sequence="EVQLVESGGGLVQPGGSLRLSCAASGFTFS...",  # Your sequence
)
```

### 2. Run Validation

```python
from etrial import ValidationPipeline

# Initialize pipeline with default configuration
pipeline = ValidationPipeline.from_config("config/default_config.yaml")

# Run validation
result = pipeline.validate(candidate)

# Check decision
print(f"Decision: {result.overall_decision}")
```

### 3. Generate Report

```python
# Generate HTML report
pipeline.generate_dossier(result, format="html")

# Generate PDF report
pipeline.generate_dossier(result, format="pdf")
```

### Command-Line Usage

```bash
# Run validation from command line
python scripts/run_validation.py \
    --name AB-001 \
    --sequence "EVQLVESGGGLVQPGGSLRLSCAASGFTFS..." \
    --target PD-L1 \
    --modality antibody \
    --format html

# Or use an example workflow
python workflows/example_antibody.py
```

## Configuration

### Configuration Files

eTrial uses two main configuration files:

1. **default_config.yaml** - Pipeline settings (which modules to run, parameters)
2. **thresholds.yaml** - Pass/Fail/Revise decision criteria

### Customizing Configuration

```yaml
# config/my_config.yaml

global:
  output_dir: "my_results"
  gpu_enabled: true

modules:
  structure: true
  specificity: true
  developability: true
  immunogenicity: true
  safety: true
  pkpd: false  # Disable PK/PD for faster runs

structure:
  md:
    enabled: true
    n_replicates: 3
    simulation_time_ns: 100
```

Load custom configuration:

```python
pipeline = ValidationPipeline.from_config("config/my_config.yaml")
```

### Adjusting Thresholds

Modify `config/thresholds.yaml` to customize decision criteria:

```yaml
structure_affinity:
  binding_affinity:
    PASS:
      predicted_kd_nm: [null, 30]  # More stringent: KD ≤ 30 nM
    REVISE:
      predicted_kd_nm: [30, 300]
    KILL:
      predicted_kd_nm: [300, null]
```

## Running Validations

### Programmatic API

```python
from etrial import TherapeuticCandidate, ValidationPipeline, Modality
from etrial.structure.affinity import AffinityValidationModule
from etrial.developability.solubility import SolubilityValidationModule

# 1. Create candidate
candidate = TherapeuticCandidate(
    name="Peptide-X",
    modality=Modality.PEPTIDE,
    target="EGFR",
    sequence="MKTAYIAKQR...",
)

# 2. Initialize pipeline
pipeline = ValidationPipeline.from_config("config/default_config.yaml")

# 3. Register specific modules (optional - can auto-register)
pipeline.register_module(AffinityValidationModule)
pipeline.register_module(SolubilityValidationModule)

# 4. Run validation
result = pipeline.validate(candidate)

# 5. Access results
print(f"Decision: {result.overall_decision.value}")

for module_name, module_result in result.module_results.items():
    print(f"{module_name}: {module_result.decision.value}")
    print(f"  {module_result.summary}")
```

### Running Specific Modules

```python
# Run only specific modules
result = pipeline.validate(
    candidate,
    modules=["structure_affinity", "developability"]
)
```

### Batch Processing

```python
# Validate multiple candidates
candidates = [
    TherapeuticCandidate(name="AB-001", ...),
    TherapeuticCandidate(name="AB-002", ...),
    TherapeuticCandidate(name="AB-003", ...),
]

results = {}
for candidate in candidates:
    results[candidate.name] = pipeline.validate(candidate)

# Compare results
for name, result in results.items():
    print(f"{name}: {result.overall_decision.value}")
```

## Understanding Results

### Decision Levels

- **PASS** - Ready for wet-lab validation
- **REVISE** - Addressable issues, recommend optimization
- **KILL** - Fundamental flaws, recommend termination
- **INFORMATIVE** - For non-gating modules

### Accessing Metrics

```python
result = pipeline.validate(candidate)

# Get module results
affinity_result = result.get_module_result("structure_affinity")

# Get specific metrics
kd = affinity_result.get_metric_value("predicted_kd")
print(f"Predicted KD: {kd} nM")

# Check all metrics
for metric in affinity_result.metrics:
    print(f"{metric.name}: {metric.value} {metric.unit or ''}")
    print(f"  Decision: {metric.decision.value}")
```

### Recommendations and Risks

```python
# Get all recommendations
for module_result in result.module_results.values():
    if module_result.recommendations:
        print(f"Recommendations from {module_result.module_name}:")
        for rec in module_result.recommendations:
            print(f"  - {rec}")

# Get all risks
all_risks = []
for module_result in result.module_results.values():
    all_risks.extend(module_result.risks)
```

### Exporting Results

```python
# Export as JSON
result.to_json("results/AB-001.json")

# Generate HTML dossier
pipeline.generate_dossier(result, format="html")

# Generate PDF dossier
pipeline.generate_dossier(result, format="pdf")
```

## Next Steps

### 1. Implement Full Modules

The current implementation includes stub modules. To get production-ready:

```bash
# Implement structure/affinity module
# - Integrate AlphaFold-Multimer
# - Add FoldX, Rosetta scoring
# - Implement OpenMM MD simulations

# Implement specificity module
# - Add BLAST/HMMER search
# - Integrate Foldseek
# - Add pathway analysis

# Implement immunogenicity module
# - Integrate NetMHCIIpan
# - Add population HLA modeling
# - Implement de-immunization suggestions
```

### 2. Customize for Your Use Case

```python
# Create custom validation module
from etrial.core.base import ValidationModule, ValidationResult, Decision

class MyCustomModule(ValidationModule):
    def validate(self, candidate):
        # Your custom validation logic
        metrics = [...]
        decision = Decision.PASS
        summary = "..."

        return self._create_result(
            candidate=candidate,
            decision=decision,
            metrics=metrics,
            summary=summary,
        )

# Register and use
pipeline.register_module(MyCustomModule())
```

### 3. Optimize for Your Hardware

```yaml
# config/gpu_config.yaml

global:
  gpu_enabled: true
  n_jobs: 4  # Number of parallel CPU jobs

structure:
  md:
    gpu_accelerated: true
    device: "cuda:0"

immunogenicity:
  tcell:
    batch_size: 128  # Optimize for your GPU memory
```

### 4. Set Up Continuous Validation

```python
# Integrate with your design pipeline
def design_and_validate(sequence_generator):
    pipeline = ValidationPipeline.from_config("config/default_config.yaml")

    for i, sequence in enumerate(sequence_generator):
        candidate = TherapeuticCandidate(
            name=f"Design-{i:04d}",
            sequence=sequence,
            ...
        )

        result = pipeline.validate(candidate)

        if result.overall_decision == "PASS":
            # Save for experimental testing
            save_for_wet_lab(candidate, result)
```

## Troubleshooting

### Common Issues

**Issue: Module not found**
```bash
# Make sure you're in the eTrial directory
cd /path/to/eTrial

# Install in development mode
pip install -e .
```

**Issue: GPU not detected**
```python
import torch
print(torch.cuda.is_available())  # Should be True

# Check CUDA version matches PyTorch
print(torch.version.cuda)
```

**Issue: Configuration not loading**
```python
from pathlib import Path

# Check paths are correct
config_path = Path("config/default_config.yaml")
print(config_path.exists())  # Should be True
print(config_path.absolute())  # Full path
```

## Support

For questions or issues:

1. Check the documentation in `docs/`
2. Review example workflows in `workflows/`
3. Examine module implementations in `etrial/`
4. Open an issue on GitHub

## What's Next?

Now that you have the framework set up, the next steps are:

1. **Implement production modules** - Replace stub implementations with real tools
2. **Validate the validation** - Run on known therapeutics to calibrate thresholds
3. **Optimize for your hardware** - Tune batch sizes and parallel execution
4. **Build your library** - Screen your designed candidates
5. **Iterate** - Use feedback to improve your design algorithm

Happy validating!
