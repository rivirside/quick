#!/usr/bin/env python3
from etrial import TherapeuticCandidate, Modality, ValidationPipeline
from etrial.structure.affinity import AffinityValidationModule
from etrial.developability.solubility import SolubilityValidationModule
from pathlib import Path

# Test Adalimumab
candidate = TherapeuticCandidate(
    name='Adalimumab',
    modality=Modality.ANTIBODY,
    target='TNF-alpha',
    heavy_chain='EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS',
    light_chain='DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQRYNRAPYTFGQGTKVEIK'
)

pipeline = ValidationPipeline.from_config(Path('config/default_config.yaml'))
pipeline.register_module(AffinityValidationModule)
pipeline.register_module(SolubilityValidationModule)

result = pipeline.validate(candidate)

# Debug: manually aggregate decisions
module_decisions = {name: res.decision.value for name, res in result.module_results.items()}
print(f'Module decisions dict: {module_decisions}')
print(f'Any KILL? {"KILL" in module_decisions.values()}')
print(f'Any REVISE? {"REVISE" in module_decisions.values()}')
print(f'Count REVISE: {list(module_decisions.values()).count("REVISE")}')

print(f'\nAdalimumab Overall: {result.overall_decision.value}')
print(f'Total modules: {len(result.module_results)}')
print(f'\nAll module decisions:')
for name, mod_result in result.module_results.items():
    print(f'  {name}: {mod_result.decision.value}')
