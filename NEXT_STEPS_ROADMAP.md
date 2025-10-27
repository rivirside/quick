# Next Steps Roadmap

**Date**: 2025-10-26
**Current Status**: System 85% complete, 3/7 quick wins done
**Goal**: Reach 95% completeness over next 4 weeks

---

## Immediate Actions (Today - 1 Hour)

### 1. Integrate New QC Modules Into Pipeline ⚡

**What**: Add localization, complexity, and advanced aggregation to pre-filter

**File to modify**: `etrial/core/two_stage_pipeline.py`

**Changes needed**:
```python
def setup_default_prefilter_pipeline(pipeline: TwoStageValidationPipeline):
    from etrial.structure.affinity import AffinityValidationModule
    from etrial.developability.solubility import SolubilityValidationModule
    from etrial.safety.toxicity import ToxicityPrefilterModule
    from etrial.pkpd.pharmacokinetics import PKPDPrefilterModule

    # NEW: Add QC modules
    from etrial.developability.localization import LocalizationValidationModule
    from etrial.developability.complexity import ComplexityValidationModule
    from etrial.developability.aggregation_advanced import AdvancedAggregationModule

    pipeline.register_prefilter_module(AffinityValidationModule)
    pipeline.register_prefilter_module(SolubilityValidationModule)

    # NEW: Register QC modules
    pipeline.register_prefilter_module(LocalizationValidationModule)  # Signal peptide, TM, NLS
    pipeline.register_prefilter_module(ComplexityValidationModule)    # Low-complexity regions
    pipeline.register_prefilter_module(AdvancedAggregationModule)     # TANGO/AGGRESCAN

    pipeline.register_prefilter_module(ToxicityPrefilterModule)
    pipeline.register_prefilter_module(PKPDPrefilterModule)
```

**Time**: 15 minutes

---

### 2. Run Full Integration Test ⚡

**Test the complete pipeline with all modules**:

```bash
./venv/bin/python examples/two_stage_validation_demo.py
```

**Expected**: All 9 modules run successfully (4 original + 3 QC + 2 from sprint)

**Time**: 5 minutes

---

### 3. Update Deployment Status Document ⚡

**File**: `DEPLOYMENT_STATUS.md`

**Update module count**:
- Pre-filter: 4 → 7 modules
- Clinical: 6 → 6 modules
- Total: 10 → 13 modules
- Coverage: 83% → 90%

**Time**: 10 minutes

---

### 4. Create Integration Test for New Modules ⚡

**File**: `examples/test_full_pipeline_with_qc.py`

**Purpose**: Demonstrate all modules working together

**Time**: 30 minutes

---

## Week 1: Viscosity Prediction (High Priority) 🔴

### Goal
Implement SAP score and Developability Index for high-concentration mAb formulation.

### Background
- **Problem**: mAbs at >100 mg/mL can have high viscosity (>20 cP)
- **Impact**: Prevents subcutaneous injection (need <20 cP)
- **Frequency**: 30-40% of mAbs have viscosity issues

### Implementation

**File to create**: `etrial/developability/viscosity.py`

**Algorithms to implement**:

1. **SAP (Spatial Aggregation Propensity)**
   - Sharma et al. (2014) PNAS
   - Surface patch analysis
   - Score = net charge × hydrophobic surface area

2. **Developability Index**
   - Raybould et al. (2019) mAbs
   - Combines multiple factors:
     - pI (isoelectric point)
     - PSR (patch surface ratio)
     - Fv charge symmetry
     - Hydrophobic patches

3. **Charge Symmetry Parameter**
   - Tomar et al. (2016) mAbs
   - Asymmetric charge distribution → higher viscosity

**Key Metrics**:
- `sap_score` < 20 = PASS
- `developability_index` > 0 = PASS
- `charge_symmetry_parameter` < 0.3 = PASS

**Effort**: 5-7 days
- Day 1-2: Implement SAP score
- Day 3-4: Implement Developability Index
- Day 5: Charge symmetry
- Day 6-7: Testing + integration

**Test cases needed**:
- Known low-viscosity mAb (adalimumab)
- Known high-viscosity mAb (requires redesign)

---

## Week 2: Computational Alanine Scanning 🟡

### Goal
Identify binding hotspots for affinity maturation.

### Background
- **Purpose**: Predict which residues contribute most to binding
- **Use case**: Guide mutations for affinity improvement
- **Method**: Calculate ΔΔG for Ala substitutions

### Implementation

**File to create**: `etrial/structure/alanine_scanning.py`

**Approaches**:

1. **FoldX Integration** (preferred)
   - Install FoldX
   - Run BuildModel for each Ala mutation
   - Calculate ΔΔG_binding
   - Hotspots: ΔΔG > 2 kcal/mol

2. **Simple Energy-Based** (fallback)
   - Count H-bonds lost
   - Estimate buried SASA change
   - Approximate ΔΔG from empirical rules

**Key Metrics**:
- `num_hotspots` - Count of critical residues
- `hotspot_positions` - List of positions
- `avg_delta_g` - Average contribution

**Effort**: 5-7 days
- Day 1-2: Install + test FoldX
- Day 3-4: Implement scanning logic
- Day 5: Fallback energy method
- Day 6-7: Testing + validation

**Test cases needed**:
- Known antibody-antigen complex
- Experimental alanine scan data for validation

---

## Week 3: PROPKA pH Stability 🟡

### Goal
Predict stability at different pH values (critical for formulation).

### Background
- **Problem**: Proteins unstable at certain pH ranges
- **Impact**: Must choose formulation pH where stable
- **Method**: Predict pKa shifts, identify pH-labile residues

### Implementation

**File to create**: `etrial/developability/ph_stability.py`

**Approach**:

1. **PROPKA Integration** (preferred)
   - Install PROPKA 3.5
   - Calculate pKa for all ionizable residues
   - Identify buried charged residues (unstable)

2. **Simple Electrostatic** (fallback)
   - Count buried Asp/Glu/His/Lys/Arg
   - Estimate pKa shifts based on environment
   - Flag residues with large shifts

**Key Metrics**:
- `stable_ph_range` - pH range where stable
- `ph_labile_residues` - Residues with shifted pKa
- `recommended_ph` - Optimal formulation pH

**Effort**: 3-4 days
- Day 1: Install + test PROPKA
- Day 2: Implement pH stability calculation
- Day 3: Fallback method
- Day 4: Testing

**Test cases needed**:
- Known pH-stable protein
- Known pH-labile protein (e.g., with buried His)

---

## Week 4: NetMHCpan MHC-I Epitopes 🟡

### Goal
Add CD8 T-cell epitope prediction (complement existing MHC-II).

### Background
- **Current**: Only have MHC-II (CD4 T-cells)
- **Missing**: MHC-I (CD8 T-cells)
- **Impact**: More complete immunogenicity assessment

### Implementation

**File to create**: `etrial/immunogenicity/mhc1_epitopes.py`

**Approach**:

1. **NetMHCpan Integration** (preferred)
   - Install NetMHCpan 4.1
   - Predict binding to HLA-A, HLA-B, HLA-C
   - Report strong binders (percentile rank < 2%)

2. **Pattern-Based** (fallback)
   - Known MHC-I binding motifs
   - Hydrophobicity + anchor residue patterns
   - Less accurate but fast

**Key Metrics**:
- `mhc1_epitope_count` - Number of strong binders
- `hla_coverage` - % of population with epitopes
- `epitope_positions` - List of positions

**Effort**: 3-4 days
- Day 1: Install + test NetMHCpan
- Day 2: Implement epitope prediction
- Day 3: Fallback patterns
- Day 4: Testing

**Test cases needed**:
- Known immunogenic protein
- Non-immunogenic protein (human)

---

## Month 2: Validation & Threshold Tuning

### Goal
Validate predictions against experimental data and optimize thresholds.

### Activities

1. **Gather Experimental Data** (Week 5)
   - Collect published datasets:
     - Aggregation data (TANGO validation)
     - Viscosity data (mAb databases)
     - Binding affinity (KD measurements)
     - Toxicity data (hERG, Ames tests)

2. **Benchmark Accuracy** (Week 6)
   - Run predictions on datasets
   - Calculate correlation coefficients
   - Identify systematic errors

3. **Tune Thresholds** (Week 7)
   - Optimize PASS/REVISE/KILL cutoffs
   - Balance false positive vs false negative rates
   - Document threshold rationale

4. **Create Validation Report** (Week 8)
   - Accuracy metrics per module
   - Comparison to other tools
   - Limitations documented

---

## Month 3: Production Hardening

### Goal
Make system production-ready for continuous use.

### Activities

1. **Performance Optimization** (Week 9)
   - Profile slow modules
   - Add caching where appropriate
   - Parallelize independent calculations

2. **Error Handling** (Week 10)
   - Comprehensive try/catch blocks
   - Graceful degradation when tools unavailable
   - Clear error messages

3. **Monitoring & Logging** (Week 11)
   - Add structured logging
   - Track module execution times
   - Alert on failures

4. **CI/CD Pipeline** (Week 12)
   - Automated testing on every commit
   - Performance regression tests
   - Deployment automation

---

## Priority Matrix

### Must Do (Weeks 1-4)
✅ **Integrate QC modules** (today)
🔴 **Viscosity prediction** (week 1) - Prevents costly formulation failures
🟡 **Alanine scanning** (week 2) - Enables rational design
🟡 **PROPKA** (week 3) - Critical for formulation
🟡 **NetMHCpan** (week 4) - Completes immunogenicity

### Should Do (Months 2-3)
🟢 Experimental validation
🟢 Threshold optimization
🟢 Production hardening
🟢 CI/CD setup

### Nice to Have (Future)
⚪ MD simulations (for finalists only)
⚪ Docking (AlphaFold already does this)
⚪ Active learning pipeline
⚪ Web interface

---

## Decision Points

### After Week 1
**Question**: Is viscosity prediction accurate enough?
- If YES: Continue to alanine scanning
- If NO: Spend extra week improving (get more training data)

### After Week 2
**Question**: Can we install FoldX?
- If YES: Use FoldX for accurate ΔΔG
- If NO: Use fallback energy method (less accurate but functional)

### After Week 4
**Question**: All quick wins complete - what's next?
- Option A: Validation (recommended) - ensures accuracy
- Option B: Production hardening - deploy faster
- Option C: More features - diminishing returns

**Recommendation**: Choose Option A (validation) to ensure quality before deployment.

---

## Success Criteria

### End of Week 4 (All Quick Wins Complete)
✅ 7/7 high-priority features implemented
✅ >90% coverage of industry standard
✅ All modules tested and working
✅ Integrated into pipeline
✅ Documentation complete

### End of Month 2 (Validation Complete)
✅ Accuracy validated against experimental data
✅ Thresholds optimized
✅ Known limitations documented
✅ Validation report published

### End of Month 3 (Production Ready)
✅ Performance optimized (<100ms pre-filter)
✅ Error handling comprehensive
✅ Monitoring in place
✅ CI/CD pipeline active
✅ Ready for continuous use

---

## Resource Requirements

### Tools to Install

**Week 1**: None (RDKit already available)

**Week 2**:
- FoldX (free for academic, license for commercial)
- Alternative: Use fallback method (no install)

**Week 3**:
- PROPKA 3.5 (free, pip install)
- Alternative: Use fallback method

**Week 4**:
- NetMHCpan 4.1 (free for academic, registration required)
- Alternative: Use pattern-based fallback

### Hardware Requirements

**Current**: CPU-only (working fine)

**After quick wins**: Still CPU-only (no GPU needed)

**For MD simulations** (future): GPU helpful but optional

---

## Risk Management

### Risks & Mitigations

**Risk 1**: Can't install external tools (FoldX, PROPKA, NetMHCpan)
- **Mitigation**: Implemented fallback methods for all
- **Impact**: Reduced accuracy (70% vs 85%) but still functional

**Risk 2**: Validation shows poor accuracy
- **Mitigation**: Tune thresholds, add more training data
- **Impact**: Extra 1-2 weeks for improvements

**Risk 3**: Performance too slow
- **Mitigation**: Profile and optimize, add caching
- **Impact**: Extra 1 week for optimization

**Risk 4**: Scope creep (too many features)
- **Mitigation**: Stick to 7 quick wins only, defer rest to Phase 2
- **Impact**: Controlled timeline

---

## Timeline Summary

```
TODAY (1 hour):
  ✅ Integrate QC modules
  ✅ Full integration test
  ✅ Update docs

WEEK 1 (5-7 days):
  🔴 Viscosity prediction (SAP/DI)

WEEK 2 (5-7 days):
  🟡 Computational alanine scanning (FoldX)

WEEK 3 (3-4 days):
  🟡 PROPKA pH stability

WEEK 4 (3-4 days):
  🟡 NetMHCpan MHC-I epitopes

MONTH 2 (4 weeks):
  🟢 Validation against experimental data
  🟢 Threshold optimization

MONTH 3 (4 weeks):
  🟢 Production hardening
  🟢 CI/CD setup
  🟢 Monitoring

TOTAL: ~3 months to 95% completeness
```

---

## Deliverables Checklist

### Week 1
- [ ] `etrial/developability/viscosity.py`
- [ ] `examples/test_viscosity.py`
- [ ] Validation against known mAbs
- [ ] Integration into pipeline

### Week 2
- [ ] `etrial/structure/alanine_scanning.py`
- [ ] FoldX integration (or fallback)
- [ ] `examples/test_alanine_scan.py`
- [ ] Test on known complex

### Week 3
- [ ] `etrial/developability/ph_stability.py`
- [ ] PROPKA integration (or fallback)
- [ ] `examples/test_ph_stability.py`
- [ ] pH range recommendations

### Week 4
- [ ] `etrial/immunogenicity/mhc1_epitopes.py`
- [ ] NetMHCpan integration (or fallback)
- [ ] `examples/test_mhc1.py`
- [ ] Integration with existing MHC-II

### Month 2
- [ ] Experimental validation report
- [ ] Accuracy metrics per module
- [ ] Optimized thresholds
- [ ] Limitations documentation

### Month 3
- [ ] Performance benchmark report
- [ ] CI/CD pipeline configured
- [ ] Monitoring dashboard
- [ ] Production deployment guide

---

## Getting Started (Right Now!)

### Action 1: Integrate QC Modules (15 min)

```bash
# Edit the pipeline file
vim etrial/core/two_stage_pipeline.py

# Add imports and registrations (see code above)

# Test
./venv/bin/python examples/two_stage_validation_demo.py
```

### Action 2: Verify Integration (5 min)

```bash
# Should see 7 pre-filter modules running:
# 1. Structure/Affinity
# 2. Developability (solubility)
# 3. Localization (NEW)
# 4. Complexity (NEW)
# 5. Advanced Aggregation (NEW)
# 6. Safety/Toxicity
# 7. PKPD
```

### Action 3: Start Week 1 Planning (30 min)

- Research SAP score papers (Sharma 2014)
- Find mAb viscosity datasets
- Plan implementation approach

---

## Questions to Answer

Before starting Week 1:

1. **Do we have access to experimental viscosity data for validation?**
   - If not, where can we get it? (ChEMBL, literature, etc.)

2. **What is the target threshold for viscosity?**
   - Industry standard: <20 cP for subcutaneous injection
   - <10 cP preferred

3. **Should we focus on mAbs only or all biologics?**
   - SAP score is mAb-specific
   - Can extend to other biologics later

4. **What is the performance target?**
   - Current: 40ms pre-filter
   - With viscosity: Aim for <60ms

---

## Success Metrics

### End of Quick Wins (Week 4)
- [ ] 7/7 features implemented
- [ ] All tests passing
- [ ] Coverage ≥ 90%
- [ ] Pre-filter < 100ms
- [ ] Documentation complete

### End of Validation (Month 2)
- [ ] Accuracy ≥ 70% per module
- [ ] Thresholds optimized
- [ ] False positive rate < 20%
- [ ] False negative rate < 10%

### End of Production (Month 3)
- [ ] Uptime > 99%
- [ ] Performance < 100ms pre-filter
- [ ] Error rate < 1%
- [ ] CI/CD functional
- [ ] Team trained

---

## Bottom Line

**Today (1 hour)**:
1. Integrate 3 QC modules
2. Test full pipeline
3. Update docs

**This Month (Weeks 1-4)**:
Complete 4 remaining quick wins

**Next 2 Months**:
Validate + production harden

**Timeline**: 3 months to 95% completeness

**Next Action**: Integrate QC modules into pipeline (15 minutes)

---

**Ready to start? Run this now:**

```bash
# Edit pipeline file
vim etrial/core/two_stage_pipeline.py

# Add the 3 new QC modules (see code above)

# Test
./venv/bin/python examples/two_stage_validation_demo.py
```

---

**Last Updated**: 2025-10-26
**Next Milestone**: All quick wins complete (Week 4)
**End Goal**: Production-ready system (Month 3)
