#!/usr/bin/env bash
# Reruns the analysis from the retrieved extracts (source_data/, ldref/, kp_*.csv) and regenerates the result tables and figures.
# Requirements: Python 3.11 (numpy, pandas, scipy, matplotlib, python-docx, openpyxl, Pillow), R 4.3 with MendelianRandomization 0.10.0,
# MRPRESSO 1.0, coloc 6.0.3, susieR 0.12.35; PLINK 1.9 on PATH as plink1.9. Expected run time about 2 hours (MR-PRESSO seeds and the B = 10^6 convergence runs dominate).
# MVMR reference functions (ref_MVMR/, MVMR GitHub commit 8be0d94) are sourced; no internet access is needed.
# GTEx, Europe PMC and the eight FinnGen R12 records were retrieved in a browser on 2026-09-29 and are archived in source_data/.
set -euo pipefail
python3 step6_unit_check.py
python3 step6_mmHg_calibration.py
python3 step6_v4_forward.py
python3 step6_phenome_screen.py
Rscript mr_v4.R main
for c in S_mvpEU D_mvpEU H_mvpEU; do PRESSO_N=2400 Rscript mr_v4.R presso "$c"; done
Rscript mr_v4.R main
python3 step6_v4_loo.py
python3 step6_mvmr_prep.py
Rscript mvmr_v4.R
python3 step6_coloc_v4_prep.py
Rscript coloc_v4.R
# ---- v5 verification analyses
mkdir -p presso_stab
for c in S_mvpEU D_mvpEU H_mvpEU; do for s in 314159 20261 20262; do Rscript presso_stability.R "$c" 2400 "$s" > "presso_stab/log_${c}_${s}" 2>&1; done; done
python3 step7_presso_stability.py
python3 step7_plink_clump.py
python3 step7_se_validation.py
Rscript mvmr_v5.R
Rscript mvmr_condF_cov.R
Rscript susie_L20_tert.R
python3 step7_finngen_audit.py
python3 step7_instrument_flow.py
python3 step7_output_consistency.py
python3 step7_dataset_register.py
python3 step7_phenome_dictionary.py
python3 -c "import numpy as np,pandas as pd;R=pd.read_csv('mr_v4_inputs/forward_r2_0.1_correlation.csv',index_col=0).values;e=np.linalg.eigvalsh(R);o=np.abs(R[np.triu_indices_from(R,1)]);pd.DataFrame([dict(n_variants=R.shape[0],min_eigenvalue=e.min(),max_eigenvalue=e.max(),condition_number=e.max()/e.min(),max_abs_offdiag_r=o.max(),n_pairs_abs_r_ge_0_1=int((o>=0.1).sum()),positive_definite=bool(e.min()>0))]).to_csv('step7_gls_conditioning.csv',index=False)"
# ---- v6 verification analyses (V5-01 to V5-30)
python3 step8_se_formula_tests.py
python3 step8_plink_both.py
python3 step8_reverse_source.py
Rscript mr_v6.R
python3 step8_presso_convergence.py
python3 step8_combine_results.py
Rscript mvmr_v6.R
python3 step8_ld_dosage.py
Rscript coloc_v6.R
python3 step8_finngen_af_fin.py
python3 step8_dataset_register.py
python3 step8_aux_outputs.py
# ---- v7 (response to the v6 evidence audit)
python3 step8_plink_both.py            # allele-aware LD record lookup; PLINK --r2 dprime; selection-chain mechanism
python3 step8_reverse_source.py; Rscript mr_v6.R; python3 step8_combine_results.py
Rscript mvmr_v7_qhet.R                 # per-exposure qhet export with bootstrap (bootstrap B = 1000)
python3 step8_aux_outputs.py
python3 step9_reverse_inputs.py        # S13 + reverse reconstruction test (fails the build on mismatch)
python3 step9_fallback_ledger.py       # S38
python3 step9_output_consistency.py    # S46 (fails the build on mismatch)
python3 step9_mvmr_primary.py          # S21
for c in S_mvpEU D_mvpEU H_mvpEU; do Rscript presso_equivalence.R "$c" 3000 20260930; done
python3 step9_presso_equivalence.py    # S70 (fails the build on mismatch)
# ---- v8 (response to the v7 review)
python3 step10_environment.py         # software versions at build time
python3 step8_dataset_register.py      # Table 1 / S1 role wording
python3 step10_ld_lookup_trail.py      # S74 LD-panel lookup correction trail
Rscript step10_qhet_objective_test.R   # S75 pointwise qhet objective test
# ---- v9 (response to the pre-submission review)
python3 step11_environment.py           # software versions (mr.raps from GitHub dd79b5b with the rsnps plotting import removed)
python3 step11_nakao_weights.py         # S81 Nakao et al. European weights for the 116 instruments (from the stored portal query)
Rscript step11_robust_mr.R              # S80 contamination mixture, cML-MA-BIC-DP, MR-RAPS; Nakao-weight IVW/WM/Egger
python3 step11_review_checks.py         # S76 Steiger-removed reverse variants, S77 rs2736100 SE, S78 exact-value checks
python3 step11_hypertension_se_pattern.py   # S79 MVP hypertension SE pattern and frequency-model SE sensitivity
Rscript coloc_v9_conditional.R          # S82 conditional colocalization (about 30 min)
# ---- v10 (response to the review of v9)
python3 step12_environment.py
Rscript coloc_v10_htn_modelSE.R        # S84 MVP hypertension coloc with allele-frequency model SEs
Rscript step12_conmix_check.R          # S85 independent contamination-mixture likelihood and phi
python3 -c "import pandas as pd; pd.read_csv('step12_conmix_check.csv').merge(pd.read_csv('step12_conmix_phi.csv'), on='code').to_csv('step12_conmix_check_S85.csv', index=False)"
python3 step12_window_check.py         # S83 window check
# ---- v12 (response to the v11 correction sheet)
python3 step14_dbgap_match.py           # S86 MVP portal records matched to dbGaP phs002453 analyses (analysis list retrieved 3 Oct 2026: step14_dbgap_raw.csv)
Rscript step14_conmix_profile.R          # S88 contamination-mixture profile versus quadratic
Rscript coloc_v12_htn_conditional.R    # S87 conditional colocalization with MVP hypertension (case fraction from S86)
# ---- figures (main-figure text checked at 180 mm print width; figures/proof_report.csv)
rm -f figures/proof_report.csv
python3 fig1_workflow_v6.py; python3 fig2_mr_v7.py; python3 fig3_locus_v7.py; python3 fig4_ancestry_v6.py
python3 figS1_mr_diag_v6.py; python3 figS2_phenome_coloc_v5.py; python3 figS3_mr_grid_v6.py; python3 figS4_coloc_grid_v5.py; python3 figS5_alignment_v5.py
python3 figS7_susie_diag_v6.py; python3 figS6_verification_v6.py
# ---- locus-specific Wald ratios (Supplemental Table S3), added after the analysis release
(cd mr_v4_inputs && python3 ../../locus_wald_ratios/locus_heterogeneity_test.py)
sha256sum step8_mr_results_v6.csv step8_mvmr_v6_estimates.csv step8_mvmr_v6_conditionalF.csv step9_mvmr_qhet.csv step8_presso_convergence_summary.csv step8_plink_both_directions_summary.csv step8_reverse_source_ledger.csv step8_susie_stability_summary.csv step8_gls_ld_sensitivity.csv step8_finngen_harmonisation_audit_v6.csv step6_coloc_v4_abf.csv step9_reverse_inputs_S13.csv step9_presso_equivalence.csv step10_ld_lookup_correction_S74.csv step11_robust_mr.csv step11_steiger_reverse_removed.csv step11_conditional_coloc.csv step11_hypertension_se_ivw.csv step12_htn_coloc_model_se.csv step12_conmix_check_S85.csv step12_window_check.csv step14_dbgap_mvp_match.csv step14_htn_conditional_coloc.csv step14_conmix_profile.csv > run_all_checksums.txt
echo "Done. Compare run_all_checksums.txt with ../expected_checksums_v15.txt (MR-PRESSO package rows at B = 2,400 can differ if the random stream differs)."
