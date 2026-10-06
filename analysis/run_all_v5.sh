#!/usr/bin/env bash
# Reruns the v5 analysis from the retrieved extracts (source_data/, ldref/, kp_*.csv) and regenerates every table, figure and document.
# Requirements: Python 3.11 (numpy, pandas, scipy, matplotlib, python-docx, openpyxl, Pillow), R 4.3 with MendelianRandomization 0.10.0,
# MRPRESSO 1.0, coloc 6.0.3, susieR 0.12.35; PLINK 1.9 on PATH as plink1.9. Expected run time about 75 minutes (MR-PRESSO seeds dominate).
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
# ---- figures and documents
python3 fig1_workflow_v5.py; python3 fig2_mr_v5.py; python3 fig3_locus_v5.py; python3 fig4_ancestry_v5.py
python3 figS1_mr_diag_v5.py; python3 figS2_phenome_coloc_v5.py; python3 figS3_mr_grid_v5.py; python3 figS4_coloc_grid_v5.py; python3 figS5_alignment_v5.py; python3 figS6_verification_v5.py
python3 build_supplement_v5.py; python3 build_manuscript_v5.py; python3 build_strobe_v5.py; python3 build_response_v5.py
sha256sum step6v4_mr_results_R.csv step7_mvmr_results.csv step6_coloc_v4_abf.csv step7_plink_clump_ivw.csv step7_se_validation.csv step7_susie_L20_TERT_coloc.csv > run_all_checksums.txt
echo "Done. Compare run_all_checksums.txt with expected_checksums_v5.txt (MR-PRESSO rows can differ if the random stream differs)."
