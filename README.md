# Genetically predicted telomere length and blood pressure: Mendelian randomization in the Million Veteran Program and colocalization at *TERC*

Analysis code, derived summary-level data and result tables for the manuscript of the same title (Ikram MF, Jaleel A, Aziz U, Uppal SS, Jaleel S, Zia A, Farrukh S, Baig S, Arora I, Yaqinuddin A, Ganguly P, Allsopp RC), submitted to *Circulation: Genomic and Precision Medicine*.

## Contents

| Folder / file | What it is |
|---|---|
| `analysis/` | All analysis scripts (Python and R), the retrieved summary-statistic extracts (`source_data/`, `kp_*.csv`, `ldref/`), harmonized inputs (`mr_v4_inputs/`) and result files |
| `analysis/run_all.sh` | Reruns the analyses and regenerates the figures from the stored extracts |
| `locus_wald_ratios/` | Locus-specific Wald ratios at *TERC* and *TERT* (Supplemental Table S3): `locus_wald_ratios.py` (ratios and sampling-error tests) and `locus_heterogeneity_test.py` (region-excluded comparators and heterogeneity-adjusted tests); `locus_wald_ratios.csv` holds the final values |
| `release_artifacts/` | MR-PRESSO shared simulation arrays, qhet bootstrap replicates and PLINK clumping outputs needed to reproduce the verification checks |
| `figures/` | Main Figures 1–4 as submitted |
| `supplemental_tables/` | Supplemental Tables S6–S88 workbook (with companion sheets for Tables S1 and S3) |
| `environment.txt` | Software versions of the analysis release |
| `release_manifest.csv`, `expected_checksums_v15.txt`, `build_checks.csv` | SHA-256 manifest of inputs, scripts and results, expected result checksums and the automated release checks (release ID v15-dddb0ef8a7b4) |

## Requirements

Python 3.11 (numpy, pandas, scipy, matplotlib, openpyxl), R 4.3 with MendelianRandomization 0.10.0, MRPRESSO 1.0, coloc 6.0.3 and susieR 0.12.35, mr.raps (GitHub qingyuanzhao/mr.raps, commit dd79b5b), MVMR reference functions (included in `analysis/ref_MVMR/`, GitHub WSpiller/MVMR commit 8be0d94) and PLINK 1.9 on the path as `plink1.9`. See `environment.txt`.

## Reproducing the results

```bash
cd analysis
cp -r ../release_artifacts/* .
bash run_all.sh            # about 2 hours; no internet access needed
sha256sum -c ../expected_checksums_v15.txt
```

MR-PRESSO rows from the R package at 2,400 simulations can differ if the R random stream differs; all other checksummed outputs should match.

## Data sources and terms

All inputs are public summary statistics, retrieved on 29 September 2026:

- UK Biobank leukocyte telomere length GWAS (Codd et al., *Nat Genet* 2021), figshare release — extract of genome-wide significant variants and the *TERC*/*TERT* regions only.
- Common Metabolic Diseases Knowledge Portal (https://hugeamp.org): Million Veteran Program (Verma et al., *Science* 2024; dbGaP phs002453.v1.p1), Biobank Japan, Genes & Health and other datasets, as β, P and N per variant.
- FinnGen release R12, endpoint I9_HYPTENS (extract of the analyzed variants).
- 1000 Genomes phase 3 (LD reference, regions and instruments only), GTEx v8 (eQTL look-ups).

Third-party extracts are redistributed here for reproducibility under each provider's terms; please cite the original data sources (see the manuscript reference list) when reusing them.

## Licence

Code: MIT licence (`LICENSE`). Derived data and result tables: CC BY 4.0. Third-party summary statistics remain under their providers' terms.

## Citation

Please cite the article (details to be added on publication) and this repository (`CITATION.cff`).
