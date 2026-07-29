# figshare deposit — Data S1 (MAP-D supplementary data workbook)

## Suggested title
**Data S1. MAP-D Atlas: integrated proteomic associations with cardiometabolic hallmarks of
type 2 diabetes, triangulated with Mendelian randomization and STEP 1/2 semaglutide trials**

## Description (abstract / body)

This dataset is the supplementary data workbook (Data S1) for the **Metabolic Atlas of the
Proteome in Diabetes (MAP-D)**. MAP-D maps the circulating proteome across the progression from
normoglycemia to type 2 diabetes (T2D) and triangulates each protein across three independent
lines of evidence — cross-sectional observational association, genetic causal inference
(bidirectional Mendelian randomization), and interventional response to the GLP-1 receptor
agonist semaglutide (STEP 1/2 trials) — before following persistent proteins to incident disease
risk and approved-drug targets.

Associations were computed for **2,923 circulating proteins** (Olink Explore 3072) in UK Biobank
participants (analytic samples of ~30,500–47,963 per protein) against three cardiometabolic
hallmarks — **adiposity (BMI), glycemia (HbA1c), and insulin resistance (triglyceride-to-HDL
cholesterol ratio, TG/HDL)** — within three mutually exclusive glycemic strata: **normoglycemia,
prediabetes, and incident T2D**. Genetic evidence comes from a split-sample bidirectional MR
design (UK Biobank cis-pQTL instruments and held-out trait GWAS; deCODE SomaScan pQTLs for the
reverse direction). Interventional evidence comes from published semaglutide-induced protein
changes in the STEP 1 (obesity without T2D) and STEP 2 (T2D) trials (SomaScan). Incident coronary
artery disease (CAD) associations and DrugBank target matching complete the resource.

The workbook is a single Microsoft Excel file (`.xlsx`) containing **8 worksheets**:

1. **Fwd_Rev_Integrated_Evidence** (8,769 rows = 2,923 proteins × 3 hallmarks) — the master
   integrated-evidence table. For each protein × hallmark it reports bidirectional observational
   associations (`Fwd_Obs_*` = protein-first, protein → hallmark; `Rev_Obs_*` = hallmark-first,
   hallmark → protein, each with beta, SE, P), bidirectional Mendelian randomization
   (`Fwd_MR_*` protein-first and `Rev_MR_*` hallmark-first, with beta, SE, P, Bonferroni flag and
   significance indicator), an MR concordance flag, semaglutide trial effects
   (`Estimate_STEP1`/`Pval_STEP1`, `Estimate_STEP2`/`Pval_STEP2`, `STEP_Mean`), colocalization
   posterior probability (`Coloc_PP_H4`, `Coloc_Sig`), and a summary `Causal_Class`
   (Neither / Forward_Only / Reverse_Only / Bidirectional_Concordant / Bidirectional_Discordant).
   `Phenotype` ∈ {BMI, HbA1c, TG/HDL-C Ratio}.

2. **STEP1_2_VALIDATION** (5,217 rows) — stage-stratified comparison of UK Biobank reverse
   observational associations (hallmark → protein) against semaglutide-induced protein changes in
   STEP 1 (obesity, non-T2D) and STEP 2 (T2D). Columns: `Phenotype`, `Subgroup`
   (Normoglycemic / Prediabetes / T2D), `Protein_UKB`, `EntrezGeneID`, `Estimate_UKB`,
   `Significance_UKB`, the matched STEP protein, and `Estimate_STEP1/2` with significance. This
   sheet spans a broader panel of traits (BMI, HbA1c, HDL, LDL, TG/HDL, systolic and diastolic
   blood pressure) than the three focal hallmarks.

3. **Evidence_Funnel_Table** (5 rows) — the number of proteins retained at each successive
   evidence layer: (1) UK Biobank observational (Bonferroni-significant); (2) + nominal STEP trial
   change (P < 0.05); (3) + nominal reverse MR (IVW P < 0.05); (4) + reverse-MR sign concordance;
   (5) CAD DrugBank hits — reported for BMI and HbA1c across the three glycemic stages.

4. **INCIDENT CAD REGRESSION OUTPUT** (1,008 rows) — per-protein logistic-regression associations
   with incident coronary artery disease in the UK Biobank Olink subsample, for the 1,008
   trial-triangulated proteins. Columns include `estimate` (log odds ratio), `std.error`,
   `p.value`, Benjamini–Hochberg `FDR`, model `AUC`, and `SampleSize`; 367 proteins reach
   FDR < 0.05.

5. **DRUGBANK_HITS** (79 rows) — approved drug–protein matches for CAD-associated proteins with a
   directionally appropriate action (inhibitors/antagonists for higher-risk proteins, OR > 1;
   agonists/activators for protective proteins, OR < 1). Columns include `gene_symbol`, `OR`,
   `drug_name`, `drugbank_id`, `actions`, and `action_needed`; 79 drug–protein pairs spanning
   62 approved drugs and 27 proteins.

6. **MAP-D_GTEX** (2,667 rows) — tissue-specificity annotation (GTEx) of hallmark-associated
   proteins by glycemic stage: for each protein × hallmark × subgroup, the GTEx tissue(s)
   (`Term`) in which the encoding gene is expression-enriched (≥ 4-fold above the cross-tissue
   mean), alongside the UK Biobank effect estimate.

7. **LDSC_genetic_corrs** (3 rows) — pairwise LD-score-regression genetic correlations (rg) among
   the three hallmarks, estimated in the held-out UK Biobank non-proteomics GWAS, with SE, Z, P,
   95% CI, SNP count, and per-trait observed-scale SNP heritabilities.

8. **Analysis_Layer_Model_Specs** (5 rows) — a methods key describing each analysis layer
   (reverse observational PWAS; STEP 1 and STEP 2 trial triangulation; split-sample forward MR;
   two-sample reverse MR via deCODE): direction, inclusion subsample, model/estimand, outcome,
   exposure, covariates, approximate N, and instruments/notes.

**Usage notes.** Observational models were adjusted for age, sex, assessment centre, household
income, fasting time, and 40 genetic principal components; the observational Bonferroni threshold
is p < 1.9 × 10⁻⁶ (0.05 / 26,307 tests). Cross-platform (Olink ↔ SomaScan) comparisons are
limited to the 2,142 proteins with gene-symbol matches. Individual-level UK Biobank data are not
included and are available only via application to UK Biobank.

**Related resources.**
- Manuscript: Tangirala, Isaac, Gehad, … Tierney & Patel (Nature Metabolism).
- Interactive atlas (web application): https://stejat98-map-d.share.connect.posit.cloud/
- Analysis code and app source: https://github.com/stejat98/MAP-D
- Hallmark GWAS summary statistics (BMI, HbA1c, TG/HDL): https://doi.org/10.6084/m9.figshare.32048550
- Code archive (Zenodo): https://doi.org/10.5281/zenodo.17071087

## Keywords
proteomics; type 2 diabetes; UK Biobank; Olink; Mendelian randomization; GLP-1 receptor agonist;
semaglutide; STEP trial; proteome-wide association study; cardiometabolic; coronary artery disease;
drug repurposing; DrugBank

## Suggested citation
Tangirala S, Isaac S, Gehad Y, et al. Data S1. MAP-D Atlas: integrated proteomic associations with
cardiometabolic hallmarks of type 2 diabetes, triangulated with Mendelian randomization and
STEP 1/2 semaglutide trials. figshare. Dataset. https://doi.org/10.6084/m9.figshare.30007306

## License
CC BY 4.0 (recommended for open data reuse).
