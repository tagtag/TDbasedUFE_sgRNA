# TDbasedUFE_sgRNA
This is a sample R code to perform the analyses in "Gene and cell line efficiency of CRISPR computed by tensor decomposition in genome-wide CRISPR-Cas9 knockout screens". https://doi.org/10.1101/2025.06.12.659265

The following files must be downlaoded

* msb145216-sup-0001-datasets1.xlsx

* msb145216-sup-0004-datasets4.xlsx

https://doi.org/10.15252/msb.20145216

* yusa_raw_v10.tab
* tko_counts.txt
* Wang2015_2017_merged_counts.txt
* Achilles_raw_GeckoV2.tab
* Avana_sgrna_raw_readcounts_matched.csv
* Avana_sgrnamapping.csv



from Data.zip
https://www.doi.org/10.6084/m9.figshare.6002438

#------

In addition to above, two additional R files are added to address reviwers' comments.

jacks_auc.R

To comute AUC for JACKS using the standard procedure

It requires JACKS.ZIP from https://www.doi.org/10.6084/m9.figshare.6002438, essential.txt, nonessential.txt.


make_model_ids_from_jacks_avana_v2.R

depmap_chronos_mean_auc_fixed.R

To compute AUC for DepMap for Avana cell line.

It requires CRISPRGeneEffect.csv and Model,csv from https://depmap.org/portal/data_page/?tab=allData (25Q3), essential.txt, nonessential.txt.
