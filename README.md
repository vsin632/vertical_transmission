# vertical_transmission
Scripts associated with manuscript "Extensive inheritance of gut microbial communities in a superorganismal termite". Authors: Justinn Renelies-Hamilton and Sergio Andreu-Sanchez. Author contact: JR-H (https://orcid.org/0000-0003-3136-6171; claxon71@gmail.com), SA-S (https://orcid.org/0000-0002-3503-9971). Author contributions: JR-H wrote script to analyze 16S amplicon sequencing data (number 1-14b, 16-17), SA-S wrote script for Random Forest model (number 15). Sequencing data available at SRA archive in GenBank (BioProject PRJNA860052). Step 1 refers to transmission to founding reproductive. Step 2 refers to transmission to offspring colonies

# 1. Dada_setup 
from raw files to dada2 objects. Utilizes files in SRA BioProject, Batch_1_new_metadata.csv, Dict_db.fasta, batch_2_metadata.csv, seqtab.nochim.rds, taxa_ditchdb_cp_2.rds

# 2. dada_wrapup
putting dada2 objects together.

# 2b decontamination

# 2c. mock communities

# 2d. rarefaction

# 3. vertical_transmission: 
defining vertically transmitted taxa.

# 4. Fig_1: 
information pertaining to figure 1 including VT strains by lineage etc.

# 5. Venn_figure_2: 
information pertaining to figure 2.

# 6. Correlation_transmission_diversity: 
Correlations between ASV-level diversity, their transmission and their abundance.

# 7. Step_1_emitter_automated: 
What role does abundance play in maternal colony worker strain transmission during step 1? (Fig 3c top panel)

# 8. Step_2_emitter_automated: 
What role does abundance play in maternal colony worker strain transmission during step 2? (Fig 3c bottom panel)

# 9. Abundance_emitter_stats: 
Stats on previous two scripts.

# 10. Step_1_receiver_automateds_stats: 
What role does transmission play in alate strain abundance? (Fig. 4a top panel)

# 11. Step_2_receiver_automated_stats: 
What role does transmission play in offspring colony caste strain abundance? (Fig 4a bottom panel)

# 12. ALDEx2_model: 
Aldex2 models for step 1 and 2.

# 13. Simulations_Abu_dri_null_model.

# 14a. Network_calculation: 
Calculate networks

# 14b. Network_analysis:
Analyze networks.

# 15. RF_Analysis: 
Random Forest Model- utilizes dataframe df_serio.rds and df_sergio_20210210.rds

# 16. RF_Interpretation: 
Value transformation and SHAP value plots.

# 17. RF_Plots: 
AUC Plots.
