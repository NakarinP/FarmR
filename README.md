# FarmR
Estimating farm-to-farm disease transmission using pathogen genetic and animal movement data

# Input data; example file
1. Time resolved phylogenetic tree (NEXUS format from BEAST MCC tree) of pathogen's genetic sequences; beastmcc.tree
2. Metadata of the sequence samples including at least ID of isolates, sampling date, and farm premises; metadata.csv
3. Animal movement data comprising origin (tail), destination (head), and date of shipment; movement.csv

# Optional data for ERGMs analysis
1. Metadata spreadsheet may include farm production type, herd size, or farm's geographical coordinates
2. All other farms' coordinates in the area where the samples were collected may be used for farm's spatial point density calculation; locationall.csv

# Getting started
1. Open "Restimation-ergm.R"
2. Install all required packages (first 10 lines)
3. At line 12 of the script: change "x" to "path to your working directory" that keeps the input files
4. Run them all!

# Primary outputs
1. Transphylo's transmission tree: MCMC diagnostics (trace & ESS) and colored transmission tree
2. Infection chain length (ICL; animal level) vs Movement path-length (MPL; farm-level): Normality test, Pearson's correlation, and linear regression plot
3. Farm-level reproduction number (R) estimation: Boxplot of ICL ~ MPL with the cut-off value at MPL = 1 and statistic summary of the farm-level R

# Farm-to-farm transmission network (Additional outputs for ERGMs analysis)
1. Dyads: dyadic relationship between all samples including transmission link (1,0), movement path-length (steps), sampling interval (days), and distance between farms (km)
2. Node attributes: same as the metadata file plus sampling season, farm's spatial point density, and in-degree and out-degree of the 6-month period (according to sampling date) movement network 
