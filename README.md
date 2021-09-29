# FarmR
Estimating farm-to-farm disease transmission using pathogen genetic and animal movement data

# Prerequisites; example file
1. Time resolved phylogenetic tree (NEXUS format from BEAST MCC tree) of pathogen's genetic sequences; beastmcc.tree
2. Metadata of the sequence samples including at least ID of isolates, sampling date, and farm premises; metadata.csv
3. Animal movement data comprising origin (tail), destination (head), and date of shipment; movement.csv
4. Optional data for ERGMs analysis:  
   4.1 Metadata spreadsheet may include farm production type, herd size, or farm's geographical coordinates
   4.2 All other farms' coordinates in the area where the samples were collected may be used for farm's spatial density calculation; locationall.csv
