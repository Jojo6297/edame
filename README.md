# EDAME: Edge Detection via Attractor Mismatch Evaluation

## Description
The EDAME algorithm detects discrepancies in Boolean attractors (stable states) between two networks, using these variations to identify edge differences. This approach enables the prediction of microbial interaction networks based on their Boolean stable states.

EDAME builds on previous research [1], which employed Boolean network modeling to systematically analyze interactions among low-abundance species in the human gut microbiome.

### Repository Structure
- **[`example`](https://github.com/Jojo6297/edame/tree/main/example)**: Contains a Python script to run EDAME on a simulated Boolean network, offering a demonstration of the algorithm’s functionality.
- **[`datasets`](https://github.com/Jojo6297/edame/tree/main/datasets)**: Holds datasets for comparative analysis across three health states: Healthy, IBD (Inflammatory Bowel Disease), and CRC (Colorectal Cancer). These datasets were generated with the R package **curatedMetagenomicData** [2] using the provided R script.
- **[`applications`](https://github.com/Jojo6297/edame/tree/main/applications/)**: Contains scripts and instructions to apply EDAME on datasets from [2], along with primary results published in the paper *"Evaluating Changes in Attractor Sets Under Small Network Perturbations to Infer Reliable Microbial Interaction Networks from Abundance Patterns"* (citation pending).

## Requirements
Scripts are compatible with Python >= 3.8.
Package versions are listed in [`requirements.txt`](https://github.com/Jojo6297/edame/blob/main/requirements.txt).

### References
1. Claussen JC, Skiecevičienė J, Wang J, Rausch P, Karlsen TH, et al. (2017). Boolean analysis reveals systematic interactions among low-abundance species in the human gut microbiome. *PLOS Computational Biology*, 13(6): e1005361. https://doi.org/10.1371/journal.pcbi.1005361

2. Pasolli E, Schiffer L, Manghi P, et al. (2017). Accessible, curated metagenomic data through ExperimentHub. *Nat Methods*, 14(11):1023-1024. https://doi.org/10.1038/nmeth.4468
