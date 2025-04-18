# Visualization of EDAME Network Results

## Description
This directory contains the filtered networks generated by running EDAME on microbial abundance data from the [`datasets`](https://github.com/Jojo6297/edame/tree/main/datasets) folder, categorized by the three health states (healthy, CRC, IBD).

- **[1_EDAME_Heatmaps.py](https://github.com/Jojo6297/edame/blob/main/applications/combined/1_EDAME_Heatmaps.py)**: This script generates a heatmap for each `pos_edge/neg_edge` parameter combination (shown here: `pos_edge` = 5, `neg_edge` = 2) using the filtered networks as input (see [`README.md`](https://github.com/Jojo6297/edame/blob/main/applications/README.md) for parameter details). Each heatmap displays network similarity across all combinations, measured by the Jaccard Index, and is saved in the [`Heatmaps`](https://github.com/Jojo6297/edame/tree/main/applications/combined/Heatmaps) folder.

- **[1_EDAME_Networks.py](https://github.com/Jojo6297/edame/blob/main/applications/combined/1_EDAME_Networks.py)**: This script produces cumulative networks for each health state and parameter combination. The resulting cumulative networks are stored in the [`CumulativeNetworks_AllStudies`](https://github.com/Jojo6297/edame/tree/main/applications/combined/CumulativeNetworks_AllStudies) directory.

Since the EDAME algorithm has stochastic elements with some randomized steps, the output networks can vary a little bit.