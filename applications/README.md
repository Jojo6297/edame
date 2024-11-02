
# Applications

## Description
This folder contains data and scripts required to generate networks for each health state (healthy, IBD, CRC) using the data available in the [`datasets`](https://github.com/Jojo6297/edame/tree/main/datasets) folder.

## Instructions for Running the Scripts
Each health state folder contains cohort-specific datasets and a Python script (e.g., [`0_healthy.py`](https://github.com/Jojo6297/edame/blob/main/applications/healthy/0_healthy.py)) to execute the EDAME algorithm on these datasets. To run a script, use the following command in the command prompt:

```bash
python <Script.py> <edge_info> <pos_edge> <neg_edge> <input_abundance> <fail_limit>
```

**Parameters:**
- **`Script.py`** – Can be any of the three: [`0_healthy.py`](https://github.com/Jojo6297/edame/blob/main/applications/healthy/0_healthy.py), [`0_CRC.py`](https://github.com/Jojo6297/edame/blob/main/applications/CRC/0_CRC.py), or [`0_IBD.py`](https://github.com/Jojo6297/edame/blob/main/applications/IBD/0_IBD.py).
- **`edge_info`** – Set to either `"yesedge"` or `"noedge"`, depending on whether you want to initialize EDAME with a predefined number of edges from the ESABO network.
- **`pos_edge`** – If `"yesedge"` is selected, this takes an integer (e.g., 4, 5, or 6, as described in the paper [citation pending]). If `"noedge"` is selected, set this to `"noedge"`.
- **`neg_edge`** – If `"yesedge"` is selected, this takes an integer (e.g., 2 as used in the paper [citation pending]). If `"noedge"` is selected, set this to `"noedge"`.
- **`input_abundance`** – Specifies the abundance dataset, with rows representing *phyla* and columns representing *samples*.
- **`fail_limit`** – An integer specifying the maximum number of attempts for the EDAME algorithm to find a network that matches the input data’s attractors. In the paper [citation pending], we set this to 50.

Additionally, each health state folder includes a pre-configured bash script (e.g., [`0_healthy_bash.sh`](https://github.com/Jojo6297/edame/blob/main/applications/healthy/0_healthy_bash.sh)) to run EDAME on all datasets within that cohort. The bash script already includes the command arguments outlined above for convenient execution.


#### Output
Each run of an `input_abundance` dataset with specified `pos_edge/neg_edge` values produces the following outputs:
- **Final Networks**: All networks (*`___FinalNetworkIntMat.csv`*) predicted with the highest *Jaccard Index* between their attractors and the original attractors.
- **Iteration Results**: The Jaccard Index of attractors obtained in each iteration (*`___Iterations.csv`*).

Once networks are generated for all conditions, you can filter for the most frequently occurring network solutions by running the script [`FilterResults.py`](https://github.com/Jojo6297/edame/blob/main/applications/FilterResults.py) with a single argument:
- `<health state>` – Specify as `healthy`, `CRC`, or `IBD`. <br><br>


The filtered networks are transferred to the [`combined`](https://github.com/Jojo6297/edame/tree/main/applications/combined) folder, which contains scripts for visualizing the results. For detailed instructions, please refer to the [`README.md`](https://github.com/Jojo6297/edame/blob/main/applications/combined/README.md) file.