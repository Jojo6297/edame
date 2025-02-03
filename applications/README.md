# Applications

## Description
This folder contains data and scripts required to generate networks for each health state (healthy, T2D, CRC, IBD) using the data available in the [`datasets`](https://github.com/Jojo6297/edame/tree/main/datasets) folder.

## Instructions for Running the Scripts
Each health state folder contains cohort-specific datasets and a Python script (e.g., [`0_healthy.py`](https://github.com/Jojo6297/edame/blob/main/applications/healthy/0_healthy.py)) to execute the EDAME algorithm on these datasets. To run a script, use the following command in the command prompt:

```bash
python <Script.py> <edge_info> <pos_edge> <neg_edge> <input_abundance> <fail_limit>
```

**Parameters:**
- **`Script.py`** – Can be any of the four: [`0_healthy.py`](https://github.com/Jojo6297/edame/blob/main/applications/healthy/0_healthy.py), [`0_T2D.py`](https://github.com/Jojo6297/edame/blob/main/applications/T2D/0_T2D.py), [`0_CRC.py`](https://github.com/Jojo6297/edame/blob/main/applications/CRC/0_CRC.py), or [`0_IBD.py`](https://github.com/Jojo6297/edame/blob/main/applications/IBD/0_IBD.py).
- **`edge_info`** – We first infer a partial network from ESABO and then run EDAME on it. While using ESABO, one can decide a threshold (e.g. selecting top 6 edges with the highest positive zscores) for edges that should come out of ESABO inference. Set this argument to either `"yesedge"` or `"noedge"`, depending on whether you want to initialize EDAME with a predefined number of edges from the ESABO network or not.
- **`pos_edge`** – If `"yesedge"` is selected, this takes an integer (e.g., 4, 5, or 6, as described in the paper [citation pending]). If `"noedge"` is selected, set this to `"noedge"`.
- **`neg_edge`** – If `"yesedge"` is selected, this takes an integer (e.g., 2 as used in the paper [citation pending]). If `"noedge"` is selected, set this to `"noedge"`.
- **`input_abundance`** – Specifies the abundance dataset, with rows representing *phyla* and columns representing *samples*. (e.g. [`1_LifeLinesDeep_2016_healthy.relative.abundance.csv`](https://github.com/Jojo6297/edame/blob/main/applications/healthy/1_LifeLinesDeep_2016_healthy.relative.abundance.csv))
- **`fail_limit`** – An integer specifying the maximum number of attempts for the EDAME algorithm to find a network that matches the input data’s attractors. In the paper [citation pending], we set this to 50.

Additionally, each health state folder includes a pre-configured bash script (e.g., [`0_healthy_bash.sh`](https://github.com/Jojo6297/edame/blob/main/applications/healthy/0_healthy_bash.sh)) to run EDAME on all datasets within that cohort. The bash script already includes the command arguments outlined above for convenient execution.


**Example command in folder "healthy":**
```bash
python 0_healthy.py yesedge 6 2 1_LifeLinesDeep_2016_healthy.relative.abundance.csv 50
```

#### Output
Each run of an `input_abundance` dataset with specified `pos_edge/neg_edge` values produces the following outputs:
- **Final Networks**: All networks (*`*_FinalNetworkIntMat.csv`*) predicted with the highest *Jaccard Index* between their attractors and the original attractors.
- **Iteration Results**: The Jaccard Index of attractors obtained in each iteration (*`*_Iterations.csv`*).
- **Final Attractors**: Attractors (*`*_FinalAttractors.csv`*) of the final network.

Once networks are generated for all conditions, you can filter for the most frequently occurring network solutions by running the script [`FilterResults.py`](https://github.com/Jojo6297/edame/blob/main/applications/FilterResults.py) with a single argument:
- `<health state>` – Specify as `healthy`, `T2D`, `CRC`, or `IBD`. <br><br>

**Example command in folder "applications":**
```bash
python FilterResults.py healthy
```

The filtered networks are transferred to the [`combined`](https://github.com/Jojo6297/edame/tree/main/applications/combined) folder, which contains scripts for visualizing the results. For detailed instructions, please refer to the [`README.md`](https://github.com/Jojo6297/edame/blob/main/applications/combined/README.md) file.