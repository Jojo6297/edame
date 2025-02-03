# Example

## Enhancement of initial networks by EDAME
Here we illustrate with an example how sensitivity of attractors with respect to small changes in the network can be utilized to predict a network from its partially predicted network version (see Methods for the underlying algorithm). In an iterative application of the algorithm EDAME, we evolve a partially predicted network to a network that has the same set of attractors as those of the original network.

We simulate a random undirected network G0 with N nodes, M+ positive edges, M− negative edges, then find attractors of the network, A0, using synchronous update rules given in Eq (1). The initial network is inferred using the ESABO ([19,20]) technique. ESABO (Entropy Shifts of Abundance Vectors under Boolean Operations) is a network inference method that uses the co-occurrence of two nodes in the attractors to
calculate entropy shift scores assigned to each possible edge in the network. Edges of the network are deduced based on these scores. Edges with a score above a certain positive threshold are positive edges and below a negative threshold are considered negative edges in the inferred network. The network inferred by ESABO is denoted by G∗, with attractor set A∗. In most cases, ESABO is only able to predict the network partially, hence G∗ is not equal to G0. We apply EDAME on G∗ to reach G0, as shown in Fig 4. Even though the initial network overlap at 0th iteration is above 80%, attractor overlap is at 50%.


## Description

The python script: [`EDAME_example.py`](https://github.com/Jojo6297/edame/blob/main/example/EDAME_example.py) depicts the functionality of *EDAME* algorithm on a random signed ($+, -$) undirected Boolean (simulated) network with *n(nodes) = 10*, *n(edges) = 20*. We first use the [ESABO](https://doi.org/10.1371/journal.pcbi.1005361) algorithm to partially predict the network, and then use the iterative technique, *EDAME*, to completely predict the network.

The script outputs a figure: [`Network_Evolution_from_ESABO_to_Original.png`](https://github.com/Jojo6297/edame/blob/main/example/Network_Evolution_from_ESABO_to_Original.png) that depicts the steps of the network deduction process.


![Network_Evolution_from_ESABO_to_Original](https://github.com/Jojo6297/edame/blob/main/example/Network_Evolution_from_ESABO_to_Original.png?raw=true)

### Figure Description
Evolution of partially inferred *ESABO* network, G\*, to a network with the same attractor set as A<sub>0</sub>, which gives us the original network, G<sub>0</sub>. Network parameters are N (nodes) = M<sub>+</sub> (positive edges) = M<sub>-</sub> (negative edges) = 10. Shown below the networks are the attractors of the networks with green: attractors matching to A<sub>0</sub>; black: attractors absent in A<sub>0</sub>; grey: attractors present in A<sub>0</sub> but missing in the corresponding attractor set. The algorithm corrects for false positive and false negative edges shown as dashed edges in the network G\*.