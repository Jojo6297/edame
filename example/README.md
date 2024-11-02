# Example

## Description

The python script: [`EDAME_example.py`](https://github.com/Jojo6297/edame/blob/main/example/EDAME_example.py) depicts the functionality of *EDAME* algorithm on a random signed ($+, -$) undirected Boolean (simulated) network with *n(nodes) = 10*, *n(edges) = 20*. We first use the [ESABO](https://doi.org/10.1371/journal.pcbi.1005361) algorithm to partially predict the network, and then use the iterative technique, *EDAME*, to completely predict the network.

The script outputs a figure: [`Network_Evolution_from_ESABO_to_Original.png`](https://github.com/Jojo6297/edame/blob/main/example/Network_Evolution_from_ESABO_to_Original.png) that depicts the steps of the network deduction process.


![Network_Evolution_from_ESABO_to_Original](https://github.com/Jojo6297/edame/blob/main/example/Network_Evolution_from_ESABO_to_Original.png?raw=true)

### Figure Description
Evolution of partially inferred *ESABO* network, G\*, to a network with the same attractor set as A<sub>0</sub>, which gives us the original network, G<sub>0</sub>. Network parameters are N (nodes) = M<sub>+</sub> (positive edges) = M<sub>-</sub> (negative edges) = 10. Shown below the networks are the attractors of the networks with green: attractors matching to A<sub>0</sub>; black: attractors absent in A<sub>0</sub>; grey: attractors present in A<sub>0</sub> but missing in the corresponding attractor set. The algorithm corrects for false positive and false negative edges shown as dashed edges in the network G\*.