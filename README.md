# Network inference algorithms applied to M1 model

## Description
We evaluated the Bayesian stochastic block modeling method (using the graph-tool package) on the same motor cortex model described above (approx 0.5M synapses). The complex local microcircuit connectivity depended on cell class (3 excitatory and 2 inhibitory) and cortical depth (with bins of 100 or 150 μm for excitatory cells). See bioRxiv paper for details (Figure 2): https://www.biorxiv.org/content/early/2018/03/27/201707.full.pdf+html

The algorithm inferred neuronal groupings and connection strengths between them that closely matched those used in the original rules. 
This was done based solely on individual synapses without prior knowledge of cell populations or locations. 
The algorithm also grouped cells at other scales, providing interesting insights into network structure, e.g., clearly separating superficial from deeper layers (see Results below). 

Execution time on a single core was 890 seconds, suggesting it is well-suited for larger-scale problems, given the algorithm implementation is parallelizable. 

## Results

- The 32 inferred neuronal groupings and its composition in terms of the M1 network populations. Grouping matchly close original populations and subdivisions based on cortical depth:
![m1_infer_level_0_groupPiecharts](https://github.com/salvadord/network_inference/blob/master/data/m1_infer_level_0_groupPiecharts.png)

- The 32 inferred groupings represented by colors using the original cell locations (note cell locations were not known by the inference algorithm):
![m1_infer_state_sfdp](https://github.com/salvadord/network_inference/blob/master/data/m1_infer_state_sfdp.png)

- Connectivity matrix of the 32 inferred groupings. Matrix closely matches original connectivity rules:
![m1_infer_level_0_mat](https://github.com/salvadord/network_inference/blob/master/data/m1_infer_level_0_mat.png)

- Nested Stochastic Block Model of the M1 model. The 32 groupings are labeled showing the corresponding M1 population that composes it. The nested model shows interesting insights into the network structure at multiple scales. For example, L5 populations (IT5A, IT5B, PT5B) are subdivided into subpopulations, similar to the original connectivity cortical-depth based rules. Also, at the larger scale, superficial and deeper populations are clearly grouped together:
![m1_infer_state_labeled](https://github.com/salvadord/network_inference/blob/master/data/m1_infer_state_labeled.png)
