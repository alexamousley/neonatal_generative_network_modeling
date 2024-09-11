**Data key**

Data available at: https://osf.io/ng43c/

Mousley, A., Akarca, D., & Astle, D.E. (2023). Premature birth changes wiring constraints in neonatal structural brain networks. PREPRINT available at Research Square [https://doi.org/10.21203/rs.3.rs-3062369/v1]

For any questions regarding the use of this repository, please get in touch at alexa.mousley@mrc-cbu.cam.ac.uk

**additional_analyses**
- **multiple_density_networks**
  	- Multiple_densities.mat is a struct that contains the thresholds, densities, and consensus networks used in the additional density analysis (Supplementary Fig. 2)
- **multiple_minimum_fiber_lengths**
	- This folder contains 4 versions of binarized networks - each that was tracked using a different minimum fiber length. The number at the end of the title indicates the fiber length used (e.g., binarized_connectomes5.mat was tracked using 5mm minimum fiber length).
- **propensity_matched_analysis**
	- Propensity_unthresholded_connectomes.mat are the unthresholded networks (total streamline count per node).
	- Propensity_binarized_connectomes.mat are the thresholded and binarized networks used in the graph theory analysis.
	- Propensity_generative_models.mat are the GNMs that were used to asses parameter differences between term and preterm groups.
	- **demographics**
		- This folder contains all the same data as the whole-sample demographics folder but here ther is only the data associated with infants used in the propensity matched analysis.
	- **derived_data**
		- Propensity_global_statistics.mat contains the graph theory metrics used in the propensity matched analysis.
		- Top_GNM_parameters.mat contains the top parameters for each infant in the propensity matched analysis.
**atlas**
- Euclidean distance and region labels for the AAL90 neonatal atlas. The atlas is publicly available at: https://www.nitrc.org/projects/pediatricatlas/
- Atlas publication: Shi, F., Yap, P. T., Wu, G., Jia, H., Gilmore, J. H., Lin, W., & Shen, D. (2011). Infant brain atlases from neonates to 1-and 2-year-olds. PloS one, 6(4), e18746.
**demographics**
- This folder contains all demographic data needed to run the generalized additive model statistics.
**derived_data**
- Observed_density.mat is the density of all the thresholded networks used in the manuscript.
- Observed_global_statistics are the global graph theory metrics (Fig. 2).
- Observed_local_statistics are the local graph theory metrics used to asses fit of generative models (Figure 4).
- Rich_club_nodes.mat are the nodes that were identified as rich clubs in the 1.C_identify_rich_club_nodes.m script.
	- **density-controlled_data**
		- This folder contains the same data as the parent folder but has been calculated using the 10% density networks for the density-controlled analysis.
	- **gnm_data**
		- Consensus_Eall_sorted.mat is the top 10 performing models (lowest energies) for the 13 different generative models fit to the consensus network.
		- Consensus_simulated_global_statistics.mat and consensus_simulated_local_statistics.mat contain the graph theory metrics used to assess fit of the simulated consensus networks to the observed consensus networks.
		- Individual_energy_sorted.mat are all the energies for the individually-fit generative models.
		- Individual_mean_top_parameters.mat are the parameters for the best-fitting generative models for each infant.
		- Individual_simulated_global_statistics.mat and individual_simulated_local_statistics.mat contain the graph theory metrics used to assess the fit of simulated networks to observed networks for each infant.
		- Individual_tfdissimilarity.mat contains the TF dissimilarity scores for for each infant.
		- Mean_preterm_counts.mat and mean_preterm_lengths.mat are the averaged connection counts and lengths for the 'preterm' group of generative models (Fig. 6). The same goes for the mean_term_counts.mat and mean_term_lengths.mat for the 'term' group.
		- Total_preterm_lengths.mat and total_term_lengths.mat are the average connection lengths across the whole network (no matter what type of connection) for each group (Fig. 6C).
**observed_networks**
- Binarized_connectomes.mat are the thresholded, binarized networks that have an AVERAGE density of 10%.
- Consensus_network.mat is the consensus network used to fit the original 13 different generative models.
- Density_controlled_binarized_connectomes.mat are the thresholded, binarized networks that are all exactly 10% density (for the density-controlled analysis)
- Density_controlled_consensus_network.mat is the consensus network created from the 10% dense networks.
- Unthresholded_connectomes.mat are the raw, unthresholded networks used in the analysis (total streamlines per node)
**simulated_networks**
NOTE: these are large files that will take a few minutes to download as the generative models have thousands of networks for each observed network.
- Developmental_gnms.mat is a struct that contains the developmental generative models used to represent preterm and term group parameters. Note that these are 'developmental' in that they contain every step of the generative model's growth and thus contain every step of the development of the simulation.
- Individual_generative_models.mat is a struct that contains all the individually-fit generative models (Fig. 4).
- Initial_generative_models.mat is a struct that contains the 13 generative models that were fit to the consensus network to determine which 'value' rule to use in the individually fit networks (Fig. 4A).
