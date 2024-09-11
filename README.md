**Code repository for:**
Mousley, A., Akarca, D., & Astle, D.E. (2023). Premature birth changes wiring constraints in neonatal structural brain networks. PREPRINT available at Research Square [https://doi.org/10.21203/rs.3.rs-3062369/v1]

For any questions regarding the use of this repository, please get in touch at alexa.mousley@mrc-cbu.cam.ac.uk

**Requirements**
Software
1) MATLAB 2020B (installation: https://uk.mathworks.com/help/install/install-products.html)
2) RStudio 4.1.2 (installation: https://rstudio.com/products/rstudio/download/)

Toolboxes/Functions
1) Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/home?authuser=0)
Rubinov, M., & Sporns, O. (2010). Complex network measures of brain connectivity: uses and 	interpretations. Neuroimage, 52(3), 1059-1069.
2) Consensus network creation (https://www.brainnetworkslab.com/coderesources)
Betzel, R. F., Griffa, A., Hagmann, P., & Mišić, B. (2019). 
Distance-dependent consensus thresholds for generating group-representative
structural brain networks. Network neuroscience, 3(2), 475-496.

**Data availability**
All derived, anonymized data used in this publication is available at: https://osf.io/ng43c/. With this data downloaded, all published code can be run (please see "data_key" for further description). The reconstructed images used are also freely available at: https://brain.labsolver.org/hcp_d2.html. The original data is available at: https://www.developingconnectome.org/data-release/third-data-release/.

**Script outline:**

set_paths.m - This script is where you set all the paths to the data and will be called for all MATLAB scripts.

1) Observed topological analysis
   A network_thresholding.m
   B create_consensus_network.m
   C identify_rich_club_nodes.m
   D calculate_organizational_measures.m
   E organizational_analysis.R  

3) Initial generative network model selection 
	A run_initial_generative_models.m
	B analyze_initial_generative_models.m  
	C compare_initial_models_plot.R

4) Individually fit generative network models 
	A run_individual_generative_models.m
	B organize_individual_model_output.m 
	C analyze_individual_generative_models.m 
	D visualize_individual_generative_models.R 

5) Group-representative developmental models 
	A run_developmental_generative_models.m
	B analyze_developmental_generative_models.m 
	C visualize_developmental_generative_models.R 

Propensity matched analysis - This code is slightly adapted versions of the code above to be run with the propensity matched sample
	A calculate_organizational_properties.m
	B organizational_analysis

