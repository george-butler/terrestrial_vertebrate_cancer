# No evidence for Peto’s Paradox in terrestrial vertebrates

### George Butler<sup>1,2</sup>, Joanna Baker<sup>3</sup>, Sarah R. Amend<sup>2</sup>, Kenneth J. Pienta<sup>2</sup>, Chris Venditti<sup>3</sup>

<sup><sup>1</sup>University College London Cancer Institute, University College London, London, UK</sup>

<sup><sup>2</sup>Cancer Ecology Center, The Brady Urological Institute, Johns Hopkins School of Medicine, Baltimore, USA</sup>

<sup><sup>3</sup>School of Biological Sciences, University of Reading, Reading, UK</sup>

<p align="center">
  <img width="380" height="355" src="/example_image/exceptional_birds_and_mammals.png">
</p>

Our paper is available [here](https://www.pnas.org/doi/10.1073/pnas.2422861122)

All code is written in R with the majority of the scripts utilizing the [MCMCglmm](https://cran.r-project.org/web/packages/MCMCglmm/index.html) library to fit Multivariate Phylogenetic Generalized Linear Mixed Models (MPGLMMs). 

--------------------------------------------------------------------------------------------------------------------------------------------------
The necessary code to reproduce the analysis in "No evidence for Peto's Paradox in terrestrial vertebrates" can be found in the directories: "body_size_only", "body_size_and_path_length", and "exceptional_species". 

To reproduce the analysis please first clone the entire repository and then extract the contents of the repository to the home directory. The scripts contained within the 3 directories mentioned above assume that the "data" directory containing the species trait data, "species_data.csv", and the phylogenetic trees, "terrestrial_vertebrate_tree.nexus.trees", is present in the home directory and that you are using a UNIX machine. If not, please change the lines within each script pointing to the location of the species trait data and the phylogenetic tree. Similarly, the scripts assume that you have a directory entitled "results" in the home directory. If not, please make a directory entitled "results". Once the location of the data and phylogenetic tree have been defined, the scripts should automatically produce the results outlined in "No evidence for Peto's Paradox in terrestrial vertebrates".

The assumptions of the fitted models can be checked by using the "model_diagnostic_output" script within the "model_diagnostics" directory. 

-------------------------------------------------------------------------------------------------------------------------------------------------
Body size and path-wise rate information for each species is contained with the "species_data.csv" file within the "data" directory. Raw body size data (body mass for birds and mammals and snout-vent length (SVL) for amphibians and squamate reptiles) are available from [Cooney and Thomas​](https://www.nature.com/articles/s41559-020-01321-y)​, as are the posterior distributions of the body size rate-scaled phylogenetic trees for each of the four vertebrate classes. The "posterior_path_length_measure" script within the "measuring_path_wise_rate" directory is designed to work with the posterior distributions from [Cooney and Thomas​](https://www.nature.com/articles/s41559-020-01321-y) and to calculate the median average path-wise rate for each species. Due to the size of the posterior distributions, the script is designed to run in parallel on a UNIX machine (be advised, the script is both computationally and memory intensive).
