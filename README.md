# Resolving Peto's Paradox

### George Butler<sup>1,</sup>, Joanna Baker<sup>2</sup>, Sarah R. Amend<sup>1</sup>, Kenneth J. Pienta<sup>1</sup>, Chris Venditti<sup>2</sup>

<sup><sup>1</sup>Cancer Ecology Center, The Brady Urological Institute, Johns Hopkins School of Medicine, Baltimore, USA</sup>

<sup><sup>2</sup>School of Biological Sciences, University of Reading, Reading, UK</sup>

<p align="center">
  <img width="380" height="355" src="/example_image/exceptional_birds_and_mammals.png">
</p>


All code is written in R with the majority of the scripts utilizing the [MCMCglmm](https://cran.r-project.org/web/packages/MCMCglmm/index.html) library to fit Multivariate Phylogenetic Generalized Linear Mixed Models (MPGLMMs). Large portions of code are designed to run on a UNIX machine and take advantage of multiple threads where possible. If you wish to run this code on a Windows machine please change the "pbmclapply" function to "lapply" with the understanding that this will increase the run time considerably.  


