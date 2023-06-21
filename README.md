
# R package for "Identifying Disease-Associated Biomarker Network Features Through Conditional Graphical Model"

<img src="https://img.shields.io/badge/Study%20Status-Results%20Available-yellow.svg" alt="Study Status: Results Available"> 

Time series data collected from a network of random variables are useful for identifying temporal pathways among the network nodes. Observed measurements may contain multiple sources of signals and noises, including Gaussian signals of interest and non-Gaussian noises including artifacts, structured noise, and other unobserved factors (e.g., genetic risk factors, disease susceptibility). Existing methods including vector autoregression (VAR) and dynamic causal modeling do not account for unobserved non-Gaussian components. Furthermore, existing methods cannot effectively distinguish contemporaneous relationships from temporal relations. In this work, we propose a novel method to identify latent temporal pathways using time series biomarker data collected from multiple subjects. The model adjusts for the non-Gaussian components and separates the temporal network from the contemporaneous network. Specifically, an independent component analysis (ICA) is used to extract the unobserved non-Gaussian components, and residuals are used to estimate the contemporaneous and temporal networks among the node variables based on method of moments. The algorithm is fast and can easily scale up. 

- Title: **Identifying Temporal Pathways Using Biomarkers in the Presence of Latent Non-Gaussian Components**

- Authors: **Shanghong Xie<sup>a,b</sup> (shanghongxie@gmail.com), Donglin Zeng<sup>c</sup>, and Yuanjia Wang<sup>b</sup>**

- Affiliations:
   + 1. **School of Statistics, Southwestern University of Finance and Economics, Chengdu, China**
   + 2. **Department of Biostatistics, Mailman School of Public Health, Columbia University, New York, USA**
   + 3. **Department of Biostatistics, University of North Carolina, Chapel Hill, North Carolina, USA**
  



## Setup Requirements
- R


## Code Instructions

The code for the proposed methodology is included in cNetworkC.cpp and cNetworkR.R. The arguments are described inside the code.

The main function for first stage is LmNetwork_1st and the main function for second stage is EnetLm.

SimGenerate.R includes the code to generate simulation data.

SimSample.R provides an example of simulation study.
