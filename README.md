
# R package for "Identifying Temporal Pathways Using Biomarkers in the Presence of Latent Non-Gaussian Components"

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

- The code for the proposed methodology is included in **ICATemporalNetwork** folder. Please download all the files in the folder to implement the method.
  + The main function for the method is **ICATemporalNet.R** and **ICATemporalNetBoot.R** which allows bootstraps.
  + To use **ICATemporalNet.R**, it requires **estICA.R** which estimates the non-Gaussian signals and then removes them from raw data, and **temporalNet.R** which estimates the temporal network. 


 
- **Examples** folder contains examples.
   + **genData.R**: generate simulated data
   + **example.R**: an example to implement the method
   + **Sim_Scenario1.R**: simulations of Scenario 1
   + **Sim_Scenario2.R**: simulations of Scenario 2

### Main function: ICATemporalNet
#### Arguments
+ `Yts`: input data, the user must supply a list of Yts, where each element is a N*K data matrix at time t. N is sample size, K is the number of nodes.
+ `N`: sample size
+ `Ntime`: total number of time points
+ `ncomp`:  maximum number of independent components to be chosen
+  `Ta`: use t<=Ta time points used for estimating temporal network A
+  `Tc`: ues t>Tc time points used for estimating contemporaneous network Gamma

#### Value
+ `estIC`: results from non-Gaussian estimation step. output from estICA.R
+ `estRts`: residuals after removing non-Gaussian signals
+ `estS`: independent components S
+ `estUts`: non-Gaussian signals U(t)
+ `estWt`: weight matrix w(t)
+  `nIC`: number of independent components
+  `A`: temporal network
+  `Gamma`: contemporaneous network
+  `Omega`: covariance matrix of e(t), inverse of Gamma

### The arguments of other functions are described within R files.
 
