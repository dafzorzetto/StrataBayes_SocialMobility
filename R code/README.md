## Code files for StrataBayes:

### General functions (in **`src`** folder):
 - **`Gibbs_Sampler.R`**:

   Gibbs sampler to estimate our proposed model: StrataBayes

 - **`Competitor.R`**:

   Gibbs sampler to estimate the benchmark model by Shwartz Li and Mealli(2011)
 
### Reproduce the results in the Simulation Study Section (in **`simulation study`** folder):
 - **`Generating_Mechanism.R`**:

   generate different settings
 - **`Estimation_model.R`**:

   estimate proposed model for the different simulated settings

 - **`Estimation_model_parallel.R`**:

   estimate proposed model for the different simulated settings using library **`parallel`**
   
  - **`Estimation_model_parallel_windows.R`**:

   estimate proposed model for the different simulated settings using library **`doparallel`**

  - **`Analysis_Results.R`**:

   results analysis of the proposed model against benchmark models

### Envrironmental and socioeconomic application (in **`application`** folder):
 - **`application.R`**:

   application with social mobility, PM2.5 exposure, educational attainment and confounders data
   
 - **`plots_paper.R`**:

   analysis and visualization of the results
