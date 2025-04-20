# Simulation and Real-World Model Evaluation

## Overview

This project performs simulations and evaluates model performance on both synthetic and real-world data. Follow the steps below to run the analysis.

## Steps

1. **Install Dependencies**  
   Run `install_dependencies.R` to install all the required libraries.

2. **Run Simulation**  
   Execute `simulation.R` to:
   - Estimate model coefficients (betas)
   - Compute mean squared error (MSE) for both training and test datasets
   - Save the results for further analysis

3. **Generate Plots**  
   Run `plots.R` to visualize:
   - Boxplots of the mean squared error (MSE)
   - Estimated covariates across models

4. **Evaluate Model Performance**  
   Run `table_TPR_NPR.R` to calculate the following performance metrics:
   - True Positive Rate (TPR)
   - True Negative Rate (TNR)
   - Positive Predictive Value (PPV)
   - Negative Predictive Value (NPV)

   The results will be saved as `table_TPR_NPR.RData`.

5. **Apply to Real-World Data**  
   Run `real-life-study.R` to apply the analysis on the real dataset `Clean_Dataset.csv`.
