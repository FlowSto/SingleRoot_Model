# SingleRoot_Model
Single Root water uptake model for four different maize root types (brace, crown, seminal. lateral)
# Root–Soil Conductance Model

This repository contains MATLAB code for estimating radial and axial conductance and running a single-root soil interaction model including root hairs.

## Requirements

- MATLAB (tested with version 2025)

## Getting Started

### 1. Estimate Radial and Axial Conductance

Run the following script: Estimation_radialConductance.m

This script computes and saves the radial and axial conductance values required for parametrizing the main model.

### 2. Run the Main Model

After obtaining the conductance parameters, run: single_root_with_soil_domain_and_root_hairs.m

This script executes the main model and produces and saves all relevant outputs.

## Additional Script

NonLinear_SingleRoot.m  
  Calculates rhizosphere conductance of a single root segment, as shown in Figure 2 of the associated publication.

## Notes

- Ensure all scripts are in the same working directory or added to the MATLAB path.
- Output files are saved automatically during execution.
