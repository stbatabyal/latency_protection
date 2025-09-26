## Title: MTBI treatment program
## author: Saikat Batabyal and Vitaly Ganusov
## date: 2025-09-24

## Overview
Stopping a program of community-wide prophylactic treatment of Mtb infection may result in increase in TB cases.
Benefits of MTBI treatment program depend on TB prevalence and immunological details. Here we investigate the net difference between the benefits of reduction during MTBI treatment program (AOC)
and the increase of TB individuals after stopping the treatment program (AUC).

## Files
*   `[epdm_dynamics_one_stage_high-low_short.R]`: Stopping policy to treat MTBI may result in surge of TB cases (Figure 3,4).
*   `[epdm_net_effect_one_stage_high-low_short.R]`: Net beneficial results are observed in MTBI treatment. program (Figure 5).

## Usage
Explain how to run the R code.

1.  Set the working directory to the project folder.
2.  Install necessary packages by running `install.packages()`.
3.  Source the main R script using `source("[epdm_net_effect_one_stage_high-low_short.R]")`.
4.  You can use the same method to install any other packages that are reported as missing package during execution.
5.  To execute the code line by line, highlight the desired section and press Ctrl + Enter (for Windows/Linux) or Cmd + Enter (for macOS). Ensure that the cursor is placed within the script editor pane when performing this action. Alternatively, to run the entire script at once, you may press Ctrl + A to select all code, and     then click the "Run" button located at the top-right corner of the editor interface or press Ctrl + Enter.

## Dependencies
*   R version [e.g., 4.1.2]
*   List of R packages (installation) required to run the code.
    *   `deSolve`
    *   `rootSolve`
    *   `bayestestR`
    *   `dplyr`
    *   `ggplot2`
    *   `viridis`
    *   `scales`
    
    
## Overview
This study aims to investigate the influence of prior Mtb infection on both the timing and pathological severity of progression to active tuberculosis.

## Files
*   `[fitting_models_to_progression_data-Badger-1949.nb]`: Mtb infection (MTBI) delays TB development upon re-exposure (Figure 2).

## Data Availability
*  `[badger_art49-all-edited.csv]`: Published data from doi: 10.1164/art.1949.60.3.305

## Dependencies
*   Mathematica version [12]
    
