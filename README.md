# depth-profile-inversion
MCMC inversion code for cosmogenic nuclide bedrock depth profiles

This MCMC inversion code generates the model output presented and discussed in the manuscript: Ice-sheet burial and erosion inferred from cosmogenic nuclide bedrock depth profiles: implications for the glaciation history of northeastern Fennoscandia, Andersen et al. Quaternary Science Reviews 344 (2024) 
https://doi.org/10.1016/j.quascirev.2024.109010.

Note 1: At present the code is set up to handle 10Be and 26Al (both are mandatory).

Note 2: Production rates are presently calculated using the 'St' scaling scheme, which does not take time-variable production into account.

Note 3: files with filenames containing '_joined' handles depth profiles of different lengths (e.g. joined inversion of the two tor profiles)

Note 4: start by calling: addpath('Functions'); addpath('data'); in the MATLAB command window.

Primary functions

1. scenarios.m: function that returns variables for each of the 16 scenarios for ice-sheet burial and erosion that is discussed in the manuscript.

2a. compile_synthetic_data.m: Calculates production paramters and 'true' concentrations for the synthetic scenarios and saves data structure in .mat format in the 'data' folder.

2b. compile_data.m: Reads data from excelfile for one or more sites, calculates production and saves data structure in .mat format in the 'data' folder.

3. bedrockMCvJ4_true.m: Performs inversion of depth profiles or surface samples. For a quick test of performance reduce number of walkers and number of models per walker in the code. To invert synthetic scenarios used in main manuscript call this function in a for-loop (consider parallellising):

for j=1:16, %loop scenarios
    bedrockMCvJ4_true(1,j,1);
end

Output is saved to the 'models' folder 

4. forward_bedrock.m: Forward calculation of nuclide concentrations in a sample given the exposure and exhumation history parameters.

5. compile_results_mn5.m: Compiles inversion model results for all samples from each site for use in plotting scripts.

6. compareInversions_synthetic.m and compareInversions_data.m: Creates figures equivalent to subplots in Fig. 3 and Fig. 9 in manuscript respectively. Shows results of inversion, true values for synthetic scenarios in red.

8. SummarizeInversions.m: create figure 4 in manuscript - requires that all 16 scenarios have been inverted.

9. makereport_syn_dp.m: Plots information about walkers, parameter distributions, model-data fit, and saves report as pdf in models/reports folder.

Fig1.m and Fig2.m creates figure 1 and figure 2 of the manuscript respectively.

The data folder contains: 
- Naakakarhakka.xlsx and Lamuvaara.xlsx and stores mat-files generated from compile_data and compile_synthetic_data functions.
- d18Ocurves.mat (Lisiecki & Raymo 2005), Data downloaded from https://lorraine-lisiecki.com/stack.html
- mat-files with production parameters and data for synthetic and real data, these are overwritten if 'compile_synthetic_data_fw_dp_vJ4.m' or 'compile_data_vJ4.m' are run.
- N10all.xlsx and N26all.xlsx which are used to find the likely error on a synthetic datapoint given the history of AMS measurements at The AARAMS Centre.

The Function folder contains
- codes used for calculation of production rates (stone2000.m)
- Muon production code by Balco et al., 2017
- export_fig package from https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig
- Sub-functions for support and visualisation


References

- Balco, G. (2017) Production rate calculations for cosmic-ray-muon-produced 10Be and 26Al benchmarked against geological calibration data. Quaternary Geochronology, 39, 150-173.
- Lisiecki, L. E., & Raymo, M. E. (2005). A Pliocene‐Pleistocene stack of 57 globally distributed benthic δ18O records. Paleoceanography, 20(1).

/Jane Lund Andersen, June 2024
