% This script loads the volumetric data extracted from FreeSurfer,
% calculates some summary metrics, and saves the data again as a .csv file
% for use with PGLS analysis in R.

%% Load data
clear; clc;
load('data/volumetricData.mat');

%% Calculate some summary metrics
cort_surf = sum(aparc_area{:,2:end}, 2);
cort_thickness = mean([aparc_thickness.lh_MeanThickness_thickness, aparc_thickness.rh_MeanThickness_thickness], 2);
total_cc_area = sum(cc_area{:,2:end},2);

%% Save as .csv file for PGLS in R
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol', 'CortexVol', 'CerebralWhiteMatterVol', 'SubcortGrayVol'}));
outTable.CerebralGrayMatterVol = aseg.CortexVol + aseg.SubCortGrayVol;
outTable.cort_surf = cort_surf;
outTable.cort_thickness = cort_thickness;
outTable.total_cc_area = total_cc_area;
writetable(outTable, 'pgls/volumetricData.csv');
