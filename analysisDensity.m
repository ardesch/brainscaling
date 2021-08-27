% This script loads the connectivity data, calculates the network density
% per subject, and correlates this with cerebral volume.

%% Load data
clear; clc; rng('default');

% Set some parameters
nRegions = 50; % 25, 50, or 100 per hemisphere

% Load the volumetric and connectivity data
load('data/volumetricData.mat');
load(sprintf('data/connectivityData%i.mat', nRegions));

% Calculate some overall metrics
cort_surf = sum(aparc_area{:,2:end}, 2);
cort_surf = cort_surf(1:13); % remove gibbon (no dwi)
aseg = aseg(1:13,:); % remove gibbon (no dwi)

% Find outliers per species
% Macaque
tmp_macaque_idx = find(strcmp(dwi.species, 'macaque'));
tmp_macaque_dens = getDensity(dwi.connectivity(:,:,1,tmp_macaque_idx));
loCutoff = median(tmp_macaque_dens) - 1.5*iqr(tmp_macaque_dens);
hiCutoff = median(tmp_macaque_dens) + 1.5*iqr(tmp_macaque_dens);
outliers = tmp_macaque_dens < loCutoff | tmp_macaque_dens > hiCutoff;
tmp_macaque_outliers = tmp_macaque_idx(outliers);

% Chimp
tmp_chimp_idx = find(strcmp(dwi.species, 'chimpanzee'));
tmp_chimp_dens = getDensity(dwi.connectivity(:,:,1,tmp_chimp_idx));
loCutoff = median(tmp_chimp_dens) - 1.5*iqr(tmp_chimp_dens);
hiCutoff = median(tmp_chimp_dens) + 1.5*iqr(tmp_chimp_dens);
outliers = tmp_chimp_dens < loCutoff | tmp_chimp_dens > hiCutoff;
tmp_chimp_outliers = tmp_chimp_idx(outliers);

% Human
tmp_human_idx = find(strcmp(dwi.species, 'human'));
tmp_human_dens = getDensity(dwi.connectivity(:,:,1,tmp_human_idx));
loCutoff = median(tmp_human_dens) - 1.5*iqr(tmp_human_dens);
hiCutoff = median(tmp_human_dens) + 1.5*iqr(tmp_human_dens);
outliers = tmp_human_dens < loCutoff | tmp_human_dens > hiCutoff;
tmp_human_outliers = tmp_human_idx(outliers);

% Remove outliers
outliers = [tmp_macaque_outliers; tmp_chimp_outliers; tmp_human_outliers];
dwi.connectivity(:,:,:,outliers) = [];
dwi.species(outliers) = [];
dwi.regionProperties(:,:,outliers) = [];

%% Calculate density of intrahemispheric connectivity

% Get the left and right hemispheres
n = size(dwi.connectivity, 1);
dwi.lh = dwi.connectivity(1:n/2,1:n/2,:,:);
dwi.rh = dwi.connectivity(n/2+1:end,n/2+1:end,:,:);

% Left hemisphere
N = size(dwi.lh, 4);
densLH = nan(N, 1);

for i = 1:N
   
    A = dwi.lh(:,:,1,i) > 0;
    densLH(i) = density_und(A);
  
end

% Average across species with multiple subjects
densLHavg(1:7,:) = densLH(1:7,:);
densLHavg(8,:) = nanmean(densLH(strcmp(dwi.species, 'macaque'), :), 1);
densLHavg(9,:) = nanmean(densLH(strcmp(dwi.species, 'chimpanzee'), :), 1);
densLHavg(10,:) = nanmean(densLH(strcmp(dwi.species, 'human'), :), 1);
densLHavg(11,:) = nanmean(densLH(strcmp(dwi.species, 'bonobo'), :), 1);
densLHavg(12,:) = nanmean(densLH(strcmp(dwi.species, 'gorilla'), :), 1);
densLHavg(13,:) = nanmean(densLH(strcmp(dwi.species, 'orangutan'), :), 1);

% Right hemisphere
N = size(dwi.rh, 4);
densRH = nan(N, 1);

for i = 1:N
   
    A = dwi.rh(:,:,1,i) > 0;
    densRH(i) = density_und(A);
  
end

densRHavg(1:7,:) = densRH(1:7,:);
densRHavg(8,:) = nanmean(densRH(strcmp(dwi.species, 'macaque'), :), 1);
densRHavg(9,:) = nanmean(densRH(strcmp(dwi.species, 'chimpanzee'), :), 1);
densRHavg(10,:) = nanmean(densRH(strcmp(dwi.species, 'human'), :), 1);
densRHavg(11,:) = nanmean(densRH(strcmp(dwi.species, 'bonobo'), :), 1);
densRHavg(12,:) = nanmean(densRH(strcmp(dwi.species, 'gorilla'), :), 1);
densRHavg(13,:) = nanmean(densRH(strcmp(dwi.species, 'orangutan'), :), 1);

%% Save results for PGLS

% Unnormalized (for figures)
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol'}));
outTable.cort_surf = cort_surf;
outTable.SupraTentorialVol = aseg.SupraTentorialVol;
outTable.densityLH = densLHavg;
outTable.densityRH = densRHavg;
writetable(outTable, 'pgls/density.csv');

% Normalized (for regression coefficients)
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol'}));
outTable.cort_surf = zscore(log(cort_surf));
outTable.SupraTentorialVol = zscore(log(aseg.SupraTentorialVol));
outTable.densityLH = zscore(densLHavg);
outTable.densityRH = zscore(densRHavg);
writetable(outTable, 'pgls/densityNormalized.csv');

%% Figures

% Load PGLS slope estimates
pgls = readtable('pgls/pglsDensity.csv');
pgls(:,1) = []; % remove first column

% Supratentorial volume and density (LH)
y = densLHavg;
x = aseg.SupraTentorialVol;
mdl = fitlm(log(x), y);
f = figureLinearLog(y, x, pgls.slope_estimate(1), pgls.intercept_estimate(1), ...
'Network density', 'Cerebral volume (mm^3)', ...
'confidenceBands', mdl);
ylim([0,1]);

% Supratentorial volume and density (RH)
y = densRHavg;
x = aseg.SupraTentorialVol;
mdl = fitlm(log(x), y);
f = figureLinearLog(y, x, pgls.slope_estimate(2), pgls.intercept_estimate(2), ...
'Network density', 'Cerebral volume (mm^3)', ...
'confidenceBands', mdl);
ylim([0,1]);

