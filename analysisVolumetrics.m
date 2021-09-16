% This script loads the volumetric data extracted from FreeSurfer and the
% model estimates calculated by PGLS to produce the data panels in Figures
% 2 and 3.

%% Load data
clear; clc;
load('data/volumetricData.mat');

% Calculate some overall metrics
cort_surf = sum(aparc_area{:,2:end}, 2);
cort_thickness = mean([aparc_thickness.lh_MeanThickness_thickness, aparc_thickness.rh_MeanThickness_thickness], 2);
total_cc_area = sum(cc_area{:,2:end},2);

% Import PGLS model estimates into MATLAB to create figures
pgls = readtable('pgls/pglsVolumetrics.csv');
pgls(:,1) = []; % remove first column

%% Figure 2B: Cortical surface vs supratentorial volume
fB = figureLogLog(cort_surf, aseg.SupraTentorialVol, pgls.slope_estimate(1), pgls.intercept_estimate(1), ...
    'Cortical surface area (mm^2)', 'Cerebral volume (mm^3)', ...
    'addLine', 2/3, 'confidenceBands', fitlm(log(aseg.SupraTentorialVol), log(cort_surf)));

%% Figure 2C: White matter volume vs gray matter volume
% Pagel's lambda is not 0 so need to calculate confidence bands using the
% adjusted beta from PGLS analysis
x = aseg.CortexVol;
y = aseg.CerebralWhiteMatterVol;
logxmin = min(log(x));
logxmax = max(log(x));
logxrange = logxmax - logxmin;
xlogci = [logxmin:logxrange/1000:logxmax]'; % define some points for x axis on log scale
beta = [pgls.intercept_estimate(4); pgls.slope_estimate(4)];
mdl = fitlm(log(x), log(y));
sigma = mdl.CoefficientCovariance;
dfe = length(x) - 1; % df = N - 1
ylogci = predictCI(xlogci, beta, sigma, dfe);
yci = exp(ylogci);
xci = exp(xlogci);
datapoints(:,1) = xci;
datapoints(:,2:3) = yci;

% Create the figure
fD = figureLogLog(aseg.CerebralWhiteMatterVol, aseg.CortexVol, pgls.slope_estimate(4), pgls.intercept_estimate(4), ...
    'White matter volume (mm^3)', 'Cortical gray matter volume (mm^3)', ...
    'addLine', 1, 'plotCustomBand', datapoints);

%% Figure 2D: Cortical surface vs white matter volume
% Pagel's lambda is not 0 so need to calculate confidence bands using the
% adjusted beta from PGLS analysis
x = aseg.CerebralWhiteMatterVol;
y = cort_surf;
logxmin = min(log(x));
logxmax = max(log(x));
logxrange = logxmax - logxmin;
xlogci = [logxmin:logxrange/1000:logxmax]'; % define some points for x axis on log scale
beta = [pgls.intercept_estimate(5); pgls.slope_estimate(5)];
mdl = fitlm(log(x), log(y));
sigma = mdl.CoefficientCovariance;
dfe = length(x) - 1; % df = N - 1
ylogci = predictCI(xlogci, beta, sigma, dfe);
yci = exp(ylogci);
xci = exp(xlogci);
datapoints(:,1) = xci;
datapoints(:,2:3) = yci;
fE = figureLogLog(cort_surf, aseg.CerebralWhiteMatterVol, pgls.slope_estimate(5), pgls.intercept_estimate(5), ...
    'Cortical surface area (mm^2)', 'White matter volume (mm^3)', ...
    'addLine', 2/3, 'plotCustomBand', datapoints);

%% Figure 3C: Corpus callosum surface area vs cortical surface area
f = figureLogLog(total_cc_area, cort_surf, pgls.slope_estimate(7), pgls.intercept_estimate(7), ...
    'CC cross-sectional area (mm^2)', 'Cortical surface area (mm^2)', ...
    'addLine', 1, 'confidenceBands', fitlm(log(cort_surf), log(total_cc_area)));
