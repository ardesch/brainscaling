% This script checks whether the resolution (brain volume divided by
% voxel size, i.e. number of voxels in the brain) is correlated with 
% brain size or connectivity asymmetry.

%% Load data
load('../data/voxelSize.mat');
load('../data/volumetricData.mat');
conn = readtable('../pgls/connectivityAsymmetry.csv');

% Match order of subjects between voxelSizes and aseg/connectivity data
[lia, locb] = ismember(aseg.subject, subjects);
voxelSizeStructural = voxelSizeStructural(locb);
voxelSizeDWI = voxelSizeDWI(locb);

%% Structural voxel size
% First check whether voxel size is correlated with cerebral volume
mdl = fitlm(log(aseg.SupraTentorialVol), voxelSizeStructural);
figureLinearLog(voxelSizeStructural, aseg.SupraTentorialVol, ...
    mdl.Coefficients.Estimate(2), mdl.Coefficients.Estimate(1),...
    'Voxel size T1 (mm)', 'Cerebral volume (mm^3)');  
[r, p] = corr(log(aseg.SupraTentorialVol), voxelSizeStructural);

% Resolution: cerebral volume divided by voxel size^3
resStructural = aseg.SupraTentorialVol ./ (voxelSizeStructural.^3);
mdl = fitlm(log(aseg.SupraTentorialVol), resStructural);
figureLinearLog(resStructural, aseg.SupraTentorialVol, ...
    mdl.Coefficients.Estimate(2), mdl.Coefficients.Estimate(1),...
    'Voxel resolution T1', 'Cerebral volume (mm^3)');   
[r, p] = corr(log(aseg.SupraTentorialVol), resStructural);

%% DWI voxel size
aseg(14,:) = []; % remove gibbon (no DWI)
voxelSizeDWI(14) = []; % remove gibbon (no DWI)

% First check whether voxel size is correlated with cerebral volume
mdl = fitlm(log(aseg.SupraTentorialVol), voxelSizeDWI);
figureLinearLog(voxelSizeDWI, aseg.SupraTentorialVol, ...
    mdl.Coefficients.Estimate(2), mdl.Coefficients.Estimate(1),...
    'Voxel size DWI (mm)', 'Cerebral volume');
[r, p] = corr(log(aseg.SupraTentorialVol), voxelSizeDWI);

% Resolution: cerebral volume divided by voxel size^3
resDWI = aseg.SupraTentorialVol ./ (voxelSizeDWI.^3);
mdl = fitlm(log(aseg.SupraTentorialVol), resDWI);
figureLinearLog(resDWI, aseg.SupraTentorialVol, ...
    mdl.Coefficients.Estimate(2), mdl.Coefficients.Estimate(1),...
    'Voxel resolution DWI', 'Cerebral volume', 'confidenceBands', mdl);
[r, p] = corr(log(aseg.SupraTentorialVol), resDWI);

%% Connectivity asymmetry
f = figure('Color', 'white', 'position', [200,200,290,290], 'Renderer', 'painters');
scatter(conn.conn_asymmetry, resDWI, 75, 'k', 'filled');
l = lsline; l.LineStyle = '--'; l.Color = [0.8, 0.8, 0.8];
xlabel('Connectivity asymmetry');
ylabel('Voxel resolution DWI');
set(gca, 'FontSize', 14);
set(gca, 'FontName', 'Helvetica'); 
[r, p] = corr(conn.conn_asymmetry, resDWI);