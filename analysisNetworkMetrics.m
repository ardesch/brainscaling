% This script loads the connectivity data and computes characteristic path
% length and clustering coefficient as network metrics, and correlates them
% with cerebral volume. It also produces the scatter plots for Figures 4B
% and 4C.

%% Load data
clear; clc; rng('default');

% Set some parameters
nRegions = 50; % 25, 50, or 100 per hemisphere
hemi = 'lh'; % lh, rh
nrand = 1000; % number of randomized reference networks

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

% Select hemisphere
n = size(dwi.connectivity, 1);
if strcmp(hemi, 'lh')
    
    dwi.connectivity = dwi.connectivity(1:n/2,1:n/2,:,:);
    
elseif strcmp(hemi, 'rh')
    
    dwi.connectivity = dwi.connectivity(n/2+1:end,n/2+1:end,:,:);
    
else
    
    error("value for 'hemi' not recognized");

end

% Set density equal
minD = min(getDensity(dwi.connectivity));

for i = 1:size(dwi.connectivity, 4)

    T = thresholdDensity(dwi.connectivity(:,:,1,i), double(dwi.connectivity(:,:,1,i) > 0), minD);
    dwi.connectivity(:,:,:,i) = dwi.connectivity(:,:,:,i) .* double(T > 0);

end

%% Compute metrics
N = size(dwi.connectivity, 4);
L = nan(N, 1);
C = nan(N, 1);
Lrand = nan(N, nrand);
Crand = nan(N, nrand);

for i = 1:N
   
    fprintf('%i/%i\n', i, N);
    
    A = dwi.connectivity(:,:,1,i) > 0; % binary
    
    L(i) = mean(squareform(distance_bin(A))); % path length
    C(i) = mean(clustering_coef_bu(A)); % clustering coefficient
    
    if L(i) == Inf
        
        L(i) = NaN; % set to NaN
        fprintf('L infinite for %s...\n', dwi.species{i});
        % two macaques have unconnected nodes in lh, and
        % 1 in rh (discarding these)
        
    end
        
    % Randomized references networks
    for j = 1:nrand
        
        R = randmio_und(A, 5);
        Lrand(i,j) = mean(squareform(distance_bin(R)));
        Crand(i,j) = mean(clustering_coef_bu(R));
        
        if Lrand(i,j) == Inf
            
            Lrand(i,j) = NaN; % set to NaN            
            
        end
    
    end
    
end

Lnorm = L ./ nanmean(Lrand, 2);
Cnorm = C ./ nanmean(Crand, 2);

% Average by species
Lavg(1:7,:) = L(1:7,:);
Lavg(8,:) = nanmean(L(strcmp(dwi.species, 'macaque'), :), 1);
Lavg(9,:) = nanmean(L(strcmp(dwi.species, 'chimpanzee'), :), 1);
Lavg(10,:) = nanmean(L(strcmp(dwi.species, 'human'), :), 1);
Lavg(11,:) = nanmean(L(strcmp(dwi.species, 'bonobo'), :), 1);
Lavg(12,:) = nanmean(L(strcmp(dwi.species, 'gorilla'), :), 1);
Lavg(13,:) = nanmean(L(strcmp(dwi.species, 'orangutan'), :), 1);

Cavg(1:7,:) = C(1:7,:);
Cavg(8,:) = mean(C(strcmp(dwi.species, 'macaque'), :), 1);
Cavg(9,:) = mean(C(strcmp(dwi.species, 'chimpanzee'), :), 1);
Cavg(10,:) = mean(C(strcmp(dwi.species, 'human'), :), 1);
Cavg(11,:) = mean(C(strcmp(dwi.species, 'bonobo'), :), 1);
Cavg(12,:) = mean(C(strcmp(dwi.species, 'gorilla'), :), 1);
Cavg(13,:) = mean(C(strcmp(dwi.species, 'orangutan'), :), 1);

LnormAvg(1:7,:) = Lnorm(1:7,:);
LnormAvg(8,:) = nanmean(Lnorm(strcmp(dwi.species, 'macaque'), :), 1);
LnormAvg(9,:) = nanmean(Lnorm(strcmp(dwi.species, 'chimpanzee'), :), 1);
LnormAvg(10,:) = nanmean(Lnorm(strcmp(dwi.species, 'human'), :), 1);
LnormAvg(11,:) = nanmean(Lnorm(strcmp(dwi.species, 'bonobo'), :), 1);
LnormAvg(12,:) = nanmean(Lnorm(strcmp(dwi.species, 'gorilla'), :), 1);
LnormAvg(13,:) = nanmean(Lnorm(strcmp(dwi.species, 'orangutan'), :), 1);

CnormAvg(1:7,:) = Cnorm(1:7,:);
CnormAvg(8,:) = mean(Cnorm(strcmp(dwi.species, 'macaque'), :), 1);
CnormAvg(9,:) = mean(Cnorm(strcmp(dwi.species, 'chimpanzee'), :), 1);
CnormAvg(10,:) = mean(Cnorm(strcmp(dwi.species, 'human'), :), 1);
CnormAvg(11,:) = mean(Cnorm(strcmp(dwi.species, 'bonobo'), :), 1);
CnormAvg(12,:) = mean(Cnorm(strcmp(dwi.species, 'gorilla'), :), 1);
CnormAvg(13,:) = mean(Cnorm(strcmp(dwi.species, 'orangutan'), :), 1);

%% Correlate network metrics to cerebral volume

% Save results for PGLS
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol'}));
outTable.cort_surf = cort_surf;
outTable.Lnorm = LnormAvg;
outTable.Cnorm = CnormAvg;
writetable(outTable, sprintf('pgls/networkMetrics_%s.csv', hemi));

% Also save normalized results (easier interpretation of PGLS coefficients)
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol'}));
outTable.cort_surf = zscore(log(cort_surf));
outTable.Cnorm = zscore(CnormAvg);
outTable.Lnorm = zscore(LnormAvg);
outTable.SupraTentorialVol = zscore(log(aseg.SupraTentorialVol));
writetable(outTable, sprintf('pgls/networkMetricsNormalized_%s.csv', hemi));

% Load PGLS results
pgls = readtable('pgls/pglsNetworkMetrics.csv');
pgls(:,1) = []; % remove first column

if strcmp(hemi, 'lh')

    % Path length (LH)
    y = LnormAvg;
    x = aseg.SupraTentorialVol;
    mdl = fitlm(log(x), y);
    fB = figureLinearLog(y, x, pgls.slope_estimate(2), pgls.intercept_estimate(2), ...
    'Normalized path length', 'Cerebral volume (mm^3)', ...
    'confidenceBands', mdl);
    
    % Clustering (LH)
    y = CnormAvg;
    x = aseg.SupraTentorialVol;
    mdl = fitlm(log(x), y);
    fC = figureLinearLog(y, x, pgls.slope_estimate(1), pgls.intercept_estimate(1), ...
    'Normalized clustering coefficient', 'Cerebral volume (mm^3)', ...
    'confidenceBands', mdl);

elseif strcmp(hemi, 'rh')

    % Path length (RH)
    y = LnormAvg;
    x = aseg.SupraTentorialVol;
    mdl = fitlm(log(x), y);
    fB = figureLinearLog(y, x, pgls.slope_estimate(4), pgls.intercept_estimate(4), ...
    'Normalized path length', 'Cerebral volume (mm^3)', ...
    'confidenceBands', mdl);

    % Clustering (RH)
    y = CnormAvg;
    x = aseg.SupraTentorialVol;
    mdl = fitlm(log(x), y);
    fC = figureLinearLog(y, x, pgls.slope_estimate(3), pgls.intercept_estimate(3), ...
    'Normalized clustering coefficient', 'Cerebral volume (mm^3)', ...
    'confidenceBands', mdl);

end