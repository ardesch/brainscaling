% This script loads the connectivity data, computes connectivity
% asymmetry based on binary topology, and correlates it to cerebral
% volume.

%% Load data
clear; clc; rng('default');

% Set some parameters
nRegions = 50; % 25, 50, or 100 per hemisphere
nrand = 1000; % number of permuted networks

% Load the volumetric and connectivity data
load('../data/volumetricData.mat');
load(sprintf('../data/connectivityData%i.mat', nRegions));

% Calculate some overall metrics
cort_surf = sum(aparc_area{:,2:end}, 2);
cort_surf = cort_surf(1:13); % remove gibbon (no dwi)
aseg = aseg(1:13,:); % remove gibbon (no dwi)

% Choose x (supratentorial volume or cortical surface area)
x = log(aseg.SupraTentorialVol);
% x = log(cort_surf); % alternative: cortical surface area

%% Prepare connectivity data
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

% Set equal density
for i = 1:size(dwi.connectivity, 4)

    T = thresholdDensity(dwi.connectivity(:,:,1,i), dwi.connectivity(:,:,1,i), 0.1);
    dwi.connectivity(:,:,:,i) = dwi.connectivity(:,:,:,i) .* double(T > 0);

end

%% Main analysis
% Initialize variables
N = size(dwi.connectivity, 1);
n = size(dwi.connectivity, 4);
asymmetry = nan(n, 1);

% Calculate binary connectivity divergence
for i = 1:n
    
    W = double(dwi.connectivity(:,:,1,i) > 0); % binary connectivity   
    W(1:N/2, N/2+1:end) = 0; % remove interhemispheric connections
    W(N/2+1:end, 1:N/2) = 0; % remove interhemispheric connections
    
    asymmetry(i) = 1 - nanmean(calcOverlapBin(W, ...
        dwi.regionDescriptions, 'includeInter', 0, 'presenceOnly', 1)); % excluding interhemispheric connections     
    
end

% Average values for species with more than one subject
avg_asymmetry(1:7) = asymmetry(1:7);
avg_asymmetry(8) = mean(asymmetry(strcmp(dwi.species, 'macaque')), 1);
avg_asymmetry(9) = mean(asymmetry(strcmp(dwi.species, 'chimpanzee')), 1);
avg_asymmetry(10) = mean(asymmetry(strcmp(dwi.species, 'human')), 1);
avg_asymmetry(11) = mean(asymmetry(strcmp(dwi.species, 'bonobo')), 1);
avg_asymmetry(12) = mean(asymmetry(strcmp(dwi.species, 'gorilla')), 1);
avg_asymmetry(13) = mean(asymmetry(strcmp(dwi.species, 'orangutan')), 1);

% Standard correlation (for PGLS see below)
[beta, pval] = corr(x, avg_asymmetry');
mdl = fitlm(x, avg_asymmetry');
b = mdl.Coefficients.Estimate(2);
    
%% Null model 
betas_null = nan(nrand, 1);
pvals_null = nan(nrand, 1);
bs_null = nan(nrand, 1);
asymmetry_null = nan(n, nrand);
y_null = nan(length(avg_asymmetry), nrand);

for k = 1:nrand
    
    fprintf('Permutation %i/%i...\n', k, nrand);

    for l = 1:n      
    
        W = double(dwi.connectivity(:,:,1,l) > 0);    
        W(1:N/2, N/2+1:end) = 0; % remove interhemispheric connections
        W(N/2+1:end, 1:N/2) = 0; % remove interhemispheric connections
        WR = W;
        WR(1:N/2, 1:N/2) = randmio_und(W(1:N/2, 1:N/2), 5); % rewire LH
        WR(N/2+1:end, N/2+1:end) = randmio_und(W(N/2+1:end, N/2+1:end), 5); % rewire RH
        
        asymmetry_null(l,k) = 1 - nanmean(calcOverlapBin(WR, ...
            dwi.regionDescriptions, 'includeInter', 0, 'presenceOnly', 1)); % excluding interhemispheric connections     
  
    end

    % Average values for species with more than one subject
    k_avg_asymmetry_null(1:7) = asymmetry_null(1:7, k);
    k_avg_asymmetry_null(8) = mean(asymmetry_null(strcmp(dwi.species, 'macaque'), k), 1);
    k_avg_asymmetry_null(9) = mean(asymmetry_null(strcmp(dwi.species, 'chimpanzee'), k), 1);
    k_avg_asymmetry_null(10) = mean(asymmetry_null(strcmp(dwi.species, 'human'), k), 1); 
    k_avg_asymmetry_null(11) = mean(asymmetry_null(strcmp(dwi.species, 'bonobo'), k), 1); 
    k_avg_asymmetry_null(12) = mean(asymmetry_null(strcmp(dwi.species, 'gorilla'), k), 1); 
    k_avg_asymmetry_null(13) = mean(asymmetry_null(strcmp(dwi.species, 'orangutan'), k), 1); 
    y_null(:,k) = k_avg_asymmetry_null;
    
    % Correlation
    [betas_null(k), pvals_null(k)] = corr(x, k_avg_asymmetry_null');
    mdl = fitlm(x, k_avg_asymmetry_null');
    bs_null(k) = mdl.Coefficients.Estimate(2);    

end

%% Save for PGLS
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol'}));
outTable.cort_surf = cort_surf;
outTable.conn_asymmetry = avg_asymmetry';
writetable(outTable, 'connectivityAsymmetryBinary.csv');

% Also save the null results
outTable = array2table(y_null);
outTable.subject = aseg.subject;
outTable = [outTable(:,end), outTable(:,1:end-1)]; % reorder
writetable(outTable, 'connectivityAsymmetryBinaryNull.csv');

% Also save normalized results (easier interpretation of PGLS coefficients)
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol'}));
outTable.cort_surf = zscore(log(cort_surf));
outTable.conn_asymmetry = zscore(avg_asymmetry)';
outTable.SupraTentorialVol = zscore(log(aseg.SupraTentorialVol));
writetable(outTable, 'connectivityAsymmetryBinaryNormalized.csv');

%% Figures
% Load PGLS slope estimates
pgls = readtable('pglsConnectivityAsymmetryBinary.csv');
pgls(:,1) = []; % remove first column

% Supratentorial volume
y = avg_asymmetry;
x = aseg.SupraTentorialVol;
mdl = fitlm(log(x), y);
fA = figureLinearLog(y, x, pgls.slope_estimate(1), pgls.intercept_estimate(1), ...
'Connectivity asymmetry', 'Cerebral volume (mm^3)', ...
'confidenceBands', mdl);
ylim([0.6*min(y), 1.05*max(y)]);

% Cortical surface area
y = avg_asymmetry;
x = cort_surf;
mdl = fitlm(log(x), y);
fB = figureLinearLog(y, x, pgls.slope_estimate(2), pgls.intercept_estimate(2), ...
'Connectivity asymmetry', 'Cortical surface area (mm^2)', ...
'confidenceBands', mdl);
xlim([fB.CurrentAxes.XLim(1), 2.5*10^5]);
ylim([0.6*min(y), 1.05*max(y)]);

% Add null model (for the x chosen in first section of this script)
figure('Color', 'white', 'position', [200,200,100,100]);
[counts,edges] = histcounts(betas_null);
edges = edges(2:end) - (edges(2)-edges(1))/2;
h = plot(edges, counts);
box('Off');
hold on;
h2 = fill([min(edges), edges, max(edges)], [0, counts, 0], [0.9, 0.9, 0.9]);
l2 = xline(beta);
l2.LineStyle = '--';
l2.Color = 'r';
ylabel('Count');
xlabel("Pearson's {\it r}");
xlim([-1,1]);
