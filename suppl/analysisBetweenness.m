% This script loads the connectivity data and computes nodal betweenness
% centrality and correlates it with cerebral volume.

%% Load data
clear; clc; rng('default');

% Set some parameters
nRegions = 50; % 25, 50, or 100 per hemisphere
hemi = 'lh'; % lh, rh
nrand = 1000; % number of randomized reference networks

% Load the volumetric and connectivity data
load('../data/volumetricData.mat');
load(sprintf('../data/connectivityData%i.mat', nRegions));

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
BC = nan(N, 1);
BCskew = nan(N, 1);
BCrand = nan(N, nrand);
BCskewRand = nan(N, nrand);

for i = 1:N
   
    fprintf('%i/%i\n', i, N);
    
    A = dwi.connectivity(:,:,1,i) > 0; % binary
    
    iBC = betweenness_bin(A);
    BC(i) = mean(iBC); % average betweeness centrality
    BCskew(i) = skewness(iBC);
    
    if BC(i) == Inf
        
        BC(i) = NaN; % set to NaN
        fprintf('BC infinite for %s...\n', dwi.species{i});
        % two macaques have unconnected nodes in lh, and
        % 1 in rh (discarding these)
        
    end
        
    % Randomized references networks 
    for j = 1:nrand

        R = randmio_und(A, 5);
        ijBC = betweenness_bin(R);
        BCrand(i,j) = mean(ijBC);
        BCskewRand(i,j) = skewness(ijBC);

        if BCrand(i,j) == Inf

            BCrand(i,j) = NaN; % set to NaN            

        end

    end
    
end

% Normalize the metrics   
BC = BC ./ nanmean(BCrand, 2);
BCskew = BCskew ./ nanmean(BCskewRand, 2);

% Average by species
BCavg(1:7,:) = BC(1:7,:);
BCavg(8,:) = mean(BC(strcmp(dwi.species, 'macaque'), :), 1);
BCavg(9,:) = mean(BC(strcmp(dwi.species, 'chimpanzee'), :), 1);
BCavg(10,:) = mean(BC(strcmp(dwi.species, 'human'), :), 1);
BCavg(11,:) = mean(BC(strcmp(dwi.species, 'bonobo'), :), 1);
BCavg(12,:) = mean(BC(strcmp(dwi.species, 'gorilla'), :), 1);
BCavg(13,:) = mean(BC(strcmp(dwi.species, 'orangutan'), :), 1);

BCskewAvg(1:7,:) = BCskew(1:7,:);
BCskewAvg(8,:) = mean(BCskew(strcmp(dwi.species, 'macaque'), :), 1);
BCskewAvg(9,:) = mean(BCskew(strcmp(dwi.species, 'chimpanzee'), :), 1);
BCskewAvg(10,:) = mean(BCskew(strcmp(dwi.species, 'human'), :), 1);
BCskewAvg(11,:) = mean(BCskew(strcmp(dwi.species, 'bonobo'), :), 1);
BCskewAvg(12,:) = mean(BCskew(strcmp(dwi.species, 'gorilla'), :), 1);
BCskewAvg(13,:) = mean(BCskew(strcmp(dwi.species, 'orangutan'), :), 1);

%% Correlate network metrics to cerebral volume
% Save normalized results (easier interpretation of PGLS coefficients)
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol'}));
outTable.avg = zscore(BCavg);
outTable.skew = zscore(BCskewAvg);
outTable.SupraTentorialVol = zscore(log(aseg.SupraTentorialVol));
writetable(outTable, sprintf('betweennessNormalized_%s.csv', hemi));
