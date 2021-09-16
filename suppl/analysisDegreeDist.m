% This script loads the connectivity data and computes the degree 
% distribution and skewness, and correlates it with cerebral volume.

%% Load data
clear; clc; rng('default');

% Set some parameters
nRegions = 50; % 25, 50, or 100 per hemisphere
hemi = 'lh'; % lh, rh

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

%% Calculate degree distribution and skewness
% Prepare some variables
N = size(dwi.connectivity, 4);
skew = nan(N, 1);

% Get degrees and skewness
for i = 1:N
   
    A = dwi.connectivity(:,:,1,i) > 0;
    k = sum(A); % degrees
    skew(i) = skewness(k, 0);
       
end

% Average by species
skewAvg(1:7,:) = skew(1:7,:);
skewAvg(8,:) = nanmean(skew(strcmp(dwi.species, 'macaque'), :), 1);
skewAvg(9,:) = nanmean(skew(strcmp(dwi.species, 'chimpanzee'), :), 1);
skewAvg(10,:) = nanmean(skew(strcmp(dwi.species, 'human'), :), 1);
skewAvg(11,:) = nanmean(skew(strcmp(dwi.species, 'bonobo'), :), 1);
skewAvg(12,:) = nanmean(skew(strcmp(dwi.species, 'gorilla'), :), 1);
skewAvg(13,:) = nanmean(skew(strcmp(dwi.species, 'orangutan'), :), 1);

%% Save for PGLS
% Save normalized results (easier interpretation of PGLS coefficients)
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject', 'SupraTentorialVol'}));
outTable.skew = zscore(skewAvg);
outTable.SupraTentorialVol = zscore(log(aseg.SupraTentorialVol));
writetable(outTable, sprintf('degreeDistNormalized_%s.csv', hemi));
