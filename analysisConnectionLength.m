% This script loads the connectivity data, divides the connections into
% bins based on their connection length relative to the anterior-posterior
% distance of the corresponding species' brain, and correlates this
% normalized connection length with cerebral volume. It also produces the
% heatmap used in Figure 4A.

%% Load data
clear; clc; rng('default');

% Set some parameters
nRegions = 50; % 25, 50, or 100 per hemisphere
nBins = 10; % number of connectin lenght bins

% Load the volumetric and connectivity data
load('data/volumetricData.mat');
load(sprintf('data/connectivityData%i.mat', nRegions));

% Calculate some overall metrics
cort_surf = sum(aparc_area{:,2:end}, 2);
cort_surf = cort_surf(1:13); % remove gibbon (no dwi)
aseg = aseg(1:13,:); % remove gibbon (no dwi)

% Load anterior-posterior distance data
dataAP = readtable('data/anteriorPosteriorDistance.csv');
AP = dataAP.AP;
[lia, locb] = ismember(aseg.subject, dataAP.subject); % match order of subjects
AP = AP(locb);
dataAP.subject = dataAP.subject(locb);

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

% Remove interhemispheric connections
n = size(dwi.connectivity, 1);
mask = zeros(n);
mask(1:n/2,1:n/2) = 1;
mask(n/2+1:end, n/2+1:end) = 1;
dwi.connectivity = dwi.connectivity .* mask;

%% Correlate AP distance to proportion of long fibers

% Make bins based on AP length for each species
% Store the proportion of fibers that belong to each bin
N = size(dwi.connectivity, 4);
connLengthProps = nan(N, nBins);
propLongerThanAP = nan(N,1);

for i = 1:N
    
    A = dwi.connectivity(:,:,2,i); % fiber length
    idxSpecies = strcmp(dataAP.subject, dwi.species{i});
    iAP = AP(idxSpecies); % AP of current species    
    binWidth = iAP/nBins;
    
    for j = 1:nBins-1
         
        % Proportion of all connections inside the current fiber length bin
        connLengthProps(i,j) = nnz(A > binWidth*(j-1) & A <= binWidth*j)/nnz(A);
                
    end

    % For the highest bin: include any fibers longer than AP as well
    connLengthProps(i,nBins) = nnz(A > binWidth*(nBins-1))/nnz(A);
    propLongerThanAP(i) = nnz(A > iAP)/nnz(A);
    
    
end

% Average across species with multiple subjects
connLengthPropsAvg(1:7,:) = connLengthProps(1:7,:);
connLengthPropsAvg(8,:) = nanmean(connLengthProps(strcmp(dwi.species, 'macaque'), :), 1);
connLengthPropsAvg(9,:) = nanmean(connLengthProps(strcmp(dwi.species, 'chimpanzee'), :), 1);
connLengthPropsAvg(10,:) = nanmean(connLengthProps(strcmp(dwi.species, 'human'), :), 1);
connLengthPropsAvg(11,:) = nanmean(connLengthProps(strcmp(dwi.species, 'bonobo'), :), 1);
connLengthPropsAvg(12,:) = nanmean(connLengthProps(strcmp(dwi.species, 'gorilla'), :), 1);
connLengthPropsAvg(13,:) = nanmean(connLengthProps(strcmp(dwi.species, 'orangutan'), :), 1);

% Make a heatmap of connection length proportion for each bin and species
[~, I] = sort(aseg.SupraTentorialVol);
y = 1:length(aseg.subject);
fA = figure('Color', 'white', 'Position', [0,0,400,200]);
imagesc(1:nBins,y(end:-1:1),flipud(connLengthPropsAvg)); colorbar;
xlabel('Connection length bin');
axis('xy');
colormap(makeColorMap([1,1,1], [0.9, 0.1, 0], 1000));
yticks([1:length(aseg.subject)]); yticklabels(aseg.subject(I));
axis equal tight;

% Group the connections into short and long connections based on some
% cutoff (took bins 1-2 and bins 3-10 here)
propLong = sum(connLengthPropsAvg(:,3:end), 2);
% figure; scatter(log(aseg.SupraTentorialVol), propLong);
propShort = sum(connLengthPropsAvg(:,1:2),2);
% figure; scatter(log(aseg.SupraTentorialVol), propShort);

% Save normalized values for PGLS
outTable = aseg(:, ismember(aseg.Properties.VariableNames, ...
    {'subject'}));
outTable.propShort = zscore(propShort);
outTable.propLong = zscore(propLong);
outTable.SupraTentorialVol = zscore(log(aseg.SupraTentorialVol));
writetable(outTable, 'pgls/connectionLengthNormalized.csv');
