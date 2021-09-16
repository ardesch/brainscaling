% This script loads the connectivity data and computes the rich-club 
% coefficient and compares it to a null distribution of randomized
% reference networks.

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

%% Compute rich-club coefficient
N = size(dwi.connectivity, 4);
n = size(dwi.connectivity, 1);
RC = nan(N, n);
RCrand = nan(N, n, nrand);

for i = 1:N
   
    fprintf('%i/%i\n', i, N);
    
    A = dwi.connectivity(:,:,1,i) > 0; % binary
    
    iRC = rich_club_bu(A); % rich-club coefficient for each k-level
    RC(i,1:length(iRC)) = iRC;
   
        
    % Randomized references networks
    for j = 1:nrand
        
        R = randmio_und(A, 5);
        jRC = rich_club_bu(R);
        RCrand(i,1:length(jRC),j) = jRC;
        
    end
    
end

% Assess significance of rich-club coefficient
pRC = nan(N, n);
for i = 1:N
    
    % For each k-level, get the non-exceedance probability of the observed
    % rich-club coefficient at that k vs the distribution of randomized
    % reference networks at that same k.
    iP = diag(invprctile(squeeze(RCrand(i,:,:))', squeeze(RC(i,:))'));
    pRC(i,:) = (100-iP)/100*2; % two-tailed p-values

end

% Set p-values of exactly 0 to NaN (happens in the absence of a nice
% distribution of rich-club coefficients in the reference networks)
pRC(pRC == 0) = NaN;
pRC(pRC > 1) = 1; % cap p-values to 1 (can be > 1 due to two-tailed testing)

% Select an arbitrary subject to display (e.g. the last one)
toPlot = [
    [1:7]'; % first 7 species only have 1 subject
    find(strcmp(dwi.species, 'macaque'), 1, 'last');
    find(strcmp(dwi.species, 'chimpanzee'), 1, 'last');
    find(strcmp(dwi.species, 'human'), 1, 'last');
    find(strcmp(dwi.species, 'bonobo'), 1, 'last');
    find(strcmp(dwi.species, 'gorilla'), 1, 'last');
    find(strcmp(dwi.species, 'orangutan'), 1, 'last')
];

% Sort on increasing brain size
[~, I] = sort(aseg.SupraTentorialVol, 'ascend');
toPlot = toPlot(I);

% FDR correct p-values
pRCplot = pRC(toPlot,:);
[~, ~, ~, pAdj] = fdr_bh(pRCplot);

%% Figure
for i = 1:length(toPlot)
   
    % Plot rich-club curve
    iidx = toPlot(i);
    f = figure('Color','white','Position',[100,100,100,150]);
    plot(RC(iidx,:), 'r'); hold on; 
    
    % Plot mean and 95% CI rich-club curves of reference networks
    mvals = mean(squeeze(RCrand(iidx,:,:)), 2);
    plot(mvals, '--', 'Color', [0.4,0.4,0.4]);
    box('off');
    x = squeeze(RCrand(iidx,:,:));
    CIFn = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); % 95% CI, see https://nl.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
    yvals = CIFn(x', 95)';
    toRemove = sum(isnan(yvals), 2);
    yvals = yvals(toRemove == 0, :);
    xvals = find(toRemove == 0);
    patchx = [xvals; xvals(end:-1:1)];
    patchy = [yvals(:,1); yvals(end:-1:1,2)];
    patch('XData', patchx, 'YData', patchy, ...
        'FaceColor', [0.5,0.5,0.5], 'FaceAlpha', 0.1, 'EdgeAlpha', 0);
    ylabel('phi(k)');
    xlabel('Degree (k)');
    title(sprintf('%s', dwi.species{iidx}));
    xlim([0,floor(0.6*n)]);
    ylim([0,1.2]);
    
    % Add asterisks for significance
    xpos = find(pAdj(i,:) < 0.05); 
    hold on;
    plot(xpos, repmat(1.1,length(xpos),1), '*k');
        
end