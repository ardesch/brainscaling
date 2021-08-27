function [overlap, regions] = calcOverlapBin(W, regionDescriptions, varargin)
% This function takes a connectivity matrix and regionDescriptions as
% input, and returns the binary overlap between each homologous region pair
% in the left and right hemisphere (for the submatrix of regions that 
% appear in both hemispheres).

connLimit = size(W,1)/2;
includeInter = 1;
presenceOnly = 0;

% Parse optional arguments
while ~isempty(varargin)
    if numel(varargin) == 1
        error('lscatter:missing_option', ...
            'Optional arguments must come in pairs.');
    end
    
    switch lower(varargin{1})
        case 'connlimit'
            assert(isnumeric(varargin{2}))
            connLimit = varargin{2};
        case 'includeinter'
            assert(isnumeric(varargin{2}))
            includeInter = varargin{2};
        case 'presenceonly'
            assert(isnumeric(varargin{2}))
            presenceOnly = varargin{2};
        otherwise
            error('option %s unknown', varargin{1});
    end

    varargin(1:2) = [];

end

assert(nnz(W < 0) == 0, ...
    'Error: weighted matrix W should only have non-negative values');

% Extract submatrix of regions that exist in both hemispheres
lh = startsWith(regionDescriptions, 'lh-');
lhRegions = erase(regionDescriptions(lh), 'lh-');
rh = startsWith(regionDescriptions, 'rh-');
rhRegions = erase(regionDescriptions(rh), 'rh-');
overlap = intersect(lhRegions, rhRegions);
sub_ind_lh = ismember(lhRegions, overlap);
sub_ind_rh = ismember(rhRegions, overlap);
regions = lhRegions(sub_ind_lh);
sub = W([sub_ind_lh; sub_ind_rh], [sub_ind_lh; sub_ind_rh]);

N = length(regions)*2;

if includeInter

    % including interhemispheric connections
    sub_lh = sub(1:N/2,:);
    sub_rh(1:N/2,1:N/2) = sub(N/2+1:end, N/2+1:end); % intra first
    sub_rh(1:N/2,N/2+1:N) = sub(N/2+1:end, 1:N/2); % add inter

else
    
    % excluding interhemispheric connections
    sub_lh = sub(1:N/2,1:N/2);
    sub_rh = sub(N/2+1:end, N/2+1:end);
    
end

newN = size(sub_lh, 2);

% Optionally filter on strongest connections
if ~isempty(connLimit)
    
    % Pick the strongest connLimit connections from sub_lh and sub_rh
    % independently, for each region
    nconn = nan(newN, 1);
    for i = 1:newN
        
        [~, I] = sort(sub_lh(i,:), 'descend');
        sub_lh(i, I(connLimit+1:end)) = 0; % only keep the top connLimit strongest connections

        [~, I] = sort(sub_rh(i,:), 'descend');
        sub_rh(i, I(connLimit+1:end)) = 0; % only keep the top connLimit strongest connections
        
        nconn(i) = max([nnz(sub_lh(i,:)), nnz(sub_rh(i,:))]);
        
    end
    
else
   
    nconn = sum((sub_lh > 0)|(sub_rh > 0), 2);
    
end

% Binarize and calculate overlap
sub_lh = double(sub_lh > 0);
sub_rh = double(sub_rh > 0);

if presenceOnly
    
    % Overlap is number of present connections in LH and RH divided by
    % total number of present connections
    overlap = sub_lh == 1 & sub_rh == 1;
    overlap = sum(overlap, 2) ./ nconn;
    
else
    
    % Overlap is number of matching connections (both 1 or both 0) in LH
    % and RH divided by total number of possible connections
    overlap = sub_lh == sub_rh;
    overlap = sum(overlap, 2) ./ newN;
    

end
    
end