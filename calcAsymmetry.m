function [asymmetry, regions] = calcAsymmetry(W, regionDescriptions, varargin)
% This function takes a connectivity matrix and regionDescriptions as
% input, and returns the connectivity asymmetry between each homologous 
% region pair in the left and right hemisphere (for the submatrix of 
% regions that appear in both hemispheres).

connLimit = size(W,1)/2;
includeInter = 1;

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
        otherwise
            error('option %s unknown', varargin{1});
    end

    varargin(1:2) = [];

end

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

    % Including interhemispheric connections
    sub_lh = sub(1:N/2,:);
    sub_rh(1:N/2,1:N/2) = sub(N/2+1:end, N/2+1:end); % intra first
    sub_rh(1:N/2,N/2+1:N) = sub(N/2+1:end, 1:N/2); % add inter

else
    
    % Excluding interhemispheric connections
    sub_lh = sub(1:N/2,1:N/2);
    sub_rh = sub(N/2+1:end, N/2+1:end);
    
end

% Calculate connectivity profile distances
D = calcProfile(sub_lh, sub_rh, 'connLimit', connLimit);

% Keep the diagonal (each region with homologous region in contralateral
% hemisphere)
asymmetry = diag(D);

end