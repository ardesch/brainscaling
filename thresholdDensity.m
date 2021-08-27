function [T, minP] = thresholdDensity(A, P, thresh, minWeight)
% Extract a certain number of most prevalent connections of an
% individual subjects' connectivity matrix resulting in a thresholded
% matrix of predefined density. If P = A then instead of the most 
% prevalent connections, the strongest connections in weighted connectivity
% matrix A are extracted to arrive at a density of thresh. If P = A > 0
% then connections are randomly extracted to arrive at a density of thresh.
%
% Input:        A,          weighted or binary connectivity matrix
%               P,          group-level prevalence matrix
%               thresh,     proportional threshold (e.g. 0.4 for top 40%)
%               minWeight,  minimum weight to be included in case of
%                           weighted matrix
% Output:       T,          thresholded matrix
%               minP,       value of the lowest-prevalence connection in T


if thresh >= 1 || thresh <= 0
    error('Density treshold (thresh) should be between 0-1.')
end

if isa(A, 'logical')
    minWeight = 1;
elseif nargin == 3 && isa(A, 'double')
    minWeight = min(A(A > 0)); % minimum positive value in weighted matrix
end

B = A >= minWeight; % binarize matrix in case A is weighted
Pi = P .* B; % individual prevalence matrix
PiV = squareform(Pi);
n = size(B, 1); % number of nodes
V = squareform(zeros(n, n)); % empty vector

total = sum(squareform(ones(size(A)) - eye(size(A))));
conn = floor(total * thresh); % number of connections to get the required density
[p, I] = sort(PiV, 'descend');
minP = min(p(1:conn)); % minimum prevalence to be included

VA = squareform(A);
I = I(VA(I) >= minWeight); % only include entries in I of which VA is at least minWeight
PiV = PiV(I); % idem for PiV

V(I(PiV > minP)) = 1; % set all indices with PiV > minP to 1
ties = I(PiV == minP); % get all indices with PiV == minP

% Randomly draw from this pool to fill the remaining 'open spots'
remaining = conn - sum(V); % number of remaining spots
selection = randsample(ties, remaining);
V(selection) = 1; % set all indices in the selected sample to 1

T = squareform(V) > 0; % transform back to logical binary matrix
T = A .* double(T);

end