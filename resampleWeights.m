function Wresampled = resampleWeights(W)
% This function resamples a weighted connectivity matrix W to a normal
% distribution by matching each ranked entry in W to the corresponding
% ranked entry in a normal distribution of the same number of elements. The
% normal distribution has mean 1 and standard deviation 0.2 by default.

% Summary statistics for Drandom:
originalWeights = nonzeros(squareform(W));
mu = 1; % arbitrary value
sigma = 0.2; % abitrary value
NoE = size(originalWeights,1);

Drandom = normrnd(mu,sigma,1,NoE);
DrandomSorted = sort(Drandom);

% Align both datasets based on rank order
V = squareform(W);
Vresampled = V.*0;
indices = 1:numel(V);
rankOrder = [V',indices'];

nonZeroIndices = indices(V~=0)';
nonZeroWeights = rankOrder(nonZeroIndices,1);

% Get a list of indices ordered by increasing weights
[~, sortedIndices] = sort(nonZeroWeights);

Vnonzerosresampled = zeros(size(DrandomSorted));
Vnonzerosresampled(sortedIndices) = DrandomSorted;
Vresampled(nonZeroIndices) = Vnonzerosresampled;

Wresampled = squareform(Vresampled); % convert vector back to matrix

end