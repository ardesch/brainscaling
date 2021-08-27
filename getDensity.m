function D = getDensity(connectivity)
% This function returns an array of connection densities based on a 4D connectivity
% matrix dataset

n = size(connectivity, 4);
D = zeros(n, 1);

for i = 1:n
    
    A = connectivity(:,:,1,i);
    D(i) = density_und(A);

end

end