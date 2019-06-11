function A = vb_create_affinity_matrix(neighborhood)


% Create a mean centered neighborhood
neighborhood_mean = mean(neighborhood,2);
neighborhood_mc = bsxfun(@minus, neighborhood, neighborhood_mean);

neighborhood_mc(neighborhood_mc==0)=eps;

% Normalise the mean centered neighborhood
neighborhood = sqrt(sum(neighborhood_mc.^2 , 2));
neighborhood = bsxfun(@rdivide, neighborhood_mc , neighborhood);

% Dot product
AFFINITY = neighborhood*neighborhood';

% avoid "correlation" values of above 1 - numerical error that
% leads to complex numbers in AFFINITY
AFFINITY(AFFINITY>.9999999)=1;

% "Linearlise" the AFFINITY
A = 1 - (acosd(AFFINITY) / 90);
A(A<=0) = eps;