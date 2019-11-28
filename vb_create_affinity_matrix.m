function A = vb_create_affinity_matrix(neighborhood)


% Create a mean centered neighborhood
neighborhood_mean = mean(neighborhood,2);
neighborhood_mc = bsxfun(@minus, neighborhood, neighborhood_mean);

neighborhood_mc(neighborhood_mc==0)=eps;

% Normalise the mean centered neighborhood
neighborhood_w = sqrt(sum(neighborhood_mc.^2 , 2));
neighborhood = bsxfun(@rdivide, neighborhood_mc , neighborhood_w);

% Dot product
AFFINITY = neighborhood*neighborhood';

% avoid "correlation" values of above 1 - numerical error that
% leads to complex numbers in AFFINITY
AFFINITY(AFFINITY>.9999999)=1;

% "Linearlise" the AFFINITY 
% ensure positive correlations are between 0 to 1
% what i am doing here is to change cosines to angles but to ensure that
% this remains a similarity rather than dissimilarity I am treating the
% values as the sine rather than cosine. I.e. identical signals will be 90
% rather than 0. 90/90 == 1.
A = (asind(AFFINITY) / 90);

% remove negative correlations
A(A<=0) = eps;