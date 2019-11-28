function [ RESULT_EigenValues , RESULT_EigenVectors ]= vb_cluster(DATA, NORM, CORT_INDEX, CLUST_INDEX, OUTPUT, nthreads)

% claude.bajada@um.edu.mt

% create default threads
if ~exist('nthreads' , 'var')
    nthreads = 0;
end

if ~(strcmpi(NORM , 'unnorm') || strcmpi(NORM , 'geig'))
    error('the value of NORM should be unnorm or geig')
end

RESULT_EigenValues = zeros(size(DATA,1),1);
RESULT_EigenVectors = zeros(size(DATA,1),max(CLUST_INDEX));

parfor (i = 1:max(CLUST_INDEX), nthreads)
   
    % Find the neighborhood (i.e. cluseters)
    neighborhood = DATA(CLUST_INDEX==i , :); 
    
    % Create the affinity matrix
    AFFINITY = vb_create_affinity_matrix(neighborhood);
    
    % Find the algebraic connectivity (second smallest eigenvalue of laplacian)
    [~, ~, ~, ~, sortedEigenValues, sortedEigenVectors] = spectral_reorder(AFFINITY, NORM);
    
    % Normalise eigenvalues for value to be between zero and one.
    normalisation_factor = sum(sortedEigenValues) / (numel(sortedEigenValues)-1);
    EigenValues(i) = sortedEigenValues(2)/normalisation_factor;
    EigenVectors{i} = sortedEigenVectors;
    
end

for i = 1:max(CLUST_INDEX)
    
    RESULT_EigenValues(CLUST_INDEX==i) = EigenValues(i);
    RESULT_EigenVectors(CLUST_INDEX==i,i) = EigenVectors{i}(:,2);
    
end
    
    
% ensure that medial wall vertices are NaNs.
RESULT_EigenValues(logical(1-CORT_INDEX),:) = nan;
RESULT_EigenVectors(logical(1-CORT_INDEX),:) = nan;

% save results if given an OUTPUT string
if ~isempty(OUTPUT)
    save(gifti(RESULT_EigenValues), [OUTPUT , '.vb-cluster.value.shape.gii']);
    save(gifti(RESULT_EigenVectors), [OUTPUT , '.vb-cluster-vector.shape.gii']);
end