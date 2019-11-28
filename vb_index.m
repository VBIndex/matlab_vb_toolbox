function RESULT = vb_index(GIFTI_SURF, DATA, NORM, CORT_INDEX, OUTPUT, nthreads)

% Script to produce the VB index
%
% GIFTI_SURF == cortical surface GIFTI object that matches the data
% DATA == the data that matches the GIFTI surface
% CORT_INDEX == an index of cortical vertices
% OUTPUT == a string for the basename of your output file, put in empty
% vector if no output to save "[]"
% nthreads == number of threads for this command (default == 0)
%
% RESULT = vb_index(GIFTI_SURF, DATA, NORM, CORT_INDEX, OUTPUT, nthreads)
%
% claude.bajada@um.edu.mt

% create default threads
if ~exist('nthreads' , 'var')
    nthreads = 0;
end

if ~(strcmpi(NORM , 'unnorm') || strcmpi(NORM , 'geig'))
    error('the value of NORM should be unnorm or geig')
end

% initialise RESULT vector for speed optimisation
RESULT = zeros(size(DATA,1),1);

parfor (i = 1:size(GIFTI_SURF.vertices,1), nthreads)
    
    % Find the neighborhood
    % finding the rows where our vertex of interest appears
    % collapse across rows using sum() then binarise using logical
    neighborhood_idx = logical(sum(GIFTI_SURF.faces == i,2));
    I = unique(GIFTI_SURF.faces(neighborhood_idx,:));
    neighborhood = DATA(I,:);
    
    % Create the affinity matrix
    AFFINITY = vb_create_affinity_matrix(neighborhood);
    
    % Find the algebraic connectivity (second smallest eigenvalue of laplacian)
    [~, ~, ~, ~, sortedEigenValues] = spectral_reorder(AFFINITY , NORM);
    
    % Normalise result for connectivity to be between zero and one.
%     normalisation_factor = sum(sortedEigenValues) / (numel(sortedEigenValues)-1);
    normalisation_factor = mean(sortedEigenValues(2:end));
    RESULT(i,1) = sortedEigenValues(2)/normalisation_factor;
    
end

% ensure that medial wall vertices are NaNs.
RESULT(logical(1-CORT_INDEX),:) = nan;

% save results if given an OUTPUT string
if ~isempty(OUTPUT)
    save(gifti(RESULT), [OUTPUT , '.vbi.shape.gii']);
end