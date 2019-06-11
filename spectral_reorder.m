function [sortedB, p, v2, v3D, sortedEigenValues , sortedEigenVectors] = spectral_reorder(B,method)
  % Function that spectrally reorders a similarity matrix and outputs the
  % spectrum of eigenvalues and eigenvectors
  %
  % Function supplemental to the text of: "The Spectral Transformation in
  % Neuroimaging: A review and MATLAB tool" Bajada et al (in prep)
  %
  % Matlab implementation partially based on an algorithm published in supplemetary text of:
  % Johansen-Berg et al. (2004) PNAS 10.1073/pnas.0403743101.
  % For in depth mathematical background, see:
  % Barnard, S.T., Pothen, A., & Simon, H.D. (1995) Numer. Linear Algebra Appl. 2,317-334
  % For a tutorial on the different approaches please see:
  % von Luxburg, U. (2007) Statistics and Computing 17 (4), 2007
  %
  % Here we implement 4 different methods:
  % method = 'geig' Based on the generalised spectral decomposition of the
  % un-normalised Laplacian (default approach)
  % method = 'sym' Based on the spectral decomposition of the normalised
  % symmetric Laplacian - results should match geig
  % method = 'rw' Based on the spectral decomposition of the normalised
  % non-symmetric (random-walk) Laplacian
  % method = 'unnorm' Based on the spectral decomposition of the
  % un-normalised Laplacian
  %
  %
  % Function takes the form:
  % [sortedB, p, v2, v3D, sortedEigenValues , sortedEigenVectors] = spectral_reorder(B,method);
  %
  % B = similarity matrix (accepts values between -1 and 1)
  % sortedB = spectrally reordered matrix
  % p = the reordering vector that allows you to reorder your subjects in the
  % same way as your matrix (sortedB)
  % v2 = Fiedler vector (eigen vector corresponding to the second smallest
  % eigenvalue)
  % v3D = First 3 smaller eigenvectors (the first one is the Fiedler vector)
  % sortedEigenValues = the complete set of sorted eigenvalues
  % sortedEigenVectors = the complete set of sorted eigenvectors
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %
  % If you only assign one output
  % (i.e. of the form sortedB = extract_spectral_reorder(B), only the spectrally reordered matrix
  % - not the reordering vector or the other outputs - will be produced)
  %
  % Created by:
  % Claude J. Bajada - University of Manchester
  % Nelson Trujillo-Barreto - University of Manchester
  %
  % Edited by:
  % Owen Falzon - University of Malta
  %
  % Neuroscience and Aphasia Research Unit (NARU)
  % University of Manchester
  % 2014
  %
  % Updates:
  % 28/07/2016: changed code to include the generalised eigenvalue
  % problem
  % 04/08/2016: edited to incorporate 3 different methods (see above)
  % 09/08/2016: edited to include unnormalised version

  %% Initialisation
  if nargin < 2
    method = 'geig';
  end
  szB=size(B);
  if szB(1) ~= szB(2)
    error('The similarity matrix must be square');
  end
  % remove any NaNs or Infs from the matrix
  if sum(sum((isnan(B)))) > 0
    B(isnan(B)) = 0;
  end

  if sum(sum((isinf(B)))) > 0
    B(isinf(B)) = 0;
  end

  if min(B(:)) < 0 && min(B(:)) >= -1
    C = B + 1; % compute C which removes the negative values from the correlation matrix (ie values from 0 - 2; 0 = anticorrelated; 2 = correlated)
    clc
    disp('The value 1 is being added to your similarity matrix to ensure positivity.')
    disp('This may cause problems with interpretation, consider inputting a positive matrix');
    finish = 100;
    my_input = input('press enter to continue with 1 added to your matrix or type FINISH to stop the programme: ');
    if ~isempty(my_input) && my_input == finish
      sortedB = []; p = []; v2 = []; v3D = []; sortedEigenValues = []; sortedEigenVectors = [];
      return
    end

  elseif min(B(:)) < -1
    error('This function only accepts matrices with a maximum negative value of -1');
  else
    C = B; % if there are no negative values (eg in euclidian distance matrices) do not add one to Matrix B
  end

  %% Spectral decomposition

  % create the laplacian matrix (Q).
  % For all non diagonal elements, populate matrix Q with the negative
  % value of matrix C
  % For all the diagonal element, sum across the rows (excluding the
  % diagonal element) and populate the diagonal of Q with that sum

  triuC = triu(C,1); % Extract upper triangular elements;
  C = triuC + triuC'; % Reconstruct a symmetric weighted adjacency matrix eliminating possible small errors in off-diagonal elements
  D = diag(sum(C));% Compute the Degree matrix

  Q = D - C; %Compute un-normalised Laplacian

  % note the forcing of real values below is due to numerical errors
  % occasionaly creating complex eigenvalues and eigenvectors, this is
  % particularly problematic with the non-symmetric random walk normalised
  % laplacian

  switch lower(method)
    case 'geig'
      % Method using generalised spectral decomposition of the
      % un-normalised Laplacian (see Shi and Malik, 2000)
      [V,s] = eig(Q,D);
      [sortedEigenValues,b] = sort(diag(real(s)),'ascend'); % Sorteigenvalues
      if sum(isnan(sortedEigenValues)) > 0
          error('NaNs present in solution')
      end
      v2 = real(V(:,b(2))); % Get Fiedler vector. Here renormalisation is not necessary
      [~,p] = sort(v2); % Get reordering operator.
      sortedB = B(p,p); % Apply reordering.
      v3D = real(V(:,b(2:4))) ; % Extract first three eigenvectors (non-zero)
      sortedEigenVectors = real(V(:,b(1:length(b)))); % Extract the full set of eigenvectors

    case 'sym'
      % Method using eigen decomposition of Symmetric Normalised Laplacian
      % Note results should be exactly the same as geig
      t = sqrt(D);
      L = (t \ Q) / t; % Compute the normalised Laplacian Note, use of "/" and "\" is faster and more accurate than multiplying by the inverse
      L = force_symmetric(L); %force exact symmetry
      [V,s] = eig(L);
      [sortedEigenValues,b] = sort(diag(real(s)),'ascend'); % Sort eigenvalues
      if sum(isnan(sortedEigenValues)) > 0
          error('NaNs present in solution')
      end
      v2 = real(t \ V(:,b(2))); % Get Fiedler vector and renormalise
      [~,p] = sort(v2); % Get reordering operator.
      sortedB = B(p,p); % Apply reordering.
      v3D = t \ V(:,b(2:4)); % Extract first three eigenvectors and re-normalise to "stretch" the axis in the eigenspace
      sortedEigenVectors = real(t \ V(:,b(1:length(b)))); % Extract the full set of eigenvectors and re-normalise

    case 'rw'
      % Method using eigen decomposition of Random Walk Normalised Laplacian
      % This method has not been rigorously tested yet
      t = D;
      L = t\Q; % this laplacian is asymmetric
      [V,s] = eig(L);
      [sortedEigenValues,b] = sort(diag(real(s)),'ascend'); % Sort eigenvalues
      if sum(isnan(sortedEigenValues)) > 0
          error('NaNs present in solution')
      end
      v2 = real(V(:,b(2))); % Get Fiedler vector. Here renormalisation is not necessary
      [~,p] = sort(v2); % Get reordering operator.
      sortedB = B(p,p); % Apply reordering.
      v3D = real(V(:,b(2:4))) ; % Extract first three eigenvectors
      sortedEigenVectors = real(V(:,b(1:length(b)))); % Extract the full set of eigenvectors

    case 'unnorm'
      % Method using spectral decomposition of the unnormalised Laplacian
      [V,s] = eig(Q);
      [sortedEigenValues,b] = sort(diag(real(s)),'ascend'); % Sort eigenvalues
      v2 = real(V(:,b(2))); % Get Fiedler vector. Here renormalisation is not necessary
      [~,p] = sort(v2); % Get reordering operator.
      sortedB = B(p,p); % Apply reordering.
      v3D = real(V(:,b(2:4))) ; % Extract first three eigenvectors
      sortedEigenVectors = real(V(:,b(1:length(b)))); % Extract the full set of eigenvectors

    otherwise
      error('This method is not allowed: the input for method must be: "sym" , "rw" , "geig" or "unnorm"')

  end

end

  %% internal function for forced symmetry
function M = force_symmetric(M)
  DiagM = diag(diag(M));
  triuM = triu(M,1); % Extract upper triangular elements;

  % Reconstruct a symmetric L eliminating possible small errors in off-diagonal elements
  M = triuM + DiagM + triuM'; 
end
