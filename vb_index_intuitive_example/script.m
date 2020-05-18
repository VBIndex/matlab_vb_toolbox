%% Script for testing the VBIndex on "2D" colour photographs (or other images)

% The script will prompt to select an image. To recreate the image in the
% publication please select "mdina.jpg" however you may experiment with
% alternative images such as "sunset.jpg" or your own custom image.

% Part of matlab_vb_toolbox
% claude.bajada@um.edu.mt


%% Setup
clear
clc

% Load Image

[FILENAME, PATHNAME, FILTERINDEX] = uigetfile({'*'});

MY_IMAGE = imread([PATHNAME filesep FILENAME]);

EDGE_IM_unnorm = zeros(size(MY_IMAGE,1),size(MY_IMAGE,2));
EDGE_IM_geig = zeros(size(MY_IMAGE,1),size(MY_IMAGE,2));
EDGE_IM_rw = zeros(size(MY_IMAGE,1),size(MY_IMAGE,2));


%% Computations

% start from the second pixel to avoid edge effects - this won't be a
% problem on brain surfaces since they are continuous surfaces

f = waitbar(0,'Please wait...');

for i = 2:(size(MY_IMAGE,1)-1)
    
    waitbar(i / (size(MY_IMAGE,1)-1) ,f, 'Processing...');
    
    parfor (j = 2:(size(MY_IMAGE,2)-1) , 24)
    % NB: the paralelisation was not tested due to lack of license.
    % parfor works (non-paralel without the license and gives correct 
    % results. Kindly report any results with paralel toolbox so that
    % I can update this statement.
        
        neighborhood = zeros(9,3);
        
        % find the immediate neighborhood
        neighborhood(1,:) = MY_IMAGE(i-1,j-1,:); % Upper left pixel.
        neighborhood(2,:) = MY_IMAGE(i-1,j,:); % Upper middle pixel.
        neighborhood(3,:) = MY_IMAGE(i-1,j+1,:); % Upper right pixel.
        neighborhood(4,:) = MY_IMAGE(i,j-1,:); % left middle pixel.
        neighborhood(5,:) = MY_IMAGE(i,j,:); % pixel of interest (middle)
        neighborhood(6,:) = MY_IMAGE(i,j+1,:); % right middle pixel.
        neighborhood(7,:) = MY_IMAGE(i+1,j+1,:); % Lowerleft pixel.
        neighborhood(8,:) = MY_IMAGE(i+1,j,:); % lower middle pixel.
        neighborhood(9,:) = MY_IMAGE(i+1,j-1,:); % Lower left pixel.
        
        
        % ensure that the data is stored as double, othewise we will run
        % into problems later
        neighborhood = double(neighborhood);
        neighborhood(neighborhood==0) = eps;
        neighborhood = rgb2hsv(neighborhood/255); % change from rgb to hsv since it is more intuitive to find edges of this representation
        
        % Create the affinity matrix (between rgb vectors)
        AFFINITY = vb_create_affinity_matrix(neighborhood);
        
        % Solve the generalised eigenvalue problem for the
        % laplacian of the affinity matrix
        
        % note that these conditionals ensure that if there are ANY
        % warnings generated in the code, the pixel value for that
        % computation will be a NaN.
        msg=lastwarn('');
        [~, ~, ~, ~, sortedEigenValues_unnorm] = spectral_reorder(AFFINITY , 'unnorm');
        if ~isempty(msg)
            sortedEigenValues_unnorm = [nan nan];
        end
        
        msg=lastwarn('');
        [~, ~, ~, ~, sortedEigenValues_geig] = spectral_reorder(AFFINITY , 'geig');
        if ~isempty(msg)
            sortedEigenValues_geig = [nan nan];
        end
        
        msg=lastwarn('');
        [~, ~, ~, ~, sortedEigenValues_rw] = spectral_reorder(AFFINITY , 'rw');
        if ~isempty(msg)
            sortedEigenValues_geig = [nan nan];
        end
        
        % Normalise result for connectivity to be between zero and one.
        normalisation_factor_unnorm = mean(sortedEigenValues_unnorm(2:end));
        normalisation_factor_geig = mean(sortedEigenValues_geig(2:end));
        normalisation_factor_rw = mean(sortedEigenValues_rw(2:end));
        
        % Save the new pixel - second smallest eigenvalue (lin cor)
        new_pixel_unnorm = sortedEigenValues_unnorm(2)/normalisation_factor_unnorm;
        new_pixel_geig = sortedEigenValues_geig(2)/normalisation_factor_geig;
        new_pixel_rw = sortedEigenValues_rw(2)/normalisation_factor_rw;
        
        EDGE_IM_unnorm(i,j) = new_pixel_unnorm;
        EDGE_IM_geig(i,j) = new_pixel_geig;
        EDGE_IM_rw(i,j) = new_pixel_rw;
        
    end
end

% removing the edge
EDGE_IM_unnorm_no_boarder = EDGE_IM_unnorm(2:end-1,2:end-1);
EDGE_IM_geig_no_boarder = EDGE_IM_geig(2:end-1,2:end-1);
EDGE_IM_rw_no_boarder = EDGE_IM_rw(2:end-1,2:end-1);

close(f)

% display the images - NOTE I adjusted the colour bar in the images I sent
% you to better distinguish the boarders.

%% Visualisation

% note that the published image underwent some minor post processing to
% improve the final resolution, background colour and spacing.

figure;

subplot(2,2,1); image(MY_IMAGE); axis equal, axis off; title('Original Image')
subplot(2,2,2); imagesc(EDGE_IM_unnorm_no_boarder); axis equal; colormap gray, axis off; title('Unnormalised');
caxis([0 1])
subplot(2,2,3); imagesc(EDGE_IM_geig_no_boarder); axis equal; colormap gray, axis off; title('Geig');
caxis([0 1])
subplot(2,2,4); imagesc(EDGE_IM_rw_no_boarder); axis equal; colormap gray, axis off; title('rw');
caxis([0 1])