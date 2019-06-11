function [CORT_INDEX , DATA_norm] = vb_format_data(CIFTI_DATA, SIDE, GIFTI_ATLAS, CIFTI_ATLAS)

% load gifti atlas file - this file indexes the vertices that contain
% cortical data
CORT_INDEX = GIFTI_ATLAS.cdata;

% isolate cortical data from cifti atlas
CIFTI_ATLAS = logical(CIFTI_ATLAS.cdata == 1);

% Only include cortical vertices
CIFTI_DATA = CIFTI_DATA.cdata(CIFTI_ATLAS,:);

% Initialise DATA matrix to ensure correct size
DATA = zeros(size(CORT_INDEX,1) , size(CIFTI_DATA,2));

% Extract left or right sided data
if strcmpi(SIDE , 'left')    
    DATA(logical(CORT_INDEX),:) = CIFTI_DATA(1:sum(CORT_INDEX),:);
elseif strcmpi(SIDE , 'right') 
    DATA(logical(CORT_INDEX),:) = CIFTI_DATA(end-sum(CORT_INDEX)+1:end,:);
else
    error('the variable SIDE must be left or right')
end

% Create a mean centered DATA
DATA_mean = mean(DATA,2);
DATA_mc = bsxfun(@minus, DATA, DATA_mean);

DATA_mc(DATA_mc==0)=eps;

% Normalise the mean centered data
DATA_norm = sqrt(sum(DATA_mc.^2 , 2));
DATA_norm = bsxfun(@rdivide, DATA_mc , DATA_norm);
