# matlab_vb_index
Vogt-Bailey index toolbox in Matlab

## Note on the MATLAB version of the VB Toolbox
While we recommend that users use the Python version of the toolbox
(https://github.com/VBIndex/py_vb_toolbox) which can be installed as
a standalone command-line application this readme file will help users
navigate through the MATLAB version of the code.

## Using the MATLAB VB Toolbox

The two starting functions are `vb_index.m` and `vb_cluster.m`

### vb_index.m

The function `vb_index.m` is used to perform the searchlight VB Index.

`RESULT = vb_index(GIFTI_SURF, DATA, NORM, CORT_INDEX, OUTPUT, nthreads)`

The input variables should already be loaded into MATLAB.

GIFTI_SURF is a variable that contains the GIFTI structure of the surface file 
containing the faces and vertices of a surface

DATA is an array that contains the timeseries that are associated with the GIFTI
surface (above)

CORT_INDEX is an array that contains a mask of all cortical vertices (excluding
the midline vertices)

OUTPUT is a string for the basename of your output file

nthreads is an option to use more than one thread for parallelising parts of 
the process. The default is using a single thread.

### vb_cluster.m

This function is used to perform the full brain gradient analysis and the 
clustered gradient and VB Index. Note that the full brain gradient analysis is
considered to be a special case of the clustered analysis.

`[ RESULT_EigenValues , RESULT_EigenVectors ]= vb_cluster(DATA, NORM, CORT_INDEX, CLUST_INDEX, OUTPUT, nthreads)`

Most variables from this function are identical to the one above. The principal
difference is the addition of the CLUST_INDEX. This file defines the clusters /
parcellation of the brain. Should one wish to carry out a whole brain gradient
analysis, CLUST_INDEX should simply be the same variable as CORT_INDEX (the brain
is considered to be one big parcel).
