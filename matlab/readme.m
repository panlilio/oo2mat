%% OpenOrganelle MATLAB wrappers - README
% Guide to using MATLAB to access data from Janelia-COSEM's OpenOrganelle 
% Project as stored on their AWS S3 bucket: s3://janelia-cosem/.
%
% Author: Mia Panlilio (2021-12-01)
%
% Consult class and method helpstrings for more information:
% 
%   >> help openorganelle_filebrowser
%   >> help openorganelle_dataset
%   >> help openorganelle_dataset.get_organelle

%% Requirements
% The two MATLAB classes described below retrieve OpenOrganelle data by wrapping 
% the oo2mat Python module. Its main dependencies are:
%
% * numpy
% * s3fs
% * scikit-image
% * fibsem-tools
%
% Note that oo2mat relies on Janelia's fibsem-tools package (currently 
% available through pip installation), which requires an old version of the 
% hdf5 and h5py packages to read N5 data. Due to version conflicts in these 
% packages between MATLAB and fibsem-tools, _MATLAB must run Python out-of-process_. 
% To check this, run:
%%
%   >> pe = pyenv;
%   >> disp(pe.ExecutionMode)
%%
% If the output is not 'OutOfProcess', you must restart MATLAB and then change 
% the default Python environment settings via:
%%
%   >> pyenv('ExecutionMode','OutOfProcess')

%% Class: openorganelle_filebrowser
% The *openorganelle_filebrowser* class allows users to query general 
% information about the OpenOrganelle datasets, e.g. the names of all available 
% datasets and their associated segmentation masks. It should be initialized 
% without any input arguments, as the default parameters defining the locations 
% of the AWS S3 bucket, Python modules, and organelle catalogs should be stable 
% at least for use on the St. Jude network. 
fb = openorganelle_filebrowser;
disp(fb)

%% Example: list datasets stored on AWS
fb.ls_datasets;

%% Example: list all segmentation (label) layers available for jrc_macrophage-2
% The following command displays a list of all segmentation types stored in the
% jrc_hela-2 N5. In principle, any of these can be catalogued into instance
% bounding boxes (see *openorganelle_dataset.do_catalog* below).
fb.ls_segmentations('jrc_macrophage-2');

%% Class: openorganelle_dataset
% The *openorganelle_dataset* class is a subclass of openorganelle_filebrowser.
% This class is used to retrieve organelle volumes, and general image volumes, 
% belonging to a specified biological cell. The class should be initialized with 
% the cell of interest's COSEM name as the first input argument. 
%
% Let's say we want to probe jrc_hela-2:
ds = openorganelle_dataset('jrc_hela-2');
disp(ds)

%%
% In addition to returning arbitrary volumes from both FIB-SEM (grayscale) and 
% segmentation (binary label) layers, this class uses and can be used to 
% generate bounding box catalogs that retrieve sparse matrices with masks 
% corresponding to specific instances of a given organelle.
%
% Additionally, openorganelle_dataset inherits the usual properties and methods 
% of its superclass with some minor modifications to make certain queries 
% specific to the dataset of interest, e.g. openorganelle_dataset.ls_segmentations
% does not take any input arguments and simply returns the list of segmentation 
% layers associated with the initialized dataset.

%% Example: retrieve the mask for a specific mitochondrion instance
% Let's say we want the binary mask of the 25-th catalogued mitochrondrion.
%%
S = ds.get_organelle('mito_seg',25);
%%
disp(S)
%%
% Notice that the organelle_idx field here is returned from Python and
% therefore 0-indexed (24), rather than MATLAB's 1-index (25). 

%% Example: retrieve multiple masks for specific mitochondria instances
% If multiple masks are needed, it is generally much faster to submit multiple 
% indices to get_organelle at once rather than call the function
% repeatedly. Say we want the mitochrondria catalogued at
% [1,2,3,10,19,25,60] with their FIB-SEM grayscale values:
S = ds.get_organelle('mito_seg',[1:3,10,19,25,60],'with_fibsem',true);
%%
disp(S)
%%
% The image/mask volumes are stored in a cell array in the order that they
% were submitted. Then the mask and FIB-SEM data for the 25-th mitochrondrion 
% are at the 6-th index of the returned arrays. However, when 'with_fibsem' is 
% set to true, the resultant volumes are returned at a lower resolution:
disp(size(S.data_mask{6}))
%%
% This is because the label volumes were upscaled at Janelia (2x in all 
% dimensions) for segmentation purposes and get_organelle returns
% volumes at the highest possible resolution while retaining a consistent
% coordinate system with FIB-SEM data if needed.

%% Example: retrieve an arbitrary FIB-SEM volume
% Note that all bounding boxes stored in either a catalog or returned by 
% get_organelle or get_volume functions obey the scikit-image convention of
% [x_min, y_min, z_min, x_max, y_max, z_max]. Say we want to look at the 
% FIB-SEM volume defined by the bounding box [100,200,300,250,400,350]:
V = ds.get_volume('fibsem-uint8',[100,200,300,250,400,350]);
whos V
%%
disp(V(1:10,1:10))
%%
% Again, by default the highest resolution is returned. Resolution levels
% can be explicitly requested, with input bounding boxes assumed to be on the 
% coordinate system defined by the requested resolution level. The parameter
% 'resolution' should be an integer between 0 and 4 (or 5) in order of
% decreasing resolution.
V = ds.get_volume('fibsem-uint8',[100,200,300,250,400,350],'resolution',0);
disp(V(1:10,1:10))

%% Example: retrieve an arbitrary label volume
V = ds.get_volume('lyso_seg',[100,200,300,250,400,350]);
whos V

%% Example: run catalog of centrioles
% Let's say we want to look at instances of segmented centrioles using 
% get_organelle but the label volume has not been catalogued yet. To
% catalog:
[catalog_filename,J] = ds.do_catalog('cent_seg');
%%
% By default, this function will not run if the catalog already exists but
% it can be forced to perform an overwrite (see help). Cataloguing returns the 
% location of the resultant catalog:
disp(catalog_filename)
%%
% and a struct corresponding to the JSON:
disp(J)

