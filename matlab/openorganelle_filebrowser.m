classdef openorganelle_filebrowser < handle
    % OPENORGANELLE_FILEBROWSER     Handle class to browse datasets in the
    % COSEM-Janelia OpenOrganelle project.
    %
    % OPENORGANELLE_FILEBROWSER properties:
    %   s3bucket - S3 address to Janelia-COSEM bucket
    %   bounding_boxes_path - path to organelle instance bounding box catalogs
    %   oo2mat_path - path to directory containing oo2mat Python package
    %   oo2mat_mod - oo2mat Python module
    %
    % OPENORGANELLE_FILEBROWSER methods:
    %   openorganelle_filebrowser - class constructor
    %   ls_datasets - list datasets available in s3bucket
    %   ls_segmentations - list segmentations (labels) available for a
    %                      given dataset
    %   ls_fibsem - list FIB-SEM image volumes available for a given
    %               dataset
    
    properties (SetAccess = immutable)
        % s3bucket - property of openorganelle_filebrowser
        % Address of the AWS S3 bucket containing OpenOrganelle N5
        % datasets. Default is 's3://janelia-cosem'
        s3bucket = 's3://janelia-cosem/'
        % bounding_boxes_path - property of openorganelle_filebrowser
        % Path to organelle catalogs. Default is
        % '/research/sharedresources/cbi/common/mia/OpenOrganelle/lib/organelle_bounding_boxes/'
        bounding_boxes_path = '/research/sharedresources/cbi/common/mia/OpenOrganelle/lib/organelle_bounding_boxes/'
        % oo2mat_path - property of openorganelle_filebrowser
        % Path to directory containing oo2mat, the underlying Python package 
        % used for data retrieval and AWS queries. Default is
        % '/research/sharedresources/cbi/common/mia/OpenOrganelle/python/'
        oo2mat_path = '/research/sharedresources/cbi/common/mia/OpenOrganelle/python/'
        % oo2mat_mod - property of openorganelle_filebrowser
        % oo2mat Python module loaded from its location in oo2mat_path.
        oo2mat_mod = []
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                 CLASS CONSTRUCTOR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = openorganelle_filebrowser(varargin)
        %obj = openorganelle_filebrowser(Name,Value)
        %Class for browsing data directories and files associated with the 
        %OpenOrganelle project.
        %
        %INPUTS: Optional NAME, VALUE pairs
        %
        %'s3bucket'         path to Janelia COSEM's AWS bucket. 
        %                   Default = 's3://janelia-cosem/'. 
        %                   It is unlikely that users will need to change
        %                   this parameter anytime in the near future.
        %
        %'bounding_boxes_path' path to directory containing catalogued
        %                   bounding boxes corresponding to different
        %                   organelle instances for different datasets.
        %                   Default = '/research/sharedresources/cbi/common/mia/OpenOrganelle/lib/organelle_bounding_boxes/'.        
        %
        %'oo2mat_path'      Path to directory containing oo2mat package 
        %                   that uses Janelia's fibsem_tools package to 
        %                   retrieve N5 data from AWS buckets. 
        %                   Default = '/research/sharedresources/cbi/common/mia/MATLAB/OpenOrganelle/'. 
        %
        %Author: Mia Panlilio (2021-11-11)
        %
        %***************
        %*N.B. PLEASE ENSURE THAT MATLAB IS RUNNING PYTHON OUT-OF-PROCESS. 
        %You can check this with pyenv, under its ExecutionMode property. 
        %If it is set to "InProcess", you must restart MATLAB, then run 
        %pyenv("ExecutionMode","OutOfProcess").
        %
        %This class relies on Janelia's fibsem_tools, which in turn relies
        %on old versions of the h5py and hdf5 packages to access N5 volume
        %data. MATLAB on the other hand depends on the newer versions. If
        %you try to run this wrapper in-process, MATLAB will likely 
        %terminate your session because of the version conflict or throw 
        %an error. 
        %
        
        p = inputParser;
        p.addParameter('s3bucket',obj.s3bucket,@ischar);
        p.addParameter('bounding_boxes_path',obj.bounding_boxes_path,@ischar);
        p.addParameter('oo2mat_path',obj.oo2mat_path,@ischar);
        p.parse(varargin{:});

        flds = fields(p.Results);
        for f = 1:numel(flds)
            obj.(flds{f}) = p.Results.(flds{f});
        end        
        
        %Add user module to python path
        pe = pyenv;
        if strcmpi(pe.ExecutionMode,'InProcess')
            error('ERROR (openorganelle_filebrowser): Python must be run out-of-process. Please restart MATLAB and run pyenv(''ExecutionMode'',''OutOfProcess'')')
        end
        addpath(obj.oo2mat_path)
        insert(py.sys.path, int32(0), obj.oo2mat_path);
        obj.oo2mat_mod = py.importlib.import_module('oo2mat');
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%         OO2MAT PYTHON WRAPPERS: DATASET INFO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function L = ls_datasets(obj)
        %L = ls_datasets(obj)
        %Displays a list of all datasets in the Janelia-COSEM AWS bucket.
        %If output is requested, L will be a cell array of strings
        %corresponding to dataset names and nothing will be printed to the
        %command line.
        %
        
        L = obj.oo2mat_mod.ls_datasets(pyargs('s3path',obj.s3bucket));
        L = cellfun(@char,cell(L),'uniformoutput',false);
                
        if nargout==0 
            fprintf('\n----- DATASETS in %s -----\n',obj.s3bucket);
            fprintf('%s\n',L{:});
            fprintf('\n')
            clear L; 
        end
                
        end
        
        
        function L = ls_segmentations(obj,dataset_name)
        %L = ls_segmentations(obj,dataset_name)
        %Displays a list of all segmentations available for a given
        %dataset. If output requested, L will be a cell array of strings
        %corresponding to the segmentation names and nothing will be
        %printed to the command line.
        %
        %INPUT
        %
        %dataset_name   string corresponding to the name of the 
        %               OpenOrganelle dataset of interest, 
        %                   e.g. 'jrc_hela-2'.
        %               
        %               Run .ls_datasets to see the full list of datasets
        %               available in the AWS bucket.
        %
        
        L = obj.oo2mat_mod.ls_segmentations(pyargs('dataset_name',dataset_name,'s3path',obj.s3bucket));
        L = cellfun(@char,cell(L),'uniformoutput',false);
        
        if nargout==0 
            fprintf('\n----- SEGMENTATIONS for %s -----\n',dataset_name);
        	fprintf('%s\n',L{:});
            fprintf('\n')
            clear L;            
        end
        
        end
        
        
        function L = ls_fibsem(obj,dataset_name)
        %L = ls_fibsem(obj,dataset_name)
        %Displays a list of all fibsem volumes available for a given
        %dataset. If output requested, L will be a cell array of strings
        %corresponding to the fibsem dataset names and nothing will be
        %printed to the command line.
        
        L = obj.oo2mat_mod.ls_fibsem(pyargs('dataset_name',dataset_name,'s3path',obj.s3bucket));
        L = cellfun(@char,cell(L),'uniformoutput',false);
        fprintf('\n----- FIBSEM VOLUMES for %s -----\n',dataset_name);
        fprintf('%s\n',L{:});
        fprintf('\n')
        
        if nargout==0
            fprintf('\n----- FIBSEM VOLUMES for %s -----\n',dataset_name);
            fprintf('%s\n',L{:});
            fprintf('\n')
            clear L; 
        end
        
        end
        
        
        
    end
end