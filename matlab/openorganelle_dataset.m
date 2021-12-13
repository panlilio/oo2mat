classdef openorganelle_dataset < openorganelle_filebrowser
    % OPENORGANELLE_DATASET  Class to handle individual cell image volume 
    % datasets belonging to Janelia-COSEM's OpenOrganelle project. A subclass of 
    % the openorganelle_filebrowser class. Run "help openorganelle_filebrowser"
    % to see inherited properties and methods. 
    %
    % OPENORGANELLE_DATASET properties:
    %   dataset_name - name of the cell following COSEM conventions, e.g.
    %                   'jrc_hela-2'
    %
    % OPENORGANELLE_DATASET methods:
    %   openorganelle_dataset - class constructor
    %   ls_segmentations - list segmentation layers 
    %   ls_fibsem - list FIB-SEM image volumes
    %   ls_organelles - list catalogued organelles
    %   n_organelle - number of organelles of a given type, i.e. instance count
    %   get_bb_json - load JSON catalog with organelle bounding boxes
    %   get_organelle - retrieve mask volume corresponding to specific
    %                   organelle instance(s)
    %   get_volume - retrieve image volume defined by a bounding box
    %   do_catalog - catalog a specific organelle type, e.g. for instance
    %                retrieval later
    %
    
    properties
        % dataset_name - property of openorganelle_dataset
        % Cell name according to COSEM convention. To list all datasets 
        % available in the S3 bucket, run:
        %   >> fb = openorganelle_filebrowser;
        %   >> fb.ls_datasets     
        dataset_name = 'jrc_hela-2'        
    end
        
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                   CLASS CONSTRUCTOR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = openorganelle_dataset(dataset_name,varargin)
        %OBJ = OPENORGANELLE_DATASET(dataset_name,Name,Value)
        %Constructor for the handle class for OpenOrganelle datasets. 
        %
        %INPUTS
        %
        %dataset_name       N5 dataset name. Default = 'jrc_hela-2'.
        %
        %Optional NAME,VALUE arguments
        %
        %'s3bucket'         path to Janelia COSEM's AWS bucket. Default = 
        %                   's3://janelia-cosem/'. It is unlikely that users 
        %                   will need to change this parameter anytime in the 
        %                   near future.
        %
        %'bounding_boxes_path' path to directory containing catalogued bounding 
        %                   boxes corresponding to different organelle instances
        %                   for different datasets. Default = 
        %                   '/research/sharedresources/cbi/common/mia/OpenOrganelle/lib/organelle_bounding_boxes/'.        
        %
        %'oo2mat_path'      Path to directory containing oo2mat package that 
        %                   uses Janelia's fibsem_tools package to retrieve N5 
        %                   data from AWS buckets. Default = 
        %                   '/research/sharedresources/cbi/common/mia/MATLAB/OpenOrganelle/'. 
        %
        %
        %Author: Mia Panlilio (2021-11-11)
        %
        %***************
        %*N.B. PLEASE ENSURE THAT MATLAB IS RUNNING PYTHON OUT-OF-PROCESS. 
        %You can check this with pyenv, under its ExecutionMode property. If it 
        %is set to "InProcess", you must restart MATLAB, then run 
        %pyenv("ExecutionMode","OutOfProcess").
        %
        %This class relies on Janelia's fibsem_tools, which in turn relies on 
        %old versions of the h5py and hdf5 packages to access N5 volume data. 
        %MATLAB on the other hand depends on the newer versions. If you try to 
        %run this wrapper in-process, MATLAB will likely terminate your session 
        %because of the version conflict or throw an error. 
        %
        
        %Call openorganelle_filebrowser constructor
        obj@openorganelle_filebrowser(varargin{:})
        if nargin~=0
            obj.dataset_name = dataset_name;
        end
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                   DATASET INFO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function L = ls_segmentations(obj)
        %L = ls_segmentations(obj)
        %Lists available segmentations in AWS bucket for the current dataset. If
        %output argument is requested, L will be a cell array of strings with 
        %the segmentation names and nothing will be printed to the command
        %line.
                       
        if nargout==0
            ls_segmentations@openorganelle_filebrowser(obj,obj.dataset_name);
        else
            L = ls_segmentations@openorganelle_filebrowser(obj,obj.dataset_name);
        end
        
        end
        
        
        function L = ls_fibsem(obj)
        %L = ls_fibsem(obj)
        %Lists available FIBSEM volumes in AWS bucket for the current dataset. 
        %If output argument is requested, L will be a cell array of strings with
        %the FIBSEM volume names and nothing will be printed to the command
        %line.
        
        if nargout==0 
            ls_fibsem@openorganelle_filebrowser(obj,obj.dataset_name);
        else
            L = ls_fibsem@openorganelle_filebrowser(obj,obj.dataset_name);
        end
        
        end
        
        
        function L = ls_organelles(obj)
        %L = ls_organelles(obj)
        %Lists available organelle types for query based on the bounding box 
        %library located at 'bounding_boxes_path'.
        
        D0 = dir(fullfile(obj.bounding_boxes_path,obj.dataset_name,'*.json'));
        L = cell(numel(D0),1);
        
        if nargout==0, fprintf('\n----- CATALOGUED ORGANELLES for %s -----\n',obj.dataset_name); end
        
        for d = 1:numel(D0)
            [~,fname,~] = fileparts(D0(d).name);
            L{d} = fname;            
            if nargout==0, fprintf('%s\n',fname); end
        end
        
        if nargout==0, fprintf('\n'); clear L; end
        
        end
        
        
        function n = n_organelle(obj,organelle_type)
        %n = n_organelle(obj,organelle_type)
        %Returns the number of unique instances for organelle_type, based on the
        %number of bounding boxes catalogued.
        
        J = obj.load_bb_json(organelle_type);
        n = size(J.bounding_boxes,1);
        
        end
       
        
        function J = get_bb_json(obj,organelle_type)
        %J = get_bb_json(obj,organelle_type)
        %Loads the JSON catalog containing bounding boxes and additional info 
        %associated with the organelle organelle_type.
        %
        %Returns J, a struct constructed identically to the JSON.
        
        json_file = fullfile(obj.bounding_boxes_path,obj.dataset_name,sprintf('%s.json',organelle_type));
        str = fileread(json_file);
        J = jsondecode(str);
        
        end
                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                   PYTHON WRAPPERS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function S = get_organelle(obj,varargin)
        %S = get_organelle(obj,organelle_type,organelle_idx,Name,Value)
        %Retrieves volumes of instances of a given organelle type. Volumes 
        %are returned in (x,y,z)-dimensional order.
        %
        %
        %INPUTS
        %
        %organelle_type     label of the organelle of interest. Run
        %                   ls_organelles to retrieve the full list of
        %                   catalogued organelles.
        %
        %organelle_idx      organelle indices of interest as catalogued 
        %                   during instance segmentation. Expects 1-based 
        %                   indexing. Run ls_instances(organelle_type) to 
        %                   retrieve total number of instances. Can be a
        %                   single element, vector array, or 'all'. 
        %
        %Optional NAME,VALUE arguments
        %
        %'with_fibsem'      boolean flag to also retrieve the fibsem-uint16
        %                   intensity associated with each volume. Default
        %                   = false.
        %                   
        %                   ***N.B. Datasets were upscaled for organelle
        %                   annotation so that the highest resolution of
        %                   organelle label volumes is 2x higher than the
        %                   original fibsem resolution. To maintain a
        %                   consistent coordinate system, if with_fibsem =
        %                   true both the mask and fibsem data will be
        %                   returned at the same, lower resolution.
        %
        %'get_sparse'       boolean flag to retrieve the sparse label
        %                   volume, excluding all other objects that may 
        %                   overlap with the bounding box of interest.
        %                   Default = true.
        %
        %OUTPUT
        %
        %S struct with fields:
        %   
        %   data_mask       array containing mask volumes of each indexed
        %                   organelle. S.data_mask{j} corresponds to
        %                   organelle_idx(j). Note that if with_fibsem is false,
        %                   resolution will be 2x upscaled.
        %   
        %   data_fibsem     array containing fibsem intensity volumes of
        %                   each indexed organelle. S.data_fibsem{j}
        %                   corresponds to organelle_idx(j). If with_fibsem
        %                   = false, an empty array will be returned.
        %
        %   cosem_label     element or vector corresponding to the COSEM
        %                   label(s) of the requested organelle instance(s).
        %
        %   organelle_type  organelle_type requested.
        %
        %   organelle_idx   organelle_idx requested.
        %
        %   pixelResolution resolution at which the masks and/or fibsem
        %                   volumes are returned.
        %
        %   urlpath_mask    S3 URL used to retrieve masks.
        %
        %   urlpath_fibsem  S3 URL used to retrieve fibsem data. This will
        %                   be an empty string if with_fibsem = False.
        %
        %   bounding_boxes  array of bounding box coordinates that were 
        %                   used to retrieve image volumes. The format is 
        %                   bounding_boxes{j} = 
        %                       [xmin_j, ymin_j, zmin_j, xmax_j, ymax_j, zmax_j]
        %                   corresponding to organelle_idx(j). Note that
        %                   some of the catalogued bounding boxes are
        %                   padded so that the object does not necessarily
        %                   touch the edges of the volume.
        %
        %   catalog_path    path to the JSON catalog used to determine
        %                   bounding boxes of each organelle instance.
        %
        %   catalog_time    time at which the catalog was made. Given in 
        %                   seconds since epoch.
        
        %Parse inputs
        p = inputParser;
        p.addRequired('organelle_type',@ischar);
        p.addRequired('organelle_idx',@(x) strcmpi(x,'all') || (isnumeric(x) && all(x > 0)));
        p.addParameter('with_fibsem',false,@(x) numel(x)==1);
        p.addParameter('get_sparse',true,@(x) numel(x)==1);
        p.parse(varargin{:});
        
        organelle_type = p.Results.organelle_type;
        if strcmpi(p.Results.organelle_idx,'all')
            organelle_idx = int16(1:obj.n_organelle(organelle_type));
        else
            organelle_idx = int16(p.Results.organelle_idx);
        end
        with_fibsem = logical(p.Results.with_fibsem);
        get_sparse = logical(p.Results.get_sparse);
        
        %Run python function
        S = obj.oo2mat_mod.get_organelle(pyargs('dataset_name', obj.dataset_name, ...
            'organelle_type', organelle_type, ...
            'organelle_idx', organelle_idx, ...
            's3path', obj.s3bucket, ...
            'bbpath', obj.bounding_boxes_path, ...
            'with_fibsem', with_fibsem, ...
            'get_sparse', get_sparse));
        
        %Explicitly convert to MATLAB data types
        S = struct(S);
        S.data_mask = cellfun(@(x) logical(squeeze(x)), cell(S.data_mask), 'uniformoutput',false);
        S.data_fibsem = cellfun(@(x) double(squeeze(x)), cell(S.data_fibsem), 'uniformoutput',false);
        S.cosem_label = double(S.cosem_label);
        S.organelle_type = char(S.organelle_type);
        S.organelle_idx = double(S.organelle_idx);
        S.pixelResolution = struct(S.pixelResolution);
        S.pixelResolution.dimensions = cell2mat(cell(S.pixelResolution.dimensions));
        S.pixelResolution.unit = char(S.pixelResolution.unit);
        S.urlpath_mask = char(S.urlpath_mask);
        S.urlpath_fibsem = char(S.urlpath_fibsem);
        S.bounding_boxes = double(squeeze(S.bounding_boxes));
        S.catalog_path = char(S.catalog_path);
        S.catalog_time = double(S.catalog_time);
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [V,meta] = get_volume(obj,varargin)
        %[V,meta] = get_volume(obj,layer,bounding_box,Name,Value)
        %Retrieves the image volume associated with a given layer and
        %bounding box.Volumes are returned in (x,y,z)-dimensional order.
        %
        %INPUTS
        %
        %layer              string describing the dataset layer to
        %                   retrieve. Examples:
        %                       - 'fibsem-uint8'
        %                       - 'fibsem-uint16'
        %                       - 'mito_seg'
        %                       - 'lyso_seg'
        %
        %bounding_box       bounding box coordinates for the volume
        %                   requested, should be in the format
        %                   [xmin, ymin, zmin, xmax, ymax, zmax],
        %                   assuming 1-based indexing.
        %
        %Optional NAME,VALUE arguments
        %
        %'resolution'       hierarchical resolution level to retrieve. 
        %                   An integer from 0 to 4 (or 5 depending on the
        %                   dataset) where 0 is the highest resolution.
        %                   Default = 0.
        %                   
        %                   ***N.B. label (segmentation) volumes were
        %                   upscaled from the fibsem acquisition
        %                   resolution. To compare label and fibsem data 
        %                   under equivalent spatial dimensions, both in 
        %                   physical space and matrix size, set the the 
        %                   label retrieval resolution to j+1 where j is
        %                   the fibsem retrieval resolution. e.g.,
        %               
        %                   V_fibsem = obj.load_volume('fibsem-uint8',[1,1,1,100,100,100],'resolution',0);
        %                   V_mito = obj.load_volume('mito_seg',[1,1,1,100,100,100],'resolution',1);
        %
        %OUTPUTS
        %
        %V                  M-by-N-by-P matrix corresponding to the image
        %                   type and bounding box requested.
        %                   
        %meta               metadata related to the retrieved volume.
        %
        
        p = inputParser;
        p.addRequired('layer',@ischar);
        p.addRequired('bounding_box',@(x) numel(x)==6)
        p.addParameter('resolution',0,@(x) round(x)==x)
        p.parse(varargin{:})
        
        layer = p.Results.layer;
        bounding_box_py = int16(reshape(p.Results.bounding_box,1,6)-[1,1,1,0,0,0]);
        resolution = int8(p.Results.resolution);
        
        S = obj.oo2mat_mod.get_volume(pyargs('dataset_name',obj.dataset_name,...
            'layer',layer,...
            'bounding_box',bounding_box_py,...
            'resolution',resolution,...
            's3path',obj.s3bucket));
        S = struct(S);
        
        %Convert to MATLAB datatype
        switch layer
            case 'fibsem-uint8'
                V = uint8(S.V);
            case 'fibsem-uint16'
                V = uint16(S.V);
            otherwise %Label matrix
                V = uint64(S.V);
        end
        
        %Store metadata
        meta = struct('layer',layer,...
            'bounding_box_matlab',p.Results.bounding_box,...
            'bounding_box_py',bounding_box_py,...
            'resolution',resolution,...
            'urlpath',char(S.urlpath),...
            'retrieval_time',char(S.retrieval_time));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [catalog_filename,J] = do_catalog(obj,varargin)
        %[catalog_filename,J] = DO_CATALOG(obj,layer,Name,Value)
        %Detects instances of organelles in the layer label and saves their
        %bounding boxes to a JSON catalog file. Bounding boxes correspond 
        %to the highest resolution image volume available.  
        %
        %INPUT
        %
        %layer          string corresponding to the name of the 
        %               segmentation layer to catalog. e.g.
        %                   - 'mito_seg'
        %                   - 'cent_seg'
        %
        %               Run .ls_segmentations to see the full list of 
        %               available segmentation layers.
        %
        %Optional NAME,VALUE arguments
        %
        %'overwrite_if_exists'
        %               boolean flag to overwrite the existing catalog if 
        %               it exists. If set to false, no cataloguing is 
        %               performed and the existing JSON is returned.
        %               Default = false.
        %
        %'resolution_to_use'
        %               integer corresponding to the resolution at which
        %               the volume should be downsampled for cataloguing.
        %               Set to py.None (default) or an empty object to 
        %               select the resolution automatically.
        %
        %'skip_check'   boolean flag to skip process that checks 
        %               whether separate object instances were split (or 
        %               merged) into a single bounding box (or different
        %               bounding boxes) due to the initial detection at a
        %               lower resolution. Default = false. 
        %
        %               This check is generally the slowest step in the
        %               cataloguing procedure. Without it however, some
        %               instances may be omitted or from the catalog and
        %               subsequently cannot be retrieved using the
        %               .get_organelle class method. 
        %
        %               See end of function documentation for more details.
        %
        %
        %OUTPUTS
        %
        %catalog_filename   full path to JSON catalog.
        %
        %J                  struct whose contents correspond to the JSON 
        %                   catalog. 
        %
        %
        %***WARNING REGARDING THE MERGE CHECK: At present, the cataloguing 
        %function cannot be keyboard-interrupted from MATLAB. Therefore if 
        %a very large number of bounding boxes need to be checked for 
        %merged instances, your MATLAB session may be blocked for a while.
        %
        %For example, there are 2.9 million bounding boxes for the
        %downscaled ribo_seg: in its current state, this would take a ~4
        %days (on workstation-02) to check each bounding box. For cases 
        %like this, consider submitting the job to the cluster or at least
        %running a separate worker for the task. 
        %
        %A parallelized solution is still in progress.
        
        p = inputParser;
        p.addRequired('layer',@(x) any(strcmp(x,obj.ls_segmentations)));
        p.addParameter('overwrite_if_exists',false,@islogical);
        p.addParameter('resolution_to_use',py.None,@(x) isempty(x) || (x==py.None) || (round(x)==x));
        p.addParameter('skip_merge_check',false,@islogical);
        p.parse(varargin{:});
        
        layer = p.Results.layer;
        overwrite_if_exists = logical(p.Results.overwrite_if_exists);
        resolution_to_use = p.Results.resolution_to_use;
        if isnumeric(resolution_to_use)
            resolution_to_use = int8(resolution_to_use);
        else
            resolution_to_use = py.None;
        end
        skip_merge_check = logical(p.Results.skip_merge_check);
        
        %Run oo2mat cataloguing function
        catalog_filename = obj.oo2mat_mod.do_bb_catalog(pyargs(...
            'dataset_name',obj.dataset_name,...
            'layer',layer,...
            'dest_dir',obj.bounding_boxes_path,...
            's3path',obj.s3bucket,...
            'use_res',resolution_to_use,...
            'skip_check',skip_merge_check,...
            'allow_overwrite',overwrite_if_exists));
            
        catalog_filename = char(catalog_filename);
        if nargout==2
            J = obj.get_bb_json(layer);
        end
        
        end
        
    end
    
end