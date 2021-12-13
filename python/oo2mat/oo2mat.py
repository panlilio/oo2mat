import numpy as np
import os
import json
from fibsem_tools.io import read_xarray
from datetime import datetime
import s3fs
from skimage.measure import label, regionprops
import multiprocessing as mp

"""
Package to retrieve image volumes from the OpenOrganelle project with output formatted to be the most MATLAB-friendly.

Author: Mia Panlilio (2021-11-15)
"""

def get_organelle(dataset_name='jrc_hela-2',
                   organelle_type='mito_seg',
                   organelle_idx=1,
                   s3path='s3://janelia-cosem/',
                   bbpath='/research/sharedresources/cbi/common/mia/OpenOrganelle/lib/organelle_bounding_boxes/',
                   with_fibsem=False,
                   get_sparse=True):

    #Load bounding box catalog file
    catalog_path = os.path.join(bbpath, dataset_name, f'{organelle_type}.json')
    with open(catalog_path) as catalog_json:
        J = json.load(catalog_json)

    #Load dataset into xarray
    if with_fibsem:
        #Note that the label matrices were upscaled 2x, so there is a difference in equivalent resolution levels as
        #organized in the N5 between the fibsem and label datasets
        urlpath_mask = os.path.join(s3path,dataset_name,f'{dataset_name}.n5','labels',organelle_type,'s1')
        urlpath_fibsem = os.path.join(s3path,dataset_name,f'{dataset_name}.n5','em','fibsem-uint16','s0')
        X_fibsem = read_xarray(urlpath_fibsem)
    else:
        urlpath_mask = os.path.join(s3path,dataset_name,f'{dataset_name}.n5','labels',organelle_type,'s0')
        urlpath_fibsem = ''
        X_fibsem = []

    if isinstance(organelle_idx,int):
        organelle_idx = [organelle_idx]

    if not(isinstance(organelle_idx,list)):
        organelle_idx = organelle_idx.tolist()

    #Translate from 1-indexed to 0-indexed
    organelle_idx_1 = organelle_idx.copy()
    organelle_idx = [idx-1 for idx in organelle_idx]

    #Load xarray associated with
    X_mask = read_xarray(urlpath_mask)
    BB = get_bounding_boxes(J,organelle_idx,resolution=X_mask.pixelResolution["dimensions"])
    L = [J["cosem_labels"][idx] for idx in organelle_idx]

    print(f'{datetime.now()} [INFO] User requested {organelle_type} at index = {organelle_idx_1}.')
    M = list()
    F = list()
    for idx,bb in enumerate(BB):
        print(f'{datetime.now()} [INFO] Retrieving index {idx+1}/{len(BB)}...')
        m = np.array(X_mask[bb[0]:bb[3],bb[1]:bb[4],bb[2]:bb[5]],dtype=int)
        if get_sparse:
            m = m*np.equal(m,L[idx])
        M.append(m)

        if len(X_fibsem) > 0:
            f = np.array(X_fibsem[bb[0]:bb[3],bb[1]:bb[4],bb[2]:bb[5]],dtype=np.uint16)
            F.append(f)

        print(f'{datetime.now()} [INFO] Retrieval of index {idx+1}/{len(BB)} done.')

    data = {"data_mask" : M,
            "data_fibsem" : F,
            "cosem_label" : np.array(L,dtype=int),
            "organelle_type" : organelle_type,
            "organelle_idx" : np.array(organelle_idx,dtype=int),
            "pixelResolution" : X_mask.pixelResolution,
            "urlpath_mask" : urlpath_mask,
            "urlpath_fibsem" : urlpath_fibsem,
            "bounding_boxes" : np.array(BB),
            "catalog_path" : catalog_path,
            "catalog_time" : J["catalog_time"],
            }

    return data

def get_bounding_boxes(catalog, idx, resolution):
    J = catalog
    if isinstance(idx,int):
        idx = [idx]
    BB = [J['bounding_boxes'][ix] for ix in idx]

    catalog_res = J["pixelResolution"]["dimensions"]
    requested_res = resolution
    scaling_factor = np.array(catalog_res)/np.array(requested_res)
    BB = np.array([bb * np.tile(scaling_factor, 2) for bb in BB], dtype=np.int16)
    BB = BB.tolist()

    return BB

def get_volume(dataset_name='jrc_hela-2',
               layer='fibsem-uint16',
               bounding_box=(0,0,0,100,100,100),
               resolution=0,
               s3path='s3://janelia-cosem/'):

    if layer.startswith('fibsem'):
        subdir = 'em'
    else:
        subdir = 'labels'

    urlpath = os.path.join(s3path,dataset_name,'{}.n5'.format(dataset_name),subdir,layer,'s{}'.format(resolution))
    X = read_xarray(urlpath)
    V = np.array(X[bounding_box[0]:bounding_box[3],bounding_box[1]:bounding_box[4],bounding_box[2]:bounding_box[5]])

    data = {"V" : V, "urlpath" : urlpath, "retrieval_time" : f"{datetime.now()}"}

    return data

def get_catalog(dataset_name='jrc_hela-2',
                organelle_type='mito_seg',
                bbpath='/research/sharedresources/cbi/common/mia/OpenOrganelle/lib/organelle_bounding_boxes/'):
    catalog_path = os.path.join(bbpath, dataset_name, f'{organelle_type}.json')
    with open(catalog_path) as catalog_json:
        J = json.load(catalog_json)
    return J

def ls_datasets(s3path='s3://janelia-cosem/'):
    fs = s3fs.S3FileSystem(anon=True)
    L = fs.ls(s3path)
    L = [os.path.basename(ll) for ll in L if not(ll.endswith(('.md','.json')))]
    return L

def ls_segmentations(dataset_name='jrc_hela-2',s3path='s3://janelia-cosem/'):
    fs = s3fs.S3FileSystem(anon=True)
    L = fs.ls(os.path.join(s3path,dataset_name,'{}.n5'.format(dataset_name),'labels'))
    L = [os.path.basename(ll) for ll in L if not (ll.endswith(('.md','.json')))]
    return L

def ls_fibsem(dataset_name='jrc_hela-2',s3path='s3://janelia-cosem/'):
    fs = s3fs.S3FileSystem(anon=True)
    L = fs.ls(os.path.join(s3path, dataset_name, '{}.n5'.format(dataset_name), 'em'))
    L = [os.path.basename(ll) for ll in L if not (ll.endswith(('.md', '.json')))]
    return L

def do_bb_catalog(dataset_name='jrc_hela-2',
                  layer='mito_seg',
                  dest_dir='/research/sharedresources/cbi/common/mia/OpenOrganelle/lib/organelle_bounding_boxes/',
                  s3path='s3://janelia-cosem/',
                  use_res=None,
                  skip_check=False,
                  allow_overwrite=False):

    t0 = datetime.now()
    catalog_name = os.path.join(dest_dir,dataset_name,"{}.json".format(layer))
    if os.path.exists(catalog_name) and not(allow_overwrite):
        print(f"{datetime.now()} [ERROR] Catalog for {dataset_name}{os.path.sep}{layer} already exists: {catalog_name}. Nothing to do.")
        return catalog_name
    elif os.path.exists(catalog_name):
        print(f"{datetime.now()} [WARNING] Catalog for {dataset_name}{os.path.sep}{layer} already exists: {catalog_name}. Existing file will be overwritten.")

    layerpath = os.path.join(s3path,dataset_name,f"{dataset_name}.n5",'labels',layer)

    #Read in dataset at increasing resolution until objects are detected
    X0 = read_xarray(os.path.join(layerpath,'s0'))

    if use_res is not None:
        urlpath = os.path.join(layerpath, f"s{use_res}")
        X = read_xarray(urlpath)
        BB,LL = find_bb_at_res(X,X0,skip_check=skip_check)

        if BB is None:
            print(f"{datetime.now()} [ERROR] No objects found in layer {layer} at resolution s{use_res}. No catalog created.")
            return
    else:
        res = 3
        while res >= 0:
            urlpath = os.path.join(layerpath,f"s{res}")
            X = read_xarray(urlpath)
            BB,LL = find_bb_at_res(X,X0,skip_check=skip_check)

            if BB is not None:
                res = -1
            else:
                print(f"{datetime.now()} [INFO] No objects found in layer {layer} at resolution s{res}. Trying s{res-1}.")
                res -= 1

    t1 = datetime.now()
    print(f"{datetime.now()} [INFO] Total time to catalog = {t1-t0}.")
    data = {"bounding_boxes": BB.tolist(),
            "cosem_labels": LL.tolist(),
            "name": X0.name,
            "pixelResolution": X0.pixelResolution,
            "merge_check_skipped" : skip_check,
            "urlpath_bb": X0.urlpath,
            "urlpath_downsample": X.urlpath,
            "catalog_time": datetime.timestamp(datetime.now()),
            "time_to_catalog" : f"{t1-t0}"
            }

    with open(f"{catalog_name}", "w") as write_file:
        json.dump(data, write_file, indent=2)

    return catalog_name

def find_bb_at_res(X,X0,skip_check=False):
    [B, n_obj] = label(np.array(X), return_num=True)
    if n_obj==0:
        return None,None

    scaling_factor = np.array(X.pixelResolution['dimensions']) / np.array(X0.pixelResolution['dimensions'])
    R = regionprops(B, intensity_image=np.array(X))
    L = np.array([r.max_intensity for r in R], dtype=np.uint64)  # COSEM label associated with each region

    #Transform bounding box to highest resolution coordinates, with a  +/- 1 pixel (downscaled) cushion in case it was
    #omitted during the N5 downscaling
    R = np.array([r.bbox + np.array([-1, -1, -1, 1, 1, 1]) for r in R], dtype=np.int16)  # pixel cushion
    R = np.array([r * np.tile(scaling_factor, 2) for r in R], dtype=np.int16)  # upscale bounding box
    R = np.array([np.maximum(0, np.minimum(np.tile(X0.shape, 2), r)) for r in R],
    dtype=np.uint16)  # Ensure that cushioned window is still within data bounds

    #Iteratively adjust bounding boxes if lower resolution labels merged separate objects
    BB = R.copy()
    LL = L.copy()

    if not(skip_check):
        R,L = split_check(R,L)
        BB, LL = merge_check(X0,R,L)
    else:
        print(f'{datetime.now()} [WARNING] No check performed to ensure that no objects were merged at a lower resolution. Some organelles (instance labels) may be omitted from the catalog.')

    print(f'{datetime.now()} [INFO] Total number of organelles found = {BB.shape[0]}.')
    return BB, LL

def split_check(B,L):

    print(f'{datetime.now()} [INFO] Checking bounding box regions for objects split at lower resolution.')
    n_obj = B.shape[0]

    BB = {ll : np.zeros((0,6)) for ll in L}
    for idx,ll in enumerate(L):
        BB[ll] = np.vstack((BB[ll], B[idx]))
    BB_min = [np.min(BB[ll], axis=0) for ll in BB.keys()]
    BB_max = [np.max(BB[ll], axis=0) for ll in BB.keys()]

    B = np.array([np.append(mn[:3], mx[3:]) for mn, mx in zip(BB_min, BB_max)],dtype=np.uint16)
    L = np.array([key for key in BB.keys()],dtype=np.uint64)

    if n_obj != B.shape[0]:
        print(f'{datetime.now()} [INFO] Merged split objects. Total number of unique objects found at downscaled resolution: {n_obj} -> {B.shape[0]}')

    print(f'{datetime.now()} [INFO] Split check complete.')

    return B,L


def merge_check(X0=None,B=None,L=None):

    print(f'{datetime.now()} [INFO] Checking bounding box regions for objects merged at lower resolution.')

    # Create bounding box dictionary using COSEM label as key
    BB = {ll : bb for (ll,bb) in zip(L,B)}
    is_cushioned = {ll : [True] for ll in L} #Flag to keep track of whether the bounding box corresponds to the cushioned version (and is therefore not necessarily the minimal bounding box)
    LL = L.copy()

    for idx,ll in enumerate(L):
        print(f'{datetime.now()} [INFO] Checking region {idx+1}/{len(L)}.')
        bb = B[idx]
        x = np.array(X0[bb[0]:bb[3], bb[1]:bb[4], bb[2]:bb[5]])
        b = label(x)
        Rx = regionprops(b, intensity_image=x)
        Lx = np.array([r.max_intensity for r in Rx], dtype=np.uint64)
        Rx = np.array([r.bbox for r in Rx], dtype=np.int16)
        for lx, rx in zip(Lx, Rx):
            #Translate to global volume coords
            rx_bbox = np.array(rx) + np.tile(bb[0:3], 2)
            if lx in LL:
                # Labelled object exists in another bounding box, add to list of potential bounding boxes
                BB[ll] = np.vstack((BB[ll],rx_bbox))
                is_cushioned[ll].append(False)
            else:
                print(f'{datetime.now()} [INFO] Additional object found.')
                LL = np.append(LL,lx)
                BB[lx] = [rx_bbox]
                is_cushioned[lx] = [False]

    #Remove cushioned bounding boxes
    for ll in BB.keys():
        BB[ll] = [bb for (bb,ic) in zip(BB[ll],is_cushioned[ll]) if not ic]

    #Get the minimum bounding box for each label
    BB_min = np.array([np.min(BB[ll],axis=0) for ll in BB.keys()])
    BB_max = np.array([np.max(BB[ll],axis=0) for ll in BB.keys()])

    B = np.array([np.append(mn[:3],mx[3:]) for mn,mx in zip(BB_min,BB_max)],dtype=np.uint16)
    L = np.array([key for key in BB.keys()],dtype=np.uint64)

    return B,L

def merge_check_multi(X0=None,B=None,L=None):
    BB = {ll: bb for (ll, bb) in zip(L, B)}
    is_cushioned = {ll: [True] for ll in L}  # Flag to keep track of whether the bounding box corresponds to the cushioned version (and is therefore not necessarily the minimal bounding box)
    LL = L.copy()

    def merge_check_to_apply(bb,X0):
        x = np.array(X0[bb[0]:bb[3], bb[1]:bb[4], bb[2]:bb[5]])
        b = label(x)
        Rx = regionprops(b, intensity_image=x)
        Lx = np.array([r.max_intensity for r in Rx], dtype=np.uint64)
        Rx = np.array([r.bbox for r in Rx], dtype=np.int16)
        return Rx,Lx





def test(var):
    print(type(var))
    if var is None:
        print('hi')
    print(var)
    return var