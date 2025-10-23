import os
import pandas as pd
from scipy.sparse import csr_matrix
import scanpy as sc

def castCounts(transcs):
    transcs = transcs.groupby(['cell_id','feature_name'])['transcript_id']
    transcs = transcs.aggregate(func='count').reset_index()
    
    row = transcs.cell_id.cat.codes
    col = transcs.feature_name.cat.codes
    
    sparse_matrix = csr_matrix((transcs["transcript_id"], (row, col)),
                               shape=(transcs.cell_id.cat.categories.size, transcs.feature_name.cat.categories.size))
    return sparse_matrix

def loadXeniumCounts(cells,xenium_dir,qv_thr=20):
    cells.cell_id = pd.Categorical(cells.cell_id)
    cells.cell_id = cells.cell_id.cat.remove_unused_categories()
    transcs = pd.read_parquet(xenium_dir+"/transcripts.parquet",filters=[("cell_id","in",list(cells.cell_id)),('is_gene','=',True),('qv','>=',qv_thr)])
    transcs.cell_id = pd.Categorical(transcs.cell_id,categories=cells.cell_id)
    transcs.feature_name = pd.Categorical(transcs.feature_name)
    cyto = castCounts(transcs[transcs.overlaps_nucleus==0])
    nucl = castCounts(transcs[transcs.overlaps_nucleus==1])
    
    adata=sc.read_10x_h5(xenium_dir+"/cell_feature_matrix.h5")
    adata = adata[transcs.cell_id.cat.categories, transcs.feature_name.cat.categories].copy()
    adata.obs['cell_id'] = adata.obs_names
    adata.layers["spliced"] = cyto
    adata.layers["unspliced"] = nucl
    adata.obsm['spatial'] = cells[['x_centroid','y_centroid']].to_numpy()
    sc.pp.calculate_qc_metrics(adata,inplace=True)
    return adata

def find_folder(folder_name, search_path):
    for root, dirs, files in os.walk(search_path):
        if folder_name in dirs:
            return os.path.join(root, folder_name)
    return None