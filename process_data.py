# Import dependencies
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import pandas as pd
import decoupler as dc
import gc as py_gc

# Data reading
input_file = "SjS/data/training/06e91d29-15e6-4590-9462-7ac166c445e9.h5ad"
adata = sc.read_h5ad(input_file)

gene_map = dict(zip(adata.var_names, adata.var['feature_name'].astype(str)))
X = adata.raw.X
obs = pd.DataFrame()
obs['Status'] = adata.obs['disease'].tolist()
obs['ind_cov'] = adata.obs['donor_id'].tolist()
obs['cg_cov'] = adata.obs['cell_type'].tolist()

var_names = adata.raw.var_names.tolist()
var = pd.DataFrame(index=var_names)
tdata = ad.AnnData(X, obs=obs, var=var, dtype='int32')
tdata.obs_names = adata.obs_names

tdata.var_names = [gene_map.get(g, g) for g in tdata.var_names]
tdata.var_names_make_unique()

del(adata)
py_gc.collect()

# Quality control
tdata.var['mt'] = tdata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(tdata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
print(tdata.var['mt'].sum())

# Cells and genes filtering
ncells = tdata.shape[0]
sc.pp.filter_cells(tdata, min_genes=100)
sc.pp.filter_cells(tdata, max_genes=4000)
sc.pp.filter_genes(tdata, min_cells=ncells/1000)

# Save pseudobulk data
pseudobulk = dc.get_pseudobulk(tdata, sample_col = "ind_cov", groups_col = "cg_cov", use_raw = False, mode='sum', min_cells=1, min_counts = 1)
pseudobulk_df = pseudobulk.to_df()
pseudobulk_df.to_csv("SjS/data/training/pseudobulk.tsv", sep="\t")

# Normalization
sc.pp.normalize_total(tdata)
sc.pp.log1p(tdata)
py_gc.collect()

# Select highly variable genes from training dataset
sc.pp.highly_variable_genes(tdata, n_top_genes=2000, inplace=True, subset=True)

# Regress by total and mitochondrial expression
sc.pp.regress_out(tdata, ['total_counts', 'pct_counts_mt'])

# Scale each dataset independently
sc.pp.scale(tdata)

# Write training dataset
tdata.write_h5ad("SjS/data/training/training_processed.h5ad")
