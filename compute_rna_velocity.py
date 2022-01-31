import anndata
import igraph
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2.robjects as robj
import scanpy as sc
import scipy
import scipy.optimize
import scvelo as scv
import sklearn
import velocyto as vcy
import os

from collections import Counter
from IPython.core.display import display, HTML
from numpy_groupies import aggregate, aggregate_np
from rpy2.robjects.packages import importr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import pdist, squareform


fn1 = "/stanley/levin_dr/kwanho/projects/Amanda/Atlas/velocity_1.5m/analysis/combined_adata_for_velocity.h5ad"
if not os.path.isfile(fn1): 
    samples = ['Pm1_5_1','Pm1_5_2','Pm1_5_3','Mm1_5_1','Mm1_5_2','Mm1_5_3']
    adata_list = [anndata.read("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/velocity_1.5m/kb_counts/" + x +"/counts_unfiltered/adata.h5ad") for x in samples]
    # label samples - trace each row back to its original anndata object.
    for i in range(len(adata_list)):
        temp_adata = adata_list[i].copy()
        temp_adata.obs['run'] = samples[i]
        temp_adata.obs['bcs'] = temp_adata.obs.index
        temp_adata.obs.index = temp_adata.obs['bcs'] + '.' + temp_adata.obs['run']
        #temp_adata.var['gene_id'] = temp_adata.var.index
        #bad_name = temp_adata.var.gene_name.str.find('/')
        #bad_name[bad_name>0]
        temp_adata.var.gene_name = temp_adata.var.gene_name.str.replace("/", "_")
        temp_adata.var.reset_index(drop=False, inplace=True)
        adata_list[i] = temp_adata
    
    # concatenate adata
    #adata_with_var = anndata.concat(adata_list, merge='same')
    adata_full = anndata.concat(adata_list, merge='same')
    adata_full.var.set_index('gene_name', inplace=True, drop=False)
    adata_full.obs
    adata_full.var
    adata_full.var.index = adata_full.var.index.str.split('.').str[0]
    adata_full.var_names_make_unique()
    adata_full.var
    # take list of cell names from fully processed Seurat object
    df = pd.read_csv('/stanley/levin_dr/kwanho/projects/Amanda/Atlas/velocity_1.5m/seurat_data/CellID_and_FinalName.tsv', sep='\t', header=None)
    cell_id = df[0].to_numpy()
    cell_type = df[1].to_numpy()
    # subset adata
    adata = adata_full[cell_id]
    # add metadata
    adata.layers["ambiguous"] = scipy.sparse.csr_matrix(np.zeros(adata.X.shape))
    adata.obs["FinalName"] = cell_type
    # additional metadata
    metadf = pd.read_csv("/stanley/levin_dr/kwanho/projects/Amanda/Atlas/velocity_1.5m/seurat_data/table_additional_metadata.tsv", sep='\t', header=0, index_col=0)
    # combine with existing metadata
    adata.obs = adata.obs.merge(metadf, left_index=True, right_index=True)
    adata.obs = adata.obs.astype({'MetName': 'category', 'Phase': 'category'})
    # apply umap embeddings computed in the Seurat object
    emb = pd.read_csv('/stanley/levin_dr/kwanho/projects/Amanda/Atlas/velocity_1.5m/seurat_data/umap_embedding_from_Seurat.tsv', sep='\t', header=0, index_col=0)
    adata.obsm['X_umap'] = emb.to_numpy()
    adata.var.index.name=None
    adata.write(fn1)
    scv.pl.scatter(adata, color='MetName', basis='umap', fontsize=10, legend_loc='on data', legend_fontsize=7, save="figures/umap_MetName.png")
else:
    print("prepped anndata found!")

#####################
## scVelo pipeline
####################3

scv.set_figure_params(dpi=300, dpi_save=300)

if not os.path.isdir('figures'):
    os.mkdir("figures")

fn2 = "/stanley/levin_dr/kwanho/projects/Amanda/Atlas/velocity_1.5m/analysis/reg_out_cc/adata_reg_out_cc_dynamical_model.h5ad"
if not os.path.isfile(fn2):
    adata = scv.read(fn1, cache=True)
    
    #scv.pp.filter_genes(adata, min_shared_counts=20)
    #scv.pp.normalize_per_cell(adata)
    #scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)  # keep 2000 highly variable genes
    #scv.pp.log1p(adata)
    # above block is summarized into this function
    print("normalize!")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)  
    
    #print("remove cycling genes!")
    #s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
    #s_mk = np.in1d( adata.var_names.values.astype(str), s_genes)
    #adata = adata[:,~s_mk]
    #g2m_mk = np.in1d( adata.var_names.values.astype(str), g2m_genes)
    #adata = adata[:,~g2m_mk]
    
    print("reg out cc")
    sc.pp.regress_out(adata, keys=["S.Score", "G2M.Score"])
    
    # compute first and second order moments (means and uncentered variances) among NNs in PCA space
    print("compute moments!")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    
    # dynamical model to learn the full transcriptional dynamics of splicing kinetics
    # It is solved in a likelihood-based expectation-maximization framework, by iteratively estimating the parameters of reaction rates and latent cell-specific variables, i.e. transcriptional state and cell-internal latent time. It thereby aims to learn the unspliced/spliced phase trajectory for each gene.
    print("recover dynamics!")
    scv.tl.recover_dynamics(adata, n_jobs=8)
    
    # save
    print("saving!")
    adata.write('adata_reg_out_cc_dynamical_model.h5ad', compression='gzip')
else:
    print("dynamical model found!")




adata = scv.read(fn2, cache=True)

counts_s = scv.utils.sum_var(adata.layers['spliced'])
counts_u = scv.utils.sum_var(adata.layers['unspliced'])
fractions_u = counts_u / (counts_s + counts_u)
scv.pl.scatter(adata, color=fractions_u, smooth=True, save="figures/frac_unspliced.png")

# The dynamical model recovers the latent time of the underlying cellular processes, based only on its transcriptional dynamics
print("latent time!")
adata.obs['root_cells'] = adata.obs.MetName=="aRG"
adata.obs['end_cells'] = adata.obs.MetName=="CFuPN"
scv.tl.latent_time(adata, root_key='root_cells', end_key='end_cells')
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="figures/umap_latent_time.png")

# use latent time as a regularization for velocity estimation
scv.tl.velocity(adata, mode='dynamical', use_latent_time=True)
scv.tl.velocity_graph(adata)

# Kindetic rate parameters
# The rates of RNA transcription, splicing and degradation are estimated without the need of any experimental data.
# They can be useful to better understand the cell identity and phenotypic heterogeneity.
print("rates!")
df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
df.to_csv("df_rates.csv")

#The estimated gene-specific parameters comprise rates of transription (fit_alpha), splicing (fit_beta), degradation (fit_gamma), switching time point (fit_t_), a scaling parameter to adjust for under-represented unspliced reads (fit_scaling), standard deviation of unspliced and spliced reads (fit_std_u, fit_std_s), the gene likelihood (fit_likelihood), inferred steady-state levels (fit_steady_u, fit_steady_s) with their corresponding p-values (fit_pval_steady_u, fit_pval_steady_s), the overall model variance (fit_variance), and a scaling factor to align the gene-wise latent times to a universal, gene-shared latent time (fit_alignment_scaling).

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

plt.savefig("figures/rates.png")
scv.get_df(adata, 'fit*', dropna=True).head()

# save
print("saving!")
adata.write('adata_with_velocity.h5ad', compression='gzip')

# plot
scv.pl.velocity_embedding_stream(adata, basis='umap', color='MetName', save="figures/umap_velocity_embedding_stream.png")
scv.pl.velocity_embedding_grid(adata, basis='umap', color='MetName', save='figures/umap_velocity_embedding_grid.png', scale=0.2)
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="figures/velocity_confidence.png")
scv.pl.velocity_graph(adata, threshold=.1, color='MetName', save="figures/velocity_graph.png")

print("top likelihood genes!")
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='MetName', n_convolve=100, save="figures/heatmap_gene_exp_in_latent_time.png")

# Driver genes display pronounced dynamic behavior and are systematically detected via their characterization by high likelihoods in the dynamic model.
print("driver genes!")
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], color='MetName', ncols=5, frameon=False, save="figures/driver_genes.png")
scv.pl.scatter(adata, color='MetName', x='latent_time', y=top_genes[:15], ncols=5, frameon=False, save="figures/driver_genes_latent_time.png")

#scv.tl.velocity_pseudotime(adata, root_key='root_cells')
#scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')
#scv.tl.paga(adata, groups='MetName', root_key='root_cells')
#scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save='figures/paga.png')

scv.pl.velocity(adata, top_genes[:6], ncols=2, add_outline=True, save="figures/top_genes_velocity_plot.png")

# Cluster-specific top-likelihood genes
# Moreover, partial gene likelihoods can be computed for a each cluster of cells to enable cluster-specific identification of potential drivers.
scv.tl.rank_dynamical_genes(adata, groupby='MetName')
adata.write('adata_final_no_ccg.h5ad', compression='gzip')

df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(5)
df.to_pickle("df_rank_dynamical_genes.pkl")
df.to_csv("df_rank_dynamical_genes.csv")
# pd.read_pickle()

scv.pl.velocity(adata, df.head(20)["aRG II"], ncols=2, add_outline=True, save="figures/top_aRG2_genes_velocity_plot.png")

nc = 5
clusters = adata.obs['MetName'].cat.categories
genes = df.head(nc)[clusters].to_numpy().flatten(order='F')

fig, big_axes = plt.subplots(nrows=len(clusters), ncols=1, sharey=True) 
for cluster, big_ax in zip(clusters, big_axes):
    big_ax.set_title(cluster, fontsize=16, pad=20)
    big_ax.tick_params(axis='both', which='both', top='off', bottom='off', left='off', right='off')
    big_ax.set_xticks([])
    big_ax.set_yticks([])
    big_ax._frameon = False

i=1
for cluster, gene in zip(np.repeat(clusters, nc), genes):
    ax = fig.add_subplot(len(clusters), nc, i)
    i += 1
    scv.pl.scatter(adata, gene, color='MetName', legend_loc='none', frameon=False, ax=ax)
    ax.get_legend().remove()


# x-axis = spliced, y-axis = unspliced
fig.set_size_inches(20, 40)
plt.savefig("figures/cluster_specific_driver_genes.png", bbox_inches='tight')




