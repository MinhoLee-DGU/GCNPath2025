import numpy as np
import pandas as pd
import os
import csv
import scipy
import torch
import torch.nn as nn
from torch_geometric.data import Data, Batch
from torch_geometric.nn import graclus, max_pool
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer


def get_genes_graph(genes_path, save_path, method='pearson', thresh=0.95, p_value=False):
    """
    determining adjaceny matrix based on correlation
    :param genes_exp_path:
    :return:
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    genes_exp_df = pd.read_csv(os.path.join(genes_path, 'EXP.csv'), index_col=0)

    # calculate correlation matrix
    genes_exp_corr = genes_exp_df.corr(method=method)
    genes_exp_corr = genes_exp_corr.apply(lambda x: abs(x))
    n = genes_exp_df.shape[0]

    # binarize
    if p_value == True:
        dist = scipy.stats.beta(n / 2 - 1, n / 2 - 1, loc=-1, scale=2)
        thresh = dist.isf(0.05)

    adj = np.where(genes_exp_corr > thresh, 1, 0)
    adj = adj - np.eye(genes_exp_corr.shape[0], dtype=np.int)
    edge_index = np.nonzero(adj)
    np.save(os.path.join(save_path, 'edge_index_{}_{}.npy').format(method, thresh), edge_index)

    return n, edge_index


def ensp_to_hugo_map():
    with open('./data/9606.protein.info.v11.0.txt') as csv_file:
        next(csv_file)  # Skip first line
        csv_reader = csv.reader(csv_file, delimiter='\t')
        ensp_map = {row[0]: row[1] for row in csv_reader if row[0] != ""}

    return ensp_map


def hugo_to_ncbi_map():
    with open('./data/enterez_NCBI_to_hugo_gene_symbol_march_2019.txt') as csv_file:
        next(csv_file)  # Skip first line
        csv_reader = csv.reader(csv_file, delimiter='\t')
        hugo_map = {row[0]: int(row[1]) for row in csv_reader if row[1] != ""}

    return hugo_map


def add_rand_noise(cell_data, std=1.0, noise_seed=2021) :
    if noise_seed==-1 :
        pass
    else :
        rng = np.random.default_rng(seed=noise_seed)
        noise = rng.normal(0, std, size=cell_data.shape)
        cell_data = cell_data + noise
        print(f"# Adding noise N(0, {std}) to Cell Data...")
        print(f"{noise}\n")
    return cell_data


# def save_cell_graph_exp(genes_path, save_path) :
#     
#     if not os.path.exists(save_path):
#         os.makedirs(save_path)
#     
#     exp = pd.read_csv(os.path.join(genes_path, 'exp.csv'), index_col=0)
# 
#     index = exp.index
#     columns = exp.columns
# 
#     scaler_exp = StandardScaler()
#     exp = scaler_exp.fit_transform(exp)
#     
#     imp_mean = SimpleImputer()
#     exp = imp_mean.fit_transform(exp)
#     exp = pd.DataFrame(exp, index=index, columns=columns)
#     
#     cell_dict = {}
#     cell_names = exp.index
#     
#     for i in cell_names:
#         cell_dict[i] = Data(x=torch.tensor([exp.loc[i]], dtype=torch.float).T)
#     
#     np.save(os.path.join(save_path, 'cell_feature_exp.npy'), cell_dict)


def save_cell_graph(genes_path, save_path) :
  
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    exp = pd.read_csv(os.path.join(genes_path, 'EXP.csv'), index_col=0)
    mu = pd.read_csv(os.path.join(genes_path, 'MUT.csv'), index_col=0)
    
    cncat = "tcga" in genes_path
    if not cncat:
        # CNV in log2(CN-Ratio + 1)
        cn = pd.read_csv(os.path.join(genes_path, 'CNV.csv'), index_col=0)
    else:
        # CNV categorized for prediction of TCGA dataset
        cn = pd.read_csv(os.path.join(genes_path, 'CNV_Category.csv'), index_col=0)
    
    noise = "noise" in save_path
    if noise : 
        exp = add_rand_noise(exp)
        mu = add_rand_noise(mu)
        cn = add_rand_noise(cn)
    
    index = exp.index
    columns = exp.columns

    scaler_exp = StandardScaler()
    scaler_cnv = StandardScaler()
    
    exp = scaler_exp.fit_transform(exp)
    cn = scaler_cnv.fit_transform(cn)
    
    if genes_path == "./_data" and cncat:
        file_scaler_exp = os.path.join(save_path, 'scaler_exp.pickle')
        file_scaler_cnv = os.path.join(save_path, 'scaler_cncat.pickle')
        
        import pickle
        with open(file_scaler_exp, "wb") as f:
            pickle.dump(scaler_exp, f)
        with open(file_scaler_cnv, "wb") as f:
            pickle.dump(scaler_cnv, f)
    
    elif genes_path == "./_data" and not cncat:
        file_scaler_cnv = os.path.join(save_path, 'scaler_cn.pickle')
        
        import pickle
        with open(file_scaler_cnv, "wb") as f:
            pickle.dump(scaler_cnv, f)
    
    else: pass
    
    imp_mean = SimpleImputer()
    exp = imp_mean.fit_transform(exp)

    exp = pd.DataFrame(exp, index=index, columns=columns)
    cn = pd.DataFrame(cn, index=index, columns=columns)
    # mu = pd.DataFrame(mu, index=index, columns=columns)
    mu = pd.DataFrame(mu.values, index=index, columns=columns)
    
    cell_dict = {}
    cell_names = exp.index
    
    for i in cell_names:
        cell_dict[i] = Data(x=torch.tensor([exp.loc[i], cn.loc[i], mu.loc[i]], dtype=torch.float).T)
    
    if not cncat:
        np.save(os.path.join(save_path, 'cell_feature_all.npy'), cell_dict)
    else:
        np.save(os.path.join(save_path, 'cell_feature_all_cncat.npy'), cell_dict)


def get_STRING_graph(genes_path, thresh=0.95):
    save_path = os.path.join(genes_path, 'edge_index_PPI_{}.npy'.format(thresh))

    if not os.path.exists(save_path):
        # gene_list
        print("# [STRING] The file does not exists...")
        exp = pd.read_csv(os.path.join(genes_path, 'EXP.csv'), index_col=0)
        gene_list = exp.columns.to_list()
        gene_list = [int(gene) for gene in gene_list]
        # gene_list = [int(gene[1:-1]) for gene in gene_list]

        # load STRING
        ensp_map = ensp_to_hugo_map()
        hugo_map = hugo_to_ncbi_map()
        edges = pd.read_csv('./data/9606.protein.links.v11.0.txt', sep=' ')

        # edge_index
        selected_edges = edges['combined_score'] > (thresh * 1000)
        edge_list = edges[selected_edges][["protein1", "protein2"]].values.tolist()

        edge_list = [[ensp_map[edge[0]], ensp_map[edge[1]]] for edge in edge_list if
                     edge[0] in ensp_map.keys() and edge[1] in ensp_map.keys()]

        edge_list = [[hugo_map[edge[0]], hugo_map[edge[1]]] for edge in edge_list if
                     edge[0] in hugo_map.keys() and edge[1] in hugo_map.keys()]
        edge_index = []
        for i in edge_list:
            if (i[0] in gene_list) & (i[1] in gene_list):
                edge_index.append((gene_list.index(i[0]), gene_list.index(i[1])))
                edge_index.append((gene_list.index(i[1]), gene_list.index(i[0])))
        edge_index = list(set(edge_index))
        edge_index = np.array(edge_index, dtype=np.int64).T

        np.save(save_path, edge_index)
    else:
        print("# [STRING] The file already exists...")
        edge_index = np.load(save_path)

    return edge_index


def get_predefine_cluster(edge_index, genes_path, thresh):
    save_path = os.path.join(genes_path, 'cluster_predefine_PPI_{}.npy'.format(thresh))
    if not os.path.exists(save_path):
        print("# [Graclus] The file does not exists...")
        g = Data(edge_index=torch.tensor(edge_index, dtype=torch.long), x=torch.zeros(706, 1))
        g = Batch.from_data_list([g])
        cluster_predefine = {}
        for i in range(5):
            cluster = graclus(g.edge_index, None, g.x.size(0))
            print(len(cluster.unique()))
            g = max_pool(cluster, g, transform=None)
            cluster_predefine[i] = cluster
        np.save(save_path, cluster_predefine)
        # cluster_predefine = {i: j.to(device) for i, j in cluster_predefine.items()}
    else:
        print("# [Graclus] The file already exists...")
        cluster_predefine = np.load(save_path, allow_pickle=True).item()
        # cluster_predefine = {i: j.to(device) for i, j in cluster_predefine.items()}

    return cluster_predefine


if __name__ == '__main__':
    
    # genes_path = "./_data"
    # genes_path = "./_data_gdsc"
    # genes_path = "./_data_ccle"
    # genes_path = "./_data_tcga"
    # genes_path = "./_data_liu24_invivo"
    genes_path = "./_data_noise"
    save_path = genes_path
    
    # edge_index_850 = get_STRING_graph(genes_path, thresh=0.85)
    # edge_index_900 = get_STRING_graph(genes_path, thresh=0.90)
    # edge_index_950 = get_STRING_graph(genes_path, thresh=0.95)
    # 
    # cluster_850 = get_predefine_cluster(edge_index_850, genes_path, thresh=0.85)
    # cluster_900 = get_predefine_cluster(edge_index_900, genes_path, thresh=0.90)
    # cluster_950 = get_predefine_cluster(edge_index_950, genes_path, thresh=0.95)
    
    save_cell_graph(genes_path, save_path=save_path)
    
    
    # EXP for Similarity Graph
    if genes_path not in ["./_data_tcga", "./_data_liu24_invivo"]:
        file = "./{}/cell_feature_normalized".format(save_path)
        if os.path.exists(file):
            print("Normalized Feature Exist...")
        else:
            import pickle
            file1 = os.path.join(genes_path, 'EXP.csv')
            file2 = os.path.join(genes_path, 'EXP_Total_Scaled.csv')
            
            exp = pd.read_csv(file1, index_col=0)
            exp_scaled = pd.read_csv(file2, index_col=0)
            
            noise = "noise" in save_path
            if noise : 
                exp = add_rand_noise(exp)
                exp_scaled = add_rand_noise(exp_scaled)
            
            sum(exp_scaled.index!=exp.index)   # 0
            sum(exp_scaled.index==exp.index)   # 1399
            
            cell_feature_normalized = {k:v for k,v in zip(exp_scaled.index, exp_scaled.values)}
            with open(file, "wb") as f : 
                pickle.dump(cell_feature_normalized, f)
