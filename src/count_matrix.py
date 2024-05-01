import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, load_npz
import os
import preprocessing as pp
from joblib import Parallel, delayed
from tqdm import tqdm
import pickle
import re
from preprocessing import load_pickle

#bam_info='/Volumes/G-DRIVE USB-C/singlecell_file/M4/M4.filtered.bam.Info.csv'
#CompatibleMatrixPath='/Volumes/G-DRIVE USB-C/singlecell_file/LH/compatible_matrix'

#CompatibleMatrixPath = '/scr1/users/xu3/singlecell/project_singlecell/sample8_R9/compatible_matrix'
#gene_pkl='/scr1/users/xu3/singlecell/project_singlecell/LH/reference/geneStructureInformation.pkl'

def generate_read_df(f, geneStructureInformation):
    df = pd.read_csv(f)
    f = os.path.basename(f)
    gene = f[:-4]
    geneID = gene[-15:]
    geneName = gene[:-16]
    df.columns = ['Reads'] + df.columns.tolist()[1:]
    df = df.set_index('Reads')
    triple = df_to_triple(df)
    reads_type = []
    for tri in triple:
        _, reads, type = tri
        if 'novel' in type:
            type = 'novel'
        if 'ENST' in type:
            type = 'existing'
        reads_type.append((reads, type))
    read_df = pd.DataFrame(reads_type)
    read_df = read_df.drop_duplicates()
    read_df.columns = ['reads', 'type']
    read_df['gene'] = gene
    read_df['geneName'] = geneName
    read_df['geneID'] = geneID
    read_df['geneLength'] = geneStructureInformation[geneID][0]['geneEnd'] - geneStructureInformation[geneID][0][
        'geneStart']
    return read_df


#filter reads: filter reads mapped to multiple genes
def filter_reads(CompatibleMatrixPath, gene_pkl, num_cores):
    if isinstance(gene_pkl,str):
        geneStructureInformation = load_pickle(gene_pkl)
    else:
        geneStructureInformation = gene_pkl
    files = [os.path.join(CompatibleMatrixPath,f) for f in os.listdir(CompatibleMatrixPath) if '.csv' in f]
    DF = Parallel(n_jobs=num_cores)(delayed(generate_read_df)(f, geneStructureInformation) for f in files)
    DF = pd.concat(DF)
    DF = DF.reset_index(drop = True)
    type_priority_map = {'existing': 1,'novel': 2,'uncategorized': 3,'intronReads': 4}
    DF['priority'] = DF['type'].map(type_priority_map)
    #determin if a read maps to multiple genes
    DF['multimap'] = DF.groupby('reads')['geneID'].transform('nunique').gt(1).astype(int)
    DF['drop'] = 0
    keep_gene = DF.sort_values(by=['reads', 'priority', 'geneLength'], ascending=[True, True, False]).groupby(
        'reads').first().reset_index()
    DF = DF.merge(keep_gene[['reads', 'gene']], on='reads', how='left', suffixes=('', '_to_keep'))
    DF['drop'] = (DF['multimap'] == 1) & (DF['gene'] != DF['gene_to_keep']).astype(int)
    DF = DF.drop(columns=['priority', 'gene_to_keep'])
    return DF



def df_to_triple(df):
    rowNames = df.index.tolist()
    colNmes = df.columns.tolist()
    mat = np.array(df)
    x = mat[np.nonzero(mat)]
    rowIndex, colIndex = np.nonzero(mat)
    ij_list = list(zip(x, rowIndex,colIndex))
    triple = [(x, rowNames[i], colNmes[j]) for x, i, j in ij_list]
    return triple


def flatten_list(nested_list):
    flattened = []
    for item in nested_list:
        if isinstance(item, list):
            flattened.extend(flatten_list(item))
        else:
            flattened.append(item)
    return flattened

def deduplicate_col(df):
    if df.shape[1]<2:
        return df
    else:
        arr = df.to_numpy()
        r_ind = np.where(arr > 0)[0]
        c_ind = np.where(arr > 0)[1]
        ind_dict = {}
        for i in range(arr.shape[1]):
            r_ind_ = r_ind[np.where(c_ind == i)[0].tolist()].tolist()
            ind_dict[i] = r_ind_
        duplicated_cols = []
        for i in range(arr.shape[1] - 1):
            for j in range(i + 1, arr.shape[1]):
                if len(ind_dict[i]) == len(ind_dict[j]):
                    diff = np.array(ind_dict[i]) - np.array(ind_dict[j])
                    if sum(diff != 0) == 0:
                        duplicated_cols.append(j)
        if len(duplicated_cols)>0:
            df = df.drop(df.columns[duplicated_cols], axis=1)
        return df

def generate_count_matrix_by_gene(CompatibleMatrixPath, gene, novel_read_n = 10, output_folder=None):
    #CompatibleMatrixPath: path to compatible matrix directory
    #CompatibleMatrixPath = '/scr1/users/xu3/singlecell/project_singlecell/LH/compatible_matrix'
    #cells
    pattern = re.compile(r'_ENS.+\.csv')
    files = [i for i in os.listdir(CompatibleMatrixPath) if '.csv' in i and pattern.sub('',i)==gene]
    df_list = []
    for f in files:
        df = pd.read_csv(os.path.join(CompatibleMatrixPath, f))
        df.columns = ['Cell'] + df.columns.tolist()[1:]
        df.Cell = df['Cell'].str.split('_', expand=True).iloc[:, 0].tolist()  # change cell names
        df = df.set_index('Cell')
        df_list.append(df)
    if len(df_list)==0:
        return
    result_df = pd.concat(df_list, axis=0)
    df = result_df.fillna(0).astype(int)
    #--------delete isoforms without reads
    df_isoform = df.sum(axis=0) > 0 #cols have read
    isoformNames = df_isoform[df_isoform].index.tolist()
    df = df.loc[:, isoformNames]
    #filter uncategorized reads
    df_uncategorized = pd.DataFrame()
    if df.shape[1] > 0:
        df_uncategorized = df.filter(['uncategorized'])
        if df_uncategorized.shape[1] > 0:
            df = df.drop(columns=df_uncategorized.columns.tolist())
    if df.shape[1] > 0:
        #--------deal with multiple mappings
        df = deduplicate_col(df) # delete same mapping isoforms
        # use unique mappings to decide multiple mappings
        multiple_bool = df.sum(axis=1) > 1
        multiple_index = [i for i,k in enumerate(multiple_bool.tolist()) if k]
        unique_index = [i for i, k in enumerate(multiple_bool.tolist()) if k==False]
        df_multiple = df.iloc[multiple_index,:]
        df_unique = df.iloc[unique_index, :]
        if df_multiple.shape[0] > 0:
            if df_unique.shape[0]>0:
                column_percentages = df_unique.sum(axis=0).div(df_unique.sum(axis=0).sum(axis=0))
            else:
                column_percentages = df_multiple.sum(axis=0).div(df_multiple.sum(axis=0).sum(axis=0))
            for ii in range(df_multiple.shape[0]):
                isoforms_mapping = df_multiple.columns[(df_multiple.iloc[ii, :] == 1).values]
                if column_percentages[isoforms_mapping].sum() == 0:
                    isoforms_mapping_max = column_percentages[isoforms_mapping].index[np.random.multinomial(1,[1/len(isoforms_mapping)]*len(isoforms_mapping))==1].tolist()
                else:
                    isoforms_mapping_prob = column_percentages[isoforms_mapping]/column_percentages[isoforms_mapping].sum()
                    isoforms_mapping_max = isoforms_mapping_prob.index[np.random.multinomial(1,isoforms_mapping_prob)==1].tolist()
                for iso in list(isoforms_mapping):
                    if iso not in isoforms_mapping_max:
                        df_multiple.iloc[ii,:][iso] = 0
            df = pd.concat([df_multiple, df_unique])
        if df_uncategorized.shape[1]>0:
            df_uncategorized = df_uncategorized.iloc[multiple_index+unique_index,:]
            df = pd.concat([df, df_uncategorized],axis=1)
    else:
        df = df_uncategorized
    # filter novel isoform by supporting reads
    if df.shape[1]>0:
        df_novel = df.filter(like='novel')
        if df_novel.shape[1]>0:
            novel_isoform_drop = df_novel.sum(axis=0)<novel_read_n
            novel_isoform_drop = novel_isoform_drop[novel_isoform_drop].index.tolist()
            df_drop=df.loc[:, novel_isoform_drop].sum(axis=1).tolist()
            df = df.drop(columns = novel_isoform_drop)
            df['uncategorized_novel']=df_drop #novel  isoforms but less than supporting reads
    if df.shape[1] ==0:
        return
    df_all, df_filtered = df.copy(), df.copy()
    df_all.columns = [gene + '_' + iso for iso in df_all.columns.tolist()]
    df_filtered = df_filtered[df_filtered.columns[~df_filtered.columns.str.contains('uncategorized')]]
    df_filtered.columns = [gene + '_' + iso for iso in df_filtered.columns.tolist()]
    #confirm again------nonfiltered
    if df_all.shape[1] > 0:
        #aggregate by cells
        df_all = df_all.groupby(df_all.index).sum()
        df_gene = pd.DataFrame({'Cell':df_all.index.tolist(),gene:df_all.sum(axis=1).tolist()}).set_index('Cell')
        triple_transcript = df_to_triple(df_all)
        triple_gene = df_to_triple(df_gene)
        if output_folder is None:
            return (triple_gene, triple_transcript)
        else:
            with open(os.path.join(output_folder,str(gene)+'_unfiltered_count.pickle'),'wb') as f:
                pickle.dump((triple_gene, triple_transcript), f)
    if df_filtered.shape[1]>0:
        df_filtered = df_filtered.groupby(df_filtered.index).sum()
        df_gene = pd.DataFrame({'Cell': df_filtered.index.tolist(), gene: df_filtered.sum(axis=1).tolist()}).set_index('Cell')
        triple_transcript = df_to_triple(df_filtered)
        triple_gene = df_to_triple(df_gene)
        if output_folder is None:
            return (triple_gene, triple_transcript)
        else:
            with open(os.path.join(output_folder,str(gene)+'_filtered_count.pickle'),'wb') as f:
                pickle.dump((triple_gene, triple_transcript), f)
    else:
        if output_folder is not None:
            log_file = os.path.join(output_folder, 'log.txt')
            if os.path.isfile(log_file):
                with open(log_file, 'a') as file:
                    file.write(str(gene) + '\n')
            else:
                with open(log_file, 'w') as file:
                    file.write(str(gene) + '\n')
        return None

def generate_count_matrix_by_gene_parse(library_path_list, sample_name, gene, novel_read_n = 10, output_folder=None):
    #CompatibleMatrixPath: path to compatible matrix directory
    #CompatibleMatrixPath = '/scr1/users/xu3/singlecell/project_singlecell/LH/compatible_matrix'
    #cells
    CompatibleMatrixPath_list = [os.path.join(p, sample_name, 'compatible_matrix') for p in library_path_list]
    pattern = re.compile(r'_ENS.+\.csv')
    files_list=[]
    for comp_path in CompatibleMatrixPath_list:
        files = [i for i in os.listdir(comp_path) if '.csv' in i and pattern.sub('', i) == gene]
        files = [os.path.join(comp_path, f) for f in files]
        files_list.append(files)
    files = flatten_list(files_list)
    df_list = []
    for f in files:
        df = pd.read_csv(f)
        df.columns = ['Cell'] + df.columns.tolist()[1:]
        df.Cell = df['Cell'].str.rsplit('_', n=1).str[0].tolist()  # change cell names
        df = df.set_index('Cell')
        df_list.append(df)
    if len(df_list)==0:
        return
    result_df = pd.concat(df_list, axis=0)
    df = result_df.fillna(0).astype(int)
    #--------delete isoforms without reads
    df_isoform = df.sum(axis=0) > 0 #cols have read
    isoformNames = df_isoform[df_isoform].index.tolist()
    df = df.loc[:, isoformNames]
    #filter uncategorized reads
    df_uncategorized = pd.DataFrame()
    if df.shape[1] > 0:
        df_uncategorized = df.filter(['uncategorized'])
        if df_uncategorized.shape[1] > 0:
            df = df.drop(columns=df_uncategorized.columns.tolist())
    if df.shape[1] > 0:
        #--------deal with multiple mappings
        df = deduplicate_col(df) # delete same mapping isoforms
        # use unique mappings to decide multiple mappings
        multiple_bool = df.sum(axis=1) > 1
        multiple_index = [i for i,k in enumerate(multiple_bool.tolist()) if k]
        unique_index = [i for i, k in enumerate(multiple_bool.tolist()) if k==False]
        df_multiple = df.iloc[multiple_index,:]
        df_unique = df.iloc[unique_index, :]
        if df_multiple.shape[0] > 0:
            if df_unique.shape[0]>0:
                column_percentages = df_unique.sum(axis=0).div(df_unique.sum(axis=0).sum(axis=0))
            else:
                column_percentages = df_multiple.sum(axis=0).div(df_multiple.sum(axis=0).sum(axis=0))
            for ii in range(df_multiple.shape[0]):
                isoforms_mapping = df_multiple.columns[(df_multiple.iloc[ii, :] == 1).values]
                if column_percentages[isoforms_mapping].sum() == 0:
                    isoforms_mapping_max = column_percentages[isoforms_mapping].index[np.random.multinomial(1,[1/len(isoforms_mapping)]*len(isoforms_mapping))==1].tolist()
                else:
                    isoforms_mapping_prob = column_percentages[isoforms_mapping]/column_percentages[isoforms_mapping].sum()
                    isoforms_mapping_max = isoforms_mapping_prob.index[np.random.multinomial(1,isoforms_mapping_prob)==1].tolist()
                for iso in list(isoforms_mapping):
                    if iso not in isoforms_mapping_max:
                        df_multiple.iloc[ii,:][iso] = 0
            df = pd.concat([df_multiple, df_unique])
        if df_uncategorized.shape[1]>0:
            df_uncategorized = df_uncategorized.iloc[multiple_index+unique_index,:]
            df = pd.concat([df, df_uncategorized],axis=1)
    else:
        df = df_uncategorized
    # filter novel isoform by supporting reads
    if df.shape[1]>0:
        df_novel = df.filter(like='novel')
        if df_novel.shape[1]>0:
            novel_isoform_drop = df_novel.sum(axis=0)<novel_read_n
            novel_isoform_drop = novel_isoform_drop[novel_isoform_drop].index.tolist()
            df_drop=df.loc[:, novel_isoform_drop].sum(axis=1).tolist()
            df = df.drop(columns = novel_isoform_drop)
            df['uncategorized_novel']=df_drop #novel  isoforms but less than supporting reads
    if df.shape[1] ==0:
        return
    df_all, df_filtered = df.copy(), df.copy()
    df_all.columns = [gene + '_' + iso for iso in df_all.columns.tolist()]
    df_filtered = df_filtered[df_filtered.columns[~df_filtered.columns.str.contains('uncategorized')]]
    df_filtered.columns = [gene + '_' + iso for iso in df_filtered.columns.tolist()]
    #confirm again------nonfiltered
    if df_all.shape[1] > 0:
        #aggregate by cells
        df_all = df_all.groupby(df_all.index).sum()
        df_gene = pd.DataFrame({'Cell':df_all.index.tolist(),gene:df_all.sum(axis=1).tolist()}).set_index('Cell')
        triple_transcript = df_to_triple(df_all)
        triple_gene = df_to_triple(df_gene)
        if output_folder is None:
            return (triple_gene, triple_transcript)
        else:
            with open(os.path.join(output_folder,str(gene)+'_unfiltered_count.pickle'),'wb') as f:
                pickle.dump((triple_gene, triple_transcript), f)
    if df_filtered.shape[1]>0:
        df_filtered = df_filtered.groupby(df_filtered.index).sum()
        df_gene = pd.DataFrame({'Cell': df_filtered.index.tolist(), gene: df_filtered.sum(axis=1).tolist()}).set_index('Cell')
        triple_transcript = df_to_triple(df_filtered)
        triple_gene = df_to_triple(df_gene)
        if output_folder is None:
            return (triple_gene, triple_transcript)
        else:
            with open(os.path.join(output_folder,str(gene)+'_filtered_count.pickle'),'wb') as f:
                pickle.dump((triple_gene, triple_transcript), f)
    else:
        if output_folder is not None:
            log_file = os.path.join(output_folder, 'log.txt')
            if os.path.isfile(log_file):
                with open(log_file, 'a') as file:
                    file.write(str(gene) + '\n')
            else:
                with open(log_file, 'w') as file:
                    file.write(str(gene) + '\n')
        return None
def generate_adata(triple_list):
    cells_dict = {}
    features_dict = {}
    data = []
    cells = []
    features = []
    for x, cell, feature in tqdm(triple_list):
        if cell not in cells_dict:
            cells_dict[cell] = len(cells_dict)
            cells.append(cell)
        if feature not in features_dict:
            features_dict[feature] = len(features_dict)
            features.append(feature)
        data.append((x, cells_dict[cell], features_dict[feature]))
    x, cells_ind, features_ind = zip(*data)
    sparse_matrix = csr_matrix((x, (cells_ind, features_ind)))
    adata = ad.AnnData(sparse_matrix)
    adata.obs_names = cells
    adata.var_names = features
    return adata


def read_adata(npz_path):
    mat = load_npz(npz_path)
    pickle_path = npz_path[:-4]+'.pickle'
    meta = pp.load_pickle(pickle_path)
    adata = ad.AnnData(mat)
    adata.obs_names = meta['obs']
    adata.var_names = meta['var']
    return adata


def csv_to_adata(csv_path):
    count = pd.read_csv(csv_path)
    cells = count.iloc[:, 0].tolist()
    features = count.columns.tolist()[1:]
    mat = np.array(count.iloc[:, 1:])
    mat = csr_matrix(mat)
    adata = ad.AnnData(mat)
    adata.obs_names = cells
    adata.var_names = features
    return adata

#main function
def generate_count_matrix(path='/scr1/users/xu3/singlecell/project_singlecell/',novel_read_n=1, num_cores = 1):
    CompatibleMatrixPath = os.path.join(path, 'compatible_matrix')
    pattern = re.compile(r'_ENS.+\.csv')
    Genes = [g for g in os.listdir(CompatibleMatrixPath) if 'csv' in g]
    Genes = [pattern.sub('', g) for g in Genes]
    count_path = os.path.join(path, 'count_matrix')
    print('generating count matrix pickles')
    Parallel(n_jobs=num_cores)(delayed(generate_count_matrix_by_gene)(CompatibleMatrixPath, gene, novel_read_n, count_path)
                               for gene in tqdm(Genes))
    out_paths_unfiltered = [os.path.join(count_path, f) for f in os.listdir(count_path) if
                          f.endswith('_unfiltered_count.pickle')]
    out_paths_filtered = [os.path.join(count_path, f) for f in os.listdir(count_path) if
                          f.endswith('_filtered_count.pickle')]
    print('reading pickles for unfiltered counts')
    Out_unfiltered = []
    for op in out_paths_unfiltered:
        Out_unfiltered.append(pp.load_pickle(op))
    print('reading pickles for filtered counts')
    Out_filtered = []
    for op in out_paths_filtered:
        Out_filtered.append(pp.load_pickle(op))
    print('generating count matrix')
    triple_gene_list, triple_transcript_list = list(zip(*Out_unfiltered))
    triple_gene_list = pp.unpack_list(list(triple_gene_list))
    triple_transcript_list = pp.unpack_list(list(triple_transcript_list))
    adata_gene_unfiltered = generate_adata(triple_gene_list)
    adata_transcript_unfiltered = generate_adata(triple_transcript_list)
    triple_gene_list, triple_transcript_list = list(zip(*Out_filtered))
    triple_gene_list = pp.unpack_list(list(triple_gene_list))
    triple_transcript_list = pp.unpack_list(list(triple_transcript_list))
    adata_gene_filtered = generate_adata(triple_gene_list)
    adata_transcript_filtered = generate_adata(triple_transcript_list)
    for op in out_paths_unfiltered:
        os.remove(op)
    for op in out_paths_filtered:
        os.remove(op)
    return adata_gene_unfiltered, adata_transcript_unfiltered, adata_gene_filtered, adata_transcript_filtered

def generate_count_matrix_parse(path='/scr1/users/xu3/singlecell/project_singlecell/',novel_read_n=1, num_cores = 1,
                                key_words='parse_sub', sample_name='Sample7'):
    library_path_list = [os.path.join(path,f) for f in os.listdir(path) if key_words in f]
    pattern = re.compile(r'_ENS.+\.csv')
    CompatibleMatrixPath_list = [os.path.join(p, sample_name, 'compatible_matrix') for p in library_path_list]
    Genes=[]
    for p in CompatibleMatrixPath_list:
        Genes_ = [g for g in os.listdir(p) if 'csv' in g]
        Genes_ = [pattern.sub('', g) for g in Genes_]
        Genes.append(Genes_)
    Genes=flatten_list(Genes)
    Genes=list(set(Genes))
    count_path = os.path.join(path, 'parse_count_matrix_'+str(sample_name))
    if os.path.exists(count_path) is False:
        os.makedirs(count_path)
    print('output folder is: '+str(count_path))
    print('generating count matrix pickles')
    Parallel(n_jobs=num_cores)(delayed(generate_count_matrix_by_gene_parse)(library_path_list,sample_name, gene, novel_read_n, count_path)
                               for gene in tqdm(Genes))
    out_paths_unfiltered = [os.path.join(count_path, f) for f in os.listdir(count_path) if
                          f.endswith('_unfiltered_count.pickle')]
    out_paths_filtered = [os.path.join(count_path, f) for f in os.listdir(count_path) if
                          f.endswith('_filtered_count.pickle')]
    print('reading pickles for unfiltered counts')
    Out_unfiltered = []
    for op in out_paths_unfiltered:
        Out_unfiltered.append(pp.load_pickle(op))
    print('reading pickles for filtered counts')
    Out_filtered = []
    for op in out_paths_filtered:
        Out_filtered.append(pp.load_pickle(op))
    print('generating count matrix')
    triple_gene_list, triple_transcript_list = list(zip(*Out_unfiltered))
    triple_gene_list = pp.unpack_list(list(triple_gene_list))
    triple_transcript_list = pp.unpack_list(list(triple_transcript_list))
    adata_gene_unfiltered = generate_adata(triple_gene_list)
    adata_transcript_unfiltered = generate_adata(triple_transcript_list)
    triple_gene_list, triple_transcript_list = list(zip(*Out_filtered))
    triple_gene_list = pp.unpack_list(list(triple_gene_list))
    triple_transcript_list = pp.unpack_list(list(triple_transcript_list))
    adata_gene_filtered = generate_adata(triple_gene_list)
    adata_transcript_filtered = generate_adata(triple_transcript_list)
    for op in out_paths_unfiltered:
        os.remove(op)
    for op in out_paths_filtered:
        os.remove(op)
    return adata_gene_unfiltered, adata_transcript_unfiltered, adata_gene_filtered, adata_transcript_filtered