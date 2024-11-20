import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, load_npz
import os
import preprocessing as pp
from joblib import Parallel, delayed, Memory
from tqdm import tqdm
import pickle
import re
from scipy.io import mmwrite
from preprocessing import load_pickle
import shutil

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


def split_list(lst, n):
    chunk_size = len(lst) // n
    remainder = len(lst) % n
    result = []
    start = 0
    for i in range(n):
        end = start + chunk_size + (1 if remainder > 0 else 0)
        result.append(lst[start:end])
        start = end
        remainder -= 1
    return result

class CountMatrix:
    def __init__(self, target:list, novel_read_n: int, group_novel = True, platform = '10x-ont', workers:int = 1,
                 csv = True, mtx =True, logger = None):
        self.logger = logger
        self.target = target
        self.workers = workers
        self.novel_read_n = novel_read_n
        self.platform = platform
        self.parse = self.platform == 'parse-ont'
        self.pacbio = self.platform == '10x-pacbio'
        self.group_novel = group_novel
        self.csv = csv
        self.mtx = mtx
        self.annotation_path_meta_gene_novel = os.path.join(target[0],"reference/metageneStructureInformationwNovel.pkl")
        self.novel_isoform_del_path = os.path.join(target[0],'reference/novel_isoform_del_' + str(novel_read_n) + '.pkl')
        if platform=='parse-ont':
            self.sample_names = os.listdir(os.path.join(self.target[0], 'samples'))
            self.n_samples = len(self.sample_names)
            self.samples_folder_path = os.path.join(self.target[0], 'samples')
            self.compatible_matrix_folder_path_list = [os.path.join(self.samples_folder_path, sample_name, 'compatible_matrix') for
                                         sample_name in self.sample_names]
            self.count_matrix_folder_path_list = [os.path.join(self.samples_folder_path, sample_name, 'count_matrix') for
                sample_name in self.sample_names]
            self.read_selection_pkl_path_list = [os.path.join(self.samples_folder_path, sample_name, 'auxillary/read_selection.pkl') for sample_name in self.sample_names]
        else:
            self.n_samples = len(target)
            self.compatible_matrix_folder_path_list = [os.path.join(target_, 'compatible_matrix') for target_ in target]
            self.count_matrix_folder_path_list = [os.path.join(target_, 'count_matrix') for target_ in target]
            self.read_selection_pkl_path_list = [os.path.join(target_, 'auxillary/read_selection.pkl') for target_ in target]
    def generate_count_matrix_by_gene(self, gene, read_selection_pkl):
        # CompatibleMatrixPaths = '/scr1/users/xu3/singlecell/project_singlecell/sample7_8_ont/sample7/compatible_matrix'
        # read_selection_pkl_paths = '/scr1/users/xu3/singlecell/project_singlecell/sample7_8_ont/sample7/auxillary/read_selection.pkl'
        # read_selection_pkl: keys must add sample index
        for folder in self.count_matrix_folder_path_list:
            os.makedirs(folder, exist_ok=True)
        pattern = re.compile(r'_ENS.+\.csv')
        files_with_indicators = [
            (os.path.join(CompatibleMatrixPath, i), idx)  # Tuple with file path and index
            for idx, CompatibleMatrixPath in enumerate(self.compatible_matrix_folder_path_list)
            for i in os.listdir(CompatibleMatrixPath)
            if '.csv' in i and pattern.sub('', i) == gene
        ]
        df_list = []
        for f, idx in files_with_indicators:
            df = pd.read_csv(f)
            df.columns = ['Cell'] + df.columns.tolist()[1:]
            if self.parse:
                df.Cell = df['Cell'].str.rsplit('_', n=1).str[0].tolist()
            df['Cell'] = df['Cell'] + f':sample{idx}'
            df = df[df['Cell'].isin([cell for cell in df['Cell'] if read_selection_pkl.get(cell) == 1])]  # filtering df
            if df.shape[0] > 0:
                df['Cell'] = df['Cell'].str.rsplit('_', n=1).str[0].tolist()
                df['Cell'] = df['Cell'] + f':sample{idx}'
                df = df.set_index('Cell')
                df_list.append(df)
        if len(df_list) == 0:
            return {gene: []}
        result_df = pd.concat(df_list, axis=0)
        df = result_df.fillna(0).astype(int)
        if self.group_novel:
            df, novel_isoform_name_mapping = pp.group_novel_isoform(df, geneStrand=self.annotation_pkl[gene][0]['geneStrand'],
                                                                    parse=self.parse)
        else:
            novel_isoform_name_mapping = None
        if novel_isoform_name_mapping is not None:
            novel_isoform_del = [key for key, value in novel_isoform_name_mapping.items() if key != value]
        else:
            novel_isoform_del = []
        # split df into samples
        df['sample_id'] = df.index.str.split(':').str[1]
        df.index = df.index.str.split(':').str[0]
        df_list = [group.drop(columns='sample_id') for _, group in df.groupby('sample_id')]
        # deal each sample separately
        for i, df in enumerate(df_list):
            # --------delete isoforms without reads
            df_isoform = df.sum(axis=0) > 0  # cols have read
            isoformNames = df_isoform[df_isoform].index.tolist()
            df = df.loc[:, isoformNames]
            # filter uncategorized reads
            df_uncategorized = pd.DataFrame()
            if df.shape[1] > 0:
                df_uncategorized = df.filter(['uncategorized'])
                if df_uncategorized.shape[1] > 0:
                    df = df.drop(columns=df_uncategorized.columns.tolist())
            if df.shape[1] > 0:
                # --------deal with multiple mappings
                df = deduplicate_col(df)  # delete same mapping isoforms
                # use unique mappings to decide multiple mappings
                multiple_bool = df.sum(axis=1) > 1
                multiple_index = [i for i, k in enumerate(multiple_bool.tolist()) if k]
                unique_index = [i for i, k in enumerate(multiple_bool.tolist()) if k == False]
                df_multiple = df.iloc[multiple_index, :]
                df_unique = df.iloc[unique_index, :]
                if df_multiple.shape[0] > 0:
                    if df_unique.shape[0] > 0:
                        column_percentages = df_unique.sum(axis=0).div(df_unique.sum(axis=0).sum(axis=0))
                    else:
                        column_percentages = df_multiple.sum(axis=0).div(df_multiple.sum(axis=0).sum(axis=0))
                    for ii in range(df_multiple.shape[0]):
                        isoforms_mapping = df_multiple.columns[(df_multiple.iloc[ii, :] == 1).values]
                        if column_percentages[isoforms_mapping].sum() == 0:
                            isoforms_mapping_max = column_percentages[isoforms_mapping].index[np.random.multinomial(1, [
                                1 / len(isoforms_mapping)] * len(isoforms_mapping)) == 1].tolist()
                        else:
                            isoforms_mapping_prob = column_percentages[isoforms_mapping] / column_percentages[
                                isoforms_mapping].sum()
                            isoforms_mapping_max = isoforms_mapping_prob.index[
                                np.random.multinomial(1, isoforms_mapping_prob) == 1].tolist()
                        for iso in list(isoforms_mapping):
                            if iso not in isoforms_mapping_max:
                                df_multiple.iloc[ii, :][iso] = 0
                    df = pd.concat([df_multiple, df_unique])
                if df_uncategorized.shape[1] > 0:
                    df_uncategorized = df_uncategorized.iloc[multiple_index + unique_index, :]
                    df = pd.concat([df, df_uncategorized], axis=1)
            else:
                df = df_uncategorized
            # filter novel isoform by supporting reads
            if df.shape[1] > 0:
                df_novel = df.filter(like='novel')
                if df_novel.shape[1] > 0:
                    novel_isoform_drop = df_novel.sum(axis=0) < self.novel_read_n
                    novel_isoform_drop = novel_isoform_drop[novel_isoform_drop].index.tolist()
                    df_drop = df.loc[:, novel_isoform_drop].sum(axis=1).tolist()
                    df = df.drop(columns=novel_isoform_drop)
                    novel_isoform_del = novel_isoform_del + novel_isoform_drop
                    df['uncategorized_novel'] = df_drop  # novel  isoforms but less than supporting reads
            if df.shape[1] == 0:
                return {gene: novel_isoform_del}
            df_all, df_filtered = df.copy(), df.copy()
            df_all.columns = [gene + '_' + iso for iso in df_all.columns.tolist()]
            df_filtered = df_filtered[df_filtered.columns[~df_filtered.columns.str.contains('uncategorized')]]
            df_filtered.columns = [gene + '_' + iso for iso in df_filtered.columns.tolist()]
            # confirm again------nonfiltered
            if df_all.shape[1] > 0:
                # aggregate by cells
                df_all = df_all.groupby(df_all.index).sum()
                df_gene = pd.DataFrame({'Cell': df_all.index.tolist(), gene: df_all.sum(axis=1).tolist()}).set_index(
                    'Cell')
                triple_transcript = df_to_triple(df_all)
                triple_gene = df_to_triple(df_gene)
                with open(os.path.join(self.count_matrix_folder_path_list[i], str(gene) + '_unfiltered_count.pickle'), 'wb') as f:
                    pickle.dump((triple_gene, triple_transcript), f)
        return {gene: novel_isoform_del}

    def generate_count_matrix_by_gene_list(self, gene_list, read_selection_pkl):
        novel_isoform_del_dict = {}
        for gene in gene_list:
            novel_isoform_del_dict_gene = self.generate_count_matrix_by_gene(gene, read_selection_pkl)
            novel_isoform_del_dict.update(novel_isoform_del_dict_gene)
        return novel_isoform_del_dict
    def read_filter(self):
        read_selection_pkl = {}
        for i, path in enumerate(self.read_selection_pkl_path_list):
            read_selection_pkl_ = pp.load_pickle(path)
            read_selection_pkl_updated = {key + f':sample{i}': value for key, value in read_selection_pkl_.items()}
            read_selection_pkl.update(read_selection_pkl_updated)
        return read_selection_pkl
    def generate_multiple_samples(self):
        pattern = re.compile(r'_ENS.+\.csv')
        Genes = []
        for p in self.compatible_matrix_folder_path_list:
            Genes_ = [g for g in os.listdir(p) if 'csv' in g]
            Genes_ = [pattern.sub('', g) for g in Genes_]
            Genes.append(Genes_)
        Genes = flatten_list(Genes)
        Genes = list(set(Genes))
        for count_path in self.count_matrix_folder_path_list:
            os.makedirs(count_path, exist_ok=True)
        adata_gene_unfiltered_list, adata_transcript_unfiltered_list = [], []
        annotation_pkl = None
        if self.group_novel:
            annotation_pkl = {}
            annotation_pkl_meta = pp.load_pickle(self.annotation_path_meta_gene_novel)
            metagenes = list(annotation_pkl_meta.keys())
            for metagene in metagenes:
                multi_gene_info = annotation_pkl_meta[metagene]
                for gene_info in multi_gene_info:
                    genename = re.sub(r'[\/\\\:\*\?\"\<\>\|]', '.', gene_info[0]['geneName'])
                    annotation_pkl[genename] = gene_info
        self.annotation_pkl = annotation_pkl
        self.logger.info(f'generating read filter')
        read_selection_pkl = self.read_filter()
        self.logger.info(f'generating count matrix pickles at: {self.count_matrix_folder_path_list}')
        Genes_list = split_list(Genes, self.workers)
        novel_isoform_del_dict = Parallel(n_jobs=self.workers)(delayed(self.generate_count_matrix_by_gene_list)(gene_list, read_selection_pkl) for gene_list in Genes_list)
        novel_isoform_del = {}
        for d in novel_isoform_del_dict:
            novel_isoform_del.update(d)
        self.logger.info('count matrix generated')
        #save novel_isoform_del
        self.novel_isoform_del_dict = novel_isoform_del
        with open(self.novel_isoform_del_path, 'wb') as f:
            pickle.dump(self.novel_isoform_del_dict, f)
        if len(self.target)>1:
            for additional_target in self.target[1:]:
                dest_path = os.path.join(additional_target, 'reference/novel_isoform_del_' + str(self.novel_read_n) + '.pkl')
                shutil.copyfile(self.novel_isoform_del_path, dest_path)
        for i, count_path in enumerate(self.count_matrix_folder_path_list):
            out_paths_unfiltered = [os.path.join(count_path, f) for f in os.listdir(count_path) if f.endswith('_unfiltered_count.pickle')]
            self.logger.info(f'reading {len(out_paths_unfiltered)} count pickles')
            Out_unfiltered = []
            for op in out_paths_unfiltered:
                Out_unfiltered.append(pp.load_pickle(op))
            self.logger.info('generating count matrix')
            triple_gene_list, triple_transcript_list = list(zip(*Out_unfiltered))
            triple_gene_list = pp.unpack_list(list(triple_gene_list))
            triple_transcript_list = pp.unpack_list(list(triple_transcript_list))
            adata_gene_unfiltered = generate_adata(triple_gene_list)
            adata_transcript_unfiltered = generate_adata(triple_transcript_list)
            self.logger.info('removing pickles')
            for op in out_paths_unfiltered:
                os.remove(op)
            adata_gene_unfiltered_list.append(adata_gene_unfiltered)
            adata_transcript_unfiltered_list.append(adata_transcript_unfiltered)
        self.adata_gene_unfiltered_list = adata_gene_unfiltered_list
        self.adata_transcript_unfiltered_list = adata_transcript_unfiltered_list
    def save_multiple_samples(self):
        if self.mtx:
            for i in range(self.n_samples):
                # save gene
                self.logger.info('saving count matrix in mtx format')
                print('saving count matrix on gene level ')
                gene_meta_unfiltered = {'obs': self.adata_gene_unfiltered_list[i].obs.index.tolist(),
                                        "var": self.adata_gene_unfiltered_list[i].var.index.tolist()}
                with open(os.path.join(self.count_matrix_folder_path_list[i],'adata_gene_unfiltered' + str(self.novel_read_n) + '.pickle'),'wb') as f:
                    pickle.dump(gene_meta_unfiltered, f)
                mmwrite(os.path.join(self.count_matrix_folder_path_list[i],'adata_gene_unfiltered' + str(self.novel_read_n)+ '.mtx'),self.adata_gene_unfiltered_list[i].X)
                # save transcript
                print('saving count matrix on transcript level ')
                transcript_meta_unfiltered = {'obs': self.adata_transcript_unfiltered_list[i].obs.index.tolist(),
                                              "var": self.adata_transcript_unfiltered_list[i].var.index.tolist()}
                mmwrite(os.path.join(self.count_matrix_folder_path_list[i], 'adata_transcript_unfiltered' + str(self.novel_read_n) + '.mtx'),
                        self.adata_transcript_unfiltered_list[i].X)
                with open(os.path.join(self.count_matrix_folder_path_list[i],
                                       'adata_transcript_unfiltered' + str(self.novel_read_n) + '.pickle'),'wb') as f:
                    pickle.dump(transcript_meta_unfiltered, f)
        if self.csv:
            self.logger.info('saving count matrix in csv format')
            for i in range(self.n_samples):
                output_gene_unfiltered = os.path.join(self.count_matrix_folder_path_list[i],'adata_gene_unfiltered' + str(self.novel_read_n) + '.csv')
                output_transcript_unfiltered = os.path.join(self.count_matrix_folder_path_list[i], 'adata_transcript_unfiltered' + str(self.novel_read_n) + '.csv')
                # save gene
                print('saving count matrix on gene level ')
                adata_gene_unfiltered_df = self.adata_gene_unfiltered_list[i].to_df()
                adata_gene_unfiltered_df.to_csv(output_gene_unfiltered)
                # save transcript
                print('saving count matrix on transcript level ')
                adata_transcript_unfiltered_df = self.adata_transcript_unfiltered_list[i].to_df()
                adata_transcript_unfiltered_df.to_csv(output_transcript_unfiltered)
    def _extract_attribute(self, attributes, attribute_name):
        try:
            return attributes.split(f'{attribute_name} "')[1].split('"')[0]
        except IndexError:
            return None
    def filter_gtf(self):
        self.logger.info('updating gtf annotation file')
        for target in self.target:
            input_gtf = os.path.join(target, 'reference/SCOTCH_updated_annotation.gtf')
            output_gtf = os.path.join(target, 'reference/SCOTCH_updated_annotation_filtered.gtf')
            with open(input_gtf, "r") as infile, open(output_gtf, "w") as outfile:
                for line in infile:
                    if line.startswith("#"):
                        outfile.write(line)
                        continue
                    columns = line.split("\t")
                    attributes = columns[8]
                    gene_name = self._extract_attribute(attributes, 'gene_name')
                    transcript_id = self._extract_attribute(attributes, 'transcript_id')
                    if gene_name in self.novel_isoform_del_dict and transcript_id in self.novel_isoform_del_dict[gene_name]:
                        continue
                    outfile.write(line)
