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
from scipy.io import mmwrite, mmread
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
        return df, []
    else:
        arr = df.to_numpy()
        r_ind, c_ind = np.where(arr > 0)
        ind_dict = {}
        for i in range(arr.shape[1]):
            r_ind_ = r_ind[np.where(c_ind == i)[0].tolist()].tolist()
            ind_dict[i] = r_ind_
        duplicated_cols = []
        for i in range(arr.shape[1] - 1):
            for j in range(i + 1, arr.shape[1]):
                if ind_dict[i] == ind_dict[j]:
                    duplicated_cols.append(j)
        dropped_col_names = []
        if duplicated_cols:
            dropped_col_names = df.columns[duplicated_cols].tolist()
            df = df.drop(columns=dropped_col_names)
        dropped_col_names_novel = [cn for cn in dropped_col_names if cn.startswith('novel')]
        return df, dropped_col_names_novel



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
    def __init__(self, target:list, novel_read_n: int, novel_read_pct: float = 0,
                 group_novel = True, platform = '10x-ont', workers:int = 1,
                 csv = True, mtx =True, logger = None, gene_subset = None):
        self.logger = logger
        self.target = target
        self.workers = workers
        self.novel_read_n = novel_read_n
        self.novel_read_pct = novel_read_pct
        self.platform = platform
        self.parse = self.platform == 'parse-ont'
        self.pacbio = self.platform == '10x-pacbio'
        self.group_novel = group_novel
        self.novel_name_substitution = []
        self.csv = csv
        self.mtx = mtx
        self.gene_subset = gene_subset
        self.annotation_path_meta_gene_novel = os.path.join(target[0],"reference/metageneStructureInformationwNovel.pkl")
        self.novel_isoform_del_path = os.path.join(target[0],f'reference/novel_isoform_del_{str(self.novel_read_n)}_{str(self.novel_read_pct)}.pkl')
        self.novel_name_substitution_path = os.path.join(target[0],'reference/novel_name_substitutions.pkl')
        if platform=='parse-ont':
            self.sample_names = os.listdir(os.path.join(self.target[0], 'samples'))
            self.n_samples = len(self.sample_names)
            self.samples_folder_path = os.path.join(self.target[0], 'samples')
            self.compatible_matrix_folder_path_list = [os.path.join(self.samples_folder_path, sample_name, 'compatible_matrix') for
                                         sample_name in self.sample_names]
            self.count_matrix_folder_path_list = [os.path.join(self.samples_folder_path, sample_name, 'count_matrix') for
                sample_name in self.sample_names]
            self.read_selection_pkl_path_list = [os.path.join(self.samples_folder_path, sample_name, 'auxillary/read_selection.pkl') for sample_name in self.sample_names]
            self.spliced_compatible_matrix_folder_path_list = [
                os.path.join(self.samples_folder_path, sample_name, 'spliced_compatible_matrix') for
                sample_name in self.sample_names]
            self.unspliced_compatible_matrix_folder_path_list = [
                os.path.join(self.samples_folder_path, sample_name, 'unspliced_compatible_matrix') for
                sample_name in self.sample_names]
            self.count_matrix_spliced_folder_path_list = [os.path.join(self.samples_folder_path, sample_name, 'count_matrix', 'spliced')
                                                  for
                                                  sample_name in self.sample_names]
            self.count_matrix_unspliced_folder_path_list = [os.path.join(self.samples_folder_path, sample_name, 'count_matrix', 'unspliced')
                                                  for
                                                  sample_name in self.sample_names]

        else:
            self.n_samples = len(target)
            self.compatible_matrix_folder_path_list = [os.path.join(target_, 'compatible_matrix') for target_ in target]
            self.count_matrix_folder_path_list = [os.path.join(target_, 'count_matrix') for target_ in target]
            self.read_selection_pkl_path_list = [os.path.join(target_, 'auxillary/read_selection.pkl') for target_ in target]
            self.spliced_compatible_matrix_folder_path_list = [os.path.join(target_, 'spliced_compatible_matrix') for target_ in target]
            self.unspliced_compatible_matrix_folder_path_list = [os.path.join(target_, 'unspliced_compatible_matrix') for target_ in target]
            self.count_matrix_spliced_folder_path_list = [os.path.join(target_, 'count_matrix', 'spliced') for target_ in target]
            self.count_matrix_unspliced_folder_path_list = [os.path.join(target_, 'count_matrix', 'unspliced') for target_ in target]

    def _load_annotation_pkl(self):
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

    def _get_count_output_paths(self, folder_path, level, splicing=None):
        if splicing is None:
            base_name = f'adata_{level}_{self.novel_read_n}_{self.novel_read_pct}'
        else:
            base_name = f'adata_{level}_unfiltered{self.novel_read_n}'
        return {
            'csv': os.path.join(folder_path, base_name + '.csv'),
            'mtx': os.path.join(folder_path, base_name + '.mtx'),
            'pickle': os.path.join(folder_path, base_name + '.pickle')
        }

    def _load_saved_matrix_df(self, folder_path, level, splicing=None):
        paths = self._get_count_output_paths(folder_path, level, splicing=splicing)
        if os.path.exists(paths['csv']):
            return pd.read_csv(paths['csv'], index_col=0)
        if os.path.exists(paths['mtx']) and os.path.exists(paths['pickle']):
            with open(paths['pickle'], 'rb') as handle:
                meta = pickle.load(handle)
            matrix = mmread(paths['mtx']).tocsr()
            return pd.DataFrame.sparse.from_spmatrix(matrix, index=meta['obs'], columns=meta['var'])
        return pd.DataFrame()

    def _save_matrix_df(self, df, folder_path, level, splicing=None):
        paths = self._get_count_output_paths(folder_path, level, splicing=splicing)
        save_csv = self.csv or os.path.exists(paths['csv'])
        save_mtx = self.mtx or os.path.exists(paths['mtx']) or os.path.exists(paths['pickle'])
        dense_df = df if isinstance(df, pd.DataFrame) else pd.DataFrame(df)
        if save_csv:
            dense_df.to_csv(paths['csv'])
        if save_mtx:
            with open(paths['pickle'], 'wb') as handle:
                pickle.dump({'obs': dense_df.index.tolist(), 'var': dense_df.columns.tolist()}, handle)
            if any(isinstance(dtype, pd.SparseDtype) for dtype in dense_df.dtypes):
                matrix = dense_df.sparse.to_coo().tocsr()
            else:
                matrix = csr_matrix(dense_df.to_numpy())
            mmwrite(paths['mtx'], matrix)

    def _load_subset_pickles_df(self, folder_path, subset_genes):
        suffix = '_unfiltered_count.pickle'
        out_paths = [
            os.path.join(folder_path, f) for f in os.listdir(folder_path)
            if f.endswith(suffix) and f[:-len(suffix)] in subset_genes
        ]
        if len(out_paths) == 0:
            return pd.DataFrame(), pd.DataFrame(), []
        out_unfiltered = [pp.load_pickle(path) for path in out_paths]
        triple_gene_list, triple_transcript_list = list(zip(*out_unfiltered))
        triple_gene_list = pp.unpack_list(list(triple_gene_list))
        triple_transcript_list = pp.unpack_list(list(triple_transcript_list))
        gene_df = generate_adata(triple_gene_list).to_df() if len(triple_gene_list) > 0 else pd.DataFrame()
        transcript_df = generate_adata(triple_transcript_list).to_df() if len(triple_transcript_list) > 0 else pd.DataFrame()
        return gene_df, transcript_df, out_paths

    def _splice_matrix_df(self, existing_df, subset_df, subset_genes, level):
        if level == 'gene':
            base_df = existing_df.drop(columns=[col for col in existing_df.columns if col in subset_genes], errors='ignore')
        else:
            subset_prefixes = tuple(f'{gene}_' for gene in subset_genes)
            base_df = existing_df.loc[:, [col for col in existing_df.columns if not col.startswith(subset_prefixes)]]
        if base_df.empty and subset_df.empty:
            return pd.DataFrame()
        all_index = base_df.index.union(subset_df.index)
        base_df = base_df.reindex(all_index, fill_value=0)
        subset_df = subset_df.reindex(all_index, fill_value=0)
        merged_df = pd.concat([base_df, subset_df], axis=1)
        if merged_df.shape[1] > 0:
            merged_df = merged_df.loc[:, ~merged_df.columns.duplicated(keep='last')]
        return merged_df.sort_index()

    def _merge_novel_metadata(self, subset_novel_isoform_del, subset_novel_name_substitution):
        if os.path.exists(self.novel_isoform_del_path):
            novel_isoform_del = pp.load_pickle(self.novel_isoform_del_path)
        else:
            novel_isoform_del = {}
        if os.path.exists(self.novel_name_substitution_path):
            novel_name_substitution = pp.load_pickle(self.novel_name_substitution_path)
        else:
            novel_name_substitution = {}
        for gene in self.gene_subset:
            novel_isoform_del.pop(gene, None)
            novel_name_substitution.pop(gene, None)
        novel_isoform_del.update(subset_novel_isoform_del)
        novel_name_substitution.update(subset_novel_name_substitution)
        self.novel_isoform_del_dict = novel_isoform_del
        self.novel_name_substitution_dict = novel_name_substitution
        with open(self.novel_isoform_del_path, 'wb') as handle:
            pickle.dump(self.novel_isoform_del_dict, handle)
        with open(self.novel_name_substitution_path, 'wb') as handle:
            pickle.dump(self.novel_name_substitution_dict, handle)
        if len(self.target)>1:
            for additional_target in self.target[1:]:
                dest_path = os.path.join(additional_target, f'reference/novel_isoform_del_{str(self.novel_read_n)}_{str(self.novel_read_pct)}.pkl')
                shutil.copyfile(self.novel_isoform_del_path, dest_path)
                dest_path = os.path.join(additional_target, 'reference/novel_name_substitutions.pkl')
                shutil.copyfile(self.novel_name_substitution_path, dest_path)

    def _update_saved_matrices_for_mode(self, subset_genes, folder_path_list, splicing=None):
        generated_paths = []
        gene_dfs, transcript_dfs = [], []
        for folder_path in folder_path_list:
            gene_df, transcript_df, out_paths = self._load_subset_pickles_df(folder_path, subset_genes)
            gene_dfs.append(gene_df)
            transcript_dfs.append(transcript_df)
            generated_paths.extend(out_paths)
        for i, folder_path in enumerate(folder_path_list):
            existing_gene_df = self._load_saved_matrix_df(folder_path, 'gene', splicing=splicing)
            existing_transcript_df = self._load_saved_matrix_df(folder_path, 'transcript', splicing=splicing)
            merged_gene_df = self._splice_matrix_df(existing_gene_df, gene_dfs[i], subset_genes, level='gene')
            merged_transcript_df = self._splice_matrix_df(existing_transcript_df, transcript_dfs[i], subset_genes, level='transcript')
            self._save_matrix_df(merged_gene_df, folder_path, 'gene', splicing=splicing)
            self._save_matrix_df(merged_transcript_df, folder_path, 'transcript', splicing=splicing)
        for path in generated_paths:
            os.remove(path)


    def generate_count_matrix_by_gene(self, gene, read_selection_pkl, splicing = None):
        # CompatibleMatrixPaths = '/scr1/users/xu3/singlecell/project_singlecell/sample7_8_ont/sample7/compatible_matrix'
        # read_selection_pkl_paths = '/scr1/users/xu3/singlecell/project_singlecell/sample7_8_ont/sample7/auxillary/read_selection.pkl'
        # read_selection_pkl: keys must add sample index
        if splicing=='spliced':
            count_matrix_folder_path_list = self.count_matrix_spliced_folder_path_list
            compatible_matrix_folder_path_list = self.spliced_compatible_matrix_folder_path_list
        elif splicing=='unspliced':
            count_matrix_folder_path_list = self.count_matrix_unspliced_folder_path_list
            compatible_matrix_folder_path_list = self.unspliced_compatible_matrix_folder_path_list
        else:
            count_matrix_folder_path_list = self.count_matrix_folder_path_list
            compatible_matrix_folder_path_list = self.compatible_matrix_folder_path_list
        for folder in count_matrix_folder_path_list:
            os.makedirs(folder, exist_ok=True)
        pattern = re.compile(r'_ENS.+\.csv')
        files_with_indicators = [
            (os.path.join(CompatibleMatrixPath, i), idx)  # Tuple with file path and index
            for idx, CompatibleMatrixPath in enumerate(compatible_matrix_folder_path_list)
            for i in os.listdir(CompatibleMatrixPath)
            if '.csv' in i and pattern.sub('', i) == gene]
        df_list = []
        for f, idx in files_with_indicators:
            df = pd.read_csv(f)
            df.columns = ['Cell'] + df.columns.tolist()[1:]
            if self.parse:
                df.Cell = df['Cell'].str.rsplit('_', n=1).str[0].tolist()
            df['Cell'] = df['Cell'] + f':sample{idx}'
            df = df[df['Cell'].isin([cell for cell in df['Cell'] if read_selection_pkl.get(cell) == 1])]  # filtering out reads not keep
            if df.shape[0] > 0:
                df['Cell'] = df['Cell'].str.rsplit('_', n=1).str[0].tolist()
                df['Cell'] = df['Cell'] + f':sample{idx}'
                df = df.set_index('Cell')
                df_list.append(df)
        if len(df_list) == 0:
            return {gene: []}, {gene: []}
        result_df = pd.concat(df_list, axis=0)
        df = result_df.fillna(0).astype(int)
        #novel isoform in annotation but not compatible matrix
        novel_isoform_del = [iso_name for iso_name in list(self.annotation_pkl[gene][2].keys()) if iso_name.startswith('novelIsoform_') and iso_name not in df.columns.tolist()]
        novel_name_substitution = []
        if self.group_novel:
            df, novel_isoform_name_mapping = pp.group_novel_isoform(df, geneStrand=self.annotation_pkl[gene][0]['geneStrand'],
                                                                    parse=self.parse)
        else:
            novel_isoform_name_mapping = None
        if novel_isoform_name_mapping is not None:
            #delete isoform being subsituted
            novel_isoform_del += [delete for delete, keep in novel_isoform_name_mapping.items() if delete != keep] #{isoform0 (delete): isoform1}
            #record name subsitution
            novel_name_substitution += [(delete, keep) for delete, keep in novel_isoform_name_mapping.items() if delete != keep]
        # split df into samples
        df['sample_id'] = df.index.str.split(':').str[1]
        df.index = df.index.str.split(':').str[0]
        df_list = [group.drop(columns='sample_id') for _, group in df.groupby('sample_id')]
        # deal each sample separately
        for i, df in enumerate(df_list):
            # --------delete isoforms without reads
            df_isoform = df.sum(axis=0) > 0  # cols have read
            novel_isoform_del += [isoname for isoname in df_isoform[df_isoform==False].index.tolist() if isoname.startswith('novelIsoform_')]
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
                df, dup_iso_name_novel = deduplicate_col(df)  # delete same mapping isoforms
                novel_isoform_del += dup_iso_name_novel
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
                    novel_read_n = 0 if splicing is not None else self.novel_read_n
                    novel_read_pct = 0 if splicing is not None else self.novel_read_pct
                    #drop based on absolute read support
                    novel_isoform_drop0 = df_novel.sum(axis=0) < novel_read_n
                    novel_isoform_drop0 = novel_isoform_drop0[novel_isoform_drop0].index.tolist()
                    # drop based on read pct
                    novel_isoform_drop1 = df.sum(axis=0)/df.sum(axis=0).sum() < novel_read_pct
                    novel_isoform_drop1 = [n_iso for n_iso in novel_isoform_drop1[novel_isoform_drop1].index.tolist() if n_iso.startswith('novel')]
                    novel_isoform_drop = list(set(novel_isoform_drop0 + novel_isoform_drop1))
                    novel_isoform_del += novel_isoform_drop
                    if novel_isoform_drop:  # only proceed if there are columns to drop
                        df_drop = df.loc[:, novel_isoform_drop].sum(axis=1).tolist()
                        df = df.drop(columns=novel_isoform_drop)
                        df['uncategorized_novel'] = df_drop
            if df.shape[1] == 0:
                return {gene: novel_isoform_del}, {gene: novel_name_substitution}
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
                with open(os.path.join(count_matrix_folder_path_list[i], str(gene) + '_unfiltered_count.pickle'), 'wb') as f:
                    pickle.dump((triple_gene, triple_transcript), f)
        return {gene: novel_isoform_del}, {gene: novel_name_substitution}

    def generate_count_matrix_by_gene_list(self, gene_list, read_selection_pkl, splicing = None):
        novel_isoform_del_dict , novel_name_substitution_dict = {}, {}
        for gene in gene_list:
            novel_isoform_del_dict_gene, novel_name_substitution_dict_gene = self.generate_count_matrix_by_gene(gene, read_selection_pkl, splicing = splicing)
            novel_isoform_del_dict.update(novel_isoform_del_dict_gene)
            novel_name_substitution_dict.update(novel_name_substitution_dict_gene)
        return novel_isoform_del_dict, novel_name_substitution_dict
    def read_filter(self):
        read_selection_pkl = {}
        for i, path in enumerate(self.read_selection_pkl_path_list):
            read_selection_pkl_ = pp.load_pickle(path)
            read_selection_pkl_updated = {key + f':sample{i}': value for key, value in read_selection_pkl_.items()}
            read_selection_pkl.update(read_selection_pkl_updated)
        return read_selection_pkl

    def update_multiple_samples_incremental(self, generate_splicing=False):
        if self.gene_subset is None or len(self.gene_subset) == 0:
            raise ValueError('Incremental count-matrix update requires gene_subset.')
        subset_genes = sorted(set(self.gene_subset))
        self._load_annotation_pkl()
        read_selection_pkl = self.read_filter()
        self.logger.info(f'updating count matrices for {len(subset_genes)} genes')
        genes_list = [gene_list for gene_list in split_list(subset_genes, self.workers) if len(gene_list) > 0]
        results = Parallel(n_jobs=self.workers)(
            delayed(self.generate_count_matrix_by_gene_list)(gene_list, read_selection_pkl) for gene_list in genes_list
        )
        subset_novel_isoform_del, subset_novel_name_substitution = {}, {}
        for d1, d2 in results:
            subset_novel_isoform_del.update(d1)
            subset_novel_name_substitution.update(d2)
        self._merge_novel_metadata(subset_novel_isoform_del, subset_novel_name_substitution)
        self._update_saved_matrices_for_mode(subset_genes, self.count_matrix_folder_path_list, splicing=None)
        if generate_splicing:
            Parallel(n_jobs=self.workers)(
                delayed(self.generate_count_matrix_by_gene_list)(gene_list, read_selection_pkl, splicing='spliced')
                for gene_list in genes_list
            )
            self._update_saved_matrices_for_mode(subset_genes, self.count_matrix_spliced_folder_path_list, splicing='spliced')
            Parallel(n_jobs=self.workers)(
                delayed(self.generate_count_matrix_by_gene_list)(gene_list, read_selection_pkl, splicing='unspliced')
                for gene_list in genes_list
            )
            self._update_saved_matrices_for_mode(subset_genes, self.count_matrix_unspliced_folder_path_list, splicing='unspliced')

    def generate_multiple_samples(self, generate_splicing = False):
        pattern = re.compile(r'_ENS.+\.csv')
        Genes = []
        for p in self.compatible_matrix_folder_path_list:
            Genes_ = [g for g in os.listdir(p) if 'csv' in g]
            Genes_ = [pattern.sub('', g) for g in Genes_]
            Genes.append(Genes_)
        Genes = flatten_list(Genes)
        Genes = list(set(Genes))
        if self.gene_subset is not None:
            gene_names_set = set(self.gene_subset)
            Genes = [g for g in Genes if g in gene_names_set]
            self.logger.info(f'Gene subset applied: {len(Genes)} genes match the {len(gene_names_set)} requested gene names')
        for count_path in self.count_matrix_folder_path_list:
            os.makedirs(count_path, exist_ok=True)
        if generate_splicing:
            for count_path in self.count_matrix_spliced_folder_path_list:
                os.makedirs(count_path, exist_ok=True)
            for count_path in self.count_matrix_unspliced_folder_path_list:
                os.makedirs(count_path, exist_ok=True)
        adata_gene_unfiltered_list, adata_transcript_unfiltered_list = [], []
        adata_gene_unfiltered_list_spliced, adata_transcript_unfiltered_list_spliced = [], []
        adata_gene_unfiltered_list_unspliced, adata_transcript_unfiltered_list_unspliced = [], []
        self._load_annotation_pkl()
        self.logger.info(f'generating read filter')
        read_selection_pkl = self.read_filter()
        self.logger.info(f'generating count matrix pickles at: {self.count_matrix_folder_path_list}')
        Genes_list = split_list(Genes, self.workers)
        # ---- generate overall count matrix
        results = Parallel(n_jobs=self.workers)(delayed(self.generate_count_matrix_by_gene_list)(gene_list, read_selection_pkl) for gene_list in Genes_list)
        novel_isoform_del, novel_name_substitution = {}, {}
        for d1, d2 in results:
            novel_isoform_del.update(d1)
            novel_name_substitution.update(d2)
        self.logger.info('count matrix generated')
        #save novel_isoform_del
        self.novel_isoform_del_dict = novel_isoform_del
        self.novel_name_substitution_dict = novel_name_substitution
        # ---- generate spliced/unspliced count matrix
        if generate_splicing:
            self.logger.info('generate count matrix for spliced')
            Parallel(n_jobs=self.workers)(
                delayed(self.generate_count_matrix_by_gene_list)(gene_list, read_selection_pkl, splicing = 'spliced') for gene_list in Genes_list)
            self.logger.info('generate count matrix for unspliced')
            Parallel(n_jobs=self.workers)(
                delayed(self.generate_count_matrix_by_gene_list)(gene_list, read_selection_pkl, splicing = 'unspliced') for gene_list in Genes_list)
        with open(self.novel_isoform_del_path, 'wb') as f:
            pickle.dump(self.novel_isoform_del_dict, f)
        with open(self.novel_name_substitution_path, 'wb') as f:
            pickle.dump(self.novel_name_substitution_dict, f)
        if len(self.target)>1:
            for additional_target in self.target[1:]:
                dest_path = os.path.join(additional_target, f'reference/novel_isoform_del_{str(self.novel_read_n)}_{str(self.novel_read_pct)}.pkl')
                shutil.copyfile(self.novel_isoform_del_path, dest_path)
                dest_path = os.path.join(additional_target, 'reference/novel_name_substitutions.pkl')
                shutil.copyfile(self.novel_name_substitution_path, dest_path)
        #--------overall count mat--------#
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
        # --------splicing count mat--------#
        if generate_splicing:
            #spliced
            for i, count_path in enumerate(self.count_matrix_spliced_folder_path_list):
                out_paths_unfiltered = [os.path.join(count_path, f) for f in os.listdir(count_path) if
                                        f.endswith('_unfiltered_count.pickle')]
                self.logger.info(f'reading {len(out_paths_unfiltered)} count pickles')
                Out_unfiltered = []
                for op in out_paths_unfiltered:
                    Out_unfiltered.append(pp.load_pickle(op))
                self.logger.info('generating spliced count matrix')
                triple_gene_list, triple_transcript_list = list(zip(*Out_unfiltered))
                triple_gene_list = pp.unpack_list(list(triple_gene_list))
                triple_transcript_list = pp.unpack_list(list(triple_transcript_list))
                adata_gene_unfiltered = generate_adata(triple_gene_list)
                adata_transcript_unfiltered = generate_adata(triple_transcript_list)
                self.logger.info('removing pickles')
                for op in out_paths_unfiltered:
                    os.remove(op)
                adata_gene_unfiltered_list_spliced.append(adata_gene_unfiltered)
                adata_transcript_unfiltered_list_spliced.append(adata_transcript_unfiltered)
            self.adata_gene_unfiltered_list_spliced = adata_gene_unfiltered_list_spliced
            self.adata_transcript_unfiltered_list_spliced = adata_transcript_unfiltered_list_spliced
            # unspliced
            for i, count_path in enumerate(self.count_matrix_unspliced_folder_path_list):
                out_paths_unfiltered = [os.path.join(count_path, f) for f in os.listdir(count_path) if
                                        f.endswith('_unfiltered_count.pickle')]
                self.logger.info(f'reading {len(out_paths_unfiltered)} count pickles')
                Out_unfiltered = []
                for op in out_paths_unfiltered:
                    Out_unfiltered.append(pp.load_pickle(op))
                self.logger.info('generating unspliced count matrix')
                triple_gene_list, triple_transcript_list = list(zip(*Out_unfiltered))
                triple_gene_list = pp.unpack_list(list(triple_gene_list))
                triple_transcript_list = pp.unpack_list(list(triple_transcript_list))
                adata_gene_unfiltered = generate_adata(triple_gene_list)
                adata_transcript_unfiltered = generate_adata(triple_transcript_list)
                self.logger.info('removing pickles')
                for op in out_paths_unfiltered:
                    os.remove(op)
                adata_gene_unfiltered_list_unspliced.append(adata_gene_unfiltered)
                adata_transcript_unfiltered_list_unspliced.append(adata_transcript_unfiltered)
            self.adata_gene_unfiltered_list_unspliced = adata_gene_unfiltered_list_unspliced
            self.adata_transcript_unfiltered_list_unspliced = adata_transcript_unfiltered_list_unspliced
    def save_multiple_samples(self, generate_splicing = False):
        if self.mtx:
            for i in range(self.n_samples):
                # save gene
                self.logger.info('saving count matrix in mtx format')
                print('saving count matrix on gene level ')
                gene_meta_unfiltered = {'obs': self.adata_gene_unfiltered_list[i].obs.index.tolist(),
                                        "var": self.adata_gene_unfiltered_list[i].var.index.tolist()}
                fn_pickle = os.path.join(self.count_matrix_folder_path_list[i],
                             f'adata_gene_{self.novel_read_n}_{self.novel_read_pct}.pickle')
                fn_mtx = os.path.join(self.count_matrix_folder_path_list[i],
                                         f'adata_gene_{self.novel_read_n}_{self.novel_read_pct}.mtx')
                with open(fn_pickle,'wb') as f:
                    pickle.dump(gene_meta_unfiltered, f)
                mmwrite(fn_mtx,self.adata_gene_unfiltered_list[i].X)
                # save transcript
                print('saving count matrix on transcript level ')
                transcript_meta_unfiltered = {'obs': self.adata_transcript_unfiltered_list[i].obs.index.tolist(),
                                              "var": self.adata_transcript_unfiltered_list[i].var.index.tolist()}
                fn_pickle = os.path.join(self.count_matrix_folder_path_list[i],
                                         f'adata_transcript_{self.novel_read_n}_{self.novel_read_pct}.pickle')
                fn_mtx = os.path.join(self.count_matrix_folder_path_list[i],
                                      f'adata_transcript_{self.novel_read_n}_{self.novel_read_pct}.mtx')
                mmwrite(fn_mtx, self.adata_transcript_unfiltered_list[i].X)
                with open(fn_pickle,'wb') as f:
                    pickle.dump(transcript_meta_unfiltered, f)
                if generate_splicing:
                    #---------spliced
                    self.logger.info('saving spliced count matrix in mtx format')
                    print('saving spliced count matrix on gene level ')
                    gene_meta_unfiltered = {'obs': self.adata_gene_unfiltered_list_spliced[i].obs.index.tolist(),
                                            "var": self.adata_gene_unfiltered_list_spliced[i].var.index.tolist()}
                    with open(os.path.join(self.count_matrix_spliced_folder_path_list[i],
                                           'adata_gene_unfiltered' + str(self.novel_read_n) + '.pickle'), 'wb') as f:
                        pickle.dump(gene_meta_unfiltered, f)
                    mmwrite(os.path.join(self.count_matrix_spliced_folder_path_list[i],
                                         'adata_gene_unfiltered' + str(self.novel_read_n) + '.mtx'),
                            self.adata_gene_unfiltered_list_spliced[i].X)
                    # save transcript
                    print('saving spliced count matrix on transcript level ')
                    transcript_meta_unfiltered = {'obs': self.adata_transcript_unfiltered_list_spliced[i].obs.index.tolist(),
                                                  "var": self.adata_transcript_unfiltered_list_spliced[i].var.index.tolist()}
                    mmwrite(os.path.join(self.count_matrix_spliced_folder_path_list[i],
                                         'adata_transcript_unfiltered' + str(self.novel_read_n) + '.mtx'),
                            self.adata_transcript_unfiltered_list_spliced[i].X)
                    with open(os.path.join(self.count_matrix_spliced_folder_path_list[i],
                                           'adata_transcript_unfiltered' + str(self.novel_read_n) + '.pickle'), 'wb') as f:
                        pickle.dump(transcript_meta_unfiltered, f)
                    # ---------unspliced
                    self.logger.info('saving unspliced count matrix in mtx format')
                    print('saving unspliced count matrix on gene level ')
                    gene_meta_unfiltered = {'obs': self.adata_gene_unfiltered_list_unspliced[i].obs.index.tolist(),
                                            "var": self.adata_gene_unfiltered_list_unspliced[i].var.index.tolist()}
                    with open(os.path.join(self.count_matrix_unspliced_folder_path_list[i],
                                           'adata_gene_unfiltered' + str(self.novel_read_n) + '.pickle'), 'wb') as f:
                        pickle.dump(gene_meta_unfiltered, f)
                    mmwrite(os.path.join(self.count_matrix_unspliced_folder_path_list[i],
                                         'adata_gene_unfiltered' + str(self.novel_read_n) + '.mtx'),
                            self.adata_gene_unfiltered_list_unspliced[i].X)
                    # save transcript
                    print('saving unspliced count matrix on transcript level ')
                    transcript_meta_unfiltered = {
                        'obs': self.adata_transcript_unfiltered_list_unspliced[i].obs.index.tolist(),
                        "var": self.adata_transcript_unfiltered_list_unspliced[i].var.index.tolist()}
                    mmwrite(os.path.join(self.count_matrix_unspliced_folder_path_list[i],
                                         'adata_transcript_unfiltered' + str(self.novel_read_n) + '.mtx'),
                            self.adata_transcript_unfiltered_list_unspliced[i].X)
                    with open(os.path.join(self.count_matrix_unspliced_folder_path_list[i],
                                           'adata_transcript_unfiltered' + str(self.novel_read_n) + '.pickle'),
                              'wb') as f:
                        pickle.dump(transcript_meta_unfiltered, f)


        if self.csv:
            self.logger.info('saving count matrix in csv format')
            for i in range(self.n_samples):
                output_gene_unfiltered = os.path.join(self.count_matrix_folder_path_list[i],
                                         f'adata_gene_{self.novel_read_n}_{self.novel_read_pct}.csv')
                output_transcript_unfiltered = os.path.join(self.count_matrix_folder_path_list[i],
                                         f'adata_transcript_{self.novel_read_n}_{self.novel_read_pct}.csv')
                # save gene
                print('saving count matrix on gene level ')
                adata_gene_unfiltered_df = self.adata_gene_unfiltered_list[i].to_df()
                adata_gene_unfiltered_df.to_csv(output_gene_unfiltered)
                # save transcript
                print('saving count matrix on transcript level ')
                adata_transcript_unfiltered_df = self.adata_transcript_unfiltered_list[i].to_df()
                adata_transcript_unfiltered_df.to_csv(output_transcript_unfiltered)

                if generate_splicing:
                    #---------spliced
                    self.logger.info('saving spliced count matrix in csv format')
                    output_gene_unfiltered = os.path.join(self.count_matrix_spliced_folder_path_list[i],
                                                          'adata_gene_unfiltered' + str(self.novel_read_n) + '.csv')
                    output_transcript_unfiltered = os.path.join(self.count_matrix_spliced_folder_path_list[i],
                                                                'adata_transcript_unfiltered' + str(
                                                                    self.novel_read_n) + '.csv')
                    # save gene
                    print('saving count matrix on gene level ')
                    adata_gene_unfiltered_df = self.adata_gene_unfiltered_list_spliced[i].to_df()
                    adata_gene_unfiltered_df.to_csv(output_gene_unfiltered)
                    # save transcript
                    print('saving count matrix on transcript level ')
                    adata_transcript_unfiltered_df = self.adata_transcript_unfiltered_list_spliced[i].to_df()
                    adata_transcript_unfiltered_df.to_csv(output_transcript_unfiltered)

                    # ---------unspliced
                    self.logger.info('saving unspliced count matrix in csv format')
                    output_gene_unfiltered = os.path.join(self.count_matrix_unspliced_folder_path_list[i],
                                                          'adata_gene_unfiltered' + str(self.novel_read_n) + '.csv')
                    output_transcript_unfiltered = os.path.join(self.count_matrix_unspliced_folder_path_list[i],
                                                                'adata_transcript_unfiltered' + str(
                                                                    self.novel_read_n) + '.csv')
                    # save gene
                    print('saving count matrix on gene level ')
                    adata_gene_unfiltered_df = self.adata_gene_unfiltered_list_unspliced[i].to_df()
                    adata_gene_unfiltered_df.to_csv(output_gene_unfiltered)
                    # save transcript
                    print('saving count matrix on transcript level ')
                    adata_transcript_unfiltered_df = self.adata_transcript_unfiltered_list_unspliced[i].to_df()
                    adata_transcript_unfiltered_df.to_csv(output_transcript_unfiltered)

    def _extract_attribute(self, attributes, attribute_name):
        try:
            return attributes.split(f'{attribute_name} "')[1].split('"')[0]
        except IndexError:
            return None
    def filter_gtf(self):
        self.logger.info('updating gtf annotation file')
        for target in self.target:
            input_gtf = os.path.join(target, f'reference/SCOTCH_updated_annotation.gtf')
            output_gtf = os.path.join(target, f'reference/SCOTCH_updated_annotation_filtered_{self.novel_read_n}_{self.novel_read_pct}.gtf')
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

    def finalize_read_isoform_mapping(self):
        """Filter and rename novel isoforms in the read-isoform mapping TSV.

        Streams the large all_read_isoform_exon_mapping.tsv line by line to:
        1. Remove reads mapped to novel isoforms that were filtered out (below threshold)
        2. Rename novel isoforms using the grouped name convention
        3. Update exon index and coordinates for renamed isoforms
        """
        import csv

        # Build per-gene rename and drop lookups
        rename_by_gene = {}
        for gene, substitutions in self.novel_name_substitution_dict.items():
            if substitutions:
                rename_by_gene[gene] = dict(substitutions)

        drop_by_gene = {}
        for gene, del_list in self.novel_isoform_del_dict.items():
            if del_list:
                drop_by_gene[gene] = set(del_list)

        for target in self.target:
            auxillary_dir = os.path.join(target, 'auxillary')
            tsv_path = os.path.join(auxillary_dir, 'all_read_isoform_exon_mapping.tsv')
            if not os.path.exists(tsv_path):
                self.logger.warning(f'Read-isoform mapping TSV not found: {tsv_path}')
                continue

            output_path = os.path.join(auxillary_dir,
                f'all_read_isoform_exon_mapping_filtered_{self.novel_read_n}_{self.novel_read_pct}.tsv')
            tmp_path = output_path + '.tmp'

            self.logger.info(f'Filtering read-isoform mapping: {tsv_path}')
            n_kept = 0
            n_dropped = 0
            n_renamed = 0

            with open(tsv_path, 'r', newline='') as fin, open(tmp_path, 'w', newline='') as fout:
                reader = csv.DictReader(fin, delimiter='\t')
                writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter='\t')
                writer.writeheader()

                for row in reader:
                    iso = row['Isoform']
                    gene_name_raw = row['geneName']
                    gene = re.sub(r'[\/\\:\*\?"<>|]', '.', gene_name_raw)

                    # Non-novel isoforms: check if in drop set (shouldn't be, but safe)
                    if not iso.startswith('novelIsoform_'):
                        writer.writerow(row)
                        n_kept += 1
                        continue

                    # Rename using substitution dict
                    rename_map = rename_by_gene.get(gene, {})
                    final_iso = rename_map.get(iso, iso)

                    # Check if the final isoform should be dropped
                    if final_iso in drop_by_gene.get(gene, set()):
                        n_dropped += 1
                        continue

                    # Also drop if original name is in drop set (before renaming)
                    if iso in drop_by_gene.get(gene, set()):
                        n_dropped += 1
                        continue

                    # Update row if renamed
                    if final_iso != iso:
                        n_renamed += 1
                        row['Isoform'] = final_iso
                        # Update exon info from annotation if available
                        if self.annotation_pkl and gene in self.annotation_pkl:
                            gene_info = self.annotation_pkl[gene]
                            exon_info = gene_info[1]
                            isoform_info = gene_info[2]
                            if final_iso in isoform_info:
                                exon_idx = isoform_info[final_iso]
                                exon_coords = pp.merge_exons(
                                    [(exon_info[ind][0], exon_info[ind][1]) for ind in exon_idx])
                                row['Exon Index'] = ','.join(map(str, exon_idx))
                                row['Exon Coordinates'] = ','.join(
                                    f"{exon}" for exon in exon_coords)

                    writer.writerow(row)
                    n_kept += 1

            os.replace(tmp_path, output_path)
            self.logger.info(
                f'Read-isoform mapping filtered: {n_kept} kept, {n_dropped} dropped, '
                f'{n_renamed} renamed -> {output_path}')
