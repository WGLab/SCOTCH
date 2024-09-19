import pysam
import pandas as pd
import pickle
from itertools import chain
import numpy as np
import os
from Bio.Seq import Seq
import networkx as nx
import community.community_louvain as community_louvain
from collections import defaultdict

#bam="/scr1/users/xu3/singlecell/project_singlecell/sample8_R10/bam/sample8_R10.filtered.bam"
#gene_pkl="/scr1/users/xu3/singlecell/project_singlecell/M4/reference/geneStructureInformation.pkl"



#---------------some utility functions----------------------#
# Function to convert the dictionary to GTF


def sort_multigeneInfo(Info_multigenes):
    Info_multigenes_sort = []
    for info_singlegene in Info_multigenes:
        isoformInfo = info_singlegene[2]
        if len(isoformInfo)>0:
            isoformInfo = dict(sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))
            info_singlegene[0]['isoformNames'] = list(isoformInfo.keys())
            info_singlegene[2] = isoformInfo
        Info_multigenes_sort.append(info_singlegene)
    return Info_multigenes_sort

def summarise_metagene(Info_multigenes):
    n_genes = len(Info_multigenes)
    geneChr = Info_multigenes[0][0]['geneChr']
    start, end, isoforms = [], [],[]
    for i in range(n_genes):
        start.append(Info_multigenes[i][0]['geneStart'])
        end.append(Info_multigenes[i][0]['geneEnd'])
    start = min(start)
    end = max(end)
    return geneChr, start, end

def load_pickle(file):
    if os.path.exists(file):
        with open(file,'rb') as file:
            data=pickle.load(file)
    else:
        data = None
    return data

def unpack_list(x):
    return list(chain(*x))

def iter_length(iter):
    return sum(1 for _ in iter)

def check_values_in_list(list_a, list_b, rule='all'):
    out=None
    if rule=='all':
        out = all(value in list_b for value in list_a)
    if rule=='any':
        out = any(value in list_b for value in list_a)
    return out

def remove_dict_items(dictionary,threshold):
    keys_to_remove = []
    exons = dictionary.keys()
    for exon in exons:
        cutoff = threshold[exon]
        if dictionary[exon] < cutoff:
            keys_to_remove.append(exon)
    for key in keys_to_remove:
        del dictionary[key]
    return dictionary

def isoformInfo_to_df(isoformInfo, geneInfo):
    all_exons = set(range(geneInfo['numofExons']))
    df = pd.DataFrame(-1, index=isoformInfo.keys(), columns=sorted(all_exons))
    for isoform, exons in isoformInfo.items():
        df.loc[isoform, exons] = 1
    return df

def isoformInfo_to_onehot(isoformInfo, geneInfo):
    all_exons = set(range(geneInfo['numofExons']))
    assignment_vectors = []
    for isoform, exons in isoformInfo.items():
        vector = [1 if exon in exons else -1 for exon in all_exons]
        assignment_vectors.append(vector)
    return assignment_vectors





#---------------some functions for novel reads----------------------#
def read_similarity(df_assign):
    data = df_assign.to_numpy()
    ones = (data == 1).astype(int)
    minus_ones = (data == -1).astype(int)
    conflict_matrix = np.tensordot(ones, minus_ones, axes=([1], [1]))
    similarity_matrix = (conflict_matrix == 0).astype(int)
    consistent_matrix = np.tensordot(ones, ones, axes=([1], [1]))+np.tensordot(minus_ones, minus_ones, axes=([1], [1]))
    similarity_matrix = similarity_matrix*consistent_matrix
    np.fill_diagonal(similarity_matrix, 0)
    return similarity_matrix


def novelid_to_exonid(novelid):
    exonid = []
    position = 0
    while novelid > 0:
        if novelid & 1:
            exonid.append(position)
        novelid >>= 1
        position += 1
    return exonid


def find_novel(df_assign, df_pct):
    similarity_matrix_assign = read_similarity(df_assign)
    G = nx.Graph()
    for i in range(len(df_pct)):
        G.add_node(i)
    for i in range(len(similarity_matrix_assign)):
        for j in range(i + 1, len(similarity_matrix_assign)):
            if similarity_matrix_assign[i][j] > 0:
                G.add_edge(i, j, weight=similarity_matrix_assign[i][j])
    partition = community_louvain.best_partition(G)
    read_groups = defaultdict(list)
    for read, group in partition.items():
        read_groups[group].append(read)
    novelisoform_dict, assigns = {}, []  ##this is what i want
    for k in range(len(read_groups)):
        if len(read_groups[k]) > 1:
            assign = df_assign.iloc[read_groups[k]].mean(axis=0).round().tolist()
            pct = df_pct.iloc[read_groups[k]].mean(axis=0).round().tolist()
            for i, j in enumerate(pct):
                if assign[i] > -1 and j > 0:
                    assign[i] = 1
                else:
                    assign[i] = -1
            isoform_index = [i for i, p in enumerate(pct) if assign[i] > -1 and p > 0]
            isoform_id = 'novelIsoform_' + str(sum([2 ** e for e in isoform_index]))
            if len(assign) > -sum(assign) and isoform_id not in novelisoform_dict.keys():
                assigns.append(assign)
                novelisoform_dict[isoform_id] = (isoform_index)
    return novelisoform_dict, assigns

def find_novel_by_chunk(df_pct, df_assign, chunk_size=1500):
    row_count = df_pct.shape[0]
    num_chunks = max(1, np.ceil(row_count / chunk_size).astype(int))
    df_pct_list = np.array_split(df_pct, num_chunks)
    df_assign_list = np.array_split(df_assign, num_chunks)
    novelisoform_dict_list, novel_df_list, novel_df_empty = [], [], []
    i = 0
    print(str(len(df_assign_list)) + ' chunks in total')
    while i < len(df_assign_list):
        print('chunk ' + str(i) + ' out of ' + str(len(df_assign_list)))
        novelisoform_dict, assigns = find_novel(df_assign_list[i], df_pct_list[i])
        novelisoform_dict, novel_df, novel_df_empty = map_read_to_isoform(novelisoform_dict, assigns, df_assign)
        if novel_df is not None:
            if novel_df.shape[0] > 10:
                novelisoform_dict_list.append(novelisoform_dict)
                novel_df_list.append(novel_df)
                # update split
                df_assign = df_assign.loc[novel_df_empty]
                df_pct = df_pct.loc[novel_df_empty]
                row_count = df_pct.shape[0]
                num_chunks = max(1, np.ceil(row_count / chunk_size).astype(int))
                df_pct_list = np.array_split(df_pct, num_chunks)
                df_assign_list = np.array_split(df_assign, num_chunks)
            else:
                i += 1
        else:
            i += 1
    novelisoform_dict = {}
    for single_dict in novelisoform_dict_list:
        novelisoform_dict.update(single_dict)
    return novelisoform_dict, novel_df_list, novel_df_empty

def group_novel_isoform(df, geneStrand):
    df_novel = df.filter(like='novelIsoform')
    novel_isoform_name_mapping = {}
    if df_novel.shape[1]>1:
        df_uncategorized = df.filter(like='uncategorized')
        drop_c = df_novel.columns.tolist()+df_uncategorized.columns.tolist()
        df_existing = df.drop(columns=drop_c)
        #sort novel isoform id
        novel_isoform_id = [int(ni.split('_')[1]) for ni in df_novel.columns.tolist()]
        novel_isoform_id_len = [len(novelid_to_exonid(i)) for i in novel_isoform_id]
        sorted_pairs = sorted(zip(novel_isoform_id_len, novel_isoform_id))
        novel_isoform_id = [id for length, id in sorted_pairs]
        #group novel isoforms
        novel_isoform_group_list = []
        parent_id = None
        while len(novel_isoform_id) > 0:
            query_id = None if parent_id is None else novel_isoform_group_list[-1][-1]
            novel_isoform_group_list, novel_isoform_id, parent_id = pair_isoform(novel_isoform_group_list, novel_isoform_id,
                                                                                 query_id, geneStrand)
        novel_isoform_group_name = ['novelIsoformGroup_'+str(i[-1]) for i in novel_isoform_group_list]
        df_novel_group = []
        for i in range(len(novel_isoform_group_list)):
            isoform_group_name = ['novelIsoform_' + str(ii) for ii in novel_isoform_group_list[i]]
            for isoform_name in isoform_group_name:
                novel_isoform_name_mapping[isoform_name] = novel_isoform_group_name[i]
            df_novel_ = df_novel[isoform_group_name].max(axis=1)
            df_novel_group.append(df_novel_)
        df_novel_grouped = pd.concat(df_novel_group, axis=1)
        df_novel_grouped.columns = novel_isoform_group_name
        df = pd.concat([df_existing, df_novel_grouped, df_uncategorized], axis=1)
    return df, novel_isoform_name_mapping


def pair_isoform(novel_isoform_group_list, novel_isoform_id, query_id, geneStrand):
    novel_isoform_id_list = novel_isoform_id.copy()
    # create new isoform group
    if query_id is None:
        child_id = novel_isoform_id_list[0]
        novel_isoform_id_list.remove(child_id)
        parent_id = find_parent_isoform(novel_isoform_id_list, child_id, geneStrand)
        if parent_id is not None:
            novel_isoform_id_list.remove(parent_id)
            novel_isoform_group_list.append([child_id, parent_id])
        else:
            novel_isoform_group_list.append([child_id])
    # update existing isoform group
    else:
        child_id = query_id
        parent_id = find_parent_isoform(novel_isoform_id_list, child_id, geneStrand)
        if parent_id is not None:
            novel_isoform_id_list.remove(parent_id)
            novel_isoform_group_list[-1].append(parent_id)
    return novel_isoform_group_list, novel_isoform_id_list,parent_id


def find_parent_isoform(novel_isoform_id_list, child_id, geneStrand):
    for i in range(len(novel_isoform_id_list)):
        child_exon = novelid_to_exonid(child_id)
        parent_exon = novelid_to_exonid(novel_isoform_id_list[i])
        if(set(child_exon).issubset(set(parent_exon))):
            diff_elements = set(parent_exon)-set(child_exon)
            if (geneStrand=='+' and len([ele for ele in diff_elements if ele<min(child_exon)])==len(diff_elements)) or (geneStrand == '-' and len([ele for ele in diff_elements if ele > max(child_exon)]) == len(
                    diff_elements)):
                return novel_isoform_id_list[i]


def generate_gene_annotation_novel_isoform(geneID, isoformID, geneStructureInformation_path, gtf_path,output_gtf_path,filter_existing=None):
    geneStructureInformation = load_pickle(geneStructureInformation_path)
    geneName = geneStructureInformation[geneID][0]['geneName']
    novel_dict_list = []
    for isoformID_ in isoformID:
        exon_ind_list = novelid_to_exonid(isoformID_)
        isoform_annotation = [geneStructureInformation[geneID][1][ind] for ind in exon_ind_list]
        isoform_annotation = merge_exons(isoform_annotation)
        isoform_annotation = [(a + 1, b) for (a, b) in isoform_annotation]# change to 1 based
        new_transcript_id = geneName+'-novelIsoformGroup0'+str(isoformID_)
        novel_dict_list.append((isoform_annotation,new_transcript_id))
    add_novel_isoform_to_gtf(gtf_path, geneID, geneName, novel_dict_list, output_gtf_path,
                                 only_this_gene=True, filter_existing=filter_existing)

def get_transcript_id(attribute):
    parts = attribute.split(';')
    transcript_name_part = next((part for part in parts if 'transcript_id' in part), None)
    if transcript_name_part:
        return transcript_name_part.split('"')[1]
    return None

def add_novel_isoform_to_gtf(gtf_path, gene_id, gene_name, novel_dict_list, output_gtf_path,
                             only_this_gene=False, filter_existing=None):
    # Load the GTF file into a DataFrame
    columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf_data = pd.read_csv(gtf_path, sep='\t', comment='#', names=columns)
    # Find the existing annotations for the gene to get the source, score, strand, and frame
    gene_data = gtf_data[gtf_data['attribute'].str.contains(f'gene_id "{gene_id}";')]
    if gene_data.empty:
        raise ValueError(f'Gene {gene_id} not found in the GTF file.')
    if filter_existing is not None:
        gene_data['transcript_id'] = gene_data['attribute'].apply(get_transcript_id)
        gene_data = gene_data[gene_data['transcript_id'].isin(filter_existing)]
        gene_data.drop('transcript_id', axis=1, inplace=True)
    # Use the first entry as a template
    template = gene_data.iloc[0]
    if only_this_gene:
        gtf_data = gene_data
    new_entry = template.copy()
    if len(novel_dict_list)>0:
        for novel_dict in novel_dict_list:
            novel_exon_annotation, new_transcript_id = novel_dict
            for start, end in novel_exon_annotation:
                new_entry['feature'] = 'exon'
                new_entry['start'] = start
                new_entry['end'] = end
                new_entry[
                    'attribute'] = f'gene_id "{gene_id}"; transcript_id "{new_transcript_id}"; gene_name "{gene_name}"; '
                gtf_data = gtf_data.append(new_entry, ignore_index=True)
            new_entry['feature'] = 'transcript'
            new_entry['start'] = novel_exon_annotation[0][0]
            gtf_data = gtf_data.append(new_entry, ignore_index=True)
    # Write the updated DataFrame to a new GTF file
    gtf_data.to_csv(output_gtf_path, sep='\t', index=False, header=False, quoting=3)

#---------------some functions for read-exon operations----------------------#
def merge_exons(exons):
    if len(exons)<2:
        merged_exons = exons
    else:
        merged_exons = []
        current_start, current_end = exons[0]
        for start, end in exons[1:]:
            if start == current_end:
                current_end = max(current_end, end)
            else:
                merged_exons.append((current_start, current_end))
                current_start, current_end = start, end
        merged_exons.append((current_start, current_end))
    return merged_exons


def read_exon_match(read, Info_singlegene):
    readName, readStart, readEnd = read.qname, read.reference_start, read.reference_end
    referencePositions = read.get_reference_positions(full_length=False)
    geneInfo, exonInfo, isoformInfo = Info_singlegene
    gene_name =geneInfo['geneName']
    gene_length = geneInfo['geneEnd'] - geneInfo['geneStart']
    if geneInfo['geneStrand']=="+":
        readDist = readStart - geneInfo["geneStart"]
    else:
        readDist = readEnd - geneInfo["geneEnd"]
    if readStart >= geneInfo['geneStart'] and readEnd <= geneInfo['geneEnd']:
        mapExons = exon_hit(referencePositions, exonInfo)
        n_mapExons=sum(mapExons)
    else:
        n_mapExons = -1
    return gene_name, readDist, gene_length, n_mapExons


def exon_hit(mapPositions, exonInfo):
    mapExons = []
    exon_index = 0
    for i in mapPositions:
        while i >= exonInfo[exon_index][1] and exon_index < len(exonInfo) - 1:
            exon_index += 1
        if int(i) >= int(exonInfo[exon_index][0]) and int(i) < int(exonInfo[exon_index][1]):
            mapExons.append(exon_index)
    return mapExons

#---------------some functions for read-isoform operations-----------------#


##TODO: this function is for 3' kit, should consider 5'kit
def detect_poly(read, window = 20, n = 15):
    query_sequence = read.query_sequence
    umi = read.get_tag('UR')
    res = False
    nT,nA=0,0
    seq = Seq(query_sequence)
    umi = Seq(umi)
    umi_rc = umi.reverse_complement()
    pos1, pos2 = seq.find(umi),seq.find(umi_rc)
    poly=None
    if pos1>0:
        nT = seq[pos1 +len(umi): pos1 +len(umi)+window].count('T')
        poly = 'T'
    if pos2>0:
        nA = seq[pos2-window:pos2].count('A')
        poly = 'A'
    if nT>=n or nA>=n:
        res = True
    return res, poly

def detect_poly_parse(read, window = 20, n = 10):
    def detect_poly_parse_(query_sequence, window, n, AorT):
        if len(query_sequence) == 0:
            return False
        poly_bool = False
        seq = Seq(query_sequence)
        for i in range(len(seq) - window + 1):
            window_seq = seq[i:i + window]
            poly_bool = None
            if AorT == 'T':
                if window_seq.count('T') > n:
                    poly_bool = True
                    break
            else:
                if window_seq.count('A') > n:
                    poly_bool = True
                    break
        return poly_bool
    query_sequence0 = read.query_sequence[:(read.query_alignment_start + 1)] #soft clip 0
    query_sequence1 = read.query_sequence[read.query_alignment_end:] #soft clip 1
    poly_bool = False
    if detect_poly_parse_(query_sequence0,window, n,'T'):
        poly = 'T'
        poly_bool = True
    elif detect_poly_parse_(query_sequence1, window, n, 'A'):
        poly = 'A'
        poly_bool = True
    else:
        poly = None
    return poly_bool, poly

def choose_gene_from_meta(read, Info_multigenes, lowest_match=0.05, pacbio = False):
    geneName, readDist, geneLength, nMapExons,index = [], [], [], [],[]
    for i, info_singlegene in enumerate(Info_multigenes):
        gene_name, read_dist, gene_length, n_mapExons = read_exon_match(read, info_singlegene)
        geneName.append(gene_name)
        geneLength.append(gene_length)
        nMapExons.append(n_mapExons)
        readDist.append(read_dist)
        index.append(i)
    df = pd.DataFrame({'geneName': geneName, 'readDist': readDist, 'geneLength': geneLength, 'nMapExons': nMapExons, 'index': index})
    df_exon = df[df['nMapExons'] > 0].reset_index(drop=True)
    df_intron = df[df['nMapExons'] == 0].reset_index(drop=True)
    if df_exon.shape[0] == 0 and df_intron.shape[0] > 0:
        df_intron = df_intron.sort_values(by=['geneLength'], ascending=[False])
        ind = df_intron['index'][0]
        read_novelisoform_tuple, read_isoform_compatibleVector_tuple = Map_read_to_gene(read,
                                                                                        Info_multigenes[ind],
                                                                                        lowest_match, geneInfo=None,
                                                                                        exonInfo=None, isoformInfo=None, pacbio = pacbio)
    elif (df_exon.shape[0] > 0):  # the read locates within at least one gene region
        df_exon = df_exon.sort_values(by=['nMapExons', 'readDist', 'geneLength'], ascending=[False, True, False])
        if sum(df_exon.nMapExons) >= 2:
            df_exon = df_exon[df_exon['nMapExons'] > 0].reset_index(drop=True)
            # sort by existing/novel
            inds = df_exon['index'].tolist()
            read_novelisoform_tuple_dict, read_isoform_compatibleVector_tuple_dict, readType = {}, {}, []
            for ind in inds:
                read_novelisoform_tuple_, read_isoform_compatibleVector_tuple_ = Map_read_to_gene(read,
                                                                                                  Info_multigenes[ind],
                                                                                                  lowest_match,
                                                                                                  geneInfo=None,
                                                                                                  exonInfo=None,
                                                                                                  isoformInfo=None, pacbio = pacbio)
                read_novelisoform_tuple_dict[ind] = read_novelisoform_tuple_
                read_isoform_compatibleVector_tuple_dict[ind] = read_isoform_compatibleVector_tuple_
                if read_novelisoform_tuple_ is not None:
                    readType.append('1')  # novel
                elif sum(read_isoform_compatibleVector_tuple_[1]) > 0:
                    readType.append('0')  # existing
                else:
                    readType.append('2')  # uncharacterized
            df_exon['readType'] = readType
            df_exon = df_exon.sort_values(by=['readType', 'nMapExons', 'readDist', 'geneLength'],
                                          ascending=[True, False, True, False])
            ind = df_exon['index'][0]
            read_novelisoform_tuple, read_isoform_compatibleVector_tuple = read_novelisoform_tuple_dict[ind], \
            read_isoform_compatibleVector_tuple_dict[ind]
        else:  # map to only one gene exon
            ind = df_exon['index'].tolist()[0]
            read_novelisoform_tuple, read_isoform_compatibleVector_tuple = Map_read_to_gene(read,
                                                                                            Info_multigenes[ind],
                                                                                            lowest_match, geneInfo=None,
                                                                                            exonInfo=None,
                                                                                            isoformInfo=None, pacbio=pacbio)
    else:
        read_novelisoform_tuple, read_isoform_compatibleVector_tuple, ind = None, None, -1
    return ind, read_novelisoform_tuple, read_isoform_compatibleVector_tuple

def choose_gene_from_meta_parse(read, Info_multigenes, lowest_match=0.05, poly=False):
    geneName, readDist, geneLength, nMapExons,index = [], [], [], [],[]
    for i, info_singlegene in enumerate(Info_multigenes):
        gene_name, read_dist, gene_length, n_mapExons = read_exon_match(read, info_singlegene)
        geneName.append(gene_name)
        geneLength.append(gene_length)
        nMapExons.append(n_mapExons)
        readDist.append(read_dist)
        index.append(i)
    df = pd.DataFrame({'geneName': geneName, 'readDist': readDist, 'geneLength': geneLength, 'nMapExons': nMapExons, 'index': index})
    df_exon = df[df['nMapExons'] > 0].reset_index(drop=True)
    df_intron = df[df['nMapExons'] == 0].reset_index(drop=True)
    if df_exon.shape[0] == 0 and df_intron.shape[0] > 0:
        df_intron = df_intron.sort_values(by=['geneLength'], ascending=[False])
        ind = df_intron['index'][0]
        read_novelisoform_tuple, read_isoform_compatibleVector_tuple = Map_read_to_gene_parse(read,
                                                                                        Info_multigenes[ind],
                                                                                        lowest_match, geneInfo=None,
                                                                                        exonInfo=None, isoformInfo=None,
                                                                                              poly=poly)
    elif (df_exon.shape[0] > 0):  # the read locates within at least one gene region
        df_exon = df_exon.sort_values(by=['nMapExons', 'readDist', 'geneLength'], ascending=[False, True, False])
        if sum(df_exon.nMapExons) >= 2:
            df_exon = df_exon[df_exon['nMapExons'] > 0].reset_index(drop=True)
            # sort by existing/novel
            inds = df_exon['index'].tolist()
            read_novelisoform_tuple_dict, read_isoform_compatibleVector_tuple_dict, readType = {}, {}, []
            for ind in inds:
                read_novelisoform_tuple_, read_isoform_compatibleVector_tuple_ = Map_read_to_gene_parse(read,
                                                                                                  Info_multigenes[ind],
                                                                                                  lowest_match,
                                                                                                  geneInfo=None,
                                                                                                  exonInfo=None,
                                                                                                  isoformInfo=None,
                                                                                                  poly=poly)
                read_novelisoform_tuple_dict[ind] = read_novelisoform_tuple_
                read_isoform_compatibleVector_tuple_dict[ind] = read_isoform_compatibleVector_tuple_
                if read_novelisoform_tuple_ is not None:
                    readType.append('1')  # novel
                elif sum(read_isoform_compatibleVector_tuple_[1]) > 0:
                    readType.append('0')  # existing
                else:
                    readType.append('2')  # uncharacterized
            df_exon['readType'] = readType
            df_exon = df_exon.sort_values(by=['readType', 'nMapExons', 'readDist', 'geneLength'],
                                          ascending=[True, False, True, False])
            ind = df_exon['index'][0]
            read_novelisoform_tuple, read_isoform_compatibleVector_tuple = read_novelisoform_tuple_dict[ind], \
            read_isoform_compatibleVector_tuple_dict[ind]
        else:  # map to only one gene exon
            ind = df_exon['index'].tolist()[0]
            read_novelisoform_tuple, read_isoform_compatibleVector_tuple = Map_read_to_gene_parse(read,
                                                                                            Info_multigenes[ind],
                                                                                            lowest_match, geneInfo=None,
                                                                                            exonInfo=None,
                                                                                            isoformInfo=None,
                                                                                                  poly=poly)
    else:
        read_novelisoform_tuple, read_isoform_compatibleVector_tuple, ind = None, None, -1
    return ind, read_novelisoform_tuple, read_isoform_compatibleVector_tuple

def map_read_to_gene_parse(read, Info_singlegene, lowest_match=0.2, geneInfo = None,exonInfo= None, isoformInfo= None, qualifyExon= None,
                     exonMatch1= None, exonMatch2= None,poly = False):
    #optional: geneInfo,exonInfo, isoformInfo, qualifyExon, exonMatch
    readName, readStart, readEnd = read.qname, read.reference_start, read.reference_end
    referencePositions = read.get_reference_positions(full_length=False)
    if geneInfo is None:
        geneInfo, exonInfo, isoformInfo = Info_singlegene
        isoformInfo = dict(sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))  # long to short
        geneInfo['isoformNames'] = list(isoformInfo.keys())
        qualifyExon = [i for i, (a, b) in enumerate(exonInfo) if b - a >= 20]
        exonMatch1 = [max(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]  # exon length threshold high
        exonMatch2 = [min(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]  # exon length threshold low
    read_isoform_compatibleVector = [1] * len(geneInfo['isoformNames'])#initialize read isoform compativle vector
    #read-exon mapping percentage
    mapExons = exon_hit(referencePositions, exonInfo)
    exon_map_pct, exon_map_vector_trunct = [], []
    if len(mapExons) > 0:
        exonhitbases = dict(pd.Series(mapExons).value_counts())
        exon_map_base = [exonhitbases.get(i, 0) for i in range(len(exonInfo))]
        exon_map_pct = [round(base / (exonInfo[i][1] - exonInfo[i][0]), 2) if base else 0 for i, base in enumerate(exon_map_base)]
        #read-exon mapping vector
        exon_map_vector_mapped = [1 if exon_map_p>=exon_map_up else 0 for exon_map_p, exon_map_up in list(zip(exon_map_base,exonMatch1))]
        exon_map_vector_skipped= [-1 if exon_map_p<=exon_map_low else 0 for exon_map_p, exon_map_low in list(zip(exon_map_base,exonMatch2))]
        exon_map_vector_notrunct = [mapped + skipped for mapped, skipped in list(zip(exon_map_vector_mapped, exon_map_vector_skipped))]
        if poly:
            if geneInfo['geneStrand']=="+":
                index_of_one = exon_map_vector_notrunct.index(1) if 1 in exon_map_vector_notrunct else len(exon_map_vector_notrunct)
                exon_map_vector_trunct = [0 if i < index_of_one else exon_map_vector_notrunct[i] for i in range(len(exon_map_vector_notrunct))]
            else: #"-"
                last_index_of_one = len(exon_map_vector_notrunct) - 1 - exon_map_vector_notrunct[::-1].index(1) if 1 in exon_map_vector_notrunct else -1
                exon_map_vector_trunct = [exon_map_vector_notrunct[i] if i <= last_index_of_one else 0 for i in range(len(exon_map_vector_notrunct))]
        else:
            if 1 in exon_map_vector_notrunct:
                index_of_one = exon_map_vector_notrunct.index(1)
                last_index_of_one = len(exon_map_vector_notrunct) - 1 - exon_map_vector_notrunct[::-1].index(1)
                exon_map_vector_trunct = [
                    exon_map_vector_notrunct[i] if index_of_one <= i <= last_index_of_one else 0
                    for i in range(len(exon_map_vector_notrunct))
                ]
            else:
                # If there are no '1's in the vector, set the whole vector to 0
                exon_map_vector_trunct = [0] * len(exon_map_vector_notrunct)
        #do not consider non-qualified exons
        #exon_map_vector_notrunct = [vec if i in qualifyExon else 0 for i,vec in enumerate(exon_map_vector_notrunct)]
        exon_map_vector_trunct = [vec if i in qualifyExon else 0 for i,vec in enumerate(exon_map_vector_trunct)]
        isoform_exon_map_list = isoformInfo_to_onehot(isoformInfo, geneInfo)
        read_exon_map_list = exon_map_vector_trunct
        read_isoform_compatibleVector = map_read_to_isoform_known(isoform_exon_map_list, read_exon_map_list, exon_map_pct,
                                  max(lowest_match,1-lowest_match), min(lowest_match,1-lowest_match), geneInfo['geneStrand'],
                                                                  qualifyExon, poly)
    else:
        # uncharacterized type3: not map to any exons---intron
        read_isoform_compatibleVector = [0] * len(geneInfo['isoformNames'])
    #---novel/uncharacterized+existing compatible vector
    read_novelisoform_tuple, read_isoform_compatibleVector_tuple = None, None
    # novel isoform
    if sum(read_isoform_compatibleVector) == 0 and len(mapExons) > 0:
        read_novelisoform_tuple = (readName, exon_map_pct, exon_map_vector_trunct)
    else:  # existing isoform/uncategorized
        read_isoform_compatibleVector_tuple = (readName, read_isoform_compatibleVector)
    return read_novelisoform_tuple, read_isoform_compatibleVector_tuple


def Map_read_to_gene_parse(read, Info_singlegene, lowest_match=0.2, geneInfo=None, exonInfo=None, isoformInfo=None, poly=False):
    if geneInfo is None:
        geneInfo, exonInfo, isoformInfo = Info_singlegene
        isoformInfo = dict(
            sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))  # long to short
        geneInfo['isoformNames'] = list(isoformInfo.keys())
    qualifyExon0 = [i for i, (a, b) in enumerate(exonInfo) if b - a >= 20]  # only exons longer than 20 will be used as proof
    exonMatch0_1 = [max(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]
    exonMatch0_2 = [min(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]
    qualifyExon1 = list(range(len(exonInfo)))
    exonMatch1_1 = [max(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]
    exonMatch1_2 = [min(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]
    # conservative match
    read_novelisoform_tuple0, read_isoform_compatibleVector_tuple0 = map_read_to_gene_parse(read, Info_singlegene,
                                                                                      lowest_match, geneInfo,
                                                                                      exonInfo, isoformInfo,
                                                                                      qualifyExon0, exonMatch0_1, exonMatch0_2, poly)
    # take account of small exons
    read_novelisoform_tuple1, read_isoform_compatibleVector_tuple1 = map_read_to_gene_parse(read, Info_singlegene,
                                                                                      lowest_match, geneInfo,
                                                                                      exonInfo,
                                                                                      isoformInfo, qualifyExon1,
                                                                                      exonMatch1_1, exonMatch1_2, poly)
    read_novelisoform_tuple = None
    if read_isoform_compatibleVector_tuple0 is None and read_isoform_compatibleVector_tuple1 is not None:
        read_isoform_compatibleVector_tuple = read_isoform_compatibleVector_tuple1
    elif read_isoform_compatibleVector_tuple0 is not None and read_isoform_compatibleVector_tuple1 is None:
        read_isoform_compatibleVector_tuple = read_isoform_compatibleVector_tuple0
    elif read_isoform_compatibleVector_tuple0 is not None and read_isoform_compatibleVector_tuple1 is not None:
        if sum(read_isoform_compatibleVector_tuple0[1]) == 0 or sum(read_isoform_compatibleVector_tuple1[1]) == 0:
            new_map = [max(i, j) for i, j in
                       zip(read_isoform_compatibleVector_tuple0[1], read_isoform_compatibleVector_tuple1[1])]
        else:
            new_map = [1 if i == 1 and j == 1 else 0 for (i, j) in
                       zip(read_isoform_compatibleVector_tuple0[1], read_isoform_compatibleVector_tuple1[1])]
        read_isoform_compatibleVector_tuple = (read_isoform_compatibleVector_tuple0[0], new_map)
    else:  # both are none
        read_novelisoform_tuple, read_isoform_compatibleVector_tuple = read_novelisoform_tuple0, read_isoform_compatibleVector_tuple0
    return read_novelisoform_tuple, read_isoform_compatibleVector_tuple

def map_read_to_gene(read, Info_singlegene, lowest_match=0.2, geneInfo = None,exonInfo= None, isoformInfo= None, qualifyExon= None,
                     exonMatch1= None, exonMatch2= None, pacbio = False):
    #optional: geneInfo,exonInfo, isoformInfo, qualifyExon, exonMatch
    readName, readStart, readEnd = read.qname, read.reference_start, read.reference_end
    if pacbio:
        readName = readName + '_' + str(readEnd - readStart)
    referencePositions = read.get_reference_positions(full_length=False)
    if geneInfo is None:
        geneInfo, exonInfo, isoformInfo = Info_singlegene
        isoformInfo = dict(sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))  # long to short
        geneInfo['isoformNames'] = list(isoformInfo.keys())
        qualifyExon = [i for i, (a, b) in enumerate(exonInfo) if b - a >= 20]
        exonMatch1 = [max(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]  # exon length threshold high
        exonMatch2 = [min(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]  # exon length threshold low
    read_isoform_compatibleVector = [1] * len(geneInfo['isoformNames'])#initialize read isoform compativle vector
    #read-exon mapping percentage
    mapExons = exon_hit(referencePositions, exonInfo)
    exon_map_pct, exon_map_vector_trunct = [], []
    if len(mapExons) > 0:
        exonhitbases = dict(pd.Series(mapExons).value_counts())
        exon_map_base = [exonhitbases.get(i, 0) for i in range(len(exonInfo))]
        exon_map_pct = [round(base / (exonInfo[i][1] - exonInfo[i][0]), 2) if base else 0 for i, base in enumerate(exon_map_base)]
        #read-exon mapping vector
        exon_map_vector_mapped = [1 if exon_map_p>=exon_map_up else 0 for exon_map_p, exon_map_up in list(zip(exon_map_base,exonMatch1))]
        exon_map_vector_skipped= [-1 if exon_map_p<=exon_map_low else 0 for exon_map_p, exon_map_low in list(zip(exon_map_base,exonMatch2))]
        exon_map_vector_notrunct = [mapped + skipped for mapped, skipped in list(zip(exon_map_vector_mapped, exon_map_vector_skipped))]
        if geneInfo['geneStrand']=="+":
            index_of_one = exon_map_vector_notrunct.index(1) if 1 in exon_map_vector_notrunct else len(exon_map_vector_notrunct)
            exon_map_vector_trunct = [0 if i < index_of_one else exon_map_vector_notrunct[i] for i in range(len(exon_map_vector_notrunct))]
        else: #"-"
            last_index_of_one = len(exon_map_vector_notrunct) - 1 - exon_map_vector_notrunct[::-1].index(1) if 1 in exon_map_vector_notrunct else -1
            exon_map_vector_trunct = [exon_map_vector_notrunct[i] if i <= last_index_of_one else 0 for i in range(len(exon_map_vector_notrunct))]
        #do not consider non-qualified exons
        #exon_map_vector_notrunct = [vec if i in qualifyExon else 0 for i,vec in enumerate(exon_map_vector_notrunct)]
        exon_map_vector_trunct = [vec if i in qualifyExon else 0 for i,vec in enumerate(exon_map_vector_trunct)]
        isoform_exon_map_list = isoformInfo_to_onehot(isoformInfo, geneInfo)
        read_exon_map_list = exon_map_vector_trunct
        read_isoform_compatibleVector = map_read_to_isoform_known(isoform_exon_map_list, read_exon_map_list, exon_map_pct,
                                  max(lowest_match,1-lowest_match), min(lowest_match,1-lowest_match), geneInfo['geneStrand'], qualifyExon)
    else:
        # uncharacterized type3: not map to any exons---intron
        read_isoform_compatibleVector = [0] * len(geneInfo['isoformNames'])
    #---novel/uncharacterized+existing compatible vector
    read_novelisoform_tuple, read_isoform_compatibleVector_tuple = None, None
    # novel isoform
    if sum(read_isoform_compatibleVector) == 0 and len(mapExons) > 0:
        read_novelisoform_tuple = (readName, exon_map_pct, exon_map_vector_trunct)
    else:  # existing isoform/uncategorized
        read_isoform_compatibleVector_tuple = (readName, read_isoform_compatibleVector)
    return read_novelisoform_tuple, read_isoform_compatibleVector_tuple

def Map_read_to_gene(read, Info_singlegene, lowest_match=0.05, geneInfo=None, exonInfo=None, isoformInfo=None, pacbio = False):
    if geneInfo is None:
        geneInfo, exonInfo, isoformInfo = Info_singlegene
        isoformInfo = dict(
            sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))  # long to short
        geneInfo['isoformNames'] = list(isoformInfo.keys())
    qualifyExon0 = [i for i, (a, b) in enumerate(exonInfo) if b - a >= 20]  # only exons longer than 20 will be used as proof
    exonMatch0_1 = [max(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]
    exonMatch0_2 = [min(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]
    qualifyExon1 = list(range(len(exonInfo)))
    exonMatch1_1 = [max(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]
    exonMatch1_2 = [min(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]
    # conservative match
    read_novelisoform_tuple0, read_isoform_compatibleVector_tuple0 = map_read_to_gene(read, Info_singlegene,
                                                                                      lowest_match, geneInfo,
                                                                                      exonInfo, isoformInfo,
                                                                                      qualifyExon0, exonMatch0_1, exonMatch0_2, pacbio)
    # take account of small exons
    read_novelisoform_tuple1, read_isoform_compatibleVector_tuple1 = map_read_to_gene(read, Info_singlegene,
                                                                                      lowest_match, geneInfo,
                                                                                      exonInfo,
                                                                                      isoformInfo, qualifyExon1,
                                                                                      exonMatch1_1, exonMatch1_2, pacbio)
    read_novelisoform_tuple = None
    if read_isoform_compatibleVector_tuple0 is None and read_isoform_compatibleVector_tuple1 is not None:
        read_isoform_compatibleVector_tuple = read_isoform_compatibleVector_tuple1
    elif read_isoform_compatibleVector_tuple0 is not None and read_isoform_compatibleVector_tuple1 is None:
        read_isoform_compatibleVector_tuple = read_isoform_compatibleVector_tuple0
    elif read_isoform_compatibleVector_tuple0 is not None and read_isoform_compatibleVector_tuple1 is not None:
        if sum(read_isoform_compatibleVector_tuple0[1]) == 0 or sum(read_isoform_compatibleVector_tuple1[1]) == 0:
            new_map = [max(i, j) for i, j in
                       zip(read_isoform_compatibleVector_tuple0[1], read_isoform_compatibleVector_tuple1[1])]
        else:
            new_map = [1 if i == 1 and j == 1 else 0 for (i, j) in
                       zip(read_isoform_compatibleVector_tuple0[1], read_isoform_compatibleVector_tuple1[1])]
        read_isoform_compatibleVector_tuple = (read_isoform_compatibleVector_tuple0[0], new_map)
    else:  # both are none
        read_novelisoform_tuple, read_isoform_compatibleVector_tuple = read_novelisoform_tuple0, read_isoform_compatibleVector_tuple0
    return read_novelisoform_tuple, read_isoform_compatibleVector_tuple


def map_read_to_isoform(isoformInfo_dict, isoform_assignment_vector_list, read_assignment_df):
    #isoformInfo_dict: isoform - exon annotation
    #isoform_assignment_vector_list: isoform - exon
    #read_assignment_df: read - exon
    #output: reduced isoformInfo_dict annotation; read - isoform compatible matrix, unmapped reads
    novel_isoform_assignment = np.array(isoform_assignment_vector_list)
    data = read_assignment_df.to_numpy()
    if novel_isoform_assignment.ndim == 1 or data.ndim == 1:
        return None, None, read_assignment_df.index.tolist()
    novel_ones = (novel_isoform_assignment == 1).astype(int)
    novel_minus_ones = (novel_isoform_assignment == -1).astype(int)
    ones = (data == 1).astype(int)
    minus_ones = (data == -1).astype(int)
    conflict_matrix = np.tensordot(novel_ones, minus_ones, axes=([1], [1])) + np.tensordot(novel_minus_ones, ones,axes=([1], [1]))
    compatible_matrix = (conflict_matrix == 0).astype(int).transpose()
    novel_df = pd.DataFrame(compatible_matrix)
    novel_df.index = read_assignment_df.index
    novel_df.columns = list(isoformInfo_dict.keys())
    novel_df_empty = novel_df[novel_df.sum(axis=1) == 0]
    novel_df = novel_df[novel_df.sum(axis=1) > 0]
    read_novelisoform_tuples = [(row, col) for (row, col), value in novel_df.stack().items() if value == 1]
    novel_isoform_names = list(set([b for a, b in read_novelisoform_tuples]))
    novelisoform_dict = {key: value for key, value in isoformInfo_dict.items() if key in novel_isoform_names}
    novelisoform_dict = dict(sorted(novelisoform_dict.items(), key=lambda item: -len(item[1])))
    novel_df = novel_df[novelisoform_dict.keys()]
    return novelisoform_dict, novel_df, novel_df_empty.index.tolist()

def map_read_to_isoform_known(isoform_exon_onehot_list, read_exon_map_list, exon_map_pct, threshold_high, threshold_low,
                              geneStrand, qualifyExon, poly = True):
    # isoform_exon_onehot_list: isoform - exon, one hot encoding [[-1,1,1],[1,1,1]]
    #read_exon_map_list: [0,0,1]
    #read_assignment_df: row: read - col: exon, values: 0,1,-1
    #output: reduced isoformInfo_dict annotation; read - isoform compatible matrix, unmapped reads
    def map_read_to_isoform_known_(isoform_exon_onehot_list, read_exon_map_list):
        read_exon_df = pd.DataFrame([read_exon_map_list])
        read_exon_array = read_exon_df.to_numpy()
        isoform_exon_array = np.array(isoform_exon_onehot_list)
        if isoform_exon_array.ndim == 1 or read_exon_array.ndim == 1:
            return None
        isoform_ones = (isoform_exon_array == 1).astype(int)
        isoform_minus_ones = (isoform_exon_array == -1).astype(int)
        reads_ones = (read_exon_array == 1).astype(int)
        reads_minus_ones = (read_exon_array == -1).astype(int)
        conflict_matrix = np.tensordot(isoform_ones, reads_minus_ones, axes=([1], [1])) + \
                          np.tensordot(isoform_minus_ones, reads_ones, axes=([1], [1]))
        compatible_matrix = (conflict_matrix == 0).astype(int).transpose()
        return compatible_matrix, isoform_ones
    compatible_matrix, isoform_ones = map_read_to_isoform_known_(isoform_exon_onehot_list, read_exon_map_list)
    #multiple mapping
    ####method1: treat the exon connecting poly tail as mapped
    if compatible_matrix.sum()>1:
        exon_map_pct_update = exon_map_pct
        if poly: ######nanopore and parse poly true
            if geneStrand=="+":
                #change the first non-zero element pct = 1
                for i, pct in enumerate(exon_map_pct_update):
                    if pct>0:
                        exon_map_pct_update[i] = 1
                        break
                exon_map_vector_mapped = [1 if exon_map_p >= threshold_high else 0 for exon_map_p in exon_map_pct_update]
                exon_map_vector_skipped = [-1 if exon_map_p <= threshold_low else 0 for exon_map_p in exon_map_pct_update]
                exon_map_vector_notrunct = [mapped + skipped for mapped, skipped in list(zip(exon_map_vector_mapped, exon_map_vector_skipped))]
                index_of_one = exon_map_vector_notrunct.index(1) if 1 in exon_map_vector_notrunct else len(exon_map_vector_notrunct)
                exon_map_vector_trunct = [0 if i < index_of_one else exon_map_vector_notrunct[i] for i in range(len(exon_map_vector_notrunct))]
                exon_map_vector_trunct = [vec if i in qualifyExon else 0 for i, vec in enumerate(exon_map_vector_trunct)]
            else: # geneStrand=="-":
                # change the last non-zero element pct = 1
                for i in range(len(exon_map_pct_update) - 1, -1, -1):
                    if exon_map_pct_update[i] > 0:
                        exon_map_pct_update[i] = 1
                        break
                exon_map_vector_mapped = [1 if exon_map_p >= threshold_high else 0 for exon_map_p in exon_map_pct_update]
                exon_map_vector_skipped = [-1 if exon_map_p <= threshold_low else 0 for exon_map_p in exon_map_pct_update]
                exon_map_vector_notrunct = [mapped + skipped for mapped, skipped in
                                            list(zip(exon_map_vector_mapped, exon_map_vector_skipped))]
                last_index_of_one = len(exon_map_vector_notrunct) - 1 - exon_map_vector_notrunct[::-1].index(
                    1) if 1 in exon_map_vector_notrunct else -1
                exon_map_vector_trunct = [exon_map_vector_notrunct[i] if i <= last_index_of_one else 0 for i in
                                          range(len(exon_map_vector_notrunct))]
                exon_map_vector_trunct = [vec if i in qualifyExon else 0 for i, vec in enumerate(exon_map_vector_trunct)]
            #use exon_map_vector_trunct to map again
            compatible_matrix_new, _ = map_read_to_isoform_known_(isoform_exon_onehot_list, exon_map_vector_trunct)
            if compatible_matrix_new.sum() > 0 and compatible_matrix_new.sum()<compatible_matrix.sum():
                compatible_matrix = compatible_matrix_new
                exon_map_pct = exon_map_pct_update
    ####method2: assign the cloest isoform to read
    if compatible_matrix.sum() > 1:
        read_pct = np.array(exon_map_pct)
        transformed_read_pct = np.where(read_pct > threshold_high, 1,
                                     np.where(read_pct < threshold_low, 0,
                                              (read_pct - threshold_low) / (threshold_high-threshold_low)))
        distances = abs(isoform_ones-transformed_read_pct).sum(axis = 1)
        compatible_ones = np.where(compatible_matrix == 1)[1]
        min_index = compatible_ones[np.argmin(distances[compatible_ones])]
        compatible_matrix[compatible_matrix==1] = 0
        compatible_matrix[0,min_index] = 1
    compatible_vector = compatible_matrix.tolist()[0]
    return compatible_vector


def polish_compatible_vectors(Read_novelIsoform, Read_Isoform_compatibleVector, n_isoforms):
    if len(Read_novelIsoform)==1:
        Read_Isoform_compatibleVector.append((Read_novelIsoform[0][0], [0] * n_isoforms))
        read_novelisoform_tuples = []
        novelisoform_dict ={}
    else:
        novel_dict_pct = {readname: readpct for readname, readpct, readassignment in Read_novelIsoform}
        novel_dict_assign = {readname: readassignment for readname, readpct, readassignment in Read_novelIsoform}
        df_pct = pd.DataFrame.from_dict(novel_dict_pct, orient='index')
        df_assign = pd.DataFrame.from_dict(novel_dict_assign, orient='index')
        novelisoform_dict, novel_df_list, novel_df_empty = find_novel_by_chunk(df_pct, df_assign, chunk_size=1500)
        if len(novelisoform_dict) > 0:
            read_novelisoform_tuples = [
                (row, col)
                for novel_df in novel_df_list
                for (row, col), value in novel_df.stack().items()
                if value == 1]
            for rd in novel_df_empty:
                Read_Isoform_compatibleVector.append((rd,[0]* n_isoforms))
            novelisoform_dict = dict(sorted(novelisoform_dict.items(), key=lambda item: -len(item[1])))
        else:
            for ii in range(len(Read_novelIsoform)):
                Read_Isoform_compatibleVector.append((Read_novelIsoform[ii][0], [0] * n_isoforms))
            read_novelisoform_tuples = []
            novelisoform_dict = {}
    return read_novelisoform_tuples, novelisoform_dict, Read_Isoform_compatibleVector

def compile_compatible_vectors(Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished,
                               lowest_match, geneInfo, exonInfo, Read_novelIsoform, poly=True):
    read_novelisoform_df = None
    #make novel isoform df
    if len(Read_novelIsoform_polished) > 0:
        read_poly_dict = {}
        for i, (readname, readpct, readmapping) in enumerate(Read_novelIsoform):
            if isinstance(poly, list):
                read_poly_dict[readname] = poly[i]
            else:
                read_poly_dict[readname] = poly
        read_novelisoform_df = pd.DataFrame(Read_novelIsoform_polished, columns=['read', 'isoform'])
        #identify multiple/unique mapped reads
        duplicated_reads_df = read_novelisoform_df[read_novelisoform_df.duplicated(subset=['read'], keep=False)]
        unique_reads_df = read_novelisoform_df.drop_duplicates(subset=['read'], keep=False)
        #remap duplicated reads
        read_novelisoform_df_duplicated, read_novelisoform_df_unique = None, None
        novel_isoformInfo = dict(sorted(novel_isoformInfo_polished.items(), key=lambda x: len(x[1]), reverse=True))
        if len(duplicated_reads_df)>0:
            qualifyExon = [i for i, (a, b) in enumerate(exonInfo) if b - a >= 20]
            geneStrand = geneInfo['geneStrand']
            novel_isoform_exon_map_list = isoformInfo_to_onehot(novel_isoformInfo_polished, geneInfo)
            duplicated_reads_list = list(set(duplicated_reads_df.read.tolist()))
            read_exon_map_list_all = [readmapping for readname, readpct, readmapping in Read_novelIsoform if
                              readname in duplicated_reads_list]
            exon_map_pct_all = [readpct for readname, readpct, readmapping in Read_novelIsoform if
                        readname in duplicated_reads_list]
            rows = []
            for i in range(len(duplicated_reads_list)):
                read_exon_map_list, exon_map_pct = read_exon_map_list_all[i], exon_map_pct_all[i]
                row = map_read_to_isoform_known(novel_isoform_exon_map_list, read_exon_map_list, exon_map_pct,
                                                max(lowest_match,1-lowest_match),min(lowest_match,1-lowest_match),
                                                geneStrand, qualifyExon, duplicated_reads_list[i])
                rows.append(row)
            read_novelisoform_df_duplicated = pd.DataFrame(rows)
            read_novelisoform_df_duplicated.index = duplicated_reads_list
            read_novelisoform_df_duplicated.columns = list(novel_isoformInfo_polished.keys())
            read_novelisoform_df_duplicated = read_novelisoform_df_duplicated[
                [col for col in novel_isoformInfo.keys() if col in read_novelisoform_df_duplicated.columns]]
        if len(unique_reads_df)>0:
            # manipulate unique reads
            read_novelisoform_df_unique_indicator = pd.get_dummies(unique_reads_df['isoform'])
            read_novelisoform_df_unique = pd.concat([unique_reads_df['read'], read_novelisoform_df_unique_indicator], axis=1)
            read_novelisoform_df_unique.set_index('read', inplace=True, drop=True)
            read_novelisoform_df_unique = read_novelisoform_df_unique[[col for col in novel_isoformInfo.keys() if col in read_novelisoform_df_unique.columns]]
        if read_novelisoform_df_duplicated is not None and read_novelisoform_df_unique is not None:
            read_novelisoform_df = pd.concat([read_novelisoform_df_unique,read_novelisoform_df_duplicated])
        if read_novelisoform_df_duplicated is not None and read_novelisoform_df_unique is None:
            read_novelisoform_df = read_novelisoform_df_duplicated
        if read_novelisoform_df_unique is not None and read_novelisoform_df_duplicated is None:
            read_novelisoform_df = read_novelisoform_df_unique
        #read_novelisoform_df = read_novelisoform_df.groupby(read_novelisoform_df.index).sum()
    #combine existing/uncharacterized isoforms with novel
    read_annoisoform_df = pd.DataFrame.from_dict(dict(Read_knownIsoform_polished), orient='index',columns=geneInfo['isoformNames'])
    if read_novelisoform_df is not None:
        Read_Isoform_df = read_annoisoform_df.merge(read_novelisoform_df, how='outer', left_index=True,right_index=True).fillna(0).astype(int)
    else:
        Read_Isoform_df = read_annoisoform_df
    #label uncategorized reads
    if (Read_Isoform_df.sum(axis=1) == 0).any():
        Read_Isoform_df['uncategorized'] = (Read_Isoform_df.sum(axis=1) == 0).astype(int)
    colNames = Read_Isoform_df.columns.tolist()
    Read_Isoform_compatibleVector = Read_Isoform_df.T.to_dict(orient='list')
    geneName = geneInfo['geneName']
    geneID = geneInfo['geneID']
    geneStrand = geneInfo['geneStrand']
    return geneName, geneID,geneStrand, colNames, Read_Isoform_compatibleVector



##TODO: change functions used this function: exonInfo,isoformInfo
def save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector,
                                  qname_cbumi_dict,exonInfo,isoformInfo,output_folder=None):
    #save compatible vector
    geneName = geneName.replace('/', '.')
    geneName = geneName + "_" + geneID
    print('gene ' + str(geneName) + ' processing')
    readmapping_data = []
    if Read_Isoform_compatibleVector is not None:
        # save read-isoform mapping
        for readname, indicators in Read_Isoform_compatibleVector.items():
            for idx, value in enumerate(indicators):
                if value == 1:  # Only consider mappings where the indicator is 1
                    if qname_cbumi_dict is not None:
                        cb, umi = qname_cbumi_dict[readname].split('_')[0], qname_cbumi_dict[readname].split('_')[0]
                    else:
                        cb, umi = '.', '.'
                    isoform_name = colNames[idx]
                    if isoform_name != 'uncategorized':
                        exon_indices = isoformInfo[isoform_name]
                        exon_coords = ",".join([f"({exonInfo[i][0]},{exonInfo[i][1]})" for i in exon_indices])
                        readmapping_data.append([readname, isoform_name, ','.join(map(str, exon_indices)), exon_coords, cb, umi])
                    else:
                        readmapping_data.append([readname, isoform_name, "-", "-",cb, umi])
        mat = np.array(list(dict(Read_Isoform_compatibleVector).values()))
        rowNames = list(dict(Read_Isoform_compatibleVector).keys())
        if qname_cbumi_dict is not None:
            rowNames = [qname_cbumi_dict[rn] for rn in rowNames]
        output = dict({'Gene': geneName, 'compatibleMatrix': mat, 'rowNames_cbumi': rowNames, 'colNames_isoforms': colNames})
    else:
        rowNames=[]
        output = None
    if (output_folder is not None and len(rowNames)>0):
        data_df = pd.DataFrame(output['compatibleMatrix'], index=output['rowNames_cbumi'],columns=output['colNames_isoforms'])
        data_df, novel_isoform_name_mapping = group_novel_isoform(data_df, geneStrand)
        # Save read-isoform mappings to a TSV file
        if len(readmapping_data)>0:
            output_folder0 = os.path.join(output_folder, 'auxillary')
            if not os.path.exists(output_folder0):
                os.makedirs(output_folder0)
            readmapping = pd.DataFrame(readmapping_data,
                                       columns=['Read', 'Isoform', 'Exon Index', 'Exon Coordinates', 'Cell', 'Umi'])
            readmapping_filename = os.path.join(output_folder0, geneName + '_read_isoform_exon_mapping.tsv')
            readmapping['Isoform_Name_Final'] = readmapping['Isoform'].apply(lambda isoform: novel_isoform_name_mapping.get(isoform, isoform))
            readmapping.to_csv(readmapping_filename, sep='\t', index=False)
        #save compatible matrix
        output_folder = os.path.join(output_folder,'compatible_matrix')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        file_name = os.path.join(output_folder, str(geneName) + '.csv')
        data_df.to_csv(file_name)
        print('gene ' + str(geneName) + ' saved')
    elif (output_folder is not None and len(rowNames)==0): #log the gene because no reads mapped
        output_folder = os.path.join(output_folder, 'compatible_matrix')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        log_file = os.path.join(output_folder, 'log.txt')
        if os.path.isfile(log_file):
            with open(log_file, 'a') as file:
                file.write(str(geneName)+'\n')
        else:
            with open(log_file, 'w') as file:
                file.write(str(geneName)+'\n')
        print('gene ' + str(geneName) + ' logged')
    else:
        return output


def process_read(read, geneInfo, exonInfo, isoformInfo, qname_dict, lowest_match, Info_multigenes, parse=False, pacbio = False):
    readName, readStart, readEnd = read.qname, read.reference_start, read.reference_end
    if pacbio:
        readName = readName + '_' + str(readEnd - readStart)
    novelIsoformResults = None
    isoformCompatibleVectorResults = None
    if qname_dict is None: #for bulk
        qname_dict = {}
        qname_dict[readName] = readName
    if parse:
        poly_bool, poly = detect_poly_parse(read, window=20, n=10)
        if (readStart >= geneInfo['geneStart'] and readEnd < geneInfo['geneEnd'] and readName == qname_dict[readName]):
            read_novelisoform_tuple, read_isoform_compatibleVector_tuple = Map_read_to_gene_parse(read, Info_multigenes[0], lowest_match, geneInfo, exonInfo, isoformInfo, poly_bool)
            if read_novelisoform_tuple is not None:
                novelIsoformResults = read_novelisoform_tuple
            if read_isoform_compatibleVector_tuple is not None:
                isoformCompatibleVectorResults = read_isoform_compatibleVector_tuple
    elif pacbio:
        if (readStart >= geneInfo['geneStart'] and readEnd < geneInfo['geneEnd'] and readName ==
                qname_dict[readName]):
            read_novelisoform_tuple, read_isoform_compatibleVector_tuple = Map_read_to_gene(read, Info_multigenes[0],
                                                                                            lowest_match, geneInfo,
                                                                                            exonInfo, isoformInfo, True)
            if read_novelisoform_tuple is not None:
                novelIsoformResults = read_novelisoform_tuple
            if read_isoform_compatibleVector_tuple is not None:
                isoformCompatibleVectorResults = read_isoform_compatibleVector_tuple
    else:
        poly_bool, poly = detect_poly(read, window=20, n=15)
        if (readStart >= geneInfo['geneStart'] and readEnd < geneInfo['geneEnd'] and poly_bool and readName == qname_dict[readName]):
            read_novelisoform_tuple, read_isoform_compatibleVector_tuple = Map_read_to_gene(read, Info_multigenes[0], lowest_match, geneInfo, exonInfo, isoformInfo, False)
            if read_novelisoform_tuple is not None:
                novelIsoformResults = read_novelisoform_tuple
            if read_isoform_compatibleVector_tuple is not None:
                isoformCompatibleVectorResults = read_isoform_compatibleVector_tuple
    return novelIsoformResults, isoformCompatibleVectorResults
def process_read_metagene(read, start, end, qname_dict, Info_multigenes, lowest_match, parse=False, pacbio = False):
    readName, readStart, readEnd = read.qname, read.reference_start, read.reference_end
    if qname_dict is None: #for bulk
        qname_dict = {}
        qname_dict[readName] = readName
    if pacbio:
        readName = readName + '_' + str(readEnd - readStart)
    if parse:
        poly_bool, poly = detect_poly_parse(read, window=20, n=10)
        if (readStart >= start and readEnd < end and readName == qname_dict[readName]):
            ind, read_novelisoform_tuple, read_isoform_compatibleVector_tuple = choose_gene_from_meta_parse(read,
                                                                                                      Info_multigenes,
                                                                                                      lowest_match,poly_bool)
            if ind >= 0:
                return ind, read_novelisoform_tuple, read_isoform_compatibleVector_tuple
    elif pacbio:
        if (readStart >= start and readEnd < end and readName == qname_dict[readName]):
            ind, read_novelisoform_tuple, read_isoform_compatibleVector_tuple = choose_gene_from_meta(read,
                                                                                                      Info_multigenes,
                                                                                                      lowest_match, pacbio=True)
            if ind >= 0:
                return ind, read_novelisoform_tuple, read_isoform_compatibleVector_tuple
    else:
        poly_bool, poly = detect_poly(read, window=20, n=15)
        if (readStart >= start and readEnd < end and poly_bool and readName == qname_dict[readName]):
            ind, read_novelisoform_tuple, read_isoform_compatibleVector_tuple = choose_gene_from_meta(read, Info_multigenes, lowest_match,pacbio=False)
            if ind >= 0:
                return ind, read_novelisoform_tuple, read_isoform_compatibleVector_tuple
    return None



#-------------------------------end----------------------------#





#####some functions to delete ########










