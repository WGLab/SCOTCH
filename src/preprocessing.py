import pysam
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
import pickle
from itertools import chain
import numpy as np
import os
from Bio.Seq import Seq
import math
import reference as ref
import networkx as nx
import community.community_louvain as community_louvain
from collections import defaultdict


#bam="/scr1/users/xu3/singlecell/project_singlecell/sample8_R10/bam/sample8_R10.filtered.bam"
#gene_pkl="/scr1/users/xu3/singlecell/project_singlecell/M4/reference/geneStructureInformation.pkl"



#extract read information
def read_to_dict(aligned_segment):
    # Create a dictionary of the aligned segment data you want to keep
    data = {
        'qname': aligned_segment.qname,
        'qstart': aligned_segment.qstart,
        'qend': aligned_segment.qend,
        'reference_start': aligned_segment.reference_start,
        'reference_end': aligned_segment.reference_end,
        'get_reference_positions': aligned_segment.get_reference_positions(full_length=False),
        'query_sequence': aligned_segment.query_sequence,
        'umi': aligned_segment.get_tag('UR')
    }
    return data

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

def load_pickle(file):
    with open(file,'rb') as file:
        data=pickle.load(file)
    return data

def unpack_list(x):
    return list(chain(*x))

def iter_length(iter):
    return sum(1 for _ in iter)

def novelid_to_exonid(novelid):
    exonid = []
    position = 0
    while novelid > 0:
        if novelid & 1:
            exonid.append(position)
        novelid >>= 1
        position += 1
    return exonid

def merge_exons(exons):
    merged_exons = []
    current_start, current_end = exons[0]
    for start, end in exons[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged_exons.append((current_start, current_end))
            current_start, current_end = start, end
    merged_exons.append((current_start, current_end))
    return merged_exons

def extract_bam_info_folder(bam_folder, num_cores, parse=False, pacbio = False):
    files = os.listdir(bam_folder)
    bamfiles = [os.path.join(bam_folder,f) for f in files if f.endswith('.bam')]
    if parse:
        df = Parallel(n_jobs=num_cores)(delayed(extract_bam_info_parse)(bam) for bam in bamfiles)
    elif pacbio:
        df = Parallel(n_jobs=num_cores)(delayed(extract_bam_info_pacbio)(bam) for bam in bamfiles)
    else:
        df = Parallel(n_jobs=num_cores)(delayed(extract_bam_info)(bam) for bam in bamfiles)
    ReadTagsDF = pd.concat(df).reset_index(drop=True)
    return ReadTagsDF

def extract_bam_info(bam):
    #extract readname, cb, umi from bam file
    # bam: path to bam file
    bamFilePysam = pysam.Samfile(bam, "rb")
    ReadTags = [(read.qname, read.get_tag('CB'), read.get_tag('UB'), read.qend-read.qstart) for read in bamFilePysam]
    ReadTagsDF = pd.DataFrame(ReadTags)
    if ReadTagsDF.shape[0]>0:
        ReadTagsDF.columns = ['QNAME', 'CB', 'UMI', 'LENGTH']
        ReadTagsDF = ReadTagsDF.sort_values(by=['CB','UMI','LENGTH'],ascending=[True, True, False]).reset_index(drop=True)
        ReadTagsDF['CBUMI']=ReadTagsDF.CB.astype(str)+'_'+ReadTagsDF.UMI.astype(str)
    else:
        ReadTagsDF = None
    return ReadTagsDF

def extract_bam_info_parse(bam):
    bamFilePysam = pysam.Samfile(bam, "rb")
    #qname cb umi cbumi length
    ReadTags = [(read.qname, read.qname.split('_')[-5]+'_'+read.qname.split('_')[-4]+'_'+read.qname.split('_')[-3], read.qname.split('_')[-1] , len(read.query_alignment_sequence), read.get_tag('pS')) for read in bamFilePysam]
    ReadTagsDF = pd.DataFrame(ReadTags)
    if ReadTagsDF.shape[0] > 0:
        ReadTagsDF.columns = ['QNAME', 'CB', 'UMI', 'LENGTH','SAMPLE']
        ReadTagsDF = ReadTagsDF.sort_values(by=['SAMPLE','CB', 'UMI', 'LENGTH'], ascending=[True, True, True, False]).reset_index(drop=True)
        ReadTagsDF['CBUMI'] = ReadTagsDF.CB.astype(str) + '_' + ReadTagsDF.UMI.astype(str)
    else:
        ReadTagsDF = None
    return ReadTagsDF

def extract_bam_info_pacbio(bam):
    bamFilePysam = pysam.Samfile(bam, "rb")
    #qname cb umi cbumi length
    ReadTags = [(read.qname, read.get_tag('CB'), read.get_tag('XM'), read.reference_end-read.reference_start) for read in bamFilePysam]
    ReadTagsDF = pd.DataFrame(ReadTags)
    if ReadTagsDF.shape[0] > 0:
        ReadTagsDF.columns = ['QNAME', 'CB', 'UMI', 'LENGTH']
        ReadTagsDF = ReadTagsDF.sort_values(by=['CB','UMI','LENGTH'],ascending=[True, True, False]).reset_index(drop=True)
        ReadTagsDF['CBUMI'] = ReadTagsDF.CB.astype(str) + '_' + ReadTagsDF.UMI.astype(str)
        ReadTagsDF['QNAME'] = ReadTagsDF.QNAME.astype(str) + '_' + ReadTagsDF.LENGTH.astype(str)
    else:
        ReadTagsDF = None
    return ReadTagsDF



def bam_info_to_dict(bam_info, parse=False):
    #input bam_info is a dataframe(bam_info is sorted already when generating)
    #the output dict qname_dict: key: qname, value: the right qname to keep; qname_CBUMI_dict: qname_cbumi
    #generate qname_dict
    max_length_df = bam_info.drop_duplicates(subset='CBUMI', keep='first')
    cbumi_to_max_qname = pd.Series(max_length_df.QNAME.values, index=max_length_df.CBUMI).to_dict()
    qname_dict = {row['QNAME']: cbumi_to_max_qname[row['CBUMI']] for index, row in bam_info.iterrows()}
    #generate qname_cbumi_dict
    # Assuming 'bam_info' is your pandas DataFrame
    qname_cbumi_dict = dict(zip(bam_info['QNAME'], bam_info['CBUMI']))
    if parse:
        qname_sample_dict = dict(zip(bam_info['QNAME'], bam_info['SAMPLE']))
        return qname_dict, qname_cbumi_dict, qname_sample_dict
    else:
        return qname_dict, qname_cbumi_dict

def process_gene(geneID, geneName, genes ,exons, build=None):
    #gtf file is 1 based, change annotation to 0 based
    GeneDf = genes[genes.iloc[:, 3] == geneID].reset_index(drop=True)
    ExonsDf = exons[exons.iloc[:, 3] == geneID].reset_index(drop=True)
    exonsdf = ExonsDf.iloc[:, :3].drop_duplicates().reset_index(drop=True)
    exonsdf['exon_num'] = list(range(exonsdf.shape[0]))
    geneChr, geneStart, geneEnd = GeneDf.iloc[0, 0], GeneDf.iloc[0, 1], GeneDf.iloc[0, 2]
    if build is not None:
        geneChr = str(build) + '_'+str(geneChr)
    geneStrand = GeneDf.iloc[0, 6]
    isoformNames = ExonsDf.TRANSCRIPT.unique().tolist()
    GeneInfo = {'geneName': geneName, 'geneID': geneID, 'geneChr': geneChr, 'geneStart': geneStart -1 , 'geneEnd': geneEnd,
                'geneStrand':geneStrand, 'numofExons': exonsdf.shape[0], 'numofIsoforms': len(isoformNames),
                'isoformNames': isoformNames}
    ExonPositions = list(zip(exonsdf.iloc[:, 1] -1 , exonsdf.iloc[:, 2]))
    ExonsDf = ExonsDf.merge(exonsdf, how='left')
    ExonIsoformDict = dict()
    for isoform in isoformNames:
        ExonIsoformDict[isoform] = ExonsDf[ExonsDf.TRANSCRIPT == isoform].exon_num.tolist()
    return geneID, [GeneInfo, ExonPositions, ExonIsoformDict]

def extract_annotation_info(refGeneFile, num_cores=8, output="geneStructureInformation.pkl", build=None):
    metageneStructureInformation = None
    meta_output = os.path.join(os.path.dirname(output),'meta'+os.path.basename(output))
    genes, exons = ref.generate_reference_df(gtf_path = refGeneFile)
    Genes = list(zip(genes.iloc[:,3].tolist(), genes.iloc[:,4].tolist()))#id, name
    if os.path.isfile(output)==False:
        geneStructureInformation = Parallel(n_jobs=num_cores)(delayed(process_gene)(geneID, geneName,genes ,exons, build) for geneID, geneName in tqdm(Genes))
        geneStructureInformation = dict(geneStructureInformation)
        if output is not None:
            with open(output, 'wb') as file:
                pickle.dump(geneStructureInformation, file)
    else:
        geneStructureInformation = load_pickle(output)
    if os.path.isfile(meta_output)==False:
        #group genes
        genes.columns = ['CHR', 'START', 'END', 'GENE_ID', 'GENE_NAME', 'GENE_TYPE', 'STRAND' ,'META_GENE']
        grouped = genes.groupby("META_GENE")
        grouped_dict = {key: group.iloc[:, 3].tolist() for key, group in grouped}  # META_GENE STARTS FROM 1
        metageneStructureInformation = {}
        for i in range(1, 1 + len(grouped_dict)):
            meta_gene = 'meta_gene_' + str(i)
            gene_ids = grouped_dict[i]
            meta_gene_info = []
            for id in gene_ids:
                meta_gene_info.append(geneStructureInformation[id])
            metageneStructureInformation[meta_gene] = meta_gene_info
        with open(meta_output, 'wb') as file:
            pickle.dump(metageneStructureInformation, file)
    return metageneStructureInformation

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

def read_exon_match(read, Info_singlegene):
    if isinstance(read, dict):
        readName, readStart, readEnd = read['qname'], read['reference_start'], read['reference_end']
        referencePositions = read['get_reference_positions']
    else:
        readName, readStart, readEnd = read.qname, read.reference_start, read.reference_end
        referencePositions = read.get_reference_positions(full_length=False)
    geneInfo, exonInfo, isoformInfo = Info_singlegene
    gene_name =geneInfo['geneName']
    gene_length = geneInfo['geneEnd'] - geneInfo['geneStart']
    if geneInfo['geneStrand']=="+":
        readDist = readStart - geneInfo["geneStart"]
    else:
        readDist = readEnd - geneInfo["geneEnd"]
    if (readStart >= geneInfo['geneStart'] and readEnd < geneInfo['geneEnd']):
        mapExons = exon_hit(referencePositions, exonInfo)
        n_mapExons=sum(mapExons)
    else:
        n_mapExons = -1
    return gene_name, readDist, gene_length, n_mapExons


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

def map_read_to_gene_parse(read, Info_singlegene, lowest_match=0.05, geneInfo = None,exonInfo= None, isoformInfo= None,
                           qualifyExon= None, exonMatch1= None, exonMatch2= None, poly=False):
    #optional: geneInfo,exonInfo, isoformInfo, qualifyExon, exonMatch
    if isinstance(read, dict):
        readName = read['qname']
        referencePositions = read['get_reference_positions']
    else:
        readName = read.qname
        referencePositions = read.get_reference_positions(full_length=False)
    if geneInfo is None:
        geneInfo, exonInfo, isoformInfo = Info_singlegene
        isoformInfo = dict(sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))  # long to short
        geneInfo['isoformNames'] = list(isoformInfo.keys())
        qualifyExon = [i for i, (a, b) in enumerate(exonInfo) if b - a >= 20]
        exonMatch1 = [max(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]  # exon length threshold
        exonMatch2 = [min(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]  # exon length threshold
    read_isoform_compatibleVector = [1] * len(geneInfo['isoformNames'])#initialize read isoform compativle vector
    mapExons = exon_hit(referencePositions, exonInfo)
    #INITIALIZATION
    #---existing annotation compatible vector
    exonCovered_qualify1,exonCovered_qualify2=[],[]
    if len(mapExons) > 0:
        #high
        exonhitbases1 = dict(pd.Series(mapExons).value_counts())
        exonhitbases1 = remove_dict_items(exonhitbases1, exonMatch1)
        exonCovered1 = list(exonhitbases1.keys()) # exon index the read mapped to
        #low
        exonhitbases2 = dict(pd.Series(mapExons).value_counts())
        exonhitbases2 = remove_dict_items(exonhitbases2, exonMatch2)
        exonCovered2 = list(exonhitbases2.keys())  # exon index the read mapped to
        rule_out = []
        isoform_exons = list(isoformInfo.values())
        exonCovered_qualify1 = [ec for ec in exonCovered1 if ec in qualifyExon]#used for infer covered exons
        exonCovered_qualify2 = [ec for ec in exonCovered2 if ec in qualifyExon]#used for infer skipped exons
        ##----rule1: reads covers exonA, rule out isoforms not covering the exonA based on annotation
        if len(exonCovered_qualify1) > 0:
            for ind, isoform_exons_ in enumerate(isoform_exons):
                # rule out the isoforms do not contain exons that covered
                if check_values_in_list(exonCovered_qualify1, isoform_exons_, 'all') == False:
                    rule_out.append(ind)
        ##----rule2: reads skips exonA, rule out isoforms covering the exonA based on annotation
        if len(exonCovered_qualify2) >0:
            exonSkipped = [v for v in list(range(min(exonCovered2), max(exonCovered2) + 1)) if v not in exonCovered2]
            if poly:
                if geneInfo['geneStrand']=="+":# 5`truncation
                    exonSkipped = exonSkipped + list(range(max(exonCovered2) + 1, geneInfo['numofExons']))
                    exonSkipped_qualify = [ec for ec in exonSkipped if ec in qualifyExon]
                else:  # 3`truncation; - strand
                    exonSkipped = exonSkipped + list(range(min(exonCovered2)))
                    exonSkipped_qualify = [ec for ec in exonSkipped if ec in qualifyExon]
            else:
                exonSkipped_qualify = [ec for ec in exonSkipped if ec in qualifyExon]
            if len(exonSkipped_qualify) > 0:
                for ind, isoform_exons_ in enumerate(isoform_exons):
                    if check_values_in_list(exonSkipped_qualify, isoform_exons_, 'any') == True:
                        rule_out.append(ind)
        else:  # uncharacterized type2:  exonCovered_qualify==0;mapped to exons but less than the threshold
            read_isoform_compatibleVector = [0] * len(geneInfo['isoformNames'])
        #for reads covered qualified exons, update read_isoform_compatibleVector based on rule_out
        if len(rule_out) > 0:
            for ii in rule_out:
                read_isoform_compatibleVector[ii] = 0
    else:
        # uncharacterized type3: not map to any exons---intron
        read_isoform_compatibleVector = [0] * len(geneInfo['isoformNames'])
    #---novel/uncharacterized+existing compatible vector
    read_novelisoform_tuple, read_isoform_compatibleVector_tuple = None, None
    if sum(read_isoform_compatibleVector) == 0 and len(exonCovered_qualify2)>0: # novel isoform
        exonCovered_novel = [0] * len(exonInfo)
        for i in exonSkipped:
            exonCovered_novel[i] = -1
        for i in exonCovered_qualify1:
            exonCovered_novel[i] = 1
        exonhitbases = dict(pd.Series(mapExons).value_counts())
        P = []
        for k in range(len(exonInfo)):
            if k in list(exonhitbases.keys()):
                p = round(exonhitbases[k] / (exonInfo[k][1] - exonInfo[k][0]), 2)
            else:
                p = 0
            P.append(p)
        read_novelisoform_tuple = (readName, P, exonCovered_novel)
    else:  # existing isoform/uncategorized
        read_isoform_compatibleVector_tuple = (readName, read_isoform_compatibleVector)
    return read_novelisoform_tuple, read_isoform_compatibleVector_tuple

def Map_read_to_gene_parse(read, Info_singlegene, lowest_match=0.05, geneInfo=None, exonInfo=None, isoformInfo=None, poly=False):
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

def map_read_to_gene(read, Info_singlegene, lowest_match=0.05, geneInfo = None,exonInfo= None, isoformInfo= None, qualifyExon= None, exonMatch1= None, exonMatch2= None, pacbio = False):
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
        exonMatch1 = [max(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]  # exon length threshold
        exonMatch2 = [min(lowest_match,1-lowest_match) * (b - a) for a, b in exonInfo]  # exon length threshold
    read_isoform_compatibleVector = [1] * len(geneInfo['isoformNames'])#initialize read isoform compativle vector
    mapExons = exon_hit(referencePositions, exonInfo)
    #INITIALIZATION
    #---existing annotation compatible vector
    exonCovered_qualify1,exonCovered_qualify2=[],[]
    if len(mapExons) > 0:
        #high
        exonhitbases1 = dict(pd.Series(mapExons).value_counts())
        exonhitbases1 = remove_dict_items(exonhitbases1, exonMatch1)
        exonCovered1 = list(exonhitbases1.keys()) # exon index the read mapped to
        #low
        exonhitbases2 = dict(pd.Series(mapExons).value_counts())
        exonhitbases2 = remove_dict_items(exonhitbases2, exonMatch2)
        exonCovered2 = list(exonhitbases2.keys())  # exon index the read mapped to
        rule_out = []
        isoform_exons = list(isoformInfo.values())
        exonCovered_qualify1 = [ec for ec in exonCovered1 if ec in qualifyExon]#used for infer covered exons
        exonCovered_qualify2 = [ec for ec in exonCovered2 if ec in qualifyExon]#used for infer skipped exons
        ##----rule1: reads covers exonA, rule out isoforms not covering the exonA based on annotation
        if len(exonCovered_qualify1) > 0:
            for ind, isoform_exons_ in enumerate(isoform_exons):
                # rule out the isoforms do not contain exons that covered
                if check_values_in_list(exonCovered_qualify1, isoform_exons_, 'all') == False:
                    rule_out.append(ind)
        ##----rule2: reads skips exonA, rule out isoforms covering the exonA based on annotation
        if len(exonCovered_qualify2) >0:
            exonSkipped = [v for v in list(range(min(exonCovered2), max(exonCovered2) + 1)) if v not in exonCovered2]
            if geneInfo['geneStrand']=="+":# 5`truncation
                exonSkipped = exonSkipped + list(range(max(exonCovered2) + 1, geneInfo['numofExons']))
                exonSkipped_qualify = [ec for ec in exonSkipped if ec in qualifyExon]
            else:  # 3`truncation; - strand
                exonSkipped = exonSkipped + list(range(min(exonCovered2)))
                exonSkipped_qualify = [ec for ec in exonSkipped if ec in qualifyExon]
            if len(exonSkipped_qualify) > 0:
                for ind, isoform_exons_ in enumerate(isoform_exons):
                    if check_values_in_list(exonSkipped_qualify, isoform_exons_, 'any') == True:
                        rule_out.append(ind)
        else:  # uncharacterized type2:  exonCovered_qualify==0;mapped to exons but less than the threshold
            read_isoform_compatibleVector = [0] * len(geneInfo['isoformNames'])
        #for reads covered qualified exons, update read_isoform_compatibleVector based on rule_out
        if len(rule_out) > 0:
            for ii in rule_out:
                read_isoform_compatibleVector[ii] = 0
    else:
        # uncharacterized type3: not map to any exons---intron
        read_isoform_compatibleVector = [0] * len(geneInfo['isoformNames'])
    #---novel/uncharacterized+existing compatible vector
    read_novelisoform_tuple, read_isoform_compatibleVector_tuple = None, None
    if sum(read_isoform_compatibleVector) == 0 and len(exonCovered_qualify2)>0: # novel isoform
        exonCovered_novel = [0] * len(exonInfo)
        for i in exonSkipped:
            exonCovered_novel[i] = -1
        for i in exonCovered_qualify1:
            exonCovered_novel[i] = 1
        exonhitbases = dict(pd.Series(mapExons).value_counts())
        P = []
        for k in range(len(exonInfo)):
            if k in list(exonhitbases.keys()):
                p = round(exonhitbases[k] / (exonInfo[k][1] - exonInfo[k][0]), 2)
            else:
                p = 0
            P.append(p)
        read_novelisoform_tuple = (readName, P, exonCovered_novel)
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



def map_read_to_isoform(isoformInfo_dict, isoform_assignment_vector_list, read_assignment_df):
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
    #for rd in novel_df_empty.index.tolist():
    #    Read_Isoform_compatibleVector.append((rd, [0] * n_isoforms))
    return novelisoform_dict, novel_df, novel_df_empty.index.tolist()

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



def compile_compatible_vectors(Read_novelIsoform, novel_isoformInfo, Read_Isoform_compatibleVector, geneInfo):
    #this function is for all reads mapped to one gene
    #input are output of map_read_to_gene
    #novel isoforms
    read_novelisoform_df = None
    #make novel isoform df
    if len(Read_novelIsoform) > 0:
        read_novelisoform_df = pd.DataFrame(Read_novelIsoform, columns=['read', 'isoform'])
        read_novelisoform_df_indicator = pd.get_dummies(read_novelisoform_df['isoform'])
        read_novelisoform_df = pd.concat([read_novelisoform_df['read'], read_novelisoform_df_indicator], axis=1)
        read_novelisoform_df.set_index('read', inplace=True, drop=True)
        novel_isoformInfo = dict(sorted(novel_isoformInfo.items(), key=lambda x: len(x[1]), reverse=True))
        read_novelisoform_df = read_novelisoform_df[list(novel_isoformInfo.keys())]
        read_novelisoform_df = read_novelisoform_df.groupby(read_novelisoform_df.index).sum()
    #combine existing/uncharacterized isoforms with novel
    read_annoisoform_df = pd.DataFrame.from_dict(dict(Read_Isoform_compatibleVector), orient='index',columns=geneInfo['isoformNames'])
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


def save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector,qname_cbumi_dict, output_folder=None):
    geneName = geneName.replace('/', '.')
    geneName = geneName + "_" + geneID
    print('gene ' + str(geneName) + ' processing')
    if Read_Isoform_compatibleVector is not None:
        mat = np.array(list(dict(Read_Isoform_compatibleVector).values()))
        rowNames = list(dict(Read_Isoform_compatibleVector).keys())
        rowNames = [qname_cbumi_dict[rn] for rn in rowNames]
        output = dict({'Gene': geneName, 'compatibleMatrix': mat, 'rowNames_cbumi': rowNames, 'colNames_isoforms': colNames})
    else:
        rowNames=[]
        output = None
    if (output_folder is not None and len(rowNames)>0):
        data_df = pd.DataFrame(output['compatibleMatrix'], index=output['rowNames_cbumi'],columns=output['colNames_isoforms'])
        data_df = group_novel_isoform(data_df, geneStrand)
        output_folder = os.path.join(output_folder,'compatible_matrix')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        file_name = os.path.join(output_folder, str(geneName) + '.csv')
        data_df.to_csv(file_name)
        print('gene ' + str(geneName) + ' saved')
    elif (output_folder is not None and len(rowNames)==0):
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


def extract_compatible_matrix_by_metagene(bamFile, meta_gene_pkl, meta_gene, qname_dict, qname_cbumi_dict, lowest_match=0.2, output_folder=None):
    #extract compatible matrix by metagene: output gene name and compatible matrix
    #bamFile: path to bam file
    #gene_pkl: path to gene annotation pkl file, or a geneStructureInformation object
    print(meta_gene)
    if isinstance(meta_gene_pkl,str):
        metageneStructureInformation = load_pickle(meta_gene_pkl)
    else:
        metageneStructureInformation = meta_gene_pkl
    Info_multigenes = metageneStructureInformation[meta_gene]#list: length is gene number
    Info_multigenes_sort = []
    for info in Info_multigenes:
        isoformInfo = info[2]
        isoformInfo = dict(sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))
        info[0]['isoformNames'] = list(isoformInfo.keys())
        info[2] = isoformInfo
        Info_multigenes_sort.append(info)
    Info_multigenes = Info_multigenes_sort
    if os.path.isfile(bamFile)==False: #bamFile is a folder
        bamFile_name = [f for f in os.listdir(bamFile) if f.endswith('.bam') and '.'+Info_multigenes[0][0]['geneChr']+'.' in f]
        bamFile = os.path.join(bamFile,bamFile_name[0])
    bamFilePysam = pysam.Samfile(bamFile, "rb")
    ##########-----only contain one gene for metagene
    if len(Info_multigenes)==1:
        geneInfo, exonInfo, isoformInfo = Info_multigenes[0]
        n_isoforms = len(isoformInfo)
        reads = bamFilePysam.fetch(geneInfo['geneChr'], geneInfo['geneStart'], geneInfo['geneEnd'])
        Read_Isoform_compatibleVector, Read_novelIsoform, novel_isoformInfo = [], [], {}
        for read in reads:
            result = process_read(read, geneInfo, exonInfo, isoformInfo, qname_dict, lowest_match, Info_multigenes)
            novelIsoformResults, isoformCompatibleVectorResults = result
            if novelIsoformResults is not None:
                Read_novelIsoform.append(novelIsoformResults)
            if isoformCompatibleVectorResults is not None:
                Read_Isoform_compatibleVector.append(isoformCompatibleVectorResults)
        if len(Read_novelIsoform) > 0:
            Read_novelIsoform, novel_isoformInfo, Read_Isoform_compatibleVector = polish_compatible_vectors(
                Read_novelIsoform, Read_Isoform_compatibleVector, n_isoforms)
        #loop finished, start compiling into compatible matrix
        geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(Read_novelIsoform, novel_isoformInfo, Read_Isoform_compatibleVector, geneInfo)
        save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector, qname_cbumi_dict, output_folder)
    ############-----multiple genes for one metagene
    else:
        geneChr, start, end = summarise_metagene(Info_multigenes)  # meta gene location
        reads = bamFilePysam.fetch(geneChr, start, end)  # fetch reads within meta gene region
        #process reads metagene
        #read_dicts = [read_to_dict(read) for read in reads]
        results = []
        for read in reads:
            out = process_read_metagene(read, start, end, qname_dict, Info_multigenes, lowest_match)
            if out is not None:
                results.append(out)
        Ind, read_novelisoform_tuple_metagene, read_isoform_compatibleVector_tuple_metagene = [], [], []
        for result in results:
            if result is not None:
                ind, novelisoform_tuple, isoformCompatibleVector_tuple = result
                Ind.append(ind)
                read_novelisoform_tuple_metagene.append(novelisoform_tuple)
                read_isoform_compatibleVector_tuple_metagene.append(isoformCompatibleVector_tuple)
        unique_ind = list(set(Ind))
        log_ind = [ind for ind in range(len(Info_multigenes)) if ind not in unique_ind]
        #logging genes without any reads
        for index in log_ind:
            save_compatibleVector_by_gene(Info_multigenes[index][0]['geneName'], Info_multigenes[index][0]['geneID'], Info_multigenes[index][0]['geneStrand'], colNames=None,
                                          Read_Isoform_compatibleVector=None, qname_cbumi_dict=None,output_folder=output_folder)
        for index in unique_ind:
            print('processing gene'+str(index))
            #loop over genes within metagene; for one single gene:
            Read_novelIsoform, Read_Isoform_compatibleVector, novel_isoformInfo=[],[], {}
            for j, i in enumerate(Ind):
                    #loop for reads
                if i==index and read_novelisoform_tuple_metagene[j] is not None:
                    Read_novelIsoform.append(read_novelisoform_tuple_metagene[j])
                if i==index and read_isoform_compatibleVector_tuple_metagene[j] is not None:
                    Read_Isoform_compatibleVector.append(read_isoform_compatibleVector_tuple_metagene[j])
            if len(Read_novelIsoform) > 0:
                Read_novelIsoform, novel_isoformInfo, Read_Isoform_compatibleVector = polish_compatible_vectors(
                    Read_novelIsoform, Read_Isoform_compatibleVector, len(Info_multigenes[index][2]))
            #loop over reads for one gene
            geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(Read_novelIsoform,
                                                                                               novel_isoformInfo,
                                                                                               Read_Isoform_compatibleVector,
                                                                                               Info_multigenes[index][0])
            save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector, qname_cbumi_dict, output_folder)
    #finish mapping reads to genes
    bamFilePysam.close()

def extract_compatible_matrix_by_metagene_pacbio(bamFile, meta_gene_pkl, meta_gene, qname_dict, qname_cbumi_dict, lowest_match=0.2, output_folder=None):
    #extract compatible matrix by metagene: output gene name and compatible matrix
    #bamFile: path to bam file
    #gene_pkl: path to gene annotation pkl file, or a geneStructureInformation object
    print(meta_gene)
    if isinstance(meta_gene_pkl,str):
        metageneStructureInformation = load_pickle(meta_gene_pkl)
    else:
        metageneStructureInformation = meta_gene_pkl
    Info_multigenes = metageneStructureInformation[meta_gene]#list: length is gene number
    Info_multigenes_sort = []
    for info in Info_multigenes:
        isoformInfo = info[2]
        isoformInfo = dict(sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))
        info[0]['isoformNames'] = list(isoformInfo.keys())
        info[2] = isoformInfo
        Info_multigenes_sort.append(info)
    Info_multigenes = Info_multigenes_sort
    if os.path.isfile(bamFile)==False: #bamFile is a folder
        bamFile_name = [f for f in os.listdir(bamFile) if f.endswith('.bam') and '.'+Info_multigenes[0][0]['geneChr']+'.' in f]
        bamFile = os.path.join(bamFile,bamFile_name[0])
    bamFilePysam = pysam.Samfile(bamFile, "rb")
    ##########-----only contain one gene for metagene
    if len(Info_multigenes)==1:
        geneInfo, exonInfo, isoformInfo = Info_multigenes[0]
        n_isoforms = len(isoformInfo)
        reads = bamFilePysam.fetch(geneInfo['geneChr'], geneInfo['geneStart'], geneInfo['geneEnd'])
        Read_Isoform_compatibleVector, Read_novelIsoform, novel_isoformInfo = [], [], {}
        for read in reads:
            result = process_read(read, geneInfo, exonInfo, isoformInfo, qname_dict, lowest_match, Info_multigenes, pacbio=True)
            novelIsoformResults, isoformCompatibleVectorResults = result
            if novelIsoformResults is not None:
                Read_novelIsoform.append(novelIsoformResults)
            if isoformCompatibleVectorResults is not None:
                Read_Isoform_compatibleVector.append(isoformCompatibleVectorResults)
        if len(Read_novelIsoform) > 0:
            Read_novelIsoform, novel_isoformInfo, Read_Isoform_compatibleVector = polish_compatible_vectors(
                Read_novelIsoform, Read_Isoform_compatibleVector, n_isoforms)
        #loop finished, start compiling into compatible matrix
        geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(Read_novelIsoform, novel_isoformInfo, Read_Isoform_compatibleVector, geneInfo)
        save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector, qname_cbumi_dict, output_folder)
    ############-----multiple genes for one metagene
    else:
        geneChr, start, end = summarise_metagene(Info_multigenes)  # meta gene location
        reads = bamFilePysam.fetch(geneChr, start, end)  # fetch reads within meta gene region
        #process reads metagene
        #read_dicts = [read_to_dict(read) for read in reads]
        results = []
        for read in reads:
            out = process_read_metagene(read, start, end, qname_dict, Info_multigenes, lowest_match, pacbio = True)
            if out is not None:
                results.append(out)
        Ind, read_novelisoform_tuple_metagene, read_isoform_compatibleVector_tuple_metagene = [], [], []
        for result in results:
            if result is not None:
                ind, novelisoform_tuple, isoformCompatibleVector_tuple = result
                Ind.append(ind)
                read_novelisoform_tuple_metagene.append(novelisoform_tuple)
                read_isoform_compatibleVector_tuple_metagene.append(isoformCompatibleVector_tuple)
        unique_ind = list(set(Ind))
        log_ind = [ind for ind in range(len(Info_multigenes)) if ind not in unique_ind]
        #logging genes without any reads
        for index in log_ind:
            save_compatibleVector_by_gene(Info_multigenes[index][0]['geneName'], Info_multigenes[index][0]['geneID'], Info_multigenes[index][0]['geneStrand'], colNames=None,
                                          Read_Isoform_compatibleVector=None, qname_cbumi_dict=None,output_folder=output_folder)
        for index in unique_ind:
            print('processing gene'+str(index))
            #loop over genes within metagene; for one single gene:
            Read_novelIsoform, Read_Isoform_compatibleVector, novel_isoformInfo=[],[], {}
            for j, i in enumerate(Ind):
                    #loop for reads
                if i==index and read_novelisoform_tuple_metagene[j] is not None:
                    Read_novelIsoform.append(read_novelisoform_tuple_metagene[j])
                if i==index and read_isoform_compatibleVector_tuple_metagene[j] is not None:
                    Read_Isoform_compatibleVector.append(read_isoform_compatibleVector_tuple_metagene[j])
            if len(Read_novelIsoform) > 0:
                Read_novelIsoform, novel_isoformInfo, Read_Isoform_compatibleVector = polish_compatible_vectors(
                    Read_novelIsoform, Read_Isoform_compatibleVector, len(Info_multigenes[index][2]))
            #loop over reads for one gene
            geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(Read_novelIsoform,
                                                                                               novel_isoformInfo,
                                                                                               Read_Isoform_compatibleVector,
                                                                                               Info_multigenes[index][0])
            save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector, qname_cbumi_dict, output_folder)
    #finish mapping reads to genes
    bamFilePysam.close()

def extract_compatible_matrix_by_metagene_parse(bamFile, meta_gene_pkl, meta_gene, qname_dict, qname_cbumi_dict, qname_sample_dict, lowest_match=0.2, output_folder=None):
    #extract compatible matrix by metagene: output gene name and compatible matrix
    #bamFile: path to bam file
    #gene_pkl: path to gene annotation pkl file, or a geneStructureInformation object
    print(meta_gene)
    if isinstance(meta_gene_pkl,str):
        metageneStructureInformation = load_pickle(meta_gene_pkl)
    else:
        metageneStructureInformation = meta_gene_pkl
    Info_multigenes = metageneStructureInformation[meta_gene]#list: length is gene number
    Info_multigenes_sort = []
    for info in Info_multigenes:
        isoformInfo = info[2]
        isoformInfo = dict(sorted(isoformInfo.items(), key=lambda item: len(item[1]), reverse=True))
        info[0]['isoformNames'] = list(isoformInfo.keys())
        info[2] = isoformInfo
        Info_multigenes_sort.append(info)
    Info_multigenes = Info_multigenes_sort
    if os.path.isfile(bamFile)==False: #bamFile is a folder
        bamFile_name = [f for f in os.listdir(bamFile) if f.endswith('.bam') and '.'+Info_multigenes[0][0]['geneChr']+'.' in f]
        bamFile = os.path.join(bamFile,bamFile_name[0])
    bamFilePysam = pysam.Samfile(bamFile, "rb")
    ##########-----only contain one gene for metagene
    if len(Info_multigenes)==1:
        geneInfo, exonInfo, isoformInfo = Info_multigenes[0]
        n_isoforms = len(isoformInfo)
        reads = bamFilePysam.fetch(geneInfo['geneChr'], geneInfo['geneStart'], geneInfo['geneEnd'])
        Read_Isoform_compatibleVector, Read_novelIsoform, novel_isoformInfo = [], [], {}
        samples_novel, samples_compatible=[],[]
        for read in reads:
            result = process_read(read, geneInfo, exonInfo, isoformInfo, qname_dict, lowest_match, Info_multigenes, True)
            novelIsoformResults, isoformCompatibleVectorResults = result
            if novelIsoformResults is not None:
                Read_novelIsoform.append(novelIsoformResults)
                samples_novel.append(qname_sample_dict[read.qname])
            if isoformCompatibleVectorResults is not None:
                Read_Isoform_compatibleVector.append(isoformCompatibleVectorResults)
                samples_compatible.append(qname_sample_dict[read.qname])
        unique_sample = list(set(samples_novel+samples_compatible))
        for sample in unique_sample:
            Read_novelIsoform_sample, Read_Isoform_compatibleVector_sample=[],[]
            output_folder_sample = os.path.join(output_folder, sample)
            sample_index_novel = [i for i, s in enumerate(samples_novel) if s==sample]
            sample_index_compatible = [i for i, s in enumerate(samples_compatible) if s == sample]
            if len(sample_index_novel)>0:
                Read_novelIsoform_sample = [Read_novelIsoform[i] for i in sample_index_novel]
            if len(sample_index_compatible)>0:
                Read_Isoform_compatibleVector_sample = [Read_Isoform_compatibleVector[i] for i in sample_index_compatible]
            if len(Read_novelIsoform_sample) > 0:
                Read_novelIsoform_sample, novel_isoformInfo, Read_Isoform_compatibleVector_sample = polish_compatible_vectors(Read_novelIsoform_sample,
                                                                                                                Read_Isoform_compatibleVector_sample, n_isoforms)
            geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector_sample = compile_compatible_vectors(Read_novelIsoform_sample, novel_isoformInfo,
                                                                                                                      Read_Isoform_compatibleVector_sample, geneInfo)
            save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector_sample, qname_cbumi_dict, output_folder_sample)
    ############-----multiple genes for one metagene
    else:
        geneChr, start, end = summarise_metagene(Info_multigenes)  # meta gene location
        reads = bamFilePysam.fetch(geneChr, start, end)  # fetch reads within meta gene region
        #process reads metagene
        results, samples = [], []
        for read in reads:
            out = process_read_metagene(read, start, end, qname_dict, Info_multigenes, lowest_match, True)
            if out is not None:
                results.append(out)
                samples.append(qname_sample_dict[read.qname])
        unique_sample = list(set(samples))
        for sample in unique_sample:
            output_folder_sample = os.path.join(output_folder, sample)
            sample_index = [i for i, s in enumerate(samples) if s == sample]
            result_sample = [results[i] for i in sample_index]
            Ind, read_novelisoform_tuple_metagene, read_isoform_compatibleVector_tuple_metagene = [], [], []
            for result in result_sample:
                if result is not None:
                    ind, novelisoform_tuple, isoformCompatibleVector_tuple = result
                    Ind.append(ind)
                    read_novelisoform_tuple_metagene.append(novelisoform_tuple)
                    read_isoform_compatibleVector_tuple_metagene.append(isoformCompatibleVector_tuple)
            unique_ind = list(set(Ind))
            log_ind = [ind for ind in range(len(Info_multigenes)) if ind not in unique_ind]
            #logging genes without any reads
            for index in log_ind:
                save_compatibleVector_by_gene(Info_multigenes[index][0]['geneName'], Info_multigenes[index][0]['geneID'], Info_multigenes[index][0]['geneStrand'], colNames=None,
                                              Read_Isoform_compatibleVector=None, qname_cbumi_dict=None,output_folder=output_folder_sample)
            for index in unique_ind:
                print('processing gene'+str(index))
                #loop over genes within metagene; for one single gene:
                Read_novelIsoform, Read_Isoform_compatibleVector, novel_isoformInfo=[],[], {}
                for j, i in enumerate(Ind):
                        #loop for reads
                    if i==index and read_novelisoform_tuple_metagene[j] is not None:
                        Read_novelIsoform.append(read_novelisoform_tuple_metagene[j])
                    if i==index and read_isoform_compatibleVector_tuple_metagene[j] is not None:
                        Read_Isoform_compatibleVector.append(read_isoform_compatibleVector_tuple_metagene[j])
                if len(Read_novelIsoform) > 0:
                    Read_novelIsoform, novel_isoformInfo, Read_Isoform_compatibleVector = polish_compatible_vectors(
                        Read_novelIsoform, Read_Isoform_compatibleVector, len(Info_multigenes[index][2]))
                #loop over reads for one gene
                geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(Read_novelIsoform,
                                                                                                   novel_isoformInfo,
                                                                                                   Read_Isoform_compatibleVector,
                                                                                                   Info_multigenes[index][0])
                save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector, qname_cbumi_dict, output_folder_sample)
    #finish mapping reads to genes
    bamFilePysam.close()


#main function
def generate_compatibleVector(bamFile, qname_dict, qname_cbumi_dict, metagene_pkl, lowest_match, output_folder, job_index=0,
                              total_jobs=1, cover_existing = True):
    if isinstance(metagene_pkl, str):
        metageneStructureInformation = load_pickle(metagene_pkl)
    else:
        metageneStructureInformation = metagene_pkl
    target = os.path.join(output_folder, 'compatible_matrix')
    if cover_existing:
        print('If there are existing compatible matrix files, SCOTCH will overwrite them')
        genes_existing = []
    else:
        print('If there are existing compatible matrix files, SCOTCH will not overwrite them')
        genes_existing = [g[:-4] for g in os.listdir(target)]
        if os.path.isfile(os.path.join(target, 'log.txt')):
            gene_df = pd.read_csv(os.path.join(target, 'log.txt'), header=None)
            genes_existing = genes_existing + gene_df.iloc[:, 0].tolist()
    MetaGene_Gene_dict = {}
    for key, values in metageneStructureInformation.items():
        genes_ = []
        for value in values:
            gene = str(value[0]['geneName'])+'_'+str(value[0]['geneID'])
            if gene not in genes_existing:
                genes_.append(gene)
        if len(genes_) > 0:
            MetaGene_Gene_dict[key] = genes_
    MetaGenes = list(MetaGene_Gene_dict.keys())
    print('total metagene number is: ' + str(len(MetaGenes)))
    if total_jobs > 1:
        step_size = math.ceil(len(MetaGenes) / total_jobs)
        s = int(list(range(0, len(MetaGenes), step_size))[job_index])
        e = int(s + step_size)
        MetaGenes_ = MetaGenes[s:e]
        print('processing: ' + str(len(MetaGenes_)) + ' metagenes')
    else:
        MetaGenes_ = MetaGenes
    for meta_gene in MetaGenes_:
        extract_compatible_matrix_by_metagene(bamFile, metagene_pkl, meta_gene, qname_dict,qname_cbumi_dict, lowest_match, output_folder)



def generate_compatibleVector_pacbio(bamFile, qname_dict, qname_cbumi_dict, metagene_pkl, lowest_match, output_folder, job_index=0,
                              total_jobs=1, cover_existing = True):
    if isinstance(metagene_pkl, str):
        metageneStructureInformation = load_pickle(metagene_pkl)
    else:
        metageneStructureInformation = metagene_pkl
    target = os.path.join(output_folder, 'compatible_matrix')
    if cover_existing:
        print('If there are existing compatible matrix files, SCOTCH will overwrite them')
        genes_existing = []
    else:
        print('If there are existing compatible matrix files, SCOTCH will not overwrite them')
        genes_existing = [g[:-4] for g in os.listdir(target)]
        if os.path.isfile(os.path.join(target, 'log.txt')):
            gene_df = pd.read_csv(os.path.join(target, 'log.txt'), header=None)
            genes_existing = genes_existing + gene_df.iloc[:, 0].tolist()
    MetaGene_Gene_dict = {}
    for key, values in metageneStructureInformation.items():
        genes_ = []
        for value in values:
            gene = str(value[0]['geneName'])+'_'+str(value[0]['geneID'])
            if gene not in genes_existing:
                genes_.append(gene)
        if len(genes_) > 0:
            MetaGene_Gene_dict[key] = genes_
    MetaGenes = list(MetaGene_Gene_dict.keys())
    print('total metagene number is: ' + str(len(MetaGenes)))
    if total_jobs > 1:
        step_size = math.ceil(len(MetaGenes) / total_jobs)
        s = int(list(range(0, len(MetaGenes), step_size))[job_index])
        e = int(s + step_size)
        MetaGenes_ = MetaGenes[s:e]
        print('processing: ' + str(len(MetaGenes_)) + ' metagenes')
    else:
        MetaGenes_ = MetaGenes
    for meta_gene in MetaGenes_:
        extract_compatible_matrix_by_metagene_pacbio(bamFile, metagene_pkl, meta_gene, qname_dict,qname_cbumi_dict, lowest_match, output_folder)


def generate_compatibleVector_parse(bamFile, qname_dict, qname_cbumi_dict, qname_sample_dict, metagene_pkl, lowest_match, output_folder, job_index=0,
                              total_jobs=1, cover_existing = True):
    if isinstance(metagene_pkl, str):
        metageneStructureInformation = load_pickle(metagene_pkl)
    else:
        metageneStructureInformation = metagene_pkl
    target = os.path.join(output_folder, 'compatible_matrix')
    if cover_existing:
        print('If there are existing compatible matrix files, SCOTCH will overwrite them')
        genes_existing = []
    else:
        print('If there are existing compatible matrix files, SCOTCH will not overwrite them')
        genes_existing = [g[:-4] for g in os.listdir(target)]
        if os.path.isfile(os.path.join(target, 'log.txt')):
            gene_df = pd.read_csv(os.path.join(target, 'log.txt'), header=None)
            genes_existing = genes_existing + gene_df.iloc[:, 0].tolist()
    MetaGene_Gene_dict = {}
    for key, values in metageneStructureInformation.items():
        genes_ = []
        for value in values:
            gene = str(value[0]['geneName'])+'_'+str(value[0]['geneID'])
            if gene not in genes_existing:
                genes_.append(gene)
        if len(genes_) > 0:
            MetaGene_Gene_dict[key] = genes_
    MetaGenes = list(MetaGene_Gene_dict.keys())
    print('total metagene number is: ' + str(len(MetaGenes)))
    if total_jobs > 1:
        step_size = math.ceil(len(MetaGenes) / total_jobs)
        s = int(list(range(0, len(MetaGenes), step_size))[job_index])
        e = int(s + step_size)
        MetaGenes_ = MetaGenes[s:e]
        print('processing: ' + str(len(MetaGenes_)) + ' metagenes')
    else:
        MetaGenes_ = MetaGenes
    for meta_gene in MetaGenes_:
        extract_compatible_matrix_by_metagene_parse(bamFile, metagene_pkl, meta_gene, qname_dict, qname_cbumi_dict,qname_sample_dict,
                                                  lowest_match, output_folder)


def detect_poly_parse(read, window = 20, n = 10):
    query_sequence = read.query_sequence
    poly_bool = False
    seq = Seq(query_sequence)
    nT = max([seq[:window].count('T'),seq[-window].count('T')])
    nA = max([seq[:window].count('A'),seq[-window].count('A')])
    poly=None
    if nT>=n:
        poly = 'T'
    if nA >= n:
        poly = 'A'
    if nT>=n or nA>=n:
        poly_bool = True
    return poly_bool, poly

##TODO: this function is for 3' kit, should consider 5'kit
def detect_poly(read, window = 20, n = 15):
    if isinstance(read, dict):
        query_sequence = read['query_sequence']
        umi = read['umi']
    else:
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


def exon_hit(mapPositions,exonInfo):
     mapExons=[]
     exon_index = 0
     for i in mapPositions:
         while i >= exonInfo[exon_index][1] and exon_index<len(exonInfo)-1:
             exon_index += 1
         if int(i)>=int(exonInfo[exon_index][0]) and int(i) < int(exonInfo[exon_index][1]):
             mapExons.append(exon_index)
     return mapExons

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


############group novel isoforms#############
def group_novel_isoform(df, geneStrand):
    df_novel = df.filter(like='novelIsoform')
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
            isoform_group_name = ['novelIsoform_' + str(i) for i in novel_isoform_group_list[i]]
            df_novel_ = df_novel[isoform_group_name].max(axis=1)
            df_novel_group.append(df_novel_)
        df_novel_grouped = pd.concat(df_novel_group, axis=1)
        df_novel_grouped.columns = novel_isoform_group_name
        df = pd.concat([df_existing, df_novel_grouped, df_uncategorized], axis=1)
    return df


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

#bamFile='/scr1/users/xu3/singlecell/project_singlecell/LH/bam/LH.filtered.bam'
#gene_pkl='/scr1/users/xu3/singlecell/project_singlecell/LH/reference/geneStructureInformation.pkl'



