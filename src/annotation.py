from collections import defaultdict
import pysam
from joblib import Parallel, delayed
import reference as ref
import os
import pandas as pd
import re
import pickle
from preprocessing import load_pickle, merge_exons
from scipy.ndimage import gaussian_filter1d
import numpy as np
import shutil


######################################################################
##############################annotation##############################
######################################################################

def extract_bam_info_folder(bam_folder, num_cores, parse=False, pacbio = False):
    files = os.listdir(bam_folder)
    bamfiles = [os.path.join(bam_folder,f) for f in files if f.endswith('.bam')]
    if parse:
        bamfiles = [bam for bam in bamfiles if bam+'.bai' in files]
        df = Parallel(n_jobs=num_cores)(delayed(extract_bam_info_parse)(bam) for bam in bamfiles)
    elif pacbio:
        df = Parallel(n_jobs=num_cores)(delayed(extract_bam_info_pacbio)(bam) for bam in bamfiles)
    else:
        df = Parallel(n_jobs=num_cores)(delayed(extract_bam_info)(bam) for bam in bamfiles)
    ReadTagsDF = pd.concat(df).reset_index(drop=True)
    return ReadTagsDF

def extract_bam_info(bam, barcode_cell = 'CB', barcode_umi = 'UB'):
    #extract readname, cb, umi from bam file
    # bam: path to bam file
    bamFilePysam = pysam.Samfile(bam, "rb")
    ReadTags = [(read.qname, read.get_tag(barcode_cell), read.get_tag(barcode_umi), read.qend-read.qstart) for read in bamFilePysam]
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
    match = re.search(r'sublibrary(\d+)', bam)
    if match:
        sublib = match.group(1)
    else:
        sublib = '1'
    #qname cb umi cbumi length
    ReadTags = [(read.qname, read.qname.split('_')[-5]+'_'+read.qname.split('_')[-4]+'_'+read.qname.split('_')[-3],
                 read.qname.split('_')[-1] , len(read.query_alignment_sequence),
                 read.get_tag('pS'), sublib) for read in bamFilePysam]
    ReadTagsDF = pd.DataFrame(ReadTags)
    if ReadTagsDF.shape[0] > 0:
        ReadTagsDF.columns = ['QNAME', 'CB', 'UMI', 'LENGTH','SAMPLE', 'SUBLIB']
        ReadTagsDF = ReadTagsDF.sort_values(by=['SAMPLE','CB', 'UMI', 'LENGTH'], ascending=[True, True, True, False]).reset_index(drop=True)
        ReadTagsDF['CBUMI'] = ReadTagsDF.CB.astype(str) + '_' + ReadTagsDF.UMI.astype(str) + '_' + ReadTagsDF.SUBLIB.astype(str)
    else:
        ReadTagsDF = None
    return ReadTagsDF

def extract_bam_info_pacbio(bam, barcode_cell = 'XC', barcode_umi = 'XM'):
    bamFilePysam = pysam.Samfile(bam, "rb")
    #qname cb umi cbumi length
    ReadTags = [(read.qname, read.get_tag(barcode_cell), read.get_tag(barcode_umi), read.reference_length) for read in bamFilePysam]
    ReadTagsDF = pd.DataFrame(ReadTags)
    if ReadTagsDF.shape[0] > 0:
        ReadTagsDF.columns = ['QNAME', 'CB', 'UMI', 'LENGTH']
        ReadTagsDF.dropna(inplace=True)
        ReadTagsDF = ReadTagsDF.sort_values(by=['CB','UMI','LENGTH'],ascending=[True, True, False]).reset_index(drop=True)
        ReadTagsDF['CBUMI'] = ReadTagsDF.CB.astype(str) + '_' + ReadTagsDF.UMI.astype(str)
        ReadTagsDF['QNAME'] = ReadTagsDF.QNAME.astype(str) + '_' + ReadTagsDF.LENGTH.astype(int).astype(str)
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
    qname_sample_dict = None
    if parse:
        qname_sample_dict = dict(zip(bam_info['QNAME'], bam_info['SAMPLE']))
    return qname_dict, qname_cbumi_dict, qname_sample_dict

def process_gene(geneID, geneName, genes ,exons, build=None):
    GeneDf = genes[genes.iloc[:, 3] == geneID].reset_index(drop=True)
    ExonsDf = exons[exons.iloc[:, 3] == geneID].reset_index(drop=True)
    exonsdf = ExonsDf.iloc[:, :3].drop_duplicates().reset_index(drop=True)
    exonsdf['exon_num'] = list(range(exonsdf.shape[0]))
    geneChr, geneStart, geneEnd = GeneDf.iloc[0, 0], GeneDf.iloc[0, 1], GeneDf.iloc[0, 2]
    if build is not None:
        geneChr = str(build) + '_'+str(geneChr)
    geneStrand = GeneDf.iloc[0, 6]
    isoformNames = ExonsDf.TRANSCRIPT.unique().tolist()
    GeneInfo = {'geneName': geneName, 'geneID': geneID, 'geneChr': geneChr, 'geneStart': geneStart, 'geneEnd': geneEnd,
                'geneStrand':geneStrand, 'numofExons': exonsdf.shape[0], 'numofIsoforms': len(isoformNames),
                'isoformNames': isoformNames}
    ExonPositions = list(zip(exonsdf.iloc[:, 1], exonsdf.iloc[:, 2]))
    ExonsDf = ExonsDf.merge(exonsdf, how='left')
    ExonIsoformDict = dict()
    for isoform in isoformNames:
        ExonIsoformDict[isoform] = ExonsDf[ExonsDf.TRANSCRIPT == isoform].exon_num.tolist()
    return geneID, [GeneInfo, ExonPositions, ExonIsoformDict]


######################################################################

def get_continuous_blocks(data):
    blocks = []
    if len(data) == 0:
        return blocks
    if isinstance(data, dict):
        sorted_positions = sorted(data.keys())
    elif isinstance(data, list):
        sorted_positions = sorted(data)
    else:
        raise TypeError("Input must be a dictionary or a list.")
    current_start = sorted_positions[0]
    # for a single base, pysam uses 0-based
    # we will use the same 0-based half open system as in bed file: e.g. (0,3): 1st,2nd,3rd bases
    current_end = sorted_positions[0] + 1
    for pos in sorted_positions[1:]:
        if pos == current_end:
            current_end = pos + 1
        else:
            blocks.append((current_start, current_end))
            current_start = pos
            current_end = pos + 1
    blocks.append((current_start, current_end))
    return blocks

def count_and_sort_tuple_elements(tuple_list):
    freq_dict = defaultdict(int)
    for a, b in tuple_list:
        freq_dict[a] += 1
        freq_dict[b] += 1
    sorted_freq_dict = dict(sorted(freq_dict.items()))
    return sorted_freq_dict


def merge_boundaries_by_evidence(boundaries, merge_distance=10):
    # Sort boundaries by value (evidence) in descending order
    sorted_boundaries = sorted(boundaries.items(), key=lambda item: item[1], reverse=True)
    merged_boundaries = {}
    while sorted_boundaries:
        # Take the boundary with the highest evidence
        base_position, base_value = sorted_boundaries.pop(0)
        merged_value = base_value
        to_remove = []
        # Check nearby positions within the merge_distance
        for i, (position, value) in enumerate(sorted_boundaries):
            if abs(position - base_position) <= merge_distance:
                merged_value += value
                to_remove.append(i)
        # Remove merged positions from sorted_boundaries
        for index in sorted(to_remove, reverse=True):
            sorted_boundaries.pop(index)
        # Add the merged boundary to the result
        merged_boundaries[base_position] = merged_value
    return merged_boundaries



def get_splicejuction_from_annotation(exonInfo, isoformInfo):
    isoform_names = list(isoformInfo.keys())
    sjInfo = {}
    for i in range(len(isoform_names)):
        exon_indexs = isoformInfo[isoform_names[i]]
        exons = [exonInfo[ind] for ind in exon_indexs]
        exons = merge_exons(exons)
        splice_junctions = []
        if len(exons)>1:
            a, b = exons[0]
            for s, e in exons[1:]:
                splice_junctions.append((b, s))
                a, b = s, e
        sjInfo[isoform_names[i]] = splice_junctions
    return sjInfo


def get_splicejuction_from_read(read):
    junctions = []
    ref_pos = read.reference_start
    for cigartype, cigarlen in read.cigartuples:
        if cigartype == 3:  # 'N' operation indicates a splice junction
            junction_start = ref_pos
            junction_end = ref_pos + cigarlen
            junctions.append((junction_start, junction_end))
        if cigartype in (0, 2, 3):
            ref_pos += cigarlen
    return junctions #this returns introns e.g. (103, 105) is 0 based containing 2 bases

def cut_exons_by_derivative(exons, coverage, sigma = 1, z_score_threshold = 3):
    #exon_dict: just one exon
    #coverage: the whole coverage object {1000:100, 1001:101}
    derivatives, positions = [], []
    for exon in exons:
        exon_start, exon_end = exon
        for pos in range(exon_start, exon_end):
            if pos not in coverage:
                coverage[pos] = 0
        coverage_values = np.array([coverage[i] for i in range(exon_start, exon_end)])
        smoothed_coverage = gaussian_filter1d(coverage_values, sigma)
        derivative = np.diff(smoothed_coverage)
        derivatives.append(derivative.tolist())
        positions.extend(range(exon_start + 1, exon_end))
    derivatives = np.asarray([ii for i in derivatives for ii in i])
    positions = np.asarray(positions)
    mean_derivative = np.mean(derivatives)
    std_derivative = np.std(derivatives)
    if std_derivative==0:
        cut_positions = []
    else:
        z_scores = (derivatives - mean_derivative) / std_derivative
        cut_positions = positions[np.abs(z_scores) > z_score_threshold]
    if len(cut_positions)==0:
        return exons
    cut_positions_valid = [cut_positions[0]]
    if len(cut_positions)>1:
        for i in range(1, len(cut_positions)):
            if cut_positions[i] - cut_positions[i - 1] > 10:
                cut_positions_valid.append(cut_positions[i])
    new_exons = []
    for exon in exons:
        exon_start, exon_end = exon
        cuts = [pos for pos in cut_positions_valid if pos - exon_start > 10 and exon_end-pos > 10]
        if len(cuts)==0:
            new_exons.append(exon)
        else:
            exon_positions = [exon_start]+cuts +[exon_end]
            for i in range(len(exon_positions)-1):
                new_exons.append((exon_positions[i], exon_positions[i+1]))
    return new_exons


def get_non_overlapping_exons(bam_file, chrom, gene_start, gene_end,
                              coverage_threshold=0.02, boundary_threshold=0.02,
                              z_score_threshold = 10):
    """
    get non overlapping exons based on bam files given a gene region
    :param bam_file: a single bam file path/oject or a list of multiple bam file paths/objects that will be used for concensus result
    :param chrom: chromosome
    :param gene_start: start position of gene
    :param gene_end: end position of gene
    :param coverage_threshold: threshold of coverage to filter out exons, percentage to maximum coverage
    :param boundary_threshold: threshold of read number to filter out splicing positions, percentage to maximum read supports of splicing positions
    :param z_score_threshold: threshold to filter out within exon partitions based on sharp changes of derivatives
    :return: a list of tuples representing coordinates of exons
    """
    if isinstance(bam_file,str):
        # a single path
        bams = [pysam.AlignmentFile(bam_file, "rb")]
    elif isinstance(bam_file,list):
        # a list of path
        if isinstance(bam_file[0],str):
            bams = [pysam.AlignmentFile(bam_, "rb") for bam_ in bam_file]
        else:
            # a list of objects
            bams = bam_file
    else:
        # a single object
        bams = [bam_file] #a pysam object
    #--> bams is a list of pysam objects
    exons, coverage, junctions = [], {}, []
    #get read coverage
    for bam in bams:
        for read in bam.fetch(chrom, gene_start, gene_end):
            if read.reference_start >= gene_start and read.reference_end < gene_end:
                junctions += get_splicejuction_from_read(read)
                for read_start, read_end in read.get_blocks():
                    for pos in range(read_start, read_end):
                        if pos in coverage:
                            coverage[pos] += 1
                        else:
                            coverage[pos] = 1
    if len(coverage)==0:
        return exons
    coverage_threshold_absolute = max(max(coverage.values())*coverage_threshold,20)
    #only keep positions that meet the coverage threshold
    coverage = {pos: cov for pos, cov in coverage.items() if cov > coverage_threshold_absolute}
    coverage_blocks = get_continuous_blocks(coverage)
    coverage_blocks = [(a,b) for a, b in coverage_blocks if b-a >= 20]#delete less than 20bp exons
    if len(coverage_blocks)==0:
        return exons
    # fill up holes less than 20bp
    if len(coverage_blocks) >1:
        coverage_blocks_ = []
        current_start, current_end = coverage_blocks[0]
        for start, end in coverage_blocks[1:]:
            if start - current_end > 20:
                coverage_blocks_.append((current_start, current_end))
                current_start, current_end = start, end
            else:
                current_end = end
        coverage_blocks_.append((current_start, current_end))
        coverage_blocks = coverage_blocks_
    #update exons to meta-exons
    exons = coverage_blocks
    #partition1: use splicing junctions to find sub-exons
    boundaries = count_and_sort_tuple_elements(junctions)
    #filter out splicing positions with few support
    if len(boundaries)>1:
        boundary_threshold_absolute = max(max(list(boundaries.values())) * boundary_threshold, 20)
        #group boundaries by meta-exons
        filtered_boundaries = [{} for _ in coverage_blocks]
        for boundary, freq in boundaries.items():
            if freq >= boundary_threshold_absolute:
                for i, (start, end) in enumerate(coverage_blocks):
                    if start < boundary < end:
                        filtered_boundaries[i][boundary] = freq
                        break
        Filtered_Boundaries= []
        for filtered_boundaries_ in filtered_boundaries:
            merged_boundaries = merge_boundaries_by_evidence(filtered_boundaries_, merge_distance=10)
            sorted_boundaries = dict(sorted(merged_boundaries.items()))
            Filtered_Boundaries.append(sorted_boundaries)
        #futher filter the boundaries that are too close to a boundary of an exon_block
        for i, block_boundaries in enumerate(Filtered_Boundaries):
            for boundary in list(block_boundaries.keys()):
                if abs(boundary - coverage_blocks[i][0])<=10 or abs(boundary - coverage_blocks[i][1])<=10:
                    del Filtered_Boundaries[i][boundary]
        #partition exon blocks into sub exons based on splicing positions
        exons = [] #----re-initiate exons
        for i in range(len(coverage_blocks)):
            start, end = coverage_blocks[i]
            positions = [start]+list(Filtered_Boundaries[i].keys())+[end]
            for ii in range(len(positions)-1):
                exons.append((positions[ii], positions[ii+1]))
    #partition2: partition meta-exons using read coverage derivatives
    exons = cut_exons_by_derivative(exons, coverage, sigma=1, z_score_threshold=z_score_threshold)
    return exons


def get_genes_from_bam(input_bam_path, coverage_threshold = 5, min_region_size=50):
    def process_bam_file(bam_file, coverage_threshold, min_region_size):
        #bam_file: a single path or a list of paths for consensus results
        if isinstance(bam_file, str):
            bams = pysam.AlignmentFile(bam_file, "rb")
        else: #list
            bams = [ pysam.AlignmentFile(bam, "rb") for bam in bam_file]
        coverage = defaultdict(lambda: defaultdict(int))
        chromosomes = list(set([bam.references for bam in bams]))
        for chrom in chromosomes:
            for bam in bams:
                if chrom in bam.references:
                    for read in bam.fetch(chrom):
                        if not read.is_unmapped:
                            for pos in range(read.reference_start, read.reference_end):
                                coverage[chrom][pos] += 1
        genes = {}
        for chrom, cov_dict in coverage.items():
            cov_dict = {pos: cov for pos, cov in cov_dict.items() if cov > coverage_threshold}
            coverage_blocks = get_continuous_blocks(cov_dict)
            genes[chrom] = [(a, b) for a, b in coverage_blocks if b - a > min_region_size]
        return genes
    if isinstance(input_bam_path, list):
        if os.path.isdir(input_bam_path[0]): # a list of folders
            bam_files = [os.path.join(folder, file) for folder in input_bam_path for file in folder]
            #bam_files = [os.path.join(input_bam_path, f) for f in os.listdir(input_bam_path) if f.endswith('.bam')]
        else: # a list of paths
            bam_files = input_bam_path
    else:#single sample
        if os.path.isdir(input_bam_path):
            bam_files = [os.path.join(input_bam_path, f) for f in os.listdir(input_bam_path) if f.endswith('.bam')]
        else:
            bam_files = [input_bam_path]
    all_genes = process_bam_file(bam_files, coverage_threshold, min_region_size)
    sorted_all_genes = {}
    for chrom in sorted(all_genes.keys()):
        sorted_all_genes[chrom] = sorted(all_genes[chrom], key=lambda x: (x[0], x[1]))
    return sorted_all_genes



def update_exons(A, B, distance_threshold = 20):
    ##B is reference
    ##A is based on bam
    def correct_point(point, reference_points, threshold = distance_threshold):
        closest_point = None
        min_distance = float('inf')
        for ref_point in reference_points:
            distance = abs(point - ref_point)
            if distance < threshold and distance < min_distance:
                closest_point = ref_point
                min_distance = distance
        return closest_point if closest_point is not None else point
    def partition_intervals(A, B):
        all_intervals = A + B
        all_intervals.sort()
        partitions = []
        current_start, current_end = all_intervals[0]
        for start, end in all_intervals[1:]:
            if start >= current_end:
                if current_start < current_end:
                    partitions.append((current_start, current_end))
                current_start, current_end = start, end
            else:
                temp = [current_start, current_end, start, end]
                temp.sort()
                for i in [0, 1]:
                    if temp[i]<temp[i + 1]:
                        partitions.append((temp[i], temp[i + 1]))
                current_start, current_end = temp[2], temp[3]
        if current_start<current_end:
            partitions.append((current_start, current_end))
        return partitions
    if len(A)==0:
        return B
    if len(B)==1:
        return B
    reference_points = [point for segment in B for point in segment]
    corrected_A = []
    for start, end in A:
        corrected_start = correct_point(start, reference_points)
        corrected_end = correct_point(end, reference_points)
        corrected_A.append((corrected_start, corrected_end))
    corrected_A = [(a,b) for a,b in corrected_A if a<b]
    return partition_intervals(corrected_A,B)




#modes involving bam file
def annotate_genes(geneStructureInformation, bamfile_path,
                   coverage_threshold_gene=20, coverage_threshold_exon=0.02,coverage_threshold_splicing=0.02,z_score_threshold = 10,
                   min_gene_size=50, workers=1):
    """
    generate geneStructureInformation either using bamfile alone (leave geneStructureInformation blank) or update existing annotation file using bam file
    :param geneStructureInformation: pickle file object of existing gene annotation, not path
    :param bamfile_path: can be (1) a single bam file from one sample, (2) a folder path of multiple bams from one sample, (3) a list of bam files from multiple samples, (4) a list of bam folder paths from multiple samples
    :param coverage_threshold_gene: minimal read coverage for gene
    :param coverage_threshold_exon: minimal read coverage for exon, percentage to the maximum coverage
    :param coverage_threshold_splicing: minimal read coverage to support splicing, percentage to the maximum splicing
    :param z_score_threshold: threshold to identify sharpe changes of coverage derivative
    :param min_gene_size: minimal length for novel gene discovery
    :return: updated geneStructureInformation annotation
    """
    def novel_gene_annotation(chrom, all_genes_dict, bamfile, coverage_threshold_exon=0.01, coverage_threshold_splicing=0.01,z_score_threshold = 10):
        #all_genes: {'chr1':[(100,200),(300,400)]}, bamfile: path to one bam file
        bamfiles = None
        #----from one sample
        if isinstance(bamfile, str):
            if os.path.isfile(bamfile) == False:  # bamfile is a folder path
                bamFile_name = [f for f in os.listdir(bamfile) if
                                f.endswith('.bam') and '.' + str(chrom) + '.' in f]
                bamfiles = [pysam.AlignmentFile(os.path.join(bamfile, bamFile_name[0]))]# bamfile is a file path now
            else:
                bamfiles = [pysam.AlignmentFile(bamfile, "rb")] # bamfile is a file path
        #-----from multiple samples
        elif isinstance(bamfile, list):
            if os.path.isfile(bamfile[0]) == False:  # bamfile = [folderpath1, folderpath2, ...]
                bamFile_name = [] #bamFile_name = [folderpath1/bamname1, folderpath2/bamname2, ...]
                for folder in bamfile:
                    for file in os.listdir(folder):
                        if file.endswith('.bam') and '.' + str(chrom) + '.' in file:
                            bamFile_name.append(os.path.join(folder,file))
                bamfiles = [pysam.AlignmentFile(bfn, "rb") for bfn in bamFile_name]
            else:
                bamfiles = [pysam.AlignmentFile(bfn, "rb") for bfn in bamfile]
        else:
            print('bamfile must be a list or str')
        #bamfiles is a list of pysam ojects
        genes = all_genes_dict[chrom]  # [(100,200),(300,400)]
        gene_annotations = {}
        for i, (gene_start, gene_end) in enumerate(genes):
            exonInfo = get_non_overlapping_exons(bamfiles, chrom, gene_start, gene_end, coverage_threshold_exon, coverage_threshold_splicing,
                                                 z_score_threshold)
            geneInfo = {'geneName': 'gene_'+ str(i+1)+'_'+chrom, 'geneID': 'gene_'+ str(i+1)+'_'+chrom,
                        'geneChr': chrom, 'geneStart':gene_start, 'geneEnd':gene_end, 'geneStrand': '.',
                        'numofExons': len(exonInfo), 'numofIsoforms': 0, 'isoformNames':[]}
            isoformInfo = {}
            gene_annotations['gene_'+ str(i+1)+'_'+chrom] = [geneInfo, exonInfo, isoformInfo]
        return gene_annotations
    def update_isoform_info(original_exons, updated_exons, isoformInfo):
        exon_mapping = {}
        for original_idx, (original_start, original_end) in enumerate(original_exons):
            for updated_idx, (updated_start, updated_end) in enumerate(updated_exons):
                if original_start <= updated_start and original_end >= updated_end:
                    if original_idx not in exon_mapping:
                        exon_mapping[original_idx] = []
                    exon_mapping[original_idx].append(updated_idx)
        updated_isoform_info = {}
        for isoform, exon_indices in isoformInfo.items():
            updated_indices = [exon_mapping[idx] for idx in exon_indices if idx in exon_mapping]
            updated_indices = [ii for i in updated_indices for ii in i]
            updated_isoform_info[isoform] = updated_indices
        return updated_isoform_info
    def update_annotation(geneStructureInformation, geneID, bamfile_path,coverage_threshold_exon, coverage_threshold_splicing, z_score_threshold):
        chrom, gene_start, gene_end = geneStructureInformation[geneID][0]['geneChr'], \
        geneStructureInformation[geneID][0]['geneStart'], geneStructureInformation[geneID][0]['geneEnd']
        bamfiles = None
        # ----from one sample
        if isinstance(bamfile_path, str):
            if os.path.isfile(bamfile_path) == False:  # bamfile is a folder path
                bamFile_name = [f for f in os.listdir(bamfile_path) if
                                f.endswith('.bam') and '.' + str(chrom) + '.' in f]
                bamfiles = [pysam.AlignmentFile(os.path.join(bamfile_path, bamFile_name[0]))]  # bamfile is a file path now
            else:
                bamfiles = [pysam.AlignmentFile(bamfile_path, "rb")]  # bamfile is a file path
        # -----from multiple samples
        elif isinstance(bamfile_path, list):
            if os.path.isfile(bamfile_path[0]) == False:  # bamfile = [folderpath1, folderpath2, ...]
                bamFile_name = []  # bamFile_name = [folderpath1/bamname1, folderpath2/bamname2, ...]
                for folder in bamfile_path:
                    for file in os.listdir(folder):
                        if file.endswith('.bam') and '.' + str(chrom) + '.' in file:
                            bamFile_name.append(os.path.join(folder, file))
                bamfiles = [pysam.AlignmentFile(bfn, "rb") for bfn in bamFile_name]
            else:
                bamfiles = [pysam.AlignmentFile(bfn, "rb") for bfn in bamfile_path]
        else:
            print('bamfile must be a list or str')
        exons_bam = get_non_overlapping_exons(bamfiles, chrom, gene_start, gene_end, coverage_threshold_exon, coverage_threshold_splicing, z_score_threshold)
        original_exons = geneStructureInformation[geneID][1]
        updated_exons = update_exons(exons_bam, original_exons)
        geneStructureInformation_copy = geneStructureInformation.copy()
        # geneInfo
        geneInfo = geneStructureInformation_copy[geneID][0]
        geneInfo['numofExons'] = len(updated_exons)
        # exonInfo
        exonInfo = updated_exons
        # isoformInfo
        isoformInfo = update_isoform_info(original_exons, updated_exons, geneStructureInformation_copy[geneID][2])
        return {geneID:[geneInfo, exonInfo, isoformInfo]}
    #generate gene annotation solely based on bam file
    if geneStructureInformation is None:
        all_genes = get_genes_from_bam(bamfile_path, coverage_threshold_gene, min_gene_size) #{'chr1':[(100,200),(300,400)]}
        chroms = list(all_genes.keys())
        results = Parallel(n_jobs=workers)(
            delayed(novel_gene_annotation)(chrom, all_genes, bamfile_path, coverage_threshold_exon, coverage_threshold_splicing, z_score_threshold) for chrom in chroms)
        annotations = {k: v for result in results for k, v in result.items()}
    #update existing gene annotation using bam file
    else:
        geneIDs = list(geneStructureInformation.keys())
        #re-order geneID
        #chunks = np.array_split(geneIDs, workers)
        #geneIDs = [item for sublist in zip(*chunks) for item in sublist]
        results = Parallel(n_jobs=workers)(
            delayed(update_annotation)(geneStructureInformation, geneID, bamfile_path,coverage_threshold_exon, coverage_threshold_splicing, z_score_threshold) for geneID in geneIDs)
        annotations = {k: v for result in results for k, v in result.items()}
    return annotations #{geneID:[geneInfo, exonInfo, isoformInfo]}


def extract_annotation_info(refGeneFile_gtf_path, refGeneFile_pkl_path, bamfile_path, num_cores=8,
                            output="geneStructureInformation.pkl", build=None,
                            coverage_threshold_gene=5, coverage_threshold_exon=0.01,coverage_threshold_splicing=0.01, z_score_threshold = 10,
                            min_gene_size=50):
    """
    wrapper function to extract gene annotation information including exons and isoforms
    :param refGeneFile_path: path to gene annotation gtf file, or pickle file that can be directly used
    :param bamfile_path: path to bam file or the folder for bam files
    :param num_cores: number of workers for parallel computing
    :param output: default: geneStructureInformation.pkl
    :param build: used for parse platform naming
    :return: metageneStructureInformation
    """
    geneStructureInformation = None
    meta_output = os.path.join(os.path.dirname(output), 'meta' + os.path.basename(output))
    if refGeneFile_gtf_path is not None:
        genes, exons = ref.generate_reference_df(gtf_path=refGeneFile_gtf_path)
    else:
        genes = None
    #####################################################
    #option1: ---------rely on bam file alone---------###
    #####################################################
    if refGeneFile_gtf_path is None and refGeneFile_pkl_path is None and bamfile_path is not None:
        print('rely on bam file alone to generate gene annotations')
        geneStructureInformation = annotate_genes(geneStructureInformation = None,
                                                  bamfile_path = bamfile_path,
                                                  coverage_threshold_gene=coverage_threshold_gene,
                                                  coverage_threshold_exon=coverage_threshold_exon,
                                                  coverage_threshold_splicing = coverage_threshold_splicing,
                                                  z_score_threshold = z_score_threshold,
                                                  min_gene_size=min_gene_size,
                                                  workers=num_cores)#no need to add build
        if output is not None:
            with open(output, 'wb') as file:
                pickle.dump(geneStructureInformation, file)
    #####################################################
    # option2: --------rely on existing annotation alone#
    #####################################################
    if refGeneFile_pkl_path is None and refGeneFile_gtf_path is not None: ###use gtf
        print('use the existing gtf file for gene annotations')
        Genes = list(zip(genes.iloc[:, 3].tolist(), genes.iloc[:, 4].tolist()))  # id, name
        #generate single gene annotations if not existing
        if os.path.isfile(output) == False:
            geneStructureInformation = Parallel(n_jobs=num_cores)(delayed(process_gene)(geneID, geneName, genes, exons, build) for geneID, geneName in Genes)
            geneStructureInformation = dict(geneStructureInformation)
            geneStructureInformation = add_build(geneStructureInformation, build)
            print('finish generating geneStructureInformation.pkl')
            #save to output, single gene
            if output is not None:
                with open(output, 'wb') as file:
                    pickle.dump(geneStructureInformation, file)
        else:#there exist pre-computate annotation file
            print('load existing annotation pickle file of each single gene at: '+str(output))
            geneStructureInformation = load_pickle(output)
    if refGeneFile_pkl_path is not None: ###use pickle
        assert refGeneFile_gtf_path is not None, 'gtf reference file is still needed! please input one'
        print('load existing annotation pickle file of each single gene at: ' + str(refGeneFile_gtf_path))
        geneStructureInformation = load_pickle(refGeneFile_pkl_path)
        #check if geneStructureInformation contains build
        if build is not None:
            if not geneStructureInformation[list(geneStructureInformation.keys())[0]][0]['geneChr'].startswith(build):
                geneStructureInformation = add_build(geneStructureInformation, build)
        shutil.copy(refGeneFile_pkl_path, output) #copy existing pickle file over
        ##############################################################
        #option3: ---------update existing annotation using bam file##
        ##############################################################
        if bamfile_path is not None:
            print('rely on bam file to update existing gene annotations')
            output_update = output[:-4]+'updated.pkl'
            if os.path.isfile(output_update):
                print('enhanced gene annotation already exist')
                geneStructureInformation = load_pickle(output_update)
            else:
                geneStructureInformation = annotate_genes(geneStructureInformation=geneStructureInformation,
                                                      bamfile_path=bamfile_path,
                                                      coverage_threshold_gene=coverage_threshold_gene,
                                                      coverage_threshold_exon=coverage_threshold_exon,
                                                      coverage_threshold_splicing=coverage_threshold_splicing,
                                                      z_score_threshold = z_score_threshold,
                                                      min_gene_size=min_gene_size, workers=num_cores)
                if output is not None:
                    with open(output[:-4]+'updated.pkl', 'wb') as file:
                        pickle.dump(geneStructureInformation, file)
    #########group genes into meta-genes########
    if os.path.isfile(meta_output) == False:
        print('meta gene information does not exist, will generate.')
        if genes is None: #bam generated reference
            geneIDs = list(geneStructureInformation.keys())
            rows = []
            for i, geneID in enumerate(geneIDs):
                row = [geneStructureInformation[geneID][0]['geneChr'], geneStructureInformation[geneID][0]['geneStart'],
                geneStructureInformation[geneID][0]['geneEnd'], geneStructureInformation[geneID][0]['geneID'],
                geneStructureInformation[geneID][0]['geneName'], '.','.', i]
                rows.append(row)
            genes = pd.DataFrame(rows)
        genes.columns = ['CHR', 'START', 'END', 'GENE_ID', 'GENE_NAME', 'GENE_TYPE', 'STRAND', 'META_GENE']
        grouped = genes.groupby("META_GENE")
        grouped_dict = {key: group.iloc[:, 3].tolist() for key, group in grouped}  # META_GENE STARTS FROM 1
        metageneStructureInformation = {}
        for i in range(len(grouped_dict)):
            meta_gene = 'meta_gene_' + str(i+1)
            gene_ids = grouped_dict[i+1]
            meta_gene_info = []
            for id in gene_ids:
                meta_gene_info.append(geneStructureInformation[id])
            metageneStructureInformation[meta_gene] = meta_gene_info
        # save to output, meta gene
        with open(meta_output, 'wb') as file:
            pickle.dump(metageneStructureInformation, file)
    else:
        print('load existing annotation pickle file of meta-gene')
        metageneStructureInformation = load_pickle(meta_output)
    return metageneStructureInformation

def add_build(geneStructureInformation, build):
    if build is not None:
        keys = list(geneStructureInformation.keys())
        for key in keys:
            geneStructureInformation[key][0]['geneChr'] = build + '_' + geneStructureInformation[key][0]['geneChr']
    return geneStructureInformation

class Annotator:
    def __init__(self, target:list, reference_gtf_path:str, reference_pkl_path:str, bam_path:list, update_gtf, workers,
                 coverage_threshold_gene, coverage_threshold_exon, coverage_threshold_splicing,z_score_threshold,
                 min_gene_size, build = None, platform = '10x', logger = None):
        """
        target: root path to save annotation files. SCOTCH will automatically create sub folders
        reference_gtf: path to gtf annotation (optional)
        bam_path: path to bam file, or path to the bam file folder; can be a str for a single sample, or a list for multiple samples
        update_gtf: whether to update gtf annotation using bam file
        build: parse parameter
        """
        self.logger = logger
        self.multiple_bam = True if len(bam_path)>1 else False
        self.multiple_samples = True if len(target)>1 else False
        self.workers = workers
        self.target = target
        self.reference_gtf_path = reference_gtf_path
        self.reference_pkl_path = reference_pkl_path
        self.bam_path = bam_path
        self.update_gtf = update_gtf
        self.build = build
        self.platform = platform
        self.parse = self.platform == 'parse'
        self.pacbio = self.platform == 'pacbio'
        # gene annotation information
        self.annotation_folder_path = [os.path.join(t, "reference") for t in target]
        self.annotation_path_single_gene = [os.path.join(t, "reference/geneStructureInformation.pkl") for t in target]
        self.annotation_path_meta_gene = [os.path.join(t, "reference/metageneStructureInformation.pkl") for t in target]
        # bam information
        self.bamInfo_folder_path = [os.path.join(t, "bam") for t in target]
        self.bamInfo_pkl_path = [os.path.join(t, 'bam/bam.Info.pkl') for t in target]
        self.bamInfo2_pkl_path = [os.path.join(t, 'bam/bam.Info2.pkl') for t in target]
        self.bamInfo3_pkl_path = [os.path.join(t, 'bam/bam.Info3.pkl') for t in target]  # only available for parse
        self.bamInfo_csv_path = [os.path.join(t, 'bam/bam.Info.csv') for t in target]
        # some parameters
        self.coverage_threshold_gene = coverage_threshold_gene
        self.coverage_threshold_exon = coverage_threshold_exon
        self.coverage_threshold_splicing = coverage_threshold_splicing
        self.z_score_threshold = z_score_threshold
        self.min_gene_size = min_gene_size
    def annotate_genes(self):
        for i in range(len(self.target)):
            if not os.path.exists(self.annotation_folder_path[i]):
                os.makedirs(self.annotation_folder_path[i])
            if os.path.isfile(self.annotation_path_single_gene[i]) and os.path.isfile(self.annotation_path_meta_gene[i]):
                self.logger.info(f'Complete gene annotation information exist at {self.annotation_path_single_gene[i]} and {self.annotation_path_meta_gene[i]}')
            else:
                if i==0:
                    self.logger.info(
                        f'Complete gene annotation information does not exist, we will generate...')
                    #annotation free mode
                    if self.reference_gtf_path is None and self.reference_pkl_path is None:
                        self.logger.info(f'Annotation-free Mode: we will rely on given bam files to generate gene annotations')
                        _ = extract_annotation_info(self.reference_gtf_path, self.reference_pkl_path,
                                                    self.bam_path, self.workers,
                                                    self.annotation_path_single_gene[i], self.build,
                                                    self.coverage_threshold_gene, self.coverage_threshold_exon,
                                                    self.coverage_threshold_splicing,self.z_score_threshold,
                                                    self.min_gene_size)
                    if self.update_gtf:
                        self.logger.info(
                            f'Enhanced-annotation Mode: we will update existing gene annotations using given bam files')
                        _ = extract_annotation_info(self.reference_gtf_path, self.reference_pkl_path,
                                                    self.bam_path, self.workers,
                                                    self.annotation_path_single_gene[i], self.build,
                                                    self.coverage_threshold_gene, self.coverage_threshold_exon,
                                                    self.coverage_threshold_splicing,self.z_score_threshold,
                                                    self.min_gene_size)
                    else:
                        self.logger.info(
                            f'Annotation-only Mode: we will only use existing gene annotations')
                        _ = extract_annotation_info(self.reference_gtf_path, self.reference_pkl_path,
                                                    None, self.workers,
                                                    self.annotation_path_single_gene[i], self.build,
                                                    self.coverage_threshold_gene, self.coverage_threshold_exon,
                                                    self.coverage_threshold_splicing,self.z_score_threshold,
                                                    self.min_gene_size)
                else:
                    # Copy the generated files from the first target to the current target
                    self.logger.info(f'Copying generated files from {self.annotation_path_single_gene[0]} to {self.annotation_path_single_gene[i]}')
                    shutil.copy(self.annotation_path_meta_gene[0], self.annotation_path_meta_gene[i])

    def annotation_bam(self):
        for i in range(len(self.target)):
            if not os.path.exists(self.bamInfo_folder_path[i]):
                os.makedirs(self.bamInfo_folder_path[i])
            if os.path.isfile(self.bamInfo_pkl_path[i]) == True and os.path.isfile(self.bamInfo_csv_path[i]) == True:
                self.logger.info(f'bam file information exist at {self.bamInfo_pkl_path[i]} and {self.bamInfo_csv_path[i]}')
            if os.path.isfile(self.bamInfo_pkl_path[i]) == False and os.path.isfile(self.bamInfo_csv_path[i]) == True:
                self.logger.info('Extracting bam file pickle information')
                bam_info = pd.read_csv(self.bamInfo_csv_path[i])
                qname_dict, qname_cbumi_dict, qname_sample_dict = bam_info_to_dict(bam_info, self.parse)
                with open(self.bamInfo_pkl_path[i], 'wb') as file:
                    pickle.dump(qname_dict, file)
                with open(self.bamInfo2_pkl_path[i], 'wb') as file:
                    pickle.dump(qname_cbumi_dict, file)
                if qname_sample_dict is not None:
                    with open(self.bamInfo3_pkl_path[i], 'wb') as file:
                        pickle.dump(qname_sample_dict, file)
            if os.path.isfile(self.bamInfo_pkl_path[i]) == False and os.path.isfile(self.bamInfo_csv_path[i]) == False:
                self.logger.info('Extracting bam file information')
                if os.path.isfile(self.bam_path[i])==False:
                    bam_info = extract_bam_info_folder(self.bam_path[i], self.workers, self.parse, self.pacbio)
                else:
                    if self.parse:
                        bam_info = extract_bam_info_parse(self.bam_path[i])
                    elif self.pacbio:
                        bam_info = extract_bam_info_pacbio(self.bam_path[i])
                    else:
                        bam_info = extract_bam_info(self.bam_path[i])
                bam_info.to_csv(self.bamInfo_csv_path[i])
                self.logger.info('Generating bam file pickle information')
                qname_dict, qname_cbumi_dict, qname_sample_dict = bam_info_to_dict(bam_info, self.parse)
                with open(self.bamInfo_pkl_path[i], 'wb') as file:
                    pickle.dump(qname_dict, file)
                with open(self.bamInfo2_pkl_path[i], 'wb') as file:
                    pickle.dump(qname_cbumi_dict, file)
                if qname_sample_dict is not None:
                    with open(self.bamInfo3_pkl_path[i], 'wb') as file:
                        pickle.dump(qname_sample_dict, file)

























