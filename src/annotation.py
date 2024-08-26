from collections import defaultdict
import pysam
from joblib import Parallel, delayed
import reference as ref
import os
import pandas as pd
import pickle
from preprocessing import load_pickle
import tqdm

######################################################################
##############################annotation##############################
######################################################################

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



###main function
##TODO: update existing annotation, generate annotation, no repeat read bam file
def get_non_overlapping_exons(bam_file, chrom, gene_start, gene_end, coverage_threshold=20):
    bam = pysam.AlignmentFile(bam_file, "rb")
    exons, coverage, junctions = [], {}, []
    #get read coverage
    for read in bam.fetch(chrom, gene_start, gene_end):
        if read.reference_start >= gene_start and read.reference_end < gene_end:
            junctions += get_splicejuction_from_read(read)
            for read_start, read_end in read.get_blocks():
                for pos in range(read_start, read_end):
                    if pos in coverage:
                        coverage[pos] += 1
                    else:
                        coverage[pos] = 1
    #only keep positions that meet the coverage threshold
    coverage = {pos: cov for pos, cov in coverage.items() if cov > coverage_threshold}
    coverage_blocks = get_continuous_blocks(coverage)
    coverage_blocks = [(a,b) for a, b in coverage_blocks if b-a >= 20]#delete less than 20bp exons
    if len(coverage_blocks)==0:
        return exons
    #use splicing junctions to find sub-exons
    boundaries = count_and_sort_tuple_elements(junctions)
    filtered_boundaries = [{} for _ in coverage_blocks]
    for boundary, freq in boundaries.items():
        for i, (start, end) in enumerate(coverage_blocks):
            if start < boundary < end:
                filtered_boundaries[i][boundary] = freq
                break
    Filtered_Boundaries = []
    for filtered_boundaries_ in filtered_boundaries:
        merged_boundaries = merge_boundaries_by_evidence(filtered_boundaries_, merge_distance=10)
        sorted_boundaries = dict(sorted(merged_boundaries.items()))
        filtered_sorted_boundaries = {k: v for k, v in sorted_boundaries.items() if v > 20}
        Filtered_Boundaries.append(filtered_sorted_boundaries)
    for i in range(len(coverage_blocks)):
        start, end = coverage_blocks[i]
        positions = [start]+list(Filtered_Boundaries[i].keys())+[end]
        for ii in range(len(positions)-1):
            exons.append((positions[ii], positions[ii+1]))
    bam.close()
    return exons


def get_genes_from_bam(input_bam_path, coverage_threshold = 5, min_region_size=50, workers = 1):
    def process_bam_file(bam_file, coverage_threshold, min_region_size):
        bam = pysam.AlignmentFile(bam_file, "rb")
        coverage = defaultdict(lambda: defaultdict(int))
        chromosomes = bam.references
        for chrom in chromosomes:
            for read in bam.fetch(chrom):
                if not read.is_unmapped:
                    for pos in range(read.reference_start, read.reference_end):
                        coverage[chrom][pos] += 1
        bam.close()
        genes = {}
        for chrom, cov_dict in coverage.items():
            cov_dict = {pos: cov for pos, cov in cov_dict.items() if cov > coverage_threshold}
            coverage_blocks = get_continuous_blocks(cov_dict)
            genes[chrom] = [(a, b) for a, b in coverage_blocks if b - a > min_region_size]
        return genes
    if os.path.isdir(input_bam_path):
        bam_files = [os.path.join(input_bam_path, f) for f in os.listdir(input_bam_path) if f.endswith('.bam')]
    else:
        bam_files = [input_bam_path]
    results = Parallel(n_jobs=workers)(delayed(process_bam_file)(bam_file, coverage_threshold, min_region_size) for bam_file in bam_files)
    all_genes={}
    for result in results:
        for chrom, gene_list in result.items():
            if chrom not in all_genes:
                all_genes[chrom] = gene_list
            else:
                all_genes[chrom].extend(gene_list)
    sorted_all_genes = {}
    for chrom in sorted(all_genes.keys()):
        sorted_all_genes[chrom] = sorted(all_genes[chrom],
                                         key=lambda x: (x[0], x[1]))
    return sorted_all_genes


def update_exons(exonInfo_annotation, exonInfo_annotation_free):
    if len(exonInfo_annotation)==1:
        return exonInfo_annotation
    updated_exonInfo = list(exonInfo_annotation)
    for i in range(len(exonInfo_annotation) - 1):
        end_of_current = exonInfo_annotation[i][1]
        start_of_next = exonInfo_annotation[i + 1][0]
        if end_of_current < start_of_next:
            gap_start = end_of_current
            gap_end = start_of_next
            for exon in exonInfo_annotation_free:
                exon_start, exon_end = exon
                # Find overlap with the gap
                if (exon_start>=gap_start and exon_start <= gap_end) or (exon_end >= gap_start and exon_end <= gap_end):
                    segment_start = max(gap_start, exon_start)
                    segment_end = min(gap_end, exon_end)
                    if segment_start < segment_end:
                        updated_exonInfo.append((segment_start, segment_end))
    updated_exonInfo.sort()
    return updated_exonInfo



def annotate_genes(geneStructureInformation, bamfile_path,
                   coverage_threshold_gene = 5, coverage_threshold_exon = 20,
                   min_gene_size=50, workers=1):
    """
    generate geneStructureInformation either using bamfile alone (leave geneStructureInformation blank) or update existing annotation file using bam file
    :param geneStructureInformation: pickle file object of existing gene annotation, not path
    :param bamfile_path: can be a folder path of multiple bams or a single bam file
    :param coverage_threshold_gene: minimal read coverage for gene
    :param coverage_threshold_exon: minimal read coverage for exon
    :param min_gene_size: minimal length for novel gene discovery
    :return: updated geneStructureInformation annotation
    """
    def novel_gene_annotation(chrom, all_genes_dict, bamfile, coverage_threshold_exon=20):
        #all_genes: {'chr1':[(100,200),(300,400)]}, bamfile: path to one bam file
        if os.path.isfile(bamfile) == False:  # bamfile is a folder
            bamFile_name = [f for f in os.listdir(bamfile) if
                            f.endswith('.bam') and '.' + str(chrom) + '.' in f]
            bamfile = os.path.join(bamfile, bamFile_name[0])
        genes = all_genes_dict[chrom]  # [(100,200),(300,400)]
        gene_annotations = {}
        for i, (gene_start, gene_end) in enumerate(genes):
            exonInfo = get_non_overlapping_exons(bamfile, chrom, gene_start, gene_end, coverage_threshold_exon)
            geneInfo = {'geneName': 'gene_'+ str(i+1)+'_'+chrom, 'geneID': 'gene_'+ str(i+1)+'_'+chrom,
                        'geneChr': chrom, 'geneStart':gene_start, 'geneEnd':gene_end, 'geneStrand': '.',
                        'numofExons': len(exonInfo), 'numofIsoforms': 0, 'isoformNames':[]}
            isoformInfo = {}
            gene_annotations['gene_'+ str(i+1)+'_'+chrom] = [geneInfo, exonInfo, isoformInfo]
        return gene_annotations
    def update_isoform_info(original_exons, updated_exons, isoformInfo):
        exon_mapping = {}
        for original_idx, original_exon in enumerate(original_exons):
            for updated_idx, updated_exon in enumerate(updated_exons):
                if original_exon == updated_exon:
                    exon_mapping[original_idx] = updated_idx
        updated_isoform_info = {}
        for isoform, exon_indices in isoformInfo.items():
            updated_indices = [exon_mapping[idx] for idx in exon_indices if idx in exon_mapping]
            updated_isoform_info[isoform] = updated_indices
        return updated_isoform_info
    def update_annotation(geneStructureInformation, geneID, bamfile_path,coverage_threshold_exon):
        chrom, gene_start, gene_end = geneStructureInformation[geneID][0]['geneChr'], \
        geneStructureInformation[geneID][0]['geneStart'], geneStructureInformation[geneID][0]['geneEnd']
        exons_bam = get_non_overlapping_exons(bamfile_path, chrom, gene_start, gene_end, coverage_threshold_exon)
        original_exons = geneStructureInformation[geneID][1]
        updated_exons = update_exons(original_exons, exons_bam)
        # geneInfo
        geneInfo = geneStructureInformation[geneID][0]
        geneInfo['numofExons'] = len(updated_exons)
        # exonInfo
        exonInfo = updated_exons
        # isoformInfo
        isoformInfo = update_isoform_info(original_exons, updated_exons, geneStructureInformation[geneID][2])
        return {geneID:[geneInfo, exonInfo, isoformInfo]}
    #generate gene annotation solely based on bam file
    if geneStructureInformation is None:
        all_genes = get_genes_from_bam(bamfile_path, coverage_threshold_gene, min_gene_size) #{'chr1':[(100,200),(300,400)]}
        chroms = list(all_genes.keys())
        results = Parallel(n_jobs=workers)(
            delayed(novel_gene_annotation)(chrom, all_genes, bamfile_path, coverage_threshold_exon) for chrom in chroms)
        annotations = {k: v for result in results for k, v in result.items()}
    #update existing gene annotation using bam file
    else:
        geneIDs = list(geneStructureInformation.keys())
        results = Parallel(n_jobs=workers)(
            delayed(update_annotation)(geneStructureInformation, geneID, bamfile_path,coverage_threshold_exon) for geneID in geneIDs)
        annotations = {k: v for result in results for k, v in result.items()}
    return annotations

def extract_annotation_info(refGeneFile_path, bamfile_path, num_cores=8,
                            output="geneStructureInformation.pkl", build=None,
                            coverage_threshold_gene=5, coverage_threshold_exon=20,min_gene_size=50):
    """
    extract gene annotation information
    :param refGeneFile_path: path to gene annotation gtf file
    :param bamfile_path: path to bam file or the folder for bam files
    :param num_cores: number of workers for parallel computing
    :param output: default: geneStructureInformation.pkl
    :param build: used for parse platform namingss
    :return: metageneStructureInformation
    """
    geneStructureInformation = None
    meta_output = os.path.join(os.path.dirname(output), 'meta' + os.path.basename(output))
    #####################################################
    #option1: ---------rely on bam file alone---------###
    #####################################################
    print('rely on bam file alone to generate gene annotations')
    if refGeneFile_path is None and bamfile_path is not None:
        geneStructureInformation = annotate_genes(geneStructureInformation = None, bamfile_path = bamfile_path,
                       coverage_threshold_gene=coverage_threshold_gene,
                       coverage_threshold_exon=coverage_threshold_exon,
                       min_gene_size=min_gene_size, workers=num_cores)
    #####################################################
    # option2: --------rely on existing annotation alone#
    #####################################################
    print('rely on existing gene annotations')
    if refGeneFile_path is not None:
        genes, exons = ref.generate_reference_df(gtf_path=refGeneFile_path)
        Genes = list(zip(genes.iloc[:, 3].tolist(), genes.iloc[:, 4].tolist()))  # id, name
        #generate single gene annotations if not existing
        if os.path.isfile(output) == False:
            geneStructureInformation = Parallel(n_jobs=num_cores)(
                delayed(process_gene)(geneID, geneName, genes, exons, build) for geneID, geneName in tqdm(Genes))
            geneStructureInformation = dict(geneStructureInformation)
            #save to output, single gene
            if output is not None:
                with open(output, 'wb') as file:
                    pickle.dump(geneStructureInformation, file)
        else:
            print('load existing annotation pickle file of each single gene')
            geneStructureInformation = load_pickle(output)
        ##############################################################
        #option3: ---------update existing annotation using bam file##
        ##############################################################
        print('rely on bam file to update existing gene annotations')
        if bamfile_path is not None:
            geneStructureInformation = annotate_genes(geneStructureInformation=geneStructureInformation,
                                                      bamfile_path=bamfile_path,
                                                      coverage_threshold_gene=coverage_threshold_gene,
                                                      coverage_threshold_exon=coverage_threshold_exon,
                                                      min_gene_size=min_gene_size, workers=num_cores)
    #########group genes into meta-genes########
    if os.path.isfile(meta_output) == False:
        # group genes
        genes.columns = ['CHR', 'START', 'END', 'GENE_ID', 'GENE_NAME', 'GENE_TYPE', 'STRAND', 'META_GENE']
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
        # save to output, meta gene
        with open(meta_output, 'wb') as file:
            pickle.dump(metageneStructureInformation, file)
    else:
        print('load existing annotation pickle file of meta-gene')
        metageneStructureInformation = load_pickle(meta_output)
    return metageneStructureInformation


class Annotator:
    def __init__(self, target, reference_gtf_path, bam_path, update_gtf, workers,
                 coverage_threshold_gene, coverage_threshold_exon, min_gene_size):
        """
        target: root path to save annotation files. SCOTCH will automatically create sub folders
        reference_gtf: path to gtf annotation (optional)
        bam_path: path to bam file, or path to the bam file folder
        update_gtf: whether to update gtf annotation using bam file
        """
        self.workers = workers
        self.target = target
        self.reference_gtf_path = reference_gtf_path
        self.bam_path = bam_path
        self.update_gtf = update_gtf
        # gene annotation information
        self.annotation_folder_path = os.path.join(target, "reference")
        self.annotation_path_single_gene = os.path.join(target, "reference/geneStructureInformation.pkl")
        self.annotation_path_meta_gene = os.path.join(target, "reference/metageneStructureInformation.pkl")
        # bam information
        self.bamInfo_folder_path = os.path.join(target, "bam")
        self.bamInfo_pkl_path = os.path.join(target, 'bam/bam.Info.pkl')
        self.bamInfo2_pkl_path = os.path.join(target, 'bam/bam.Info2.pkl')
        self.bamInfo_csv_path = os.path.join(target, 'bam/bam.Info.csv')
        # some parameters
        self.coverage_threshold_gene = coverage_threshold_gene
        self.coverage_threshold_exon = coverage_threshold_exon
        self.min_gene_size = min_gene_size
    def annotate_genes(self):
        if not os.path.exists(self.annotation_folder_path):
            os.makedirs(self.annotation_folder_path)
        if os.path.isfile(self.annotation_path_single_gene) and os.path.isfile(self.annotation_path_meta_gene):
            print('gene annotation information exist')
        else:
            print('gene annotation information does not exist, we will generate')
            #annotation free mode
            if self.reference_gtf_path is None:
                print('annotation free mode')
                _ = extract_annotation_info(None, self.bam_path, self.workers,
                                "geneStructureInformation.pkl", None,
                                self.coverage_threshold_gene, self.coverage_threshold_exon,
                                            self.min_gene_size)
            if self.update_gtf:
                print('update existing annotation using bam file')
                _ = extract_annotation_info(self.reference_gtf_path, self.bam_path, self.workers,
                                            "geneStructureInformation.pkl", None,
                                            self.coverage_threshold_gene, self.coverage_threshold_exon,
                                            self.min_gene_size)
            else:
                print('using existing annotation file')
                _ = extract_annotation_info(self.reference_gtf_path, None, self.workers,
                                            "geneStructureInformation.pkl", None,
                                            self.coverage_threshold_gene, self.coverage_threshold_exon,
                                            self.min_gene_size)
    def annotation_bam(self):
        if not os.path.exists(self.bamInfo_folder_path):
            os.makedirs(self.bamInfo_folder_path)
        if os.path.isfile(self.bamInfo_pkl_path) == True and os.path.isfile(self.bamInfo_csv_path) == True:
            print('bam file information exist')
        if os.path.isfile(self.bamInfo_pkl_path) == False and os.path.isfile(self.bamInfo_csv_path) == True:
            print('extracting bam file pickle information')
            bam_info = pd.read_csv(self.bamInfo_csv_path)
            qname_dict, qname_cbumi_dict = bam_info_to_dict(bam_info)
            with open(self.bamInfo_pkl_path, 'wb') as file:
                pickle.dump(qname_dict, file)
            with open(self.bamInfo2_pkl_path, 'wb') as file:
                pickle.dump(qname_cbumi_dict, file)
        if os.path.isfile(self.bamInfo_pkl_path) == False and os.path.isfile(self.bamInfo_csv_path) == False:
            print('extracting bam file information')
            if os.path.isfile(self.bam_path)==False:
                bam_info = extract_bam_info_folder(self.bam_path, self.workers)
            else:
                bam_info = extract_bam_info(self.bam_path)
            bam_info.to_csv(self.bamInfo_csv_path)
            print('generating bam file pickle information')
            qname_dict, qname_cbumi_dict = bam_info_to_dict(bam_info)
            with open(self.bamInfo_pkl_path, 'wb') as file:
                pickle.dump(qname_dict, file)
            with open(self.bamInfo2_pkl_path, 'wb') as file:
                pickle.dump(qname_cbumi_dict, file)

























