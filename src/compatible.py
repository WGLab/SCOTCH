import pickle

from preprocessing import *
import pysam
import re
import csv
import math
import copy
from joblib import Parallel, delayed
def convert_to_gtf(metageneStructureInformationNovel, output_file, gtf_df = None, num_cores=1):
    def update_annotation_gene(geneID, gtf_df, geneStructureInformationwNovel):
        gtf_df_sub = gtf_df[gtf_df['attribute'].str.contains(f'gene_id "{geneID}"', regex=False)].reset_index(drop=True)
        geneInfo, exonInfo, isoformInfo = geneStructureInformationwNovel[geneID]
        if gtf_df_sub.shape[0] > 0:
            new_row_gene = []
        else:
            new_row_gene = [[geneInfo['geneChr'], 'SCOTCH', 'gene', geneInfo['geneStart']+1, geneInfo['geneEnd'], '.',
                            geneInfo['geneStrand'], '.', f'gene_id "{geneID}"; gene_name "{geneInfo["geneName"]}"']]
        new_rows_isoforms = []
        for isoform_name, exon_indices in isoformInfo.items():
            isoform_start = exonInfo[exon_indices[0]][0]
            isoform_end = exonInfo[exon_indices[-1]][1]
            if gtf_df_sub[gtf_df_sub['attribute'].str.contains(f'transcript_id "{isoform_name}"', regex=False)].shape[0] == 0:
                new_rows_isoform = [geneInfo['geneChr'], 'SCOTCH', 'transcript', isoform_start+1, isoform_end, '.',
                                    geneInfo['geneStrand'], '.',
                                    f'gene_id "{geneID}"; gene_name "{geneInfo["geneName"]}"; transcript_id "{isoform_name}"; transcript_name "{isoform_name}"']
                new_rows_isoforms.append(new_rows_isoform)
                exons_isoform = []
                for exon_index in exon_indices:
                    exons_isoform.append(exonInfo[exon_index])
                merged_exons_isoform = merge_exons(exons_isoform)
                for exon_num, (exon_start, exon_end) in enumerate(merged_exons_isoform, start=1):
                    new_rows_exon = [geneInfo['geneChr'], 'SCOTCH', 'exon', exon_start+1, exon_end, '.',
                                     geneInfo['geneStrand'], '.',
                                     f'gene_id "{geneID}"; gene_name "{geneInfo["geneName"]}"; transcript_id "{isoform_name}"; transcript_name "{isoform_name}"; exon_number "{exon_num}"']
                    new_rows_isoforms.append(new_rows_exon)
        new_row = new_row_gene + new_rows_isoforms
        gtf_df_sub_new = pd.DataFrame(new_row, columns=column_names)
        gtf_df_gene = pd.concat([gtf_df_sub, gtf_df_sub_new], ignore_index=True)
        return gtf_df_gene
    def partition_metagene(metageneStructureInformation):
        meta_gene_names = list(metageneStructureInformation.keys())
        geneStructureInformation = {}
        for meta_gene in meta_gene_names:
            genes_info_list = metageneStructureInformation[meta_gene]
            for gene_info in genes_info_list:
                geneInfo, exonInfo, isoformInfo = gene_info
                geneID = geneInfo['geneID']
                geneStructureInformation[geneID] = [geneInfo, exonInfo, isoformInfo]
        return geneStructureInformation
    geneStructureInformationwNovel = partition_metagene(metageneStructureInformationNovel)
    geneIDs = list(geneStructureInformationwNovel.keys())
    column_names = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    if gtf_df is None:
        gtf_df = pd.DataFrame(columns=column_names)
    gtf_df_gene_list = Parallel(n_jobs=num_cores)(delayed(update_annotation_gene)(geneID, gtf_df, geneStructureInformationwNovel) for geneID in geneIDs)
    gtf_df_geness = pd.concat(gtf_df_gene_list, ignore_index=True)
    gtf_df_geness.to_csv(output_file, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)



def summarise_annotation(target):
    reference_folders = []
    for root, dirs, files in os.walk(target):
        if 'reference' in dirs:
            reference_folders.append(os.path.join(root, 'reference'))
    def get_numeric_key(key):
        return int(key.split('_')[-1])
    for reference_folder in reference_folders:
        output_pkl = os.path.join(reference_folder, "metageneStructureInformationwNovel.pkl")
        output_gtf = os.path.join(reference_folder, "SCOTCH_updated_annotation.gtf")
        if os.path.isfile(output_gtf):
            continue
        file_names_pkl = [os.path.join(reference_folder, f) for f in os.listdir(reference_folder) if
                          re.match(r'metageneStructureInformationwNovel_\d+\.pkl', f)]
        file_names_gtf = [os.path.join(reference_folder, f) for f in os.listdir(reference_folder) if
                          re.match(r'gene_annotations_scotch_\d+\.gtf', f)]
        if len(file_names_pkl)>0 and len(file_names_gtf)>0:
            # merge pkl annotation file
            print('merging new isoform annotations')
            metageneStructureInformationwNovel = {}
            for file_name_pkl in file_names_pkl:
                metageneStructureInformation = load_pickle(file_name_pkl)
                metageneStructureInformationwNovel.update(metageneStructureInformation)
            metageneStructureInformationwNovel = dict(
                sorted(metageneStructureInformationwNovel.items(), key=lambda item: get_numeric_key(item[0]))
            )
            with open(output_pkl, 'wb') as file:
                pickle.dump(metageneStructureInformationwNovel, file)
            for file_name_pkl in file_names_pkl:
                os.remove(file_name_pkl)
            print('mergered new isoform annotation saved at: '+str(output_pkl))
            # merge gtf annotation file
            print('Merging new GTF annotations...')
            gtf_lines = []
            for file_name_gtf in file_names_gtf:
                with open(file_name_gtf, 'r') as gtf_file:
                    for line in gtf_file:
                        if not line.startswith('#'):  # Skip any header lines
                            gtf_lines.append(line.strip())
            with open(output_gtf, 'w') as output_gtf_file:
                for line in gtf_lines:
                    output_gtf_file.write(line + '\n')
            for file_name_gtf in file_names_gtf:
                os.remove(file_name_gtf)
            print('Merged GTF annotations saved at: ' + output_gtf)
        else:
            print('novel isoform annotations does not exist!')

def summarise_auxillary(target):
    def process_group(group):
        highest_priority = group['priority'].max()
        highest_priority_rows = group[group['priority'] == highest_priority]
        group['GeneMapping'] = 'delete'
        group['Keep'] = 0
        if len(highest_priority_rows) == 1:
            group.loc[highest_priority_rows.index, 'GeneMapping'] = 'unique'
            group.loc[highest_priority_rows.index, 'Keep'] = 1
        else:
            group.loc[highest_priority_rows.index, 'GeneMapping'] = 'ambiguous'
            random_index = highest_priority_rows.sample(n=1).index
            group.loc[random_index, 'Keep'] = 1
        return group
    def read_file(file_path):
        df = pd.read_csv(file_path, sep='\t')
        df['gene'] = df['geneName'] + '_' + df['geneID']
        return df
    # Collect all 'auxillary' directories
    auxillary_folders = []
    for root, dirs, files in os.walk(target):
        if 'auxillary' in dirs:
            auxillary_folders.append(os.path.join(root, 'auxillary'))
    for auxillary_folder in auxillary_folders:
        print('summarising read-isoform mapping files at: ' + str(auxillary_folder))
        file_paths = [os.path.join(auxillary_folder, f) for f in os.listdir(auxillary_folder) if 'ENSG' in f]
        df_list = Parallel(n_jobs=-1)(delayed(read_file)(file_path) for file_path in file_paths)
        DF = pd.concat(df_list, axis=0, ignore_index=True).reset_index(drop=True)
        DF['MappingScore'] = DF['MappingScore'].fillna(-1)
        conditions = [
            DF['Isoform'].str.startswith('ENST'),  # 1 Highest priority
            DF['Isoform'].str.startswith('novel'),  # 2 Medium priority
            DF['Isoform'] == 'uncategorized'  # 3 Lowest priority
        ]
        choices = [1, 2, 3]
        DF['priority'] = np.select(conditions, choices, default=3)
        DF['priority'] *= DF['MappingScore']
        DF['GeneMapping'] = 'delete'  # Initialize as 'delete'
        DF = DF.sort_values(by=['geneChr', 'Read', 'priority'], ascending=[True, True, False])
        #Split DF into DF_unique and DF_multiple
        unique_mask = ~DF['Read'].duplicated(keep=False)  # Reads that appear only once
        multiple_mask = DF['Read'].duplicated(keep=False)  # Reads that appear more than once
        DF_unique = DF[unique_mask].copy()
        DF_unique['GeneMapping'] = 'unique'
        DF_unique['Keep'] = 1
        DF_multiple = DF[multiple_mask].copy()
        grouped = DF_multiple.groupby(['Read'], group_keys=False)
        # Process groups in parallel and concatenate results
        print('Processing groups for multiple reads...')
        processed_groups = Parallel(n_jobs=-1)(delayed(process_group)(group) for _, group in grouped)
        DF_multiple_processed = pd.concat(processed_groups).reset_index(drop=True)
        DF_final = pd.concat([DF_unique, DF_multiple_processed], ignore_index=True)
        output_file_tsv = os.path.join(auxillary_folder, 'all_read_isoform_exon_mapping.tsv')
        print('saving read-isoform mapping file: '+str(output_file_tsv))
        DF_final.to_csv(output_file_tsv, sep='\t', index=False)
        print('removing temporary files in: '+ str(auxillary_folder))
        for file in file_paths:
            os.remove(file)
        output_file_pkl = os.path.join(auxillary_folder, 'read_selection.pkl')
        cbumi_keep_dict = DF_final.set_index('CBUMI')['Keep'].to_dict()
        print('saving read filtering file: ' + str(output_file_pkl))
        with open(output_file_pkl, 'wb') as pickle_file:
            pickle.dump(cbumi_keep_dict, pickle_file)







class ReadMapper:
    def __init__(self, target:list, bam_path:list, lowest_match=0.2, lowest_match1 = 0.6, small_exon_threshold = 0,
                 small_exon_threshold1=80, truncation_match=0.4, platform = '10x-ont',
                 reference_gtf_path = None, logger = None):
        self.logger = logger
        self.target = target
        self.bam_path = bam_path
        column_names = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        if reference_gtf_path is not None or reference_gtf_path=='None':
            self.gtf_df = pd.read_csv(reference_gtf_path, sep='\t', comment='#', header=None, names=column_names)
        else:
            self.gtf_df = pd.DataFrame(columns=column_names)
        # gene annotation information
        self.annotation_folder_path_list = [os.path.join(target_, 'reference') for target_ in target]
        self.annotation_path_single_gene_list = [os.path.join(target_, 'reference/geneStructureInformation.pkl') for target_ in target]
        self.annotation_path_meta_gene_list = [os.path.join(target_, "reference/metageneStructureInformation.pkl") for target_ in target]
        self.annotation_path_meta_gene_novel_list = [os.path.join(target_, "reference/metageneStructureInformationwNovel.pkl") for target_ in target]
        self.annotation_path_gtf_novel_list = [os.path.join(target_, "reference/gene_annotations_scotch.gtf") for target_ in target]
        # bam information path
        self.bamInfo_folder_path_list = [os.path.join(target_, "bam") for target_ in target]
        self.bamInfo_pkl_path_list = [os.path.join(target_, 'bam/bam.Info.pkl') for target_ in target]#bamInfo_pkl_file
        self.bamInfo2_pkl_path_list = [os.path.join(target_, 'bam/bam.Info2.pkl') for target_ in target]#bamInfo2_pkl_file
        self.bamInfo3_pkl_path_list = [os.path.join(target_, 'bam/bam.Info3.pkl') for target_ in target] # bamInfo2_pkl_file
        self.bamInfo_csv_path_list = [os.path.join(target_, 'bam/bam.Info.csv') for target_ in target]
        # parameters
        self.small_exon_threshold = small_exon_threshold
        self.small_exon_threshold1 = small_exon_threshold1
        self.truncation_match = truncation_match
        self.lowest_match = lowest_match
        self.lowest_match1 = lowest_match1
        self.platform = platform
        self.parse = self.platform == 'parse-ont'
        self.pacbio = self.platform == '10x-pacbio'
        # some paths
        if platform != 'parse-ont':
            self.nsamples = len(self.target)
            self.compatible_matrix_folder_path_list = [os.path.join(target_, "compatible_matrix") for target_ in target] #not for parse
            self.read_mapping_path_list = [os.path.join(target_, "auxillary") for target_ in target] #not for parse
            self.qname_dict_list = [load_pickle(bamInfo_pkl_path) for bamInfo_pkl_path in self.bamInfo_pkl_path_list]
            self.qname_cbumi_dict_list = [load_pickle(bamInfo2_pkl_path) for bamInfo2_pkl_path in self.bamInfo2_pkl_path_list]
            self.sorted_bam_path_list = None
            self.qname_sample_dict_list = [load_pickle(bamInfo3_pkl_path) for bamInfo3_pkl_path in self.bamInfo3_pkl_path_list]
        else:
            self.qname_dict = load_pickle(self.bamInfo_pkl_path_list[0])
            self.qname_cbumi_dict = load_pickle(self.bamInfo2_pkl_path_list[0])
            self.sorted_bam_path = None
            self.qname_sample_dict = load_pickle(self.bamInfo3_pkl_path_list[0])
        self.metageneStructureInformation = load_pickle(self.annotation_path_meta_gene_list[0])
        self.metageneStructureInformationwNovel = self.metageneStructureInformation.copy()
    def read_bam(self, chrom = None):
        #if parse: the input length is 1
        # bam_path is a folder
        bamFilePysam_list = []
        for bam_path in self.bam_path:
            if os.path.isfile(bam_path) == False: # If it's a folder, find the BAM file based on chrom
                if chrom is not None:
                    bamFile_name = [f for f in os.listdir(bam_path) if f.endswith('.bam') and '.' + chrom + '.' in f]
                    if bamFile_name:  # Ensure a matching BAM file is found
                        bamFile = os.path.join(bam_path, bamFile_name[0])  # Not .bai
                        bamFilePysam = pysam.Samfile(bamFile, "rb")
                        bamFilePysam_list.append(bamFilePysam)
                else:
                    #read the merged bam file, has to run merge_bam first
                    bamFilePysam = pysam.Samfile(self.sorted_bam_path, "rb")
                    bamFilePysam_list.append(bamFilePysam)
            else:
                # If it's a BAM file path, read it directly
                bamFilePysam = pysam.Samfile(bam_path, "rb")
                bamFilePysam_list.append(bamFilePysam)
        if self.parse:
            bamFilePysam_list = bamFilePysam_list[0] #not a list
        return bamFilePysam_list
    def map_reads(self, meta_gene, save = True):
        Info_multigenes = copy.deepcopy(self.metageneStructureInformation[meta_gene])
        Info_multigenes = sort_multigeneInfo(Info_multigenes)
        bamFilePysams = self.read_bam(chrom=Info_multigenes[0][0]['geneChr'])
        if len(Info_multigenes)==1:
            Info_singlegene = Info_multigenes[0]
            geneInfo, exonInfo, isoformInfo = Info_singlegene
            n_isoforms = len(isoformInfo)
            Read_novelIsoform = []  # [('read name',[read-exon percentage],[read-exon mapping])]
            Read_knownIsoform = []  # [('read name',[read-isoform mapping])]
            Read_knownIsoform_scores = {}  # [readname: read-isoform mapping scores]
            novel_isoformInfo = {}  # {'novelIsoform_1234':[2,3,4]}
            qname_sample_dict = {}
            for i, bamFilePysam in enumerate(bamFilePysams):
                sample = 'sample' + str(i)
                reads = bamFilePysam.fetch(geneInfo['geneChr'], geneInfo['geneStart'], geneInfo['geneEnd'])
                for read in reads:
                    readName, readStart, readEnd = read.qname, read.reference_start, read.reference_end
                    result = process_read(read, self.qname_dict_list[i], self.lowest_match, self.lowest_match1,
                                          self.small_exon_threshold,
                                          self.small_exon_threshold1, self.truncation_match, Info_singlegene,
                                          self.parse, self.pacbio)
                    if result is not None:
                        if self.pacbio:
                            readName = readName + '_' + str(readEnd - readStart)
                        qname_sample_dict[readName] = sample
                    result_novel, result_known, result_known_scores = result
                    if result_novel is not None:
                        Read_novelIsoform.append(result_novel)
                    if result_known is not None:
                        Read_knownIsoform.append(result_known)
                        Read_knownIsoform_scores[result_known[0]] = result_known_scores
            #expand uncategorized novel reads into Read_knownIsoform
            if len(Read_novelIsoform) > 0:
                #novel_isoformInfo_polished: novel isoform annotation: {'novelIsoform_7':[0,1,2]}
                Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished = polish_compatible_vectors(
                    Read_novelIsoform, Read_knownIsoform, n_isoforms, exonInfo, self.small_exon_threshold, self.small_exon_threshold1)
            else:
                Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished = (Read_novelIsoform, novel_isoformInfo,
                                                                                                      Read_knownIsoform)
            #compile output into compatible matrix
            geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(
                    Read_novelIsoform_polished, Read_knownIsoform_polished, geneInfo)
            #update annotation information in self
            self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'] = \
                self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames']+list(novel_isoformInfo_polished.keys())
            self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] = \
                self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] + len(list(
                    novel_isoformInfo_polished.keys()))
            self.metageneStructureInformationwNovel[meta_gene][0][2].update(novel_isoformInfo_polished)
            for key in Read_knownIsoform_scores:
                Read_knownIsoform_scores[key].extend([-1] * len(list(novel_isoformInfo_polished.keys())))
            sample_list = []
            for sample_ind in range(self.nsamples):
                sample_target = self.target[sample_ind]
                Read_Isoform_compatibleVector_sample, Read_knownIsoform_scores_sample = {}, {}
                for readname in list(Read_Isoform_compatibleVector.keys()):
                    if qname_sample_dict[readname] == 'sample'+str(sample_ind):
                        Read_Isoform_compatibleVector_sample[readname] = Read_Isoform_compatibleVector[readname]
                        if readname in Read_knownIsoform_scores.keys():
                            Read_knownIsoform_scores_sample[readname] = Read_knownIsoform_scores[readname]
                if save:
                    # save compatible matrix of each gene, save read-isoform mappings
                    save_compatibleVector_by_gene(geneName, geneID, geneChr, colNames,
                                                  Read_Isoform_compatibleVector_sample, Read_knownIsoform_scores_sample,
                                                  self.qname_cbumi_dict_list[sample_ind], self.metageneStructureInformationwNovel[meta_gene][0][1],
                                                  self.metageneStructureInformationwNovel[meta_gene][0][2], sample_target,
                                                  self.parse)
                else:
                    sample_list.append({'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector_sample, 'isoforms': colNames,
                     'exonInfo': self.metageneStructureInformationwNovel[meta_gene][0][1],
                    'isoformInfo':self.metageneStructureInformationwNovel[meta_gene][0][2]})
            if not save:
                return sample_list
        else:
            geneChr, start, end = summarise_metagene(Info_multigenes)  # geneChr, start, end
            qname_sample_dict={}
            results = []
            for i, bamFilePysam in enumerate(bamFilePysams):
                reads = bamFilePysam.fetch(geneChr, start, end)  # fetch reads within meta gene region
                # process reads metagene
                for read in reads:
                    readName, readStart, readEnd = read.qname, read.reference_start, read.reference_end
                    out = process_read_metagene(read,self.qname_dict_list[i], Info_multigenes, self.lowest_match,self.lowest_match1,
                                                self.small_exon_threshold,self.small_exon_threshold1,
                                                self.truncation_match, self.parse, self.pacbio)
                    if out is not None: #may not within this meta gene region
                        results.append(out)
                        if self.pacbio:
                            readName = readName + '_' + str(readEnd - readStart)
                        qname_sample_dict[readName] = 'sample' + str(i)
            Ind, Read_novelIsoform_metagene, Read_knownIsoform_metagene, Read_knownIsoform_metagene_scores = [], [], [], []
            for result in results:
                if result is not None:
                    ind, novelisoform, knownisoform, knownisoformscores = result
                    Ind.append(ind)
                    Read_novelIsoform_metagene.append(novelisoform)
                    Read_knownIsoform_metagene.append(knownisoform)
                    Read_knownIsoform_metagene_scores.append(knownisoformscores)
            unique_ind = list(set(Ind))
            # logging genes without any reads
            log_ind = [ind for ind in range(len(Info_multigenes)) if ind not in unique_ind]
            for index in log_ind:
                for sample_ind in range(self.nsamples):
                    sample_target = self.target[sample_ind]
                    save_compatibleVector_by_gene(geneName=Info_multigenes[index][0]['geneName'],
                                              geneID=Info_multigenes[index][0]['geneID'],
                                              geneChr=Info_multigenes[index][0]['geneChr'],
                                              colNames=None,Read_Isoform_compatibleVector=None, #set this to None for log
                                              Read_knownIsoform_scores = None,
                                              qname_cbumi_dict=None, exonInfo=None,isoformInfo=None,
                                              output_folder=sample_target, parse=self.parse)
            #save compatible matrix by genes
            return_list = []
            for index in unique_ind:
                print('processing gene' + str(index))
                # loop over genes within metagene; for one single gene:
                Read_novelIsoform, Read_knownIsoform, Read_knownIsoform_scores, novel_isoformInfo = [], [], {}, {}
                for j, i in enumerate(Ind):#i: gene index; j: index of index---# loop for reads
                    if i == index and Read_novelIsoform_metagene[j] is not None:
                        Read_novelIsoform.append(Read_novelIsoform_metagene[j])
                    if i == index and Read_knownIsoform_metagene[j] is not None:
                        Read_knownIsoform.append(Read_knownIsoform_metagene[j])
                        Read_knownIsoform_scores[Read_knownIsoform_metagene[j][0]] = Read_knownIsoform_metagene_scores[j]
                if len(Read_novelIsoform) > 0:
                    Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished = polish_compatible_vectors(
                        Read_novelIsoform, Read_knownIsoform, len(Info_multigenes[index][2]),
                    Info_multigenes[index][1], self.small_exon_threshold, self.small_exon_threshold1)
                else:
                    Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished = (
                    Read_novelIsoform, novel_isoformInfo, Read_knownIsoform)
                geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(
                    Read_novelIsoform_polished, Read_knownIsoform_polished, Info_multigenes[index][0])
                for key in Read_knownIsoform_scores:
                    Read_knownIsoform_scores[key].extend([-1] * len(list(novel_isoformInfo_polished.keys())))
                # update annotation information in self
                self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames'] = \
                    self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames'] + list(
                        novel_isoformInfo_polished.keys())
                self.metageneStructureInformationwNovel[meta_gene][index][0]['numofIsoforms'] = \
                    self.metageneStructureInformationwNovel[meta_gene][index][0]['numofIsoforms'] + len(list(
                        novel_isoformInfo_polished.keys()))
                self.metageneStructureInformationwNovel[meta_gene][index][2].update(novel_isoformInfo_polished)
                for sample_ind in range(self.nsamples):
                    sample_target = self.target[sample_ind]
                    sample = 'sample'+str(sample_ind)
                    Read_Isoform_compatibleVector_sample, Read_knownIsoform_scores_sample = {}, {}
                    for readname in list(Read_Isoform_compatibleVector.keys()):
                        if qname_sample_dict[readname] == sample:
                            Read_Isoform_compatibleVector_sample[readname] = Read_Isoform_compatibleVector[readname]
                            if readname in Read_knownIsoform_scores.keys():
                                Read_knownIsoform_scores_sample[readname] = Read_knownIsoform_scores[readname]
                    if save:
                        save_compatibleVector_by_gene(geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector_sample,
                                                      Read_knownIsoform_scores_sample,self.qname_cbumi_dict_list[sample_ind],
                                                      self.metageneStructureInformationwNovel[meta_gene][index][1],
                                                      self.metageneStructureInformationwNovel[meta_gene][index][2],
                                                      sample_target, self.parse)
                    else:
                        return_list.append({'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector_sample,
                                            'isoforms': colNames,
                                'exonInfo': self.metageneStructureInformationwNovel[meta_gene][index][1],
                                'isoformInfo': self.metageneStructureInformationwNovel[meta_gene][index][2]})
            if save==False:
                return return_list
    def map_reads_parse(self, meta_gene, save = True):
        Info_multigenes = copy.deepcopy(self.metageneStructureInformation[meta_gene])
        Info_multigenes = sort_multigeneInfo(Info_multigenes)
        bamFilePysam = self.read_bam()
        if len(Info_multigenes)==1:
            Info_singlegene = Info_multigenes[0]
            geneInfo, exonInfo, isoformInfo = Info_singlegene
            n_isoforms = len(isoformInfo)
            reads = bamFilePysam.fetch(geneInfo['geneChr'], geneInfo['geneStart'], geneInfo['geneEnd'])
            Read_novelIsoform = [] #[('read name',[read-exon percentage],[read-exon mapping])]
            Read_knownIsoform = [] #[('read name',[read-isoform mapping])]
            Read_knownIsoform_scores = {}
            novel_isoformInfo = {} #{'novelIsoform_1234':[2,3,4]}
            samples_list = []
            for read in reads:
                poly, _ = detect_poly_parse(read, window=20, n=15)
                result = process_read(read, self.qname_dict, self.lowest_match,self.lowest_match1, self.small_exon_threshold,self.small_exon_threshold1,
                                      self.truncation_match, Info_singlegene, self.parse, self.pacbio)
                result_novel, result_known, result_known_scores = result
                samples_list.append(self.qname_sample_dict[read.qname])
                if result_novel is not None:
                    Read_novelIsoform.append(result_novel)
                if result_known is not None:
                    Read_knownIsoform.append(result_known)
                    Read_knownIsoform_scores[result_known[0]] = result_known_scores
            unique_samples = list(set(samples_list))
            return_samples = []
            if len(Read_novelIsoform) > 0:
                # novel_isoformInfo_polished: novel isoform annotation: {'novelIsoform_7':[0,1,2]}
                Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished = polish_compatible_vectors(
                    Read_novelIsoform, Read_knownIsoform, n_isoforms, exonInfo, self.small_exon_threshold,
                    self.small_exon_threshold1)
            else:
                Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished = (
                Read_novelIsoform, novel_isoformInfo,
                Read_knownIsoform)
            # compile output into compatible matrix
            geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(
                Read_novelIsoform_polished, Read_knownIsoform_polished, geneInfo)
            # update annotation information in self
            self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'] = \
                self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'] + list(
                    novel_isoformInfo_polished.keys())
            self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] = \
                self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] + len(list(
                    novel_isoformInfo_polished.keys()))
            self.metageneStructureInformationwNovel[meta_gene][0][2].update(novel_isoformInfo_polished)
            for key in Read_knownIsoform_scores:
                Read_knownIsoform_scores[key].extend([-1] * len(list(novel_isoformInfo_polished.keys())))
            for sample in unique_samples:
                sample_target = os.path.join(self.target[0], 'samples/' + sample)
                if not os.path.exists(sample_target):
                    os.makedirs(sample_target)
                Read_Isoform_compatibleVector_sample, Read_knownIsoform_scores_sample = {}, {}
                for readname in list(Read_Isoform_compatibleVector.keys()):
                    if self.qname_sample_dict[readname]==sample:
                        Read_Isoform_compatibleVector_sample[readname] = Read_Isoform_compatibleVector[readname]
                        if readname in Read_knownIsoform_scores.keys():
                            Read_knownIsoform_scores_sample[readname] = Read_knownIsoform_scores[readname]
                if save:
                    save_compatibleVector_by_gene(geneName, geneID, geneChr, colNames,
                                                  Read_Isoform_compatibleVector_sample, Read_knownIsoform_scores_sample,
                                                  self.qname_cbumi_dict,
                                                  self.metageneStructureInformationwNovel[meta_gene][0][1],
                                                  self.metageneStructureInformationwNovel[meta_gene][0][2],
                                                  sample_target, self.parse)
                else:
                    return_sample = {'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector_sample, 'isoforms': colNames,
                             'exonInfo': self.metageneStructureInformationwNovel[meta_gene][0][1],
                             'isoformInfo': self.metageneStructureInformationwNovel[meta_gene][0][2]}
                    return_samples.append(return_sample)
            if save==False:
                return return_samples
        else:
            geneChr, start, end = summarise_metagene(Info_multigenes)  # geneChr, start, end
            reads = bamFilePysam.fetch(geneChr, start, end)  # fetch reads within meta gene region
            # process reads metagene
            results, samples, polies = [], [], []
            for read in reads:
                poly, _ = detect_poly_parse(read, window=20, n=15)
                out = process_read_metagene(read, self.qname_dict, Info_multigenes, self.lowest_match, self.lowest_match1,self.small_exon_threshold,self.small_exon_threshold1,
                                            self.truncation_match, self.parse, self.pacbio)
                if out is not None: #may not within this meta gene region
                    polies.append(poly)
                    results.append(out)
                    samples.append(self.qname_sample_dict[read.qname])
            unique_samples = list(set(samples))
            Ind, Read_novelIsoform_metagene, Read_knownIsoform_metagene, Read_knownIsoform_metagene_scores = [], [], [], []
            for result in results:
                if result is not None:
                    ind, novelisoform, knownisoform, knownisoformscores = result
                    Ind.append(ind)
                    Read_novelIsoform_metagene.append(novelisoform)
                    Read_knownIsoform_metagene.append(knownisoform)
                    Read_knownIsoform_metagene_scores.append(knownisoformscores)
            unique_ind = list(set(Ind))
            # logging genes without any reads
            log_ind = [ind for ind in range(len(Info_multigenes)) if ind not in unique_ind]
            for index in log_ind:
                for sample in unique_samples:
                    sample_target = os.path.join(self.target[0], 'samples/' + sample)
                    if not os.path.exists(sample_target):
                        os.makedirs(sample_target)
                    save_compatibleVector_by_gene(geneName=Info_multigenes[index][0]['geneName'],
                                                  geneID=Info_multigenes[index][0]['geneID'],
                                                  geneChr=Info_multigenes[index][0]['geneChr'],
                                                  colNames=None, Read_Isoform_compatibleVector=None,
                                                  # set this to None for log
                                                  Read_knownIsoform_scores=None,
                                                  qname_cbumi_dict=None, exonInfo=None, isoformInfo=None,
                                                  output_folder=sample_target, parse=self.parse)
            return_samples = []
            for index in unique_ind: #loop gene
                print('processing gene' + str(index))
                # loop over genes within metagene; for one single gene:
                Read_novelIsoform, Read_knownIsoform, Read_knownIsoform_scores, novel_isoformInfo = [], [], {}, {}
                for j, i in enumerate(Ind):  # i: gene index; j: index of index---# loop for reads
                    if i == index and Read_novelIsoform_metagene[j] is not None:
                        Read_novelIsoform.append(Read_novelIsoform_metagene[j])
                    if i == index and Read_knownIsoform_metagene[j] is not None:
                        Read_knownIsoform.append(Read_knownIsoform_metagene[j])
                        Read_knownIsoform_scores[Read_knownIsoform_metagene[j][0]] = Read_knownIsoform_metagene_scores[
                            j]
                if len(Read_novelIsoform) > 0:
                    Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished = polish_compatible_vectors(
                        Read_novelIsoform, Read_knownIsoform, len(Info_multigenes[index][2]),
                        Info_multigenes[index][1], self.small_exon_threshold, self.small_exon_threshold1)
                else:
                    Read_novelIsoform_polished, novel_isoformInfo_polished, Read_knownIsoform_polished = (
                        Read_novelIsoform, novel_isoformInfo, Read_knownIsoform)
                geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(
                    Read_novelIsoform_polished, Read_knownIsoform_polished, Info_multigenes[index][0])
                for key in Read_knownIsoform_scores:
                    Read_knownIsoform_scores[key].extend([-1] * len(list(novel_isoformInfo_polished.keys())))
                # update annotation information in self
                self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames'] = \
                    self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames'] + list(
                        novel_isoformInfo_polished.keys())
                self.metageneStructureInformationwNovel[meta_gene][index][0]['numofIsoforms'] = \
                    self.metageneStructureInformationwNovel[meta_gene][index][0]['numofIsoforms'] + len(list(
                        novel_isoformInfo_polished.keys()))
                self.metageneStructureInformationwNovel[meta_gene][index][2].update(novel_isoformInfo_polished)
                if save:
                    for sample in unique_samples:
                        sample_target = os.path.join(self.target[0], 'samples/' + sample)
                        if not os.path.exists(sample_target):
                            os.makedirs(sample_target)
                        Read_Isoform_compatibleVector_sample, Read_knownIsoform_scores_sample = {}, {}
                        for readname in list(Read_Isoform_compatibleVector.keys()):
                            if self.qname_sample_dict[readname] == sample:
                                Read_Isoform_compatibleVector_sample[readname] = Read_Isoform_compatibleVector[readname]
                                if readname in Read_knownIsoform_scores.keys():
                                    Read_knownIsoform_scores_sample[readname] = Read_knownIsoform_scores[readname]
                        save_compatibleVector_by_gene(geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector_sample,
                                                  Read_knownIsoform_scores_sample, self.qname_cbumi_dict,
                                                  self.metageneStructureInformationwNovel[meta_gene][index][1],
                                                  self.metageneStructureInformationwNovel[meta_gene][index][2],
                                                  sample_target, self.parse)
                else:
                    for sample in unique_samples:
                        Read_Isoform_compatibleVector_sample, Read_knownIsoform_scores_sample = {}, {}
                        for readname in list(Read_Isoform_compatibleVector.keys()):
                            if self.qname_sample_dict[readname] == sample:
                                Read_Isoform_compatibleVector_sample[readname] = Read_Isoform_compatibleVector[readname]
                                if readname in Read_knownIsoform_scores.keys():
                                    Read_knownIsoform_scores_sample[readname] = Read_knownIsoform_scores[readname]
                        return_samples.append(
                            {'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector_sample, 'isoforms': colNames,
                             'exonInfo': self.metageneStructureInformationwNovel[meta_gene][index][1],
                             'isoformInfo': self.metageneStructureInformationwNovel[meta_gene][index][2]})
            if save == False:
                return return_samples
    def map_reads_allgenes(self, cover_existing = True, total_jobs = 1, current_job_index = 0):
        if self.parse==False:
            for compatible_matrix_folder_path in self.compatible_matrix_folder_path_list:
                if not os.path.exists(compatible_matrix_folder_path):
                    os.makedirs(compatible_matrix_folder_path, exist_ok=True)
        MetaGenes = list(self.metageneStructureInformation.keys()) #all meta genes
        if total_jobs > 1:
            step_size = math.ceil(len(MetaGenes) / total_jobs)
            s = int(list(range(0, len(MetaGenes), step_size))[current_job_index])
            e = int(s + step_size)
            MetaGenes_job = MetaGenes[s:e]
        else:#total_jobs = 1
            MetaGenes_job = MetaGenes
        if cover_existing:
            print('If there are existing compatible matrix files, SCOTCH will overwrite them')
            genes_existing = []
        else:
            print('If there are existing compatible matrix files, SCOTCH will not overwrite them')
            if self.parse:
                self.compatible_matrix_folder_paths = find_subfolder(self.target[0], subfolder='compatible_matrix')
                self.read_mapping_paths = find_subfolder(self.target[0], subfolder='auxillary')
                genes_existing = [file[:-4] for folder_path in self.compatible_matrix_folder_paths
                                  for file in os.listdir(folder_path) if file.endswith('.csv')]
                for folder_path in self.compatible_matrix_folder_paths:
                    log_file_path = os.path.join(folder_path, 'log.txt')
                    if os.path.isfile(log_file_path):
                        gene_df = pd.read_csv(log_file_path, header=None)
                        genes_existing += gene_df.iloc[:, 0].tolist()
            else:
                genes_existing = [file[:-4] for folder_path in self.compatible_matrix_folder_path_list
                                  for file in os.listdir(folder_path) if file.endswith('.csv')]
                for folder_path in self.compatible_matrix_folder_path_list:
                    log_file_path = os.path.join(folder_path, 'log.txt')
                    if os.path.isfile(log_file_path):
                        gene_df = pd.read_csv(log_file_path, header=None)
                        genes_existing += gene_df.iloc[:, 0].tolist()
            print('there exist ' + str(len(set(genes_existing))) + ' genes')
        MetaGene_Gene_dict = {}
        for metagene_name, genes_info in self.metageneStructureInformation.items():
            if metagene_name in MetaGenes_job:
                genes_ = []
                for gene_info in genes_info:
                    gene = str(gene_info[0]['geneName']) + '_' + str(gene_info[0]['geneID'])
                    if gene not in genes_existing:
                        genes_.append(gene)
                if len(genes_) > 0:
                    MetaGene_Gene_dict[metagene_name] = genes_
        MetaGenes_job = list(MetaGene_Gene_dict.keys())
        self.logger.info(f'{str(len(MetaGenes_job))} metagenes for job {current_job_index}')
        #print('processing ' + str(len(MetaGenes_job)) + ' metagenes for this job')
        if self.parse:
            for meta_gene in MetaGenes_job:
                print(meta_gene)
                self.map_reads_parse(meta_gene, save=True)
        else:
            for meta_gene in MetaGenes_job:
                print(meta_gene)
                self.map_reads(meta_gene, save=True)
        for key in MetaGenes:
            if key not in MetaGenes_job:
                del self.metageneStructureInformationwNovel[key]
        gene_ids = [g_name_id.split('_')[1] for g_name_ids in list(MetaGene_Gene_dict.values()) for g_name_id in g_name_ids]
        gene_ids_pattern = '|'.join([f'gene_id "{gene_id}"' for gene_id in gene_ids])
        self.gtf_df_job = self.gtf_df[self.gtf_df['attribute'].str.contains(gene_ids_pattern, regex=True)].reset_index(drop=True)
    def save_annotation_w_novel_isoform(self, total_jobs = 1, current_job_index = 0):
        self.logger.info(f'Saving annotation file...')
        for i in range(len(self.annotation_path_meta_gene_list)):
            if total_jobs>1:
                file_name_pkl = self.annotation_path_meta_gene_novel_list[i][:-4] + '_' + str(current_job_index) +'.pkl'
                file_name_gtf = self.annotation_path_gtf_novel_list[i][:-4] + '_' + str(current_job_index) +'.gtf'
            else:
                file_name_pkl = self.annotation_path_meta_gene_novel_list[i]
                file_name_gtf = self.annotation_path_gtf_novel_list[i]
            #save pickle file
            with open(file_name_pkl, 'wb') as file:
                pickle.dump(self.metageneStructureInformationwNovel, file)
            #save gtf file
            convert_to_gtf(self.metageneStructureInformationwNovel, file_name_gtf, self.gtf_df_job, num_cores=1)




