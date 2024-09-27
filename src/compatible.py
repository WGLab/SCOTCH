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
            new_row_gene = [[geneInfo['geneChr'], 'SCOTCH', 'gene', geneInfo['geneStart'], geneInfo['geneEnd'], '.',
                            geneInfo['geneStrand'], '.', f'gene_id "{geneID}"; gene_name "{geneInfo["geneName"]}"']]
        new_rows_isoforms = []
        for isoform_name, exon_indices in isoformInfo.items():
            isoform_start = exonInfo[exon_indices[0]][0]
            isoform_end = exonInfo[exon_indices[-1]][1]
            if gtf_df_sub[gtf_df_sub['attribute'].str.contains(f'transcript_id "{isoform_name}"', regex=False)].shape[0] == 0:
                new_rows_isoform = [geneInfo['geneChr'], 'SCOTCH', 'transcript', isoform_start, isoform_end, '.',
                                    geneInfo['geneStrand'], '.',
                                    f'gene_id "{geneID}"; gene_name "{geneInfo["geneName"]}"; transcript_id "{isoform_name}"; transcript_name "{isoform_name}"']
                new_rows_isoforms.append(new_rows_isoform)
                exons_isoform = []
                for exon_index in exon_indices:
                    exons_isoform.append(exonInfo[exon_index])
                merged_exons_isoform = merge_exons(exons_isoform)
                for exon_num, (exon_start, exon_end) in enumerate(merged_exons_isoform, start=1):
                    new_rows_exon = [geneInfo['geneChr'], 'SCOTCH', 'exon', exon_start, exon_end, '.',
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
        output_gtf = os.path.join(reference_folder, "metageneStructureInformationwNovel.gtf")
        pattern_pkl = re.compile(r".*_\d+\.pkl$")
        pattern_gtf = re.compile(r".*_\d+\.gtf$")
        file_names_pkl = [os.path.join(target, 'reference', f) for f in os.listdir(os.path.join(target, "reference")) if pattern_pkl.match(f)]
        file_names_gtf = [os.path.join(target, 'reference', f) for f in os.listdir(os.path.join(target, "reference")) if
                          pattern_gtf.match(f)]
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
            #for file_name_pkl in file_names_pkl:
            #    os.remove(file_name_pkl)
            print('mergered new isoform annotation saved at: '+str(file_names_pkl))
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
            #for file_name_gtf in file_names_gtf:
            #    os.remove(file_name_gtf)
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
        file_paths = [os.path.join(auxillary_folder, f) for f in os.listdir(auxillary_folder)]
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
        DF['priority'] = DF['MappingScore']*DF['priority']
        DF['GeneMapping'] = 'delete'  # Initialize as 'delete'
        DF = DF.sort_values(by=['geneChr', 'Read', 'priority'], ascending=[True, True, False])
        grouped = DF.groupby(['Read'], group_keys=False)
        # Process groups in parallel and concatenate results
        processed_groups = Parallel(n_jobs=-1)(delayed(process_group)(group) for _, group in grouped)
        DF = pd.concat(processed_groups).reset_index(drop=True)
        DF.to_csv(os.path.join(auxillary_folder, 'all_read_isoform_exon_mapping.tsv'), sep='\t', index=False)






class ReadMapper:
    def __init__(self, target, bam_path, lowest_match=0.2, small_exon_threshold = 20, small_exon_threshold1=100, truncation_match=0.5, platform = '10x',
                 reference_gtf_path = None):
        self.target = target
        self.bam_path = bam_path
        self.small_exon_threshold = small_exon_threshold
        self.small_exon_threshold1 = small_exon_threshold1
        self.truncation_match = truncation_match
        column_names = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        if reference_gtf_path is not None or reference_gtf_path=='None':
            self.gtf_df = pd.read_csv(reference_gtf_path, sep='\t', comment='#', header=None, names=column_names)
        else:
            self.gtf_df = pd.DataFrame(columns=column_names)
        # gene annotation information
        self.annotation_folder_path = os.path.join(target, "reference")
        self.annotation_path_single_gene = os.path.join(target, "reference/geneStructureInformation.pkl")
        self.annotation_path_meta_gene = os.path.join(target, "reference/metageneStructureInformation.pkl")
        self.annotation_path_meta_gene_novel = os.path.join(target, "reference/metageneStructureInformationwNovel.pkl")
        self.annotation_path_gtf_novel = os.path.join(target, "reference/gene_annotations_scotch.gtf")
        # bam information path
        self.bamInfo_folder_path = os.path.join(target, "bam")
        self.bamInfo_pkl_path = os.path.join(target, 'bam/bam.Info.pkl')#bamInfo_pkl_file
        self.bamInfo2_pkl_path = os.path.join(target, 'bam/bam.Info2.pkl')#bamInfo2_pkl_file
        self.bamInfo3_pkl_path = os.path.join(target, 'bam/bam.Info3.pkl') # bamInfo2_pkl_file
        self.bamInfo_csv_path = os.path.join(target, 'bam/bam.Info.csv')
        # parameters
        self.lowest_match = lowest_match
        self.platform = platform
        self.parse = self.platform == 'parse'
        self.pacbio = self.platform == 'pacbio'
        # some paths
        if platform != 'parse':
            self.compatible_matrix_folder_path = os.path.join(target, "compatible_matrix")#not for parse
            self.read_mapping_path = os.path.join(target, "auxillary")#not for parse
            #self.count_matrix_folder_path = os.path.join(target, "count_matrix")#not for parse
        # bam information file
        self.qname_dict = load_pickle(self.bamInfo_pkl_path)
        self.qname_cbumi_dict = load_pickle(self.bamInfo2_pkl_path)
        self.sorted_bam_path = None
        self.qname_sample_dict = load_pickle(self.bamInfo3_pkl_path)
        self.metageneStructureInformation = load_pickle(self.annotation_path_meta_gene)
        self.metageneStructureInformationwNovel = self.metageneStructureInformation.copy()
    def merge_bam(self):
        merged_folder = os.path.join(self.bam_path, 'merged')
        os.makedirs(merged_folder, exist_ok=True)
        bamFile_name = [f for f in os.listdir(self.bam_path) if
                        f.endswith('.bam') and f + '.bai' in os.listdir(self.bam_path)]
        print('merging bam files, usually sublibraries')
        self.sorted_bam_path = os.path.join(merged_folder, 'merged.sorted.bam')
        pysam.merge(os.path.join(merged_folder, 'merged.bam'), bamFile_name)
        print('sorting bam files, usually sublibraries')
        pysam.sort("-o", self.sorted_bam_path, os.path.join(merged_folder, 'merged.bam'))
        pysam.index(self.sorted_bam_path)
    def read_bam(self, chrom = None):
        # parse: input a folder, find merged bam file to read; read the single bam file if input is a file path
        # bam_path is a folder
        if os.path.isfile(self.bam_path) == False:
            # find the bam file based on chrom
            if chrom is not None:
                bamFile_name = [f for f in os.listdir(self.bam_path) if
                            f.endswith('.bam') and '.' + chrom + '.' in f]
                bamFile = os.path.join(self.bam_path, bamFile_name[0]) #not bai
                bamFilePysam = pysam.Samfile(bamFile, "rb")
            else:
                #read the merged bam file, has to run merge_bam first
                bamFilePysam = pysam.Samfile(self.sorted_bam_path, "rb")
        else:
            bamFilePysam = pysam.Samfile(self.bam_path, "rb")
        return bamFilePysam
    def map_reads(self, meta_gene, save = True):
        #used for ont and pacbio
        Info_multigenes = copy.deepcopy(self.metageneStructureInformation[meta_gene])
        Info_multigenes = sort_multigeneInfo(Info_multigenes)
        bamFilePysam = self.read_bam(chrom=Info_multigenes[0][0]['geneChr'])
        if len(Info_multigenes)==1:
            Info_singlegene = Info_multigenes[0]
            geneInfo, exonInfo, isoformInfo = Info_singlegene
            n_isoforms = len(isoformInfo)
            reads = bamFilePysam.fetch(geneInfo['geneChr'], geneInfo['geneStart'], geneInfo['geneEnd'])
            Read_novelIsoform = [] #[('read name',[read-exon percentage],[read-exon mapping])]
            Read_knownIsoform = [] #[('read name',[read-isoform mapping])]
            Read_knownIsoform_scores = {}#[readname: read-isoform mapping scores]
            novel_isoformInfo = {} #{'novelIsoform_1234':[2,3,4]}
            for read in reads:
                result = process_read(read, self.qname_dict, self.lowest_match,self.small_exon_threshold,self.small_exon_threshold1,
                                      self.truncation_match,Info_singlegene, self.parse, self.pacbio)
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
            if save:
                # save compatible matrix of each gene, save read-isoform mappings
                save_compatibleVector_by_gene(geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector, Read_knownIsoform_scores,
                                              self.qname_cbumi_dict, self.metageneStructureInformationwNovel[meta_gene][0][1],
                                              self.metageneStructureInformationwNovel[meta_gene][0][2], self.target,
                                              self.parse)
            else:
                return [{'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector, 'isoforms': colNames,
                 'exonInfo': self.metageneStructureInformationwNovel[meta_gene][0][1],
                'isoformInfo':self.metageneStructureInformationwNovel[meta_gene][0][2]}]
        else:
            geneChr, start, end = summarise_metagene(Info_multigenes)  # geneChr, start, end
            reads = bamFilePysam.fetch(geneChr, start, end)  # fetch reads within meta gene region
            # process reads metagene
            results = []
            for read in reads:
                out = process_read_metagene(read,self.qname_dict, Info_multigenes, self.lowest_match,
                                            self.small_exon_threshold,self.small_exon_threshold1,
                                            self.truncation_match, self.parse, self.pacbio)
                if out is not None: #may not within this meta gene region
                    results.append(out)
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
                save_compatibleVector_by_gene(geneName=Info_multigenes[index][0]['geneName'],
                                              geneID=Info_multigenes[index][0]['geneID'],
                                              geneChr=Info_multigenes[index][0]['geneChr'],
                                              colNames=None,Read_Isoform_compatibleVector=None, #set this to None for log
                                              Read_knownIsoform_scores = None,
                                              qname_cbumi_dict=None, exonInfo=None,isoformInfo=None,
                                              output_folder=self.target, parse=self.parse)
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
                if save:
                    save_compatibleVector_by_gene(geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector,
                                                  Read_knownIsoform_scores,self.qname_cbumi_dict,
                                                  self.metageneStructureInformationwNovel[meta_gene][index][1],
                                                  self.metageneStructureInformationwNovel[meta_gene][index][2],
                                                  self.target, self.parse)
                else:
                    return_list.append({'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector, 'isoforms': colNames,
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
            samples_novel, samples_known = [], []
            Read_novelIsoform_poly = []
            for read in reads:
                poly, _ = detect_poly_parse(read, window=20, n=15)
                result = process_read(read, self.qname_dict, self.lowest_match,self.small_exon_threshold,self.small_exon_threshold1,
                                      self.truncation_match, Info_singlegene, self.parse, self.pacbio)
                result_novel, result_known, result_known_scores = result
                if result_novel is not None:
                    Read_novelIsoform.append(result_novel)
                    samples_novel.append(self.qname_sample_dict[read.qname])
                    Read_novelIsoform_poly.append(poly)
                if result_known is not None:
                    Read_knownIsoform.append(result_known)
                    samples_known.append(self.qname_sample_dict[read.qname])
                    Read_knownIsoform_scores[result_known[0]] = result_known_scores
            unique_samples = list(set(samples_novel+samples_known))
            return_samples = []
            for sample in unique_samples:
                Read_novelIsoform_sample, Read_knownIsoform_sample, Read_knownIsoform_scores_sample, Read_novelIsoform_poly_sample = [],[], {}, []
                sample_target = os.path.join(self.target, 'samples/'+sample)
                if not os.path.exists(sample_target):
                    os.makedirs(sample_target)
                sample_index_novel = [i for i, s in enumerate(samples_novel) if s == sample]
                sample_index_known = [i for i, s in enumerate(samples_known) if s == sample]
                if len(sample_index_novel) > 0:
                    Read_novelIsoform_sample = [Read_novelIsoform[i] for i in sample_index_novel]
                    #Read_novelIsoform_poly_sample = [Read_novelIsoform_poly[i] for i in sample_index_novel]
                if len(sample_index_known) > 0:
                    Read_knownIsoform_sample = [Read_knownIsoform[i] for i in sample_index_known]
                    Read_knownIsoform_scores_sample = {Read_knownIsoform[i][0]:Read_knownIsoform_scores[Read_knownIsoform[i][0]] for i in sample_index_known}
                if len(Read_novelIsoform_sample) > 0:
                    Read_novelIsoform_sample_polished, novel_isoformInfo_polished, Read_knownIsoform_sample_polished = polish_compatible_vectors(
                        Read_novelIsoform_sample,Read_knownIsoform_sample, n_isoforms, exonInfo, self.small_exon_threshold, self.small_exon_threshold1)
                else:
                    Read_novelIsoform_sample_polished, novel_isoformInfo_polished, Read_knownIsoform_sample_polished = (
                    Read_novelIsoform_sample, novel_isoformInfo, Read_knownIsoform_sample)
                geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector_sample = compile_compatible_vectors(
                    Read_novelIsoform_sample_polished, Read_knownIsoform_sample_polished, geneInfo)
                # update annotation information in self
                for novel_isoform_name in list(novel_isoformInfo_polished.keys()):
                    if novel_isoform_name not in self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames']:
                        self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'].append(novel_isoform_name)
                self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] = len(self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'])
                self.metageneStructureInformationwNovel[meta_gene][0][2].update(novel_isoformInfo_polished)
                n = self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms']-n_isoforms
                for key in Read_knownIsoform_scores_sample:
                    Read_knownIsoform_scores_sample[key].extend([-1] * n)
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
                out = process_read_metagene(read, self.qname_dict, Info_multigenes, self.lowest_match, self.small_exon_threshold,self.small_exon_threshold1,
                                            self.truncation_match, self.parse, self.pacbio)
                if out is not None: #may not within this meta gene region
                    polies.append(poly)
                    results.append(out)
                    samples.append(self.qname_sample_dict[read.qname])
            unique_samples = list(set(samples))
            return_samples = []
            for sample in unique_samples:
                sample_target = os.path.join(self.target, 'samples/'+sample)
                if not os.path.exists(sample_target):
                    os.makedirs(sample_target)
                sample_index = [i for i, s in enumerate(samples) if s == sample]
                result_sample = [results[i] for i in sample_index]
                poly_sample = [polies[i] for i in sample_index]
                Ind, Read_novelIsoform_metagene, Read_knownIsoform_metagene, Read_knownIsoform_metagene_scores, poly_sample_filtered = [], [], [], [], []
                for ri, result in enumerate(result_sample):
                    if result is not None:
                        ind, novelisoform, knownisoform, knownisoformscores = result
                        Ind.append(ind)
                        Read_novelIsoform_metagene.append(novelisoform)
                        Read_knownIsoform_metagene.append(knownisoform)
                        Read_knownIsoform_metagene_scores.append(knownisoformscores)
                        poly_sample_filtered.append(poly_sample[ri])
                unique_ind = list(set(Ind))
                # logging genes without any reads
                log_ind = [ind for ind in range(len(Info_multigenes)) if ind not in unique_ind]
                for index in log_ind:
                    save_compatibleVector_by_gene(geneName=Info_multigenes[index][0]['geneName'],
                                                  geneID=Info_multigenes[index][0]['geneID'],
                                                  geneChr=Info_multigenes[index][0]['geneChr'],
                                                  colNames=None,Read_Isoform_compatibleVector=None, #set this to None for log
                                                  Read_knownIsoform_scores = None,
                                                  qname_cbumi_dict=None, exonInfo=None,isoformInfo=None,
                                                  output_folder=sample_target, parse = self.parse)
                #save compatible matrix by genes
                for index in unique_ind:
                    print('processing gene' + str(index))
                    # loop over genes within metagene; for one single gene:
                    Read_novelIsoform, Read_knownIsoform, Read_knownIsoform_scores, novel_isoformInfo, Read_novelIsoform_poly = [], [],{}, {}, []
                    for j, i in enumerate(Ind):#i: gene index; j: index of index---# loop for reads
                        if i == index and Read_novelIsoform_metagene[j] is not None:
                            Read_novelIsoform.append(Read_novelIsoform_metagene[j])
                            Read_novelIsoform_poly.append(poly_sample_filtered[j])
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
                        Read_novelIsoform_polished,Read_knownIsoform_polished, Info_multigenes[index][0])
                    # update annotation information in self
                    for novel_isoform_name in list(novel_isoformInfo_polished.keys()):
                        if novel_isoform_name not in self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames']:
                            self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames'].append(novel_isoform_name)
                    self.metageneStructureInformationwNovel[meta_gene][index][0]['numofIsoforms'] = len(self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames'])
                    self.metageneStructureInformationwNovel[meta_gene][index][2].update(novel_isoformInfo_polished)
                    n = self.metageneStructureInformationwNovel[meta_gene][index][0]['numofIsoforms'] - Info_multigenes[index][0]['numofIsoforms']
                    for key in Read_knownIsoform_scores:
                        Read_knownIsoform_scores[key].extend([-1] * n)
                    if save:
                        save_compatibleVector_by_gene(geneName, geneID, geneChr, colNames, Read_Isoform_compatibleVector,Read_knownIsoform_scores,
                                                      qname_cbumi_dict = self.qname_cbumi_dict,
                                                      exonInfo = self.metageneStructureInformationwNovel[meta_gene][index][1],
                                                      isoformInfo=self.metageneStructureInformationwNovel[meta_gene][index][2],
                                                      output_folder = sample_target, parse=self.parse)
                    else:
                        return_samples.append({'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector, 'isoforms': colNames,
                                'exonInfo': self.metageneStructureInformationwNovel[meta_gene][index][1],
                                'isoformInfo': self.metageneStructureInformationwNovel[meta_gene][index][2]})
                if save==False:
                    return return_samples
    def map_reads_allgenes(self, cover_existing = True, total_jobs = 1, current_job_index = 0):
        if self.parse==False:
            if not os.path.exists(self.compatible_matrix_folder_path):
                os.makedirs(self.compatible_matrix_folder_path)
        MetaGenes = list(self.metageneStructureInformation.keys()) #all meta genes
        if total_jobs > 1:
            step_size = math.ceil(len(MetaGenes) / total_jobs)
            s = int(list(range(0, len(MetaGenes), step_size))[current_job_index])
            e = int(s + step_size)
            MetaGenes_job = MetaGenes[s:e]
        else:#total_jobs = 1
            MetaGenes_job = MetaGenes
        print(str(len(MetaGenes_job)) + ' metagenes for this job')
        if cover_existing:
            print('If there are existing compatible matrix files, SCOTCH will overwrite them')
            genes_existing = []
        else:
            print('If there are existing compatible matrix files, SCOTCH will not overwrite them')
            if self.parse:
                self.compatible_matrix_folder_paths = find_subfolder(self.target, subfolder='compatible_matrix')
                self.read_mapping_paths = find_subfolder(self.target, subfolder='auxillary')
                genes_existing = [file[:-4] for folder_path in self.compatible_matrix_folder_paths
                                  for file in os.listdir(folder_path) if file.endswith('.csv')]
                for folder_path in self.compatible_matrix_folder_paths:
                    log_file_path = os.path.join(folder_path, 'log.txt')
                    if os.path.isfile(log_file_path):
                        gene_df = pd.read_csv(log_file_path, header=None)
                        genes_existing += gene_df.iloc[:, 0].tolist()
            else:
                genes_existing = [g[:-4] for g in os.listdir(self.compatible_matrix_folder_path)]
                if os.path.isfile(os.path.join(self.compatible_matrix_folder_path, 'log.txt')):
                    gene_df = pd.read_csv(os.path.join(self.compatible_matrix_folder_path, 'log.txt'), header=None)
                    genes_existing = genes_existing + gene_df.iloc[:, 0].tolist()
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
        print('processing ' + str(len(MetaGenes_job)) + ' metagenes for this job')
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
        if total_jobs>1:
            file_name_pkl = self.annotation_path_meta_gene_novel[:-4] + '_' + str(current_job_index) +'.pkl'
            file_name_gtf = self.annotation_path_gtf_novel[:-4] + '_' + str(current_job_index) +'.gtf'
        else:
            file_name_pkl = self.annotation_path_meta_gene_novel
            file_name_gtf = self.annotation_path_gtf_novel
        #save pickle file
        with open(file_name_pkl, 'wb') as file:
            pickle.dump(self.metageneStructureInformationwNovel, file)
        #save gtf file
        convert_to_gtf(self.metageneStructureInformationwNovel, file_name_gtf, self.gtf_df_job, num_cores=1)








