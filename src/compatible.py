from preprocessing import *
import pysam


class ReadMapper:
    def __init__(self, target, bam_path, lowest_match=0.2, platform = '10x'):
        self.target = target
        self.bam_path = bam_path
        # gene annotation information
        self.annotation_folder_path = os.path.join(target, "reference")
        self.annotation_path_single_gene = os.path.join(target, "reference/geneStructureInformation.pkl")
        self.annotation_path_meta_gene = os.path.join(target, "reference/metageneStructureInformation.pkl")
        # bam information path
        self.bamInfo_folder_path = os.path.join(target, "bam")
        self.bamInfo_pkl_path = os.path.join(target, 'bam/bam.Info.pkl')#bamInfo_pkl_file
        self.bamInfo2_pkl_path = os.path.join(target, 'bam/bam.Info2.pkl')#bamInfo2_pkl_file
        self.bamInfo3_pkl_path = os.path.join(target, 'bam/bam.Info3.pkl')  # bamInfo2_pkl_file
        self.bamInfo_csv_path = os.path.join(target, 'bam/bam.Info.csv')
        # parameters
        self.lowest_match = lowest_match
        self.platform = platform
        self.parse = self.platform == 'parse'
        self.pacbio = self.platform == 'pacbio'
        # some paths
        if platform!='parse':
            self.compatible_matrix_folder_path = os.path.join(target, "compatible_matrix")#not for parse
            self.read_mapping_path = os.path.join(target, "auxillary")#not for parse
            self.count_matrix_folder_path = os.path.join(target, "count_matrix")#not for parse
        # bam information file
        self.qname_dict = load_pickle(self.bamInfo_pkl_path)
        self.qname_cbumi_dict = load_pickle(self.bamInfo2_pkl_path)
        self.sorted_bam_path = None
        if self.parse:
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
        Info_multigenes = self.metageneStructureInformation[meta_gene]
        Info_multigenes = sort_multigeneInfo(Info_multigenes)
        bamFilePysam = self.read_bam(chrom=Info_multigenes[0][0]['geneChr'])
        if len(Info_multigenes)==1:
            geneInfo, exonInfo, isoformInfo = Info_multigenes[0]
            n_isoforms = len(isoformInfo)
            reads = bamFilePysam.fetch(geneInfo['geneChr'], geneInfo['geneStart'], geneInfo['geneEnd'])
            Read_novelIsoform = [] #[('read name',[read-exon percentage],[read-exon mapping])]
            Read_knownIsoform = [] #[('read name',[read-isoform mapping])]
            novel_isoformInfo = {} #{'novelIsoform_1234':[2,3,4]}
            for read in reads:
                result = process_read(read, geneInfo, exonInfo, isoformInfo, self.qname_dict, self.lowest_match,
                                      Info_multigenes, self.parse, self.pacbio)
                result_novel, result_known = result
                if result_novel is not None:
                    Read_novelIsoform.append(result_novel)
                if result_known is not None:
                    Read_knownIsoform.append(result_known)
            #expand uncategorized novel reads into Read_knownIsoform
            if len(Read_novelIsoform) > 0:
                Read_novelIsoform, novel_isoformInfo, Read_knownIsoform = polish_compatible_vectors(
                    Read_novelIsoform, Read_knownIsoform, n_isoforms)
            #compile output into compatible matrix
            geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(
                    Read_novelIsoform, novel_isoformInfo, Read_knownIsoform, geneInfo)
            #update annotation information in self
            self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'] = \
                self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames']+list(novel_isoformInfo.keys())
            self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] = \
                self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] + len(list(
                    novel_isoformInfo.keys()))
            self.metageneStructureInformationwNovel[meta_gene][0][2].update(novel_isoformInfo)
            if save:
                # save compatible matrix of each gene, save read-isoform mappings
                save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector,
                                      self.qname_cbumi_dict, self.metageneStructureInformationwNovel[meta_gene][0][1],
                                          self.metageneStructureInformationwNovel[meta_gene][0][2], self.target)
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
                out = process_read_metagene(read, start, end, self.qname_dict, Info_multigenes, self.lowest_match, self.parse, self.pacbio)
                if out is not None: #may not within this meta gene region
                    results.append(out)
            Ind, Read_novelIsoform_metagene, Read_knownIsoform_metagene = map(list, zip(*results))
            unique_ind = list(set(Ind))
            # logging genes without any reads
            log_ind = [ind for ind in range(len(Info_multigenes)) if ind not in unique_ind]
            for index in log_ind:
                save_compatibleVector_by_gene(geneName=Info_multigenes[index][0]['geneName'],
                                              geneID=Info_multigenes[index][0]['geneID'],
                                              geneStrand=Info_multigenes[index][0]['geneStrand'],
                                              colNames=None,Read_Isoform_compatibleVector=None, #set this to None for log
                                              qname_cbumi_dict=None, exonInfo=None,isoformInfo=None,
                                              output_folder=self.target)
            #save compatible matrix by genes
            return_list = []
            for index in unique_ind:
                print('processing gene' + str(index))
                # loop over genes within metagene; for one single gene:
                Read_novelIsoform, Read_knownIsoform, novel_isoformInfo = [], [], {}
                for j, i in enumerate(Ind):#i: gene index; j: index of index---# loop for reads
                    if i == index and Read_novelIsoform_metagene[j] is not None:
                        Read_novelIsoform.append(Read_novelIsoform_metagene[j])
                    if i == index and Read_knownIsoform_metagene[j] is not None:
                        Read_knownIsoform.append(Read_knownIsoform_metagene[j])
                if len(Read_novelIsoform) > 0:
                    Read_novelIsoform, novel_isoformInfo, Read_knownIsoform = polish_compatible_vectors(
                        Read_novelIsoform, Read_knownIsoform, len(Info_multigenes[index][2]))
                geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(
                    Read_novelIsoform, novel_isoformInfo, Read_knownIsoform, Info_multigenes[index][0])
                # update annotation information in self
                self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames'] = \
                    self.metageneStructureInformationwNovel[meta_gene][index][0]['isoformNames'] + list(
                        novel_isoformInfo.keys())
                self.metageneStructureInformationwNovel[meta_gene][index][0]['numofIsoforms'] = \
                    self.metageneStructureInformationwNovel[meta_gene][index][0]['numofIsoforms'] + len(list(
                        novel_isoformInfo.keys()))
                self.metageneStructureInformationwNovel[meta_gene][index][2].update(novel_isoformInfo)
                if save:
                    save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector,
                                                  self.qname_cbumi_dict,
                                                  self.metageneStructureInformationwNovel[meta_gene][index][1],
                                                  self.metageneStructureInformationwNovel[meta_gene][index][2],
                                                  self.target)
                else:
                    return_list.append({'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector, 'isoforms': colNames,
                            'exonInfo': self.metageneStructureInformationwNovel[meta_gene][index][1],
                            'isoformInfo': self.metageneStructureInformationwNovel[meta_gene][index][2]})
            if save==False:
                return return_list
    def map_reads_parse(self, meta_gene, save = True):
        Info_multigenes = self.metageneStructureInformation[meta_gene]
        Info_multigenes = sort_multigeneInfo(Info_multigenes)
        bamFilePysam = self.read_bam()
        if len(Info_multigenes)==1:
            geneInfo, exonInfo, isoformInfo = Info_multigenes[0]
            n_isoforms = len(isoformInfo)
            reads = bamFilePysam.fetch(geneInfo['geneChr'], geneInfo['geneStart'], geneInfo['geneEnd'])
            Read_novelIsoform = [] #[('read name',[read-exon percentage],[read-exon mapping])]
            Read_knownIsoform = [] #[('read name',[read-isoform mapping])]
            novel_isoformInfo = {} #{'novelIsoform_1234':[2,3,4]}
            samples_novel, samples_known = [], []
            for read in reads:
                result = process_read(read, geneInfo, exonInfo, isoformInfo, self.qname_dict, self.lowest_match,
                                      Info_multigenes, self.parse, self.pacbio)
                result_novel, result_known = result
                if result_novel is not None:
                    Read_novelIsoform.append(result_novel)
                    samples_novel.append(self.qname_sample_dict[read.qname])
                if result_known is not None:
                    Read_knownIsoform.append(result_known)
                    samples_known.append(self.qname_sample_dict[read.qname])
            unique_samples = list(set(samples_novel+samples_known))
            return_samples = []
            for sample in unique_samples:
                Read_novelIsoform_sample, Read_knownIsoform_sample = [], []
                sample_target = os.path.join(self.target, 'samples/'+sample)
                sample_index_novel = [i for i, s in enumerate(samples_novel) if s == sample]
                sample_index_known = [i for i, s in enumerate(samples_known) if s == sample]
                if len(sample_index_novel) > 0:
                    Read_novelIsoform_sample = [Read_novelIsoform[i] for i in sample_index_novel]
                if len(sample_index_known) > 0:
                    Read_knownIsoform_sample = [Read_knownIsoform[i] for i in sample_index_known]
                if len(Read_novelIsoform_sample) > 0:
                    Read_novelIsoform_sample, novel_isoformInfo, Read_knownIsoform_sample = polish_compatible_vectors(Read_novelIsoform_sample,Read_knownIsoform_sample, n_isoforms)
                geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector_sample = compile_compatible_vectors(Read_novelIsoform_sample, novel_isoformInfo,Read_knownIsoform_sample, geneInfo)
                # update annotation information in self
                for novel_isoform_name in list(novel_isoformInfo.keys()):
                    if novel_isoform_name not in self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames']:
                        self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'].append(novel_isoform_name)
                self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] = len(self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'])
                self.metageneStructureInformationwNovel[meta_gene][0][2].update(novel_isoformInfo)
                if save:
                    save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames,
                                                  Read_Isoform_compatibleVector_sample, self.qname_cbumi_dict,
                                                  self.metageneStructureInformationwNovel[meta_gene][0][1],
                                                  self.metageneStructureInformationwNovel[meta_gene][0][2],
                                                  sample_target)
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
            results, samples = [], []
            for read in reads:
                out = process_read_metagene(read, start, end, self.qname_dict, Info_multigenes, self.lowest_match, self.parse, self.pacbio)
                if out is not None: #may not within this meta gene region
                    results.append(out)
                    samples.append(self.qname_sample_dict[read.qname])
            unique_samples = list(set(samples))
            return_samples = []
            for sample in unique_samples:
                sample_target = os.path.join(self.target, 'samples/'+sample)
                sample_index = [i for i, s in enumerate(samples) if s == sample]
                result_sample = [results[i] for i in sample_index]
                Ind, Read_novelIsoform_metagene, Read_knownIsoform_metagene = [], [], []
                for result in result_sample:
                    if result is not None:
                        ind, novelisoform, knownisoform = result
                        Ind.append(ind)
                        Read_novelIsoform_metagene.append(novelisoform)
                        Read_knownIsoform_metagene.append(knownisoform)
                unique_ind = list(set(Ind))
                # logging genes without any reads
                log_ind = [ind for ind in range(len(Info_multigenes)) if ind not in unique_ind]
                for index in log_ind:
                    save_compatibleVector_by_gene(geneName=Info_multigenes[index][0]['geneName'],
                                                  geneID=Info_multigenes[index][0]['geneID'],
                                                  geneStrand=Info_multigenes[index][0]['geneStrand'],
                                                  colNames=None,Read_Isoform_compatibleVector=None, #set this to None for log
                                                  qname_cbumi_dict=None, exonInfo=None,isoformInfo=None,
                                                  output_folder=sample_target)
                #save compatible matrix by genes
                for index in unique_ind:
                    print('processing gene' + str(index))
                    # loop over genes within metagene; for one single gene:
                    Read_novelIsoform, Read_knownIsoform, novel_isoformInfo = [], [], {}
                    for j, i in enumerate(Ind):#i: gene index; j: index of index---# loop for reads
                        if i == index and Read_novelIsoform_metagene[j] is not None:
                            Read_novelIsoform.append(Read_novelIsoform_metagene[j])
                        if i == index and Read_knownIsoform_metagene[j] is not None:
                            Read_knownIsoform.append(Read_knownIsoform_metagene[j])
                    if len(Read_novelIsoform) > 0:
                        Read_novelIsoform, novel_isoformInfo, Read_knownIsoform = polish_compatible_vectors(
                            Read_novelIsoform, Read_knownIsoform, len(Info_multigenes[index][2]))
                    geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector = compile_compatible_vectors(
                        Read_novelIsoform, novel_isoformInfo, Read_knownIsoform, Info_multigenes[index][0])
                    # update annotation information in self
                    for novel_isoform_name in list(novel_isoformInfo.keys()):
                        if novel_isoform_name not in self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames']:
                            self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'].append(novel_isoform_name)
                    self.metageneStructureInformationwNovel[meta_gene][0][0]['numofIsoforms'] = len(self.metageneStructureInformationwNovel[meta_gene][0][0]['isoformNames'])
                    self.metageneStructureInformationwNovel[meta_gene][0][2].update(novel_isoformInfo)
                    if save:
                        save_compatibleVector_by_gene(geneName, geneID, geneStrand, colNames, Read_Isoform_compatibleVector,
                                                      self.qname_cbumi_dict,
                                                      self.metageneStructureInformationwNovel[meta_gene][index][1],
                                                      self.metageneStructureInformationwNovel[meta_gene][index][2],
                                                      sample_target)
                    else:
                        return_samples.append({'Read_Isoform_compatibleVector': Read_Isoform_compatibleVector, 'isoforms': colNames,
                                'exonInfo': self.metageneStructureInformationwNovel[meta_gene][index][1],
                                'isoformInfo': self.metageneStructureInformationwNovel[meta_gene][index][2]})
                if save==False:
                    return return_samples
    def map_reads_allgenes(self, cover_existing = True, total_jobs = 1, current_job_index = 0):
        if not os.path.exists(self.compatible_matrix_folder_path):
            os.makedirs(self.compatible_matrix_folder_path)
        if cover_existing:
            print('If there are existing compatible matrix files, SCOTCH will overwrite them')
            genes_existing = []
        else:
            print('If there are existing compatible matrix files, SCOTCH will not overwrite them')
            genes_existing = [g[:-4] for g in os.listdir(self.compatible_matrix_folder_path)]
            if os.path.isfile(os.path.join(self.compatible_matrix_folder_path, 'log.txt')):
                gene_df = pd.read_csv(os.path.join(self.compatible_matrix_folder_path, 'log.txt'), header=None)
                genes_existing = genes_existing + gene_df.iloc[:, 0].tolist()
        MetaGene_Gene_dict = {}
        for key, values in self.metageneStructureInformation.items():
            genes_ = []
            for value in values:
                gene = str(value[0]['geneName']) + '_' + str(value[0]['geneID'])
                if gene not in genes_existing:
                    genes_.append(gene)
            if len(genes_) > 0:
                MetaGene_Gene_dict[key] = genes_
        MetaGenes = list(MetaGene_Gene_dict.keys())
        print('total metagene number is: ' + str(len(MetaGenes)))
        if total_jobs > 1:
            step_size = math.ceil(len(MetaGenes) / total_jobs)
            s = int(list(range(0, len(MetaGenes), step_size))[current_job_index])
            e = int(s + step_size)
            MetaGenes_ = MetaGenes[s:e]
            print('processing: ' + str(len(MetaGenes_)) + ' metagenes')
        else:#total_jobs = 1
            MetaGenes_ = MetaGenes
        if self.parse:
            for meta_gene in MetaGenes_:
                self.map_reads_parse(meta_gene, save=True)
        else:
            for meta_gene in MetaGenes_:
                self.map_reads(meta_gene, save=True)









