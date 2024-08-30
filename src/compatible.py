from preprocessing import *
import pysam

##TODO: add other platforms
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
        self.bamInfo_csv_path = os.path.join(target, 'bam/bam.Info.csv')
        # bam information file
        self.qname_dict = load_pickle(self.bamInfo_pkl_path)
        self.qname_cbumi_dict = load_pickle(self.bamInfo2_pkl_path)
        self.metageneStructureInformation = load_pickle(self.annotation_path_meta_gene)
        self.metageneStructureInformationwNovel = self.metageneStructureInformation.copy()
        # parameters
        self.lowest_match = lowest_match
        self.platform = platform
        # some paths
        self.compatible_matrix_folder_path = os.path.join(target, "compatible_matrix")
        self.read_mapping_path = os.path.join(target, "auxillary")
        self.count_matrix_folder_path = os.path.join(target, "count_matrix")
    def read_bam(self, chrom = None):
        # bam_path is a folder
        if os.path.isfile(self.bam_path) == False:
            # find the bam file
            bamFile_name = [f for f in os.listdir(self.bam_path) if
                            f.endswith('.bam') and '.' + chrom + '.' in f]
            bamFile = os.path.join(self.bam_path, bamFile_name[0])
            bamFilePysam = pysam.Samfile(bamFile, "rb")
        else:
            bamFilePysam = pysam.Samfile(self.bam_path, "rb")
        return bamFilePysam
    def map_reads(self, meta_gene, save = True):
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
                                      Info_multigenes)
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
                out = process_read_metagene(read, start, end, self.qname_dict, Info_multigenes, self.lowest_match)
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
        for meta_gene in MetaGenes_:
            self.map_reads(meta_gene, save=True)









