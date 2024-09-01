import argparse
import count_matrix as cm
import os
import annotation as annot
import compatible as cp

parser = argparse.ArgumentParser(description='Preprocess files')
#mandatory options
parser.add_argument('--task',type=str,help="choose task from annotation, compatible matrix, count matrix or all")
parser.add_argument('--platform',type=str,default='10x',help="platform: 10x, parse, or pacbio")
parser.add_argument('--target',type=str,nargs='+', help="path to target root folders for output files")#a list
parser.add_argument('--bam',type=str,nargs='+', help="one or multiple bam file paths or bam folder paths")#a list
#mandatory for parse
parser.add_argument('--build',type=str, default='hg38',help="genome build")
#task is annotation
parser.add_argument('--reference',type=str, default='../data/geneStructureInformation.pkl',help="Path to gene annotation file in either gtf format or the geneStructureInformation.pkl file. Leave it blank if using annotation-free mode. Default is human hg38 build")
parser.add_argument('--update_gtf', action='store_true', help='use bam file to update existing gtf annotations')
parser.add_argument('--update_gtf_off', action='store_false',dest='update_gtf',help='do NOT use bam file to update existing gtf annotations')
parser.add_argument('--coverage_threshold_gene',type=int, default= 20, help="coverage threshold to support gene discovery")
parser.add_argument('--coverage_threshold_exon',type=float, default=0.02, help="coverage threshold to support exon discovery, percentage to the maximum coverage")
parser.add_argument('--coverage_threshold_splicing',type=float, default=0.02, help="threshold to support splicing discovery, percentage to the maximum splicing junctions")
parser.add_argument('--z_score_threshold',type=int, default=10, help="threshold to support exon coverage sharp change discovery")
parser.add_argument('--min_gene_size',type=int, default=50, help="minimal length of novel discovered gene")

#task is matrix
parser.add_argument('--job_index',type=int, default=0, help="work array index")
parser.add_argument('--total_jobs',type=int, default=1, help="number of subwork")
parser.add_argument('--cover_existing',action='store_true')
parser.add_argument('--cover_existing_false', action='store_false',dest='cover_existing')

#task is count
parser.add_argument('--novel_read_n',type=int, default=0, help="filter out novel isoforms with supporting read number smaller than n")
#optional
parser.add_argument('--match',type=float,default=0.8, help="the lowest base percentage for matching an exon")
parser.add_argument('--workers',type=int,default=8, help="number of workers per work")


#bam = '/scr1/users/xu3/singlecell/project_singlecell/sample8_R9/bam/sample8_R9.filtered.bam'

def main():
    global args
    args = parser.parse_args()
    for t in args.target:
        if not os.path.exists(t):
            os.makedirs(t)
    def run_annotation():
        annotator = annot.Annotator(target=args.target, reference_gtf_path=args.reference,
                                    bam_path=args.bam, update_gtf=args.update_gtf,
                                    workers=args.workers, coverage_threshold_gene=args.coverage_threshold_gene,
                                    coverage_threshold_exon=args.coverage_threshold_exon,
                                    coverage_threshold_splicing=args.coverage_threshold_splicing,
                                    z_score_threshold=args.z_score_threshold,
                                    min_gene_size=args.min_gene_size, build=args.build, platform=args.platform)
        # generate gene annotation
        annotator.annotate_genes()
        # bam information
        annotator.annotation_bam()

    def run_compatible():
        for i in range(len(args.target)):
            print('processing the sample of: '+str(args.bam[i]))
            readmapper = cp.ReadMapper(target=args.target[i], bam_path = args.bam[i],
                                       lowest_match=args.match, platform = args.platform)
            if args.platform == 'parse':
                readmapper.merge_bam()
            print('generating compatible matrix')
            readmapper.map_reads_allgenes(cover_existing=args.cover_existing,
                                          total_jobs=args.total_jobs,current_job_index=args.job_index)
            print('saving annotations with identified novel isoforms')
            readmapper.save_annotation_w_novel_isoform(total_jobs=args.total_jobs,current_job_index=args.job_index)

    def run_count():
        for i in range(len(args.target)):
            print('generating count matrix at: ' + str(args.target[i]))
            countmatrix = cm.CountMatrix(target = args.target[i], novel_read_n = args.novel_read_n,
                           platform = args.platform, workers = args.workers)
            if args.platform=='parse':
                countmatrix.generate_multiple_samples()
                countmatrix.save_multiple_samples(csv=True, mtx=True)
            else:
                countmatrix.generate_single_sample()
                countmatrix.save_single_sample(csv=True, mtx=True)
    def run_summary():
        for i in range(len(args.target)):
            cp.summarise_annotation(args.target[i])

    if args.task=='annotation':
        run_annotation()
    if args.task=='compatible matrix':#task is to generate compatible matrix
        run_compatible()
    if args.task == 'count matrix': # task is to generate count matrix
        run_count()
    if args.task =='summary':
        run_summary()

    if args.task =='all':
        run_annotation()
        run_compatible()
        run_count()
        run_summary()




if __name__ == '__main__':
    main()