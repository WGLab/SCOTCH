import argparse
import count_matrix as cm
import os
import annotation as annot
import compatible as cp
import logging
import shutil
import visualization as vis
import subprocess

parser = argparse.ArgumentParser(description='SCOTCH preprocessing pipeline')
#mandatory options
parser.add_argument('--task',type=str,help="choose task from annotation, compatible matrix, count matrix, summary, or all; or visualization")#compatible matrix splicing
parser.add_argument('--generate_splicing', action='store_true', default= False, help="do not generate spliced/unspliced count matrix for velocity calculation")
parser.add_argument('--platform',type=str,default='10x-ont',help="platform: 10x-ont, parse-ont, or 10x-pacbio")
parser.add_argument('--target',type=str,nargs='+', help="path to target root folders for output files")#a list
parser.add_argument('--bam',type=str,nargs='+', help="one or multiple bam file paths or bam folder paths")#a list
parser.add_argument('--build',type=str,help="genome build, bam files for parse platform have contigs in the format of build_chr1")
#task is annotation
parser.add_argument('--reference',type=str, default='/scr1/users/xu3/singlecell/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf',help="Path to gene annotation file in gtf format. Set it as None if using annotation-free mode. To save time, users can use SCOTCH generated annotation file usint the argument of --reference_pkl. But the gtf file is still needed.")
parser.add_argument('--reference_pkl',type=str, default='../data/geneStructureInformation.pkl',help="Path to gene annotation file generated by SCOTCH named with geneStructureInformation.pkl file. Set this argument and --reference as None if using annotation-free mode. Set this argument to None and set --reference to use known annotation file, then SCOTCH will generate the annotation pickle file. Default is in human hg38 build")
parser.add_argument('--update_gtf', action='store_true', help='use bam file to update existing gtf annotations')
parser.add_argument('--update_gtf_off', action='store_false',dest='update_gtf',help='do NOT use bam file to update existing gtf annotations')
parser.add_argument('--coverage_threshold_gene',type=int, default= 20, help="coverage threshold to support gene discovery")
parser.add_argument('--coverage_threshold_exon',type=float, default=0.02, help="coverage threshold to support exon discovery, percentage to the maximum coverage")
parser.add_argument('--coverage_threshold_splicing',type=float, default=0.02, help="threshold to support splicing discovery, percentage to the maximum splicing junctions")
parser.add_argument('--z_score_threshold',type=int, default=10, help="threshold to support exon coverage sharp change discovery")
parser.add_argument('--min_gene_size',type=int, default=50, help="minimal length of novel discovered gene")
parser.add_argument('--barcode_cell',type=str, help="cell barcode tag in bam file")
parser.add_argument('--barcode_umi',type=str,  help="umi barcode tag in bam file")

#task is compatible matrix
parser.add_argument('--job_index',type=int, default=0, help="work array index")
parser.add_argument('--total_jobs',type=int, default=1, help="number of subwork")
parser.add_argument('--cover_existing',action='store_true')
parser.add_argument('--cover_existing_false', action='store_false',dest='cover_existing')
parser.add_argument('--small_exon_threshold',type=int,default=0, help="dynamic exon length threshold to ignore for includsion and exclusion")
parser.add_argument('--small_exon_threshold_high',type=int,default=80, help="the upper bound of dynamic exon length threshold to ignore for includsion and exclusion")
parser.add_argument('--truncation_match',type =float, default=0.4, help="higher than this threshold at the truncation end will be adjusted to 1")
parser.add_argument('--match_low',type=float,default=0.1, help="the base percentage to call a read-exon unmatched")
parser.add_argument('--match_high',type=float,default=0.6, help="the base percentage to call a read-exon matched")

#task is compatible matrix splicing
parser.add_argument('--unsplice_threshold',type=int, default=15, help="threshold for deciding unspliced reads")

#task is count
parser.add_argument('--novel_read_n',type=int, default=0, help="filter out novel isoforms with supporting read number smaller than n")
parser.add_argument('--group_novel', action='store_true', help="whether to further group novel isoforms generated in compatible matrix, default is true")
parser.add_argument('--group_novel_off', action='store_false', dest='group_novel')
parser.add_argument('--save_csv', action='store_true', help="whether to save count matrix output as csv format")
parser.add_argument('--save_csv_false', action='store_false', dest='save_csv')
parser.add_argument('--save_mtx', action='store_true', help="whether to save count matrix output as mtx format")
parser.add_argument('--save_mtx_false', action='store_false', dest='save_mtx')

#general
parser.add_argument('--workers',type=int,default=8, help="number of workers per work")
parser.add_argument('--single_cell',action='store_true',help="default setting for preprocessing single cell data")
parser.add_argument('--bulk', action='store_false',dest='single_cell',help="bulk data")


#task is visualization
parser.add_argument('--gene',type=str, help="gene name to visualize")
parser.add_argument('--target_vis',type=str, help="target path for visualization")
parser.add_argument('--sample_names',type=str,nargs='+', help="sample names for visualization")#a list
parser.add_argument('--novel_pct',type=float,default=0.1, help="only keep novel isoform annotation if expression in count matrix surpass the threshold")
parser.add_argument('--junction_num',type=int,default=10, help="only keep junctions over this number")
parser.add_argument('--width',type=int,default=12)
parser.add_argument('--height',type=float,default=1)
parser.add_argument("--annotation_scale", type=float, default=0.25, help="annotation plot scale")

def setup_logger(target, task_name):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    if logger.hasHandlers():
        logger.handlers.clear()
    log_file = os.path.join(target, f'{task_name}.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.ERROR)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger, log_file

def copy_log_to_targets(log_file, targets):
    for target in targets[1:]:
        shutil.copy(log_file, os.path.join(target, 'annotation.log'))


def main():
    global args
    args = parser.parse_args()
    if args.reference_pkl=='None':
        args.reference_pkl = None
    if args.platform=='10x-ont':
        if args.barcode_cell is None:
            args.barcode_cell = "CB"
        if args.barcode_umi is None:
            args.barcode_umi = "UB"
    elif args.platform == '10x-pacbio':
        if args.barcode_cell is None:
            args.barcode_cell = "CB"
        if args.barcode_umi is None:
            args.barcode_umi = "XM"
    for t in args.target:
        if not os.path.exists(t):
            os.makedirs(t)
    def run_annotation():
        logger, log_file = setup_logger(args.target[0], 'annotation')
        logger.info('Start annotation step for all targets.')
        logger.info(f'Target directories: {args.target}')
        logger.info(f'Reference annotation file: {args.reference}')
        logger.info(f'Reference annotation pickle file: {args.reference_pkl}')
        logger.info(f'Update GTF option: {args.update_gtf}')
        logger.info(f'BAM files: {args.bam}')
        logger.info(f'Platform: {args.platform}')
        logger.info(f'Coverage threshold (gene): {args.coverage_threshold_gene}')
        logger.info(f'Coverage threshold (exon): {args.coverage_threshold_exon}')
        logger.info(f'Coverage threshold (splicing): {args.coverage_threshold_splicing}')
        logger.info(f'Z-score threshold: {args.z_score_threshold}')
        logger.info(f'Minimum gene size: {args.min_gene_size}')
        annotator = annot.Annotator(target=args.target,
                                    reference_gtf_path=args.reference, reference_pkl_path = args.reference_pkl,
                                    bam_path=args.bam, update_gtf=args.update_gtf,
                                    workers=args.workers, coverage_threshold_gene=args.coverage_threshold_gene,
                                    coverage_threshold_exon=args.coverage_threshold_exon,
                                    coverage_threshold_splicing=args.coverage_threshold_splicing,
                                    z_score_threshold=args.z_score_threshold,
                                    min_gene_size=args.min_gene_size, build=args.build, platform=args.platform,
                                    logger=logger)
        # generate gene annotation
        logger.info('Start generating gene annotation.')
        annotator.annotate_genes()
        # bam information
        logger.info('Start processing bam file information.')
        annotator.annotation_bam(args.barcode_cell, args.barcode_umi)
        copy_log_to_targets(log_file, args.target)
    def run_compatible():
        logger, log_file = setup_logger(args.target[0], 'compatible')
        logger.info(f'Start generating compatible matrix step for all targets. Job: {args.job_index}')
        logger.info(f'Target directories: {args.target}. Job: {args.job_index}')
        logger.info(f'BAM files: {args.bam}. Job: {args.job_index}')
        logger.info(f'Lowest match: {args.match_low}. Job: {args.job_index}')
        logger.info(f'Highest match: {args.match_high}. Job: {args.job_index}')
        logger.info(f'Small exon threshold (low): {args.small_exon_threshold}. Job: {args.job_index}')
        logger.info(f'Small exon threshold (high): min(average exon length, {args.small_exon_threshold_high}). Job: {args.job_index}')
        logger.info(f'Truncation match: {args.truncation_match} or 100bps. Job: {args.job_index}')
        logger.info(f'Platform: {args.platform}. Job: {args.job_index}')
        logger.info(f'Reference GTF Path: {args.reference}. Job: {args.job_index}')
        logger.info(f'Update GTF option: {args.update_gtf}. Job: {args.job_index}')
        readmapper = cp.ReadMapper(target=args.target, bam_path = args.bam,
                                   lowest_match=args.match_low, lowest_match1=args.match_high,
                                    small_exon_threshold = args.small_exon_threshold,
                                   small_exon_threshold1=args.small_exon_threshold_high,
                                   truncation_match = args.truncation_match,
                                   platform = args.platform, reference_gtf_path=args.reference,
                                   logger = logger)
        readmapper.map_reads_allgenes(cover_existing=True,total_jobs=args.total_jobs,current_job_index=args.job_index)
        logger.info(f'saving annotations with identified novel isoforms  Job: {args.job_index}')
        readmapper.save_annotation_w_novel_isoform(total_jobs=args.total_jobs,current_job_index=args.job_index)
        logger.info(f'Completed generating compatible matrix for all targets.  Job: {args.job_index}')
        copy_log_to_targets(log_file, args.target)

    def run_compatible_splicing():
        logger, log_file = setup_logger(args.target[0], 'compatible matrix splicing')
        logger.info(f'Start generating spliced/unspliced compatible matrix for all targets. Job: {args.job_index}')
        logger.info(f'BAM files: {args.bam}. Job: {args.job_index}')
        logger.info(f'Unspliced threshold: {args.unsplice_threshold}. Job: {args.job_index}')
        cr = cp.ClassifyReadsSplice(scotch_target = args.target,
                                    bam_path = args.bam,
                                    unsplice_threshold = args.unsplice_threshold,
                                    n_jobs = args.total_jobs, job_index = args.job_index,
                                    logger=logger)
        cr.split_compatible()
        logger.info(f'Completed generating spliced/unspliced compatible matrix for all targets.  Job: {args.job_index}')

    def run_count():
        # target has to be len 1 for parse platform
        logger, log_file = setup_logger(args.target[0], 'count')
        logger.info('Start generating count matrix for all targets.')
        logger.info(f'Spliced/unspliced count matrix generation is set as {args.generate_splicing}')
        logger.info(f'Target directories: {args.target}')
        logger.info(f'Novel read threshold: {args.novel_read_n}')
        logger.info(f'Platform: {args.platform}')
        logger.info(f'Group novel isoforms: {args.group_novel}')
        logger.info(f'Workers: {args.workers}')
        logger.info(f'saving count matrix csv: {args.save_csv}')
        logger.info(f'saving count matrix mtx: {args.save_mtx}')
        countmatrix = cm.CountMatrix(target = args.target, novel_read_n = args.novel_read_n,
                                        platform = args.platform, workers = args.workers,
                                     group_novel = args.group_novel, logger = logger,
                                     csv = args.save_csv, mtx = args.save_mtx)
        if args.platform=='parse':
            assert len(args.target) == 1, "Error: The length of target must be 1 when platform is 'parse'."
        countmatrix.generate_multiple_samples(generate_splicing=args.generate_splicing)
        logger.info('Saving count matrix')
        countmatrix.save_multiple_samples(generate_splicing=args.generate_splicing)
        countmatrix.filter_gtf()
        logger.info('Completed generating count matrix for all targets.')
        copy_log_to_targets(log_file, args.target)

    def run_summary():
        logger, log_file = setup_logger(args.target[0], 'summary')
        logger.info('Start summarizing read isoform mapping and annotations for all targets.')
        logger.info(f'Target directories: {args.target}')
        for i in range(len(args.target)):
            logger.info(f'Start summarizing annotation for target: {args.target[i]}')
            cp.summarise_annotation(args.target[i])
            logger.info(f'Start summarizing read mapping information for target: {args.target[i]}')
            cp.summarise_auxillary(args.target[i])
        logger.info('Completed summarizing annotations and auxiliary information for all targets.')
        copy_log_to_targets(log_file, args.target)

    if args.task=='annotation':
        run_annotation()
    if args.task=='compatible matrix':#task is to generate compatible matrix
        run_compatible()
    if args.task == 'summary':
        run_summary()
    if args.task=='compatible matrix splicing':
        run_compatible_splicing()
    if args.task == 'count matrix': # task is to generate count matrix
        run_count()

    if args.task =='all':
        run_annotation()
        run_compatible()
        run_summary()
        run_count()

    if args.task =='visualization': #currently only support 10X file structure
        vis.visualization(args.gene, args.bam, args.target, args.novel_pct, args.junction_num, args.target_vis,
                          args.sample_names, args.width, args.height, args.annotation_scale)

if __name__ == '__main__':
    main()