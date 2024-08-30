import argparse
import preprocessing as pp
import count_matrix as cm
import os
import annotation as annot
import compatible as cp
#from scipy.sparse import csr_matrix, save_npz
from scipy.io import mmwrite
import pickle

parser = argparse.ArgumentParser(description='Preprocess files')
#mandaroty options
parser.add_argument('--target',type=str,help="path to target root folder for output files")
parser.add_argument('--task',type=str,help="choose task between annotation, matrix(read * isoform), and count (count matrix)")
parser.add_argument('--bam',type=str,help="Path to bam file or bam folder")

#task is annotation
parser.add_argument('--ref',type=str, help="Path to gene annotation file in gtf format, output pickel file, leave it blank if using annotation-free mode")
parser.add_argument('--update_gtf', action='store_true', help='use bam file to update existing gtf annotations')
parser.add_argument('--update_gtf_off', action='store_false',dest='update_gtf',help='do NOT use bam file to update existing gtf annotations')
parser.add_argument('--coverage_threshold_gene',type=int, default= 5, help="coverage threshold to support gene discovery")
parser.add_argument('--coverage_threshold_exon',type=int, default=0.01, help="coverage threshold to support exon discovery, percentage to the maximum coverage")
parser.add_argument('--coverage_threshold_splicing',type=int, default=0.01, help="threshold to support splicing discovery, percentage to the maximum splicing junctions")
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
    if not os.path.exists(args.target):
        os.makedirs(args.target)
    if args.task=='annotation':
        annotator = annot.Annotator(target=args.target, reference_gtf_path = args.ref,
                                    bam_path = args.bam, update_gtf = args.update_gtf,
                                    workers = args.workers,coverage_threshold_gene = args.coverage_threshold_gene,
                                    coverage_threshold_exon = args.coverage_threshold_exon, coverage_threshold_splicing = args.coverage_threshold_splicing,
                                    min_gene_size = args.min_gene_size)
        #generate gene annotation
        annotator.annotate_genes()
        #bam information
        annotator.annotation_bam()
    elif args.task=='matrix':#task is to generate compatible matrix
        readmapper = cp.ReadMapper(target=args.target, bam_path = args.bam,
                                   lowest_match=args.match, platform = '10x')
        print('generating compatible matrix')
        readmapper.map_reads_allgenes(cover_existing=args.cover_existing, total_jobs=args.total_jobs,
                                      current_job_index=args.job_index)
    else: # task is to generate count matrix
        adata_gene_unfiltered, adata_transcript_unfiltered, adata_gene_filtered, adata_transcript_filtered = cm.generate_count_matrix(path=args.target, novel_read_n=args.novel_read_n, num_cores=args.workers)
        print('count matrix generated')
        output_gene_unfiltered = os.path.join(args.target,'count_matrix/adata_gene_unfiltered'+str(args.novel_read_n)+'.csv')
        output_transcript_unfiltered = os.path.join(args.target,  'count_matrix/adata_transcript_unfiltered'+str(args.novel_read_n)+'.csv')
        output_gene_filtered = os.path.join(args.target,'count_matrix/adata_gene_filtered' + str(args.novel_read_n) + '.csv')
        output_transcript_filtered = os.path.join(args.target, 'count_matrix/adata_transcript_filtered' + str(args.novel_read_n) + '.csv')
        # save gene
        print('saving count matrix on gene level ')
        adata_gene_unfiltered_df = adata_gene_unfiltered.to_df()
        adata_gene_filtered_df = adata_gene_filtered.to_df()
        adata_gene_unfiltered_df.to_csv(output_gene_unfiltered)
        adata_gene_filtered_df.to_csv(output_gene_filtered)
        gene_meta_unfiltered = {'obs':adata_gene_unfiltered.obs.index.tolist(), "var":adata_gene_unfiltered.var.index.tolist()}
        gene_meta_filtered = {'obs': adata_gene_filtered.obs.index.tolist(),"var": adata_gene_filtered.var.index.tolist()}
        with open(os.path.join(args.target,'count_matrix/adata_gene_unfiltered'+str(args.novel_read_n)+'.pickle'), 'wb') as f:
            pickle.dump(gene_meta_unfiltered, f)
        mmwrite(os.path.join(args.target,'count_matrix/adata_gene_unfiltered'+str(args.novel_read_n)+'.mtx'),adata_gene_unfiltered.X)
        with open(os.path.join(args.target,'count_matrix/adata_gene_filtered'+str(args.novel_read_n)+'.pickle'), 'wb') as f:
            pickle.dump(gene_meta_filtered, f)
        mmwrite(os.path.join(args.target, 'count_matrix/adata_gene_filtered' + str(args.novel_read_n) + '.mtx'),
                adata_gene_filtered.X)
        # save transcript
        print('saving count matrix on transcript level ')
        adata_transcript_unfiltered_df = adata_transcript_unfiltered.to_df()
        adata_transcript_filtered_df = adata_transcript_filtered.to_df()
        adata_transcript_unfiltered_df.to_csv(output_transcript_unfiltered)
        adata_transcript_filtered_df.to_csv(output_transcript_filtered)
        transcript_meta_unfiltered = {'obs': adata_transcript_unfiltered.obs.index.tolist(), "var": adata_transcript_unfiltered.var.index.tolist()}
        mmwrite(os.path.join(args.target, 'count_matrix/adata_transcript_unfiltered' + str(args.novel_read_n) + '.mtx'), adata_transcript_unfiltered.X)
        transcript_meta_filtered = {'obs': adata_transcript_filtered.obs.index.tolist(),"var": adata_transcript_filtered.var.index.tolist()}
        mmwrite(os.path.join(args.target, 'count_matrix/adata_transcript_filtered' + str(args.novel_read_n) + '.mtx'),
                adata_transcript_filtered.X)
        with open(os.path.join(args.target, 'count_matrix/adata_transcript_unfiltered' + str(args.novel_read_n) + '.pickle'),
                  'wb') as f:pickle.dump(transcript_meta_unfiltered, f)
        with open(os.path.join(args.target, 'count_matrix/adata_transcript_filtered' + str(args.novel_read_n) + '.pickle'),
                  'wb') as f:pickle.dump(transcript_meta_filtered, f)
    return


if __name__ == '__main__':
    main()