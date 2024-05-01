import argparse
import preprocessing as pp
import count_matrix as cm
import os
import pandas as pd
from scipy.sparse import csr_matrix, save_npz
from scipy.io import mmwrite
import pickle

parser = argparse.ArgumentParser(description='Preprocess files')
#mandaroty options
parser.add_argument('--target',type=str,help="path to target root folder for output files")
parser.add_argument('--task',type=str,help="choose task between annotation, matrix(read * isoform), and count (count matrix)")

#give one of the arguments below
parser.add_argument('--ref',type=str,default='/scr1/users/xu3/singlecell/ref/hg38_ensemble/genes',
                    help="Path to gene annotation file in gtf format, output pickel file")
parser.add_argument('--geneinfo',type=str, help="path to pickle file of gene annotations")

#task is matrix
parser.add_argument('--bam',type=str,help="Path to bam file")
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
        #gene annotation information
        print('extracting annotation information')
        if not os.path.exists(os.path.join(args.target,'reference')):
            os.makedirs(os.path.join(args.target,'reference'))
        output = os.path.join(args.target,"reference/geneStructureInformation.pkl")
        output2 = os.path.join(args.target, "reference/metageneStructureInformation.pkl")
        if os.path.isfile(output) and os.path.isfile(output2):
            print('gene annotation information exist')
        else:
            _ = pp.extract_annotation_info(args.ref, num_cores=args.workers, output=output)
        #bam information
        if not os.path.exists(os.path.join(args.target,'bam')):
            os.makedirs(os.path.join(args.target,'bam'))
        bamInfo_pkl_file = os.path.join(args.target, 'bam/bam.Info.pkl')
        bamInfo2_pkl_file = os.path.join(args.target, 'bam/bam.Info2.pkl')
        bamInfo_file = os.path.join(args.target, 'bam/bam.Info.csv')
        if os.path.isfile(bamInfo_pkl_file) == True and os.path.isfile(bamInfo_file) == True:
            print('bam file information exist')
        if os.path.isfile(bamInfo_pkl_file) == False and os.path.isfile(bamInfo_file) == True:
            print('extracting bam file pickle information')
            bam_info = pd.read_csv(bamInfo_file)
            qname_dict, qname_cbumi_dict = pp.bam_info_to_dict(bam_info)
            with open(bamInfo_pkl_file, 'wb') as file:
                pickle.dump(qname_dict, file)
            with open(bamInfo2_pkl_file, 'wb') as file:
                pickle.dump(qname_cbumi_dict, file)
        if os.path.isfile(bamInfo_pkl_file) == False and os.path.isfile(bamInfo_file) == False:
            print('extracting bam file information')
            if os.path.isfile(args.bam)==False:
                bam_info = pp.extract_bam_info_folder(args.bam, args.workers)
            else:
                bam_info = pp.extract_bam_info(args.bam)
            bam_info.to_csv(bamInfo_file)
            print('generating bam file pickle information')
            qname_dict, qname_cbumi_dict = pp.bam_info_to_dict(bam_info)
            with open(bamInfo_pkl_file, 'wb') as file:
                pickle.dump(qname_dict, file)
            with open(bamInfo2_pkl_file, 'wb') as file:
                pickle.dump(qname_cbumi_dict, file)
    elif args.task=='matrix':#task is to generate compatible matrix
        bamInfo_pkl_file = os.path.join(args.target, 'bam/bam.Info.pkl')
        bamInfo2_pkl_file = os.path.join(args.target, 'bam/bam.Info2.pkl')
        qname_dict = pp.load_pickle(bamInfo_pkl_file)
        qname_cbumi_dict = pp.load_pickle(bamInfo2_pkl_file)
        metagene_pkl = pp.load_pickle(os.path.join(args.target, "reference/metageneStructureInformation.pkl"))
        print('generating compatible matrix')
        pp.generate_compatibleVector(args.bam, qname_dict, qname_cbumi_dict, metagene_pkl, args.match, args.target,
                                     args.job_index, args.total_jobs, args.cover_existing)
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