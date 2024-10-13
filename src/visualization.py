from preprocessing import *
import os
import shutil
import pandas as pd
import pysam
import re
import subprocess
from scipy.io import mmread

#NFKBIA--535612368/535612383

def sub_gtf(gene_name, gtf_path, out_path = None):
    gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None)
    gtf.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    filtered_gtf = gtf[gtf['attribute'].str.contains(f'gene_name "{gene_name}"')]
    transcript_ids = filtered_gtf['attribute'].apply(
        lambda x: re.search(r'transcript_id "([^"]+)"', x).group(1) if re.search(r'transcript_id "([^"]+)"',
                                                                                 x) else None
    ).dropna().unique()
    gene_chr = filtered_gtf['seqname'].iloc[0]
    gene_start = filtered_gtf['start'].min()
    gene_end = filtered_gtf['end'].max()
    gene_strand = filtered_gtf['strand'].iloc[0]
    if out_path is not None:
        filtered_gtf.to_csv(out_path, sep='\t', index=False, header=False, quoting=3)
    return filtered_gtf, transcript_ids, gene_chr, gene_start, gene_end, gene_strand


def merge_gtf_files(gtf_paths: list, out_path=None):
    merged_gtf = pd.concat([pd.read_csv(gtf_path, sep='\t', comment='#', header=None) for gtf_path in gtf_paths], ignore_index=True)
    merged_gtf = merged_gtf.drop_duplicates()
    if out_path is not None:
        merged_gtf.to_csv(out_path, sep='\t', index=False, header=False,quoting=3)
    return merged_gtf


def generate_subbam_subgtf_single_sample(gene, bamFile, target, novel_pct=0.1):
    out_path = os.path.join(target, 'visualization', gene)
    gtf_path = os.path.join(target, 'reference/SCOTCH_updated_annotation_filtered.gtf')
    if not os.path.exists(out_path):
        os.makedirs(out_path, exist_ok=True)
    filtered_gtf, transcript_ids, gene_chr, gene_start, gene_end, gene_strand = sub_gtf(gene, gtf_path, out_path=None)
    #-------get major novel isoform names and known isoform names
    count_matrix_path = os.path.join(target, 'count_matrix')
    count_matrix_path_mtx = [os.path.join(count_matrix_path, i) for i in os.listdir(count_matrix_path) if
                             'adata_transcript' in i and '.mtx' in i][0]
    count_matrix_path_pkl = [os.path.join(count_matrix_path, i) for i in os.listdir(count_matrix_path) if
                             'adata_transcript' in i and '.pickle' in i][0]
    mtx = mmread(count_matrix_path_mtx)
    pkl = load_pickle(count_matrix_path_pkl)
    dense_mtx = mtx.todense()
    gene_cols_index = [i for i, col in enumerate(pkl['var']) if gene in col]
    sub_mtx = dense_mtx[:, gene_cols_index]
    col_pct = np.sum(sub_mtx, axis=0) / np.sum(sub_mtx)
    over_1_percent_cols = np.where(col_pct > novel_pct)[1].tolist()
    known_isoform_names = [tr for tr in transcript_ids if 'ENST' in tr]
    novel_isoform_names = [pkl['var'][gene_cols_index[ind]] for ind in over_1_percent_cols if 'novel' in pkl['var'][gene_cols_index[ind]]]
    novel_isoform_names = [isoform.split('_',1)[1] for isoform in novel_isoform_names]
    selected_isoform = known_isoform_names + novel_isoform_names
    transcript_ids = filtered_gtf['attribute'].apply(
        lambda x: re.search(r'transcript_id "([^"]+)"', x).group(1) if re.search(r'transcript_id "([^"]+)"',
                                                                                 x) else None)
    filtered_gtf = pd.concat([filtered_gtf.iloc[[0]], filtered_gtf[transcript_ids.isin(selected_isoform)]], ignore_index=True)
    #locate reads belonging to known and novel
    compatible_matrix_path = os.path.join(target,'compatible_matrix')
    compatible_matrix_path = os.path.join(compatible_matrix_path,[cp for cp in os.listdir(compatible_matrix_path) if cp.startswith(str(gene)+'_')][0])
    df=pd.read_csv(compatible_matrix_path)
    new_columns = df.columns.tolist()
    new_columns[0] = 'CBUMI'
    df.columns = new_columns
    mask_enst = df[known_isoform_names].gt(0).any(axis=1)
    cbumi_enst = df.loc[mask_enst, 'CBUMI'].tolist()
    ##TODO: add parse
    geneStrand = filtered_gtf.iloc[[0]]['strand'].tolist()[0]
    df_novel = df.filter(like='novelIsoform')
    child_id_list = [int(col.split('_')[1]) for col in df_novel.columns.tolist()]
    novel_isoform_id_list = [int(novel.split('_')[1]) for novel in novel_isoform_names]
    novel_cbumi_dict = {}
    for ii in range(len(novel_isoform_id_list)):
        col_index_list = []
        for i, child_id in enumerate(child_id_list):
            if find_parent_isoform([novel_isoform_id_list[ii]], child_id, geneStrand) is not None:
                col_index_list.append(i)
        if len(col_index_list)>0:
            mask_novel = df_novel.iloc[:,col_index_list].gt(0).any(axis=1)
            cbumi_novel = df.loc[mask_novel, 'CBUMI'].tolist()
            novel_cbumi_dict[novel_isoform_names[ii]] = cbumi_novel
    # write bam file
    if os.path.isfile(bamFile)==False: #bamFile is a folder
        bamFile_name = [f for f in os.listdir(bamFile) if f.endswith('.bam') and '.'+gene_chr+'.' in f]
        bamFile = os.path.join(bamFile,bamFile_name[0])
    bamFilePysam = pysam.Samfile(bamFile, "rb")
    reads = bamFilePysam.fetch(gene_chr, gene_start, gene_end)
    #-------------------write bam files------------------#
    bam_existing_path=os.path.join(out_path, gene+'_known_isoform.bam')
    bam_existing_sorted_path = os.path.join(out_path,gene+'_known_isoform_sorted.bam')
    bam_existing = pysam.AlignmentFile(bam_existing_path, "wb", template=bamFilePysam)
    bam_novel_files = {}
    bam_novel_paths = {}
    bam_novel_sorted_paths = {}
    for novel_isoform in novel_cbumi_dict.keys():
        bam_novel_path = os.path.join(out_path, f'{gene}_{novel_isoform}_isoform.bam')
        bam_novel_sorted_path = os.path.join(out_path, f'{gene}_{novel_isoform}_isoform_sorted.bam')
        bam_novel_files[novel_isoform] = pysam.AlignmentFile(bam_novel_path, "wb", template=bamFilePysam)
        # Store paths for sorting and indexing later
        bam_novel_paths[novel_isoform] = bam_novel_path
        bam_novel_sorted_paths[novel_isoform] = bam_novel_sorted_path
    for read in reads:
        CBUMI = read.get_tag('CB')+'_'+read.get_tag('UB')
        if CBUMI in cbumi_enst:
            bam_existing.write(read)
        else:
            for novel_isoform, cbumi_list in novel_cbumi_dict.items():
                if CBUMI in cbumi_list:
                    bam_novel_files[novel_isoform].write(read)
    bam_existing.close()
    for bam_file in bam_novel_files.values():
        bam_file.close()
    subprocess.run(["samtools", "sort", "-o", bam_existing_sorted_path, bam_existing_path])
    subprocess.run(["samtools", "index", bam_existing_sorted_path])
    for novel_isoform, bam_novel_path in bam_novel_paths.items():
        bam_novel_sorted_path = bam_novel_sorted_paths[novel_isoform]
        subprocess.run(["samtools", "sort", "-o", bam_novel_sorted_path, bam_novel_path])
        subprocess.run(["samtools", "index", bam_novel_sorted_path])
    filtered_gtf.to_csv(os.path.join(out_path, gene+'_SCOTCH_filtered.gtf'), sep='\t', index=False, header=False,quoting=3)
    return gene_chr, gene_start, gene_end, gene_strand

def generate_subbam_subgtf_multiple_samples(gene:str, bamFiles:list, targets:list, novel_pct=0.1,
                                            final_target = None, sample_names = None):
    gtf_paths, target_gene_folders = [], []
    gene_chr, gene_start, gene_end, gene_strand = None, None, None, None
    for bamFile, target in list(zip(bamFiles, targets)):
        gene_chr, gene_start, gene_end, gene_strand = generate_subbam_subgtf_single_sample(gene, bamFile, target, novel_pct)
        gtf_path = os.path.join(target, 'visualization', gene, gene+'_SCOTCH_filtered.gtf')
        gtf_paths.append(gtf_path)
        target_gene_folder = os.path.join(target, 'visualization', gene)
        target_gene_folders.append(target_gene_folder)
    if final_target is not None:
        final_target = os.path.join(final_target, 'visualization', gene)
        if sample_names is not None:
            if not os.path.exists(final_target):
                os.makedirs(final_target, exist_ok=True)
            for folder, sample_name in zip(target_gene_folders, sample_names):
                dest_folder = os.path.join(final_target, sample_name)
                if os.path.exists(dest_folder):
                    shutil.rmtree(dest_folder)  # Remove existing folder if any
                shutil.copytree(folder, dest_folder)
            merged_gtf_path = os.path.join(final_target, gene + '_merged.gtf')
            _ = merge_gtf_files(gtf_paths, out_path=merged_gtf_path)
    return gene_chr, gene_start, gene_end, gene_strand

#path, bam, label, color
def run_trackplot(target_vis, gene_name, gene_chr, gene_start, gene_end,gene_strand,
                  dpi=300, width=12, height=1, junction_num=100, intron_scale=1, annotation_scale=0.25):
    #target_vis, files are at target_vis/visualization/gene/sample1/
    gtf_file = os.path.join(target_vis, 'visualization', gene_name, gene_name+'_merged.gtf')
    folder_path = os.path.join(target_vis, 'visualization', gene_name)
    bamfile_tsv_path = os.path.join(folder_path, 'bamfile.tsv')
    output_pdf = os.path.join(folder_path, gene_name+'_vis.pdf')
    field1 =[]
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".bam") and 'sorted' in file:
                field1.append(os.path.join(root, file))
    field2 = ['bam'] * len(field1)
    field3_0 = [os.path.basename(os.path.dirname(path)) for path in field1]
    field3_1 = [
        (match.group(0) if (match := re.search(r'novelIsoform_\d+', os.path.basename(path))) else "known")
        for path in field1
    ]
    field3 = [f'{a}_{b}' for a, b in zip(field3_0, field3_1)]
    field4 = ['#2166ac' if ii=='known' else '#b2182b' for ii in field3_1 ]
    bamfile_tsv = pd.DataFrame({'path':field1, 'type':field2, 'label':field3, 'color':field4})
    bamfile_tsv = bamfile_tsv.sort_values(by='label').reset_index(drop=True)
    bamfile_tsv.to_csv(bamfile_tsv_path, sep='\t', index=False, header=False)
    genomic_region = f"{gene_chr}:{gene_start}-{gene_end}:{gene_strand}"
    command = [
        "trackplot",
        "-e", genomic_region,
        "-r", gtf_file,
        "--density", bamfile_tsv_path,
        "-o", output_pdf,  # Output PDF file path
        "--dpi", str(dpi),  # DPI resolution
        "--width", str(width),  # Plot width
        "--height", str(height),  # Plot height
        "--show-junction-num",  # Show junction numbers
        "-t", str(junction_num),  # threshold
        "--intron-scale", str(intron_scale),  # Intron scaling factor
        "--same-y",
        "--distance-ratio", '0.1',
        "--annotation-scale", str(annotation_scale)
    ]
    try:
        subprocess.run(command, check=True)
        print("Trackplot generated successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while generating the plot: {e}")

def visualization(gene, bamFiles, targets, novel_pct, junction_num, final_target, sample_names, width, height, annotation_scale):
    output = generate_subbam_subgtf_multiple_samples(gene=gene, bamFiles=bamFiles, targets=targets,
                                                         novel_pct=novel_pct, final_target=final_target,
                                                         sample_names=sample_names)
    gene_chr, gene_start, gene_end, gene_strand = output
    run_trackplot(target_vis=final_target, gene_name=gene, gene_chr=gene_chr, gene_start=gene_start, gene_end=gene_end,
                  gene_strand = gene_strand,dpi=300, width=width, height=height, junction_num=junction_num, intron_scale=1,
                  annotation_scale = annotation_scale)

#gene='CD74'
#bamFile='/mnt/isilon/wang_lab/xinya/projects/single_cell_pipeline/CAG_SingleCell/sample7-R10-allpass-v4/wf-single-cell-v1-output-sample7R10-allpass-ed1/reseq/bams'
#target='/scr1/users/xu3/singlecell/project_singlecell/sample7_8_ont/sample7'





#trackplot \
#-e chr19:14031761-14053218 \
#-r genes_IL27RA.gtf \
#--density bamfile.tsv \
#-o IL27RA_sashimi.pdf \
#--dpi 300 \
#--width 12 \
#--height 1 --show-junction-num -t 10 \
#--intron-scale 1


#output = generate_subbam_subgtf_multiple_samples(gene='CD74',
#                                                 bamFiles=['/mnt/isilon/wang_lab/xinya/projects/single_cell_pipeline/CAG_SingleCell/sample7-R10-allpass-v4/wf-single-cell-v1-output-sample7R10-allpass-ed1/reseq/bams',
#                                                           '/mnt/isilon/wang_lab/xinya/projects/single_cell_pipeline/CAG_SingleCell/sample8-R10-allpass-v4/wf-single-cell-v1-output-sample8R10-allpass-ed1/reseq/bams'],
#                                                 targets=[
#                                                     '/scr1/users/xu3/singlecell/project_singlecell/sample7_8_ont/sample7',
#                                                     '/scr1/users/xu3/singlecell/project_singlecell/sample7_8_ont/sample8'],
#                                                 novel_pct=0.1,
#                                                 final_target='/scr1/users/xu3/singlecell/project_singlecell/sample7_8_ont/visualization',
#                                                 sample_names=['sample7', 'sample8'])
#gene_chr, gene_start, gene_end = output