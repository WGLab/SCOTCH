import pandas as pd
import csv
import os

def read_top(path, N):
    header = []
    lines = []
    with open(path, 'rt') as file:
        for i in range(N):
            line = next(file).strip()
            if not line.startswith('#'):
                line = line.split('\t')
                lines.append(line)
            else:
                header.append(line)
    return header, lines


def parse_info(x: str,sep=';',inner_sep = ' ')->dict:
    x=x.split(sep)
    x = [x_ for x_ in x if x_!='']
    out = []
    for i in x:
        pair = [p for p in i.split(inner_sep) if p!='']
        if len(pair)==2 and (pair[0] in ['gene_name','gene_id','transcript_id','gene_type']):
            pair = (pair[0],pair[1].replace('\"',''))
            out.append(pair)
    return dict(out)

def lines2df(lines: list, sep=';',inner_sep = ' ', parse = True):
    df = pd.DataFrame(lines)
    if parse:
        info = list(df.iloc[:, -1])
        info = [parse_info(x, sep, inner_sep) for x in info]
        df2 = pd.DataFrame(info)
        df1 = df.iloc[: , :-1]
        df = pd.concat([df1, df2], axis=1)
    return df

def read_gtf(path, N, condition=None):
    header = []
    lines = []
    with open(path, 'rt') as file:
        if N is not None:
            for i in range(N):
                line = next(file).strip()
                if not line.startswith('#'):
                    line = line.split('\t')
                    if condition is not None:
                        if line[2] in condition:
                            lines.append(line)
                    else:
                        lines.append(line)
                else:
                    header.append(line)
        else:
            for line in file:
                if not line.startswith('#'):
                    line = [x for x in line.split('\t')]
                    if condition is not None:
                        if line[2] in condition:
                            lines.append(line)
                    else:
                        lines.append(line)
        df = lines2df(lines, parse=True)
    return header, lines, df


def generate_reference_df(gtf_path):
    outdir = os.path.dirname(gtf_path)
    _, _, df = read_gtf(gtf_path, None, ['gene', 'exon'])
    df = df.dropna(axis=0, subset=['gene_name']).reset_index(drop=True)
    df.columns = ['CHR', 'SOURCE', 'TYPE', 'START', 'END', '.', 'STRAND', '.', 'gene_id','gene_type','gene_name',  'transcript_id']
    df_gene = df[df.iloc[:, 2] == "gene"].loc[:, ['CHR', 'START', 'END', 'gene_id','gene_name', 'gene_type','STRAND']].reset_index(
        drop=True)
    df_exon = df[df.iloc[:, 2] == "exon"].loc[:, ['CHR', 'START', 'END','gene_id', 'gene_name', 'transcript_id']].reset_index(
        drop=True)
    outfile_gene = os.path.join(outdir,'gene_gtf.bed')
    outfile_exon = os.path.join(outdir, 'exon_gtf.bed')
    #change gtf 1-based to 0-based
    df_gene.START = [str(int(s)-1) for s in df_gene.START.tolist()]
    df_exon.START = [str(int(s) - 1) for s in df_exon.START.tolist()]
    df_gene.to_csv(outfile_gene,sep = '\t', header = False, index = False, quoting=csv.QUOTE_NONE)
    df_exon.to_csv(outfile_exon, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
    #sort
    input_path = os.path.join(outdir,'exon_gtf.bed')
    output_path = os.path.join(outdir,'exon_gtf.sort.bed')
    sort_command = f'sort -k1,1 -k2,2n -k3,3n {input_path} > {output_path}'
    os.system(sort_command)
    #partition
    output_path2 = os.path.join(outdir, 'exon_gtf.part.bed')
    part_command = f'bedops -p {output_path} > {output_path2}'
    os.system(part_command)
    #intersect
    output_path3 = os.path.join(outdir, 'exon_gtf.part2.bed')
    intersect_command = f'bedtools intersect -a {output_path2} -b {output_path} -wa -wb > {output_path3}'
    os.system(intersect_command)
    #sort gene
    output_path4 = os.path.join(outdir, 'gene_gtf.sort.bed')
    sort_gene_command = f'sort -k1,1 -k2,2n -k3,3n {outfile_gene} > {output_path4}'
    os.system(sort_gene_command)
    # meta-gene
    outfile_gene_meta = os.path.join(outdir, 'meta_gene_gtf.bed')
    meta_gene_command = f'bedtools cluster -i {output_path4}  > {outfile_gene_meta}'
    os.system(meta_gene_command)
    #read file
    exons = pd.read_csv(output_path3, sep='\t', header=None, low_memory=False)
    exons.columns = ['CHR', 'START', 'END', 'CHR2', 'START2', 'END2', 'GENE_ID','GENE', 'TRANSCRIPT']
    exons = exons.loc[:, ['CHR', 'START', 'END', 'GENE_ID','GENE', 'TRANSCRIPT']]
    genes = pd.read_csv(outfile_gene_meta, sep='\t', header=None, low_memory=False)
    #remove
    remove_command = f'rm {output_path} {output_path2} {output_path3} {output_path4}'
    os.system(remove_command)
    return genes, exons
















