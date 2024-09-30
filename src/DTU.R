library(SCOTCH)
library(argparse)

parser <- ArgumentParser(description = "SCOTCH DTU gene analysis")

# Add arguments
parser$add_argument("--gene_count1", required=TRUE, help="Path to the first gene count matrix file")
parser$add_argument("--gene_count2", required=TRUE, help="Path to the second gene count matrix file")

parser$add_argument("--transcript_count1", required=TRUE, help="Path to the first transcript count matrix file")
parser$add_argument("--transcript_count2", required=TRUE, help="Path to the second transcript count matrix file")

parser$add_argument("--nworkers", type="integer", required=FALSE, default=1, help="Number of workers")
parser$add_argument("--out", required=TRUE, help="Folder to save output")
parser$add_argument("--file_name", required=TRUE, default = 'SCOTCH_DTU.csv', help="File name for the output file")

# Parse command-line arguments
args <- parser$parse_args()



#----read gene-level count matrix-----#
print('reading in gene level count matrix:')
print(paste0('File 1: ',args.gene_count1))
print(paste0('File 2: ',args.gene_count2))
gene_mat1=as.matrix(read.csv(args.gene_count1,row.names = 'X'))
gene_mat2=as.matrix(read.csv(args.gene_count2,row.names = 'X'))


#----read transcript-level count matrix-----#
transcript_mat1=as.matrix(read.csv(args.transcript_count1,row.names = 'X'))
gene_transcript_df1 = data.frame(genes=str_remove(colnames(transcript_mat1),"[_-](ENST|novel|uncategorized).+"),
                                    transcripts=colnames(transcript_mat1))

transcript_mat2=as.matrix(read.csv(args.transcript_count2,row.names = 'X'))
gene_transcript_df2 = data.frame(genes=str_remove(colnames(transcript_mat2),"[_-](ENST|novel|uncategorized).+"),
                                 transcripts=colnames(transcript_mat2))


#----gene-level analysis-----#
print('Performing differential gene expression analysis......')
df_gene = scotch_gene(gene_mat1, gene_mat2, epsilon=0.01,ncores=10)%>%
  filter(pct1>=0.01|pct2>=0.01)

#----transcript-level analysis-----#
print('Performing differential transcript usage analysis......')
df_transcript = scotch_transcript(gene_transcript_df1,gene_transcript_df2, 
                                  transcript_mat1, transcript_mat2, ncores=args.nworkers)
df_scotch = df_gene%>%left_join(df_transcript,by=c('genes'='gene'))

print('Saving DTU results')
write.csv(df_scotch, file.path(args.out, args.file_name))


