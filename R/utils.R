

QC = function(seurat_object,organism='human', nFeature_RNA_threshold=NA,mt_threshold=NA){
  require(Seurat)
  if (organism=='mouse'){
    pattern = "^mt-"
  }else{pattern = "^MT-"}
  seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = pattern)
  if (is.na(nFeature_RNA_threshold)==T){nFeature_RNA_threshold = quantile(seurat_object@meta.data$nFeature_RNA,0.99)}
  if (is.na(mt_threshold)==T){mt_threshold = quantile(seurat_object@meta.data$percent.mt,0.99)}
  seurat_object <- subset(seurat_object, 
                          subset = nFeature_RNA > 200 & 
                            nFeature_RNA < nFeature_RNA_threshold & 
                            percent.mt < mt_threshold)
  p=VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0.05)
  return(list(object = seurat_object,plot=p))
}

preprocess = function(seurat_object,dims,resolution,seed=42){
  require(Seurat)
  #normalization
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  #scaling
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  #clustering
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object),
                          seed.use = seed)
  p1 = ElbowPlot(seurat_object)
  seurat_object <- FindNeighbors(seurat_object,dims = 1:dims)
  seurat_object <- FindClusters(seurat_object, resolution = resolution)
  seurat_object <- RunUMAP(seurat_object, dims = 1:dims,seed.use = seed)
  p2 = DimPlot(seurat_object, reduction = "umap")
  return(list(object = seurat_object,plot = list(p1,p2)))
}

seurat_pipeline = function(seurat_object_raw,pc_dims,resolution,seed=42,
                           organism='human',nFeature_RNA_threshold=NA,mt_threshold=NA){
  qc_list = QC(seurat_object_raw, organism,nFeature_RNA_threshold,mt_threshold)
  seurat_object = qc_list[[1]];vlnplot = qc_list[[2]]
  preprocess_list = preprocess(seurat_object,pc_dims,resolution,seed)
  seurat_object = preprocess_list[[1]];elbowplot = preprocess_list[[2]][[1]];umapplot= preprocess_list[[2]][[2]]
  return(list(object = seurat_object,plot=list(vlnplot,elbowplot,umapplot)))
}

read_seurat_from_mtx = function(mtx_file){
  require(Matrix)
  require(stringr)
  project = str_remove(basename(mtx_file),'.mtx')
  pickle_file = str_replace(mtx_file,'mtx','pickle')
  print('reading mtx file')
  data= readMM(mtx_file)
  print('reading pickle file')
  pickle_data <-py_load_object(pickle_file)
  data = t(data)
  rownames(data) = pickle_data[[2]];colnames(data) = pickle_data[[1]]
  print('creating seurat object')
  object = CreateSeuratObject(data,project=project,min.cells = 3)
  return(object)
}

process_transcript_from_gene = function(gene_data, transcript_data){
  require(Seurat)
  keep =which(str_replace(rownames(transcript_data@assays$RNA),"-uncategorized.*$|-EN.*$|-novel.*$", "")%in%rownames(gene_data@assays$RNA))
  transcript_data = subset(x=transcript_data, cells = Cells(gene_data),features = rownames(transcript_data@assays$RNA)[keep])
  transcript_data@meta.data$seurat_clusters = gene_data@meta.data$seurat_clusters
  transcript_data@meta.data$labeled_clusters = Idents(gene_data)
  mat = transcript_data@assays$RNA$counts
  df = data.frame(genes = str_replace(rownames(mat),"-uncategorized.*$|-EN.*$|-novel.*$",""),
                  transcripts = rownames(mat))
  return(list(transcript_data=transcript_data, gene_transcript_name = df))
}

query_count_matrix = function(gene, cluster, transcript_data_list, remove_uncategorize = T, remove_novel=F){
  transcript_data = transcript_data_list[[1]];gene_transcript_name = transcript_data_list[[2]]
  if (is.na(cluster)==F){
    mat = t(as.matrix(transcript_data@assays$RNA$counts[which(gene_transcript_name$genes==gene),
                                                        which(transcript_data@meta.data$labeled_clusters==cluster)]))
  }else{
    mat = t(as.matrix(transcript_data@assays$RNA$counts[which(gene_transcript_name$genes==gene),]))
  }
  
  #deal with uncategorized
  colindex = which(colnames(mat)%in%paste0(gene,c('-uncategorized_novel','-uncategorized')))
  #deal with novel
  colindex2 =which(str_detect(colnames(mat),paste0(gene,'-novelIsoform'))==T)
  nr = nrow(mat);nc=ncol(mat)
  if (nc-length(colindex)<2){
    mat = NULL
  }else if (remove_uncategorize==T){
    if (length(colindex)>0){
      mat = mat[,-colindex]
    }
  }else{#merge
    if (length(colindex)>1){
      mat = cbind(mat[,-colindex], rowSums(mat[,colindex]));colnames(mat)[dim(mat)[2]] = paste0(gene,'-uncategorized')
    }
    if (length(colindex)==1){
      mat = cbind(mat[,-colindex], mat[,colindex]);colnames(mat)[dim(mat)[2]] = paste0(gene,'-uncategorized')
    }
  }
  nr = nrow(mat);nc=ncol(mat)
  if (nc-length(colindex2)<2){
    mat = NULL
  }else if (remove_novel==T){
    if (length(colindex2)>0){
      mat = mat[,-colindex2]
    }
  }else{mat = mat}
  #filter dropout
  if (is.null(mat)==F){
    if (sum(rowSums(mat)==0)>0){
      #print('removing drop out cells')
      mat = mat[-which(rowSums(mat)==0),]
    }
  }
  return(mat)
}
