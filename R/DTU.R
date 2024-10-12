

perform_mle = function(X){
  require(rlang)
  require(extraDistr)
  neg_log_likelihood <- function(alpha, X) {
    if(any(alpha <= 0)) return(Inf)
    ll <- sum(ddirmnom(X, rowSums(X), alpha, log = TRUE))
    return(-ll)
  }
  start_values1 <- rep(1, ncol(X))
  start_values2 <- colSums(X)/sum(X)
  m = try(suppressWarnings({optim(
    par = start_values1, 
    fn = neg_log_likelihood, 
    X = X, 
    method = "L-BFGS-B",  
    lower = rep(0.0001, ncol(X))  
  )}))
  if (inherits(m, "try-error")) {
    m = try(suppressWarnings({optim(
      par = start_values2, 
      fn = neg_log_likelihood, 
      X = X, 
      method = "L-BFGS-B",  
      lower = rep(0.0001, ncol(X))  
    )}))
  }
  if (inherits(m, "try-error")) {
    m = NULL
    print('MLE failed to converge')
  }
  return(m)
}

bootstrap_dispersion = function(X){
  pars = list()
  for (i in 1:1000){
    pars[[i]]=perform_mle(X[sample(1:nrow(X),nrow(X),replace = T),])$par 
  }
  pars = do.call(rbind,pars)
  return(quantile(1/rowSums(pars),c(0.025,0.5,0.975)))
}


LRT_test = function(X1, X2, group_novel = FALSE, bootstrap = FALSE){
  require(dplyr)
  if (sum(X1)==0|sum(X2)==0){
    return(NULL)
  }
  if (group_novel){
    novel_cols1 <- grep("novel", colnames(X1))
    novel_cols2 <- grep("novel", colnames(X2))
    if(length(novel_cols1) > 0){
      novel_isoform_sum1 <- rowSums(X1[, novel_cols1, drop = FALSE])
      novel_isoform_name <- str_remove(colnames(X1)[novel_cols1][1], '[_-]\\d+')
      X1 <- X1[, -novel_cols1, drop = FALSE]
      X1 <- cbind(X1, novel_isoform_sum1)
      colnames(X1)[ncol(X1)] <- novel_isoform_name
    }
    if (length(novel_cols2) > 0){
      novel_isoform_sum2 <- rowSums(X2[, novel_cols2, drop = FALSE])
      novel_isoform_name <- str_remove(colnames(X2)[novel_cols2][1], '[_-]\\d+')
      X2 <- X2[, -novel_cols2, drop = FALSE]
      X2 <- cbind(X2, novel_isoform_sum2)
      colnames(X2)[ncol(X2)] <- novel_isoform_name
    }
  }
  #combine rare isoforms
  isoform_top = union(colnames(X1)[colSums(X1)/sum(X1)>=0.05],colnames(X2)[colSums(X2)/sum(X2)>=0.05])
  top=length(isoform_top)
  #add isoforms only show up in one group
  add1 = isoform_top[isoform_top%in%colnames(X1)==F]
  add2 = isoform_top[isoform_top%in%colnames(X2)==F]
  if (length(add1)>0){
    X1_add = matrix(0,ncol = length(add1),nrow = nrow(X1),dimnames = list(rownames(X1),add1))
    X1 = cbind(X1,X1_add)}
  if (length(add2)>0){
    X2_add = matrix(0,ncol = length(add2),nrow = nrow(X2),dimnames = list(rownames(X2),add2))
    X2 = cbind(X2,X2_add)}
  if (top>0){
    X1_ = X1[,colnames(X1)%in%isoform_top==F]
    
    if (is.null(dim(X1_))){
      X1 = cbind(X1[,isoform_top], other = X1_)
      colnames(X1) = c(isoform_top, 'other')
    }else if (ncol(X1_)>0){
      X1 = cbind(X1[,isoform_top], other = rowSums(X1_))
      colnames(X1) = c(isoform_top, 'other')
    }else{
      X1=X1[,isoform_top]
    }
    
    X2_ = X2[,colnames(X2)%in%isoform_top==F]
    if (is.null(dim(X2_))){
      X2 = cbind(X2[,isoform_top], other = X2_)
      colnames(X2) = c(isoform_top, 'other')
    }else if (ncol(X2_)>0){
      X2 = cbind(X2[,isoform_top], other = rowSums(X2_))
      colnames(X2) = c(isoform_top, 'other')
    }else{
      X2=X2[,isoform_top]
    }
  }
  
  if (ncol(X1)<ncol(X2)){
    X1 = cbind(X1,other=0)
  }
  if (ncol(X2)<ncol(X1)){
    X2 = cbind(X2,other=0)
  }
  if (sum(rowSums(X1))>=20 & sum(rowSums(X2))>=20){
    X = rbind(X1, X2)
    X = X[rowSums(X)>0,]
    X1 = X1[rowSums(X1)>0,]
    X2 = X2[rowSums(X2)>0,]
    #start testing-----gene level
    H0 = perform_mle(X);if (is.null(H0)){return(NULL)}
    H1 = perform_mle(X1);if (is.null(H1)){return(NULL)}
    H2 = perform_mle(X2);if (is.null(H2)){return(NULL)}
    
    chi = 2*(H0$value-H1$value-H2$value);p= pchisq(chi,ncol(X),lower.tail = F)
    #testing----transcript level
    p_transcript = c()
    if (ncol(X1)>2){
      for (i in 1:ncol(X1)){
        X1_temp = cbind(X1[,i],rowSums(X1[,-i]))
        X2_temp = cbind(X2[,i],rowSums(X2[,-i]))
        X_temp = rbind(X1_temp, X2_temp)
        H0_ = perform_mle(X_temp);if (is.null(H0_)){return(NULL)}
        H1_ = perform_mle(X1_temp);if (is.null(H1_)){return(NULL)}
        H2_ = perform_mle(X2_temp);if (is.null(H2_)){return(NULL)}
        if (is.null(H0_)|is.null(H1_)|is.null(H2_)){
          p_transcript= c(p_transcript,NA)
        }else{
          chi = 2*(H0_$value-H1_$value-H2_$value)
          p_transcript= c(p_transcript,pchisq(chi,ncol(X),lower.tail = F)) 
        }
      }
    }else{
      p_transcript = rep(p,ncol(X1))
    }
    #alpha
    df1 = as.data.frame(setNames(H1$par,colnames(X1)))%>%mutate(isoforms = rownames(.))
    colnames(df1)[1] = 'alpha1'
    df2 = as.data.frame(setNames(H2$par,colnames(X2)))%>%mutate(isoforms = rownames(.))
    colnames(df2)[1] = 'alpha2'
    df_gene = data.frame(isoforms = colnames(X1))%>%left_join(df1)%>%left_join(df2)
    alpha0_1 = sum(df_gene$alpha1);alpha0_2 = sum(df_gene$alpha2)
    df_gene = df_gene%>%mutate(TU1 = alpha1/alpha0_1,TU2 = alpha2/alpha0_2, 
                               TU_diff = abs(TU1-TU2))%>%
      mutate(TU_var1 =TU1*(1-TU1)/(1+alpha0_1) ,
             TU_var2 = TU2*(1-TU2)/(1+alpha0_2))%>%
      mutate(dispersion1 = 1/(alpha0_1+1),dispersion2 = 1/(alpha0_2+1))
    if (bootstrap==T){
      dis1 = bootstrap_dispersion(X1)
      dis2 = bootstrap_dispersion(X2)
      df_gene = df_gene%>%mutate(dispersion1_median = dis1[2], dispersion1_CI0 = dis1[1],dispersion1_CI1 = dis1[3],
                                 dispersion2_median = dis2[2], dispersion2_CI0 = dis2[1],dispersion2_CI1 = dis2[3])
    }
    df_topisoform=df_gene%>%filter(isoforms!='other')%>%filter(TU1==max(TU1)|TU2==max(TU2))
    isoform_switch = ifelse(nrow(df_topisoform)==1,F,T)
    if (isoform_switch==T){
      isoform_switch_ES =abs((df_topisoform$TU1[1]-df_topisoform$TU2[1])+(df_topisoform$TU2[2]-df_topisoform$TU1[2]))
    }else{
      isoform_switch_ES = abs(df_topisoform$TU1-df_topisoform$TU2)}
    df_gene = df_gene%>%mutate(isoform_switch=isoform_switch,isoform_switch_ES=isoform_switch_ES,
                               p_DTU_gene = p,p_transcript=p_transcript)
  }else{df_gene = NULL}
  return(df_gene)
}




#' Differential gene expression analysis using Wilcoxon Rank Sum test
#'
#' @param X1_gene cell(row) x gene(column) count matrix for cell population 1. 
#' @param X2_gene cell(row) x gene(column) count matrix for cell population 2. 
#' @param epsilon a small number to add on gene count matrix for computational purpose.
#' @param ncores number of cores to use for parallel computing
#' @return result dataframe.pct1 and pct2: cell percentage expressing the gene; logFC: log fold change; p_gene: p-value for the test.
#' @examples
#' df_gene = scotch_gene(sample8_CD4_gene, sample8_others_gene, epsilon=0.01,ncores=10)
#' @export
scotch_gene = function(X1_gene, X2_gene, epsilon=0.01,ncores=10){
  require(doParallel)
  require(foreach)
  registerDoParallel(cores=ncores)
  logFC = c();p=c();pct1=c();pct2=c()
  X1_gene=X1_gene+epsilon;X2_gene=X2_gene+epsilon
  genes = intersect(colnames(X1_gene),colnames(X2_gene))
  results = foreach (i =1:length(genes),.combine = 'rbind')%dopar%{
    wtest = wilcox.test(X1_gene[,genes[i]],X2_gene[,genes[i]])
    pct1 = sum(X1_gene[,genes[i]]>epsilon)/nrow(X1_gene)
    pct2= sum(X2_gene[,genes[i]]>epsilon)/nrow(X2_gene)
    logFC =log2(mean(X1_gene[,genes[i]])/mean(X2_gene[,genes[i]]))
    p = wtest$p.value
    c(pct1,pct2,logFC,p)
  }
  colnames(results) = c('pct1','pct2','logFC','p_gene')
  df_gene = cbind(genes = genes,as.data.frame(results))%>%
    mutate(pct_diff=pct1-pct2,p_gene_adj = p.adjust(p_gene))
  return(df_gene)
}


LRT_by_gene = function(gene_name,gene_transcript_df1,gene_transcript_df2, X1_transcript, X2_transcript, group_novel, bootstrap){
  require(dplyr)
  transcripts1_gene = gene_transcript_df1%>%filter(genes==gene_name)%>%pull(transcripts)
  transcripts2_gene = gene_transcript_df2%>%filter(genes==gene_name)%>%pull(transcripts)
  X1 = X1_transcript[,transcripts1_gene]
  X2 = X2_transcript[,transcripts2_gene]
  if (is.numeric(nrow(X1)) && is.numeric(nrow(X2)) && nrow(X1) >= 20 && nrow(X2) >= 20){
    out = LRT_test(X1,X2, group_novel, bootstrap)
  }else{out = NULL}
  if (is.null(out)==F){
    out$gene = gene_name
  }
  return(out)
}


#' Differential transcript usage analysis between two cell populations
#'
#' @param gene_transcript_df1 dataframe storing gene and transcript names for the transcript count matrix: X1_transcript. 
#' @param gene_transcript_df2 dataframe storing gene and transcript names for the transcript count matrix: X2_transcript. 
#' @param X1_transcript cell(row) x transcript(column) count matrix for cell population 1. 
#' @param X2_transcript cell(row) x transcript(column) count matrix for cell population 2. 
#' @param ncores number of cores to use for parallel computing
#' @param group_novel whether to group all novel isoforms together when testing DTU genes, default is FALSE
#' @param bootstrap whether to performe bootstraping to generate CI for alpha, default is FALSE
#' @return result dataframe.
#' @examples
#' df_transcript = scotch_transcript(gene_transcript_CD4_df,gene_transcript_CD8_df, sample8_CD4_transcript, sample8_CD8_transcript, 
#' ncores=10)
#' @export
scotch_transcript = function(gene_transcript_df1,gene_transcript_df2, X1_transcript, X2_transcript,ncores=10, group_novel=FALSE, bootstrap = FALSE){
  require(doParallel)
  require(foreach)
  registerDoParallel(cores=ncores)
  genes = union(gene_transcript_df1$genes,gene_transcript_df2$genes)
  
  print('start isoform analysis')
  results = foreach (i =1:length(genes),.combine = 'rbind',.errorhandling = 'pass')%dopar%{
    tryCatch({
      gene_name = genes[i]
      LRT_by_gene(gene_name, gene_transcript_df1,gene_transcript_df2, X1_transcript, X2_transcript, group_novel, bootstrap)
    },error = function(e){
      cat("Error at i =", i, "with gene:", gene_name, "Error message:", e$message, "\n")
      NULL
    })
  }  
  results = results[!sapply(results, is.null), ]
  if (ncol(results) > 15) {
    results = results[, c('gene', colnames(results)[1:15])]
  }
  print('completed')
  results=results%>%mutate(p_transcript_adj = p.adjust(p_transcript,'fdr'))
  results = left_join(results,results%>%dplyr::select(gene,p_DTU_gene)%>%unique()%>%
                        mutate(p_DTU_gene_adj = p.adjust(p_DTU_gene,'fdr')))
  return(results)
}


