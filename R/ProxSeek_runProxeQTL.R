

#' Run Proximity-based eQTL Analysis
#'
#' This function performs proximity-based eQTL analysis using SNP and gene expression data, 
#' integrating cell type-specific physical interaction regions (P-PIRs). It outputs significant 
#' Prox-eQTLs, considering cellular context.
#'
#' @usage ProxSeek_runProxeQTL(SNP_file_name, expression_file_name, snps_location_file_name,
#'        gene_location_file_name, covariates_file_name, cell.type = NA, chr = NA,
#'        output.file = "./test/ciseqtl.tsv", n.cluster = 4)
#'
#' @param SNP_file_name Path to the SNP genotype file.
#' @param expression_file_name Path to the gene expression file.
#' @param snps_location_file_name Path to the file with SNP locations.
#' @param gene_location_file_name Path to the file with gene locations.
#' @param covariates_file_name Path to the covariates file.
#' @param cell.type  Specific cell type to analyze.
#' @param chr integer. Chromosome to analyze.
#' @param output.file Output file path for the results.
#' @param n.cluster Number of clusters for analysis.
#' @param P_PIR_folder The folder where P-PIR dataset were downloaded.
#'
#' @details
#' This function integrates MatrixeQTL pipeline and P-PIR information to identify prox-eQTLs.
#'
#' @return A data frame with the results of the proximity-based eQTL analysis.
#'
#' @examples
#' # Example usage of ProxSeek_runProxeQTL
#' \dontrun{
#' result <- ProxSeek_runProxeQTL("path/to/SNP_file.txt", "path/to/expression_file.txt",
#'                                "path/to/snps_location_file.txt", "path/to/gene_location_file.txt",
#'                                "path/to/covariates_file.txt", cell.type = "monocyte", chr = "X")
#' }
#' @export
#' @import data.table
#' @import dplyr
#' @import plyr
#' @import stringr
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @import MatrixEQTL


ProxSeek_runProxeQTL<-function(SNP_file_name, expression_file_name,snps_location_file_name,
                               gene_location_file_name,covariates_file_name, cell.type = NA, chr = NA,
                               output.file = "output.csv", n.cluster =4, P_PIR_folder){
  if(FALSE){
    cell.type<-"Mon"
    SNP_file_name = system.file("data", "chr19_geno2.csv", package = "ProxSeek")
    expression_file_name = system.file("data", "GE2.csv", package = "ProxSeek")
    snps_location_file_name = system.file("data", "chr19_geno_loc2.csv", package = "ProxSeek")
    gene_location_file_name = system.file("data", "geneloc_TSS2.csv", package = "ProxSeek")
    covariates_file_name = system.file("data", "Covariates.csv", package = "ProxSeek")
    chr<-19
  }
  ##### 0. check chromosome id 
  chr.i<-chr
  chr.arr<-c(seq(1,22), "X", "Y")
  if(!chr.i %in% chr.arr){
    stop("Chromosome must be between 1-22 or X and Y!!!")
  }
  
  ##### 1. read p-pir file and get the furthest distance
  # cell.type<-"monocyte"
  cell.type.file<-paste0("P_PIR_", cell.type,".csv")
  cell.type.file<-system.file("data", cell.type.file, package = "ProxSeek")
  
  # Check if the file exists
  if (!file.exists(cell.type.file)) {
    stop("Cell type P-PIR file not found: ", cell.type.file)
  }
  
  # Attempt to read the file
  tryCatch({
    cell.hic <- fread(cell.type.file)
  }, error = function(e) {
    stop("Failed to read cell type P-PIR file: ", e$message)
  })
  
  #cell.hic<-as.data.frame(cell.hic)
  # Validate the file content
  requiredColumns <- c("P_chr","P_start", "P_end", "P_gene","PIR_chr","PIR_start", "PIR_end", "PIR_gene")
  if (!all(requiredColumns %in% colnames(cell.hic))) {
    stop("Cell type P-PIR file does not contain required columns.")
  }
  
  # get the furthest distance 
  cell.hic.i<-cell.hic %>% filter(P_chr == chr.i & PIR_chr == chr.i)
  dist.df<-ddply(cell.hic.i, .(P_gene), function(df){
    ### get max distance to TSS of the gene
    ### down stream 
    down.idx<-which(df$PIR_start - df$P_end >0)
    up.idx<-which(df$PIR_end - df$P_start<0)
    
    if(length(down.idx)!=0){
      down.dist<-df$PIR_end - df$P_start
      down.dist<-down.dist[down.idx]
      max.down.dist<-max(down.dist)
    }else{max.down.dist<-0}
    if(length(up.idx)!=0){
      up.dist<-df$P_end - df$PIR_start
      up.dist<-up.dist[up.idx]
      max.up.dist<-max(up.dist)
    }else{max.up.dist<-0}
    
    return(max(max.down.dist,max.up.dist))
  })
  max.cis.dist<-max(dist.df$V1)
  
  ##### 2. run MatrixeQTL 
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Genotype file name
  SNP_file_name = SNP_file_name
  snps_location_file_name = snps_location_file_name
  
  # Gene expression file name
  expression_file_name = expression_file_name
  gene_location_file_name = gene_location_file_name
  
  # Covariates file name
  # Set to character() for no covariates
  covariates_file_name = covariates_file_name
  
  # Output file name
  output_file_name_cis = tempfile();
  output_file_name_tra = tempfile();
  
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = 0.05;
  pvOutputThreshold_tra = 0;
  
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
  
  # Distance for local gene-SNP pairs
  cisDist = max.cis.dist;
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  unlink(output_file_name_tra);
  unlink(output_file_name_cis);
  
  ### add required column names c("snps", "gene","Symbol","chr","pos")
  snp.loc.df<-fread(snps_location_file_name)
  
  result.df<-me$cis$eqtls
  result.df$chr<-snp.loc.df$chr[match(result.df$snps, snp.loc.df$snp)]
  result.df$pos<-snp.loc.df$pos[match(result.df$snps, snp.loc.df$snp)]
  result.df$Symbol<-result.df$gene
  
  result.df<-result.df %>% filter(FDR<=0.05)
  
  
  if(!is.na(output.file)){
    fwrite(result.df, output.file)
  }else{
    stop("Must have an output file!")
  }
  
  result<-ProxSeek_selectProxCis(ciseqtl.file = output.file, P_PIR_folder = P_PIR_folder, cell.type = cell.type, n.cluster = n.cluster)
  return(result)
}
