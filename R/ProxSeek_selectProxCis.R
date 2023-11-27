#' ProxSeek_selectProxCis
#'
#' This function select cis-eQTLs that lie within promoter and promoter interaction region (P-PIR) based on tissue/cell type Hi-C data.
#' It reads cis-eQTL and cell type P-PIR files, performs filtering and annotation, 
#' and returns an annotated data frame.
#'
#' @param ciseqtl.file A string specifying the path to the cis-eQTL file.
#'                     Default value is NA. It's recommended to provide the actual file path for specific analyses.
#' @param cell.type A string specifying the cell type. This is used to construct the file name of the cell type-specific P-PIR file.
#' @param n.cluster An integer specifying the number of clusters to use in parallel processing. Default is 4.
#'                  This should be set according to the hardware capabilities of the system.
#'
#' @return An annotated data frame with the results of the cis-eQTL and Hi-C interaction processing.
#'
#' @examples
#' # Example usage:
#' \dontrun{
#' demo.ciseqtl.file<-system.file("data", "ciseqtl.csv", package = "ProxSeek")
#' result <- ProxSeek_selectProxCis(ciseqtl.file = demo.ciseqtl.file, P_PIR_folder = /path/to/downloaded/PIRfiles/,
#' cell.type = "monocyte", n.cluster =4)
#' }
#' 
#'
#' @import data.table
#' @import dplyr
#' @import plyr
#' @import stringr
#' @import foreach
#' @import doSNOW
#' @import parallel
#' @export



ProxSeek_selectProxCis<-function(ciseqtl.file=NA, P_PIR_folder=NA, cell.type=NA, n.cluster = 4){
  ###### 1. read cis-eqtl tsv file
  # ciseqtl.file<-"/media/xingewang/JRlab10TB/Project_Space/Proxi_HIC_Validation/Prepare_CellTypePPIR/ciseqtl.tsv"
  cis.eqtl<-fread(ciseqtl.file)
  
  if (!file.exists(ciseqtl.file)) {
    stop("Cis-eQTL file not found: ", ciseqtl.file)
  }
  
  # Attempt to read the file
  tryCatch({
    cis.eqtl <- fread(ciseqtl.file)
  }, error = function(e) {
    stop("Failed to read cis-eQTL file: ", e$message)
  })
  
  # Validate the file content
  requiredColumns <- c("snps", "gene","Symbol","chr","pos")
  if (!all(requiredColumns %in% colnames(cis.eqtl))) {
    stop("Cis-eQTL file does not contain the required columns.")
  }
  
  ###### 2. read cell.type P-PIR file
  # cell.type<-"monocyte"
  cell.type.file<-paste0("P_PIR_", cell.type,".csv")
  cell.type.file<-paste0(P_PIR_folder,"/",cell.type.file)

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
  
  # Validate the file content
  requiredColumns <- c("P_chr","P_start", "P_end", "P_gene","PIR_chr","PIR_start", "PIR_end", "PIR_gene")
  if (!all(requiredColumns %in% colnames(cell.hic))) {
    stop("Cell type P-PIR file does not contain required columns.")
  }
  
  ##### 3. report genes overlapped 
  bait.arr<-unique(as.character(str_split(cell.hic$P_gene, ";", simplify = T)))
  bait.arr<-bait.arr[bait.arr!=""]
  bait.arr<-str_split(bait.arr, "[.]", simplify = T)[,1]
  bait.arr<-unique(bait.arr)
  hic.gene.n<-length(bait.arr) 
  cat(paste0("P-PIR file reported ", hic.gene.n, " genes!"))
  
  sel.genes<-intersect(unique(cis.eqtl$Symbol), bait.arr)
  share.gene.n<-length(sel.genes)
  cat(paste0("Shared number of genes with cis-eQTL report ", share.gene.n))
  message("")
  
  # 4. Subset cis-eqtl with genes report
  eqtl.df.sub<-cis.eqtl %>% filter(Symbol %in% sel.genes)
  eqtl.df.remain<-cis.eqtl %>% filter(!Symbol %in% sel.genes)

  # 5. annotate the P PIR for each SNP
  ##### get all the gene names into each row and split gene symbol and remove version
  kk<-str_split(cell.hic$P_gene, ";", simplify = T)
  kk.new<-gsub("\\..*", "", kk)
  head(kk.new)
  
  ######## loop each gene and subset eQTLs 
  g.all<-unique(eqtl.df.sub$Symbol)
  g.all<-g.all[g.all!=""]
  
  cl <- makeCluster(n.cluster)
  registerDoSNOW(cl)
  
  iterations <- length(g.all)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  tmp<-foreach(gi=g.all,.combine = 'rbind', .verbose = T, .options.snow = opts) %dopar% {
    ### subset eqtl for this gene
    gi.eqtl<-eqtl.df.sub[eqtl.df.sub$Symbol==gi,]
    ### get promoter and PIR region
    sel.idx<-which(rowSums(kk.new == gi)!=0)
    gi.hic.df<-cell.hic[sel.idx,]
    # loop each snps for this gene
    in.pir<-sapply(gi.eqtl$pos, function(x){
      x %in% unlist(Map(`:`, gi.hic.df$PIR_start, gi.hic.df$PIR_end))})
    snp.idx.pir<-which(in.pir == 1)
    
    in.prm <-sapply(gi.eqtl$pos, function(x){
      x %in% unlist(Map(`:`, gi.hic.df$P_start, gi.hic.df$P_end))})
    snp.idx.prm<-which(in.prm == 1)
    
    #### within snp 
    gi.eqtl$Within<-"non_PPIR"
    gi.eqtl$Within[snp.idx.pir]<-"PIR"
    gi.eqtl$Within[snp.idx.prm]<-"Promoter"
    return(gi.eqtl)
  }
  
  stopCluster(cl) 
  
  eqtl.df.remain$Within<-"No_HiC_Support"
  cis.eqtl.annotate<-rbind(tmp, eqtl.df.remain)
  return(cis.eqtl.annotate)
}