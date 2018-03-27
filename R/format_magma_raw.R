
# Christopher Medway
# Requires that correlations between genes have already been calculated 
# using twas_matrix_for_magma.R 


format_twas_raw_file <- function(matrix_dir, pos_file, twas_file, sampleSize = 100000, nparam = 100, symbol2EntrezId) {
  
  matrix_files <- list.files(matrix_dir, pattern = ".csv", full.names = T)
  matrix_names <- basename(gsub(matrix_files, pattern = ".csv", replacement = ""))
  pos_df       <- read.table(pos_file, header=T, stringsAsFactors = F)
  twas_results <- read.table(twas_file, header = T, stringsAsFactors = F)
  symbol2EntrezId <- read.table(symbol2EntrezId, header = F)
  
  # filter pos file to include only genes with i) entrezid and ii) valid twas statistic
  pos_df <- filter_pos_file(pos_df, twas_results, symbol2EntrezId)
  
  # order pos file
  pos_df <- pos_df[order(pos_df$CHR, pos_df$P0, pos_df$P1),]
  
   #there should be a matrix file for every row in pos file.
  if (file.exists("./output/twagma_raw/twagma.missing")) {file.remove("./output/twagma_raw/twagma.missing")}
  if (!(all(pos_df$ID %in% matrix_names))) {
    warning(paste0("FILE NOT AVAILABLE FOR ALL TWAS GENES"))
    write.table(pos_df[!(pos_df$ID %in% matrix_names),], file = "./output/twagma_raw/twagma.missing", row.names = F, col.names = F, quote = F)
  }

  # loop over pos_file
  returned <- loop_over_files(pos=pos_df, symbol2EntrezId, twas_results, nparam, sampleSize, matrix_dir)
  
  qc <- lapply(seq(length(returned)), function(x) returned[[x]][["QC"]])
  qc <- do.call(rbind, qc)
  qc_check(qcObject = qc, pos_df = pos_df)
  
  raw <- lapply(seq(length(returned)), function(x) { raw <- returned[[x]][["RAW"]]})
  i <- unlist(lapply(raw, function(x) {length(x) > 1}))
  if (file.exists("./output/twagma_raw/twagma.raw")) {file.remove("./output/twagma_raw/twagma.raw")}
  lapply(raw[i], function(x) {
    write.table(x, file = "./output/twagma_raw/twagma.raw", append = T, quote = F, row.names = F, col.names = F, sep = " ")
  })
  
  covar <- lapply(seq(length(returned)), function(x) { covar <- returned[[x]][["COV"]]})
  i <- unlist(lapply(raw, function(x) {length(x) > 1}))
  if (file.exists("./output/twagma_raw/twagma.covar")) {file.remove("./output/twagma_raw/twagma.covar")}
  lapply(covar[i], function(x) {
    write.table(x, file = "./output/twagma_raw/twagma.covar", append = T, quote = F, row.names = F, col.names = F, sep = " ")
  })
  
  return(qc)
}

filter_pos_file <- function(pos_df, twas_results, symbol2EntrezId) {
  
  pos_df$EntrezId <- symbol2EntrezId[match(pos_df$ID, symbol2EntrezId$V2),1]
  pos_df$twas_z   <- twas_results[match(pos_df$ID, twas_results$ID),"TWAS.Z"]
  pos_df$twas_p   <- twas_results[match(pos_df$ID, twas_results$ID),"TWAS.P"]
  pos_df <- pos_df[!(is.na(pos_df$twas_z)),]
  pos_df <- pos_df[complete.cases(pos_df[,c("twas_z","EntrezId")]),]
  # remove duplicate entrezid 
  pos_df <- pos_df[!(duplicated(pos_df$EntrezId)),]
  # renormalise twas
  #pos_df$twas_z <- (pos_df$twas_z - mean(pos_df$twas_z)) / sd(pos_df$twas_z)
  return(pos_df)
}


loop_over_files <- function(pos, symbol2EntrezId, twas_results, nparam, sampleSize, matrix_dir) {
  
  # loop over rows of pos_file
  out <- lapply( seq(dim(pos)[1]), function(x) {
    
    # initialise empty df tp store qc info
    qc <- as.data.frame(matrix(nrow = 1, ncol = 6))
    
    row   <- pos[x,]
    gene  <- row[["ID"]]
    chr   <- row[["CHR"]]
    start <- row[["P0"]]
    stop  <- row[["P1"]]
    
    qc["V1"] <- row["CHR"]
    qc["V2"] <- gene
      
    entrezid    <- pos[match(gene, pos$ID),"EntrezId"]
    twas_z      <- pos[match(gene, pos$ID),"twas_z"]
    twas_p      <- pos[match(gene, pos$ID),"twas_p"]
    
    # calculate probit transformed p-value - this is how magma does it
    twas_probit <- qnorm(twas_p, lower.tail = F)
    # p-values = 1 generate "-Inf" - conver to mean
    twas_probit[twas_probit == "-Inf"] <- -3.09 # because -Inf will break MAGMA
        
    # read ld file if exists
    file <- paste0(matrix_dir,"/",gene,".csv")
    if (file.exists(file)) {
      df <- read.table(file, header = T, stringsAsFactors = F, sep = ",")
      
      nsnps <- df$GENE1_NSNPS[1]
      model <- df$GENE1_MODEL[1]
      
      # check gene has snp weights
      if (df$GENE1_NSNPS[1] > 0) {
        
        # if upstream gene(s) exist, remove any upstream genes that are not in pos file
        if ("GENE2" %in% names(df) && sum(!(is.na(df$GENE2))) > 0) {
              df <- df[df$GENE2 %in% pos$ID,]
            }
            
            # after removing invalid upstream genes, check gene has valid upstream genes remaining. This does not warrent elimination of the index gene. 
            # Just because it has no upstream genes, it maybe upstream of another gene
            if ((("COR" %in% names(df))) && sum(!(is.na(df$COR))) > 0) {
              
              validUs <- df[!(is.na(df$COR)),]
              qc$V3 <- "VALID"
              qc$V4 <- dim(validUs)[1]
              qc$V5 <- paste(validUs$GENE2, collapse = ",")
              qc$V6 <- paste(validUs$COR, collapse = ",")
              
              out <- data.frame(
                entrezid,
                #gene,
                chr,
                start,
                stop,
                nsnps,
                nparam,
                as.integer(sampleSize),
                abs(twas_z),
                #twas_probit,
                paste(rev(abs(validUs$COR)), collapse = " ") # using the absolute correlation
              )
              
            } else {
              qc$V3 <- "VALID"
              qc$V4 <- 0
              
              out <- data.frame(
                entrezid,
                #gene,
                chr,
                start,
                stop,
                nsnps,
                nparam,            
                as.integer(sampleSize),
                #twas_probit
                abs(twas_z)
              )
            }
        # covariate file
        lasso <- as.integer(model == "lasso")
        enet  <- as.integer(model == "enet")
        blup  <- as.integer(model == "blup")
        bslmm <- as.integer(model == "bslmm")
        
        cov <- data.frame(
          "ID" = entrezid,
          "NSNPS" = nsnps,
          "isLasso" = lasso,
          "isEnet" = enet,
          "isBlup" = blup,
          "isBslmm" = bslmm,
          "TWAS.Z" = twas_z,
          "ABS.TWAS.Z" = abs(twas_z),
          "PROBIT.TWAS.Z" = twas_probit,
          "TWAS.P" = twas_p,
          "GENE" = gene
          )
            
          } else {qc$V3 <- "NO_SNP_WEIGHTS"; out <- ""; cov = ""}
        } else {qc$V3 <- "NO_CORR_FILE"; out <- ""; cov = ""}
    return(list("RAW" = out, "QC" = qc, "COV" = cov))
  })
}


qc_check <- function(qcObject, pos_df) {
  
  l <- split(qcObject, qcObject$V3)
  
  lapply(seq(dim(l[["VALID"]])[1]), function(row) {
    
    r    <- l[["VALID"]][row,]
    gene <- r$V2
    n_us <- r[[4]]
    genes_us <- unlist(stringr::str_split(r[5], ","))
    
    if (n_us == 0) {
    } else {
      rows_us <- l[["VALID"]][(row - n_us):(row - 1),]
      if(!(dim(rows_us)[1] == n_us)){stop("DIFFERENT NUMBERS")}
      if(!(all(rows_us$V2 == genes_us))) {stop(paste0(gene, ": CALCULATED UPSTREAM GENES DOES NOT MATCH AVAILABLE UPSTREAM GENES!!"))}
      if(!(all(genes_us %in% pos_df$ID))) {stop(paste0(genes_us, " NOT ALL HAVE UPSTREAM GENES ARE VALID"))}
      if (!(all(rows_us$V3 == "VALID"))) {stop(paste0(gene," NOT ALL UPSTREAM GENES ARE VALID"))}
      }
    
    if(!(gene %in% pos_df$ID)) {stop(paste0(gene, "INDEX GENE NOT VALID"))}
    
    
  })
}


require("optparse")

option_list <- list(
  make_option(c("-c","--correlation_files_dir"), type = "character", default = NULL, help = "files containing gene-gene correlations"),
  make_option(c("-p","--fusion_pos_file"), type = "character", default = NULL, help = "containing gene coordinates (hg19)"),
  make_option(c("-t","--twas_results_file"), type = "character", default = NULL, help = "twas results file"),
  make_option(c("-n","--samplesize"), type = "integer", default = 10000, help = "integer givingsample number"),
  make_option(c("-u","--number_parameters"), type = "integer", default = 100, help = "integer giving parameter number"),
  make_option(c("-m","--symbol2entrez"), type = "character", default = NULL, help = "filename")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (!(dir.exists("./output/twagma_raw"))) {dir.create("./output/twagma_raw")}

x <- format_twas_raw_file(
  matrix_dir      =  opt$correlation_files_dir,
  pos_file        =  opt$fusion_pos_file,
  twas_file       =  opt$twas_results_file,
  symbol2EntrezId =  opt$symbol2entrez,
  sampleSize      =  opt$samplesize,
  nparam          =  opt$number_parameters  
)
