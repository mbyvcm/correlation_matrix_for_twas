
# Christopher Medway
# script calculates correlation between genes using
# a specific set of SNPs (whereas MAGMA uses all SNPs in the cis region)
# USEAGE:
# Rscript --vanilla twas_matrix_for_magma <chr(1-22)> 

## MUST SET FOLLOWING ARGUMENTS:

# these files / directory are included with weights
# downloaded on the FUSION website
weights_profile <- "<./*.profile file>"
weights_pos     <- "<./*.pos>" # must be ordered by chr and start coordinate 
weights_dir     <- "<./dir/ directory>"
ldref_dir       <- "~/Resources/TWAS/ldref/"

# also location of plink binary
plink           <- "plink2"


get_best_model <- function(weights_profile, weights_pos) {
# this function will fetch the best (significant with highest R2) model
# for each gene in TWAS 
  
  # contains model information per gene
  profiledf <- read.table(weights_profile, header=T, stringsAsFactors = F)
  
  # contains gene coordinates
  pos <- read.table(weights_pos, header = T, stringsAsFactors = F)

  # loop over each gene (one gene per row in .pos file)
  out <- apply(pos,1,function(row) {
    gene <- row[["ID"]]
    message(gene)
    
    # get gene info from files
    profile_row <- profiledf[match(gene,gsub(profiledf$id, pattern = "CMC.", replacement = "")),]
    bestModel   <- tools::file_path_sans_ext(names(which.min(profile_row[c("blup.pv","enet.pv","bslmm.pv","lasso.pv")])))
    coord       <- row[c("CHR","P0","P1")]
    chr         <- as.numeric(coord[["CHR"]])
    start       <- as.numeric(coord[["P0"]])
    stop        <- as.numeric(coord[["P1"]])
    
    # identify genes upstream by no more than 500kbp
    earlierGene <- pos[(pos[,"CHR"] %in% chr) & (pos[,"P0"] < start) & (pos[,"P0"] > (start - 5000000)),"ID"]
    # also include genes with same start site, but which are shorter
    earlierGene   <- c(earlierGene, pos[(pos[,"CHR"] %in% chr) & (pos[,"P0"] == start) & (pos[,"P1"] < stop),"ID"]) 
    earlierGene   <- earlierGene[!(earlierGene %in% gene)] # remove index gene
    
    l <- list("GENE" = gene, "MODEL" = bestModel, "EARLIER_GENE" = earlierGene, "COORD" = data.frame("CHR" = chr,"START" = start,"STOP" = stop))
    return(l)
  })
  
  names(out) <- unlist(lapply(out, function(x){ x[["GENE"]]}))
  return(out)
}



get_weights_for_model <- function(bestModels, weights_dir) {
# extract SNP weights for the best model for all genes  
  out <- lapply(names(bestModels), function(genename) {
    entry <- bestModels[[genename]]
    model <- entry[["MODEL"]]
    file <- paste0(weights_dir,"CMC.",genename,".wgt.RDat")
    if(file.exists(file)) {
      load(file)
      w <- wgt.matrix[,model]
      w <- w[!(is.na(w))]
      w <- w[w > 0]
      entry[["WEIGHTS"]] <- w
      return(entry)
    }
  })
  
  names(out) <- unlist(lapply(out, function(x){ x[["GENE"]]}))
  return(out)
}



get_ld_plink <- function(bestModelWeights) {
  
  genes <- names(bestModelWeights)
  
  lapply(genes, function(gene) {
    
    message(gene)
    
    gene_slot    <- bestModelWeights[[gene]]
    gene_weights <- gene_slot[["WEIGHTS"]]
    chr          <- gene_slot[["COORD"]][["CHR"]]
    
    # if there are no weights or no earlier genes - skip
    if ((length(gene_slot[["WEIGHTS"]]) > 0) & (length(gene_slot[["EARLIER_GENE"]]) > 0)) {
      
      # loop through earlier genes
      upstream_genes <- gene_slot[["EARLIER_GENE"]]
      gene_level_reults <- lapply(upstream_genes, function(usgene) {
        
        message(paste0("....upstream gene: ", usgene))
        
        usgene_weights <- bestModelWeights[[usgene]][["WEIGHTS"]]
        if (length(usgene_weights) > 0) {
          # write snps in upstream gene and query gene to file - these will be used to --extract in PLINK
        
        snps_extract <- unique(c(names(gene_weights), names(usgene_weights)))  
            
          write.table(
            snps_extract,
            file = paste0("../output/tmp_extract_file_",gene,"_",usgene,".txt"),
            row.names = F,
            quote = F,
            col.names = F
            )

          command <- paste0(
            plink, 
            " --bfile ", ldref_dir,"1000G.EUR.",chr,
            " --extract ../output/tmp_extract_file_",gene,"_",usgene,".txt",
            " --r inter-chr with-freqs",
            " --ld-snp-list ../output/tmp_extract_file_",gene,"_",usgene,".txt",
            " --out ../output/tmp_",gene,"_",usgene
            )
          
          
          system(command)
          
          
          ldfile <- read_plink_ld_file(gene = gene, usgene = usgene)
          ldfile <- gather_weights_to_matrix(ldfile = ldfile, usgene_weights = usgene_weights, gene_weights = gene_weights)
          ldfile[["COV"]] <- covariance_t1_t2(ldfile, t1 = names(gene_weights), t2 = names(usgene_weights))
          
          # calculate gene variance
          ldfile[["VAR"]] <- calculate_genes_variance(ldfile = ldfile, t1 = names(gene_weights), t2 = names(usgene_weights))
          
          # calculate corr
          ldfile[["CORR"]] <- calculate_correlation_t1_t2(ldfile)
          
          # delete PLINK files
          file.remove(paste0("../output/tmp_",gene,"_",usgene,".log"))
          file.remove(paste0("../output/tmp_",gene,"_",usgene,".ld"))
          file.remove(paste0("../output/tmp_",gene,"_",usgene,".nosex"))
          file.remove(paste0("../output/tmp_extract_file_",gene,"_",usgene,".txt"))
          
          return(
            data.frame(
              "GENE1" = gene,
              "GENE1_CHR" = bestModelWeights[[gene]][["COORD"]][["CHR"]],
              "GENE1_START" = bestModelWeights[[gene]][["COORD"]][["START"]],
              "GENE1_STOP" = bestModelWeights[[gene]][["COORD"]][["STOP"]],
              "GENE1_MODEL" = bestModelWeights[[gene]][["MODEL"]],
              "GENE1_NSNPS" = length(bestModelWeights[[gene]][["WEIGHTS"]]),
              "GENE2" = usgene,
              "GENE2_CHR" = bestModelWeights[[usgene]][["COORD"]][["CHR"]],
              "GENE2_START" = bestModelWeights[[usgene]][["COORD"]][["START"]],
              "GENE2_STOP" = bestModelWeights[[usgene]][["COORD"]][["STOP"]],
              "GENE2_MODEL" = bestModelWeights[[usgene]][["MODEL"]],
              "GENE2_NSNPS" = length(bestModelWeights[[usgene]][["WEIGHTS"]]),
              "COV" = ldfile[["COV"]],
              "VAR_1" = ldfile[["VAR"]][["T1"]],
              "VAR_2" = ldfile[["VAR"]][["T2"]],
              "COR" = ldfile[["CORR"]])
            )
        } else {
          
            df <- data.frame(
              "GENE1" = gene,
              "GENE1_CHR"   = bestModelWeights[[gene]][["COORD"]][["CHR"]],
              "GENE1_START" = bestModelWeights[[gene]][["COORD"]][["START"]],
              "GENE1_STOP"  = bestModelWeights[[gene]][["COORD"]][["STOP"]],
              "GENE1_MODEL" = bestModelWeights[[gene]][["MODEL"]],
              "GENE1_NSNPS" = length(bestModelWeights[[gene]][["WEIGHTS"]]),
              "GENE2"      = usgene,
              "GENE2_CHR" = NA,
              "GENE2_START" = NA,
              "GENE2_STOP" = NA,
              "GENE2_MODEL" = NA,
              "GENE2_NSNPS" = NA,
              "COV" = NA,
              "VAR_1" = NA,
              "VAR_2" = NA,
              "COR" = NA
                )
            return(df)
          }
      })
      
      out <- do.call(rbind, gene_level_reults)
      write.table(out, paste0("../output/",gene,".csv"), quote = F, sep = ",", row.names = F)
    } else {
      out <- data.frame(
        "GENE1" = gene,
        "GENE1_CHR"   = bestModelWeights[[gene]][["COORD"]][["CHR"]],
        "GENE1_START" = bestModelWeights[[gene]][["COORD"]][["START"]],
        "GENE1_STOP"  = bestModelWeights[[gene]][["COORD"]][["STOP"]],
        "GENE1_MODEL" = bestModelWeights[[gene]][["MODEL"]],
        "GENE1_NSNPS" = length(bestModelWeights[[gene]][["WEIGHTS"]]),
        "UPSTREAM GENES" = paste((gene_slot[["EARLIER_GENE"]]), collapse = ";")
      )
      write.table(out, paste0("../output/",gene,".csv"), quote = F, row.names = F, sep = ",")
    }
  })
  
  
  return(1)
}

calculate_correlation_t1_t2 <- function(ldfile) {
  
  # cov(T1, T2)/√(var(T1)var(T2))
  
  cor <- ldfile[["COV"]]/sqrt(ldfile[["VAR"]][["T1"]]*ldfile[["VAR"]][["T2"]])
  return(cor)
}


calculate_genes_variance <- function(ldfile, t1, t2) {
  
  out <- lapply(list(t1,t2), function(t) {
    
    # maf t= snp names
    maf_i <- ldfile[["MAF_A"]][t,t]
    maf_j <- ldfile[["MAF_B"]][t,t]
    
    var_maf_i <- 2*(maf_i*(1-maf_i))
    var_maf_j <- 2*(maf_j*(1-maf_j))
    maf <- sqrt(var_maf_i * var_maf_j)
    
    # weights
    weights <- ldfile[["WEIGHTS_A"]][t,t] * ldfile[["WEIGHTS_B"]][t,t]
    
    # LD
    ld <- ldfile[["LD"]][t,t]
    
    # only use rowSum if at least 2x2 matrix
    if (is.matrix(ld)) {
      var <- sum(rowSums(weights * ld * maf))
    } else {
      var <- sum(weights * ld * maf)
    }
  
    return(var)
  })
  
  names(out) <- c("T1","T2")
  return(out)
}


covariance_t1_t2 <- function(ldfile, t1, t2) {
  # calculate variance of MAFs
  maf_gene   <- ldfile[["MAF_A"]][t1,t2]
  maf_usgene <- ldfile[["MAF_B"]][t1,t2]
  
  # 
  var_maf_gene <- 2*(maf_gene*(1-maf_gene))
  var_maf_usgene <- 2*(maf_usgene*(1-maf_usgene))
  maf <- sqrt(var_maf_gene * var_maf_usgene)
  
  weights <- ldfile[["WEIGHTS_A"]][t1,t2] * ldfile[["WEIGHTS_B"]][t1,t2] 
  
  ld <- ldfile[["LD"]][t1,t2]
  
  # only use rowSum if at least 2x2 matrix
  if (is.matrix(ld)) {
    cov <- sum(rowSums(weights * ld * maf))
  } else {
    cov <- sum(weights * ld * maf)
  }

  return(cov)
  }


gather_weights_to_matrix <- function(ldfile, usgene_weights, gene_weights) {
  
  usgene_weights <- usgene_weights[match(colnames(ldfile[["LD"]]), names(usgene_weights))]
  gene_weights <- gene_weights[match(row.names(ldfile[["LD"]]), names(gene_weights))]
  
  weights <- c(usgene_weights, gene_weights)
  weights <- weights[match(row.names(ldfile[["LD"]]), names(weights))]
  
  m_weights_A <- matrix(rep(weights, length(weights)), byrow = F, nrow = length(weights))
  m_weights_B <- matrix(rep(weights, length(weights)), byrow = T, ncol = length(weights))
  
  colnames(m_weights_A) <- names(weights)
  colnames(m_weights_B) <- names(weights)
  
  row.names(m_weights_A) <- names(weights)
  row.names(m_weights_B) <- names(weights)
  
  ldfile[["WEIGHTS_A"]] <- m_weights_A
  ldfile[["WEIGHTS_B"]] <- m_weights_B
  return(ldfile)
}


read_plink_ld_file <- function(gene, usgene) {
  
  ldfile <- read.table(paste0("../output/tmp_",gene,"_",usgene,".ld"), header = T, stringsAsFactors = F)

  l <- split(ldfile, ldfile$SNP_A)
  m_ld <- do.call(rbind, lapply(l, function(x) x[,"R"]))
  colnames(m_ld) <- l[[1]][,"SNP_B"]
  m_ld <- m_ld[,row.names(m_ld)]
  
  # maf
  maf <- unique(ldfile[,c("SNP_A","MAF_A")])
  maf <- maf[match(row.names(m_ld),maf$SNP_A),]
  
  m_maf_A <- matrix(rep(maf$MAF_A, length(maf$SNP_A)), byrow = F, ncol = length(maf$SNP_A))
  m_maf_B <- matrix(rep(maf$MAF_A, length(maf$SNP_A)), byrow = T, ncol = length(maf$SNP_A))
  
  colnames(m_maf_A) <- colnames(m_ld)
  colnames(m_maf_B) <- colnames(m_ld)
  
  row.names(m_maf_A) <- row.names(m_ld)
  row.names(m_maf_B) <- row.names(m_ld)
  
  return(list("LD" = m_ld, "MAF_A" = m_maf_A, "MAF_B" = m_maf_B))
} 

chr <- commandArgs(trailingOnly = T)[1]

bestModels <- get_best_model(weights_profile, weights_pos = weights_pos)
bestModelWeights <- get_weights_for_model(bestModels = bestModels, weights_dir = weights_dir)
bestModelWeights <- bestModelWeights[unlist(lapply(bestModelWeights, function(x) x[["COORD"]][["CHR"]])) == chr]
get_ld_plink(bestModelWeights = bestModelWeights)
