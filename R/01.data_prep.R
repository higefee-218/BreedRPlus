#' Load PLINK files into a bigSNP object and optionally align phenotype
#'
#' @param bed_file Path to the .bed file.
#' @param bim_file Path to the .bim file.
#' @param fam_file Path to the .fam file.
#' @param pheno Optional data.frame of phenotypes (must contain id_col_name). Default NULL.
#' @param id_col_name Name of the column in pheno containing individual IDs.
#' @param impute_method The imputation method. Default is "mode".
#' @param backingfile Path for the 'bigSNP' backing file.
#'
#' @return A list containing:
#' \item{snp_obj}{The bigSNP object.}
#' \item{pheno}{Phenotype data aligned to genotypes (if provided).}
#'
#' @importFrom bigsnpr snp_readBed snp_attach snp_fastImputeSimple
#' @export
load_plink <- function(bed_file,
                       bim_file,
                       fam_file,
                       pheno = NULL,
                       id_col_name = NULL,
                       impute_method = "mode",
                       backingfile = tempfile()) {
  
  message("Starting PLINK load process...")
  
  if (!file.exists(bed_file)) stop(".bed file not found: ", bed_file)
  if (!file.exists(bim_file)) stop(".bim file not found: ", bim_file)
  if (!file.exists(fam_file)) stop(".fam file not found: ", fam_file)
  
  # --- Read .bed into bigSNP ---
  message("--> Reading .bed file...")
  rds_path <- bigsnpr::snp_readBed(bed_file, backingfile = backingfile)
  snp_obj <- bigsnpr::snp_attach(rds_path)
  
  # --- Impute missing genotypes ---
  message("--> Imputing missing genotypes using method: '", impute_method, "'...")
  bigsnpr::snp_fastImputeSimple(snp_obj$genotypes, method = impute_method)
  
  # --- Fix: store bed file path for later use ---
  snp_obj$bedfile <- bed_file
  
  # --- Optional: align phenotypes if provided ---
  aligned_pheno <- NULL
  if (!is.null(pheno)) {
    if (is.null(id_col_name)) stop("id_col_name must be provided when pheno is supplied.")
    
    geno_ids <- snp_obj$fam$sample.ID
    message("--> Aligning phenotype to genotype IDs...")
    
    # Keep only individuals present in genotype
    aligned_pheno <- pheno[pheno[[id_col_name]] %in% geno_ids, , drop = FALSE]
    
    # Match order to genotype
    geno_match <- match(aligned_pheno[[id_col_name]], geno_ids)
    aligned_pheno <- aligned_pheno[order(geno_match), , drop = FALSE]
    
    # Sanity check
    if (!all(aligned_pheno[[id_col_name]] == geno_ids[geno_match])) {
      warning("Phenotype and genotype IDs are not perfectly aligned after matching.")
    } else {
      message("Phenotype successfully aligned to genotype order.")
    }
  }
  
  message("--- PLINK Load Complete! ---")
  
  return(list(
    snp_obj = snp_obj,
    pheno = aligned_pheno
  ))
}
