#johaGL 2025

library(GSVA)
set.seed(1) 

# Running GSVA on processed metabolites matrix, as it is performed by SPATA2


############ minimal test 

run_minimal_test <- function(){
  cat("*-*-*\n DEBUG purpose: running fast example (as no options specified)\n")
  X = matrix(data=rnorm(10*5), nrow=10, ncol=5)
  rownames(X) = as.character(seq(10))
  colnames(X) = as.character(seq(5))
  sets_list = list("set1"=c("5","7","9","6"),
                   "set2"=c("8","9","10"),
                   "set3"=c("1", "10","7"))
  print(sets_list)
  cat("\n printing the fast example (a matrix of scores per set)\n")
  gsvaPar <- gsvaParam(X, sets_list, minSize=2)  
  gsva_scores_matrix <- gsva(gsvaPar, verbose=TRUE)
  print(gsva_scores_matrix)
  cat("\nend of fast example for debug  \n*-*-*\n")
}

### functions
df_to_list <- function(df) {
   # Initialize an empty list to store results
   result_list <- list()
   
   # Loop through each row of the dataframe
   for (i in seq_len(nrow(df))) {
     # Get the row name
     row_name <- rownames(df)[i]
     
     # Extract the row as a vector
     row_values <- as.character(unname(df[i, ]))
     
     # Remove empty strings and NA values
     cleaned_values <- row_values[!(is.na(row_values) | row_values == "")]
     
     # Assign to list with rowname as key
     result_list[[row_name]] <- cleaned_values
   }
   
   return(result_list)
}


# start

args <- commandArgs(trailingOnly = TRUE)     # Note : args[0] # not usable

if (is.na(args[1])){run_minimal_test()}

###############  

# variables

data_dir_abspath = args[1]     # DATADIR (absolute/path/to/data)
SAMPLE = args[2]
PW_TYPE= args[3]
min_sample_size = 2



if (PW_TYPE == "KEGG"){
  pathways_path <- file.path(data_dir_abspath, 'metabolites-mappings', 
                           'KEGG_hsa_pathways_compounds_R110.gmt')
}

if (PW_TYPE == "WIKI"){
  pathways_path <- file.path(data_dir_abspath, 'metabolites-mappings', 
                             'WIKIPATHWAYS_human_multi_id_redundant.mmt')
}

matrix_file <- paste0(SAMPLE,"-cmps-mx-",PW_TYPE, ".tsv")
matrix_path <- file.path(data_dir_abspath, SAMPLE, matrix_file)
out_file = paste0(unlist(strsplit(x=matrix_file, split="\\."))[1], "-",
                  PW_TYPE,"_gsva.tsv")  


# open files
pre_gmt <- read.delim(pathways_path,  sep="\t",  header=FALSE,
                      row.names=1, 
                stringsAsFactors=FALSE)

pre_gmt$V1 <- ""   #   not required (pathway names), rownames provide the pathway identifiers

pre_gmt <- df_to_list(pre_gmt)


## annotations
annots_path <- file.path(data_dir_abspath, SAMPLE, paste0(SAMPLE,"-", PW_TYPE, "-CleanAnnotations.tsv"))

annots <- read.table(annots_path, row.names=NULL, header=TRUE, 
                     check.names = FALSE,
                     sep="\t", comment.char="", 
                     stringsAsFactors=FALSE)

column_in_annot <- ""
if (PW_TYPE == "WIKI"){
  column_in_annot = "moleculeIds"
}

if (PW_TYPE == "KEGG"){
  column_in_annot = "KEGG"
}


# ---------- pathways scores computation
start_time <- Sys.time()

# keep only the pathways that fill the condition: 
# at least 2 metabolite ids AND  at least 2 distinct m/z
compliant_pws = c()
for (pw in names(pre_gmt)){
  tmp_ann <- annots[annots[[column_in_annot]] %in% pre_gmt[[pw]] ,]
  if (length(unique(tmp_ann$mz)) >= 2) {
    compliant_pws <- c(compliant_pws,  pw)
  }
}
gmt <- pre_gmt[compliant_pws]


# loading quantification matrix
df <- read.table(matrix_path, row.names=1, header=TRUE, check.names = FALSE,
                sep="\t", comment.char="", 
                stringsAsFactors=FALSE)

X <- as.matrix(df)

colnames(X) = as.character(colnames(df))
cat("loaded quantifications matrix of shape: ", dim(X), "\n")

# compute  
X <- t(X)
cat("Log transform matrix\n")
X <- log10(X + 1)  # +1 to avoid log zero inf

cat("Computing GSVA: sample ", SAMPLE,"...\n")

gsvaPar <- gsvaParam(X, gmt,
                     minSize=as.integer(min_sample_size))  
rm(X)
gsva_scores_matrix <- gsva(gsvaPar, verbose=FALSE)

pathways_ids <- rownames(gsva_scores_matrix) 

exitr <- t(gsva_scores_matrix)

exitr <- cbind(
  data.frame("pathway_id"=rownames(exitr)),
  data.frame(exitr, row.names=NULL)
)
saveRDS(exitr, 
        file.path(data_dir_abspath, SAMPLE, paste0(out_file,".RDS")))
end_time <- Sys.time() # time measured up to saving the rds file

if (dim(exitr)[1] <= 90000){
  write.table(exitr, 
              file=file.path(data_dir_abspath, SAMPLE,  out_file),
              sep="\t", row.names=FALSE, col.names=TRUE)
}

cat("Sucess. Output tab delimited file saved to: ", 
    file.path(data_dir_abspath, out_file) , "\n")
duration_ms <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(SAMPLE, PW_TYPE, "Execution time in seconds:", duration_ms, "\n")



