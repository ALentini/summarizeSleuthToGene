#Extracts estimated RNA-seq count data from a sleuth object and summarizes it to gene-level using the tximport algorithm
#Uses the data.table package for speed but can be changed to base::data.frame eg. using reshape2::dcast
summarizeSleuthToGene <- function(sleuth_object, tx2gene, which_df = c("obs_norm_filt", "obs_norm", "obs_raw"), ... ){
  if( class(sleuth_object) != 'sleuth' ) stop( '"sleuth_object" must be a sleuth object' )
  if( missing(sleuth_object) ) stop( 'argument "tx2gene" is missing with no default' )
  which_df <- match.arg(which_df, c("obs_norm_filt", "obs_norm", "obs_raw") )
  require(sleuth)
  require(tximport)
  require(data.table)
  dat <- data.table(sleuth_object[[which_df]])
  cnt_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "est_counts"),row.names = 1))
  tpm_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "tpm"),row.names = 1))
  len_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "eff_len"),row.names = 1))
  txi <- list(abundance=tpm_mat,counts=cnt_mat,length=len_mat)
  summarizeToGene(txi,tx2gene=tx2gene, ... )
}
