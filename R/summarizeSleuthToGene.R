#Extracts estimated RNA-seq count data from a sleuth object and summarizes it to gene-level using the tximport algorithm
#Uses the data.table package for speed but can be changed to base::data.frame using reshape2::dcast
summarizeSleuthToGene <- function(sleuth_object,tx2gene, ... ){
  require(sleuth)
  require(tximport)
  require(data.table)
  dat <- data.table(sleuth_object$obs_norm)
  cnt_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "est_counts"),row.names = 1))
  tpm_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "tpm"),row.names = 1))
  len_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "eff_len"),row.names = 1))
  txi <- list(abundance=tpm_mat,counts=cnt_mat,length=len_mat)
  summarizeToGene(txi,tx2gene=tx2gene, ... )
}
