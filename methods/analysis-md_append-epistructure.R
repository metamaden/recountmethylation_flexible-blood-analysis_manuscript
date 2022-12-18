library(HDF5Array)

# load data
save.dpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_final_results")
# get grset
gr.fname <- "gr-gseadj_h5se_hm450k-epic-merge_0-0-3" # "gr-noob_h5se_hm450k-epic-merge_0-0-3"
gr.fpath <- file.path(save.dpath, gr.fname)
gr <- loadHDF5SummarizedExperiment(gr.fpath)
colData(gr) <- colData(gr)[,c(1:56)]
# get glint results
glint.fname <- "glint_results.epistructure.pcs.txt"
epi <- read.table(file.path(save.dpath, glint.fname))
colnames(epi) <- c("sample", paste0("glint.epi.pc",seq(10)))
# append glint results to coldata
identical(as.character(epi[,1]), as.character(colnames(gr)))
for(ii in seq(10)){
  epi.var <- paste0("glint.epi.pc", ii)
  colData(gr)$epi <- epi[,epi.var]
  colnames(colData(gr))[ncol(colData(gr))] <- epi.var
}

# save updated object
quickResaveHDF5SummarizedExperiment(gr)

# save new metadata
md.fname <- "mdf-epic-hm450k-merge_normal-blood_0-0-3.rda"
md.fpath <- file.path(save.dpath, md.fname)
md <- as.data.frame(colData(gr))
save(md, file = md.fpath)