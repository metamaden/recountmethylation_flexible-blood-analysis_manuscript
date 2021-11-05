#!/usr/bin/env R

# Author: Sean Maden
#
# Validation for study Sadikovik et al 2021, a DNAm diagnostic 
# panel for rare Mendelian disorders.
#
# Blood sample studies
#
# 1. "GSE55491"
# title: "Genomewide methylation analysis in Silver Russell syndrome patients"
# whole blood from silver russell
# PBMC from normal controls
# descriptive sample titles
#
# 2. "GSE74432" 
# title: "NSD1 Mutations Generate a Genome-Wide DNA Methylation Signature"
# descriptive sample titles
# whole blood from controls, Sotos, Weaver, NSD1 var patients
# fibroblasts from control and Sotos patients
#
# 3. "GSE97362" 
# title: "CHARGE and Kabuki syndromes: Gene-specific DNA methylation signatures"
# whole blood from Kabuki, CHARGE, and control patients
#
# 4. "GSE108815" 
# title: "Methylomic Analysis reveals Epigenetic Pattern of SH2B1 in Chinese Monozygotic Twins Discordant for Autism Spectrum Disorder"
# whole blood from twin pairs with autism and non-autistic control
# desciptive sample titles
#
# Brain sample studies:
#
# 1. "GSE80017"  
# title: "Meta-Analysis of Autism Cortical Transcriptome Reveals Discordant Region-Specific Changes Linked to Radial Glia and Inhibitory Interneurons"
# prefrontal cortex samples from normal controls, autistic patients
#
# 2. "GSE131706" 
# title: "Aberrant intragenic DNA methylation in a neurogenic region in autism spectrum disorders"
# desciptive sample titles
# postmortem brain from autistic and control patients


# GSE55491, GSE74432, GSE97362, GSE108815
# GSE80017, GSE131706

gsev <- c("GSE55491", "GSE74432", "GSE97362", "GSE108815", "GSE80017", "GSE131706")

gsev.new <- c('GSE55491','GSE60274','GSE74432','GSE80017','GSE95488','GSE108815',
              'GSE116992','GSE119778','GSE128068','GSE131706','GSE166373',
              'GSE100825','GSE120558','GSE133774','GSE157252')
gsev.check <- gsev.new[!gsev.new %in% gsev]; gsev.check

# GSE95488
# title: "Wilms Tumor in Beckwith-Wiedemann Syndrome and Loss of Methylation at Imprinting Centre 2"
# snp array and dnam array
# blood sample cases Beckwith-Wiedemann Syndrome, and controls
#
# GSE116992
# title: "Genome-wide DNA methylation analysis of the whole blood of individuals with Coffin-Siris and Nicolaides-Baraitser syndromes"
# whole blood cases (Coffin-Siris syndrome and Nicolaides-Baraitser syndrome)
# 
#
# GSE119778
# title: "Genome-wide DNA methylation analysis in Williams syndrome (WS patients vs Controls)"
# whole blood/peripheral cases (Williams syndrome) and controls
# 
#
# GSE133774
# title: "Clinical spectrum of multi-locus imprinting disturbances associated with maternal-effect variants range from overt Beckwith-Wiedemann syndrome to apparently healthy phenotype"
# blood from cases (Beckwith-Wiedemann syndrome) and controls

#----------
# load data
#----------
md.fpath <- "si1_all-md-2platforms.rda"
md <- get(load(md.fpath))

#---------
# query md
#---------
dxv <- c("(T|t)halassemia", "(R|r)etardation", "(X|x)-linked",
         "(A|a)utism", "(I|i)ntellectual( |_)(D|d)isability", 
         "(A|a)taxia", "(D|d)eaf", "(N|n)arcolepsy", "CHARGE", 
         "charge", "(K|k)abuki", "(E|e)pilep", "(S|s)otos", 
         "(F|f)ragile", "(P|p)rader( |_|-)(W|w)illi", "(S|s)ilver", 
         "(R|r)ussel", "(T|t)emple", "(K|k)agami", "(W|w)eaver", 
         "(A|a)ngelman", "(O|o)gatta", "(B|b)eckwith( |_|-)(W|w)iedemann",
         "(S|s)yndrome")

dxv.str <- paste0(paste0(".*", dxv, ".*"), collapse = "|")
filt.cond <- grepl(dxv.str, md$gsm_title)
filt.cond <- filt.cond | grepl(dxv.str, md$disease)
table(filt.cond)
# filt.cond
# FALSE  TRUE 
# 68658   100

# get studies and samples
gsev <- unique(md[filt.cond,]$gse) # string matches
gsev <- c(gsev, "GSE97362") # charge/kabuki study
# exclude cancer studies and other FP matches
gsev.exclude1 <- c("GSE60274", "GSE128068", "GSE166373", "GSE100825", 
                  "GSE120558", "GSE157252")
# exclude geo studies used to train the episign model (Aref-Eshghi et al 2020)
gsev.exclude2 <- c("GSE116992", "GSE66552", "GSE74432", "GSE97362", 
                   "GSE116300", "GSE95040", "GSE104451", "GSE125367", 
                   "GSE55491", "GSE108423", "GSE116300", "GSE89353", 
                   "GSE52588", "GSE42861", "GSE85210", "GSE87571", 
                   "GSE87648", "GSE99863", "GSE35069")
gsev <- gsev[!gsev %in% c(gsev.exclude, gsev.exclude2)]
mdf <- md[md$gse %in% gsev,]
dim(mdf) # [1] 158  18
table(mdf$platform)
# epic hm450k 
# 10    148
length(unique(mdf$gse)) # [1] 6

#---------------------------------------------
# test episign model in Aref-Eshghi et al 2020
#---------------------------------------------
# get study st's with panel info
st2.fname <- "st2_episign_probes.csv"; st3.fname <- "st3_episign_dnam.csv"
st2 <- read.csv(st2.fname); st3 <- read.csv(st3.fname)
identical(st2[,1], st3[,1]) # TRUE
les <- list("cgpanel" = st2[,1], "dxscore" = st2, "dnamscore" = st3)

# episign function
episign <- function(bval, les){
  message("Filtering bval matrix..."); cgpanel <- les[["cgpanel"]]
  bf <- bval[rownames(bval) %in% cgpanel,]
  bf <- bf[order(match(rownames(bf), cgpanel)),]
  message("Getting score by dx..."); lres <- list()
  mdx <- les[["dxscore"]]; mdx <- mdx[,c(4:ncol(mdx))]
  
  which.cg <- 1
  cgi <- cgpanel[which.cg]; mdxi <- mdx[which.cg,]
  bf[cgi,1] <= mdxi
  
  bfi <- bf[,1]
  mbool <- do.call(rbind, lapply(seq(length(cgpanel)), function(cgi){
    bfi[cgi] >= mdx[cgi,]
  }))
  
  dxscore <- les[["dxscore"]]
  dxv <- colnames(dxscore)[4:ncol(dxscore)]
  
  mbool.dxv <- unlist(lapply(seq(4,ncol(dxscore),1), function(ii){
    mdxscore <- mbool[which(dxscore[,ii]),1]
    length(which(mdxscore))/length(mdxscore)}))
  names(mbool.dxv) <- dxv
  
  
  message("Getting dx predictions...")
  
  return(lres)
}

# get dnam, epic
rg.epic.fname <- "remethdb_h5se_rg_epic-hm850k_merged_1621537799-1589820348_0-0-3"
rg <- HDF5Array::loadHDF5SummarizedExperiment(rg.epic.fname)
# filter on studies, panel probes
gsev <- unique(md[filt.cond,]$gse)
gsmv <- md[md$gse %in% gsev,]$gsm
labelv <- rownames(md[md$gsm %in% gsmv,])
anno <- getAnnotation(rg); annof <- anno[cgpanel,]
addrv <- c(annof$AddressA, annof$AddressB)
rgf <- rg[rownames(rg) %in% addrv, colnames(rg) %in% labelv]
dim(rgf) # [1] 4630  185
bval <- getBeta(rgf)




