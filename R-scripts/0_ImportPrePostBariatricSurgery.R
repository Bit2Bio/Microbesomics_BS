#Pipeline Per Import Dati Microbesomics

#Import Librerie
library(Microsoft365R)
library(microbiomeMarker)
library(phyloseq)
library(xlsx)
library(reticulate)
use_condaenv("python_env")

py_run_string("
import pandas as pd
import numpy as np
import statsmodels.formula.api as sm
import pingouin as ping           
")

# Ottieni la data e l'orario correnti nel formato desiderato
data_orario <- format(Sys.time(), "%d-%m-%Y")

# Specifica il percorso della nuova cartella
percorso_cartella <- paste0("../data/", "Pre-processing", data_orario)

#Crea La Cartella
dir.create(percorso_cartella)

## Specifica il percorso dei risultati su OneDrive
onedriveres = "Progetti/Microbesomics/Microbesomics_MetS/results/Metagenomica"

# Crea la cartella
dir.create(paste0(percorso_cartella, "/","qza"))

#Import Dati da OneDrive
dy = get_business_onedrive()

#Path dei Risultati
path = "Progetti/Microbesomics/Microbesomics_MetS/Docs"

#qza_path
qza_path = paste0(percorso_cartella,"/","qza")

#Download QZA
if (!dir.exists(qza_path) || length(list.files(qza_path)) == 0) {
  dy$download_folder(paste0(path, "/", "qza"), dest = qza_path)
  print(paste0("QZA Files Downloaded in: ", qza_path))
} else {
  print(paste0("QZA Files already exist in: ", qza_path, ", skipping download."))
}

#Create Ps Object
ps<-import_qiime2(paste0(qza_path, "/", "asv_table.qza"),
                  taxa_qza = paste0(qza_path, "/", "taxonomy.qza"),
                  sam_tab = NULL,
                  tree_qza =  paste0(qza_path, "/", "fasttree_tree_rooted.qza") ,
                  refseq_qza = paste0(qza_path, "/", "rep_seq.qza"))

#Change Sample Names
sample_names(ps) = as.character(1:247)

#Import Metadata da OneDrive
SmpDt_path = paste0(percorso_cartella,"/","Copia di Microobesomics_clinical data IGA 260324.xlsx")
if (!file.exists(SmpDt_path)) {
  dy$download_file("Progetti/Microbesomics/Microbesomics_MetS/Docs/Copia di Microobesomics_clinical data IGA 260324_LA.xlsx",
                   dest = SmpDt_path)
  print(paste0("Sample Data downloaded at: ", SmpDt_path))
} else {
  print(paste0("Sample Data already exists at: ", SmpDt_path, ", skipping download."))
}


#Import Metadata da OneDrive
mtbDt_path = paste0(percorso_cartella,"/","tentativo riassunto microbesomicsLA.xlsx")
if (!file.exists(mtbDt_path)) {
  dy$download_file("Progetti/Microbesomics/Microbesomics_MetS/Docs/tentativo riassunto microbesomicsCorretto.xlsx",
                   dest = mtbDt_path)
  print(paste0("Metab Data downloaded at: ", mtbDt_path))
} else {
  print(paste0("Metab Data already exists at: ", mtbDt_path, ", skipping download."))
}

mtbDt = read.xlsx(mtbDt_path, check.names = F, row.names = 1, sheetIndex = 1)

#Import Excel
sampleDt = read.xlsx("../data/clinics/Copia di Microobesomics_clinical data IGA 260324.xlsx", check.names = F, sheetName = "group5vs6")
sampleDt =sampleDt[sampleDt$`GROUP IGA`%in% c("5", "6"),]
row.names(sampleDt) = sampleDt$`CODE IGA`
py_run_string("
sampleDt = r.sampleDt
mtbDt = r.mtbDt
mtbDt = mtbDt.dropna(subset=['IGA'])
mtbDt.index = mtbDt['IGA'].astype(int).astype(str)

smps = (set(sampleDt.index).intersection(mtbDt.index))
smps = [str(y) for y in np.sort([int(x) for x in (list(smps))])]
sampleDt = sampleDt.loc[smps]
mtbDt = mtbDt.loc[smps]      
")

sampleDt = py$sampleDt
mtbDt = py$mtbDt
row.names(mtbDt) = mtbDt$IGA
#Aggiungiamo SampleData al Phyloseq Object
ps = merge_phyloseq(ps , sample_data(sampleDt))

#Scarichiamo i dati Filtrati da OneDrive (MicrobiomeAnalyst)
filt_path = paste0(percorso_cartella,"/","filteredData.csv")
if (!file.exists(filt_path)) {
  dy$download_file("Progetti/Microbesomics/Microbesomics(n=177)/Docs/MicrobiomeAnalyst/filtering/Filtering_n177/data_filtered.csv",
                   dest = filt_path)
  print(paste0("Filtered data downloaded at: ", filt_path))
} else {
  print(paste0("Filtered data already exists at: ", filt_path, ", skipping download."))
}

#Import dati filtrati (Verranno utilizzati nella parte finale della pipeline per DESeq2)
filtOtuTab = read.csv(filt_path, row.names=1, check.names = F)

keep_taxa <- rownames(filtOtuTab)

# filtra l'oggetto phyloseq
ps_filt1vs2 <- prune_taxa(taxa_names(ps) %in% keep_taxa, ps)

#Rimuoviamo Taxa non identificati a livello di Phylum
ps_filt1vs2 = subset_taxa(ps_filt1vs2, Phylum != "NA")
mtbDt = mtbDt[,4:ncol(mtbDt)]
saveRDS(ps_filt1vs2, file = paste0("../data/metagenomics/ps_filt5vs6.rds"))
saveRDS(mtbDt, file = paste0("../data/metabolomics/mtbDt5vs6.rds"))
saveRDS(sampleDt, file = paste0("../data/clinics/sampleDt5vs6.rds"))



