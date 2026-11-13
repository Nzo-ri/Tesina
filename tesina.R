library(stringr)
library(pheatmap)

options(stringAsFactor = F)

setwd("D:/Lorenzo/Desktop/Università/Magistrale/1° anno/Bioingegneria per la genomica/Tesina")

#creiamo la directory dei risultati

dirRes = "Result/"
if(!dir.exists(dirRes)){
  dir.create(dirRes)
}else{
  print(paste(dirRes,"already exists"))
}

#creiamo la matrice dei geni del tumore

tumor = "brca"
dirTumor = paste0(dirRes,tumor,"/")
if(!dir.exists(dirTumor)){
  dir.create(dirTumor)
}else{
  print(paste(dirTumor,"already exists"))
}

filename_in = paste0(dirTumor,"matrice_",tumor,"_RNASeq.txt")

#creiamo file di output

filename_DEG = paste0(dirTumor,"DEG.txt")
filename_matrix_DEG = paste0(dirTumor,"matrix_DEG.txt")
filename_heatmap = paste0(dirTumor,"heatmap.pdf")

#configuro i valori di soglia per l'analisi dei geni
####################
prc = 0.1          #
soglia_fc = 4      #
soglia_pval = 0.05 #
####################

#importo i dati

tmp = read.table(filename_in, header = T, sep = "\t", check.names = F, quote = "", row.names = 1, nrows = 3)
classes = sapply(tmp, class)
tmp = read.table(filename_in, header = T, sep = "\t", check.names = F, quote = "", row.names = 1, colClasses = classes)

geni = rownames(tmp)
pz = colnames(tmp) #barcodes dei pazienti 

#cerchiamo pazienti sani e malati con grep()

pzN = grep("TCGA-\\w+-\\w+1\\d", pz, value = T) #pazienti sani identificati con 1
pzC = grep("TCGA-\\w+-\\w+0\\d", pz, value = T) #pazienti malati identificati con 0

#eliminiamo gli eventuali duplicati
pzN = pzN[!duplicated(str_extract(pzN,"TCGA-\\w+-\\w+-"))]
pzC_duplicati = pzC[!duplicated(str_extract(pzN,"TCGA-\\w+-\\w+-"))]
pzc_univoci = pzC[!duplicated(str_extract(pzC,"TCGA-\\w+-\\w+-"))]

nameN = str_extract(pzN, "TCGA-\\w+-\\w+-")
nameC = str_extract(pzC, "TCGA-\\w+-\\w+-")

common_pz = intersect(nameN,nameC)

pzN_com = unlist(lapply(common_pz, function(x){grep(x, pzN, value = TRUE)}))
pzC_com = unlist(lapply(common_pz, function(x){grep(x, pzN, value = TRUE)}))

dataN = tmp[,pzN_com]
dataC = tmp[, pzC_com]

data = cbind(dataN, dataC)













