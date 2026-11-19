library(stringr)
library(pheatmap)

options(stringAsFactor = FALSE)

setwd("C:/Users/Lorenzo/Desktop/Università/Magistrale/Bioing/Tesina")

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
prc = 0.1          # percentile IQR (standard)
soglia_fc = 4      # soglia fold change
soglia_pval = 0.05 # soglia p value (standard)
####################

#importo i dati

tmp = read.table(filename_in, header = T, sep = "\t", check.names = F, quote = "",
                 row.names = 1, nrows = 3)
classes = sapply(tmp, class)
tmp = read.table(filename_in, header = T, sep = "\t", check.names = F, quote = "",
                 row.names = 1, colClasses = classes)

geni = rownames(tmp)
pz = colnames(tmp) #barcodes dei pazienti 

#cerchiamo pazienti sani e malati con grep()

pzN = grep("TCGA-\\w+-\\w+-1\\d", pz, value = T) #pazienti sani identificati con 1
pzC = grep("TCGA-\\w+-\\w+-0\\d", pz, value = T) #pazienti malati identificati con 0

#rimuovo i duplicati
pzN = pzN[!duplicated(str_extract(pzN,"TCGA-\\w+-\\w+-"))]
pzC = pzC[!duplicated(str_extract(pzC,"TCGA-\\w+-\\w+-"))]

nameN = str_extract(pzN, "TCGA-\\w+-\\w+-")
nameC = str_extract(pzC, "TCGA-\\w+-\\w+-")

common_pz = intersect(nameN,nameC)

pzN_com = unlist(lapply(common_pz, function(x){grep(x, pzN, value = TRUE)}))
pzC_com = unlist(lapply(common_pz, function(x){grep(x, pzC, value = TRUE)}))

dataN = tmp[,pzN_com]
dataC = tmp[, pzC_com]

data = cbind(dataN, dataC)

pz_com = colnames(data)
rm(classes, tmp, pz, nameN, nameC, pzN, pzC, common_pz)

#STEP 2 - trasformazione logaritmica per avere una distribuzione simmetrica

dataN = log2(dataN + 1) #ridurre variabilità del dato
dataC = log2(dataC + 1)
data = log2(data + 1)

#STEP 3: pre-processing calcolo l'IQR
variation = apply(data, 1, IQR, type=5) #variabilità la voglio su tuttii dati, su tutta la matrice ma tra normale e cancro
#butto via quelli con IQR basso ma considerando tutta la matrice
#type = 5 serve affnchè venga come matlab (di default è 7)

#grafico percentile

hist(variation, 
     main = "IQR frequency distribution", 
     breaks = 100, 
     xlab = "IQR value",
     ylab = "Frequency", 
     col = "green")

abline(v= quantile(variation, prc), lty = 2, lwd = 4, col = "red")

#Filtro su IQR

ind = which(variation <= quantile(variation, prc))

#i prossimi comandi servono per buttare via roba sotto a linea tratteggiata
dataN = dataN[-ind,]
data = data[-ind,]
dataC = dataC[-ind,]
geni = geni[-ind]

rm(ind, variation)

#STEP 4 - Filtraggio

#filtro su logFC

logFC = rowMeans(dataC)-rowMeans(dataN)

#devo fare il p value ma metto un filtro, non vogliamo tutti. 
#scegliamo quindi una soglia
#ista tra 1000e 2000 -> scegliere la soglia, magari 4 è stato pure troppo cattivo

#calcolo del pvalue

N = ncol(data)/2

#eseguiamo il paired t test per ogni gene prendendo come due gruppi i malati e i sani
pval =  apply(data, 1, function(x){ 
  t.test(x[1:N],x[(N+1):(2*N)], paired = TRUE)$p.value #estraiamo solo il pvalue
  })

#calcolo del pvalue corretto

pval_adj = p.adjust(pval, method = "fdr")

# Il pval_adj viene calcolato prima del filtraggio sulla treshold del fold change
# per far rimanere le due variabili indipendenti, se calcolassi il pvalue dopo il
# filtraggio queto varierebbe in base al numero di geni

plot(logFC, -log10(pval_adj),
     main = "Volcano plot",
     ylim = c(0, 70),
     xlab = "log2 fold change",
     ylab = "-log10 p-value")

abline(h = -log10(soglia_pval), lty = 2, lwd = 4, col="red")
abline(v = c(-log2(soglia_fc), log2(soglia_fc)), lty = 2, lwd = 4, col="red")
################################################################################

#Filtro sul fold change

ind = which( abs(logFC) < log2(soglia_fc))

if(length(ind)>0){ 
  dataN = dataN[-ind,]
  dataC = dataC[-ind,]
  data = data[-ind,]
  geni = geni[-ind]
  logFC = logFC[-ind]
  pval = pval[-ind]
  pval_adj = pval_adj[-ind]
}

rm(ind)

# Filtro sul pvalue

ind = which( pval_adj > soglia_pval)

if(length(ind) > 0){
  dataN = dataN[-ind,]
  dataC = dataC[-ind,]
  data = data[-ind,]
  geni = geni[-ind]
  logFC = logFC[-ind]
  pval = pval[-ind]
  pval_adj = pval_adj[-ind]
}

rm(ind)

#grafico volcano plot

plot(logFC, -log10(pval_adj),
     main = "Volcano plot",
     ylim = c(0, 70),
     xlab = "log2 fold change",
     ylab = "-log10 p-value")

abline(h = -log10(soglia_pval), lty = 2, lwd = 4, col="red")
abline(v = c(-log2(soglia_fc), log2(soglia_fc)), lty = 2, lwd = 4, col="red")

#scrittura dei risultati

direction = ifelse(logFC > 0, "UP", "DOWN")
result = data.frame(geni,pval,pval_adj,logFC, direction)

write.table(result, file = filename_DEG, row.names = F, sep = "\t", quote = F)

write.table(data, file = filename_matrix_DEG, row.names = T,
            col.names = NA, sep = "\t", quote = F)


#box plot UPREGOLATO

gene_id = names(which.max(logFC)) #cerchiamo il gene più UP regolato

#creo un data frame con i valori del gene 
df = data.frame(normal = t(dataN[gene_id,]),
                cancer = t(dataC[gene_id,]),
                row.names = NULL)

colnames(df) = c("normal","cancer")


boxplot(df,
        main = paste0(gene_id, ", ", "p-value = ",
                      format(pval_adj[gene_id], digits=2)),
        notch = T, #il notch ci permette di confrontare la mediana di bue box più vicine
        ylab = "Gene expression value",
        xlab = "Condition",
        col = c("green","orange"),
        pars = list(boxwex = 0.3, staplewex=0.6)
        )

#box plot DOWNREGOLATO

gene_id = names(which.min(logFC)) #cerchiamo il gene più DOWN regolato

#creo un data frame con i valori del gene 
df = data.frame(normal = t(dataN[gene_id,]),
                cancer = t(dataC[gene_id,]),
                row.names = NULL)

colnames(df) = c("normal","cancer")


boxplot(df,
        main = paste0(gene_id, ", ", "p-value = ",
                      format(pval_adj[gene_id], digits=2)),
        notch = T, #il notch ci permette di confrontare la mediana di bue box più vicine
        ylab = "Gene expression value",
        xlab = "Condition",
        col = c("green","orange"),
        pars = list(boxwex = 0.3, staplewex=0.6)
)


#STEP 7: HEATMAP

test = grepl("TCGA-\\w+-\\w+-0\\d", colnames(data))

samples = ifelse(test, "normal", "cancer")

#creiamo un dataframe che associa ogni barcode ad un'etichetta normal o cancer
annotation_col = data.frame(SampleType = samples)
rownames(annotation_col) = colnames(data)

vect_color_col = c("green","orange")
names(vect_color_col) = c('normal','cancer')

annotation_row = data.frame(direction = result$direction)
rownames(annotation_row)  = result$geni

vect_color_row = c('magenta','cyan')
names(vect_color_row) = c('UP','DOWN')

annotation_colors = list(SampleType = vect_color_col,
                         direction = vect_color_row)

pheatmap(data, scale="row",
         border_color = NA,
         cluster_cols = T,
         cluster_rows = T,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         color = colorRampPalette(colors = c("blue","blue3","black","yellow3","yellow"))(100),
         show_rownames = F,
         show_colnames = F,
         cutree_cols=2,
         cutree=2,
         width = 10,height = 10,
         filename = filename_heatmap,
         #filename_heatmap<-paste0(dirTumor,"heatmap.pdf")
)














