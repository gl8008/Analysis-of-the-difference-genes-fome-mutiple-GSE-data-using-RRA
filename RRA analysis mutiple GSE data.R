###Aanlysis difference gene from mutiple GSE data using RRA (2021.12.28)

library(RobustRankAggreg)
library(AnnoProbe)
library(GEOmirror)
library(GEOquery) 
library(devtools)
library(tidyverse)
library(idmap2)
library(limma) 
library(dplyr)
library(tidyr)
## Download GSE118370 data
library(GEOmirror)
GSE = 'GSE118370'
eSet=geoChina(GSE)
eSet=eSet[[1]]
probes_expr <- exprs(eSet);dim(probes_expr)
probes_expr [is.na(probes_expr)]<- 0 
phenoDat <- pData(eSet) 
colnames(phenoDat)
head(phenoDat)[1,]
Tumor_expr = probes_expr[ , grep( ('tumor tissue'),   phenoDat$`tissue:ch1`)]
Normal_expr  = probes_expr[  , !(colnames(probes_expr) %in% colnames(Tumor_expr)) ]
exprSet=cbind(Normal_expr, Tumor_expr)

groupList = c(rep( 'Normal', ncol( Normal_expr ) ), 
              rep( 'Tumor', ncol( Tumor_expr ) ) )

groupList = factor(groupList)
groupList <- relevel(groupList, ref="Normal")
table(groupList)
gpl=eSet@annotation
gpl
checkGPL(gpl) 
printGPLInfo(gpl)

library(GEOquery)
GPL <-getGEO(gpl,destdir =".")

GPL_anno <- Table(GPL)
#head(GPL_anno)
colnames(GPL_anno)
probe2symbol_df <- GPL_anno %>% 
  
  dplyr::select(ID,c("Gene Symbol") ) %>%
  
  filter(c("Gene Symbol") != "")

colnames(probe2symbol_df) <-c("ID","Symbol")
probe2symbol_df=data.frame(probe2symbol_df)
head(probe2symbol_df)
genes_expr=exprSet
genes_expr[1:4,1:4]

genes_expr <- normalizeBetweenArrays(genes_expr)
ex <- genes_expr
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
genes_expr <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
genes_expr <- as.data.frame(genes_expr)
genes_expr <- genes_expr %>% 
  rownames_to_column(var="ID") %>% 
  inner_join(probe2symbol_df,by="ID") %>% 
  select(Symbol,everything()) %>% 
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% 
  filter(Symbol != "NA") %>% 
  arrange(desc(rowMean)) %>% 
  distinct(Symbol,.keep_all = T) %>% 
  select(-rowMean) %>% 
  filter(Symbol !="" ) %>% 
  column_to_rownames(var = "Symbol")
genes_expr <- genes_expr[,-1]
library("FactoMineR")
library("factoextra")
dat.pca <- PCA(t(genes_expr) , graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = groupList,
             addEllipses = TRUE,
             legend.title = "Groups"
)

library(limma)
design <- model.matrix(~0+factor(groupList))
colnames(design) <- levels(factor(groupList))
rownames(design) <- colnames(exprSet)
contrast.matrix <- makeContrasts(paste0(unique(groupList),collapse = "-"),levels = design)
contrast.matrix <- makeContrasts(Tumor - Normal,levels = design)
contrast.matrix

fit <- lmFit(genes_expr,design)

fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2,0.01) 
tempOutput <- topTable(fit2,adjust="fdr",number = Inf,sort.by = "B",coef = 1)
DEG <-  na.omit(tempOutput) 
dim(DEG)
DEG.sig <- DEG[which(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > 1),]
dim(DEG.sig)
GSE
write.csv(DEG,paste(GSE,"DEGen_results.0317.csv",sep = "_"),row.names = T)

## for volcano plot
df=DEG
df <- rt
attach(df)
df$v= -log10(P.Value)
df$g=ifelse(df$P.Value>0.05,'stable',
            ifelse( df$logFC >1,'up',
                    ifelse( df$logFC < -1,'down','stable') )
)
table(df$g)
df$name=rownames(df)
head(df)
library(ggpubr)
ggpubr::ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
                  label = "name", repel = T,
                  label.select =head(rownames(df)),
                  palette = c("#00AFBB", "#E7B800", "#FC4E07") )
detach(df)

x=DEG$logFC
names(x)=rownames(DEG)
cg=c(names(head(sort(x),100)),
     names(tail(sort(x),100)))
cg
library(pheatmap)
n=t(scale(t(genes_expr[cg,])))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
ac=data.frame(groupList=groupList)
rownames(ac)=colnames(n) 
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac)

## Analysis other eight GSE data (GSE19188,GSE27262,GSE32863,GSE101929,GSE74706,GSE75037,

## GSE18842,GSE31210)
####RRA analysis
padj=0.05
logFC=1

files=c("GSE19188_DEGen_results.0317.csv","GSE27262_DEGen_results.0317.csv",
        "GSE32863_DEGen_results.0317.csv","GSE101929_DEGen_results.0317.csv",
        "GSE74706_DEGen_results.0317.csv","GSE75037_DEGen_results.0317.csv",
        "GSE18842_DEGen_results.0317.csv","GSE31210_DEGen_results.0317.csv",
        "GSE118370_DEGen_results.0317.csv")

upList=list()
downList=list()
allFCList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.csv(paste("C:/LTM/LUAD/",inputFile,sep=""), header=T, na.strings=c("NA"))
  rt$X <- gsub(". ///.*", "", rt$X)
  rt=rt[order(rt$logFC),] 
  header=unlist(strsplit(inputFile,"_"))
  downList[[header[1]]]=as.vector(rt[,1]) 
  upList[[header[1]]]=rev(as.vector(rt[,1])) 
  fcCol=rt[,1:2]
  colnames(fcCol)=c("Gene",header[[1]])
  allFCList[[header[1]]]=fcCol
}


mergeLe=function(x,y){
  merge(x,y,by="Gene",all=T)}
newTab=Reduce(mergeLe,allFCList)
newTab[is.na(newTab)]=0
newTab <- newTab[!duplicated(newTab$Gene),]
rownames(newTab)=newTab[,1]
newTab=newTab[,2:ncol(newTab)]
colnames(newTab)



library(RobustRankAggreg)
upMatrix = rankMatrix(upList)
upAR = aggregateRanks(rmat=upMatrix)
colnames(upAR)=c("Name","Pvalue")
upAdj=p.adjust(upAR$Pvalue,method="bonferroni")
upXls=cbind(upAR,adjPvalue=upAdj)
upFC=newTab[as.vector(upXls[,1]),]
upXls=cbind(upXls,logFC=rowMeans(upFC))
upSig=upXls[(upXls$adjPvalue<padj & upXls$logFC>logFC),]
dim(upSig)

downMatrix = rankMatrix(downList)
downAR = aggregateRanks(rmat=downMatrix)
colnames(downAR)=c("Name","Pvalue")
downAdj=p.adjust(downAR$Pvalue,method="bonferroni")
downXls=cbind(downAR,adjPvalue=downAdj)
downFC=newTab[as.vector(downXls[,1]),]
downXls=cbind(downXls,logFC=rowMeans(downFC))

downSig=downXls[(downXls$adjPvalue<padj & downXls$logFC< -logFC),]
dim(downSig)

allSig = rbind(upSig,downSig)
colnames(allSig)
dim(allSig)
write.table(allSig,file = 'C:/LUAD/9_GEO_allSign.txt',sep = '\t',quote = F)
####logFC.tiff
hminput=newTab[c(as.vector(upSig[1:50,1]),as.vector(downSig[1:50,1])),] #各取上下调50个基因画图

n =1

for(i in hminput[,n])
{
  if(n<8)
  {hminput[,n] <- ifelse(hminput[,n] >= 8, 8, hminput[,n])
  hminput[,n] <- ifelse(hminput[,n] <= -8, -8, hminput[,n])
n=n+1}
  else{
    hminput[,9] <- ifelse(hminput[,n] <= -8, -8, hminput[,n]) 
  }
}

library(pheatmap)

pheatmap(hminput,display_numbers = T,
         fontsize_row=8,
         fontsize_col=10,
         color = colorRampPalette(c("green", "white", "red"))(50),
         cluster_cols = FALSE,cluster_rows = FALSE )


