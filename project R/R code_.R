if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scater")
library(scater)
if(!requireNamespace("BiocManager",quietly = TRUE))
{
  install.packages("BiocManager")
BiocManager::install("GEOquery")
}


library(GEOquery)
gset= getGEO("GSE153428",destdir="."  )
length(gset)
gset<-gset[[1]]

ex=exprs(gset)
print(ex)

View(gset)

#filter first 20 row
group1_treat<-(exprs(gset[,4:12]))
group2_untreat<-(exprs(gset[,1:3]))

pValue=c()
tstatistic=c()
gene.id=c()
counter=0 #number of gene
rowNumber=1
while(counter!=20)
{
  test.gene=t.test(x=group1_treat [rowNumber,],y=group2_untreat[rowNumber,],var.equal = TRUE)
  if(!is.nan(test.gene$p.value)&test.gene$p.value<=0.05)
  {
    pValue=c(pValue,test.gene$p.value)
    tstatistic=c(tstatistic,test.gene$statistic)
    gene.id=c(gene.id,rowNumber)
    counter=counter+1
  }
  rowNumber=rowNumber+1
}
print(gene.id)


install.packages("PairedData")
install.packages("ggpubr")
install.packages("tidyverse")
library(ggpubr)  #view visulisation
library(tidyverse)  #filter groupBy
library(PairedData)


group1_treat_updated<-exprs(gset[gene.id,4:12])
group2_untreat_updated<-exprs(gset[gene.id,1:3])


# Create a data frame
my_data <- data.frame( 
  group = rep(c("group1_treat", "group2_untreat "), each=120),
  geneExpression = c(group1_treat_updated,  group2_untreat_updated)
)
print(my_data)


group_by(my_data, group) %>% 
  summarise(
    count = n(),
    mean = mean(geneExpression, na.rm = TRUE),
    sd = sd(geneExpression, na.rm = TRUE),
  )


#visulusation
boxplot(exprs(gset[gene.id,]),outline=FALSE,ylab="gene expression",xlab = "sample",main="treat and untreat",col ="#00AFBB")
hist(exprs(gset[gene.id,]),ylab="gene expression",xlab = "sample",col = "green",border = "red", xlim = c(0,20), ylim = c(0,15),
     breaks = 11)

# Subset data before treatment
subset_group1_treat<- subset(my_data,  group == "group1_treat_updated", geneExpression,drop = TRUE)
# subset data after treatment
subset_group2_untreat<-subset(my_data,  group == "group2_untreat_updated", geneExpression,drop = TRUE)
# Plot paired data
library(PairedData)
pd <- paired(group1_treat_updated, group2_untreat_updated)
plot(pd, type = "profile") + theme_bw()

# Compute t-test
res=t.test(geneExpression~group,data=my_data,alt="two.sided",paired = TRUE,var.equal=TRUE)
print(res)

