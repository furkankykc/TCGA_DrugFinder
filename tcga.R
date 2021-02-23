#!/Library/Frameworks/R.framework/Versions/4.0/Resources/Rscript
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")
# source("http://bioconductor.org/biocLite.R")
#useDevel()
library('TCGAbiolinks')
# test if there is at least one argument: if not, return an error

path = "/Users/furkankykc/PycharmProjects/TCGA_DrugFinder"
args <- commandArgs(trailingOnly = TRUE)
#project = "TCGA-READ"
project = args[1]
print(project)
query <- GDCquery(project = project,
data.category = "Clinical",file.type = "xml")
print("Downloading")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info ="patient")
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
out <- clinical.drug[,c("bcr_patient_barcode","drug_name","therapy_types","project")]
out_2 <-clinical[,c("bcr_patient_barcode","days_to_death")]
out_2$tobacco_smoking_history <- tryCatch(
  {
   clinical[,c("tobacco_smoking_history")]  },
  error = function(e){
    NA
  }
)
out <- merge(out,out_2)
out$tobacco_smoking_history[out$tobacco_smoking_history!=1]<- "Smoking"
out$tobacco_smoking_history[out$tobacco_smoking_history==1]<- "Non-Smoking"
write.csv(out,paste0(path,project,".csv"), row.names = FALSE)
