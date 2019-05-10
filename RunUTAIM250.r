RunAdmixture <- function(sampleVCF,sampleNames
                         ,ref.genome = c("GRCh38","GRCh37"),temp, PopCode_data = "1000G_1652_PopCode.txt",Str_Path){
  setwd("~/admixture/Admixture_all_function/")
  ##Requirement tools
  #1.VCFtools version 0.1.15
  #2.R version 3.5.1
  #3.bcftools version: 1.7
  #4.STRUCTURE version
  library(vcfR)
  library(doParallel)
  library(ggplot2)
  library(missMDA)
  library(ade4)
  library(SNPlocs.Hsapiens.dbSNP.20120608)
  
  ##Modify VCF files
  Admixture_modifyVCF(inputFiles = sampleVCF,outputName = "UT-AIM250",threads = 1,temp = temp,positions = TRUE ,ref.genome = "grch37" ,chr.label = FALSE)
  
  ##Convert to Structure inport format
  Admixture_VCF2Structure(vcf1 = "UT-250AIMs_1000genome.vcf",vcf2 = "UT-AIM250_complete.recode.vcf",
                          PopCode_data = PopCode_data,GT_mergeTable = TRUE,
                          outputName = "UT-AIM250")
  
  ##Run Structure
  Admixture_InputStr(Input_Str = "UT-AIM250_str.txt" ,Str_Path = Str_Path ,Output_Str ="UT-AIM250_str_output.txt")
  
  ##Structure output modify
  ##distinguish Pop and group
  ##Population Summary Table
  Admixture_StrOut(Output_Str ="UT-AIM250_str_output.txt_f" ,PopGroup = c("AFR","EUR","EAS","Other"),
                   K_value = 3,fileName = "UT-AIM250_Summary.txt" )
  
  ##Individuals proportion Table
  Str_frq <- Admixture_StrOut_Frequency(Output_Str = "UT-AIM250_str_output.txt_f" ,PopGroup = c("AFR","EUR","EAS","Other"),K_value = 3,PopCluster = "UT-AIM250_Summary.txt",fileName = "UT-AIM250_Frequency.txt")
  
  ##Convert to PCA format
  df <- Admixture_gt2PCAformat(InputFile = "UT-AIM250_genotype.txt",PopCode_data = PopCode_data ,genotype_start_column = 6,
                               genotype_seperator = "|",outputName = "UT-AIM250")
  
  #df <- read.table("UT-AIM250_PCAFormat.txt", sep = "\t",stringsAsFactors = FALSE, header = TRUE)
  df_1000G <- df[which(!df$Pop == "Other"),]
  
  
  
  sampleList <- read.table(sampleNames, sep = "\t", stringsAsFactors = FALSE)
  cat("if you have only one sample, you will have Warning message!")
  if (!ncol(sampleList) == 1) { stop(" Wrong format : sampleNames")
  }
  
  
  Str_frq$Label[which(Str_frq$Pop == 10)] <- sampleList$V1
  Str_frq_1000G <- Str_frq[which(!Str_frq$Pop == 10),]
  if(nrow(sampleList) == 1) {
    
    
    #creat sub folder
    dirName <- paste(sampleList$V1,"AdmixtureAnalysis", sep = "_")
    dir.create(dirName, showWarnings = FALSE)
    
    outputName <- paste(dirName,sampleList$V1, sep = "/")
    
    ##PCA plot
    df_pca <- rbind(df[which(df$Individuals == sampleList$V1),], df_1000G)
    figFileName <- paste(dirName,sampleList$V1,sep = "/")
    Admixture_PCAplot_MissingValue(PCAformat = df_pca,PopSelect = "all" ,figFileName = figFileName)
    
    ##Proportion plot
    #Str_frq <- read.table("UT-AIM250_Frequency.txt", sep = "\t", header = TRUE,stringsAsFactors = FALSE)
    Admixture_TrianglePlot(Str_frq = Str_frq,figFileName = figFileName,PopCode_data = PopCode_data )
    
    
  }else{
    #creat sub folder
    dirName <- "UT-AIM250_AdmixtureAnalysis"
    dir.create(dirName, showWarnings = FALSE)
    
    
    ##PCA plot
    df_pca <- df
    figFileName <- paste(dirName,"UT-AIM250",sep = "/")
    Admixture_PCAplot_MissingValue(PCAformat = df_pca,PopSelect ="all" ,figFileName = figFileName)
    
    ##Proportion plot
    #Str_frq <- read.table("UT-AIM250_Frequency.txt", sep = "\t", header = TRUE,stringsAsFactors = FALSE)
    Admixture_TrianglePlot(Str_frq = Str_frq ,figFileName = figFileName,PopCode_data = PopCode_data )
    

    
  }
} 
