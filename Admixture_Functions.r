# Fuctions of admixture program

#Compress and modify vcf files 
Admixture_modifyVCF <- function(inputFiles, outputName, threads = 1, temp, positions = FALSE, ref.genome = c("GRCh38","GRCh37"), chr.label = FALSE){
  
  if(positions == FALSE){
    cat("Compress VCF files","\n",paste("Threads:", threads),"\n")
    CMD1 <- paste0("bgzip -c --threads ",threads," ",inputFiles," > ",paste0(outputName,".vcf.gz"))
    print(CMD1)
    system(CMD1)
    
    CMD2 <- paste0("tabix -p vcf ",paste0(outputName,".vcf.gz"))
    print(CMD2)
    system(CMD2)
    
    cat("Modify VCF files \n")
    CMD3 <- paste0("bcftools annotate -x INFO,^FORMAT/GT ",paste0(outputName,".vcf.gz")," --output-type z --output ",paste0(outputName,"_SNP_GT_noINFO.vcf.gz"))
    print(CMD3)
    system(CMD3)
    
    CMD4 <- paste0("vcftools --gzvcf ",paste0(outputName,"_SNP_GT_noINFO.vcf.gz")," --min-alleles 2 --max-alleles 2 --remove-indels --recode-INFO-all --max-missing 0.9 --temp ",temp," --recode --out ",paste0(outputName,"_complete"))
    
    print(CMD4)
    system(CMD4)
    
  }else{
    cat("Compress VCF files","\n",paste("Threads:", threads),"\n")
    CMD1 <- paste0("bgzip -c --threads ",threads," ",inputFiles," > ",paste0(outputName,".vcf.gz"))
    print(CMD1)
    system(CMD1)
    
    CMD2 <- paste0("tabix -p vcf ",paste0(outputName,".vcf.gz"))
    print(CMD2)
    system(CMD2)
    
    cat("Modify VCF files \n")
    CMD3 <- paste0("bcftools annotate -x INFO,^FORMAT/GT ",paste0(outputName,".vcf.gz")," --output-type z --output ",paste0(outputName,"_SNP_GT_noINFO.vcf.gz"))
    print(CMD3)
    system(CMD3)
    
    if(tolower(ref.genome) == "grch38"){
      if(chr.label == FALSE){
        AIM250_GRCH38 <- "AIM_250_grch38.txt"
      }else{
        AIM250_GRCH38 <- "AIM_250_chr_grch38.txt" 
      }
      
      CMD4 <- paste0("vcftools --gzvcf ",paste0(outputName,"_SNP_GT_noINFO.vcf.gz")," --positions ", AIM250_GRCH38," --min-alleles 2 --max-alleles 2 --remove-indels --recode-INFO-all --max-missing 0.9 --temp ",temp," --recode --out ",paste0(outputName,"_complete"))
      
    }else{
      if(chr.label == FALSE){
        AIM250_GRCH37 <- "AIM_250_grch37.txt"
      }else{
        AIM250_GRCH37 <- "AIM_250_chr_grch37.txt" 
      }
      
      CMD4 <- paste0("vcftools --gzvcf ",paste0(outputName,"_SNP_GT_noINFO.vcf.gz")," --positions ", AIM250_GRCH37," --min-alleles 2 --max-alleles 2 --remove-indels --recode-INFO-all --max-missing 0.9 --temp ",temp," --recode --out ",paste0(outputName,"_complete"))
    }
    
    print(CMD4)
    system(CMD4)
    
  }
}




###VCF files convert ot Struture
Admixture_VCF2Structure<-function(vcf1 = "UT-250AIMs_1000genome.vcf",vcf2=NULL,PopCode_data = "1000G_1652_PopCode",GT_mergeTable=FALSE,outputName){
  library(vcfR)
  library(SNPlocs.Hsapiens.dbSNP.20120608)
  
  if(!require(vcfR)){stop("Missing vcfR package")}
  
  #Check PopCode_data
  PopCode<-read.table(file = PopCode_data,sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  if(is.null(PopCode_data)){stop("Missing PopCode_data")}
  if(!sum(colnames(PopCode) == "sample") == 1){stop("PopCode_data Format: no sample names")}
  if(!sum(colnames(PopCode) == "pop") == 1){stop("PopCode_data Format:no populaiton")}
  
  #Check ouputName 
  if(is.null(outputName)){ouputName="vcf1"}
  
  
  # Check 250AIMs_1000genome.recode.vcf
  if(is.null(vcf1)) stop("No 250AIMs_1000genome.recode.vcf")
  
  if(is.null(vcf2)){
    vcf_1<-read.vcfR(vcf1,convertNA = FALSE)
    #genotype
    GT_1<-as.data.frame(vcf_1@gt,stringsAsFactors = FALSE)
    GT_1<-GT_1[,-1]
    
    #check genotype seperator
    if(!sum(grepl("[[:digit:]]/[[:digit:]]",GT_1))==0){
      GT_revise_1<-as.data.frame(matrix(data = NA,nrow = nrow(GT_1),ncol = ncol(GT_1)))
      colnames(GT_revise_1)<-colnames(GT_1)
      for (i in 1:nrow(GT_1)) {for (j in 1:ncol(GT_1)) {
        if(GT_1[i,j]=="."){GT_revise_1[i,j]<-"./."}else{GT_revise_1[i,j]<-as.character(GT_1[i,j])}
        GT_revise_1[i,j]<-gsub("/","|",GT_revise_1[i,j])
      }
      }
    }else{ GT_revise_1=GT_1}
    
    cat("Genotype Finished!","\n")
    
    #REF&ALT
    REF_1<-as.data.frame(vcf_1@fix,stringsAsFactors = FALSE)
    REF_1$CHROM <-as.character(REF_1$CHROM)
    if(!sum(grepl("chr[[:digit:]]",REF_1$CHROM))==0){}else{
      REF_1$CHROM<-as.factor(paste0(rep("chr",times=length(REF_1$CHROM)),REF_1$CHROM))
    }
    
    REF_GT_1<-cbind(REF_1[,1:5],GT_revise_1)
    STR<-REF_GT_1
    
    cat("REF/ALT Finished!","\n")
    
  } else{
    
    vcf_1<-read.vcfR(vcf1,convertNA = FALSE)
    vcf_2<-read.vcfR(vcf2,convertNA = FALSE)
    
    #genotype
    GT_1<-as.data.frame(vcf_1@gt,stringsAsFactors = FALSE)
    GT_1<-GT_1[,-1]
    GT_2<-as.data.frame(vcf_2@gt,stringsAsFactors = FALSE)
    GT_2_1<- as.data.frame(GT_2[,-1])
    colnames(GT_2_1) <- colnames(GT_2)[-1]
    
    #check genotype seperator
    if(!sum(grepl("[[:digit:]]/[[:digit:]]",GT_1))==0){
      GT_revise_1<-as.data.frame(matrix(data = NA,nrow = nrow(GT_1),ncol = ncol(GT_1)))
      colnames(GT_revise_1)<-colnames(GT_1)
      for (i in 1:nrow(GT_1)) {for (j in 1:ncol(GT_1)) {
        if(GT_1[i,j]=="."){GT_revise_1[i,j]<-"./."}else{GT_revise_1[i,j]<-as.character(GT_1[i,j])}
        GT_revise_1[i,j]<-gsub("/","|",GT_revise_1[i,j])
      }
      }
    }else{ GT_revise_1=GT_1}
    
    if(!sum(grepl("[[:digit:]]/[[:digit:]]",GT_2_1[,1]))==0){
      GT_revise_2<-as.data.frame(matrix(data = NA,nrow = nrow(GT_2_1),ncol = ncol(GT_2_1)))
      colnames(GT_revise_2)<-colnames(GT_2)[-1]
      for (i in 1:nrow(GT_2_1)) {for (j in 1:ncol(GT_2_1)) {
        if(GT_2_1[i,j]=="."){GT_revise_2[i,j]<-"./."}else{GT_revise_2[i,j]<-as.character(GT_2_1[i,j])}
        GT_revise_2[i,j]<-gsub("/","|",GT_revise_2[i,j])
      }
      }
    }else{GT_revise_2=GT_2_1}
    cat("Genotype Finished!","\n")
    #REF&ALT
    REF_1<-as.data.frame(vcf_1@fix,stringsAsFactors = FALSE)
    REF_1$CHROM <-as.character(REF_1$CHROM)
    if(!sum(grepl("chr[[:digit:]]",REF_1$CHROM))== 0){}else{
      REF_1$CHROM<-as.factor(paste0(rep("chr",times=length(REF_1$CHROM)),REF_1$CHROM))
    }
    REF_2<-as.data.frame(vcf_2@fix,stringsAsFactors = FALSE)
    if(!sum(grepl("chr[[:digit:]]",REF_2$CHROM))== 0){}else{
      REF_2$CHROM<-as.factor(paste0(rep("chr",times=length(REF_2$CHROM)),REF_2$CHROM))
    }
    
    REF_GT_1<-cbind(REF_1[,1:5],GT_revise_1)
    REF_GT_2<-cbind(REF_2[,1:5],GT_revise_2)
    STR<-merge(REF_GT_2,REF_GT_1,by=c("CHROM","POS","REF","ALT"),all.y = TRUE)
    colnames(STR)[which(colnames(STR)=="ID.x")]<-"ID"
    STR$ID<-STR$ID.y
    STR$ID.y<-NULL
    df<-sapply(STR[,-(1:5)],as.character)
    df[is.na(df)]<-".|."
    STR<-cbind(STR[,1:5],df)
    
    cat("REF/ALT Finished!","\n")
  }
  
  
  if(GT_mergeTable){x<-paste0(outputName,"_genotype.txt")
  write.table(STR,file = x,sep = '\t',col.names = TRUE,row.names = TRUE,quote = FALSE)
  }
  
  #Transer to Structure format
  filenames<-STR
  temp1<-as.data.frame(matrix(data = NA,ncol = (ncol(filenames)-5)*2,nrow = nrow(filenames)))
  for (i in 1:nrow(filenames)) {
    temp1[i,]<-unlist(unname(strsplit(sapply(filenames[i,-(1:5)],as.character),split="|",fixed = TRUE)))
  }
  
  temp2<-as.data.frame(matrix(data = NA,ncol = (ncol(filenames)-5)*2,nrow = nrow(filenames)))
  temp3<-as.data.frame(matrix(data = NA,ncol = (ncol(filenames)-5)*2,nrow = nrow(filenames)))
  
  #genotype convert to REF or ALT
  for (i in 1:nrow(temp1)) { for(j in 1:ncol(temp1)){
    temp2[i,j] <- sapply(switch(temp1[i,j],"1"=filenames$ALT[i],"0"=filenames$REF[i],"."="."),as.character)
  }
    if(i %%100 == 0) {print( i )}
  }
  
  #genotype convert to number
  for(i in 1:nrow(temp2)){for(j in 1:ncol(temp2)){
    
    if(nchar(temp2[i,j]) == 1){temp3[i,j] <- switch(temp2[i,j],"A"=1,"T"=2,"G"=3,"C"=4,"."=-9)}else 
      if(nchar(temp2[i,j]) > 1){temp3[i,j] <- nchar(temp2[i,j])+4 }else{}
  }
    if(i %% 100 == 0) {print( i )}
  } 
  
  if(sum(grepl("rs[[:digit:]]",filenames$ID)) == nrow(filenames)){  
    snpName <- matrix(data = filenames$ID,ncol = 1)
  }else{
    
    cat("Search snp rs number!\nThis step should be waited for a long time if more snps. \n")
    
    filenames$CHROM <- sub("chr","",filenames$CHROM)
    for(i in sort(as.numeric(unique(filenames$CHROM)))){
      chr <- paste0("ch",i)
      temp <- filenames[which(filenames$CHROM==i),]
      
      snps <- getSNPlocs(chr, as.GRanges=TRUE)
      
      for(j in 1:nrow(temp)){
        mypos <- temp$POS[j]
        idx <- match(mypos, start(snps))
        filenames$ID[which(filenames$CHROM == i)[j]] <- paste0("rs",mcols(snps)$RefSNP_id[idx])
        snpName <- matrix(data = filenames$ID, ncol = 1)
      }
      cat("chr",i,"Finished!\n")
    }
    
    
  }
  
  temp3 <- cbind(snpName,temp3)  
  cat("genotype2number Finished!","\n") 
  
  #Structure format
  
  STR_Format <- as.data.frame(t(temp3),stringsAsFactors = FALSE)
  PopGroup <- as.data.frame(matrix(data = rep(colnames(filenames[,-(1:5)]),each = 2),ncol = 1),stringsAsFactors = FALSE)
  PopGroup[,2] <- NA
  colnames(PopGroup) <- NA
  
  for (i in 1:nrow(PopGroup)) { if(sum(PopGroup[i,1] == PopCode$sample) == 1){
    PopGroup[i,2] <- sapply(switch(PopCode$super_pop[which(PopGroup[i,1] == PopCode$sample)],"AFR" = 1,"EUR" = 2,"EAS" = 3,"AMR" = 4,"SAS" = 5),as.character)
  }else{PopGroup[i,2] <- 10}
    if(i %% 100 == 0) {  print( i ) }
  }
  
  EMP <- data.frame(matrix(data = "",nrow = 1,ncol = 2),stringsAsFactors = FALSE)
  colnames(PopGroup) <- colnames(EMP)
  PopGroup <- rbind(EMP,PopGroup)
  STR_Format <- cbind(PopGroup,STR_Format)
  
  y <-paste0(outputName,"_str.txt")
  
  write.table(STR_Format,file = y,sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  cat("Population Code:AFR = 1,EUR = 2,EAS = 3,AMR = 4,SAS = 5, Other = 10 \n")
  return()
  
}


##genotype convert to PCA format
Admixture_gt2PCAformat<-function(InputFile,PopCode_data,genotype_start_column=1,genotype_seperator="|",outputName){
  if(!require(ggplot2)) stop("Missing ggplot2 package")
  if(!require(ade4)) stop("Missing ade4 package")
  
  
  #Check PopCode_data
  PopCode<-read.delim(file = PopCode_data,sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  if(is.null(PopCode_data)){stop("Missing PopCode_data")}
  if(!sum(colnames(PopCode) == "sample") == 1){stop("PopCode_data Format: no sample names")}
  if(!sum(colnames(PopCode) == "pop") == 1){stop("PopCode_data Format:no populaiton")}
  
  #check input file and its format
  STR <- read.table(InputFile,sep = '\t',stringsAsFactors = FALSE)
  if(is.null(InputFile)) stop("No inputFile")
  
  
  if(!sum(colnames(STR) == "ID"| colnames(STR) == "ID.x") == 1){
    stop("InputFiles no SNP ID.")
    
  }else {
    filenames <- STR[,-c(1:(genotype_start_column-1))]
    cat("Keep genotype format\n")}
  
  PCA_format <- as.data.frame(matrix(data = NA,ncol = ncol(filenames),nrow = nrow(filenames)),stringsAsFactors = FALSE)
  if(genotype_seperator=="|"){
    for (i in 1:nrow(filenames)) {
      for(j in 1:ncol(filenames)){
        PCA_format[i,j] <- switch(sapply(filenames[i,j],as.character),
                                  "0|0" = 0,
                                  "0|1" = 1,
                                  "1|1" = 2,
                                  "1|0" = 1,
                                  ".|." = NA)
        
      }
      if(i %%100 == 0) {print( i )}
      
    }}else{
      for (i in 1:nrow(filenames)) {for(j in 1:ncol(filenames)){
        PCA_format[i,j] <- switch(sapply(filenames[i,j],as.character),
                                  "0/0" = 0,
                                  "0/1" = 1,
                                  "1/1" = 2,
                                  "1/0" = 1,
                                  "./." = NA)
        
      }
        if(i %%100 == 0) {print( i )}
        
      }
      
      
    }
  cat("Genotype convert to number finished!","\n")
  
  
  PCA_format <- as.data.frame(t(PCA_format),stringsAsFactors = FALSE)
  colnames(PCA_format) <- STR$ID
  
  PopGroup <- as.data.frame(matrix(data = colnames(filenames),ncol = 1),stringsAsFactors = FALSE)
  PopGroup[,2] <- NA
  colnames(PopGroup) <- c("Individuals","Pop")
  
  for (i in 1:nrow(PopGroup)) { if(sum(PopGroup[i,1] == PopCode$sample) == 1){
    PopGroup[i,2] <- PopCode$super_pop[which(PopGroup[i,1] == PopCode$sample)]
  }else{PopGroup[i,2] <- "Other"}
    
    if(i %% 100 == 0) {print( i )}
    
  }
  cat("Convert to PCA_format finished!","\n")
  
  PCA_format <- cbind(PopGroup,PCA_format)
  
  write.table(PCA_format,file = paste0(outputName,"_PCAFormat.txt") ,sep = "\t", col.names = TRUE,row.names = FALSE,quote = FALSE)
  return(PCA_format)
}

Admixture_PCAplot_noMissingValue <- function(PCAformat, PopSelect = c(AFR = 1,EUR = 2,EAS = 3, AFR = 4,SAS = 5,Other = 10,"all"), figFileName){
  library(ggplot2)
  library(missMDA)
  library(ade4)
  
  if(!require(ggplot2)) stop("Missing ggplot2 package")
  if(!require(missMDA)) stop("Missing missMDA package")
  if(!require(ade4)) stop("Missing ade4 package")
  
  group.colors <- c(AFR = "#FF7300",EUR = "#00FF00",EAS = "#0000FF",Other = "black")
  
  if(sum(tolower(PopSelect) == "all") == 1){
    
    df <- read.table(PCAformat,header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df_1 <- df
    df_2 <- df_1[,-c(1,2)]
    pcal <- dudi.pca(df_2,nf = 2, scannf = FALSE) #nf:the number of kept axies
    
    scores <- pcal$li
    scores <- cbind(df_1[,1:2],scores)
    scores$Pop <- factor(df_1$Pop,levels = unique(df_1$Pop))
    
    p <- ggplot(data = scores, aes(x = Axis1, y = Axis2,color = scores$Pop)) +
      geom_point() +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      scale_colour_manual(name  = "Populations",breaks = unique(scores$Pop),values = group.colors)+
      xlab("PC1") +
      ylab("PC2") +
      ggtitle(paste0(figFileName,"_UT-250AIMs_PCA_Plot"))
    
    pdf(paste0( figFileName,"_UT-250AIMs_PCA_Plot"),paper = "letter")
    print(p)
    dev.off()
    
  }else{
    df <- read.table(PCAformat,header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df_1 <- df[which(!df$Pop == PopSelect),]
    df_2 <- df_1[,-c(1,2)]
    pcal <- dudi.pca(df_2,nf = 2, scannf = FALSE) #nf:the number of kept axies
    
    scores <- pcal$li
    scores <- cbind(df_1[,1:2],scores)
    scores$Pop <- factor(df_1$Pop,levels = unique(df_1$Pop))
    
    p <- ggplot(data = scores, aes(x = Axis1, y = Axis2,color = scores$Pop)) +
      geom_point() +
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      scale_colour_manual(name  = "Populations",breaks = unique(scores$Pop),values = group.colors)+
      xlab("PC1") +
      ylab("PC2") +
      ggtitle(paste0( figFileName,"_UT-250AIMs_PCA_Plot"))
    
    pdf(paste0( figFileName,"_UT-250AIMs_PCA_Plot"),paper = "letter")
    print(p)
    dev.off()
    
    
  }
  
  
}


Admixture_PCAplot_MissingValue <- function(PCAformat, PopSelect = c(AFR = 1, EUR = 2, EAS = 3, AFR = 4, SAS = 5,Other = 10,"all")
                                           ,  figFileName){
  library(ggplot2)
  library(missMDA)
  library(ade4)
  
  if(!require(ggplot2)) stop("Missing ggplot2 package")
  if(!require(missMDA)) stop("Missing missMDA package")
  if(!require(ade4)) stop("Missing ade4 package")
  
  group.colors <- c(AFR = "#FF7300",EUR = "#00FF00",EAS = "#0000FF",Other = "black")
  
  df <- PCAformat
  #df <- read.table(PCAformat,header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if(sum(tolower(PopSelect) == "all") == 1){
    
    #imputed missing_value
    nb <- estim_ncpPCA(df[,-c(1,2)],ncp.max=5)
    res.comp <- imputePCA(df[,-c(1,2)],ncp=2)
    res.pca <- dudi.pca(res.comp$completeObs,nf = 5, scannf = FALSE)
    
    scores<-res.pca$li
    scores<-cbind(df[,1:2],scores)
    scores$Pop<-factor(df$Pop,levels = unique(df$Pop))
    scores_1 <- scores[which(!scores$Pop == "Other"),]
    
    #plot
    p<-ggplot(data = scores_1, aes(x = Axis1, y = Axis2,color= scores_1$Pop)) +
      stat_ellipse(type = "t",linetype= 1, level = 0.95)+ 
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      scale_colour_manual(name  ="Populations",breaks = unique(scores$Pop),values = group.colors )+
      xlab("PC1")+
      ylab("PC2")+
      ggtitle(paste0( unlist(strsplit(figFileName,"/"))[2],"_UT-250AIMs_PCA_Plot"))
    
    q <- p+geom_point(data = scores[which(scores$Pop == "Other"),], color = "black" , size = 3 )+
      geom_text(data = scores[which(scores$Pop == "Other"),],aes(label = scores$Individuals[which(scores$Pop == "Other")], color = scores$Pop[which(scores$Pop == "Other")]),check_overlap = TRUE,hjust = 0, nudge_x = 0.5)
    pdf(paste0( figFileName,"_UT-250AIMs_PCA_Plot.pdf"),paper = "letter")
    print(q)
    dev.off()
    
  }else{
    
    df_1 <- df[which(!df$Pop == PopSelect),]
    #imputed missing_value
    nb <- estim_ncpPCA(df_1[,-c(1,2)],ncp.max=5)
    res.comp <- imputePCA(df_1[,-c(1,2)],ncp=2)
    res.pca <- dudi.pca(res.comp$completeObs,nf = 5, scannf = FALSE)
    
    scores<-res.pca$li
    scores<-cbind(df_1[,1:2],scores)
    scores$Pop<-factor(df$Pop,levels = unique(df_1$Pop))
    scores_1 <- scores[which(!scores$Pop == "Other"),]
    
    #plot
    p<-ggplot(data = scores_1, aes(x = Axis1, y = Axis2,color= scores_1$Pop)) +
      stat_ellipse(type = "t",linetype= 1, level = 0.95)+ 
      geom_hline(yintercept = 0, colour = "gray65") +
      geom_vline(xintercept = 0, colour = "gray65") +
      scale_colour_manual(name  ="Populations",breaks = unique(scores$Pop),values = group.colors )+
      xlab("PC1")+
      ylab("PC2")+
      ggtitle(paste0( unlist(strsplit(figFileName,"/"))[2],"_UT-250AIMs_PCA_Plot"))
    
    q <- p+geom_point(data = scores[which(scores$Pop == "Other"),], color = "black" , size = 3 )+
      geom_text(data = scores[which(scores$Pop == "Other"),],aes(label = scores$Individuals[which(scores$Pop == "Other")], color = scores$Pop[which(scores$Pop == "Other")]),check_overlap = TRUE,hjust = 0, nudge_x = 0.5)
    
    pdf(paste0( figFileName,"_UT-250AIMs_PCA_Plot.pdf"),paper = "letter")
    print(q)
    dev.off()
    
  }
  
  
}


#Structure Input CMD

Admixture_InputStr <- function(Input_Str = "str.txt",Str_Path, Output_Str = "str_output.txt"){
  path <- getwd()
  df <- read.table(Input_Str, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  People <- length(unique(df$X[-1]))
  no_AIM <- ncol(df[,-c(1,2)])
  cat("No. of people:",People,"\n")
  cat("No. of AIMs:",no_AIM,"\n")
  
  setwd(Str_Path)
  CMD<-paste0("./structure -m mainparams -e extraparams -L ",no_AIM," -N ",People," -K 3 -i ",paste(path,Input_Str, sep = "/")," -o ",paste(path,Output_Str,sep = "/"))
  system(CMD)
  
  cat(paste0(no_AIM,"AIMs"),"Structure Finished! \n")
  setwd(path)
  
}

#Structure output modify
#distinguish Pop and group

Admixture_StrOut <- function(Output_Str = "str_output.txt_f",PopGroup = c("AFR","EUR","EAS","Other"),K_value = 3, fileName = "Summary.txt"){
  StrOut_File <- read.delim(Output_Str,stringsAsFactors = FALSE)
  
  n = which(StrOut_File[,1] == " Pop       1      2      3      Individuals")
  
  Col_title = unlist(strsplit(StrOut_File[n,1],split = " "))
  
  PopCluster <- as.data.frame(matrix(data = NA,nrow = 1,ncol = length(Col_title[!nchar(Col_title) == 0])))
  for (i in 1:length(PopGroup)) {
    PopCluster[i,] <- unlist(strsplit(StrOut_File[(n+i),1],split = " "))[!nchar(unlist(strsplit(StrOut_File[(n+i),1],split = " "))) == 0]
    
  }
  
  PopCluster <- as.matrix(PopCluster[,-1])
  colnames(PopCluster)[ncol(PopCluster)] <- "Individuals"
  rownames(PopCluster) <- PopGroup
  for (i in 1:K_value) {
    colnames(PopCluster)[which(PopCluster[i,] == max(PopCluster[i,1:K_value]))] <- rownames(PopCluster)[i]
    
  }
  
  write.table(PopCluster,file = fileName, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
  return(PopCluster)
}


#seperate Frequency from structure output file
Admixture_StrOut_Frequency <- function(Output_Str,PopGroup = c("AFR","EUR","EAS","Other"),K_value = 3,PopCluster = "Summary.txt"
                                       , fileName = "Frequency.txt"){
  StrOut_File <- read.delim(Output_Str,stringsAsFactors = FALSE)
  
  m = which(StrOut_File[,1]=="        Label (%Miss) Pop:  Inferred clusters")
  
  Col_title = unlist(strsplit(StrOut_File[m,1],split = " "))[!nchar(unlist(strsplit(StrOut_File[m,1],split = " "))) == 0]
  Col_title[3] <- "Pop"
  Individual_Frequency <- as.data.frame(matrix(data = NA,nrow = 1,ncol = 5+K_value))
  
  PopCluster <- read.delim(PopCluster,stringsAsFactors = FALSE,header = TRUE)
  
  for (i in 1:sum(PopCluster$Individuals)) {
    Individual_Frequency[i,] <- unlist(strsplit(StrOut_File[(m+i),1],split = " "))[!nchar(unlist(strsplit(StrOut_File[(m+i),1],split = " "))) == 0]
    
    if(i %/% 100){print(i)}
    
  }
  
  Individual_Frequency <- Individual_Frequency[,-c(1,5)] 
  
  colnames(Individual_Frequency) <- c(Col_title[1:3],colnames(PopCluster)[1:K_value])    
  
  
  write.table(Individual_Frequency,file = fileName, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(Individual_Frequency)                     
}


##Proportion Plot
Admixture_ProportionPlot <- function(Str_frq, figSample = TRUE, sample_info,  figFileName){
  library(ggplot2)
  library(reshape)
  library(dplyr)
  group.colors <- c(AFR="#FF7300",EUR="#00FF00",EAS="#0000FF")
  
  df <- read.table(Str_frq,sep = '\t',stringsAsFactors = FALSE)
  
  if (figSample == TRUE) {
    
    #Sample proportion plot
    df_1 <- df[which(df$Pop.== 10),-c(1,3,4)]
    df_1$Label <- factor(df_1$Label,levels = df_1$Label)
    df_2 <- melt(df_1,id='Label')
    df_2$variable <- factor(df_2$variable,levels = c("AFR","EUR","EAS"))
    df_3 <- df_2[order(df_2$variable),]
    
    
    q <- ggplot(df_3, aes(fill = variable, y = value, x = Label)) + 
      geom_bar( stat ="identity", position = "fill") +
      xlab('Individuals') +
      ylab('Percentage') +
      ggtitle("UT-AIM250_ProportionPlot")+
      scale_fill_manual(breaks = c("AFR", "EUR","EAS"),values = group.colors) +
      scale_colour_discrete(name  ="Populations",breaks=c("AFR", "EUR","EAS"),labels=c("African", "European","East Asian")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    pdf(paste0(figFileName,"_ProportionPlot.pdf"),paper = "letter")
    print(q)
    dev.off() 
    
  }else{
    
    #1000Gs proportion plot
    AFR <- df[which(df$Pop.== 1),]
    AFR <- AFR[order(AFR$AFR),]
    EUR <- df[which(df$Pop.== 2),]
    EUR <- EUR[order(EUR$EUR),]
    EAS <- df[which(df$Pop.== 3),]
    EAS <- EAS[order(EAS$EAS),]
    df_1 <- rbind(AFR,EUR,EAS)
    df_1 <- df_1[,-c(1,3,4)]
    df_1$Label <- factor(df_1$Label,levels = df_1$Label)
    
    df_2 <- melt(df_1,id = 'Label')
    df_2$variable <- factor(df_2$variable,levels = c("AFR","EUR","EAS"))
    df_3 <- df_2[order(df_2$variable),]
    
    p <- ggplot(df_3, aes(fill = variable, y = value, x = Label)) + 
      geom_bar( stat ="identity", position ="fill") +
      xlab('Individuals') +
      ylab('Percentage') +
      ggtitle("UT-AIM250_1000G_ProportionPlot") +
      scale_fill_manual(breaks = c("AFR", "EUR","EAS"),values=group.colors) +
      scale_colour_discrete(name  = "Population",breaks = c("AFR", "EUR","EAS"),labels = c("African", "European","East Asian")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    pdf(paste0(figFileName,"_1000G_ProportionPlot.pdf"),paper = "letter")
    print(p)
    dev.off()
    
  }
  
}


#Triangle plot
Admixture_TrianglePlot <- function(Str_frq,figFileName = "TrianglePlot.pdf", PopCode_data){
  group.colors <- c(AFR="#FF7300",EUR="#00FF00",EAS="#0000FF",Other="black")
  library(ggtern)
  df <- Str_frq
  df_1 <- df[,-2]
  rownames(df_1) <- df_1[,1]
  df_1 <- df_1[,-1]
  df_1$Pop <- NA
  
  PopCode<-read.table(file = PopCode_data,sep = '\t',header = TRUE,stringsAsFactors = FALSE)
  
  for (i in 1:nrow(df_1)) { 
    if(sum(rownames(df_1)[i] == PopCode$sample) == 1){
      df_1$Pop[i] <- PopCode$super_pop[which(rownames(df_1)[i] == PopCode$sample)]
      
    }else{df_1$Pop[i] <- "Other"}
    if(i %% 100 == 0) {print( i)}
    
  }
  
  df_1$Pop <- factor(df_1$Pop,levels = c("AFR","EUR","EAS","Other"))
  p <- ggtern(data = df_1, aes(x = AFR, y = EAS, z = EUR)) + #define data sources
    geom_point(aes(colour = Pop)) +  #define data geometry
    scale_colour_manual(values = group.colors)+
    ggtitle(paste0( unlist(strsplit(figFileName,"/"))[2],"_trianglePlot"))
  
  
  pdf(paste0(figFileName,"_trianglePlot.pdf"),paper = "letter")
  print(p)
  dev.off()
}


