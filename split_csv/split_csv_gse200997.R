setwd("/home/r2310101419")#我的服务器的工作地址
getwd()

#提前让后台把包加载好，把文件该解压的让后台解压出来
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
#library(muscat)
library(reshape2)
library(celldex)
library(clustree)
library(SingleR)
library(tidyverse)
library(RColorBrewer)

#读取我的两个文件（annotation和UMI的raw data）
datadir = "GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv"
metadata = "GSE200997_GEO_processed_CRC_10X_cell_annotation.csv"
outdir = "./Dimention"
nfeature<-100
ncount<-200
mito<-10#线粒体
ribo<-100#核糖体
red<-10
# read gene-sampleID matrix and metadata
data<-fread(datadir, header = T, check.names = F)
#我的是Excel，为已经分开的版本；如果是txt版本，则需要注明制表符，用这句：
#data<-fread(datadir,sep = "\t", header = T, check.names = F)
data<-data.frame(data,check.names = F,row.names = 1)

metadata<-fread(metadata, header = T, check.names = F)#本句制表符同上
metadata<-data.frame(metadata,check.names = F,row.names = 1)

#创建SeuratObject
#装包 install.packages("Seurat")
library(Seurat)
pbmc = CreateSeuratObject(counts = data,metadata = metadata)
pbmc = AddMetaData(object = pbmc, metadata = metadata)
#开始数据拆分
a <-pbmc
library(tidyverse) #加载tidyverse
head(a@meta.data) #查看读取的数据前6行
#head(a@meta.data, 8) #查看前8行的命令
#正式拆分，使用tidyverse函数
#rownames(a@meta.data) -> a@meta.data$barcode 
#$用于从dataframe或list里面提取某个变量,上一句是把所有的行名赋值给一个新的列变量barcode，让所有列都有名字


R1T<- SplitObject(a,split.by='samples') #引号里具体的名字一定和自己文件里对应好
R1T#看一下RT1，拆分是否成功，再把拆出来的用$单独提取,输入$会自动给出可提取的列表，选择即可
C1RT<-R1T$T_cac1
C2LT<-R1T$T_cac2
C3RT<-R1T$T_cac3
C4LT<-R1T$T_cac4
C5RT<-R1T$T_cac5
C6LT<-R1T$T_cac6
C7RT<-R1T$T_cac7
C8RT<-R1T$T_cac8
C9LT<-R1T$T_cac9
C10LT<-R1T$T_cac10
C11RT<-R1T$T_cac11
C12RT<-R1T$T_cac12
C13RT<-R1T$T_cac13
C14LT<-R1T$T_cac14
C15LT<-R1T$T_cac15
C16LT<-R1T$T_cac16
C4LN<-R1T$B_cac4
C6LN<-R1T$B_cac6
C7RN<-R1T$B_cac7
C10LN<-R1T$B_cac10
C11RN<-R1T$B_cac11
C14LN<-R1T$B_cac14
C15LN<-R1T$B_cac15

#提取好的保存成各自的文件
saveRDS(C1RT, file = "C1RT.Rds")
saveRDS(C2LT, file = "C2LT.Rds")
saveRDS(C3RT, file = "C3RT.Rds")
saveRDS(C4LT, file = "C4LT.Rds")
saveRDS(C5RT, file = "C5RT.Rds")
saveRDS(C6LT, file = "C6LT.Rds")
saveRDS(C7RT, file = "C7RT.Rds")
saveRDS(C8RT, file = "C8RT.Rds")
saveRDS(C9LT, file = "C9LT.Rds")
saveRDS(C10LT, file = "C10LT.Rds")
saveRDS(C11RT, file = "C11RT.Rds")
saveRDS(C12RT, file = "C12RT.Rds")
saveRDS(C13RT, file = "C13RT.Rds")
saveRDS(C14LT, file = "C14LT.Rds")
saveRDS(C15LT, file = "C15LT.Rds")
saveRDS(C16LT, file = "C16LT.Rds")
saveRDS(C4LN, file = "C4LN.Rds")
saveRDS(C6LN, file = "C6LN.Rds")
saveRDS(C7RN, file = "C7RN.Rds")
saveRDS(C10LN, file = "C10LN.Rds")
saveRDS(C11RN, file = "C11RN.Rds")
saveRDS(C14LN, file = "C14LN.Rds")
saveRDS(C15LN, file = "C15LN.Rds")