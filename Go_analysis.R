library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包
library(biomaRt)

# 读取CSV文件
file_path <- "C:/Users/Lamarck/Desktop/UP_genes_ENSEMBL_ENTREZID.csv"
gene_data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

# 提取第二列 ENTREZID
gene <- gene_data[, 2]

# GO富集分析
# GO富集的三个部分：CC、BP、MF  这里的P值和Q值都没有限定，因为后续想手动筛选
ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = TRUE)

# 把得到的四个结果转换成数据框格式
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)

# 把GO富集分析得到的完整列表以及相关信息导出
write.csv(ego_result_BP,file = "C:/Users/Lamarck/Desktop/ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "C:/Users/Lamarck/Desktop/ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "C:/Users/Lamarck/Desktop/ego_result_MF.csv",row.names = T)
