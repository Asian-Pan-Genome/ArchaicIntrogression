library(geneHapR)
library(conflicted)
library(dplyr)
library(tidyverse)

# 明确指定优先使用dplyr的select函数
conflict_prefer("select", "dplyr")
conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::rename)

## work1: hapResult hapSummary hapNet ####
target_genes <- c("PRDM16")
target_genes <- c("PRKCH")
target_genes <- c("CSGALNACT2")
target_genes <- c("CH25H")
###============================================================================= parameter 循环法调用参数
popfile <- "/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/007_ASMai/000_popInfo/003_PopInfo_APG_addArchaic.txt"
# 读取群体信息文件（只需要读一次）
dfpop <- read_tsv(popfile) %>% select(Assembly_ID, Pop, Anc)

dfpara <- readxl::read_xlsx("/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/011_adaptiveIntrogression/005_geneHaplotype/000_gene_hsmap/gene_block_hsmap.xlsx")

# 定义基础路径
# popfile <- "/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/007_ASMai/000_popInfo/003_PopInfo_APG_addArchaic.txt"
base_hapmap_path <- "/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/011_adaptiveIntrogression/001_DenOri"
base_output_path <- "/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/012_haplotype_network/001_haplotype"
###============================================================================= parameter end

df_target <- dfpara %>% filter(GeneName %in% target_genes)
dfpara0 <- df_target

for(i in 1:nrow(dfpara0)) {
  # 获取当前基因的参数
  geneName <- dfpara0$GeneName[i]
  DenBlock <- str_sub(dfpara0$ChunkID[i],2) 
  archaic <- dfpara0$Archaic[i]
  
  Start_heatmap <- dfpara0$Start_heatmap
  End_heatmap <- dfpara0$End_heatmap
  
  cat("正在处理基因:", geneName, "，区块:", DenBlock, "\n")
  
  if(archaic == "Denisovan"){
    # 构建文件路径
    hapmapfile <- file.path(base_hapmap_path, paste0("Den_", DenBlock, ".merge.easy_region.filt_chr.hapmap.xlsx"))
    outfile_hapResult <- file.path(base_output_path, paste0(geneName,"_D",DenBlock, "_forhapResult.txt"))
    outfile_hapSummary <- file.path(base_output_path, paste0(geneName,"_D",DenBlock, "_forhapSummary.txt"))
    
  }
  if(archaic == "Neanderthal"){
    # 构建文件路径
    hapmapfile <- file.path(base_hapmap_path, paste0("Nean_", DenBlock, ".merge.easy_region.filt_chr.hapmap.xlsx"))
    outfile_hapResult <- file.path(base_output_path, paste0(geneName,"_N",DenBlock, "_forhapResult.txt"))
    outfile_hapSummary <- file.path(base_output_path, paste0(geneName,"_N",DenBlock, "_forhapSummary.txt"))
  }
  
  # 检查hapmap文件是否存在
  if(!file.exists(hapmapfile)) {
    cat("  ⚠️  文件不存在，跳过:", hapmapfile, "\n")
    next  # 跳过当前循环，继续下一个基因
  }
  
  cat("  ✅ 文件存在，开始处理...\n")  
  
  ### ============================== step1: 处理原始文件，获得 geno 数据框
  df0 <- read_tsv(hapmapfile) %>% 
    rename(Den = `Den-Den`, NeanAltai = `Nean-Altai`, NeanVin = `Nean-Vin`, NeanCha = `Nean-Cha`) %>% 
    arrange(POS)
  
  ### 数据过滤，第一步：根据条件过滤POS
  df1 <- df0 %>% 
    filter(str_length(REF) == 1, str_length(ALT) ==1) %>%
    select(-POS_HG38) %>% 
    pivot_longer(cols = `KOR01-Mat`:last_col(), names_to = "Assembly_ID", values_to = "Geno") %>% 
    left_join(dfpop,by=c("Assembly_ID")) %>% 
    ############# 第一种过滤方法：频率过滤
    ### 根据频率进行过滤,求每个位置在所有个体中的频率大于0.05
    group_by(POS) %>%
    summarise(geno_sum = sum(Geno), n = 2*n(), fre = geno_sum/n) %>%
    filter(0.05 < fre, fre < 0.95 ) %>%
    ungroup() %>%
    select(POS)
  
  df2 <- df0 %>% 
    inner_join(df1,by = c("POS"))
  
  ### 数据过滤，第二步：过滤样本
  df3 <- df2 %>% 
    pivot_longer(cols = `KOR01-Mat`:last_col(), names_to = "Assembly_ID", values_to = "Geno") %>% 
    left_join(dfpop,by=c("Assembly_ID")) %>% 
    # filter(Anc != "AFR") %>% 
    select(-c(Pop,Anc)) %>% 
    ### 转置回原来的摸样
    pivot_wider(names_from = Assembly_ID, values_from = Geno)
  
  ### 过滤完毕后，准备pheatmap所需的文件格式
  if(!is.na(Start_heatmap)){
    if(archaic == "Denisovan"){
      df4 <- df3 %>% 
        select(CHR = `#CHROM`, POS, REF, ALT, POS_HG38,c(Den:last_col())) %>% 
        ### POS 是CHM13的坐标
        filter(POS >= Start_heatmap, POS <= End_heatmap) 
    }else {
      
      df4 <- df3 %>% 
        select(CHR = `#CHROM`, POS, REF, ALT, POS_HG38,c(NeanAltai,`KOR01-Mat`:last_col())) %>% 
        ### POS 是CHM13的坐标
        filter(POS >= Start_heatmap, POS <= End_heatmap) 
    }
    
  }else{
    if(archaic == "Denisovan"){
      df4 <- df3 %>% 
        select(CHR = `#CHROM`, POS, REF, ALT, POS_HG38,c(Den:last_col()))
    }else{
      df4 <- df3 %>% 
        select(CHR = `#CHROM`, POS, REF, ALT, POS_HG38,c(NeanAltai,`KOR01-Mat`:last_col()))
    }
  }
  
  
  ### write_tsv(df4,"~/Documents/geon.txt" )
  ### step1: import data
  df <- df4 %>% mutate(across(6:last_col(), ~ case_when(
    .x == 0 ~ REF,
    .x == 1 ~ ALT,
    .x == 2 ~ ALT,
    TRUE ~ "N"  # # 将其他所有值都替换为字符"N"
  )))
  
  ################################################################################################ 过滤POS
  ### ==============过滤一下SNP数量，用小数量的数据来构建 haplotype network
  # dfintro <- read_tsv("/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/011_adaptiveIntrogression/005_geneHaplotype/001_gene/002_introgressionSite.txt")
  # 获取列名并创建数据框
  # df <- df %>% 
  #   filter(POS_HG38 %in% dfintro$POS_HG38)
  ### ==============end
  ### write_tsv(df,"~/Documents/geon_ATCG.txt" )
  ####============================================================================================ delimiter

  
  ### ============================== step2: 建立 gt.geno AccINFO
  gt.geno <- as.data.frame(df)
  
  #### 注释信息
  dfpop <- read_tsv("/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/007_ASMai/000_popInfo/003_PopInfo_APG_addArchaic.txt") %>% 
    select(Assembly_ID, Anc)
  
  AccINFO0 <- tibble(id = colnames(df)) %>% dplyr::slice(6:n()) %>% left_join(dfpop,by=c("id" = "Assembly_ID")) 
  AccINFO <- AccINFO0 %>% column_to_rownames("id")
  ### write_tsv(AccINFO0,"~/Documents/pop_group.txt" )
  
  # Import accession group/location information
  # AccINFO <- import_AccINFO("~/Documents/t.txt")
  ################################################################################ 构建对象 hapResult hapSummary
  ### step2: table2hap
  hapResult <- table2hap(gt.geno, hapPrefix = "H",
                         hetero_remove = TRUE, na_drop = TRUE) 
  # write_tsv(hapResult,"/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/012_haplotype_network/001_test/001_hapResult.txt")
  
  write_tsv(hapResult,outfile_hapResult)
  # ### 可临时添加其他信息 INFO
  # hapResult <- addINFO(hapResult,
  #                      tag = "PrChange",
  #                      values = c("C->D", "V->R", "G->N","C->D","C->D","C->D"),
  #                      replace = F) # replace the old INFO or not
  
  ### 查看单倍型结果
  hapSummary <- geneHapR::hap_summary(hapResult)
  # write_tsv(hapSummary,"/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/012_haplotype_network/001_test/001_hapSummary.txt")
  write_tsv(hapSummary,outfile_hapSummary)
  
  ### 获取每个单倍型的频率表格
  ### df_hap_fre <- tibble(Hap = hapSummary$Hap, Fre = hapSummary$freq) %>% filter(str_starts(Hap,"H"))
  ### 获取每个样本所属的单倍型
  df_sm_hap <- tibble(Hap = hapResult$Hap, Accession = hapResult$Accession) %>% filter(str_starts(Hap,"H"))
  write_tsv(df_sm_hap,str_replace(outfile_hapResult,".txt","_sm2hap_hsmap.txt"))
  
  ### 获取每个单倍型的频率表格
  df_sm_hap2 <- df_sm_hap %>% group_by(Hap) %>% dplyr::count()
  write_tsv(df_sm_hap2,str_replace(outfile_hapSummary,".txt","_hap2fre_hsmap.txt"))
  
  ##write_tsv(df_sm_hap,"/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/012_haplotype_network/001_haplotype/TRPV1_sm2hap_hsmap.txt")
  
  ### 可视化详细的变异信息
  # plotHapTable(hapSummary)
  
  # Plot variants details with additional information
  # plotHapTable(hapSummary,
  #              hapPrefix = "H",
  #              INFO_tag = c("PrChange"), ###这条信息也显示出来
  #              tag_name = c("Pr"), ###自己可以随便命名
  #              displayIndelSize = 1, 
  #              angle = 0, ##or 45
  #              replaceMultiAllele = TRUE,
  #              ALLELE.color = "grey90")
  
  ################################################################################ 直接构建 hapNet
  ### Generate haplonet
  hapNet <- get_hapNet(hapSummary,
                       AccINFO = AccINFO,
                       groupName = "Anc")
  
  ### 查看 hapNet 类中的内容：
  #print.default(hapNet)
  
  hap_group_table <- attr(hapNet, "hapGroup")
  df_hap_group_table <- hap_group_table %>% as.data.frame() %>% pivot_wider(names_from = group,values_from = Freq)
  write_tsv(df_hap_group_table, str_replace(outfile_hapSummary,".txt","_hap2group.txt"))
  
  # 开始重定向输出到文件
  sink(str_replace(outfile_hapSummary,".txt","_hapNet_output.txt"))
  print.default(hapNet)
  sink()  # 结束重定向
  
  
} ### for 循环结束

####============================================================================================ delimiter
## work2: 计算两两 haplotype 之间的差异-不过滤hap ####

#输出结果为：差异表格 和热图
#write_csv(diff_table, file.path(base_output_path,str_c(geneName,"_",DenBlock,"_haplotype_pairwise_differences.csv")))

library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(viridis)

###============================================================================= parameter 循环法调用参数
popfile <- "/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/007_ASMai/000_popInfo/003_PopInfo_APG_addArchaic.txt"
# 读取群体信息文件（只需要读一次）
dfpop <- read_tsv(popfile) %>% select(Assembly_ID, Pop, Anc)
dfpara <- readxl::read_xlsx("/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/011_adaptiveIntrogression/005_geneHaplotype/000_gene_hsmap/gene_block_hsmap.xlsx")

# 定义基础路径
base_hapmap_path <- "/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/012_haplotype_network/001_haplotype"
base_output_path <- "/Users/Aoyue/project_pos/APG_archaicDNA/002_dataAnalysis/012_haplotype_network/002_plot"
###============================================================================= parameter end

dfpara0 <- dfpara %>% 
  filter(GeneName %in% target_genes) %>% 
  mutate(a=1)

for(i in 1:nrow(dfpara0)) {
  # 获取当前基因的参数
  geneName <- dfpara0$GeneName[i]
  DenBlock <- str_sub(dfpara0$ChunkID[i],2) 
  DenBlockND <- dfpara0$ChunkID[i]
  archaic <- dfpara0$Archaic[i]
  
  Start_heatmap <- dfpara0$Start_heatmap
  End_heatmap <- dfpara0$End_heatmap
  
  cat("正在处理基因:", geneName, "，区块:", DenBlock, "\n")
  
  # ==================== 第一步：数据预处理 ====================
  hap_geno <- hapSummary %>% 
    dplyr::slice(-c(1:4)) %>%  # 移除CHR, POS, INFO, ALLELE行
    select(-Accession, -freq)  # 移除Accession和freq列
  
  # 设置行名为单倍型名称
  rownames(hap_geno) <- hap_geno$Hap
  hap_geno <- hap_geno %>% select(-Hap)
  
  # ==================== 第二步：计算两两差异 ====================
  calculate_pairwise_differences <- function(geno_df) {
    hap_names <- rownames(geno_df)
    n_haps <- length(hap_names)
    
    # 创建空矩阵存储差异
    diff_matrix <- matrix(0, nrow = n_haps, ncol = n_haps,
                          dimnames = list(hap_names, hap_names))
    
    # 计算每对单倍型间的差异
    for (i in 1:(n_haps-1)) {
      for (j in (i+1):n_haps) {
        hap1 <- as.character(geno_df[i, ])
        hap2 <- as.character(geno_df[j, ])
        
        # 计算差异位点数（排除NA）
        differences <- sum(hap1 != hap2, na.rm = TRUE)
        diff_matrix[i, j] <- differences
        diff_matrix[j, i] <- differences  # 对称矩阵
      }
    }
    
    return(diff_matrix)
  }
  
  # 计算差异矩阵
  diff_matrix <- calculate_pairwise_differences(hap_geno)
  
  df_diff_matrix <- as.data.frame(diff_matrix) %>% rownames_to_column(var = "Hap")
  writexl::write_xlsx(df_diff_matrix,file.path(base_output_path,str_c(DenBlockND,"_",geneName,"_haplotype_pairwise_differences_table.xlsx")))
  # ==================== 第三步：输出差异表格 ====================
  # 转换为长格式并保存
  diff_table <- diff_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("Hap1") %>% 
    pivot_longer(cols = -Hap1, names_to = "Hap2", values_to = "Differences") %>% 
    filter(Hap1 != Hap2)  # 移除对角线（自身比较）
  
  # 保存差异表格
  write_csv(diff_table, file.path(base_output_path,str_c(DenBlockND,"_",geneName,"_haplotype_pairwise_differences.csv")))
  cat("差异表格已保存为: haplotype_pairwise_differences.csv\n")
  
  # 查看摘要统计
  cat("差异统计摘要:\n")
  cat("最小差异:", min(diff_table$Differences), "\n")
  cat("最大差异:", max(diff_table$Differences), "\n")
  cat("平均差异:", mean(diff_table$Differences), "\n")
  
  
  # diff_table_Archaic2Morden <- diff_table %>% 
  #   filter(Hap1 == "H070")
  
  
  # ==================== 第四步：多种可视化方法 ====================

  pheatmap(diff_matrix,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.0f",
           number_color = "black",
           # fontsize_number = 8,
           fontsize_number = 4, ## prdm16
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = "The full matrix of pairwise differences between all the unique haplotypes",
           # filename = "haplotype_differences_clustered.pdf",
           width = 15,height = 10, ## prdm16
           # width = 10,height = 8, 
           filename = file.path(base_output_path,str_c(geneName,"_",DenBlockND,"_haplotype_differences_clustered.pdf"))
  )
  cat("聚类热图已保存为: haplotype_differences_clustered.pdf\n")
  
  
  # 网络图可视化差异关系
  library(igraph)
  
  create_difference_network <- function(diff_table, threshold = 5) {
    # 只保留差异较小的连接
    network_edges <- diff_table %>% 
      filter(Differences <= threshold) %>% 
      select(from = Hap1, to = Hap2, weight = Differences)
    
    # 创建图对象
    g <- graph_from_data_frame(network_edges, directed = FALSE)
    
    # 设置绘图参数
    plot(g, 
         vertex.size = 15,
         vertex.color = "lightblue",
         vertex.label.cex = 0.8,
         vertex.label.color = "black",
         edge.width = 2,
         edge.color = "gray",
         main = paste("Haplotype Similarity Network (Differences ≤", threshold, ")"))
  }
  
  # pdf("haplotype_similarity_network.pdf", width = 10, height = 8)
  pdf(file.path(base_output_path,str_c(geneName,"_",DenBlock,"_differences_distribution.pdf")),width = 10, height = 8)
  create_difference_network(diff_table, threshold = 5)
  dev.off()
  cat("网络图已保存为: haplotype_similarity_network.pdf\n")
  
  # ==================== 第五步：生成分析报告 ====================
  cat("\n=== 分析报告 ===\n")
  cat("分析的单倍型数量:", nrow(diff_matrix), "\n")
  cat("分析的SNP位点数量:", ncol(hap_geno), "\n")
  cat("总的两两比较数:", nrow(diff_table), "\n")
  cat("差异范围:", min(diff_table$Differences), "-", max(diff_table$Differences), "个位点\n")
  cat("平均遗传差异:", round(mean(diff_table$Differences), 2), "±", 
      round(sd(diff_table$Differences), 2), "个位点\n")
  
  # 找出差异最小和最大的单倍型对
  min_diff_pair <- diff_table[which.min(diff_table$Differences), ]
  max_diff_pair <- diff_table[which.max(diff_table$Differences), ]
  
  cat("\n最相似的单倍型对:", min_diff_pair$Hap1, "vs", min_diff_pair$Hap2, 
      "(差异 =", min_diff_pair$Differences, "个位点)\n")
  cat("最差异的单倍型对:", max_diff_pair$Hap1, "vs", max_diff_pair$Hap2, 
      "(差异 =", max_diff_pair$Differences, "个位点)\n")
  
  # 保存详细报告
  sink(file.path(base_output_path,str_c(geneName,"_",DenBlock,"_haplotype_analysis_report.txt")))
  cat("单倍型差异分析报告\n")
  cat("===================\n\n")
  cat("分析时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("基本统计:\n")
  cat("- 单倍型数量:", nrow(diff_matrix), "\n")
  cat("- SNP位点数量:", ncol(hap_geno), "\n")
  cat("- 差异范围:", min(diff_table$Differences), "-", max(diff_table$Differences), "\n")
  cat("- 平均差异:", round(mean(diff_table$Differences), 2), "\n\n")
  cat("最相似的单倍型对:", paste(min_diff_pair$Hap1, min_diff_pair$Hap2, sep = " vs "), 
      "(差异 =", min_diff_pair$Differences, ")\n")
  cat("最差异的单倍型对:", paste(max_diff_pair$Hap1, max_diff_pair$Hap2, sep = " vs "), 
      "(差异 =", max_diff_pair$Differences, ")\n")
  sink()
  
  cat("\n分析完成！生成的文件:\n")
  cat("1. haplotype_pairwise_differences.csv - 差异数据表格\n")
  cat("2. haplotype_differences_heatmap.pdf - 热图\n")
  cat("3. haplotype_differences_clustered.pdf - 聚类热图\n")
  cat("4. haplotype_differences_upper_triangle.pdf - 上三角热图\n")
  cat("5. differences_distribution.pdf - 差异分布图\n")
  cat("6. haplotype_similarity_network.pdf - 相似性网络图\n")
  cat("7. haplotype_analysis_report.txt - 分析报告\n")
  
}

####============================================================================================ delimiter
## work3: 计算两两 haplotype 之间的差异-不过滤hap ####



