

# 


# Filename="output/05_Tidy_TCM_CID_compound_Gene_zhebeimu.csv"
# Filename="output/05_Tidy_TCM_CID_compound_Gene_xiakucao.csv"
Filename="output/05_Tidy_TCM_CID_compound_Gene_prescription.csv"
DF_TCM_PDB_SMiles <- read.csv(Filename)

##############################多线程######################
# Nrow=nrow(DF_TCM_PDB_SMiles)
# for(n in 6552:Nrow){
#   # n=1
#   
#   SMILES <- DF_TCM_PDB_SMiles$CanonicalSMILES[n]
#   identifier = DF_TCM_PDB_SMiles$identifier[n]
#   CID=DF_TCM_PDB_SMiles$CID[n]
#   # 按照分号分割identifier
#   identifier <- unlist(strsplit(identifier, ";"))
#   if(length(identifier) > 1){
#     for(i in identifier ){
#       # print( glue::glue(  "{SMILES}_{i} \n"   ) )
#       print(n)
#       run_cmd=glue::glue("Rscript 04.2_docking_Vina.R  -s \"{SMILES}\" -n {CID} -t input/PDB/{i}.pdb ")
#       print(run_cmd)
#       system(run_cmd)
#     }
#   
#   
#   }
# }
library(parallel)
library(doParallel)
library(foreach)
library(glue)

# 使用tryCatch执行并行任务
re <- tryCatch({
  # 检测核心数并创建集群
  numCores <- 30
  # numCores <- detectCores() / 10 %>% integer()
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  print(glue("numCores:{numCores}"))
  
  Nrow <- nrow(DF_TCM_PDB_SMiles)
  
  foreach(n = 6560:Nrow, .packages = c("glue")) %dopar% {
    SMILES <- DF_TCM_PDB_SMiles$CanonicalSMILES[n]
    identifier <- DF_TCM_PDB_SMiles$identifier[n]
    CID <- DF_TCM_PDB_SMiles$CID[n]
    
    # 按照分号分割identifier
    identifier <- unlist(strsplit(identifier, ";"))
    
    if (length(identifier) > 1) {
      for (i in identifier) {
        print(n)
        cat(n)
        run_cmd <- glue("Rscript 04.2_docking_Vina.R -s \"{SMILES}\" -n {CID} -t input/PDB/{i}.pdb")
        print(run_cmd)
        cat(run_cmd)
        system(run_cmd)
        system(glue::glue("echo {n} >> output/doing.txt" ) )
      }
    }
  }
  
}, error = function(e) {
  message("捕获到错误: ", e$message)
  NULL  # 返回NULL或其他默认值
}, finally = {
  # 确保集群在任务完成或出现错误时被释放
  stopCluster(cl)
  message("集群已停止")
})

##############################多线程######################

dir.create("output/Docking_Vina")
system("mv output/results* output/Docking_Vina")

#################################tidy result########################
tbl_all <-
  list.files(path = "output/Docking_Vina/", recursive = T, pattern = ".*score.csv", full.names = TRUE) %>%
  purrr::set_names() %>%
  purrr::map_dfr(read.csv)
tbl_all$vina_results =gsub("REMARK VINA RESULT:    ","",tbl_all$vina_results)
tbl_all$vina_results =gsub("    .*","",tbl_all$vina_results)
tbl_all$ligand= gsub("output/Docking_Vina/results_Vina_","",tbl_all$Path)
tbl_all$ligand= gsub("_.*","",tbl_all$ligand)
tbl_all$receptor=gsub(".pdb_ligand_vina_out_score.csv","",tbl_all$Path)
tbl_all$receptor=gsub(".*_","",tbl_all$receptor)
tbl_all$vina_results =as.numeric(tbl_all$vina_results)

tbl_all=tbl_all %>% data.frame() %>% dplyr::distinct()
write.csv(tbl_all,"output/Docking_Vina/06_Docking_Vina_result_ALL_scores.csv",row.names = F,quote = T)



