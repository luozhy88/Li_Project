

# 


# Filename="output/05_Tidy_TCM_CID_compound_Gene_zhebeimu.csv"
Filename="output/05_Tidy_TCM_CID_compound_Gene_xiakucao.csv"
DF_TCM_PDB_SMiles <- read.csv(Filename)


Nrow=nrow(DF_TCM_PDB_SMiles)
for(n in 6552:Nrow){
  # n=1
  
  SMILES <- DF_TCM_PDB_SMiles$CanonicalSMILES[n]
  identifier = DF_TCM_PDB_SMiles$identifier[n]
  CID=DF_TCM_PDB_SMiles$CID[n]
  # 按照分号分割identifier
  identifier <- unlist(strsplit(identifier, ";"))
  if(length(identifier) > 1){
    for(i in identifier ){
      # print( glue::glue(  "{SMILES}_{i} \n"   ) )
      print(n)
      run_cmd=glue::glue("Rscript 04.2_docking_Vina.R  -s \"{SMILES}\" -n {CID} -t input/PDB/{i}.pdb ")
      print(run_cmd)
      system(run_cmd)
    }
  
  
  }
}
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



