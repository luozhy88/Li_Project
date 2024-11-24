

# 



DF_TCM_PDB_SMiles <- read.csv("output/05_Tidy_TCM_CID_compound_Gene_zhebeimu.csv")


Nrow=nrow(DF_TCM_PDB_SMiles)
for(n in 1:5){
  # n=347
  print(n)
  SMILES <- DF_TCM_PDB_SMiles$CanonicalSMILES[n]
  identifier = DF_TCM_PDB_SMiles$identifier[n]
  CID=DF_TCM_PDB_SMiles$CID[n]
  # 按照分号分割identifier
  identifier <- unlist(strsplit(identifier, ";")) %>% unique()
  if(length(identifier) >= 1){
    for(i in identifier ){
      # print( glue::glue(  "{SMILES}_{i} \n"   ) )
      run_cmd=glue::glue("Rscript 04_docking_Vina.R  -s \"{SMILES}\" -n {CID} -t input/PDB/{i}.pdb ")
      print(run_cmd)
      system(run_cmd)
    }
  
  
  }
} 
dir.create("output/Docking_Vina")
system("mv output/results* output/Docking_Vina")


