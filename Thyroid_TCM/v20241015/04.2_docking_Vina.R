
rm(list = ls())
####################
#如果mySMILES为非空，则优先执行mySMILES，否则执行Mol.path
# Rscript 01_docking_Vina.R -m input/test1/P06.mol2 -s "CC(C)CC(NC(=S)Nc1ccccc1)C(=O)NC(CO)CCO" -n 5hie -t input/test1/5hie_modified.pdb -M Vina
# Rscript 04.2_docking_Vina.R  -s "CC(C)CC(NC(=S)Nc1ccccc1)C(=O)NC(CO)CCO" -n 5hie -t input/test1/5hie_modified.pdb 
# Rscript 04_docking_Vina.R  -s "CCC(CCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C)C(C)C" -n 7Q38 -t input/PDB/7Q38.pdb
# Rscript 04.2_docking_Vina.R  -s "CC(C)CC(NC(=S)Nc1ccccc1)C(=O)NC(CO)CCO" -n 7Q38 -t input/PDB/7Q38.pdb
# Rscript 04_docking_Vina.R  -s "CCC(CCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C)C(C)C" -n 7Q38 -t input/test1/5hie_modified.pdb 

# conda install conda-forge::vina
# conda install hcc::adfr-suite
# conda install conda-forge::meeko
# git clone git@github.com:forlilab/scrubber.git
# cd scrubber
# pip install -e .
# conda install rdkit -c conda-forge

##################
library(optparse)
library(dplyr)

# 定义命令行参数
option_list <- list(
  make_option(c("-m", "--mol"), type = "character", default = "input/test1/P06.mol2",
              help = "Path to the mol2 file"),
  make_option(c("-s", "--smiles"), type = "character", 
              default = "CC(C)CC(NC(=S)Nc1ccccc1)C(=O)NC(CO)CCO",
              help = "SMILES string"),
  make_option(c("-n", "--name"), type = "character", default = "5hie",
              help = "Name of the SMILES"),
  make_option(c("-t", "--target"), type = "character", 
              default = "input/test1/5hie_modified.pdb",
              help = "Path to the target PDB file"),
  make_option(c("-M", "--method"), type = "character", default = "Vina",
              help = "Docking method")
)

# 解析命令行参数
opt <- parse_args(OptionParser(option_list=option_list))

# 使用解析后的参数
Mol.path <- opt$mol
mySMILES <- opt$smiles  
mySMILES_name <- opt$name
target.path <- opt$target
Method <- opt$method






####################
# 1. Preparing the ligand and target

# Mol.path="input/test1/P06.mol2"
# mySMILES="CCC(CCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C)C(C)C"
# mySMILES_name="7Q38"
# target.path="input/PDB/7Q38.pdb" #https://www.rcsb.org/structure/8KHF
# Method="Vina"

# Mol.path="input/test1/P06.mol2"
# mySMILES="CC(C)CC(NC(=S)Nc1ccccc1)C(=O)NC(CO)CCO"
# mySMILES_name="5hie"
# target.path="input/test1/5hie_modified.pdb" #https://www.rcsb.org/structure/8KHF
# Method="Vina"
################################get boxCenter######################################
#2. 计算calculate_box_center
calculate_box_center <- function(pdb_file_path) {
  # 读取PDB文件
  # pdb_file_path=target.path
  pdb_lines <- readLines(pdb_file_path)
  
  # 提取ATOM行
  atom_lines <- pdb_lines[grep("^ATOM", pdb_lines)]
  
  # 提取坐标信息
  coords <- do.call(rbind, strsplit(substr(atom_lines, 31, 54), "\\s+"))
  
  coords <- as.data.frame(apply(coords, 2, as.numeric))
  # 筛选非空值的列
  coords <- coords[, colSums(!is.na(coords)) > 0]
  # 删除含有空值的行
  # coords <- coords[rowSums(coords == "") == 0, ]
  coords=coords[,1:3]
  
  colnames(coords) <- c("x", "y", "z")
  
  coords <- as.data.frame(apply(coords, 2, as.numeric))
  coords=coords[!is.na(coords$x),]
  # 计算box_center
  box_center <- colMeans(coords)
  
  # 返回结果
  return(box_center)
}

# 使用函数
# pdb_file_path <- "input/test4/1P9U.pdb"
box_center <- calculate_box_center(target.path)
box_center= round(box_center,0)
boxCenter=paste(box_center[1],box_center[2],box_center[3],sep = "_")
receptor_vina_box= glue::glue("
center_x = {box_center[1]}
center_y = {box_center[2]}
center_z = {box_center[3]}
size_x = 20.0
size_y = 20.0
size_z = 20.0
")
print(paste0("boxCenter: ",receptor_vina_box))
dir.create("output/receptor_vina_box", showWarnings = FALSE)
system("rm -rf output/receptor_vina_box.txt")
target.receptor=basename(target.path) %>% gsub(".pdb","",.) 
writeLines(receptor_vina_box,glue::glue("output/receptor_vina_box/{target.receptor}_receptor_vina_box.txt") , sep = "\n")




######################### Preparing the sdf by mySMILES#########################
dir.create("output/04.2_SDF",recursive = TRUE, showWarnings = FALSE)
dir.create("output/04.2_PDB",recursive = TRUE, showWarnings = FALSE)
dir.create("output/Docking_Vina",recursive = TRUE, showWarnings = FALSE)
system(  glue::glue("/home/zhiyu/miniconda3/envs/vina/bin/scrub.py '{mySMILES}' -o output/04.2_SDF/{mySMILES_name}.sdf ")  )

################################### Preparing the ligand########################
system(  glue::glue("/home/zhiyu/miniconda3/envs/meeko/bin/mk_prepare_ligand.py  -i output/04.2_SDF/{mySMILES_name}.sdf -o  output/04.2_SDF/{mySMILES_name}_ligand.pdbqt ")  )

###################################### Preparing the receptor################### 
PDB.pdbqt=basename(target.path) %>% gsub(".pdb","_receptor.pdbqt",.)
system(  glue::glue("/home/zhiyu/miniconda3/envs/vina/bin/prepare_receptor -r {target.path} -o  output/04.2_PDB/{PDB.pdbqt}")  )

################################### run #######################################
if( file.exists(glue::glue("output/04.2_PDB/{PDB.pdbqt}")  )  & file.exists(glue::glue("output/04.2_SDF/{mySMILES_name}_ligand.pdbqt") ) ){
  cmd_run=glue::glue( "/home/zhiyu/miniconda3/envs/vina/bin/vina --receptor output/04.2_PDB/{PDB.pdbqt} --ligand output/04.2_SDF/{mySMILES_name}_ligand.pdbqt  --config output/receptor_vina_box/{target.receptor}_receptor_vina_box.txt --exhaustiveness=32 --out output/Docking_Vina/results_{Method}_{make.names(mySMILES_name)}_{basename(target.path)}_ligand_vina_out.pdbqt --cpu 5 "  )
  print(cmd_run)
  system( cmd_run  )
}
#####################################get scores#################################

file_path <- glue::glue("output/Docking_Vina/results_{Method}_{make.names(mySMILES_name)}_{basename(target.path)}_ligand_vina_out.pdbqt")
if (file.exists(file_path)) {
  lines <- readLines(file_path)
  
  # 提取含有"REMARK VINA RESULT"的行
  vina_results <- lines[grep("REMARK VINA RESULT", lines)]
  vina_results_df=data.frame(vina_results = vina_results)
  vina_results_df$Model=paste0("Model_",rownames(vina_results_df))
  # 打印结果
  print(vina_results_df)
  Out.name=gsub(".pdbqt","_score.csv",file_path)
  # glue::glue("output/results_{Method}_{make.names(basename(ligand))}_{basename(target.path)}_score.csv")
  vina_results_df$Path=Out.name
  write.csv(vina_results_df,Out.name,row.names = FALSE)
}







# mk_prepare_ligand.py -i {mySMILES_name}.sdf -o 1iep_ligand.pdbqt

# 
# # 3. bash
# if(mySMILES !=""){
#   print("You are using mySMILES!")
#   run1=glue::glue(" curl -g 'https://swissdock.ch:8443/preplig?mySMILES={mySMILES}&{Method}' ")
#   ligand=mySMILES_name
# } else {
#   run1=glue::glue(" curl -F 'myLig=@{Mol.path}' 'https://swissdock.ch:8443/preplig?{Method}' ")
#   ligand=Mol.path
# }
# 
# 
# output <- system(run1, intern = TRUE)# 运行命令并捕获输出
# session_number <- as.numeric(gsub(".*Session number: (\\d+).*", "\\1", paste(output, collapse = " ")))# 使用正则表达式提取session number
# print(paste0("Now,Session number: ", session_number))
# 
# run2=glue::glue(" curl -F 'myTarget=@{target.path}' 'https://swissdock.ch:8443/preptarget?sessionNumber={session_number}'  "    )
# print(run2)
# system(run2)# 运行命令并捕获输出
# 
# run3=glue::glue("curl https://swissdock.ch:8443/setparameters\\?sessionNumber\\={session_number}\\&exhaust=4\\&ric=5\\&boxCenter={boxCenter}\\&boxSize=20_20_20 ")#\\&cavity=70\\ &boxCenter=-67_2_87
# # run3=glue::glue("curl https://swissdock.ch:8443/setparameters")
# print(run3)
# system(run3)# 配置参数
# 
# run4=glue::glue("curl https://swissdock.ch:8443/checkstatus\\?sessionNumber\\={session_number} ")
# print(run4)
# system(run4)# 查看运行状态
# 
# run5=glue::glue("curl https://swissdock.ch:8443/startdock\\?sessionNumber\\={session_number} ")
# print(run5)
# # 执行命令并捕获输出
# output = system(run5, intern = TRUE) # 开始运行dock
# 
# # 检查输出是否包含ERROR
# if (any(grepl("ERROR", output))) {
#   cat("检测到错误,停止脚本运行。\n")
#   cat("错误信息:\n")
#   cat(paste(output, collapse = "\n"))
#   stop("脚本执行终止")
# } else {
#   cat("Dock操作成功启动。\n")
#   cat("输出信息:\n")
#   cat(paste(output, collapse = "\n"))
# }
# 
# 
# # status=system(run4, intern = TRUE)# 运行命令并捕获输出
# # run6=glue::glue("curl https://www.swissdock.ch:8443/retrievesession?sessionNumber={session_number} -o results_{session_number}.zip")
# # system(run6)# 运行命令并捕获输出
# dir.create("output", showWarnings = FALSE)
# repeat {
#   # 检查状态
#   status <- system(run4, intern = TRUE)
#   print(status)
#   if (any(grepl("finished", status))) {
#     # 如果完成，下载结果
#     
#     run6 <- glue::glue("curl https://www.swissdock.ch:8443/retrievesession?sessionNumber={session_number} -o output/results_{Method}_{make.names(basename(ligand))}_{basename(target.path)}.zip")
#     print(run6)
#     system(run6) #获取结果
#     cat("Docking finished. Results downloaded.\n")
#     
#     # 读取下载结果生成表格
#     # 加载必要的库
#     library(utils)
#     library(tools)
#     
#     # 创建临时目录
#     temp_dir <- tempfile(pattern = "temp_unzip_")
#     dir.create(temp_dir)
#     
#     # 解压文件到临时目录
#     unzip(glue::glue("output/results_{Method}_{make.names(basename(ligand))}_{basename(target.path)}.zip"), exdir = temp_dir)
#     # unzip("output/results_Vina_X5hie_5hie_modified.pdb.zip")
#     # 读取pdbqt文件
#     file_path <- file.path(temp_dir, "vina_dock.pdb")
#     lines <- readLines(file_path)
#     
#     # 提取含有"REMARK VINA RESULT"的行
#     vina_results <- lines[grep("REMARK VINA RESULT", lines)]
#     vina_results_df=data.frame(vina_results = vina_results)
#     vina_results_df$Model=paste0("Model_",rownames(vina_results_df))
#     # 打印结果
#     print(vina_results_df)
#     Out.name=glue::glue("output/results_{Method}_{make.names(basename(ligand))}_{basename(target.path)}_score.csv")
#     vina_results_df$Path=Out.name
#     write.csv(vina_results_df,Out.name,row.names = FALSE)
#     # 删除临时目录及其内容
#     unlink(temp_dir, recursive = TRUE)
# 
#     break
#   } else {
#     # 如果未完成，等待 30 秒
#     cat("Waiting 30 seconds...\n")
#     Sys.sleep(30)
#   }
# }
# 
# 
# 
