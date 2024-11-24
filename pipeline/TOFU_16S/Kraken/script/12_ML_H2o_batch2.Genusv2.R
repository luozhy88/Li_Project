library(LinDA)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(openxlsx)
library(dplyr)
library(kableExtra)
rm(list=ls())

# Create a directory named "LinDA"
BATCH="SRP151288"
REF="HC"
notREF="TC"
tax="Genus"


path <- glue::glue("../output/H20/",BATCH,"_",paste(REF,notREF,sep  = "_vs_"),".",tax   )
dir.create(path, recursive = TRUE)




# Read the phyloseq object from a file and filter out taxa with zero counts
phy <- readRDS("../../Kraken/input/output/Batch2_phy_meta.rds")
# Meta=meta(phy)
phy



pseq = phyloseq::subset_samples(phy,grepl(BATCH,Batch )  &grepl( paste(REF,notREF,sep  = "|")  ,GROUP ))
pseq
Meta=meta(pseq)

#sample_data(pseq)$age <- scale(sample_data(pseq)$age)
###############################################################
############### empty values for ranks replaced with "undefined"##
###############################################################
# Extract the tax_table from the phyloseq object
taxa <- as.data.frame(tax_table(pseq))

# Replace empty or NA values with "undefined"
taxa[is.na(taxa)] <- "undefined"
taxa[taxa == ""] <- "undefined"

# Convert the modified data frame back to a tax_table format
tax_table_modified <- tax_table(as.matrix(taxa))

# Assign the modified tax_table back to the phyloseq object
tax_table(pseq) <- tax_table_modified
pseq.core <- core(pseq, detection = 0.000001, prevalence = .1)
pseq.core <- tax_glom(pseq.core, taxrank = "Genus")
pseq.core_Tax= phyloseq_to_df(pseq.core) %>% as.data.frame() 
rownames(pseq.core_Tax)=paste(pseq.core_Tax$OTU,pseq.core_Tax$Genus,sep = "_")
pseq.core_Tax=pseq.core_Tax[,-c(1:9)]%>% t() %>% as.data.frame()
Meta=meta(pseq.core) %>% dplyr::select("GROUP")

pseq.core_df= merge(Meta,pseq.core_Tax,by=0) %>% tibble::column_to_rownames(var="Row.names") %>% as.data.frame()


###############################example3#######################################

library(h2o)
h2o.init()

pseq.core_df <- pseq.core_df %>% as.h2o()
# Response column


# Split into train & test
splits <- h2o.splitFrame(pseq.core_df ,ratios = 0.8, seed = 200)
train <- splits[[1]]
test <- splits[[2]]


# Identify predictors and response
y <- "GROUP"
x <- setdiff(names(train), y)

# For binary classification, response should be a factor
train[, y] <- as.factor(train[, y])
test[, y] <- as.factor(test[, y])

##################################new methods##############################################


# 创建一个空的 dataframe 来存储结果
results <<- data.frame(
  model_id = character(),
  accuracy = numeric(),
  precision = numeric(),
  auc = numeric(),
  recall = numeric(),
  F1 = numeric(),
  stringsAsFactors = FALSE
)
ROC_PLOT_LIST<<-list()
Varimp.df_LIST<<-list()


get.re.ml=function(model.name="DRF",train,test,x,y){
  # model.name="XGBoost"
  print(model.name)
  aml <- h2o.automl(x = x, y = y,
                    training_frame = train,
                    max_models = 1,
                    # include_algos=c( "XGBoost"),
                    include_algos=c( model.name),
                    distribution = "AUTO",
                    seed = 2)
  
  # View the AutoML Leaderboard
  lb <- aml@leaderboard
  print(lb, n = nrow(lb))  # Print all rows instead of default (6 rows)
  lb_df=lb %>% as.data.frame()
  model_id=lb_df$model_id[1]
  lb_df$model_name=gsub("_.*","",model_id) 
  lb_df=lb_df%>% dplyr::distinct(model_name, .keep_all = TRUE)
  
  
  model <- h2o.getModel(lb_df$model_id[1])
  performance <- h2o.performance(model, newdata  = test)
  Accuracy <- h2o.accuracy(performance)[[2]] %>% mean()
  Precision=h2o.precision(performance)[[2]] %>% mean()
  Auc=h2o.auc(performance)
  Recall=h2o.recall(performance)[[2]] %>% mean()
  F1=h2o.F1(performance)[[2]] %>% mean()
  
  TPR=h2o.tpr(performance)[[2]]
  FPR=h2o.fpr(performance)[[2]]
  Auc.round=round(Auc,3)
  out.name=gsub("_.*","",model_id)
  
  # 生产ROC混合模型图的参数
  ROC_plot=data.frame(TPR=TPR,FPR=FPR )
  ROC_plot$Model=glue::glue(out.name,"(AUC={Auc.round})" )
  ROC_PLOT_LIST[[model.name]]<<-ROC_plot
  dir.create("output", showWarnings = FALSE)
  try({
    # ROC plot
    # pdf(glue::glue("output/ROC_{out.name}.pdf"),width = 6, height = 6)
    pdf( glue::glue(path,"/H20_ROC_test_plot_{out.name}.pdf"),width = 6, height = 6)
   
    plot(performance,main =out.name)
    legend("bottomright", legend = glue::glue("AUC={Auc.round}"), bty = "l")
    dev.off()
    
    # Varimp plot
    exa <- h2o.explain(model, test)
    Varimp=exa$varimp
    
    Varimp.df=h2o.varimp(model) %>% as.data.frame()
    Varimp.df$model=model.name
    Varimp.df_LIST[[model.name]]<<-Varimp.df
    
    Varimp_plot=Varimp$plots[[1]] +labs(title = out.name)
    
    width=Varimp_plot[["data"]][["variable"]] %>% unique() %>% nchar() %>% max()
    Width=0.15*width +3
    # ggsave(glue::glue("output/varimp_{out.name}.pdf"),Varimp_plot,width = Width, height = 4)
    ggsave(glue::glue(path,"/H20_featur_importance_plot_{out.name}.pdf"),Varimp_plot,width = Width, height = 4)
    
  },silent = T)
  
  # 将结果添加到 dataframe 中
  results <<- rbind(results, data.frame(
    Model_id = model_id,
    Accuracy = Accuracy,
    Precision = Precision,
    Auc = Auc,
    Recall = Recall,
    F1 = F1,
    Model.name = out.name,
    stringsAsFactors = FALSE
  ))
  print(results)
  # plot(performance)
  cat(sprintf("模型在测试集 %s 的准确率为: %f\n", model_id, Accuracy))
  
  
}



get.re.ml(model.name="DRF",train,test,x,y)
get.re.ml(model.name="XGBoost",train,test,x,y)
get.re.ml(model.name="GBM",train,test,x,y)
get.re.ml(model.name="GLM",train,test,x,y)
get.re.ml(model.name="DeepLearning",train,test,x,y)


#ROC plot
roc_data <- dplyr::bind_rows(ROC_PLOT_LIST)
ROC_plot_models=ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linetype = "dashed", size = 1) +
  # geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curves for Multiple Models", x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(legend.title = element_blank())
# ggsave(glue::glue("output/ROC_plot_models.pdf"),ROC_plot_models,width = 8, height = 6)
ggsave(glue::glue(path,"/H20_ROC_models_{tax}.pdf") ,ROC_plot_models,width = 8, height = 6)
results$Model_id=NULL
openxlsx::write.xlsx(results,glue::glue("output/ML_results.xlsx"), sheetName = "Sheet1", row.names = F)


#Model barplot
results_long <- reshape2::melt(results, id.vars = "Model.name")
colors <- c("#d7191c", "#fdae61", "#ffffbf", "#a6d96a", "#1a9641")
model_barplot=ggplot(results_long, aes(x = Model.name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Score", y = "Model", fill = "") +
  theme_minimal()+
  scale_fill_manual(values = colors) +
  theme(legend.position = "top")
# ggsave(glue::glue("output/models_barplot.pdf"),model_barplot,width = 6, height = 6)
ggsave(paste0(path,"/df_scores_model_barplot.pdf" ),model_barplot,width = 6, height = 6)


#Var.imp.heatmap
varimp_combined <- dplyr::bind_rows(Varimp.df_LIST)
varimp_combined$Importance=varimp_combined$scaled_importance 
varimp_combined.plot=ggplot(varimp_combined, aes(y = variable, x = model, fill = Importance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#a6d96a") +
  theme_minimal() +
  labs(title = "Feature Importance Heatmap",
       x = "",
       y = ""
  )+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
height=length(varimp_combined$variable )%>% unique() /30
# ggsave(glue::glue("output/varimp_combined.plot.pdf"),varimp_combined.plot,width = 6, height = height)
ggsave(glue::glue(path,"/H20_feature_importance_heatmap_plot.pdf"),varimp_combined.plot,width = 6, height = height)

varimp_wide <- varimp_combined %>% dplyr::select(-relative_importance ,-scaled_importance ,-percentage )%>%tidyr::pivot_wider(names_from = model,values_from = Importance) %>% dplyr::select(-DeepLearning) %>% dplyr::filter(DRF>0.00 & GBM>0.00 & GLM>0.00 )
write.csv(varimp_wide, file = glue::glue(path,"/H20_feature_importance.csv"), row.names = FALSE)
# 
# 
# ###################################old###############################################
# # Run AutoML for 20 base models
# aml <- h2o.automl(x = x, y = y,
#                   training_frame = train,
#                   max_models = 20,
#                   exclude_algos = c("XGBoost"),
#                   distribution = "bernoulli",
#                   seed = 1)
# 
# # View the AutoML Leaderboard
# lb <- aml@leaderboard
# df_score=as.data.frame(lb)
# write.csv(df_score, file = glue::glue(path,"/H20_model_scores.csv"), row.names = FALSE)
# 
# 
# 
# print(lb, n = nrow(lb))  # Print all rows instead of default (6 rows)
# h2o.varimp(lb)
# 
# pdf(file = glue::glue(path,"/H20_feature_importance_heatmap_plot.pdf"), width = 8, height = 6)
# h2o.varimp_heatmap(lb)
# dev.off()
# 
# 
# 
# # # Explain leader model & compare with all AutoML models
# # exa <- h2o.explain(aml, test)
# # exa$varimp
# # exa$varimp_heatmap
# # exa$model_correlation_heatmap
# 
# # Explain a single H2O model (e.g. leader model from AutoML)
# exm <- h2o.explain(aml@leader, train)
# # exm 
# 
# varimp_df <- as.data.frame(h2o.varimp(aml@leader))
# write.csv(varimp_df, file = glue::glue(path,"/H20_feature_importance.csv"), row.names = FALSE)
# 
# pdf(file = glue::glue(path,"/H20_featur_importance_plot.pdf"), width = 8, height = 4)
# exm$varimp 
# dev.off()
# 
# 
# # ROC
# model <- aml@leader
# # dl_perf <- h2o.performance(model,train = T)
# # dl_perf
# dl_perf <- h2o.performance(model,newdata = test)
# dl_perf
# 
# pdf(file = glue::glue(path,"/H20_ROC_test_plot.pdf"), width = 6, height = 6)
# plot(dl_perf)
# dev.off()
# 
# h2o.auc(dl_perf)
# h2o.auc(dl_perf, train = TRUE, valid = TRUE, xval = FALSE)
# 
# 
# # Extract model score from data.frame
# df_scroe_top = df_score %>% dplyr::select("model_id","auc","mean_per_class_error")
# table_q=kbl(df_scroe_top , align = "c") %>%
#   kable_classic(full_width = F)  %>% 
#   collapse_rows(columns = 1:2, valign = "top")
# 
# save_kable(table_q,paste0(path,"/df_scores_model.pdf" ))
# 
