library(LinDA)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(openxlsx)
library(dplyr)
rm(list=ls())

# Create a directory named "LinDA"
BATCH="SRP388727"
REF="HC"
notREF="TN"

path <- glue::glue("../output/H20/",BATCH,"_",paste(REF,notREF,sep  = "_vs_") )
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
splits <- h2o.splitFrame(pseq.core_df ,ratios = 0.8, seed = 100)
train <- splits[[1]]
test <- splits[[2]]


# Identify predictors and response
y <- "GROUP"
x <- setdiff(names(train), y)

# For binary classification, response should be a factor
train[, y] <- as.factor(train[, y])
test[, y] <- as.factor(test[, y])

# Run AutoML for 20 base models
aml <- h2o.automl(x = x, y = y,
                  training_frame = train,
                  max_models = 20,
                  exclude_algos = c("XGBoost"),
                  distribution = "bernoulli",
                  seed = 1)

# View the AutoML Leaderboard
lb <- aml@leaderboard
print(lb, n = nrow(lb))  # Print all rows instead of default (6 rows)
h2o.varimp(lb)

pdf(file = glue::glue(path,"/H20_feature_importance_heatmap_plot.pdf"), width = 8, height = 6)
h2o.varimp_heatmap(lb)
dev.off()



# # Explain leader model & compare with all AutoML models
# exa <- h2o.explain(aml, test)
# exa$varimp
# exa$varimp_heatmap
# exa$model_correlation_heatmap

# Explain a single H2O model (e.g. leader model from AutoML)
exm <- h2o.explain(aml@leader, train)
# exm 

varimp_df <- as.data.frame(h2o.varimp(aml@leader))
write.csv(varimp_df, file = glue::glue(path,"/H20_feature_importance.csv"), row.names = FALSE)

pdf(file = glue::glue(path,"/H20_featur_importance_plot.pdf"), width = 8, height = 4)
exm$varimp 
dev.off()


# ROC
model <- aml@leader
# dl_perf <- h2o.performance(model,train = T)
# dl_perf
dl_perf <- h2o.performance(model,newdata = test)
dl_perf

pdf(file = glue::glue(path,"/H20_ROC_test_plot.pdf"), width = 6, height = 6)
plot(dl_perf)
dev.off()

h2o.auc(dl_perf)
h2o.auc(dl_perf, train = TRUE, valid = TRUE, xval = FALSE)

