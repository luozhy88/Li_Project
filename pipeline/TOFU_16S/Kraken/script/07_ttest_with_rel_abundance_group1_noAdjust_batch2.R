library(microbiomeMarker)
library(phyloseq)
library(microbiome)
library(metagMisc)
library(openxlsx)
rm(list=ls())

# Create a directory named "Mann Whitney T test"
dir.create("../output/differential/ttest_barplot", recursive = TRUE)
path <- "../output/differential/ttest_barplot/"

dir.create("../output/differential/ttest_abundance_barplot", recursive = TRUE)
path2 <- "../output/differential/ttest_abundance_barplot/"



# Read the phyloseq object from a file and filter out taxa with zero counts
phy <- readRDS("../../Kraken/input/output/Batch2_phy_meta.rds")
# Meta=meta(phy)
phy

BATCH="SRP388727"
pseq = phyloseq::subset_samples(phy,grepl(BATCH,Batch )  &grepl("HC|TC",GROUP ))
pseq
Meta=meta(pseq)

tax_table(pseq) <- tax_table(pseq)[,c(1:7)]

phyloseq::tax_table(pseq)[, "Species"] <-sprintf("%s%04d", "ASV", seq_along(taxa_names(pseq)))

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
###############################################################
###############################################################



###### group1 ---> group

# names(sample_data(pseq))[7] <- "group"

pseq <- phyloseq::prune_taxa(taxa_sums(pseq) > 0, pseq)
pseq <- phyloseq::prune_samples(sample_sums(pseq) > 0, pseq)
sample_data(pseq)[, length(sample_data(pseq)) + 1] <- sample_data(pseq)[, 1]

## GROUP_VAR <- paste0("sex+age+", "sex")
GROUP_VAR <- paste0("GROUP")

## here we refer to the main interest even though one has to add several covariates like above
MAIN <- trimws(GROUP_VAR, whitespace = ".*\\+")

### here I rename all names (unreliable at all) with ASVs. But we must know they are ASVs not Species
# ## better to make sure there are no improper names (empty) at Genus level
# taxa_names(pseq) <- phyloseq::tax_table(pseq)[, "Species"]

############# split the phyloseq into a list based on the variable of interest
#######
agegut_split <- metagMisc::phyloseq_sep_variable(pseq, MAIN)
############ make pairwise comparison
pair <- combn(names(agegut_split), 2)


for (i in 1:(length(pair) / 2)) {
  ## new phyloseq
  #i=1
  print(i)
  pseq.sub <- phyloseq::merge_phyloseq(agegut_split[[pair[1, i]]], agegut_split[[pair[2, i]]])


  pseq.temp <- core(pseq.sub, detection = 0.000001, prevalence = .1)

  ##
  for (my_rank in c("Species", "Genus", "Family", "Class", "Order", "Phylum")) {
    #my_rank="Genus"
    print(my_rank)
    pseq.new <- phyloseq::tax_glom(pseq.temp, taxrank = my_rank)

    fil.names <- paste0(unique(microbiome::meta(pseq.new)[, MAIN])[1], "_vs_", unique(microbiome::meta(pseq.new)[, MAIN])[2])

    sum_res_ttest <- run_test_two_groups(
      pseq.new,
      group = GROUP_VAR,
      taxa_rank = my_rank,
      transform = c("identity"),
      norm = "CPM",
      norm_para = list(),
      method = c("white.test"),
      p_adjust = c("none"),
      pvalue_cutoff = 0.05,
      diff_mean_cutoff = NULL,
      ratio_cutoff = NULL,
      conf_level = 0.95,
      nperm = 1000
    )


    sum_res <- data.frame(sum_res_ttest@marker_table)
    
    
    if (length(sum_res) > 0) {

    sum_res$group <- ifelse(sum_res$ef_diff_mean < 0, "REF", "notREF")

    rownames(sum_res) <- sum_res$feature
    ## we can take only the significantly differential taxa
    sig_ttest <- subset(sum_res, padj < 0.05)
    # Adding taxonomic labels
    taxa_info <- data.frame(tax_table(pseq.new))
    rownames(taxa_info) <- make.unique(taxa_info[, my_rank])
    
    sum_res <- merge(sum_res, taxa_info, by = 0)
    sig_ttest <- merge(sig_ttest, taxa_info, by = 0)

    REF <- sum_res %>%
      dplyr::filter(ef_diff_mean < 0) %>%
      dplyr::pull(enrich_group) %>%
      unique()
    notREF <- sum_res %>%
      dplyr::filter(ef_diff_mean > 0) %>%
      dplyr::pull(enrich_group) %>%
      unique()

    
    
    if(length(REF) == 0){
      REF <- setdiff(unique(as.data.frame(sample_data(pseq.new))[,"group"])$group, notREF)
    }
    
    if(length(notREF) == 0){
      notREF <- setdiff(unique(as.data.frame(sample_data(pseq.new))[,"group"])$group, REF)  
    }
    
    
    ##
    fil.names.new <- paste0(fil.names, "_***", REF, "_as_reference***",BATCH,"_")
    ## save the ../output for each comparison, reference shown on the file name
    openxlsx::write.xlsx(sum_res, file = paste0(path, fil.names.new, my_rank, "__all_ttest_cal.xlsx"))
    openxlsx::write.xlsx(sig_ttest, file = paste0(path, fil.names.new, my_rank, "__sig_ttest_cal.xlsx"))

    ## reformat the ../output for barplot
    sig_ttest$P.adj <-
      ifelse(sig_ttest$padj <= 0.05,
        "< 0.05",
        ifelse(sig_ttest$padj <= 0.1 & sig_ttest$padj > 0.05,
          "0.05 - 0.1",
          "0.1 - 0.2 "
        )
      )
    ##

    ### do not need to plot if there were no significantly differential taxa
    if (nrow(sig_ttest) > 0) {
      ##
      openxlsx::write.xlsx(sig_ttest, file = paste0(path, fil.names.new, my_rank, "__sig_ttest_cal.xlsx"))
      sig_ttest$taxa <- sig_ttest$Row.names
      ## one can also plot the results based on log fold change LFC but colored by padj range
      ggpubr::ggbarplot(sig_ttest,
        x = "taxa", y = "ef_diff_mean",
        fill = "enrich_group",
        color = "white",
        # palette = c("#DD498D", "#00A4D3"),
         palette = c("#00A4D3","#DD498D"),
        sort.val = "desc",
        sort.by.sex = FALSE,
        x.text.angle = 90,
        ylab = "effect_diff_mean",
        xlab = "",
        rotate = TRUE,
        title = paste0("[", notREF, "]", "-", "[", REF, "]"),
        ggtheme = theme_minimal()
      )

      ## save the plot
      ##
      ggsave(filename = paste0(path, fil.names.new, my_rank, "__ttest_cal_barplot.pdf"), width = 6, height = nrow(sig_ttest) * 0.1 + 2)
   
      p0 <- plot_abundance(sum_res_ttest, group = GROUP_VAR) + scale_fill_manual(values = c("#00A4D3","#DD498D") ) +
        scale_x_log10()
      ggsave(p0, filename = paste0(path2, fil.names.new, my_rank, "__ttest_abundance_barplot.pdf"), width = 6, height = nrow(sum_res_ttest@marker_table) * 0.15 + 2.5)
      
       }
    }
    
  }
}
