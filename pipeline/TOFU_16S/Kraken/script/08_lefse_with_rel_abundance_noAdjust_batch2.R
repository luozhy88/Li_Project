library(microbiomeMarker)#1.6.0
library(phyloseq)
library(microbiome)
library(metagMisc)
library(openxlsx)
library(dplyr)
rm(list=ls())

# Create a directory named "lefse"
dir.create("../output/differential/lefse_barplot", recursive = TRUE)
path <- "../output/differential/lefse_barplot/"

dir.create("../output/differential/lefse_cladogram", recursive = TRUE)
path2 <- "../output/differential/lefse_cladogram/"

dir.create("../output/differential/lefse_abundance_barplot", recursive = TRUE)
path3 <- "../output/differential/lefse_abundance_barplot/"


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

## GROUP_VAR <- paste0("sex+age+", "group")
GROUP_VAR <- paste0("GROUP")

## here we refer to the main interest even though one has to add several covariates like above
MAIN <- trimws(GROUP_VAR, whitespace = ".*\\+")

############# split the phyloseq into a list based on the variable of interest
#######
agegut_split <- metagMisc::phyloseq_sep_variable(pseq, MAIN)
############ make pairwise comparison
pair <- combn(names(agegut_split), 2)


for (i in 1: (length(pair) / 2)   ) {
  #new phyloseq
  # i =1
  print(i)
  pseq.sub <- phyloseq::merge_phyloseq(agegut_split[[pair[1, i]]], agegut_split[[pair[2, i]]])
  pseq.temp <- core(pseq.sub, detection = 0.000001, prevalence = .1)

  ##
  for (my_rank in c( "Species", "Genus", "Family" ,"Class", "Order", "Phylum")) { try({
    #my_rank="Genus"
    print(my_rank)
    pseq.new <- phyloseq::tax_glom(pseq.temp, taxrank = my_rank)
    
    tax_table(pseq.new)[,my_rank] <- make.names(tax_table(pseq.new)[,my_rank])
    
    fil.names <- paste0(unique(microbiome::meta(pseq.new)[, MAIN])[1], "_vs_", unique(microbiome::meta(pseq.new)[, MAIN])[2])

    sum_res_lefse <- run_lefse(
      pseq.new,
      group = GROUP_VAR,
      subgroup = NULL,
      taxa_rank = my_rank,
      transform = c("identity"),
      norm = "CPM",
      norm_para = list(),
      kw_cutoff = 0.05,
      lda_cutoff = 2,
      bootstrap_n = 30,
      bootstrap_fraction = 2 / 3,
      wilcoxon_cutoff = 0.05,
      multigrp_strat = FALSE,
      strict = c("0"),
      sample_min = 10,
      only_same_subgrp = FALSE,
      curv = FALSE
    )


    sum_res <- data.frame(sum_res_lefse@marker_table)

    if (length(sum_res) > 0) {
      
    
    sum_res$group <- ifelse(sum_res$enrich_group == unique(microbiome::meta(pseq.new)[, MAIN])[1], "REF", "notREF")
    sum_res$ef_lda <- ifelse(sum_res$group == "notREF", sum_res$ef_lda, sum_res$ef_lda * (-1))

    rownames(sum_res) <- sum_res$feature
    ## we can take only the significantly differential taxa
    sig_lefse <- subset(sum_res, padj < 0.05)
    # Adding taxonomic labels
    taxa_info <- data.frame(tax_table(pseq.new))
    
    rownames(taxa_info) <- make.unique(taxa_info[, my_rank])

    sum_res <- merge(sum_res, taxa_info, by = 0)
    sig_lefse <- merge(sig_lefse, taxa_info, by = 0)

    REF <- sum_res %>%
      dplyr::filter(ef_lda < 0) %>%
      dplyr::pull(enrich_group) %>%
      unique()
    notREF <- sum_res %>%
      dplyr::filter(ef_lda > 0) %>%
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
    openxlsx::write.xlsx(sum_res, file = paste0(path, fil.names.new, my_rank, "__all_lefse_cal.xlsx"))
    

    ## reformat the ../output for barplot
    sig_lefse$P.adj <-
      ifelse(sig_lefse$padj <= 0.05,
        "< 0.05",
        ifelse(sig_lefse$padj <= 0.1 & sig_lefse$padj > 0.05,
          "0.05 - 0.1",
          "0.1 - 0.2 "
        )
      )
    ##

    ### do not need to plot if there were no significantly differential taxa
    if (nrow(sig_lefse) > 0) {
      ##
      sig_lefse$taxa <- sig_lefse$Row.names
      openxlsx::write.xlsx(sig_lefse, file = paste0(path, fil.names.new, my_rank, "__sig_lefse_cal.xlsx"))
      ## one can also plot the results based on log fold change LFC but colored by padj range
      ggpubr::ggbarplot(sig_lefse,
        x = "taxa", y = "ef_lda",
        fill = "enrich_group",
        color = "white",
        palette = c("#00A4D3", "#DD498D"),
        # palette = c("#DD498D", "#00A4D3"),
        sort.val = "asc",
        sort.by.sex = FALSE,
        x.text.angle = 90,
        ylab = "LDA",
        xlab = "",
        rotate = TRUE,
        title = paste0("[", notREF, "]", "-", "[", REF, "]"),
        ggtheme = theme_minimal()
      )

      ## save the plot
      ##
      ggsave(filename = paste0(path, fil.names.new, my_rank, "__lefse_cal_barplot.pdf"), width = 5, height = nrow(sig_lefse) * 0.1 + 2)


      p0 <- plot_abundance(sum_res_lefse, group = GROUP_VAR) + scale_fill_manual(values = c("#00A4D3", "#DD498D")) +
        scale_x_continuous(labels = scales::label_number_auto()) +
        scale_x_log10()
      ggsave(p0, filename = paste0(path3, fil.names.new, my_rank, "__lefse_abundance_barplot.pdf"), width = 6, height = nrow(sum_res_lefse@marker_table) * 0.15 + 2.5)
    }
  }
  ###
  ###

  sum_res_lefse_all <- run_lefse(
    pseq.temp,
    group = GROUP_VAR,
    subgroup = NULL,
    taxa_rank = "all",
    transform = c("identity"),
    norm = "CPM",
    norm_para = list(),
    kw_cutoff = 0.05,
    lda_cutoff = 2,
    bootstrap_n = 30,
    bootstrap_fraction = 2 / 3,
    wilcoxon_cutoff = 0.05,
    multigrp_strat = FALSE,
    strict = c("0"),
    sample_min = 10,
    only_same_subgrp = FALSE,
    curv = FALSE
  )

  if(length(unique(sum_res_lefse_all@marker_table[["enrich_group"]])) > 1){
    p1 <- plot_cladogram(sum_res_lefse_all, only_marker = TRUE, color = c("#00A4D3", "#DD498D")) +
      theme(plot.margin = margin(0, 0, 0, 0))
    ggsave(p1, filename = paste0(path2, fil.names.new, my_rank, "__lefse_cladogram.pdf"), width = 12, height = 12)
    
  } else{
    
    p1 <- plot_cladogram(sum_res_lefse_all, only_marker = TRUE, color = c("#d95f02")) +
      theme(plot.margin = margin(0, 0, 0, 0))
    ggsave(p1, filename = paste0(path2, fil.names.new, my_rank, "__lefse_cladogram.pdf"), width = 12, height = 12)
    
  }
  
  openxlsx::write.xlsx(sum_res_lefse_all@marker_table, file = paste0(path2, fil.names.new, my_rank, "__all_ranks_lefse_cal.xlsx"))

  p <- plot_abundance(sum_res_lefse_all, group = GROUP_VAR) + scale_fill_manual(values = c("#00A4D3", "#DD498D")) +
    scale_x_continuous(labels = scales::label_number_auto()) +
    scale_x_log10()

  ggsave(p, filename = paste0(path3, fil.names.new, "all_ranks", "__lefse_abundance_barplot.pdf"), width = 6, height = nrow(sum_res_lefse_all@marker_table) * 0.15 + 2.5)
  }
 )
  }
}
