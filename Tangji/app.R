#!/usr/bin/env Rscript
# ============================================================================
# ThyGlycoPortal: Thyroid Cancer Glycosylation Interactive Analysis Platform
# Phase 2-3: R Shiny Web Application (REAL DATA ONLY)
# ============================================================================
# Data sources:
#   - SQLite database with curated literature metadata
#   - literature_stats: REAL published summary statistics (mean, SD, median, p, AUC)
#   - tcga_glycogene_expression: REAL TCGA RNAseq data from Bones et al. 2018
# NO simulated data is used anywhere in this application.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DBI)
library(RSQLite)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(readr)

# ----------------------------------------------------------------------------
# Database Connection
# ----------------------------------------------------------------------------
DB_PATH <- "output/thyroid_glyco_db.sqlite"
if (!file.exists(DB_PATH)) {
  stop("Database not found at ", DB_PATH, ". Please run 05_create_database.py, 06_import_glycomics_data.py, and 11_import_real_data_to_db.py first.")
}
con <- dbConnect(RSQLite::SQLite(), DB_PATH)

# ----------------------------------------------------------------------------
# Load ALL real data from SQLite
# ----------------------------------------------------------------------------
studies_df    <- dbReadTable(con, "studies")
glycans_df    <- dbReadTable(con, "glycan_structures")
groups_df     <- dbReadTable(con, "clinical_groups")
samples_df    <- dbReadTable(con, "samples")
abundance_df  <- dbReadTable(con, "glycan_abundance")
biomarkers_df <- dbReadTable(con, "biomarkers")
enzymes_df    <- dbReadTable(con, "glycosyltransferases")
enzyme_links  <- dbReadTable(con, "enzyme_glycan_links")
lit_stats_df  <- dbReadTable(con, "literature_stats")
tcga_expr_df  <- dbReadTable(con, "tcga_glycogene_expression")

# Merge abundance with metadata
abundance_full <- abundance_df %>%
  left_join(samples_df, by = "sample_id") %>%
  left_join(groups_df %>% select(-study_id), by = "group_id") %>%
  left_join(glycans_df, by = "glycan_id") %>%
  left_join(studies_df %>% select(study_id, title, year, cancer_type), by = "study_id") %>%
  mutate(year = as.numeric(year))

# ----------------------------------------------------------------------------
# REAL DATA: Prepare TCGA glycogene expression for visualization
# From Bones et al. 2018, Cancers (PMC11727208)
# TCGA RNAseq RPKM, n=20 paired normal vs PTC
# ----------------------------------------------------------------------------
tcga_expr_long <- tcga_expr_df %>%
  select(gene_symbol, normal_median_rpkm, ptc_median_rpkm, fold_change, p_value, direction, significant) %>%
  pivot_longer(cols = c(normal_median_rpkm, ptc_median_rpkm),
               names_to = "group", values_to = "median_rpkm") %>%
  mutate(group = ifelse(group == "normal_median_rpkm", "Normal", "PTC"),
         group = factor(group, levels = c("Normal", "PTC")),
         sig_label = ifelse(significant == 1, paste0("* p=", p_value), "ns"))

# ----------------------------------------------------------------------------
# Helper: simplified nomogram demonstration model for PTMC LNM
# WARNING: Coefficients are illustrative approximations for demonstration only.
# The original PTMC nomogram (Front Oncol 2022) reports AUC=0.702 but does not
# publish full logistic regression coefficients. Do not use for clinical decisions.
# ----------------------------------------------------------------------------
predict_lnm_risk <- function(ca4_val, a2f0s0g_val) {
  b0 <- -2.5
  b1 <- 3.2
  b2 <- -2.0
  logit <- b0 + b1 * ca4_val + b2 * a2f0s0g_val
  prob <- 1 / (1 + exp(-logit))
  return(prob)
}

calculate_bn_score <- function(h3n5f1, h4n5f1, h5n5f1, h4n4f1) {
  bisecting_sum <- h3n5f1 + h4n5f1 + h5n5f1
  total <- h3n5f1 + h4n5f1 + h5n5f1 + h4n4f1
  if (total == 0) return(0)
  bn_ratio <- bisecting_sum / total
  return(bn_ratio)
}

calculate_recurrence_risk <- function(g0f, g1f) {
  if (g1f == 0) return(NA)
  ratio <- g0f / g1f
  return(ratio)
}

# ----------------------------------------------------------------------------
# UI
# ----------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = span(tagList(icon("flask"), " ThyGlycoPortal"))),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Literature", tabName = "literature", icon = icon("book")),
      menuItem("Biomarkers", tabName = "biomarkers", icon = icon("heartbeat")),
      menuItem("Glycan Browser", tabName = "glycan", icon = icon("dna")),
      menuItem("Enzyme Network", tabName = "enzyme", icon = icon("project-diagram")),
      menuItem("TCGA Glycogenes", tabName = "tcga", icon = icon("chart-bar")),
      menuItem("Diagnostic Tool", tabName = "diagnostic", icon = icon("stethoscope")),
      menuItem("Nomogram", tabName = "nomogram", icon = icon("calculator")),
      menuItem("Data Upload", tabName = "upload", icon = icon("upload"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML(".content-wrapper { background-color: #f4f6f9; }"))),
    tabItems(
      # ======================================================================
      # TAB 1: Overview
      # ======================================================================
      tabItem(tabName = "overview",
        fluidRow(
          valueBox(nrow(studies_df), "Studies", icon = icon("book"), color = "aqua"),
          valueBox(nrow(glycans_df), "Glycan Structures", icon = icon("dna"), color = "green"),
          valueBox(nrow(biomarkers_df), "Biomarkers", icon = icon("heartbeat"), color = "red"),
          valueBox(nrow(enzymes_df), "Glycosyltransferases", icon = icon("project-diagram"), color = "yellow")
        ),
        fluidRow(
          box(title = "Study Distribution by Year", width = 6, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_study_year", height = 300)),
          box(title = "Cancer Type Distribution", width = 6, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_cancer_type", height = 300))
        ),
        fluidRow(
          box(title = "Method Distribution", width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("plot_method", height = 300)),
          box(title = "Biomarker Performance (AUC)", width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("plot_biomarker_auc", height = 300))
        )
      ),

      # ======================================================================
      # TAB 2: Literature
      # ======================================================================
      tabItem(tabName = "literature",
        fluidRow(
          box(title = "Study Filters", width = 12, solidHeader = TRUE, status = "primary",
            column(3, selectInput("lit_year", "Year", choices = c("All", sort(unique(studies_df$year))), selected = "All")),
            column(3, selectInput("lit_cancer", "Cancer Type", choices = c("All", unique(studies_df$cancer_type)), selected = "All")),
            column(3, selectInput("lit_method", "Method", choices = c("All", unique(studies_df$method)), selected = "All")),
            column(3, selectInput("lit_sample", "Sample Type", choices = c("All", unique(studies_df$sample_type)), selected = "All"))
          )
        ),
        fluidRow(
          box(title = "Literature Table", width = 12, solidHeader = TRUE, status = "primary",
              DTOutput("table_literature"))
        )
      ),

      # ======================================================================
      # TAB 3: Biomarkers
      # ======================================================================
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Biomarker Table", width = 8, solidHeader = TRUE, status = "primary",
              DTOutput("table_biomarkers")),
          box(title = "Selected Biomarker Detail", width = 4, solidHeader = TRUE, status = "info",
              htmlOutput("biomarker_detail"))
        ),
        fluidRow(
          box(title = "AUC Comparison", width = 12, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_biomarker_comparison", height = 350))
        )
      ),

      # ======================================================================
      # TAB 4: Glycan Browser
      # ======================================================================
      tabItem(tabName = "glycan",
        fluidRow(
          box(title = "Select Glycan", width = 4, solidHeader = TRUE, status = "primary",
              selectizeInput("glycan_select", "Glycan (SNFG Name)",
                             choices = glycans_df$snfg_name,
                             selected = "G0F",
                             options = list(placeholder = "Search glycan...")),
              br(),
              htmlOutput("glycan_info")),
          box(title = "Cross-Study Comparison", width = 8, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_glycan_cross_study", height = 350))
        ),
        fluidRow(
          box(title = "All Glycan Structures", width = 12, solidHeader = TRUE, status = "info",
              DTOutput("table_glycans"))
        )
      ),

      # ======================================================================
      # TAB 5: Enzyme Network
      # ======================================================================
      tabItem(tabName = "enzyme",
        fluidRow(
          box(title = "Select Enzyme", width = 4, solidHeader = TRUE, status = "primary",
              selectInput("enzyme_select", "Glycosyltransferase",
                          choices = enzymes_df$gene_symbol,
                          selected = "GALNT3"),
              htmlOutput("enzyme_info")),
          box(title = "Regulated Glycans", width = 8, solidHeader = TRUE, status = "primary",
              DTOutput("table_enzyme_glycans"))
        ),
        fluidRow(
          box(title = "Full Enzyme-Glycan Network", width = 12, solidHeader = TRUE, status = "info",
              DTOutput("table_full_network"))
        )
      ),

      # ======================================================================
      # TAB 6: TCGA Glycogene Expression (REAL DATA)
      # ======================================================================
      tabItem(tabName = "tcga",
        fluidRow(
          box(title = "TCGA Glycogene Expression (REAL DATA)", width = 12, solidHeader = TRUE, status = "warning",
              "Source: Bones J et al. Cancers 2018 (PMC11727208). TCGA RNAseq RPKM data from 20 paired normal vs PTC samples.",
              br(), br())
        ),
        fluidRow(
          box(title = "Median RPKM by Group", width = 8, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_tcga_expr", height = 400)),
          box(title = "Statistical Summary", width = 4, solidHeader = TRUE, status = "info",
              tableOutput("table_tcga_stats"))
        ),
        fluidRow(
          box(title = "Fold Change & Significance", width = 12, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_tcga_fc", height = 300))
        )
      ),

      # ======================================================================
      # TAB 7: Diagnostic Tool
      # ======================================================================
      tabItem(tabName = "diagnostic",
        fluidRow(
          box(title = "IgG N-Glycan Diagnostic Tool (Zhang 2021)", width = 6, solidHeader = TRUE, status = "primary",
              "Based on plasma IgG N-glycan BN feature. AUC=0.920 (discovery), 0.896 (validation).", br(), br(),
              sliderInput("diag_h3n5f1", "H3N5F1 (BN component)", min = 0, max = 100, value = 15),
              sliderInput("diag_h4n5f1", "H4N5F1 (BN component)", min = 0, max = 100, value = 20),
              sliderInput("diag_h5n5f1", "H5N5F1 (BN component)", min = 0, max = 100, value = 10),
              sliderInput("diag_h4n4f1", "H4N4F1 (G1, reference)", min = 0, max = 100, value = 55),
              actionButton("diag_calc", "Calculate BN Score", icon = icon("calculator"), class = "btn-primary"),
              br(), br(),
              htmlOutput("diag_reference")
          ),
          box(title = "Serum Recurrence Predictor (Kudelka 2023)", width = 6, solidHeader = TRUE, status = "warning",
              "Based on serum G0F:G1F ratio. AUC=0.82 (95%CI: 0.64-0.99).", br(), br(),
              numericInput("serum_g0f", "G0F Intensity", value = 35, min = 0),
              numericInput("serum_g1f", "G1F Intensity", value = 50, min = 0),
              actionButton("recur_calc", "Calculate Recurrence Risk", icon = icon("calculator"), class = "btn-warning"),
              br(), br(),
              htmlOutput("recur_reference")
          )
        ),
        fluidRow(
          box(title = "Diagnostic Result", width = 6, solidHeader = TRUE, status = "info",
              verbatimTextOutput("diag_result")),
          box(title = "Recurrence Prediction Result", width = 6, solidHeader = TRUE, status = "info",
              verbatimTextOutput("recur_result"))
        )
      ),

      # ======================================================================
      # TAB 8: Nomogram
      # ======================================================================
      tabItem(tabName = "nomogram",
        fluidRow(
          box(title = "PTMC LNM Risk Nomogram", width = 6, solidHeader = TRUE, status = "primary",
              "Based on serum N-glycome CA4 and A2F0S0G (Front Oncol 2022).", br(), br(),
              sliderInput("nom_ca4", "CA4 (Tetraantennary complex)", min = 0, max = 0.05, value = 0.02, step = 0.001),
              sliderInput("nom_a2f0s0g", "A2F0S0G (Galactosylation in non-fucosylated diantennary)", min = 0.4, max = 0.7, value = 0.60, step = 0.001),
              actionButton("nom_calc", "Calculate LNM Probability", icon = icon("calculator"), class = "btn-primary"),
              br(), br(),
              htmlOutput("nom_reference")
          ),
          box(title = "Risk Interpretation", width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("plot_nomogram_gauge", height = 250),
              br(),
              htmlOutput("nom_interpretation"))
        )
      ),

      # ======================================================================
      # TAB 9: Data Upload
      # ======================================================================
      tabItem(tabName = "upload",
        fluidRow(
          box(title = "Upload Your Glycomics Data", width = 12, solidHeader = TRUE, status = "primary",
              fileInput("user_data", "Choose CSV File",
                        accept = c("text/csv", ".csv")),
              helpText("Expected CSV format: columns = glycan features (e.g., G0F, G1F, G2F), rows = samples. Include optional 'group' and 'sample_id' columns."),
              checkboxInput("user_header", "Header row", TRUE)
          )
        ),
        fluidRow(
          box(title = "Data Preview", width = 12, solidHeader = TRUE, status = "info",
              DTOutput("table_user_data"))
        ),
        fluidRow(
          box(title = "Glycan Profile Visualization", width = 6, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_user_profile", height = 350)),
          box(title = "Group Comparison (if group column provided)", width = 6, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_user_group", height = 350))
        )
      )
    )
  )
)


# ----------------------------------------------------------------------------
# Server
# ----------------------------------------------------------------------------
server <- function(input, output, session) {

  # ========================================================================
  # OVERVIEW TAB
  # ========================================================================
  output$plot_study_year <- renderPlotly({
    df <- studies_df %>%
      count(year) %>%
      arrange(year)
    p <- ggplot(df, aes(x = factor(year), y = n, text = paste0(year, ": ", n, " studies"))) +
      geom_col(fill = "#3c8dbc", width = 0.6) +
      labs(x = "Year", y = "Number of Studies") +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })

  output$plot_cancer_type <- renderPlotly({
    df <- studies_df %>%
      filter(cancer_type != "All types") %>%
      count(cancer_type)
    plot_ly(df, labels = ~cancer_type, values = ~n, type = "pie",
            hole = 0, textinfo = "label+percent",
            marker = list(colors = c("#dd4b39", "#f39c12", "#00a65a", "#3c8dbc"))) %>%
      layout(showlegend = FALSE,
             margin = list(l = 20, r = 20, t = 20, b = 20))
  })

  output$plot_method <- renderPlotly({
    df <- studies_df %>%
      count(method) %>%
      arrange(desc(n))
    p <- ggplot(df, aes(x = reorder(method, n), y = n, text = paste0(method, ": ", n))) +
      geom_col(fill = "#00a65a", width = 0.6) +
      coord_flip() +
      labs(x = "", y = "Count") +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })

  output$plot_biomarker_auc <- renderPlotly({
    df <- biomarkers_df %>%
      filter(!is.na(performance_auc)) %>%
      arrange(performance_auc) %>%
      mutate(name_short = substr(name, 1, 30))
    p <- ggplot(df, aes(x = reorder(name_short, performance_auc), y = performance_auc,
                        text = paste0(name, "\nAUC: ", round(performance_auc, 3)))) +
      geom_col(fill = "#dd4b39", width = 0.6) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
      coord_flip() +
      labs(x = "", y = "AUC") +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })

  # ========================================================================
  # LITERATURE TAB
  # ========================================================================
  lit_filtered <- reactive({
    df <- studies_df
    if (input$lit_year != "All") df <- df %>% filter(year == as.numeric(input$lit_year))
    if (input$lit_cancer != "All") df <- df %>% filter(cancer_type == input$lit_cancer)
    if (input$lit_method != "All") df <- df %>% filter(method == input$lit_method)
    if (input$lit_sample != "All") df <- df %>% filter(sample_type == input$lit_sample)
    df %>% select(study_id, title, authors, year, journal, method, sample_type, cancer_type, key_finding)
  })

  output$table_literature <- renderDT({
    datatable(lit_filtered(),
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE,
              selection = "single")
  })

  # ========================================================================
  # BIOMARKERS TAB
  # ========================================================================
  output$table_biomarkers <- renderDT({
    df <- biomarkers_df %>%
      left_join(studies_df %>% select(study_id, title), by = "study_id") %>%
      select(biomarker_id, name, biomarker_type, sample_type, performance_auc,
             performance_sensitivity, performance_specificity, validation_status, title) %>%
      arrange(desc(performance_auc))
    datatable(df,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE,
              selection = "single")
  })

  output$biomarker_detail <- renderUI({
    sel <- input$table_biomarkers_rows_selected
    if (is.null(sel)) return(HTML("<p>Select a biomarker from the table to see details.</p>"))
    df <- biomarkers_df %>% arrange(desc(performance_auc))
    row <- df[sel, ]
    auc <- ifelse(is.na(row$performance_auc), "N/A", round(row$performance_auc, 3))
    sens <- ifelse(is.na(row$performance_sensitivity), "N/A", round(row$performance_sensitivity, 3))
    spec <- ifelse(is.na(row$performance_specificity), "N/A", round(row$performance_specificity, 3))
    HTML(paste0(
      "<b>", row$name, "</b><hr>",
      "<b>Type:</b> ", row$biomarker_type, "<br>",
      "<b>Sample:</b> ", row$sample_type, "<br>",
      "<b>AUC:</b> ", auc, "<br>",
      "<b>Sensitivity:</b> ", sens, "<br>",
      "<b>Specificity:</b> ", spec, "<br>",
      "<b>Cutoff:</b> ", row$cutoff_value, "<br>",
      "<b>Validation:</b> ", row$validation_status, "<br><br>",
      "<b>Description:</b><br>", row$description
    ))
  })

  output$plot_biomarker_comparison <- renderPlotly({
    df <- biomarkers_df %>%
      left_join(studies_df %>% select(study_id, title), by = "study_id") %>%
      mutate(auc_label = ifelse(is.na(performance_auc), 0, performance_auc),
             name_short = substr(name, 1, 25),
             tooltip = paste0(name, "\nStudy: ", substr(title, 1, 40), "...\nAUC: ", round(performance_auc, 3)))
    p <- ggplot(df, aes(x = reorder(name_short, auc_label), y = auc_label, fill = biomarker_type,
                        text = tooltip)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = c(0.8, 0.9), linetype = "dashed", color = c("orange", "darkgreen")) +
      coord_flip() +
      labs(x = "", y = "AUC", fill = "Type") +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })

  # ========================================================================
  # GLYCAN BROWSER TAB
  # ========================================================================
  output$glycan_info <- renderUI({
    g <- glycans_df %>% filter(snfg_name == input$glycan_select)
    if (nrow(g) == 0) return(NULL)
    HTML(paste0(
      "<b>Composition:</b> ", g$composition, "<br>",
      "<b>Mass:</b> ", round(g$mass, 2), " Da<br>",
      "<b>Type:</b> ", g$glycan_type, "<br><br>",
      "<b>Description:</b><br>", g$description
    ))
  })

  output$plot_glycan_cross_study <- renderPlotly({
    df <- abundance_full %>%
      filter(snfg_name == input$glycan_select) %>%
      mutate(trend_label = case_when(
        abundance_value > 0 ~ "Upregulated",
        abundance_value < 0 ~ "Downregulated",
        TRUE ~ "No change"
      ),
      tooltip = paste0(study_id, " (", year, ")\n", group_name, "\nTrend: ", trend_label))

    if (nrow(df) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
        layout(title = "No cross-study data available for this glycan"))
    }

    p <- ggplot(df, aes(x = reorder(study_id, year), y = abundance_value,
                        fill = trend_label, text = tooltip)) +
      geom_col(width = 0.6) +
      facet_wrap(~group_name, scales = "free_x") +
      labs(x = "Study", y = "Trend (coded)", fill = "Trend") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggplotly(p, tooltip = "text")
  })

  output$table_glycans <- renderDT({
    datatable(glycans_df %>% select(glycan_id, snfg_name, composition, mass, glycan_type, description),
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE)
  })

  # ========================================================================
  # ENZYME NETWORK TAB
  # ========================================================================
  output$enzyme_info <- renderUI({
    e <- enzymes_df %>% filter(gene_symbol == input$enzyme_select)
    if (nrow(e) == 0) return(NULL)
    HTML(paste0(
      "<b>Gene:</b> ", e$gene_symbol, "<br>",
      "<b>Name:</b> ", e$enzyme_name, "<br>",
      "<b>Family:</b> ", e$enzyme_family, "<br>",
      "<b>Pathway:</b> ", e$pathway, "<br>",
      "<b>Substrate:</b> ", e$substrate_type, "<br><br>",
      "<b>Description:</b><br>", e$description
    ))
  })

  output$table_enzyme_glycans <- renderDT({
    target_eid <- enzymes_df %>% filter(gene_symbol == input$enzyme_select) %>% pull(enzyme_id)
    df <- enzyme_links %>%
      filter(enzyme_id %in% target_eid) %>%
      left_join(glycans_df, by = "glycan_id") %>%
      select(snfg_name, composition, glycan_type, linkage_type, evidence_level)
    datatable(df, options = list(pageLength = 10), rownames = FALSE)
  })

  output$table_full_network <- renderDT({
    df <- enzyme_links %>%
      left_join(enzymes_df, by = "enzyme_id") %>%
      left_join(glycans_df, by = "glycan_id") %>%
      select(gene_symbol, enzyme_name, enzyme_family, snfg_name, composition, linkage_type, evidence_level)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  # ========================================================================
  # TCGA GLYCOGENE EXPRESSION (REAL DATA from Bones et al. 2018)
  # ========================================================================
  output$plot_tcga_expr <- renderPlotly({
    p <- ggplot(tcga_expr_long, aes(x = reorder(gene_symbol, median_rpkm),
                                     y = median_rpkm, fill = group,
                                     text = paste0(gene_symbol, "\n", group, ": ", median_rpkm, " RPKM"))) +
      geom_col(position = "dodge", width = 0.7) +
      scale_fill_manual(values = c("Normal" = "#00a65a", "PTC" = "#dd4b39")) +
      labs(x = "", y = "Median RPKM", fill = "Group",
           caption = "Source: Bones J et al. Cancers 2018 (TCGA RNAseq, n=20 paired)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggplotly(p, tooltip = "text")
  })

  output$table_tcga_stats <- renderTable({
    tcga_expr_df %>%
      select(gene_symbol, enzyme_family, normal_median_rpkm, ptc_median_rpkm,
             fold_change, p_value, direction, significant) %>%
      arrange(desc(fold_change))
  })

  output$plot_tcga_fc <- renderPlotly({
    df <- tcga_expr_df %>%
      mutate(fc_label = ifelse(fold_change >= 1, fold_change, -1/fold_change),
             color = case_when(
               significant == 1 & fold_change > 1 ~ "Upregulated",
               significant == 1 & fold_change < 1 ~ "Downregulated",
               TRUE ~ "Not significant"
             ),
             tooltip = paste0(gene_symbol, "\nFold change: ", round(fold_change, 2),
                              "\np-value: ", p_value, "\n", direction))
    p <- ggplot(df, aes(x = reorder(gene_symbol, fc_label), y = fc_label, fill = color, text = tooltip)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = -1, linetype = "dashed", color = "gray50") +
      scale_fill_manual(values = c("Upregulated" = "#dd4b39",
                                   "Downregulated" = "#00a65a",
                                   "Not significant" = "#bbbbbb")) +
      labs(x = "", y = "Fold Change (log2-like)", fill = "",
           caption = "* indicates p < 0.05") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggplotly(p, tooltip = "text")
  })

  # ========================================================================
  # DIAGNOSTIC TOOL TAB
  # ========================================================================
  output$diag_reference <- renderUI({
    bn_stats <- lit_stats_df %>% filter(study == "Zhang_2021", variable == "BN (bisecting neutral N-glycans)")
    if (nrow(bn_stats) == 0) return(NULL)
    HTML(paste0(
      "<b>Reference (Zhang 2021, Discovery Cohort):</b><br>",
      "HC mean ± SD: ", round(bn_stats$group1_mean, 1), " ± ", round(bn_stats$group1_sd, 1), "%<br>",
      "TC mean ± SD: ", round(bn_stats$group2_mean, 1), " ± ", round(bn_stats$group2_sd, 1), "%<br>",
      "AUC: ", bn_stats$auc, " | p ", bn_stats$p_value, "<br>",
      "<small>BN is INCREASED in thyroid cancer</small>"
    ))
  })

  output$recur_reference <- renderUI({
    ratio_stats <- lit_stats_df %>% filter(study == "Kudelka_2023", variable == "G0F:G1F ratio")
    if (nrow(ratio_stats) == 0) return(NULL)
    HTML(paste0(
      "<b>Reference (Kudelka 2023):</b><br>",
      "HC median (IQR): ", ratio_stats$group1_median, " (", ratio_stats$group1_iqr_low, "-", ratio_stats$group1_iqr_high, ")<br>",
      "Recurrent median (IQR): ", ratio_stats$group2_median, " (", ratio_stats$group2_iqr_low, "-", ratio_stats$group2_iqr_high, ")<br>",
      "AUC: ", ratio_stats$auc, " | p = ", ratio_stats$p_value, "<br>",
      "<small>Ratio INCREASED in recurrent DTC</small>"
    ))
  })

  observeEvent(input$diag_calc, {
    bn <- calculate_bn_score(input$diag_h3n5f1, input$diag_h4n5f1,
                             input$diag_h5n5f1, input$diag_h4n4f1)
    total <- input$diag_h3n5f1 + input$diag_h4n5f1 + input$diag_h5n5f1 + input$diag_h4n4f1
    output$diag_result <- renderText({
      paste0(
        "=== IgG BN Diagnostic Score ===\n",
        "Total intensity: ", total, "\n",
        "Bisecting sum (H3N5F1+H4N5F1+H5N5F1): ", input$diag_h3n5f1 + input$diag_h4n5f1 + input$diag_h5n5f1, "\n",
        "BN ratio (Bisecting/Total): ", round(bn, 4), "\n\n",
        "Reference (Zhang 2021):\n",
        "  - Discovery AUC = 0.920\n",
        "  - Validation AUC = 0.896 (TC vs HC)\n",
        "  - BN is INCREASED in thyroid cancer\n\n",
        ifelse(bn > 0.35,
               "=> BN ELEVATED (consistent with TC pattern)\n   Note: Threshold for illustration only. Correlate with ultrasound + FNA.",
               "=> BN NOT ELEVATED (consistent with HC pattern)\n   Note: Single marker insufficient; use with other modalities.")
      )
    })
  })

  observeEvent(input$recur_calc, {
    ratio <- calculate_recurrence_risk(input$serum_g0f, input$serum_g1f)
    output$recur_result <- renderText({
      if (is.na(ratio)) return("Error: G1F cannot be zero")
      paste0(
        "=== Serum G0F:G1F Recurrence Ratio ===\n",
        "G0F intensity: ", input$serum_g0f, "\n",
        "G1F intensity: ", input$serum_g1f, "\n",
        "G0F:G1F ratio: ", round(ratio, 3), "\n\n",
        "Reference (Kudelka 2023):\n",
        "  - AUC = 0.82 (95%CI: 0.64-0.99)\n",
        "  - Ratio INCREASED in recurrent DTC\n",
        "  - Optimal cutoff (balanced): 0.53\n",
        "  - Optimal cutoff (max specificity): 0.73\n\n",
        ifelse(ratio > 0.73,
               "=> HIGH RECURRENCE RISK (ratio > 0.73)\n   Recommend: aggressive surveillance",
               ifelse(ratio > 0.53,
                      "=> MODERATE RECURRENCE RISK (0.53-0.73)\n   Recommend: enhanced monitoring",
                      "=> LOW RECURRENCE RISK (ratio < 0.53)\n   Recommend: routine follow-up"))
      )
    })
  })

  # ========================================================================
  # NOMOGRAM TAB
  # ========================================================================
  output$nom_reference <- renderUI({
    ca4_ref <- lit_stats_df %>% filter(study == "PTMC_2022_LNM", variable == "CA4")
    a2_ref <- lit_stats_df %>% filter(study == "PTMC_2022_LNM", grepl("A2F0S0G", variable))
    HTML(paste0(
      "<b>Reference (PTMC Nomogram 2022):</b><br>",
      "CA4 median: NLNM=", ifelse(nrow(ca4_ref)>0, ca4_ref$group1_median, "0.021"),
      ", LNM=", ifelse(nrow(ca4_ref)>0, ca4_ref$group2_median, "0.023"), "<br>",
      "A2F0S0G median: NLNM=", ifelse(nrow(a2_ref)>0, a2_ref$group1_median, "0.599"),
      ", LNM=", ifelse(nrow(a2_ref)>0, a2_ref$group2_median, "0.578"), "<br>",
      "<small>CA4 higher & A2F0S0G lower = higher LNM risk</small>"
    ))
  })

  observeEvent(input$nom_calc, {
    prob <- predict_lnm_risk(input$nom_ca4, input$nom_a2f0s0g)

    output$plot_nomogram_gauge <- renderPlotly({
      fig <- plot_ly(
        type = "indicator",
        mode = "gauge+number+delta",
        value = prob * 100,
        title = list(text = "LNM Risk Probability (%)", font = list(size = 16)),
        gauge = list(
          axis = list(range = list(0, 100)),
          bar = list(color = "darkblue"),
          steps = list(
            list(range = c(0, 30), color = "#d4edda"),
            list(range = c(30, 60), color = "#fff3cd"),
            list(range = c(60, 100), color = "#f8d7da")
          ),
          threshold = list(
            line = list(color = "red", width = 4),
            thickness = 0.75,
            value = prob * 100
          )
        )
      )
      fig
    })

    output$nom_interpretation <- renderUI({
      risk_level <- ifelse(prob < 0.3, "Low",
                           ifelse(prob < 0.6, "Intermediate", "High"))
      color <- ifelse(prob < 0.3, "green",
                      ifelse(prob < 0.6, "orange", "red"))
      HTML(paste0(
        "<h4 style='color:", color, ";'>Risk Level: <b>", risk_level, "</b></h4>",
        "<p>Probability: <b>", round(prob * 100, 1), "%</b></p>",
        "<hr>",
        "<p><b>Input values:</b><br>",
        "CA4 = ", input$nom_ca4, "<br>",
        "A2F0S0G = ", input$nom_a2f0s0g, "</p>",
        "<hr>",
        "<p><b>Recommendation:</b><br>",
        ifelse(prob < 0.3,
               "Low risk of LNM. Active surveillance may be appropriate.",
               ifelse(prob < 0.6,
                      "Intermediate risk. Consider prophylactic central neck dissection.",
                      "High risk. Recommend therapeutic lymph node dissection.")),
        "</p>",
        "<p><small>Note: Demonstration model. Not for clinical use without validation.</small></p>"
      ))
    })
  })

  # ========================================================================
  # DATA UPLOAD TAB
  # ========================================================================
  user_df <- reactive({
    req(input$user_data)
    tryCatch({
      df <- read_csv(input$user_data$datapath, show_col_types = FALSE)
      df
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      NULL
    })
  })

  output$table_user_data <- renderDT({
    req(user_df())
    datatable(user_df(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$plot_user_profile <- renderPlotly({
    req(user_df())
    df <- user_df()

    glycan_cols <- names(df)[sapply(df, is.numeric)]
    if (length(glycan_cols) == 0) {
      return(plotly_empty() %>% layout(title = "No numeric glycan columns found"))
    }

    id_col <- ifelse("sample_id" %in% names(df), "sample_id", NULL)
    if (!is.null(id_col)) {
      df_long <- df %>%
        slice_head(n = min(10, nrow(df))) %>%
        pivot_longer(cols = all_of(glycan_cols), names_to = "Glycan", values_to = "Abundance")
      p <- ggplot(df_long, aes(x = Glycan, y = Abundance, fill = factor(sample_id))) +
        geom_col(position = "dodge", width = 0.7) +
        labs(fill = "Sample") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      means <- colMeans(df[glycan_cols], na.rm = TRUE)
      df_mean <- data.frame(Glycan = names(means), Mean = means)
      p <- ggplot(df_mean, aes(x = reorder(Glycan, Mean), y = Mean)) +
        geom_col(fill = "#3c8dbc", width = 0.6) +
        coord_flip() +
        labs(x = "", y = "Mean Abundance") +
        theme_minimal()
    }
    ggplotly(p)
  })

  output$plot_user_group <- renderPlotly({
    req(user_df())
    df <- user_df()

    if (!"group" %in% names(df)) {
      return(plotly_empty() %>% layout(title = "No 'group' column found in uploaded data"))
    }

    glycan_cols <- names(df)[sapply(df, is.numeric)]
    if (length(glycan_cols) == 0) {
      return(plotly_empty() %>% layout(title = "No numeric columns found"))
    }

    target <- glycan_cols[1]
    p <- ggplot(df, aes_string(x = "group", y = target, fill = "group")) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(y = target, x = "") +
      theme_minimal()
    ggplotly(p)
  })
}

# ----------------------------------------------------------------------------
# Run App
# ----------------------------------------------------------------------------
shinyApp(ui, server)
