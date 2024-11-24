library(DiagrammeR)

plot1=grViz("
digraph flowchart {
  # Node attributes
  node [shape = rectangle, style = filled, fillcolor = lightblue, fontname = 'Arial']
  
  # Define nodes
  herbs [label = 'Traditional Chinese Medicine\n(ZBM,\nXKC,\nXYT)']
  tcmsp [label = 'TCMSP Database\nMetabolite Information']
  pubchem [label = 'PubChem Database\nCID Information']
  genes1 [label = 'Gene Information\n(Meta_Bat/BindingDB/\nGuideToPharmacology)']
  genes2 [label = 'GeneCards Database\nThyroid Nodule Related Genes']
  intersection [label = 'Gene Intersection Analysis']
  
  subgraph cluster_analysis {
    label = 'Downstream Analysis'
    style = dashed
    
    enrichment [label = 'GO/KEGG Enrichment Analysis\n(clusterProfiler)']
    docking [label = 'Molecular Docking\n(AutoDock Vina)']
    visual [label = 'Results Visualization']
  }
  
  # Define connections
  herbs -> tcmsp
  tcmsp -> pubchem
  pubchem -> genes1
  genes2 -> intersection
  genes1 -> intersection
  intersection -> enrichment
  intersection -> docking
  enrichment -> visual
  docking -> visual
  
  # Graph attributes
  graph [rankdir = TB, splines = ortho]
}
")
dir.create("output")
# 将图形保存为PNG格式
library(rsvg)

# 使用正确的导出方法
plot1 %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg::rsvg_png("output/flowchart.png", width = 800, height = 600)



