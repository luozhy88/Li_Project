





library(DiagrammeR)
path <- glue::glue("../output/DiagrammeR/") 
dir.create(path, recursive = TRUE)

pp=grViz("
  digraph flowchart {
    # 图形属性
    graph [rankdir = TB, bgcolor = white]
    
    # 节点样式
    node [shape = rectangle, style = 'rounded,filled', fillcolor = white, color = red, penwidth = 2]
    
    # 定义节点
    A [label = 'Rawdata(fq)', fixedsize = true, width = 2]
    B [label = 'TOFU', fixedsize = true, width = 2]
    C [label = 'Kraken2', fixedsize = true, width = 2]
    D [label = 'Generalized\nLinear Model']
    E [label = 'eXtremely\nRandomized\nTrees']
    F [label = 'Distributed\nRandom Forest']
    G [label = 'DeepLearning']
    H [label = 'Feature importance', shape = box, style = filled, fillcolor = '#000080', fontcolor = white]
    I [label = 'Genus', shape = box, style = 'dashed,filled', fillcolor = lightgray]
    J [label = 'Species', shape = box, style = 'dashed,filled', fillcolor = lightgray]
    
    # 定义连接
    A -> B -> C
    C -> {D E F G}
    {D E F G} -> H
    H -> {I J} 
    
    # 布局调整
    {rank = same; D; E; F; G}
    {rank = same; I}
    {rank = same; J}
  }
")




library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(glue)

# 假设 pp 是您的 DiagrammeR 图形对象
# 如果还没有创建，请先创建它

# 设置更高的宽度和高度
width <- 1200  # 增加宽度
height <- 1200 # 增加高度
dpi <- 300     # 设置更高的DPI值

# 导出为SVG，然后转换为高分辨率PNG
png_file <- DiagrammeRsvg::export_svg(pp) %>%
  charToRaw() %>%
  rsvg::rsvg_png(glue::glue("{path}/flowchart.png"),
                 width = width,
                 height = height)




# # 将图形转换为SVG
# svg_xml <- DiagrammeRsvg::export_svg(pp)
# svg <- charToRaw(svg_xml)
# 
# # 保存为PNG文件
# rsvg::rsvg_png(svg, "flowchart.png", width = 800)

