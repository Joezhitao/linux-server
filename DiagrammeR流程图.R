rm(list = ls())

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# 创建改进的流程图（逻辑更清晰）
flow_diagram <- grViz("
digraph research_flow {

  # 设置图的整体属性
  graph [nodesep = 0.6, ranksep = 0.7, splines = true]

  # 定义节点样式
  node [shape = box, style = filled, fontname = Helvetica, fontsize = 14, margin = 0.25, height = 0.6, width = 2.2]
  
  # 主要节点
  A [label = '儿童高血压患者数据收集与处理', color = '#ADD8E6', fontsize = 16]
  
  # 训练和测试相关节点
  B1 [label = '训练集\\n(n=350-400例)', color = '#ADD8E6']
  C1 [label = '数据预处理与特征提取', color = '#E6F2FF']
  D1 [label = '深度学习模型开发', color = '#E6F2FF']
  
  # 测试组相关节点
  B2 [label = '内部测试集\\n(n=150-200例)', color = '#ADD8E6']
  C2 [label = '模型性能评估', color = '#E6F2FF']
  D2 [label = '模型优化与调整', color = '#E6F2FF']
  
  # 最终模型
  E [label = '最终预测模型', color = '#ADD8E6', fontsize = 15]
  
  # 验证组相关节点
  F [label = '前瞻性外部验证集\\n(n=200例)', color = '#ADD8E6']
  G [label = '模型外部验证', color = '#E6F2FF']
  
  # 临床应用相关节点
  H [label = '临床决策支持系统开发', color = '#ADD8E6']
  I1 [label = '风险评估工具', color = '#E6F2FF']
  I2 [label = '预警信号可视化', color = '#E6F2FF']
  I3 [label = '干预策略推荐', color = '#E6F2FF']
  J [label = '临床应用与推广', color = '#ADD8E6', fontsize = 16]
  
  # 节点间的连接
  A -> {B1 B2} [weight=2]
  
  # 训练流程
  B1 -> C1 -> D1
  
  # 测试流程
  B2 -> C2
  D1 -> C2 [label = '测试']
  C2 -> D2
  
  # 形成最终模型
  D2 -> E
  
  # 外部验证
  E -> G
  F -> G [label = '验证数据']
  
  # 验证后连接到临床应用
  G -> H
  H -> {I1 I2 I3}
  {I1 I2 I3} -> J
  
  # 设置节点排列
  {rank = same; B1; B2}
  {rank = same; C1; C2}
  {rank = same; F; G}
  {rank = same; I1; I2; I3}
}
")

# 将图转换为SVG格式
flow_diagram_svg <- DiagrammeRsvg::export_svg(flow_diagram)

# 将SVG保存为PNG（方形比例）
rsvg::rsvg_png(charToRaw(flow_diagram_svg), 
               file = "Hypertension_research_flow_improved.png", 
               width = 1200, height = 1100)  # 接近正方形的尺寸

# 如果需要保存为PDF
rsvg::rsvg_pdf(charToRaw(flow_diagram_svg), 
               file = "Hypertension_research_flow_improved.pdf")
