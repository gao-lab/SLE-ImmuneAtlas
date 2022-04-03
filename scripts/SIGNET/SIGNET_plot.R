# scatter plot
plot.scatter <- function(dataplot){
  p <- ggplot(dataplot, aes(UMAP_1, UMAP_2, color = CellType)) +
  geom_point(size =1, shape = 19, alpha = 0.5) + 
  scale_colour_brewer(palette="Paired") + 
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "TNM", face = "bold", size = 12),
        legend.text = element_text(family = "TNM",size = 10),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 
  return(p)
}

# heatmap
plot.heatmap <- function(data,meta,color){
	ann_col = data.frame(CellType = meta)
    rownames(ann_col) = colnames(data)
    color_label = color
    names(color_label) = levels(as.factor(meta))
	ph = pheatmap(data, scale = "row", 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation",
         treeheight_row = 30, treeheight_col = 0,
         #cutree_rows = 3, cutree_cols = 4,
         #cellwidth = 6, cellheight = 5, 
         #color = c("#FFFFFF","#000000"), 
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         color = colorRampPalette(c("#3c5488","white","#dc0000"))(25),
         annotation_col = ann_col,
         annotation_colors = list(CellType = color_label),
         annotation_names_col = F,
         annotation_legend = T,
         show_colnames = F, border_color = NA, legend = T,
         #fontfamily= "Helvetica Neue CE 55 Roman", #fontface="bold",
         fontsize = 9, fontsize_row = 6)
	return(ph)
}

# network
# mapping function
map = function(x, mat){
  y = numeric()
  n = length(x)
  for (i in 1:n) {
    if (x[i] %in% mat[,1]){
      j = min(which(x[i]==mat[,1]))
      y[i] = mat[j,2]
    }
    else {
      y[i] = NA
    }
  }
  return(y)
}

levelscount = function(x, output = "df"){
  if(!is.factor(x)){
    stop("input requires a factor")
  }
  lev = levels(x)
  l = length(lev)
  y = numeric()
  for (i in 1:l) {
    y[i] = length(which(x == lev[i]))
  }
  if (output == "v"){
    names(y) = lev
    return(y)
  }
  if (output == "df"){
    df = data.frame(levels = lev, counts = y)
    return(df)
  }
}


# generate edge or node list 
getNetElement = function(data, center, color, mode = "plane", ...){
  if(length(center) > length(color)){
    stop("The length of palette is shorter than the center")
  }
  n = length(center)
  colormap = data.frame(center = center, color[1:n])
  # edgelist
  row.id = data[,1] %in% center
  edgelist = data[row.id,]
  edgecolor = map(edgelist[,1], colormap)
  edgelist = data.frame(edgelist, edgecolor)
  
  # nodelist
  from = as.character(levels(as.factor(edgelist[,1])))
  to = as.character(levels(as.factor(edgelist[,2])))
  dif = setdiff(to, from)
  ntf = levelscount(as.factor(edgelist[,2]))
  to2 = ntf[ntf$counts>=2,1]
  dif2 = setdiff(to2, from)
  m = max(ntf$counts)
  # plane mode
  if(mode == "plane"){
    nodeid = c(from, dif)
    nodegroup = c(from, rep("NTF", length(dif)))
  }
  # multiple nodes highlight
  if(mode == "crosshighlight_truncated"){
    # node list
    dif = setdiff(dif, dif2)
    nodeid = c(from, dif, dif2)
    nodegroup = c(from, rep("NTF", length(dif)), rep("MNTF", length(dif2)))
  }
  # multiple nodes highlight version 2
  if(mode == "crosshighlight"){
    # node list
    dif = setdiff(dif, dif2)
    nodeid = c(from, dif)
    nodegroup = c(from, rep("NTF", length(dif)))
    if(m > 1){
      for (i in 2:m) {
        tom = ntf[ntf$counts==i,1]
        difm = setdiff(tom, from)
        nodeid = c(nodeid, difm)
        label = paste0("MNTF", i, sep = "")
        nodegroup = c(nodegroup, rep(label, length(difm)))
      }
    }
  }
  nodelist = data.frame(nodeid, nodegroup)
  
  # network object
  net = igraph::graph_from_data_frame(d=edgelist, vertices=nodelist, directed=T, ...)
  return(list(edge = edgelist, node = nodelist, cross_max = m, net = net))
}


plot.net <- function(net,...){
	p = ggnet2(net, mode = "kamadakawai", 
       shape = 19, alpha = 0.85, 
       size = "nodegroup", size.palette = nodesize, 
       color = "nodegroup", color.palette = nodecol, color.legend = "Group", 
       label = T, label.size = 1.25, label.color = "white", label.alpha = 1, 
       edge.alpha = 0.75, edge.color = "edgecolor", edge.size = 0.5,
       arrow.size = 8, arrow.gap = 0.01,
       legend.size = 8, legend.position = "bottom",...) +
    guides(size = F)
    return(p)
}

