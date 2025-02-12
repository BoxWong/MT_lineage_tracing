options(scipen = 999)
library(tidyverse)
library(SummarizedExperiment)
#library(Seurat)
library(Matrix)
library(readxl)
library(ComplexHeatmap)
library(circlize) # For heatmap colors
#library(magick) # For heatmap rasterization
library(mclust) # For adjusted rand index
library("ggpubr")
library(treeio)
library(ggtree)


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)


# Assign arguments to variables

args_1 <- args[1]
args_2 = args[2]
args_3 = args[3]
exp_idx = as.numeric(args[4])
mtCN = as.numeric(args[5])
mutRate = args[6]

color_df = read.csv("/data/wangxin/ref/discrete_color_v1.csv", header = T)
sns_col = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
cell_type_pal = c("#E9E3B4", "#AABD8C", "#F39B6D", "#3E6990", "#704E2E")
lineage_pal = c('#577590', '#e76f51', '#4d908e', '#e9c46a', '#bc6c25', '#C58FC9', '#D7C9AA','#264653')
expansion_ls = c("0.1", "0.2", "0.25", "0.5","0.75","0.8","0.9")

dif_mode = args_1
if (dif_mode %in% c("linear_bn", "linear_const")) {
    dif = "linear"
} else {
    dif = "bif"
}


if (dif_mode %in% c("linear_bn", "bifurcated_bn")) {
    bn = "mid"
} else {
    bn = "const"
}

vaf_cutoff=0.01

main_dir = paste0("/syn1/wangxin/work/GB_rev/simul/1209/mt_freq/",dif_mode,"/")
out_dir = "/syn1/wangxin/work/GB_rev/simul/1209/"
simul_id = args_2
gen = args_3
expansion = expansion_ls[exp_idx] 


#Check if the CIV and plots directory exist, if not, create them.
civ_dir = paste0(out_dir, "CIV")
plot_dir = paste0(out_dir, "plots")
dirs <- c(civ_dir, plot_dir)

for (dir in dirs) {
      if (!dir.exists(dir)) {
              dir.create(dir)
        cat(paste("Directory", dir, "has been created.\n"))
      } else {
        cat(paste("Directory", dir, "already exists.\n"))
      }
}

mat_name = paste0("mt_mut_freq_", bn, "_100_", mtCN, "_", mutRate, "_", simul_id, "_", expansion, "_", gen, ".mtx")
cell_name = paste0("mt_mut_freq_", bn, "_100_", mtCN, "_", mutRate, "_", simul_id, "_", expansion, "_", gen, "_cells.csv")
mut_name = paste0("mt_mut_freq_", bn, "_100_", mtCN, "_", mutRate, "_", simul_id, "_", expansion, "_", gen, "_muts.csv")
mat_path = file.path(main_dir, mat_name)
cell_path = file.path(main_dir, cell_name)
mut_path = file.path(main_dir, mut_name)
cell_id = (read.csv(cell_path, row.names = 1))[,"cells"]
mut_id = (read.csv(mut_path, row.names = 1))[,"muts"]

sparse_matrix <- readMM(mat_path)
mat = as.matrix(sparse_matrix)
rownames(mat) = cell_id
colnames(mat) = mut_id
af.dm = as.matrix(t(mat)*100)
rm(mat)
af.dm[af.dm<(vaf_cutoff*100)] = 0
vars.tib <- tibble(var = rownames(af.dm),
                   mean_af = rowMeans(af.dm))

start_time <- Sys.time()
vars.tib = mutate(vars.tib, n0 = apply(af.dm, 1, function(x) sum(x == 0)))
vars.tib = mutate(vars.tib, n1 = apply(af.dm, 1, function(x) sum(x >= 1)))
vars.tib = mutate(vars.tib, n5 = apply(af.dm, 1, function(x) sum(x >= 5)))
vars.tib = mutate(vars.tib, n10 = apply(af.dm, 1, function(x) sum(x >= 10)))
vars.tib = mutate(vars.tib, n20 = apply(af.dm, 1, function(x) sum(x >= 20)))
vars.tib = mutate(vars.tib, n50 = apply(af.dm, 1, function(x) sum(x >= 50)))
Sys.time() - start_time

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Test different variant selection thresholds #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Specify variant selection thresholds to test. voi = variant of interest
conditions.tib <- tibble(min_clone_size = rep(seq(10, 100, by = 10), 5),
                         min_vaf = rep(c("n1", "n5", "n10", "n20", "n50"), each = 10),
                         vois = NA,
                         n_vois = NA,
                         cells = NA,
                         n_cells = NA,
                         transitions = NA)
vois.ls <- vector(mode = "list", length = nrow(conditions.tib))
cells.ls <- vector(mode = "list", length = nrow(conditions.tib))

vars_filter.tib <- vars.tib %>% filter(n0 > 0.9*ncol(af.dm)) # change this from 0.9 to 0.6 to include more variants
# Fill in conditions.tib
for (x in 1:nrow(conditions.tib)) {
    min_clone_size <- conditions.tib$min_clone_size[x]
    min_vaf <- conditions.tib$min_vaf[x]

    # Define variants of which the number of cells exceeding min_vaf is higher than min_clone_size
    voi.ch <- vars_filter.tib$var[vars_filter.tib[[min_vaf]] >= min_clone_size]
    #print(voi.ch)
    if (length(voi.ch) == 0) {
         # Add information to summary table
        conditions.tib[x,"n_vois"] <- 0
        conditions.tib[x,"n_cells"] <- NA
        # Transitions vs. transversions
        conditions.tib[x,"transitions"] <- NA
        
        # Save variants and cells for these parameters
        vois.ls[[x]] <- NA
        #print(positive_cells)
        cells.ls[[x]] <- NA
    }
    else {
        #print(voi.ch)
        # Which cells are positive for at least one of the variants?
        af_subset.dm <- af.dm[voi.ch, ,drop = FALSE]
        positive_cells <- colnames( af_subset.dm[,colSums(af_subset.dm > (vaf_cutoff*100)) > 0 , drop = FALSE] ) # changed the original code
        
        # Add information to summary table
        conditions.tib[x,"n_vois"] <- length(voi.ch)
        conditions.tib[x,"n_cells"] <- length(positive_cells)
        # Transitions vs. transversions
        conditions.tib[x,"transitions"] <- mean( str_count(voi.ch, "G>A|A>G|C>T|T>C") )
        
        # Save variants and cells for these parameters
        vois.ls[[x]] <- voi.ch
        #print(positive_cells)
        cells.ls[[x]] <- positive_cells
        }
}
conditions.tib$vois <- vois.ls
conditions.tib$cells <- cells.ls
conditions.tib

conditions_subset.tib <- conditions.tib %>% filter(min_clone_size %in% c(10,20,50), min_vaf %in% c("n1","n10","n20","n50"))
conditions_subset.tib


# 6. Select variants present in at least 50 cells with a VAF of >20%

a=9
voi.ch <- conditions_subset.tib$vois[[a]]

# List cell IDs that are positive for each of the vois --------------------------------------------
positive_cells.ls <- list()
for (v in voi.ch) {
    # Determine cells with an appreciable VAF
    current_cells.ch <- colnames(af.dm)[af.dm[v,]>(vaf_cutoff*100)]
    # Save cell IDs for positive cells
    positive_cells.ls[[v]] <- current_cells.ch
}
# Make a tibble of cells marked by each voi
positive_cells.tib <- as_tibble(bind_rows(lapply(positive_cells.ls, function(x) data.frame(cell = x)), .id = "variant")[,2:1]) %>%
    mutate(variant = factor(variant, levels = voi.ch))


# Prepare matrix of variants of interest in cells that are positive for at least one
af_voi.mat <- af.dm[voi.ch,]

af_subset.mat = af_voi.mat[, apply(af_voi.mat, 2, function(x) sum(x>1))>0]

# Customize column order. This is different from the strategy for K562 subclones.
plot_order.mat <- af_subset.mat
for (x in rev(voi.ch)) { plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])] }
#for (x in voi.ch) { plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])] }

#options(repr.plot.width=16, repr.plot.height=12)
Heatmap(plot_order.mat,
              col = colorRamp2(seq(0, round(max(plot_order.mat)), length.out = 9),
                               c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
              show_row_names = ifelse(nrow(plot_order.mat) < 100, T, F),
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              row_names_gp = gpar(fontsize = 10),
              name = "AF",
              heatmap_legend_param = list(border = "#000000", grid_height = unit(10, "mm")),
              #top_annotation = ha,
              border = T,
              #width = unit(250, "mm"),
              #height = unit(120, "mm"),
              use_raster = T,
              raster_quality = 5)

cor.mat = cor(t(af_subset.mat))
#cor.mat
var.clust = hclust(as.dist(1 - cor.mat))

plot(var.clust$height, ylim = c(0, max(var.clust$height)))

length(var.clust$height) - sum(var.clust$height<0.8) + 1


ngroups = length(var.clust$height) - sum(var.clust$height<0.8) + 1
Heatmap(cor.mat,
               col = colorRamp2(c(-1,0,1), c("blue", "#DDDDDD", "red")),
               cluster_columns = var.clust,
               cluster_rows = var.clust,
               row_split = switch(ngroups < length(voi.ch), ngroups),
               column_split = switch(ngroups < length(voi.ch), ngroups),
               show_row_dend = F, # without this the visualizationn does not complete
               show_column_dend = F, # without this the visualizationn does not complete
               row_gap = unit(0.5, "mm"),
               column_gap = unit(0.5, "mm"),
               row_names_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 10),
               row_title_gp = gpar(fontsize = 10),
               #width = unit(150, "mm"),
               #height = unit(150, "mm"),
               column_title = ngroups)

Groups.tib <- tibble(var = names(cutree(var.clust, k = ngroups)), Cut = cutree(var.clust, k = ngroups))[var.clust$order,]
Groups.tib <- Groups.tib %>% mutate(Group = match(Cut, unique(Cut)))
Groups.tib <- Groups.tib %>% mutate(Group = match(Cut, unique(Cut)))
Groups.tib <- Groups.tib %>% group_by(Group) %>% summarize(vars = toString(var), nvar = n())
GroupIDs.ls <- lapply(str_split(Groups.tib$vars, ", "), function(x) c(sapply(x, function(y) colnames(af.dm[,af.dm[y,] > 1]))))
Groups.tib$ncells <- unlist(lapply(GroupIDs.ls, function(x) length(unique(unlist(x)))))           
Groups.tib <- Groups.tib %>% arrange(desc(ncells), desc(nvar))

identical(sort(voi.ch),sort(unlist(str_split(Groups.tib$vars, ", "))))

#~~~~~~~~~~~~~#
# VAF heatmap #
#~~~~~~~~~~~~~#

# Sort for all variants from the correlation matrix
plot_order.mat <- af_subset.mat[unlist(str_split(Groups.tib$vars, ", ")),]
# Customize column order.
for (x in rev(strsplit(Groups.tib$vars, ", "))) {
    if (length(x) == 1) {
        plot_order.mat <- plot_order.mat[,order(-plot_order.mat[x,])]
    } else {
        plot_order.mat <- plot_order.mat[,order(-colSums(plot_order.mat[x,]))]
    }
}

Heatmap(plot_order.mat, column_title = paste0(ncol(plot_order.mat)," cells"), 
               col = colorRamp2(seq(0, round(max(plot_order.mat)), length.out = 9),
                                c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               row_names_gp = gpar(fontsize = 10),
               name = "AF",
               heatmap_legend_param = list(border = "#000000", grid_height = unit(5, "mm")),
               #top_annotation = ha,
               border = T,
               #width = unit(100, "mm"),
               #height = unit(100, "mm"),
               use_raster = T,
               raster_quality = 5)

tmp = Groups.tib$vars[Groups.tib$nvar >1]

if (length(tmp) > 0) {Groups.tib$vars[Groups.tib$nvar >1] = sapply(strsplit(tmp, ","), function (x) {
  names(sort(rowMeans(af.dm[unlist(strsplit(gsub(" ", "",x), ",")),]), decreasing = T)[1])
}
)
}
#make a new matrix containing only CIV
CIV = unlist(Groups.tib$vars)
new_plot_mat = plot_order.mat[CIV,]
#sort the new CIV VAF matrix
for (x in rev(strsplit(Groups.tib$vars, ", "))) {
  if (length(x) == 1) {
    new_plot_mat <- new_plot_mat[,order(-new_plot_mat[x,])]
  } else {
    new_plot_mat <- new_plot_mat[,order(-colSums(new_plot_mat[x,]))]
  }
}


right_end = sapply(CIV, function(x) tail(which(new_plot_mat[x,]>0), 1))
right_bor = 0
CIV_idx = c() 
for (i in 1:length(right_end)) {
  if (right_end[i]>right_bor) {
    right_bor = right_end[i]
    CIV_idx = c(CIV_idx, i)
  }}
right_end = right_end[CIV_idx]
CIV = CIV[CIV_idx]
left_end = c(1, sapply(right_end[1:length(right_end)-1], function(x) x+1))
cell_num = right_end - left_end + 1
CIV = CIV[cell_num>0] #remove CIVs that are misclustered
new_plot_mat = new_plot_mat[CIV,]
right_end = sapply(CIV, function(x) tail(which(new_plot_mat[x,]>0), 1))
left_end = c(1, sapply(right_end[1:length(right_end)-1], function(x) x+1))
cell_num = right_end - left_end + 1

new1 = new_plot_mat[CIV[rev(order(cell_num))],]
my_command = paste0("order(", paste0("-new1[", 1:nrow(new1), ",]", collapse = ","), ")")
new2 = new1[, eval(parse(text = my_command))]
new2 = new2[ , colSums(new2)>0]
CIV = rownames(new2)
right_end = sapply(CIV, function(x) tail(which(new2[x,]>0), 1))
left_end = c(1, sapply(right_end[1:length(right_end)-1], function(x) x+1))
cell_num = right_end - left_end + 1                     
#remove clones smaller than 5
CIV = CIV[cell_num >= 5]
cell_num = cell_num[cell_num >= 5]
new_plot_mat = new2[CIV, 1:sum(cell_num)]
clone_info_df = data.frame(Cell = colnames(new_plot_mat), Clone = paste0("C", rep(1:length(CIV), cell_num)))

CIV_file_name = paste0(civ_dir,"/CIV_clones_", simul_id, "_", mtCN, "_", mutRate, "_", gen, "gen_", expansion, ".csv")
write.csv(clone_info_df, file = CIV_file_name, quote = F, row.names = F)
message("Clone info has been saved to: ", CIV_file_name)

#new_plot_mat = new2
new_plot_mat1 = new_plot_mat/100
cols = c("#136fb0", "#ff7f0e")
names(cols) = c("pre-existing", "de novo")
col_ready = list(MutType=cols)
pre_df = read.csv("/syn1/wangxin/work/GB_rev/simul/1209/pre_id/pre_id.csv", head=F, row.names=1)
lookup_id = paste0(dif_mode, "_", simul_id, "_", mtCN, "_", mutRate)
mut_state = ifelse(as.numeric(rownames(new_plot_mat)) < pre_df[lookup_id,], "pre-existing", "de novo")
CIV_type_df = data.frame(CIV = CIV, Type = mut_state, Simulation = rep(simul_id, length(CIV)), Model = expansion)
CIV_type_file_name = paste0(civ_dir, "/mut_type_", lookup_id, "_", gen, "gen_", expansion, ".csv")
write.csv(CIV_type_df, file = CIV_type_file_name, quote = F, row.names = F)
message("Clone info has been saved to: ", CIV_type_file_name)

                       
row_annot = rowAnnotation(MutType = mut_state, col=col_ready, annotation_legend_param = list(
        MutType = list(
                title = "MutType",
                at = c("pre-existing", "de novo"),
                labels = c("pre-existing", "de novo")
            )))

                       
Heatmap(new_plot_mat, column_title = paste0(ncol(new_plot_mat)," cells"), 
               col = colorRamp2(seq(0, 100, length.out = 9),
                                c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               row_names_gp = gpar(fontsize = 10),
               name = "AF",
               heatmap_legend_param = list(border = "#000000", grid_height = unit(5, "mm")),
               left_annotation = row_annot,
               #top_annotation = ha,
               border = T,
               #width = unit(100, "mm"),
               #height = unit(100, "mm"),
               use_raster = T,
               raster_quality = 5)

rownames(new_plot_mat) = paste0("Mut-", rownames(new_plot_mat))
p1.1 = Heatmap(new_plot_mat, column_title = paste0(ncol(new_plot_mat)," cells"), 
               col = colorRamp2(seq(0, 100, length.out = 9),
                                c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               row_names_gp = gpar(fontsize = 10),
               name = "VAF",
               heatmap_legend_param = list(border = "#000000", grid_height = unit(5, "mm")),
               left_annotation = row_annot,
               #top_annotation = ha,
               border = T,
               #width = unit(100, "mm"),
               #height = unit(100, "mm"),
               use_raster = T,
               raster_quality = 5)
print(p1.1)

rownames(new_plot_mat) = paste0("C", 1:length(CIV))
p1.2 = Heatmap(new_plot_mat, column_title = paste0(ncol(new_plot_mat)," cells"), 
               col = colorRamp2(seq(0, 100, length.out = 9),
                                c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               row_names_gp = gpar(fontsize = 10),
               name = "VAF",
               heatmap_legend_param = list(border = "#000000", grid_height = unit(5, "mm")),
               left_annotation = row_annot,
               #top_annotation = ha,
               border = T,
               #width = unit(100, "mm"),
               #height = unit(100, "mm"),
               use_raster = T,
               raster_quality = 5)
print(p1.2)

rownames(new_plot_mat) = paste0("Mut-", CIV, "/", paste0("C", 1:length(CIV)))
p1.3 = Heatmap(new_plot_mat, column_title = paste0(ncol(new_plot_mat)," cells"), 
               col = colorRamp2(seq(0, 100, length.out = 9),
                                c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
               show_row_names = T,
               show_column_names = F,
               cluster_columns = F,
               cluster_rows = F,
               row_names_gp = gpar(fontsize = 10),
               name = "VAF",
               heatmap_legend_param = list(border = "#000000", grid_height = unit(5, "mm")),
               left_annotation = row_annot,
               #top_annotation = ha,
               border = T,
               #width = unit(100, "mm"),
               #height = unit(100, "mm"),
               use_raster = T,
               raster_quality = 5)
print(p1.3)


heatmap_name = paste0(plot_dir,"/", lookup_id, "_", gen, "_gen_", expansion, "_CIV_heatmap.pdf")
pdf(heatmap_name, height = 2.645, width = 3.5)
p1.1
p1.2
p1.3
dev.off()
message(paste0("CIV heatmaps have been saved to: ", heatmap_name))


### Distribution of clones on phylogenetic trees
CIV_num = length(CIV)
color_df = read.csv("/data/wangxin/ref/discrete_color_v2.csv", header = T)
pal1 = color_df$color[1:CIV_num]
names(pal1) = paste0("C", 1:CIV_num)
  
if (expansion == "1.0") {
  tree_name = paste0(dif,"_tree_gt_", simul_id, ".nwk")
  tree = read.newick(paste0(out_dir,"trees/modified_trees/",tree_name))
  clone_df = data.frame(Clone=clone_info_df$Clone, row.names = sapply(clone_info_df$Cell, function(y) {paste0(strsplit(y,">")[[1]][1], ">")} ))
  #sample 1K cells from all cells and remove cells with no CIV
  cells_1k = sample(colnames(af.dm), 1000, replace = F)
  cells_1k = sapply(cells_1k, function(y) {paste0(strsplit(y,">")[[1]][1], ">")} )
  to_drop = setdiff(sapply(colnames(af.dm), function(y) {paste0(strsplit(y,">")[[1]][1], ">")} ), cells_1k)
  to_drop = c(to_drop, setdiff(cells_1k,rownames(clone_df)))
  new_tree = drop.tip(tree, to_drop)
  } else {
  ## e.g. clonal_expansion_tree_173076_0.95_800.nwk
  tree_name = paste0("mt_allmuts_", bn, "_100_", mtCN, "_", mutRate, "_", simul_id, "_", expansion, "_", gen, ".mnwk")
  tree = read.newick(paste0(out_dir,"trees/modified_trees/",tree_name))
  #print("read tree, done.")
  clone_df = data.frame(Clone=clone_info_df$Clone, row.names = clone_info_df$Cell)
  if (nrow(clone_df) >1000) {
      cells_to_keep = sample(rownames(clone_df), 1000, replace = F)}
  else {
    cells_to_keep = rownames(clone_df)
  }
  to_drop = setdiff(tree$tip.label, cells_to_keep)
  new_tree = drop.tip(tree, to_drop)
  
  #subset top10 CIV-defined clones if the number of CIV-defined clones exceed 10, and sample 1000 cells if cell population is bigger than 1000.
    if (length(CIV) > 10) {
        CIV_top = paste0("C",1:10)}
    else {
        CIV_top = paste0("C",1:length(CIV))}
    clone_df_top = data.frame( Clone = clone_info_df[clone_info_df$Clone %in% CIV_top, "Clone"], row.names = clone_info_df[clone_info_df$Clone %in% CIV_top, "Cell"])
    if (nrow(clone_df_top) >1000) {
        cells_to_keep = sample(rownames(clone_df_top), 1000, replace = F)}
    else {
        cells_to_keep = rownames(clone_df_top)
    }

    to_drop = setdiff(tree$tip.label, cells_to_keep)
    new_tree_top = drop.tip(tree, to_drop)
}
  #print("new tree, done.")
new_tree_top_name = paste0(civ_dir, "/top10_tree_downsampled_", bn, "_100_", mtCN, "_", mutRate, "_", simul_id, "_", expansion, "_", gen, ".nwk")
write.tree(new_tree_top, file = new_tree_top_name)

##plot downsampled tree
p =  ggtree(new_tree,  layout="circular", branch.length = "none") + geom_tree(size = 0.5)



#pal1.1 = sample(color_df$color)
p2 = gheatmap(p, clone_df, width=.1, colnames = FALSE, color=NA) + theme(legend.position ="right",legend.key.size = unit(0.5, 'cm'))+
      scale_fill_manual(values=pal1, na.translate = F, breaks =paste0("C",1:length(CIV)),  labels = paste0("C",1:length(CIV)), name = "Clone")
p2
  
tree_pdf_name = paste0(plot_dir,"/celldivision_tree_", lookup_id, "_",gen,"_gen_", expansion, ".pdf")
ggsave(tree_pdf_name, p2, height = 2.3, width = 5)


p2.1 = ggtree(new_tree_top,  layout="circular", branch.length = "none") + geom_tree(size = 0.5)
p2.2 = gheatmap(p2.1, clone_df_top, width=.1, colnames = FALSE, color=NA) + theme(legend.position ="right",legend.key.size = unit(0.5, 'cm'))+
        scale_fill_manual(values=pal1, na.translate = F, breaks =paste0("C",1:length(CIV_top)),  labels = paste0("C",1:length(CIV_top)), name = "Clone")

#write.csv(data.frame(name=1:length(pal1.1), color=pal1.1), file = "tmp_color.csv", quote=F)


top_clone_tree_pdf_name = paste0(plot_dir,"/top10_clones_celldivision_tree_", lookup_id, "_",gen,"_gen_", expansion, ".pdf")
ggsave(top_clone_tree_pdf_name, p2.2, height = 2.3, width = 5)
message(paste0("Cell divsion tree has been saved to: ", tree_pdf_name, " and ", top_clone_tree_pdf_name, "."))
#quit()




#search the MRCA of CIV-defined subpopulations
tree_name = paste0(dif,"_tree_gt_", simul_id, ".nwk")
tree = read.newick(paste0(out_dir,"/trees/modified_trees/", tree_name))
tib = as_tibble(tree)


all_mean_ls = c()
all_var_ls = c()
sub_mean_ls = c()
sub_var_ls = c()
all_mrca_ls = c()
sub_mrca_ls = c()
all_mv1_ls = c()
sub_mv1_ls = c()
#new_plot_mat1 = new_plot_mat/100
#rownames(new_plot_mat1)
#CIV
for (x in 1:length(CIV)) {
        #all positive cells
        ori_cells = colnames(new_plot_mat1)[new_plot_mat1[CIV[x],] >= vaf_cutoff]
        tmp_vaf_ls = new_plot_mat1[CIV[x],ori_cells]
        all_mean_ls = append(all_mean_ls, mean(tmp_vaf_ls))
        all_var_ls = append(all_var_ls, var(tmp_vaf_ls))
        all_mv1_ls = append(all_mv1_ls, mean(tmp_vaf_ls)/(1+var(tmp_vaf_ls)))
        uniq_cells = unique(sapply(ori_cells, function(y) {paste0(strsplit(y,">")[[1]][1], ">")} ))
        str1 = strsplit(as.character(tib[(tib$node == MRCA(tree, uniq_cells)),"label"]), "_")[[1]][1]
        cur_mrca = as.numeric(substr(str1, 2, nchar(str1)))
        all_mrca_ls = append(all_mrca_ls, cur_mrca)
        
        #sub cells
        ori_cells = clone_info_df$Cell[clone_info_df$Clone == paste0("C", x)]
        tmp_vaf_ls = new_plot_mat1[CIV[x], ori_cells]
        tmp_vaf_ls = tmp_vaf_ls[tmp_vaf_ls >= vaf_cutoff]
        ori_cells = ori_cells[tmp_vaf_ls >= vaf_cutoff]
        sub_mean_ls = append(sub_mean_ls, mean(tmp_vaf_ls))
        sub_var_ls = append(sub_var_ls, var(tmp_vaf_ls))
        sub_mv1_ls = append(sub_mv1_ls, mean(tmp_vaf_ls)/(1+var(tmp_vaf_ls)))
        uniq_cells = unique(sapply(ori_cells, function(y) {paste0(strsplit(y,">")[[1]][1], ">")} ))
        str1 = strsplit(as.character(tib[(tib$node == MRCA(tree, uniq_cells)),"label"]), "_")[[1]][1]
        cur_mrca = as.numeric(substr(str1, 2, nchar(str1)))
        sub_mrca_ls = append(sub_mrca_ls, cur_mrca)
}


meanvar_df_all = data.frame("Clone" = CIV, "Mean" = all_mean_ls, "Variance" = all_var_ls,"MV1" = all_mv1_ls, "MRCA" = all_mrca_ls, "Sim" = rep(simul_id, length(CIV)), "Type" = mut_state, "Method" = rep("All", length(CIV)))
meanvar_df_sub= data.frame("Clone" = CIV, "Mean" = sub_mean_ls, "Variance" = sub_var_ls, "MV1" = sub_mv1_ls, "MRCA" = sub_mrca_ls, "Sim" = rep(simul_id, length(CIV)), "Type" = mut_state, "Method" = rep("Sub", length(CIV)))




### Clone composition

#construct python command to trace lineages of each CIV-defined clone for All mode including noise.
message("Analyze lineage composition and cell type composition...")

original_tree_file = paste0(out_dir,"res/", dif_mode, "/", simul_id, "/", tree_name)
CIV_pos_cell_df = as.data.frame(matrix(ncol=2, nrow=0))
colnames(CIV_pos_cell_df) = c("Cell", "Clone")
for (i in 1:length(CIV)) {
    tmp_df = data.frame("Cell" = colnames(new_plot_mat1)[new_plot_mat1[i,] > vaf_cutoff], "Clone" = paste0("C",i))
    CIV_pos_cell_df =  rbind(CIV_pos_cell_df, tmp_df)                   
}
CIV_pos_cell_file = paste0(civ_dir,"/CIV_postive_cells_", lookup_id, "_", gen, "gen_", expansion, ".csv")
write.csv(CIV_pos_cell_df, file = CIV_pos_cell_file, quote = F, row.names = F)
message("Clone composition info has been saved to: ", CIV_pos_cell_file)


#lookup_id = paste0(dif_mode, "_", simul_id, "_", mtCN, "_", mutRate)
all_cc_output_file = paste0(civ_dir, "/CIV_composition_", lookup_id, "_", gen, "gen_", expansion, "_all.csv")
all_anc_output_file = paste0(civ_dir, "/anc_composition_", lookup_id, "_", gen, "gen_", expansion, "_all.csv")
cc_command = paste0("python /syn1/wangxin/work/mtsim/0421/src/ancestor_tracing.py ", simul_id, " ", gen, " ", CIV_pos_cell_file, " ", original_tree_file, " ",  all_cc_output_file, " ", all_anc_output_file)
message(paste0("Executing ",cc_command))
system(cc_command)

#construct python command to trace lineages of each CIV-defined clone for Sub mode
cc_output_file = paste0(civ_dir, "/CIV_composition_", lookup_id, "_", gen, "gen_", expansion, "_sub.csv")
anc_output_file = paste0(civ_dir, "/anc_composition_", lookup_id, "_", gen, "gen_", expansion, "_sub.csv")
cc_command = paste0("python /syn1/wangxin/work/mtsim/0421/src/ancestor_tracing.py ", simul_id, " ", gen, " ", CIV_file_name, " ", original_tree_file, " ",  cc_output_file, " ", anc_output_file)
message(paste0("Executing ",cc_command))
system(cc_command)




#combine mean, var, mv1, MRCA and 
all_tmp_df = read.csv(all_anc_output_file, check.names =F, row.names = 1) 
sub_tmp_df = read.csv(anc_output_file, check.names =F, row.names = 1)
meanvar_df_all$"AG_0.6" = all_tmp_df[, "Ancestor_gen_0.6"]
meanvar_df_all$"AC_0.6" = all_tmp_df[, "Ancestor_0.6"]
meanvar_df_sub$"AG_0.6" = sub_tmp_df[, "Ancestor_gen_0.6"]
meanvar_df_sub$"AC_0.6" = sub_tmp_df[, "Ancestor_0.6"]
meanvar_df_all$"AG_0.7" = all_tmp_df[, "Ancestor_gen_0.7"]
meanvar_df_all$"AC_0.7" = all_tmp_df[, "Ancestor_0.7"]
meanvar_df_sub$"AG_0.7" = sub_tmp_df[, "Ancestor_gen_0.7"]
meanvar_df_sub$"AC_0.7" = sub_tmp_df[, "Ancestor_0.75"]
meanvar_df_all$"AG_0.75" = all_tmp_df[, "Ancestor_gen_0.75"]
meanvar_df_all$"AC_0.75" = all_tmp_df[, "Ancestor_0.75"]
meanvar_df_sub$"AG_0.75" = sub_tmp_df[, "Ancestor_gen_0.75"]
meanvar_df_sub$"AC_0.75" = sub_tmp_df[, "Ancestor_0.75"]
meanvar_df_all$"AG_0.8" = all_tmp_df[, "Ancestor_gen_0.8"]
meanvar_df_all$"AC_0.8" = all_tmp_df[, "Ancestor_0.8"]
meanvar_df_sub$"AG_0.8" = sub_tmp_df[, "Ancestor_gen_0.8"]
meanvar_df_sub$"AC_0.8" = sub_tmp_df[, "Ancestor_0.8"]




meanvar_df = rbind(meanvar_df_all, meanvar_df_sub)
meanvar_df$Model = expansion
meanvar_file_name = paste0(civ_dir,"/meanvar_CIV_", lookup_id, "_", gen, "gen_", expansion, ".csv")
write.csv(meanvar_df, file = meanvar_file_name, quote = F, row.names = F)
message("Mean, variance and MRCA info have been saved to: ", meanvar_file_name)
#message("Success till here, done.")
#quit()

compo_df = read.csv(cc_output_file, check.names = F)
CIV_num = length(unique(compo_df$clone))
compo_df$clone = factor(compo_df$clone, levels = paste0("C",1:CIV_num))

p3 = ggplot(compo_df, aes(x = clone, y=counts, fill = ancestor)) + geom_bar(stat="identity") + theme_classic() + labs(title = paste0(lookup_id, " - ", gen, expansion)) + 
        theme(plot.title=element_text(hjust=0.5, size=18), axis.text = element_text(size = 14), axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 270, hjust = 0.1, vjust = 0.5)) + scale_fill_manual(values=lineage_pal) + xlab("")+ylab("Cell number") + theme(legend.position = "none") 
p3

clone_cell_vec = c() # create a vector for the total number of cells in each clone.
for (i in 1:CIV_num) {
    cell_sum = sum(compo_df[compo_df$clone == paste0("C",i), "counts"])
    clone_cell_vec = append(clone_cell_vec, cell_sum)
}
names(clone_cell_vec) = paste0("C",1:CIV_num)
compo_df$Proportion = apply(compo_df, 1, function (x) {as.numeric(x[2])/clone_cell_vec[x[3]]})

p4 = ggplot(compo_df, aes(x = clone, y=Proportion, fill = ancestor)) + geom_bar(stat="identity") + theme_classic() + labs(title = paste0(lookup_id, " - ", gen, expansion)) + 
        theme(plot.title=element_text(hjust=0.5, size=9), axis.text = element_text(size = 9), axis.title = element_text(size=9),
        axis.text.x = element_text(angle = 270, hjust = 0.1, vjust = 0.5)) + scale_fill_manual(values=lineage_pal) + xlab("")+ylab("Proportion") + theme(legend.position = "none") +  
        annotate("text", x = paste0("C", 1:length(CIV)), y = 1.05, label = as.character((table(compo_df$clone))[paste0("C", 1:length(CIV))]), size =6)
p4

cc_pdf_file = paste0(plot_dir, "/lineage_composition_", lookup_id, "_", gen, "gen_",expansion,".pdf")
ggsave(cc_pdf_file, p4, height = 2.2, width = 2.2)
message(paste0("PDF has been written to: ",cc_pdf_file))

meta_file_name = paste0(out_dir, "/res/", dif_mode, "/", simul_id, "/tree_origin_", dif, "_", simul_id, ".csv")
cell_type_meta = read.csv(meta_file_name, check.names = F)
meta_df = data.frame(row.names = paste0("<",cell_type_meta[,1], "_", cell_type_meta[,2], ">"), "type" = cell_type_meta[,4])
cell_type_df = data.frame("Clone"=NULL, "Type"=NULL, "Counts"=NULL)
#cell_type_pal = c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")

for (i in 1:CIV_num) {
    query_cells = clone_info_df[clone_info_df$Clone == paste0("C",i),"Cell"]
    query_cells = sapply(query_cells, function(y) {paste0(strsplit(y,">")[[1]][1], ">")} )
    cell_type_stat = table(meta_df[query_cells,"type"])
    tmp_cell_type = as.data.frame(cell_type_stat)
    tmp_cell_type$Clone = paste0("C",i)
    colnames(tmp_cell_type) = c("Type", "Counts", "Clone")
    cell_type_df = rbind(cell_type_df, tmp_cell_type)
    }

cell_type_df$Clone = factor(cell_type_df$Clone, levels = paste0("C",1:CIV_num))
cell_type_df$Type = factor(cell_type_df$Type, levels = as.character(0:4))

p5 = ggplot(cell_type_df, aes(x = Clone, y= Counts, fill = Type)) + geom_bar(stat="identity") + theme_classic() + labs(title = paste0(lookup_id, " - ", gen, expansion)) + 
        theme(plot.title=element_text(hjust=0.5, size=9), axis.text = element_text(size = 9), axis.title = element_text(size=9),
        axis.text.x = element_text(angle = 270, hjust = 0.1, vjust = 0.5)) + scale_fill_manual(values= cell_type_pal) + xlab("")+ylab("Cell number") + theme (legend.position = "none")
p5

clone_num_vec = table(clone_info_df$Clone)
cell_type_df$Proportion = apply(cell_type_df, 1, function (x) {as.numeric(x[2])/clone_num_vec[x[3]]})
p6 = ggplot(cell_type_df, aes(x = Clone, y= Proportion, fill = Type)) + geom_bar(stat="identity") + theme_classic() + labs(title = paste0(lookup_id, " - ", gen, expansion)) + 
        theme(plot.title=element_text(hjust=0.5, size=9), axis.text = element_text(size = 9), axis.title = element_text(size=9),
        axis.text.x = element_text(angle = 270, hjust = 0.1, vjust = 0.5)) + scale_fill_manual(values= cell_type_pal) + xlab("")+ylab("Proportion") + theme (legend.position = "none")
p6

ct_pdf_file = paste0(plot_dir, "/celltype_composition_", lookup_id, "_", gen, "gen_",expansion,".pdf")
ggsave(ct_pdf_file, p6, height = 2.2, width = 2.2)
message(paste0("PDF has been written to: ", ct_pdf_file))
#message("Success till here")
#quit()


### tSNE
message("Plot tSNEs...")
if (expansion != "1.0") {
    for (ds in c(10, 20, 30, 40, 50)){

        tsne_file = paste0(out_dir, "/tsne/tsne_", simul_id, "_", mtCN, "_", mutRate, "_", gen, "gen_", expansion, "_dis", ds, ".csv")
        tsne_df = read.csv(tsne_file, check.names = F, row.names =1)
        colnames(tsne_df) = c("tSNE1", "tSNE2")
        plot_ls = lapply(1:CIV_num, function (x) {
            Clone = paste0("C", x)
            pos_cells = rev(clone_info_df[clone_info_df$Clone == Clone, "Cell"])
            nega_cells = setdiff(rownames(tsne_df), pos_cells)
            new_cell_order = c(nega_cells, pos_cells)
            new_tsne_df = tsne_df[new_cell_order,]
            new_tsne_df$vaf = 0
            new_tsne_df[pos_cells,"vaf"] = new_plot_mat[x, pos_cells]
            graph = ggplot(new_tsne_df, aes(x=tSNE1, y= tSNE2, color = vaf)) + geom_point() + theme_void() + scale_color_continuous(low="gray", high="red") + 
            theme(legend.position = "none") + labs(title=paste0("C",x)) + 
            theme(plot.title=element_text(hjust=0.5, size=18))
            graph
        })
    tsne_pdf_file = paste0(plot_dir, "/eachCIV_on_tsne_",lookup_id, "_", gen, "gen_", expansion, "_dis", ds, ".pdf")
    ggexport(plotlist = plot_ls, nrow=4, ncol=4, filename = tsne_pdf_file, width = 20, height = 20)
    }
} else{
  tsne_file = paste0(out_dir, "/res/", dif_mode, "/", simul_id, "/tsne_", dif, "_", simul_id, ".csv")
  tsne_df = read.csv(tsne_file, check.names = F, row.names =1)
  colnames(tsne_df) = c("tSNE1", "tSNE2")
  plot_ls = lapply(1:CIV_num, function (x) {
      Clone = paste0("C", x)
      
      ori_pos_cells = rev(clone_info_df[clone_info_df$Clone == Clone, "Cell"])
      pos_cells = sapply(ori_pos_cells, function(y) {paste0(strsplit(y,">")[[1]][1], ">")} )
      nega_cells = setdiff(rownames(tsne_df), pos_cells)
      new_cell_order = c(nega_cells, pos_cells)
      new_tsne_df = tsne_df[new_cell_order,]
      new_tsne_df$vaf = 0
      new_tsne_df[pos_cells,"vaf"] = new_plot_mat[x, ori_pos_cells]
      graph = ggplot(new_tsne_df, aes(x=tSNE1, y= tSNE2, color = vaf)) + geom_point() + theme_void() + scale_color_continuous(low="gray", high="red") + 
          theme(legend.position = "none") + labs(title=paste0("C",x)) + 
          theme(plot.title=element_text(hjust=0.5, size=18))
      graph
      })
  tsne_pdf_file = paste0(plot_dir, "/eachCIV_on_tsne_",lookup_id, "_", gen, "gen_", expansion, ".pdf")
  ggexport(plotlist = plot_ls, nrow=4, ncol=4, filename = tsne_pdf_file, width = 20, height = 20)

}


print(paste0("Successfully done: ", lookup_id, "_", gen, "_gen_", expansion, "_expansion."))


