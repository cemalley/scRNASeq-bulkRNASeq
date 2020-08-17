library(devtools)
devtools::install_github("dynverse/dyno", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynmethods", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynplot", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynutils", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynwrap", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynguidelines", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynfeature", force=TRUE, ref = 'master')


library(dyno)
library(dynutils)
library(tidyverse)
library(data.table)
options(rgl.useNULL=TRUE)
.rs.restartR()
#library("plot3Drgl")

load("/Volumes/ncatssctl/NGS_related/ddSeq/Analysis_across_runs/R_scripts/Seurat/Aggr6.Seurat.Object.RData")

mycolanno <- data.table::fread('/Volumes/ncatssctl/NGS_related/ddSeq/Analysis_across_runs/R_scripts/ddseq-mycolanno.csv')
mycolanno <- subset(mycolanno, mycolanno$Barcode %in% names(as.data.frame(as.matrix(aggr6@assays$RNA@data))))
names(mycolanno) <- c('cell_ids', 'cell_info')

grouping <- as.data.frame(mycolanno)
row.names(grouping) <- mycolanno$cell_ids
grouping <- subset(grouping, select=c('cell_info'))
head(grouping)
grouping <- as.vector(grouping)
grouping <- setNames(as.character(grouping$cell_info), row.names(grouping))
head(grouping)

# for using full dataset
#aggr6.raw <- as.matrix(aggr6@raw.data)
#aggr6.data <- as.matrix(aggr6@data)

#subsample

# aggr6.subset <- SubsetData(aggr6, ident.use='day0', subset.raw=T)
# aggr6.subset.raw <- as.matrix(aggr6.subset@raw.data)
# aggr6.subset.data <- as.matrix(aggr6.subset@data)
# #
# aggr6.subset.raw <- aggr6.subset.raw[1:10000,]# doesn't work
# aggr6.subset.data <- aggr6.subset.data[1:10000,]
# #
# 
# aggr6.subset.raw <- as.matrix(aggr6@raw.data)[1:200,] # works
# aggr6.subset.data <- as.matrix(aggr6@data)[1:200,]
# 
# 
# subsample.genes <- sample(row.names(as.data.frame(as.matrix(aggr6@raw.data))), 5000, replace=F) # works.
# subsample.genes <- sample(row.names(as.data.frame(as.matrix(aggr6@raw.data))), 7000, replace=F) # works!
# subsample.genes <- sample(row.names(as.data.frame(as.matrix(aggr6@raw.data))), 9000, replace=F) # doesnt work
# subsample.genes <- sample(row.names(as.data.frame(as.matrix(aggr6@raw.data))), 8000, replace=F) # works
# subsample.genes <- sample(row.names(as.data.frame(as.matrix(aggr6@raw.data))), 8500, replace=F) #doesnt work
#subsample.genes <- sample(row.names(as.data.frame(as.matrix(aggr6@raw.data))), 8200, replace=F)  #works
subsample.genes <- aggr6@assays$RNA@var.features

subsample.genes <- unique(c(subsample.genes, 'ESRG', 'TPBG', 'PAX6', 'DUSP6', 'NANOG', 'POU5F1', 'PRTG'))

aggr6.subset.raw <- as.matrix(aggr6@assays$RNA@counts)[subsample.genes,]
aggr6.subset.data <- as.matrix(aggr6@assays$RNA@data)[subsample.genes,]

#
# subsample.genes.sig <- FindMarkers(object = aggr6, ident.1 = "day0", ident.2="day7", min.pct = 0.25)
# subsample.genes.all <- row.names(subsample.genes.sig) # trying with all the genes here even though not all are sig

# aggr6.subset.raw <- as.matrix(aggr6@raw.data)[subsample.genes.all,]
# aggr6.subset.data <- as.matrix(aggr6@data)[subsample.genes.all,]
# #
# subsample.genes.sig <- FindMarkers(object = aggr6, ident.1 = "day7", min.pct = 0.25)
# subsample.genes.all <- row.names(subsample.genes.sig) # trying with all the genes here even though not all are sig
# 
# aggr6.subset.raw <- as.matrix(aggr6@raw.data)[subsample.genes.all,]
# aggr6.subset.data <- as.matrix(aggr6@data)[subsample.genes.all,]
# #

#using aggr6 without G2M or S genes (13715)----
#genes.keep <- row.names(exp)
#aggr6.subset.raw <- as.matrix(aggr6@raw.data)[genes.keep,]
#aggr6.subset.data <- as.matrix(aggr6@data)[genes.keep,]


task <- wrap_expression(
  counts = t(aggr6.subset.raw),
  expression= t(aggr6.subset.data),
  cell_info = mycolanno
)



# # full
# task <- wrap_expression(
#   counts = t(as.matrix(aggr6@raw.data)),
#   expression = t(as.matrix(aggr6@data)),
#   cell_info = mycolanno
# )
# 
# raw <- t(as.matrix(aggr6@raw.data))
# data.norm <- t(as.matrix(aggr6@data))
# 
# 
# task <- wrap_expression(
#   counts = raw,
#   expression = data.norm
# )
# 
# head(mycolanno)
# str(task)

guidelines <- guidelines_shiny(task)
answers <- dynguidelines::answer_questions(
  multiple_disconnected = NULL,
  expect_topology = NULL,
  expected_topology = NULL,
  n_cells = 1385,
  n_features = 13953,
  prior_information = NULL
)
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE, 
  expect_topology = TRUE, 
  expected_topology = "linear", 
  n_cells = 1385, 
  n_features = 3404, 
  memory = "10GB", 
  prior_information = c("start_id", "end_id", "end_n", "start_n", "groups_n", "features_id"), 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers) 

#model <- infer_trajectory(task, 'projected_slingshot', start_id=c(mycolanno[cell_info=="day0",c("cell_ids")]),
#                          end_id=c(mycolanno[cell_info=="day7",c("cell_ids")]))

model.projectedslingshot <- infer_trajectory(task, 'projected_slingshot')
model.slingshot <- infer_trajectory(task, 'slingshot')
#model.projectedslingshot.full <- infer_trajectory(task, 'projected_slingshot')
model.celltree <- infer_trajectory(task, 'celltree_maptpx')
model.scorpius <- infer_trajectory(task, 'scorpius')
model.tscan <- infer_trajectory(task, 'tscan') # works but it labels day 2 as endpoint NPC
#model.projectedpaga <- infer_trajectory(task, 'projected_paga', verbose=T) #errors
model.paga <- infer_trajectory(task, 'paga', verbose=T)

model.monocle <- infer_trajectory(task, 'monocle_ddrtree', verbose=T)
model.mnclinca <- infer_trajectory(task, 'mnclica') #errors
#model.recat <- infer_trajectory(task, 'recat') #errors

#task <- add_prior_information(task, start_id = c(mycolanno[cell_info=="day0",c("cell_ids")]))
#model <- infer_trajectory(task, 'cellrouter', start_id=1)

model.embeddr <- infer_trajectory(task, 'embeddr') # ok, it is linear
model.slice <- infer_trajectory(task, 'slice') # errors
#model.waterfall <- infer_trajectory(task, 'waterfall') # linear though not as straight as embedder

model <- model.embeddr
model <- model.projectedslingshot
model <- model.slingshot
#model <- model.projectedslingshot.full
model <- model.monocle
model <- model.tscan
model <- model.slice
model <- model.projectedpaga

model <- model %>% add_root_using_expression(c("DUSP6"), task$expression)

# model <- label_milestones_markers(
#   model,
#   markers = list(
#     NPC = c("PAX6", "PRTG", "TPBG"),
#     Pluripotent=c('ESRG')
#     ),
#   task$expression
# )# works well for slingshot

# model <- label_milestones_markers(
#   model,
#   markers = list(
#     Neuroprogenitor = c("PAX6", "PRTG", "TPBG"),
#     Pluripotent=c('DUSP6','NANOG', 'ESRG')
#   ),
#   task$expression
# )
model <- label_milestones_markers(model, list(
  Neuroprogenitor = c("PAX6", "PRTG", "TPBG"),
  Pluripotent=c('DUSP6','NANOG', 'ESRG')),
  expression_source = task$expression, n_nearest_cells = 20)

# model <- add_waypoints(model, n_waypoints = 3,
#               resolution = sum(model$milestone_network$length)/3)

# Midpoint = c("TMSB4X", "POU3F1"),

task$grouping <- grouping
model$grouping <- grouping
model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = task$expression)

plot_dimred(
  model, 
  expression_source = task$expression, 
  grouping = grouping,
  groups=grouping,
  color_cells="auto"
)

dynplot::plot_topology(model)

plot_dendro(
  model, 
  expression_source = task$expression,
  grouping = grouping,
  groups=grouping,
  color_cells="auto"
)

plot_dimred(
  model, 
  expression_source = task$expression,
  grouping = grouping,
  color_cells="milestone"
)

plot_dimred(
  model, 
  expression_source = task$expression,
  grouping = grouping,
  color_cells="pseudotime"
)

plot_dimred(
  model,
  expression_source = task$expression,
  feature_oi = "POU5F1"
)

plot_heatmap(
  model,
  expression_source = task$expression,
  grouping = task$grouping,
  features_oi = 50,label_milestones = T, color_cells = "pseudotime"
)


plot_dimred(
  model, 
  expression_source = task$expression, 
  color_cells = "feature",
  feature_oi = "PAX6",
  color_density = "grouping",
  grouping = task$grouping,
  label_milestones = TRUE
)


branching_milestone <- model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()

#branch_feature_importance <- calculate_branching_point_feature_importance(model, expression_source=task$expression, milestones_oi = branching_milestone)
branch_feature_importance <- calculate_branch_feature_importance(model, expression_source=task$expression)

branching_point_features <- branch_feature_importance %>% top_n(20, importance) %>% pull(feature_id)

plot_heatmap(
  model,
  expression_source = task$expression,
  features_oi = branching_point_features, grouping= task$grouping
)

branching_point_features <- branch_feature_importance %>% top_n(20, importance) %>% pull(feature_id)

space <- dyndimred::dimred_mds(task$expression)
map(branching_point_features[1:12], function(feature_oi) {
  plot_dimred(model, dimred = space, expression_source = task$expression, feature_oi = feature_oi, label_milestones = FALSE) +
    theme(legend.position = "none") +
    ggtitle(feature_oi)
}) %>% patchwork::wrap_plots()

overall_feature_importances <- dynfeature::calculate_overall_feature_importance(model.slingshot, expression_source=task$expression)
features <- overall_feature_importances %>% 
  top_n(40, importance) %>% 
  pull(feature_id)

plot_heatmap(
  model,
  expression_source = task$expression,
  features_oi = features, grouping= task$grouping
)
