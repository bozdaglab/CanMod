library(pbapply)
library(HDCI)
library(data.table)
library(bnlearn)
library(doParallel)
library(parallel)

# Load utility functions to run CanMod ----------------------------------------------------
source("helper_functions.R")

### --- Global variables ---- #####
args = commandArgs(trailingOnly = T)
if (!exists("cancer.type")){
  cancer.type = ifelse(length(args) !=0 , args[1], "BRCA")
}

# Load required input for CanMod  ---------------------------------------------------------
load(paste(cancer.type, "variable.selection.data.rda", sep = "/" ))

load("miRNA_mRNA_interactions.rda")
load("TF.target.rda")

putative.miRNA.mRNA.interactions = putative.miRNA.mRNA.interactions
colnames(putative.miRNA.mRNA.interactions) = c("miRNA","target")

# set up regression.data
regression.data = variable.selection.data
regression.data$interactions = NULL

{
  samples = colnames(mRNA)[grep(colnames(mRNA), pattern = "Tumor")]
  regression.data$mRNA = mRNA[,samples]
  regression.data$miRNA = miRNA[,samples]
  regression.data$methyl = regression.data$methyl[,samples]
  regression.data$cna = regression.data$cna[,samples]
  regression.data$miRNA.target.interactions = miRNA.target.interactions
  regression.data$miRNA.target.interactions = regression.data$miRNA.target.interactions[regression.data$miRNA.target.interactions$miRNA %in% rownames(regression.data$miRNA),]
  regression.data$miRNA.target.interactions = regression.data$miRNA.target.interactions[regression.data$miRNA.target.interactions$target %in% rownames(regression.data$mRNA),]
  regression.data$tf.target.interactions = TF.target
  regression.data$tf.target.interactions = regression.data$tf.target.interactions[regression.data$tf.target.interactions$tf %in% rownames(regression.data$mRNA),]
  regression.data$tf.target.interactions = regression.data$tf.target.interactions[regression.data$tf.target.interactions$target %in% rownames(regression.data$mRNA),]
  rownames(regression.data$miRNA) = gsub(rownames(regression.data$miRNA), pattern = "[-]", replacement = ".")
  regression.data$miRNA.target.interactions$miRNA = gsub(regression.data$miRNA.target.interactions$miRNA, pattern = "[-]", replacement = ".")
}

# STEP 1: variable selection -----------------------------------
mRNA.targets = rownames(regression.data$mRNA)

cluster = parallel::makeCluster(parallel::detectCores()-1, outfile = "log/april14.out")
parallel::clusterExport(cl = cluster, 
                        varlist = c("mRNA.targets","regression.data", "getDF"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table, HDCI)))
coefs = pbapply::pblapply(cl = cluster, X = 1:length(mRNA.targets), FUN = function(mRNA.index){
  mRNA = mRNA.targets[mRNA.index]
  cat("Start", mRNA, "\n")
  include.lncRNA = grepl(mRNA, pattern = "lncRNA")
  regression.df = getDF(regression.data = regression.data, mRNA = mRNA,
                        onlymiR = F, include.lncRNA = include.lncRNA)
  if (is.null(regression.df)){
    return(NULL)
  }
  
  y.matrix = as.matrix(regression.df$mRNA)
  x.matrix = as.matrix(regression.df[,-1])
  colnames(x.matrix) = colnames(regression.df)[-1]
  candidate.regulators = colnames(x.matrix)
  # if only 1 miRNA
  if(ncol(x.matrix) == 1){
    if(cor(x.matrix, y.matrix) > 0){
      return(NULL)
    }else{
      if (cor.test(x.matrix, y.matrix)$p.value >= 0.05){
        return(NULL)
      }else{
        normal.lasso.dt = data.frame(regulator=candidate.regulators, target=mRNA, 
                                     coef = cor(x.matrix, y.matrix), 
                                     stringsAsFactors = F)
        normal.lasso.dt = as.data.table(normal.lasso.dt)
        normal.lasso.dt$pair = paste(normal.lasso.dt$regulator, normal.lasso.dt$target,sep="~")
        return(normal.lasso.dt)
      }
    }
  }
  
  # get selected coefs
  selected.coefs = lapply(1:100, function(iter){
    #cat(iter, "\n")
    y.matrix = as.matrix(regression.df$mRNA)
    x.matrix = as.matrix(regression.df[,-1])
    colnames(x.matrix) = colnames(regression.df)[-1]
    
    lasso.result = tryCatch({
      HDCI::Lasso(x = x.matrix, y = y.matrix,
                  fix.lambda = F,nfolds = 10,
                  cv.method = "cv1se",
                  standardize = T, parallel = F)
    }, warning = function(w) {
      NULL
    }, error = function(e) {
      NULL
    })
    
    if (is.null(lasso.result)){
      return(NULL)
    }
    
    coefs = lasso.result$beta
    names(coefs) =  colnames(x.matrix) 
    coefs = coefs[which(coefs != 0)]
    return(coefs)
  })
  
  # make consistently selected data strucutre
  normal.lasso.dt = NULL
  dt.list = lapply(1:length(selected.coefs), function(target.index){
    regulators = selected.coefs[[target.index]]
    if (is.null(regulators)){
      return(NULL)
    }
    target = mRNA
    regulator.names = names(regulators)
    regulator.coef =unname(regulators)
    if (length(regulator.names) > 0){
      df = data.frame(regulator=regulator.names,
                      target = rep(target,length(regulator.names)),
                      coef=as.vector(regulator.coef),
                      stringsAsFactors = F)
      df$pair = paste(df$regulator,df$target,sep="~")
      dt = as.data.table(df)
      return(dt)
    }else{
      return(NULL)
    }
  })
  if (is.null(dt.list)){
    normal.lasso.dt = NULL
    return(normal.lasso.dt)
  }
  
  normal.lasso.dt = tryCatch({
    rbindlist(dt.list)
  }, warning = function(w) {
    NULL
  }, error = function(e) {
    NULL
  })
  
  cat("Done", mRNA, "\n")
  return(normal.lasso.dt)
})
names(coefs) = mRNA.targets
save(coefs, file =  "codes/test/test_final/brca.coefs.rda")

#get bootstrap interval
bt.interval.list = pbapply::pblapply(1:length(mRNA.targets),FUN = function(mRNA.index){
  if (mRNA.index %% 100 == 0) cat(mRNA.index, "\n")
  mRNA = mRNA.targets[mRNA.index]
  include.lncRNA = grepl(mRNA, pattern = "lncRNA")
  regression.df = getDF(regression.data = regression.data, mRNA = mRNA,
                        onlymiR = F, include.lncRNA = F)
  
  if (is.null(regression.df)){
    return(NULL)
  }
  
  y.matrix = as.matrix(regression.df$mRNA)
  x.matrix = as.matrix(regression.df[,-1])
  colnames(x.matrix) = colnames(regression.df)[-1]
  
  # if only 1 miRNA
  if(ncol(x.matrix) == 1){
    if(cor(x.matrix, y.matrix) > 0){
      return(NULL)
    }else{
      if (cor.test(x.matrix, y.matrix)$p.value >= 0.05){
        return(NULL)
      }else{
        bt.interval = matrix(c(-1000,1000), ncol = 1)
        colnames(bt.interval) = colnames(x.matrix)
        return(bt.interval)
      }
    }
  }
  
  num.core = parallel::detectCores()
  doParallel::registerDoParallel(num.core)
  num.bt.replications = 100
  btlasso.result = tryCatch({
    bootLasso(x = x.matrix, y = y.matrix,
              B = num.bt.replications,
              cv.method = "cv1se",
              type.boot = "paired",
              standardize = T,  parallel.boot = T, ncores.boot = num.core)
  }, warning = function(w) {
    NULL
  }, error = function(e) {
    NULL
  })
  
  if (is.null(btlasso.result)){
    return(NULL)
  }
  
  bt.interval = btlasso.result$interval
  colnames(bt.interval) = colnames(x.matrix)
  return(bt.interval)
})

names(bt.interval.list) = mRNA.targets

regulator.list = lapply(1:length(coefs), function(index){
  gene.name = mRNA.targets[index]
  coef.dt = coefs[[gene.name]]
  bt.interval.dt = bt.interval.list[[gene.name]]
  list(coef.dt=coef.dt, bt.interval.dt=bt.interval.dt)
})
names(regulator.list) = mRNA.targets

regulator.target.pair.list = lapply(1:length(regulator.list), function(gene.index) {
  if(gene.index %% 100 == 0) print(gene.index)
  #print(gene.index)
  gene.pair.info.list = regulator.list[[gene.index]]
  gene.coef.dt = gene.pair.info.list$coef.dt
  if (nrow(gene.coef.dt) == 0 || is.null(gene.coef.dt)) {
    return(NULL)
  }
  gene.coef.stat.dt = gene.coef.dt[, list(count = .N, 
                                          median.coef = median(.SD$coef)),
                                   by = "pair"]
  gene.coef.stat.dt$regulator = sapply(strsplit(gene.coef.stat.dt$pair, split = "[~]"), function(v)
    v[1])
  gene.coef.stat.dt$target = sapply(strsplit(gene.coef.stat.dt$pair, split = "[~]"), function(v)
    v[2])
  gene.coef.stat.dt = gene.coef.stat.dt[, c(4, 5, 1, 2, 3), ]
  gene.bt.interval.dt = tryCatch({
    gene.pair.info.list$bt.interval.dt
  }, warning = function(w) {
    NULL
  }, error = function(e) {
    NULL
  })
  
  if (is.null(gene.bt.interval.dt) || nrow(gene.bt.interval.dt) == 0){
    return(NULL)
  }
  
  lower.percentile = sapply(1:nrow(gene.coef.stat.dt), function(row.index) {
    regulator = gene.coef.stat.dt$regulator[row.index]
    return(gene.bt.interval.dt[1, regulator])
  })
  gene.coef.stat.dt$lower.percentile = lower.percentile
  
  upper.percentile = sapply(1:nrow(gene.coef.stat.dt), function(row.index) {
    regulator = gene.coef.stat.dt$regulator[row.index]
    return(gene.bt.interval.dt[2, regulator])
  })
  gene.coef.stat.dt$upper.percentile = upper.percentile
  
  gene.coef.stat.dt$confidence =  (gene.coef.stat.dt$median.coef >= lower.percentile) &
    (gene.coef.stat.dt$median.coef <= upper.percentile)
  
  return(gene.coef.stat.dt)
})

regulator.target.pair.dt = rbindlist(regulator.target.pair.list)

# STEP 2: cluster regulators based on shared targets similarity -------------------------------------------------------------------
# load lasso result 
regulator.target.pair.dt = regulator.target.pair.dt[which(regulator.target.pair.dt$confidence == T & regulator.target.pair.dt$count >=75)]
regulator.target.pair.dt = regulator.target.pair.dt[-which(regulator.target.pair.dt$regulator %in% c("CNA","Methyl"))]
to.removed.indicies = which(grepl(regulator.target.pair.dt$regulator, pattern = "hsa") & regulator.target.pair.dt$median.coef >=0 )  # miRNA with positive coefficients
regulator.target.pair.dt = regulator.target.pair.dt[-to.removed.indicies]
length(unique(regulator.target.pair.dt$target))
lasso.df = regulator.target.pair.dt[,c(1,2)]

length(unique(lasso.df$target))
{
  mir.lasso.df = lasso.df[grepl(lasso.df$regulator, pattern = "hsa."),]
  tf.lasso.df = lasso.df[!grepl(lasso.df$regulator, pattern = "hsa."),]
  nrow(mir.lasso.df)
  length(unique(mir.lasso.df$regulator))
  length(unique(mir.lasso.df$target))
  nrow(tf.lasso.df)
  length(unique(tf.lasso.df$regulator))
  length(unique(tf.lasso.df$target))
  
  lasso.regulator.list = lapply(unique(lasso.df$target), function(target){
    regulators= lasso.df$regulator[lasso.df$target == target]
  })
  lasso.tf.list = sapply(lasso.regulator.list, function(reg){
    length(reg[!grepl(reg, pattern = "hsa.")])
  })
  lasso.mir.list = sapply(lasso.regulator.list, function(reg){
    length(reg[grepl(reg, pattern = "hsa.")])
  })
  summary(lasso.tf.list)
  summary(lasso.mir.list)
}

regulation.target.df = matrix(data = 0, 
                              nrow = length(unique(lasso.df$target)),
                              ncol = length(unique(lasso.df$regulator)),
                              dimnames = list(unique(lasso.df$target),
                                              unique(lasso.df$regulator)))
regulation.target.df = as.data.frame(regulation.target.df)

dim(regulation.target.df)
# populate regulation matrix 
for (row.index in 1:nrow(regulation.target.df)){
  if (row.index %% 100 == 0) print(row.index)
  target = rownames(regulation.target.df)[row.index]
  regulators = lasso.df$regulator[which(lasso.df$target == target)]
  regulation.target.df[target,regulators] = 1
}

# cluster regulators based on shared targets similarity
{
  df = t(regulation.target.df)
  d <- dist(df, method = "binary")
  hc <- hclust(d)
  plot(hc)
  # labels at the same level
  # plot(hc, hang = -1)
  regulator.cluster.list = cutree(hc, h = max(hc$height[which(hc$height < 1)]))
}

regulator.cluster.df = data.frame(regulator=names(regulator.cluster.list), cluster = unname(regulator.cluster.list),stringsAsFactors = F)
summary(as.vector(table(regulator.cluster.df$cluster)))

regulator.cluster.list = lapply(unique(regulator.cluster.df$cluster), function(cluster){
  regulator.cluster.df$regulator[regulator.cluster.df$cluster==cluster]
})

# remove regulator clusters with only 1 element
regulator.cluster.list = regulator.cluster.list[-which(sapply(regulator.cluster.list, function(cluster) length(cluster) == 1 ))]
summary(as.vector(sapply(regulator.cluster.list, function(cluster) length(unique(cluster)))))

# STEP 3: Get GO-based cluster --------------------------------------------------------------------------
de.genes = rownames(variable.selection.data$mRNA); length(de.genes)

# get similarity matrix
hsGO2 = godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)
sum(de.genes %in% unique(hsGO2@keys))

cat("Start computing GO-based similarity matrix \n")
time = proc.time()
de.gene.bp.sim = mgeneSim(de.genes, semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)
time = proc.time() - time
cat("running time for computing gene similarity based on GO:", time[3])

cat("Done with obtaining similarity matrix in cluster_de_genes for", cancer.type, "\n") 

# cluster de gene by GO similarity 
time = proc.time()
de.gene.bp.cluster = apcluster(s = de.gene.bp.sim, details = T, maxits = 10000, convits = 100)
time = proc.time() - time

de.gene.bp.cluster.list = lapply(1:length(de.gene.bp.cluster@clusters), function(index){
  cluster = names(de.gene.bp.cluster@clusters[[index]])
  return(cluster)
})

target.sim = de.gene.bp.sim
target.cluster.list = de.gene.bp.cluster.list

names(target.cluster.list) = as.character(1:length(target.cluster.list))
require(qdapTools)
target.cluster.df = list2df(target.cluster.list); names(target.cluster.df) = c("target","cluster")
summary(as.vector(table(target.cluster.df$cluster)))

# STEP 4: Generate candidate modules ------------------------------------------------------------------------
# for each regulator cluster, recruit targets
regulator.target.cluster.list = lapply(1:length(regulator.cluster.list), function(index){
  regulators = regulator.cluster.list[[index]]
  regulator.target.dt = lasso.df[lasso.df$regulator %in% regulators]
  regulator.target.dt = regulator.target.dt[regulator.target.dt$target %in% unique(target.cluster.df$target)]
  regulator.target.dt$cluster = as.numeric(target.cluster.df$cluster[match(regulator.target.dt$target, target.cluster.df$target)])
  #table(regulator.target.dt$regulator)
  #table(regulator.target.dt$target)
  return(regulator.target.dt)
}) 
summary((sapply(regulator.target.cluster.list, function(module) length(unique(module$target)))))
summary((sapply(regulator.target.cluster.list, function(module) length(unique(module$cluster)))))


# for each regulator cluster, create list of GO-based clusters, whose elments are seed gene 
seed.target.list = lapply(1:length(regulator.target.cluster.list), function(k){
  regulator.cluster =  regulator.target.cluster.list[[k]]
  target.clusters = unique(regulator.cluster$cluster)
  target.clusters = lapply(1:length(target.clusters), function(i){
    target.cluster.index = target.clusters[[i]]
    seed.targets = unique(regulator.cluster$target[regulator.cluster$cluster == target.cluster.index] )
    all.targets.in.cluster = unique(target.cluster.df$target[target.cluster.df$cluster == target.cluster.index]  )
    all.targets.in.cluster = setdiff(all.targets.in.cluster, seed.targets)
    return(list(seed.targets= seed.targets, all.targets.in.cluster= all.targets.in.cluster))
  })
  return(target.clusters)
})

cor.threshold = quantile(expression.values, 0.9)

selected.seed.target.list =  lapply(1:length(seed.target.list), function(index){
  overall.seed.list =  seed.target.list[[index]]
  seed.partner.list = lapply(1:length(overall.seed.list), function(i){
    seed.list = overall.seed.list[[i]]
    seed.targets = seed.list$seed.targets
    seed.partners = seed.list$all.targets.in.cluster
    seed.target.partner.list = lapply(1:length(seed.targets), function(k){
      seed.target = seed.targets[k]
      seed.target.cor = expression.cor[seed.target, seed.partners]
      if (length(seed.target) == 1 && abs(seed.target.cor) > cor.threshold){
        selected.partners = seed.target
        return(selected.partners)
      }
      selected.partners = names(which(abs(seed.target.cor) > cor.threshold))
      return(selected.partners)
    })
    names(seed.target.partner.list) = seed.targets
    return(seed.target.partner.list)
  })
})


selected.seed.target.list =  lapply(1:length(selected.seed.target.list), function(index){
  seed.list = selected.seed.target.list[[index]] 
  seed.list = lapply(1:length(seed.list), function(k){
    seed.group = seed.list[[k]] # all belong to a same GO-based cluster
    seed = names(seed.group)
    seed.group = unique(c(seed, unname(unlist(seed.group))))
  })
  names(seed.list) = sapply(seed.list, function(g) g[1])
  return(seed.list)
})

regulator.target.list =  lapply(1:length(regulator.cluster.list), function(index){
  regulators = regulator.cluster.list[[index]]
  target.list = selected.seed.target.list[[index]]
  regulator.target.list = lapply(target.list, function(targets){
    list(regulators=regulators, targets=targets)
  })
  return(regulator.target.list)
})

count = 0
simplified.regulator.target.list = list()
for (i in 1:length(regulator.target.list)){
  for (j in 1:length(regulator.target.list[[i]])){
    simplified.regulator.target.list = rlist::list.append(simplified.regulator.target.list, 
                                                          regulator.target.list[[i]][[j]])
  }
}

summary(sapply(simplified.regulator.target.list, function(module) length(module$regulators)))
summary(sapply(simplified.regulator.target.list, function(module) length(module$targets)))
sum(sapply(simplified.regulator.target.list, function(module) length(module$targets)) == 1)


# Step 5: refine module using SCCA ----------------------------------------------------------------------------------------------
module.expression.list = lapply(1:length(simplified.regulator.target.list), function(index){
  module = simplified.regulator.target.list[[index]]
  if (length(module$targets) == 1){
    target.df = t(as.data.frame(expression.df[module$targets,]))
    rownames(target.df) = module$targets
  }else{
    target.df = expression.df[module$targets,]
  }
  if (length(module$regulators) == 1){
    regulator.df = t(as.data.frame(expression.df[module$regulators,]))
    rownames(regulator.df) = module$regulators
  }else{
    regulator.df = expression.df[module$regulators,]
  }
  module = list(target.df=target.df,
                regulator.df=regulator.df)
})


cluster = parallel::makeCluster(detectCores()-1)
parallel::clusterExport(cl = cluster,
                        varlist = c("module.expression.list"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(PMA)))

cca.module.list = pbapply::pblapply(cl = cluster, X = 1:length(module.expression.list), FUN = function(index){
  cat("Module: " , index, "\n")
  module.expression = module.expression.list[[index]]  
  
  target.df = module.expression$target.df
  regulator.df = module.expression$regulator.df
  
  if (nrow(target.df) == 1 || nrow(regulator.df) == 1){
    cca.module = list(selected.targets=rownames(target.df), selected.regulators=rownames(regulator.df))
    return(cca.module)
  }
  
  perm.out = PMA::CCA.permute(t(target.df), t(regulator.df), nperms = 10)
  cca.result = PMA::CCA(t(target.df), t(regulator.df), 
                        penaltyx=perm.out$bestpenaltyx,
                        penaltyz=perm.out$bestpenaltyz)
  
  selected.targets = rownames(target.df)[which(cca.result$u[,1] != 0)]
  selected.regulators = rownames(regulator.df)[which(cca.result$v[,1] != 0)]
  cca.module = list(selected.targets=selected.targets, 
                    selected.regulators=selected.regulators,
                    cca.cor = cca.result$cors)
  return(cca.module)
}); stopCluster(cluster)

length(cca.module.list)

cca.target.size = sapply(cca.module.list, function(module){
  return(length(module$selected.targets))
})
cca.regulator.size = sapply(cca.module.list, function(module){
  return(length(module$selected.regulators))
})
# remove module whose 
plot(density(cca.target.size)); summary(cca.target.size)
plot(density(cca.regulator.size)); summary(cca.regulator.size)
to.removed.indices =  which(cca.target.size == 1)

if (length(to.removed.indices) > 0){
  filtered.cca.module.list = cca.module.list[-to.removed.indices]  
}

summary(sapply(filtered.cca.module.list, function(module) length(module$selected.regulators)))
summary(sapply(filtered.cca.module.list, function(module) length(module$selected.targets)))

filtered.cca.target.list = lapply(filtered.cca.module.list, function(module){
  module$selected.targets
})
filtered.cca.target.df = qdapTools::list2df(filtered.cca.target.list)
colnames(filtered.cca.target.df) = c("target","module")
nrow(filtered.cca.target.df)

target.module.df = matrix(data = 0, 
                          nrow = length(unique(filtered.cca.target.df$target)),
                          ncol = length(unique(filtered.cca.target.df$module)),
                          dimnames = list(unique(filtered.cca.target.df$target),
                                          unique(filtered.cca.target.df$module)))

target.module.df = as.data.frame(target.module.df)
for (row.index in 1:nrow(target.module.df)){
  if (row.index %% 100 == 0) print(row.index)
  target = rownames(target.module.df)[row.index]
  modules = filtered.cca.target.df$module[which(filtered.cca.target.df$target == target)]
  target.module.df[target,modules] = 1
}

length(which(target.module.df != 0))

save(lasso.df, regulator.cluster.list, target.cluster.list, simplified.regulator.target.list, filtered.cca.module.list, file = "code/code_test/explore_cluster/v1.results.rda")

# Step 6: cluster modules by shared target -----------------------------------------------------
{
  df = target.module.df
  d <- dist(t(df), method = "binary")
  module.distance = d 
  hc <- hclust(module.distance)
  plot(hc)
  hist(hc$height)
  module.cluster.list = cutree(hc, h = max(hc$height[which(hc$height < 1)]))
  table(module.cluster.list)
  
  #module.cluster.list = cutree(hc, h = 1)
  #length(module.cluster.list)
  
  
  module.cluster.df = data.frame(module=names(module.cluster.list), cluster = unname(module.cluster.list),stringsAsFactors = F)
  summary(as.vector(table(regulator.cluster.df$cluster)))
  
  module.cluster.list = lapply(unique(module.cluster.df$cluster), function(cluster){
    module.cluster.df$module[module.cluster.df$cluster==cluster]
  })
}
# for each cluster modules, recruit back targets 
final.module.list = lapply(1:length(module.cluster.list), function(index){
  modules = as.numeric(module.cluster.list[[index]])
  detailed.modules = filtered.cca.module.list[modules]
  combined.regulators = lapply(detailed.modules, function(module){
    module$selected.regulators
  })
  combined.regulators = unique(unlist(combined.regulators))
  combined.targets = lapply(detailed.modules, function(module){
    module$selected.targets
  }) 
  combined.targets = unique(unlist(combined.targets))
  return(list(regulators = combined.regulators, targets = combined.targets))
})

length(final.module.list)

regulator.size = sapply(final.module.list, function(module){
  length(module$regulators)
})
summary(regulator.size)

target.size = sapply(final.module.list, function(module){
  length(module$targets)
})
summary(target.size)


{
  hc <- hclust(module.distance)
  plot(hc)
  hist(hc$height)
  module.cluster.list = cutree(hc, h = max(hc$height[which(hc$height < 1)]))
  table(module.cluster.list)
  sum(table(module.cluster.list))
  
  module.cluster.list = cutree(hc, h = 1)
  table(module.cluster.list)
  sum(table(module.cluster.list))
  
  dtc.module.cluster.list = dynamicTreeCut::cutreeDynamicTree(hc, minModuleSize= 1)
  table(dtc.module.cluster.list)
  names(dtc.module.cluster.list) = 1:length(dtc.module.cluster.list)
  
  dtc.module.cluster.df = data.frame(module=names(dtc.module.cluster.list), cluster = unname(dtc.module.cluster.list),stringsAsFactors = F)
  summary(as.vector(table(dtc.module.cluster.df$cluster)))
  
  dtc.module.cluster.list = lapply(unique(dtc.module.cluster.df$cluster), function(cluster){
    dtc.module.cluster.df$module[dtc.module.cluster.df$cluster==cluster]
  })
  
  # for each cluster modules, recruit back targets 
  dtc.final.module.list = lapply(1:length(dtc.module.cluster.list), function(index){
    modules = as.numeric(dtc.module.cluster.list[[index]])
    dtc.detailed.modules = filtered.cca.module.list[modules]
    combined.regulators = lapply(dtc.detailed.modules, function(module){
      module$selected.regulators
    })
    combined.regulators = unique(unlist(combined.regulators))
    combined.targets = lapply(dtc.detailed.modules, function(module){
      module$selected.targets
    }) 
    combined.targets = unique(unlist(combined.targets))
    return(list(regulators = combined.regulators, targets = combined.targets))
  })
}
