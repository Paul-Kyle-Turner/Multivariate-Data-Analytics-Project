# Author : Paul Turner
library(MASS)
library(dplyr)
library(class)
library(rgl)
library(Rdimtools)
library(caret)
library(ROCR)
library(GGally)
library(MASS)
library(car)

add.to.each.string = function(string.vec, string){
  for (i in 1:length(string.vec)) {
    string.vec[i] = paste(string.vec[i], string)
  }
  return(string.vec)
}

seed = function(seed = 226262, func = NULL, ...){
  set.seed(seed)
  if(!is.null(func)){
    return(func(...))
  }
}

check.pca = function(tbl, pca = NULL, ...){
  if(is.null(pca)){
    tbl.pca = princomp(tbl, cor = TRUE, ...)
  }else{
    tbl.pca = pca
  }
  return(tbl.pca)
}

check.pca.scores = function(tbl, pca = NULL, ...){
  pca = check.pca(tbl, pca, ...)
  return(pca$scores)
}

determine.k.cv = function(tbl, results, k.min = 1, k.max = 10, seeds = TRUE, confusion.func = NULL, verbose = FALSE, ...){
  save.pred = list()
  k.dist = k.max - k.min
  max.acc = -1 
  best.acc = 1
  
  for (i in 1:k.dist) {
    if(seeds){
      seed()
    }
    pred = knn.cv(tbl, results, k = k.min + i - 1, ...)
    if(is.null(confusion.func)){
      temp = acc.confusion(pred, results)
      if(temp$Acc > max.acc){
        max.acc = temp$Acc
        best.acc = i + k.min - 1
      }
    }else{
      temp = confusion.func(pred, results)
    }
    save.pred[i] = list(temp)
  }
  if(is.null(confusion.func)){
    if(!verbose){
      return.list = list(c(max.acc, round(best.acc)))
      names(return.list) = c("best acc")
      return(return.list)
    }else{
      return.list = list(c(max.acc, round(best.acc)), save.pred)
      names(return.list) = c("best acc","confusion.acc lists")
      return(return.list)
    }
  }else{
    names(save.pred) = c("confusion.acc lists")
    return(save.pred)
  }
}

cv = function(tbl, func, form, ...){
  n = nrow(tbl)
  pred = vector(mode = "numeric")
  for (i in 1:n) {
    mod.fit = func(formula = form, data = tbl[-i,], ...)
    pred[i] = predict(object = mod.fit, newdata = tbl[i,], type = "response")
  }
  return(pred)
}

fold.cv = function(tbl, restults, n.cross = 10, ...){
  cross = nrow(tbl) / n.cross
  preds.list = list()
  for (i in 1:n.cross) {
    s.ind = ((i - 1) * cross) + 1
    e.ind = (i * cross)
    test.tbl = tbl[s.ind:e.ind,]
    train.tbl = tbl[-s.ind:-e.ind,]
    preds.list[i] = knn(train.tbl, test.tbl, ...)
  }
  
  
}

fac.to = function(fac, vec.change = NULL){
  fac.unq = unique(fac)
  fac.to.vec = vector(length = length(fac))
  if(is.null(vec.change)){
    vec.change = rainbow(length(fac.unq))
  }
  for (i in 1:length(fac)) {
    for (k in 1:length(fac.unq)) {
      if(fac[i] == fac.unq[k]){
        fac.to.vec[i] = vec.change[k]
      }
    }
  }
  return(fac.to.vec)
}

true.n = function(confusion){
  n = nrow(confusion)
  true.neg = vector(length = n)
  for (i in 1:nrow(confusion)) {
    temp = confusion[-i,-i]
    true.neg[i] = sum(temp)
  }
  return(true.neg)
}

acc.confusion = function(pred, real){
  confusion = table(pred, real)
  total = sum(confusion)
  confusion.colsum = colSums(confusion)
  confusion.rowsum = rowSums(confusion)
  
  true.pos = diag(confusion)
  true.neg = true.n(confusion)
  false.pos = confusion.colsum - true.pos
  false.neg = confusion.rowsum - true.pos
  
  pos.pred = round(true.pos / (true.pos + false.neg), 4) 
  neg.pred = round(true.neg / (true.neg + false.pos), 4)
  
  sens = round(true.pos / (true.pos + false.pos), 4)
  spec = round(true.neg / (true.neg + false.neg), 4)
  
  acc = round(sum(true.pos) / total, 4)
  
  return.list = list(confusion, acc, sens, spec, pos.pred, neg.pred)
  names(return.list) = c("Confusion Matrix", "Acc", "Sens", "Spec","Pos.pred", "Neg.pred")
  
  return(return.list)
}

acc.confusion.vec = function(pred, real, uniq = NULL,  percent = FALSE){
  if(!is.null(uniq)){
    pred.unq = uniq
  }else{
    pred.unq = unique(pred)
  }
  pred.prop = rep(0, length(pred.unq) + 1)
  totals = rep(0, length(pred.unq))
  total = length(pred)
  for (i in 1:total) {
    if(pred[i] == real[i]){
      for (k in 1:length(pred.unq)) {
        if(is.na(pred[k])){
          next()
        }
        else if(pred[i] == pred.unq[k]){
          pred.prop[k] = pred.prop[k] + 1
          totals[k] = totals[k] + 1
        }
      }
      pred.prop[length(pred.prop)] = pred.prop[length(pred.prop)] + 1
    }else{
      for (k in 1:length(pred.unq)) {
        if(real[i] == pred.unq[k]){
          totals[k] = totals[k] + 1
        }
      }
    }
  }
  totals = c(totals, total)
  if(percent){
    return(c(round(pred.prop / totals, 4), total))
  }else{
    return(c(pred.prop, total))
  }
}

plot.color.2tbl = function(tbl1, tbl2, col1 = "green", col2 = "orange"){
  return(c(rep(col1, nrow(tbl1)), rep(col2, nrow(tbl2))))
}

plot.pca.bubble = function(tbl, r.color, pca = NULL, scale.cirlces = 25, havetext = FALSE, ...){
  tbl.pca.scores = check.pca.scores(tbl, pca)
  tbl.pca.pos.col = ifelse(tbl.pca.scores[,3] > 0, "red", "blue")
  tbl.pca.pos = (tbl.pca.scores[,3] - min(tbl.pca.scores[,3])) / scale.cirlces # division to make readable
  symbols(tbl.pca.scores[,1], tbl.pca.scores[,2], circles = tbl.pca.pos,
          inches = FALSE, fg = tbl.pca.pos.col, ...)
  if(havetext){
    text(x = tbl.pca.scores[,1], y = tbl.pca.scores[,2], col = r.color)
  }else{
    points(x = tbl.pca.scores[,1], y = tbl.pca.scores[,2], col = r.color, pch = 15)
  }
}

plot.pca.3d = function(tbl, r.color, pca = NULL, ...){
  tbl.pca.scores = check.pca.scores(tbl, pca)
  plot3d(tbl.pca.scores[,1], tbl.pca.scores[,2], tbl.pca.scores[,3],
         col = r.color, ...)
}

plot.pca.scree = function(tbl = NULL, pca = NULL, ...){
  tbl.pca = check.pca(tbl, pca)
  prop.pc.var = (tbl.pca$sdev ^ 2) / sum(tbl.pca$sdev ^ 2)
  plot(cumsum(prop.pc.var), ...)
}

plot.samples = function(tbl1, tbl2, plot.func, r.col = TRUE, ...){
  if(r.col){
    tbl.r.color = c(rep("green", nrow(tbl1)), rep("orange", nrow(tbl2)))
    plot.func(rbind(tbl1, tbl2), r.color = tbl.r.color, ...)
  }else{
    plot.func(rbind(tbl1, tbl2), ...)
  }
}

plot.square = function(k = NULL, starting = 10, ending = 90, by = 5){
  if(is.null(k)){
    k = (ending - starting) / by
  }
  sq.k = sqrt(k)
  if(sq.k%%1==0){
    par(mfrow = c(sq.k, sq.k))
  }else{
    par(mfrow = c(1,sq.k))
  }
}

plot.box.sne.pcar = function(tbl, r.color, sne.function, starting = 0, ending = 9, by = .1, pca = TRUE, ...){
  k = (ending - starting)
  plot.square(k = k)
  sne.list = list()
  for (i in 1:k) {
    pcar = (i - 1) * by
    if(pcar == 0){
      sne = sne.function(tbl, pca = FALSE, ...)
    }else{
      sne = sne.function(tbl, pca = TRUE, ...)
    }
    if(all(apply(sne$Y, 2, is.finite))){
      sne.list[i] = list(sne)
      plot(sne$Y, col = r.color, main = paste("PC Variance cutoff ", pcar))
    }else{
      i = i - 1
    }
  }
  return(sne.list)
}

plot.box.sne.perp = function(tbl, r.color, sne.function, starting = 5, ending = 50, by = 5, pca = TRUE, pcaratio = .9, ...){
  k = (ending - starting) / by
  plot.square(starting = starting, ending = ending, by = by)
  sne.list = list()
  for (i in 1:k) {
    perp = starting + (i * by)
    sne = sne.function(tbl, perplexity = perp, pca = pca, pcaratio = pcaratio, ...)
    if(all(apply(sne$Y, 2, is.finite))){
      sne.list[i] = list(sne)
      plot(sne$Y, col = r.color, main = paste("Perplexity ", perp))
    }else{
      i = i - 1
    }
  }
  return(sne.list)
}

plot.box.sne.lr = function(tbl, r.color, sne.function, lr = c(.01, .0175, .025, .0375, .05, .0675, .075, .0825, .09), ...){
  k = 9
  plot.square(k = 9)
  sne.list = list()
  for (i in 1:k) {
    t.lr = lr[i]
    sne = sne.function(tbl, pca = TRUE, eta = t.lr, ...)
    if(all(apply(sne$Y, 2, is.finite))){
      sne.list[i] = list(sne)
      plot(sne$Y, col = r.color, main = paste("Learning rate ", t.lr))
    }else{
      i = i - 1
    }
  }
  return(sne.list)
}

plot.create.samples = function(tbl, results = NULL, plot.func, percent = .3, ...){
  if(is.null(results)){
    results = tbl[,1]
  }
  result.types = unique(results)
  tbl1 = sample.tbl.return(tbl[results == result.types[1],], percent = percent)
  tbl2 = sample.tbl.return(tbl[results == result.types[2],], percent = percent)
  plot.samples(tbl1, tbl2, plot.func, ...)
  return(rbind(tbl1, tbl2))
}

plot.dend = function(tbl, distance.method = "euclidean", nclust = 2, clust.method.vec = NULL, ...){
  if(is.null(clust.method.vec)){
    clust.method.vec = c("ward.D2", "single", "complete", "average")
  }
  k = length(clust.method.vec)
  par(mfrow = c(1,k))
  for (i in 1:k) {
    tbl.dist = dist(tbl, method = distance.method)
    clust.m = hclust(tbl.dist, method = clust.method.vec[i])
    plot(clust.m)
    rect.hclust(clust.m, k = nclust)
  }
}

plot.samples = function(tbl1, tbl2, plot.func, r.col = TRUE, ...){
  if(r.col){
    tbl.r.color = plot.color.2tbl(tbl1, tbl2)
    plot.func(rbind(tbl1, tbl2), r.color = tbl.r.color, ...)
  }else{
    plot.func(rbind(tbl1, tbl2), ...)
  }
}

plot.sne.3d = function(tbl, r.color, sne.function, perplexity = 30, ...){
  sne = sne.function(tbl, ndim = 3, perplexity = perplexity, ...)
  plot3d(sne$Y, col = r.color)
}

sample.tbl = function(tbl, n, ...){
  # just bananas
  seed()
  return(sample(nrow(tbl), n, ...))
}

sample.tbl.percent = function(tbl, percent = .3, ...){
  n = nrow(tbl) * percent
  return(sample.tbl(tbl, n, ...))
}

sample.tbl.return = function(tbl, percent = .3, ...){
  return(tbl[sample.tbl.percent(tbl, percent, ...),])
}

sample.tbl.func = function(tbl, func, percent = .3, ...){
  return(func(sample.tbl.return(tbl, percent), ...))
}

sample.tbl.col.func = function(tbl, col, func, percent = .3, proportional = TRUE, r.color = FALSE, ...){
  if(proportional){
    
    tbl.t = data.frame()
    col.unq = unique(col)
    col.t = vector()
    
    for (i in 1:length(col.unq)) {
      temp.t = as.data.frame(sample.tbl.return(tbl[col == col.unq[i],], percent))
      temp.t = na.omit(temp.t)
      col.t = c(col.t, rep(col.unq[i], nrow(temp.t)))
      tbl.t = rbind(tbl.t, temp.t)
    }
    
    if(r.color){
      return(func(tbl.t, r.color = col.t, ...))
    }else{
      return(func(tbl.t, col = col.t, ...))
    }
    
  }else{
    idx = sample.tbl.percent(tbl, percent)
    tbl.t = tbl[idx,]
    col.t = col[idx]
    if(r.color){
      return(func(tbl.t, r.color = col.t, ...))
    }else{
      return(func(tbl.t, col = col.t, ...))
    }
  }
}









