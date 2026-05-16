
## ReScale function -----
ReScale <- function(x, first, last) {
  (last - first) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (x - min(x, na.rm = TRUE)) + first
}

# partial Spearman correlation function -----
partialSpearmanJSz <- function(u, v, lowertail, uppertail, bounds, min.overlap, tail){
  
  if(!is.numeric(u) | !is.numeric(v)){
    stop("u and v must be numeric")
  }
  if(length(u)!=length(v)){
    stop("u and v have different length")
  }
  if(min(bounds)<0 | max(bounds)>1){
    stop("bounds must be between 0 and 1")
  }
  
  ## min. overlap filter and NAs removal
  tmp <- data.frame(u, v)
  tmp <- na.omit(tmp)
  tmp$year <- as.numeric(rownames(tmp))
  
  n <- length(tmp$u)
  
  tmp <- as.matrix(tmp)
  
  ## if there is less pairwise obs than min. required, return NA for synchrony
  if(n < min.overlap) {
    r <- NA
    
    return(r)
  }
  
  ## else, calculate ranks and split the observations into tails
  else {
    x<-tmp[,1]
    y<-tmp[,2]
    
    x.ranked <- rank(x)/(n+1)
    y.ranked <- rank(y)/(n+1)
    
    ranks <- data.frame(x.ranked = x.ranked, y.ranked = y.ranked,
                        year = rownames(data.frame(x.ranked)))
    
    x_b_lb <- x[x >= (min(lowertail)) & x <= (max(lowertail))]
    y_b_lb <- y[y >= (min(lowertail)) & y <= (max(lowertail))]
    
    x_b_ub <- x[x > (min(uppertail)) & x <= (max(uppertail))]
    y_b_ub <- y[y > (min(uppertail)) & y <= (max(uppertail))]
    
    tmp.x.lb <- data.frame(x = x_b_lb,
                           year = rownames(data.frame(x_b_lb)))
    tmp.y.lb <- data.frame(y = y_b_lb,
                           year = rownames(data.frame(y_b_lb)))
    
    tmp.x.ub <- data.frame(x = x_b_ub,
                           year = rownames(data.frame(x_b_ub)))
    tmp.y.ub <- data.frame(y = y_b_ub,
                           year = rownames(data.frame(y_b_ub)))

    # ## years for lb and ub
    years.lb <- c(tmp.x.lb$year, tmp.y.lb$year); years.lb <- unique(years.lb)
    years.ub <- c(tmp.x.ub$year, tmp.y.ub$year); years.ub <- unique(years.ub)

    ## filter ranks
    ranks.lb <- ranks[ranks[,"year"] %in% years.lb, ]
    ranks.ub <- ranks[ranks[,"year"] %in% years.ub, ]
    
    ## calculated partial correlations in both tails
    r.lb <- sum((ranks.lb[,1]-mean(ranks[,1]))*(ranks.lb[,2]-mean(ranks[,2])))/((n-1)*sqrt(var(ranks[,1])*var(ranks[,2])))
    r.ub <- sum((ranks.ub[,1]-mean(ranks[,1]))*(ranks.ub[,2]-mean(ranks[,2])))/((n-1)*sqrt(var(ranks[,1])*var(ranks[,2])))
    # r<-cor(x_b, y_b, method = "spearman")
    
    ## return partial correlation depending on selected tail
    if(tail == "lower") {
      return(r.lb)
    }
    else {
      return(r.ub)
    }
  }
}

## function used to calculate distance-decay of within-tail synchrony -----
## this function is copy-pasted here for convenience but has not been developed by the authors of the present study
## the original source of the code: Walter, J.A., Castorani, M.C., Bell, T.W., Sheppard, L., Cavanaugh, K.C. & Reuman, D.C. (2022). Tail-dependent spatial synchrony arises from nonlinear driver–response relationships. Ecology Letters, 25, 1189–1201
splineFit<-function(distmat,zmat,nresamp=1000,quantiles=c(0,0.01,0.025,0.05,0.1,0.5,0.9,0.95,0.975,0.99,1)){
  triang<-lower.tri(distmat)
  distmat<-distmat
  xemp<-distmat[triang]
  yemp<-zmat[triang]
  drop.NaNs<-!is.na(yemp)
  dfs=sqrt(nrow(distmat))
  out<-list()
  
  emp.spline<-smooth.spline(xemp[drop.NaNs],yemp[drop.NaNs],df=dfs)
  out$emp.spline<-emp.spline
  
  resamp.splines<-matrix(NA, nrow=nresamp, ncol=length(emp.spline$y))
  for(ii in 1:nresamp){
    shuffle<-sample(1:nrow(distmat), size=nrow(distmat), replace=TRUE)
    xres<-distmat[shuffle,shuffle][triang]
    yres<-zmat[shuffle,shuffle][triang]
    drop.NaNs<-!is.na(yres)
    xres<-xres[drop.NaNs]
    yres<-yres[drop.NaNs]
    yres<-yres[!(xres==0)]
    xres<-xres[!(xres==0)]
    res.spline<-smooth.spline(xres,yres,df=dfs)
    resamp.splines[ii,]<-predict(res.spline, x=emp.spline$x)$y
  }
  out$resamp.splines<-resamp.splines
  out$spline.quantiles<-apply(resamp.splines,2,quantile,probs=quantiles)
  return(out)
}


## extract.edges function
extract.edges <- function(x) {
  ## extract weights and edges
  edges <- x %>% as_edgelist()
  edges <- as.data.frame(edges)
  edges$weight <- as.numeric(E(x)$weight)
  colnames(edges) <- c("site1", "site2", "synchrony")
  
  ## extract edge coords and merge into a single data frame
  from <- xy.tmp
  colnames(from) <- c("site1", "long.from", "lat.from")
  to <- xy.tmp
  colnames(to) <- c("site2", "long.to", "lat.to")
  
  dd <- merge(edges, from, by = "site1")
  dd <- merge(dd, to, by = "site2")
  edges <- na.omit(dd)
  return(edges)
}


##
## matvec -----
## Function transforming similarity matrix into vector.
##

## matvec function
matvec <- function(x, var.name, id.vars) {
  data.sim <- t(combn(colnames(x), 2))
  data.sim <- data.frame(data.sim, sim = x[data.sim])
  colnames(data.sim)[1:2] <- id.vars
  colnames(data.sim)[3] <- var.name
  return(data.sim)
}
