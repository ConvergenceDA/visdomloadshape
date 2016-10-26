
#' @title calculate shannon entropy
#'
#' @description a sequence of load shape cluster assignments is a perfect candidate for calculating the Shannon Entropy of a sequence of symbols.
#'
#' @param p is frequency table or probability distribution vector
#'
#' @export
shannon.entropy=function(p){
  if (min(p) < 0 | sum(p) <= 0)  return(NA)
  p.norm = p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

#' @title calculate shannon entropy
#'
#' @description a sequence of load shape cluster assignments is a perfect candidate for calculating the Shannon Entropy of a sequence of symbols.
#'
#' @param p a sequence of dictionary encoded load days
#'
#' @export
# TODO: remove this function and replace it with an auto-build of a table if necessary in the main function
shannon.entropy2=function(p){
  ##
  return( shannon.entropy( as.numeric(table(p)) ) )
}

#' @title normalize raw smart meter data
#'
#' @description divide each day of consumption, to get rows (days) of normalized consumption
#'
#' @param A consumption data matrix, each row corresponds to a daily profile
#' @param mod mode of operation with 1: Euclidean (divide by the sqrt of the sum of squares) 2: L1 norm (divide by the sum)
#'
#' @export
quick.norm = function(A,mod=2){
  if (mod==1){
    t(apply(A,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sqrt(sum(i^2))}}))
  } else if (mod==2){
    t(apply(A,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sum(i)}}))
  }
}

## modified function which fixed a bug in kNNImpute function in 'imputation' package
## it will be an internal function called from 'impute' function
m.kNNImpute = function (x, k, verbose = T)
{
    if (k >= nrow(x))
        stop("k must be less than the number of rows in x")
    missing.matrix = is.na(x)
    numMissing = sum(missing.matrix)
    if (verbose) {
        print(paste("imputing on", numMissing, "missing values with matrix size",
            nrow(x) * ncol(x), sep = " "))
    }
    if (numMissing == 0) {
        return(x)
    }
    if (verbose)
        print("Computing distance matrix...")
    x.dist = as.matrix(dist(x, upper = T))
    if (verbose)
        print("Distance matrix complete")
    missing.rows.indices = which(apply(missing.matrix, 1, function(i) {
        any(i)
    }))
    if (length(missing.rows.indices)==1) x.missing = matrix((cbind(1:nrow(x), x))[missing.rows.indices, ],nrow=1)
    else x.missing = (cbind(1:nrow(x), x))[missing.rows.indices, ]

    x.missing.imputed = t(apply(x.missing, 1, function(i) {
        rowIndex = i[1]
        i.original = i[-1]
        if (verbose)
            print(paste("Imputing row", rowIndex, sep = " "))
        missing.cols = which(missing.matrix[rowIndex, ])
        if (length(missing.cols) == ncol(x))
            warning(paste("Row", rowIndex, "is completely missing",
                sep = " "))
        imputed.values = sapply(missing.cols, function(j) {
            neighbor.indices = which(!missing.matrix[, j])
            knn.ranks = order(x.dist[rowIndex, neighbor.indices])
            knn = neighbor.indices[(knn.ranks[1:k])]
            mean(x[knn, j])
        })
        i.original[missing.cols] = imputed.values
        i.original
    }))
    x[missing.rows.indices, ] = x.missing.imputed
    missing.matrix2 = is.na(x)
    x[missing.matrix2] = 0
    return(list(x = x, missing.matrix = missing.matrix))
}

#' @title Meter data imputation function
#'
#' @description Imputation function to fill in the blanks of an array of meter data. It assumes the input rows come from the same smart meter and are ordered by date.
#'
#' @param A input data matrix, one row per day of meter data, at least two rows should be valid
#' @param uidx usage column idx, for example, if A is just hourly consumption data (n by 24 matrix), uidx should be 1:24
#'
#' @export
impute=function(A,uidx=4:99){

  n = nrow(A);need = 0 ## need to run knn impute
  n.na = apply(A[,uidx],1,function(i){sum(is.na(i))})
  im.idx = which(n.na>0)

  for (i in im.idx){
    if (n.na[i]==length(uidx)){ ## if it's all NA
      ## put the mean of 2 near days.
      near = order(abs(1:n-i))
      A[i,uidx] = apply(A[near[!(near%in%im.idx)][1:2],uidx],2,mean)
    } else if (n.na[i]==1){ ## if only one point is missing, use linear interpolation as it is faster than doing knn-impute
      na.idx = which(is.na(A[i,uidx]))
      if (na.idx>1 && na.idx<length(uidx)) A[i,uidx[na.idx]] = (A[i,uidx[na.idx-1]]+A[i,uidx[na.idx+1]])/2
      else if (na.idx==1 && i>1 && !is.na(A[i-1,uidx[length(uidx)]])) A[i,uidx[na.idx]] = (A[i-1,uidx[length(uidx)]]+A[i,uidx[na.idx+1]])/2
      else if (na.idx==length(uidx) && i<n && !is.na(A[i+1,uidx[1]])) A[i,uidx[na.idx]] = (A[i,uidx[na.idx-1]]+A[i+1,uidx[1]])/2
      else need = 1
    } else {
      ## if NA value length is short, leave them and use imputation package.
      need = 1
    }
  }
  if (need) A[,uidx] = m.kNNImpute(A[,uidx], k=2, verbose=F)$x ## use 2-NN
  A
}

#' @title encode meter data load shapes according to an existing dictionary
#'
#' @description given a matrix of raw meter data (one day per row), encodes each day into the closest fit culster
#'
#' @param rdata raw data as a matrix with each row as a day of load data
#' @param dic dictionary of cluster centers to use for encodings
#' @param relerror whether to include relative error information in result

#' @return encoding: closest shape code, daily sum, L2 err on normalized data, estimated threshold (relative L2 error)
#'
#' @export
encode = function(rdata, dic, relerror=0){
  if ( !all(complete.cases(rdata) ) ) {
    print(paste('[encode] Warning: missing observations in rdata (NA or NULL) will likely lead to missing value errors.',
                'Please pass in only complete.cases() or fully interpolated data'))
  }

  n = nrow(rdata); p=ncol(rdata)

  ndata = t(apply(rdata,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sum(i)}}))
  encoded = matrix(0,n,4)
  knnres = class::knn(dic,ndata,1:nrow(dic))
  encoded[,1] = knnres
  encoded[,2] = apply(rdata,1,sum)
  encoded[,3] = apply((ndata - dic[knnres,])^2,1,sum)
  encoded[,4] = encoded[,3]/apply(dic[knnres,]^2,1,sum)
  if (relerror) {
    tmp = abs(ndata - dic[knnres,])/ndata
    tmp2 = apply(tmp,1,function(i){
      mean(i[is.finite(i)]) ## use is.finite for the cases that the denominator is zero
    })
    list(encoded=encoded,
         re_perday=apply(abs(ndata - dic[knnres,]),1,sum),
         avg_re_perhour=tmp2)
  } else return(encoded)
}

#' @title reduce the dictionary size hierarchically
#'
#' @description uses k nearest neighbors (\code{class::knn}) to reduce the size of the passed dictionary to a target number of entries by merging the most similar clusters.
#'
#' @param dic dictionary(=cluster centers) /
#' @param cl.size size of each cluster /
#' @param t.num target dictionary size
#' @param d.metric distance metric to use for similarity comparison == 1: Euclidean, ==2: Cosine
#'
#' @export
reduce.dictionary=function(dic,cl.size,t.num=1000,d.metric=1){

  n = nrow(dic)
  dist.mat = matrix(10000,n,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      dist.mat[i,j] = ifelse(d.metric==1,sum((dic[i,]-dic[j,])^2),1-sum(dic[i,]*dic[j,]))
      dist.mat[j,i] = dist.mat[i,j]
    }
  }

  new.cl = c()
  cnt = n
  while(1){
    if (cnt<=t.num) break
    else {
      tmp = which.min(dist.mat)
      xidx = ceiling(tmp/cnt)
      yidx = tmp%%cnt
      yidx = ifelse(yidx==0,cnt,yidx)
      new.size = cl.size[xidx]+cl.size[yidx]
      w.cen = (cl.size[xidx]*dic[xidx,]+cl.size[yidx]*dic[yidx,])/new.size
      if (d.metric==1) {new.cl = w.cen}
      else new.cl = w.cen/sqrt(sum(w.cen^2))

      dic[xidx,] = new.cl
      cl.size[xidx] = new.size
      for (i in (1:nrow(dic))[-xidx]){
        dist.mat[xidx,i] = ifelse(d.metric==1,sum((dic[xidx,]-dic[i,])^2),1-sum(dic[xidx,]*dic[i,]))
        dist.mat[i,xidx] = dist.mat[xidx,i]
      }
      dist.mat = dist.mat[-yidx,-yidx]
      dic = dic[-yidx,]
      cl.size = cl.size[-yidx]
      cnt = cnt-1
    }
  }
  list(n.center = dic, n.cl.size = cl.size)
}


## approximate algorithm to calculate the distance between two lifestyle groups: each group is represented by a lifestyle vector
## lifestyle vector: probability distribution vector of a certain feature
## the distance is calculated by a heuristic to estimate EMD. To find the real EMD, it should solve LP
lifestyle.distance=function(dist.mat,lsv1,lsv2){
  ## dist.mat: distance matrix btw lifestyle components, let's say it's p by p matrix
  ## lsv1, lsv2: lifestyle vector 1,2: the length should be the same with p
  base = pmin(lsv1,lsv2)
  rlsv1 = lsv1 - base
  rlsv2 = lsv2 - base
  nr = sum(rlsv1>0)
  nc = sum(rlsv2>0)
  rdist.mat = matrix(dist.mat[rlsv1>0,rlsv2>0],nr,nc)
  rlsv1 = rlsv1[rlsv1>0]
  rlsv2 = rlsv2[rlsv2>0]
  d = 0

  for (i in order(rdist.mat)){
    ridx = i%%nrow(rdist.mat)
    #ridx = ifelse(ridx==0,nrow(rdist.mat),ridx)
    ridx = ifelse(ridx==0,nr,ridx)
    #cidx = ceiling(i/nrow(rdist.mat))
    cidx = ceiling(i/nr)
    tmp = min(rlsv1[ridx],rlsv2[cidx])
    if (tmp>0){
      rlsv1[ridx] = rlsv1[ridx]-tmp
      rlsv2[cidx] = rlsv2[cidx]-tmp
      d = d + rdist.mat[i]*tmp
    }
  }
  d
}

#' @title plot top n shapes from encoded result
#'
#' @description plots the top n shapes from an encoded result
#' @param en encoded data / en.s: encoded daily sum / dic: dictionary
#' @param n number of top load shapes
#' @param show.avg whether to show avg consumption per load shape. if >0, en.s should have the same dimension with en
#'
#' @export
draw.top.n.shapes = function(en,en.s=0,dic, n=4, show.avg=0){
  en = as.vector(en)
  en.s = as.vector(en.s)
  ten = table(en)
  top.n = sort(ten,decreasing=TRUE)[1:min(n,length(ten))]
  top.u = as.numeric(names(top.n))

  if (show.avg){
    avg.top.n = c()
    ## avg consumption calculation for top10 shape
    for (i in 1:length(top.n)){
      avg.top.n[i] = mean(en.s[en==top.u[i]],na.rm=T)
    }
  }

  par(mfrow=c(ceiling(sqrt(length(top.n))),ceiling(sqrt(length(top.n)))))
  for (i in 1:length(top.n)){
    if (show.avg){
      plot(dic[top.u[i],],xlab='Hour',ylab='Norm.Usage',
           type='o',lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.5,
           main = paste('#',i,': ',round(top.n[i]/length(en)*100,2),'%\n Avg:',round(avg.top.n[i],2),'kWh'))
      grid(lwd=2)
    } else {
      plot(dic[top.u[i],],xlab='Hour',ylab='Norm.Usage',
           type='o',lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.5,
           main = paste('#',i,': ',round(top.n[i]/length(en)*100,2),'%'))
      grid(lwd=2)
    }
  }
  par(mfrow=c(1,1))
}

## calculate EMD btw two different load shapes
calculate.emd = function(a,b){
  ## a,b: two load shape vectors
  l = length(a)
  emd = matrix(0,l+1,1)
  for (i in 1:l){
    emd[i+1] = a[i]+emd[i]-b[i]
  }
  sum(abs(emd))
}

#' @title calculate the distance btw two dictionaries
#'
#' @description calculate the distance btw two dictionaries
#' @param A first dictionary: n1 by p matrix
#' @param B second dictionary: n2 by p matrix
#' @param pa probability distribution vector of dictionary A
#' @param pb probability distribution vector of dictionary B
#' @param emd whether to estimate the distance btw two load shapes as EMD (earth mover distance) or L1 distance/2
#' @param same whether A and B are the same dictionaries or not
#' @param tmpdist user can provide the distance matrix (n1 by n2) btw the dictionaries A and B if already calculated
#'
#' @export
dictionary.distance = function(A,B,pa = 1,pb = 1, emd=1, same=0,tmpdist=NULL){


  n1=nrow(A);p=ncol(A)
  n2=nrow(B)
  if (pa==1) pa = rep(1,n1)/n1
  if (pb==1) pb = rep(1,n2)/n2

  ## calculate the distance btw every pair of components in two dictionaries
  if (is.null(tmpdist)) {
  tmpdist = matrix(0,n1,n2)
    if (same){ ## A == B
      if (emd) { ## EMD
        for (i in 1:(n1-1)){
          a = A[i,]
          for (k in (i+1):n1){
            b = A[k,]
            tmp = rep(0,p+1)
            for (j in 1:p){
              tmp[j+1] = a[j]+tmp[j]-b[j]
            }
            tmpdist[i,k] = sum(abs(tmp))
            tmpdist[k,i] = tmpdist[i,k]
          }
        }
      } else { ## L1 distance/2
        for (i in 1:(n1-1)){
          a = A[i,]
          for (k in (i+1):n1){
            b = A[k,]
            tmpdist[i,k] = sum(abs(a-b))/2
            tmpdist[k,i] = tmpdist[i,k]
          }
        }
      }
    } else { ## A!= B
      if (emd) { ## EMD
        for (i in 1:n1){
          a = A[i,]
          for (k in 1:n2){
            b = B[k,]
            tmp = rep(0,p+1)
            for (j in 1:p){
              tmp[j+1] = a[j]+tmp[j]-b[j]
            }
            tmpdist[i,k] = sum(abs(tmp))
          }
        }
      } else { ## L1 distance/2
        for (i in 1:n1){
          a = A[i,]
          for (k in 1:n2){
            b = B[k,]
            tmpdist[i,k] = sum(abs(a-b))/2
          }
        }
      }
    }
  }

  ## calculate EMD by heuristic
  rdist.mat = tmpdist[pa>0,pb>0]
  nr = nrow(rdist.mat)
  pa = pa[pa>0];pb = pb[pb>0]
  d = 0

  for (i in order(rdist.mat)){
    ridx = i%%nrow(rdist.mat)
    ridx = ifelse(ridx==0,nr,ridx)
    cidx = ceiling(i/nr)
    tmp = min(pa[ridx],pb[cidx])
    if (tmp>0){
      pa[ridx] = pa[ridx]-tmp
      pb[cidx] = pb[cidx]-tmp
      d = d + rdist.mat[i]*tmp
    }
  }
  return(d)
}


## calculate the distances among load shapes within a dictionary with shift flavor
dictionary.distance2 = function(dic,shift.allow=0,d.metric=1){
  ## d.metric = 1: L1 distance
  ## d.metric = 2: L2 distance
  ## d.metric = 3: EMD distance

  ## when shift.allow!=0, it means 1hour shift allowed. Then, compare the avg error, not just sum of errors as it's not fair
  n=nrow(dic);p=ncol(dic)
  dist.mat = matrix(0,n,n)
  if (shift.allow) {
    dist.mat2 = matrix(0,n,n)
    dist.mat3 = matrix(0,n,n)
    dist.mat4 = matrix(0,n,n)
  }

  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (d.metric==3) dist.mat[i,j] = calculate.emd(dic[i,],dic[j,])
      else if (d.metric==2) dist.mat[i,j] = sum((dic[i,]-dic[j,])^2)
      else dist.mat[i,j] = sum(abs(dic[i,]-dic[j,]))
      dist.mat[j,i] = dist.mat[i,j]
      if (shift.allow) {
        if (d.metric==3) {
          dist.mat2[i,j] = calculate.emd(dic[i,],dic[j,c(2:24,1)])
          dist.mat3[i,j] = calculate.emd(dic[i,c(2:24,1)],dic[j,])
        }
        else if (d.metric==2) {
          dist.mat2[i,j] = sum((dic[i,]-dic[j,c(2:24,1)])^2)
          dist.mat3[i,j] = sum((dic[i,c(2:24,1)]-dic[j,])^2)
        }
        else {
          dist.mat2[i,j] = sum(abs(dic[i,]-dic[j,c(2:24,1)]))
          dist.mat3[i,j] = sum(abs(dic[i,c(2:24,1)]-dic[j,]))
        }
        dist.mat2[j,i] = dist.mat2[i,j]
        dist.mat3[j,i] = dist.mat3[i,j]
      }
    }
  }

  if (shift.allow) {
    dist.mat4 = dist.mat ## min distance
    #tmp.idx = dist.mat/24>dist.mat2/23 & dist.mat3>dist.mat2
    tmp.idx = dist.mat>dist.mat2 & dist.mat3>dist.mat2
    dist.mat4[tmp.idx] = dist.mat2[tmp.idx]
    #tmp.idx = dist.mat/24>dist.mat3/23 & dist.mat2>dist.mat3
    tmp.idx = dist.mat>dist.mat3 & dist.mat2>dist.mat3
    dist.mat4[tmp.idx] = dist.mat3[tmp.idx]
    return(dist.mat4)
  } else return(dist.mat)
}

#' @title Perform load shape clusing and return cluster centers
#'
#' @description This function takes load shape data in a standardized matrix format and creates a dictionary of resulting cluster centers
#'
#' @param sdata source data, assume it's already standardized (cleansed and n by p matrix format)
#' @param target.size target size of the dictionary (i.e. nubmer of clusters)
#' @param mode 1: use ths1, 2: use ths2, 3: use ths3, 4: use ths4
#' @param d.metric 1: use euclidean distance metric, otherwise use cosine distance metric
#' @param ths will be transferred to akmeans parameter according to mode setting
#' \code{ths1}: threshold to decide whether to increase k or not: check sum((sample-assigned center)^2) < ths1*sum(assigned center^2)
#' \code{ths2}: threshold to decide whether to increase k or not: check all components of |sample-assigned center| < ths2
#' \code{ths3}: threshold to decide whether to increase k or not: check inner product of (sample,assigned center) > ths3 , this is only for cosine distance metric
#' \code{ths4}: threshold to decide whether to increase k or not: check sum(abs(sample-assigned center)) < ths4
#' @param iter.max maximum iteration setting to be used in kmeans
#' @param n.start parameter to be transferred to kmeans
#' @param two.step.compress whether to reduce the dictionary only by hierarchical clustering or hier+use top N shapes.
#' this option gets activated only when the ratio (original dictionary size before compression/target.size) is larger than 10
#' @param verbose whether to show log or not
#'
#' @export
create_dictionary = function(sdata,target.size=1000, mode=1, d.metric=1, ths=0.2, iter.max=100,
							 nstart=1, two.step.compress=F,verbose=F){


  sdata = sdata[apply(sdata,1,sum)>0,]
  if (d.metric==1) sdata = quick.norm(sdata) ## if euclidean distance, normalize here. for cosine, akmeans will handle

  gc()
  akmres = akmeans(x = sdata, min.k = round(target.size/4), max.k = 99999,
                   mode=mode, d.metric=d.metric, ths1=ths,ths2=ths,ths3=ths,ths4=ths, iter.max=iter.max, nstart=nstart,verbose=verbose)
  if (two.step.compress & nrow(akmres$centers)/target.size>10) {
    ratio = nrow(akmres$centers)/target.size
    #print(akmres$size)
    #print(dim(akmres$centers))
    rdic1 = reduce.dictionary(akmres$centers,akmres$size,t.num=round(target.size*sqrt(ratio)),d.metric=d.metric) ## from empirical tests, root is a sweet spot
    rdic = rdic1$n.center[order(rdic1$n.cl.size,decreasing=T)[1:target.size],]
  } else {
    rdic = reduce.dictionary(akmres$centers,akmres$size,t.num=target.size,d.metric=d.metric)$n.center
  }
  return(rdic)
}

#' @title Cluster meter data into the best fit shape clusters and return clsuter assignments
#'
#' @description Useful for generation of a load shape cluster center dictionary and encoding in one shot
#'
#' @param rdata rawdata of format n by p matrix
#' @param is.clean whether to do data interpolation or not.
#' Note that the default interpolation method is very memory and CPU intensive.
#' You may be better off doing your own interpolation or passing only complete.cases().
#' @param use.all whether to use all data to generate a dictionary or not
#' @param s.size sample size to use to generate a dictionary
#' @param target.size target size of the dictionary (i.e. nubmer of clusters)
#' @param mode 1: use ths1, 2: use ths2, 3: use ths3, 4: use ths4
#' @param d.metric 1: use euclidean distance metric, otherwise use cosine distance metric
#' @param ths will be transferred to akmeans parameter according to mode setting
#' \code{ths1}: threshold to decide whether to increase k or not: check sum((sample-assigned center)^2) < ths1*sum(assigned center^2)
#' \code{ths2}: threshold to decide whether to increase k or not: check all components of |sample-assigned center| < ths2
#' \code{ths3}: threshold to decide whether to increase k or not: check inner product of (sample,assigned center) > ths3 , this is only for cosine distance metric
#' \code{ths4}: threshold to decide whether to increase k or not: check sum(abs(sample-assigned center)) < ths4
#' @param iter.max maximum iteration setting to be used in kmeans
#' @param n.start parameter to be transferred to kmeans
#' @param two.step.compress whether to reduce the dictionary only by hierarchical clustering or hier+use top N shapes.
#' this option gets activated only when the ratio (original dictionary size before compression/target.size) is larger than 10
#' @param verbose whether to show log or not
#'
#' @export
raw2encoded = function(rawdata, is.clean = T, use.all = T, s.size = 100000, target.size=1000,
                       mode=1, d.metric=1, ths=0.2, iter.max=100, nstart=1, two.step.compress=F,verbose=F){


  ## 1. cleanse first (impute all zero rows and any null values) but it doesn't detect outliers
  if (!is.clean) {
    print('[raw2encoded] Warning running data imputation because is.clean=F. This can be very memory and CPU intensive.')
    cdata = impute(rawdata,uidx=1:ncol(rawdata))
  }
  else cdata = rawdata
  if ( !all(complete.cases(cdata) ) ) {
    print(paste('[raw2encoded] Warning: missing observations in rawdata (NA or NULL) will likely lead to missing value errors.',
          'Please re-run with is.clean=F, interpolate missing values prior to passing data in, or',
          'pass in only complete.cases()'))
  }

  ## 2. create a dictionary
  if (use.all) {
    dic = create_dictionary(cdata, target.size=target.size, mode=mode, d.metric=d.metric, ths=ths, iter.max=iter.max,
                                       nstart=nstart, two.step.compress=two.step.compress, verbose=verbose)
  } else {
    sidx = sample(nrow(cdata),s.size)
    dic = create_dictionary(cdata[sidx,], target.size=target.size, mode=mode, d.metric=d.metric, ths=ths, iter.max=iter.max,
                            nstart=nstart, two.step.compress=two.step.compress, verbose=verbose)
  }

  ## 3. encode the data: 4 columns: load shape code, daily sum, square err, relative square err
  encoded = encode(cdata,dic)

  list(clean.data = cdata, dictionary = dic, encoded.data = encoded)
}


