# Copyright 2016 The Board of Trustees of the Leland Stanford Junior University.
# Direct inquiries to Sam Borgeson (sborgeson@stanford.edu)
# or professor Ram Rajagopal (ramr@stanford.edu)

getClosest=function(X,y,k,mode=1){
## find the closest row from y among X rows
## X: n by p matrix: pool
## y: input vector
## k: find till kth closest row
## mode: 1: L1 distance 2: L2 euclidean distance 3: cosine similarity
n = dim(X)[1]; p=dim(X)[2]
if (mode==1) {
  diff = apply(abs(X-matrix(rep(y,each=n),nrow=n)),1,sum)
} else if(mode==2) {
  diff = apply((X-matrix(rep(y,each=n),nrow=n))^2,1,sum)
} else if(mode==3) {
  diff = apply(X,1,function(i){1 - sum(i*y)/sqrt(sum(i^2)*sum(y^2))})
}
list(idx=which(rank(diff)%in%1:k),dist=diff[which(rank(diff)%in%1:k)])
## return the distance and the base idx
}

## will give you qkw1,qkw2,...,qkw96
q.cols = paste(paste('qkw',1:96,sep=''),collapse=',')

## will give you qkw1+qkw2+...+qkw96
q.cols.sum = paste(paste('qkw',1:96,sep=''),collapse='+')

## will give you hkw1,hkw2,...,hkw24
h.cols = paste(paste('hkw',1:24,sep=''),collapse=',')

## will give you hkw1+hkw2+...+hkw24
h.cols.sum = paste(paste('hkw',1:24,sep=''),collapse='+')


## modified function which fixed a bug in kNNImpute function in 'imputation' package
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

## imputation function
kjs.impute=function(A,uidx=4:99){
## assume the input comes for the same spid and ordered by date
## at least two rows should be valid
## by default uidx = 4:99 : usage column idx
  n = nrow(A);need = 0 ## need to run knn impute
  im.idx = which(apply(A[,uidx],1,function(i){
    if (sum(is.na(i))>0) return(1)
    else if (sum(i)==0) return(1)
    else return(0)
  })==1)

  ## return fail if all values are zero
  if (sum(A[,uidx],na.rm=T)==0) {print("all values are zero=>can't impute");return(A)}

  for (i in im.idx){
    if (any(is.na(A[i,uidx]))){ ## rows containing NA value

      if (sum(is.na(A[i,uidx]))==length(uidx)){ ## if it's all NA
        ## put the mean or mean of 2 near days
        ## A[i,uidx] = apply(A[-im.idx,uidx],2,mean) ## mean of all data, not good
        near = order(abs(1:n-i))
        A[i,uidx] = apply(A[near[!(near%in%im.idx)][1:2],uidx],2,mean)
      } else {
        ## if NA value length is short, leave them and use imputation package.
	need = 1
      }
    } else if (sum(A[i,uidx])==0){ ## all zero rows: treat as same with all NA
      near = order(abs(1:n-i))
      A[i,uidx] = apply(A[near[!(near%in%im.idx)][1:2],uidx],2,mean)
    }
  }
  if (need) A[,uidx] = m.kNNImpute(A[,uidx], k=2, verbose=F)$x ## use 2-NN
  A
}

## get 4 binned usage rank for 96points
get.rank=function(y,k=4,div=c(16,40,64,88),mode=1){
  ## mode 1 returns rank info of length 4
  ## mode 2 returns rank info of length 2
  m1label = c(1234, 1243, 1324, 1342, 1423, 1432, 2134, 2143, 2314, 2341, 2413, 2431,
              3124, 3142, 3214, 3241, 3412, 3421, 4123, 4132, 4213, 4231, 4312, 4321)
  m2label = c(12, 13, 14, 21, 23, 24, 31, 32, 34, 41, 42, 43)
  ## assume division is less than 10 slots
  s = matrix(0,1,k)
  for (i in 1:(k-1)){
    s[i] = sum(y[(div[i]+1):div[i+1]])
  }
  s[k] = sum(y[1:div[1]])
  if (div[k]<96) s[k] = s[k] + sum(y[(div[k]+1):96])
  if (mode ==1) {
    rank = as.numeric(paste(order(s,decreasing=TRUE),collapse=''))
    list(rank=rank,label=which(m1label==rank),ps=s)
  } else if (mode==2){
    rank = substr(paste(order(s,decreasing=TRUE),collapse=''),1,2)
    list(rank=rank,label=which(m2label==rank),ps=s)
  }
}

get.rank24=function(y,k=4,div=c(4,10,16,22),mode=1){
  ## mode 1 returns rank info of length 4
  ## mode 2 returns rank info of length 2
  m1label = c(1234, 1243, 1324, 1342, 1423, 1432, 2134, 2143, 2314, 2341, 2413, 2431,
              3124, 3142, 3214, 3241, 3412, 3421, 4123, 4132, 4213, 4231, 4312, 4321)
  m2label = c(12, 13, 14, 21, 23, 24, 31, 32, 34, 41, 42, 43)
  ## assume division is less than 10 slots
  s = matrix(0,1,k)
  for (i in 1:(k-1)){
    s[i] = sum(y[(div[i]+1):div[i+1]])
  }
  s[k] = sum(y[1:div[1]])
  if (div[k]<24) s[k] = s[k] + sum(y[(div[k]+1):24])
  if (mode ==1) {
    rank = as.numeric(paste(order(s,decreasing=TRUE),collapse=''))
    list(rank=rank,label=which(m1label==rank),ps=s)
  } else if (mode==2){
    rank = substr(paste(order(s,decreasing=TRUE),collapse=''),1,2)
    list(rank=rank,label=which(m2label==rank),ps=s)
  }
}

## get date in the format wanted, and return the day of week info: Monday is 1, and....
get.day = function(date, format = '%Y-%m-%d'){
  days = c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')
  day.ind = c()
  for (i in 1:length(date)){
    day.ind[i] = which(days==weekdays(as.Date(as.character(date[i]),format=format)))
  }
  day.ind
}

## normalization function
quick.norm = function(A,mod=2){
  if (mod==1){ ## L2 normalization
    t(apply(A,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sqrt(sum(i^2))}}))
  } else if (mod==2){ ## L1 normalization
    t(apply(A,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sum(i)}}))
  }
}

ch.24 = function(A){
  if (dim(A)[2]==24) return(A)
  else {
    t(apply(A,1,function(i){
      apply(matrix(i,nrow=4),2,sum)
    }))
  }
}

## encoding function (load shape code, daily consumption, L2 error, relative error
kjs.encode = function(rdata, dic, relerror=0){
  ## rdata: raw data
  ## dic: dictionary
  ## encoding: closest shape code, daily sum, err on normalized data, est ths
  n = nrow(rdata); p=ncol(rdata)
  ndata = t(apply(rdata,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sum(i)}}))
  encoded = matrix(0,n,4)
  require(class)
  knnres = knn(dic,ndata,1:nrow(dic))
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

## reduce the dictionary size
reduce.dic=function(dic,cl.size,t.num=1000,d.metric=1){
  ## d.metric == 1: Euclidean, ==2: Cosine
  ## dic: dictionary(=cluster centers) / cl.size: size of each cluster / t.num: target dictionary size
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

## load profile emulator version 1
## get the data of individual and emulate its load profiles
## independent from daily sum emulator
## this version needs 1 year data to make the seasonal effect
########################################################
## model description ##
## Y_i: 24 by 1 matrix: normalized load shape for day i
## Y_hat_i: estimated Y on day i
## minimize sum of |Y_i-Y_hat_i|^2 for i=2,..,365
## Y_hat_i = w1 * sum of (P(Cj|yesterday shape)*Cj)
##         + w2 * sum of (P(Cj|weekday)*Cj)
##         + w3 * sum of (P(Cj|season)*Cj)
##         + w4 * sum of (P(Cj|special day)*Cj)
## subject to: sum of wi = 1, wi>=0
## then based on w1,w2,w3,w4, we know P(Cj) on day i=> randomly populate the load shape
kjs.lp.emul=function(dic,encoded,s.date='2010-08-01',e.date='2011-07-31',emul.sdate='2011-08-01',emul.edate='2012-07-31'){
  ## dic: dictionary
  ## encoded: encoded shape code
  ## s.date: starting date
  ## e.date: ending date

  ## training process
  wd = c("Monday","Tuesday" ,  "Wednesday", "Thursday" , "Friday"  ,  "Saturday" , "Sunday")
  len = length(encoded) ## length of train days
  wdidx = (matrix(0:6,1,len) + which(wd==weekdays(as.Date(s.date))))%%7
  wdidx[wdidx==0]=7 ## to make Sunday as 7
  ## P(Cj|weekday) calculation
  p.wd = matrix(0,7,nrow(dic))
  for (i in 1:7){
    p.wd[i,sort(unique(encoded[wdidx==i]))] = table(encoded[wdidx==i])/sum(wdidx==i)
  }

  ## P(Cj|yesterday) calculation
  p.ye = matrix(0,nrow(dic),nrow(dic))
  for (i in 1:(len-1)){
    p.ye[encoded[i],encoded[i+1]] = p.ye[encoded[i],encoded[i+1]] + 1
  }
  p.ye = t(apply(p.ye,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sum(i)}}))

  ## P(Cj|season) calculation
  ## how to set the season
  ## divide into 2 or 4 fixed seasons?
  ## Or, set sliding window and calculate P(Cj|season) for every day?
  ## for now, fix 4 seasons and calculate!!
  ## practically, the date would be 2010-08-01 ~ 2011-07-31
  ## take the california season division from web:
  ## http://answers.yahoo.com/question/index?qid=20110227141200AAYDfcn
  ## Summer: June 22 to September 21
  ## Fall: September 22 to December 21
  ## Winter: December 22 to March 21
  ## Spring: March 22 to June 21
  p.se = matrix(0,4,nrow(dic))
  for (i in 1:365){
    if (i<53 | i>325) p.se[1,encoded[i]] = p.se[1,encoded[i]] + 1/92 ## summer
    else if (i>52 & i<144) p.se[2,encoded[i]] = p.se[2,encoded[i]] + 1/91 ## fall
    else if (i>143 & i<235) p.se[3,encoded[i]] = p.se[3,encoded[i]] + 1/91 ## fall
    else p.se[4,encoded[i]] = p.se[4,encoded[i]] + 1/91 ## fall
  }

  ## P(Cj|specific day) calculation
  ## # of classes: for now 2=> 1.not holidays 2.holidays
  ## holiday except weekends
  require(timeDate)
  hd=c("2010-09-06", "2010-11-25", "2010-12-24" ,"2011-01-17","2011-02-21","2011-04-22" ,"2011-05-30" ,"2011-07-04") ## in training period
  hdidx = as.numeric(as.Date(hd)-as.Date(s.date))+1
  #weekdays(as.Date(hd))
  ## "Monday"   "Thursday" "Friday"   "Monday"   "Monday"   "Friday"   "Monday"  "Monday"
  p.sp = matrix(0,2,nrow(dic))
  p.sp[1,sort(unique(encoded[-hdidx]))] = table(encoded[-hdidx])/(len-length(hd))
  p.sp[2,sort(unique(encoded[hdidx]))] = table(encoded[hdidx])/(length(hd))

  ##########################################################
  ## calculate w1,w2,w3,w4 by fitting minimize sum of |Y_i-Y_hat_i|^2 for i=2,..,365
  ## Y_hat = X*b
  ## make Y
  Y = matrix(t(dic[encoded[2:len],]))
  ## make X
  se.idx = 0
  X = c()
  for (i in 2:365){
    if (i<53 | i>325) se.idx=1 ## summer
    else if (i>52 & i<144) se.idx=2 ## fall
    else if (i>143 & i<235) se.idx=3 ## winter
    else se.idx=4 ## spr

    X = rbind(X,cbind(t(dic)%*%p.ye[encoded[i-1],], ## sum of P(Cj|yesterday)*Cj
          t(dic)%*%p.wd[wdidx[i],], ## sum of (P(Cj|weekday)*Cj)
          t(dic)%*%p.se[se.idx,], ## sum of (P(Cj|season)*Cj)
          t(dic)%*%p.sp[ifelse(i%in%hdidx,2,1),] ## sum of (P(Cj|special day)*Cj)
    ))
  }

  ## solve quadratic programming
  require(quadprog)
  Dmat = t(X)%*%X
  dvec = t(X)%*%Y
  Amat = t(rbind(rep(1,4),diag(4)))
  bvec = c(1,0,0,0,0)
  sol = solve.QP(Dmat, dvec, Amat, bvec, meq=1)

  w = sol$solution ## now got the weight

  ## start emulation by given emul.len
  emul.len = as.numeric(as.Date(emul.edate)-as.Date(emul.sdate))+1
  emul.encode = matrix(0,1,emul.len+1)

  ## starting point
  tmp = as.numeric(as.Date(e.date)-as.Date(emul.sdate))
  emul.encode[1] = ifelse(tmp>-2,encoded[365-tmp-1], ## overlap=> use real starting point
                          sort(unique(encoded))[which.max(table(encoded))]) ## not overlap => choose most frequent shapes

  ## holiday check in emulation period
  ## as.Date(holidayNYSE(2011))
  chd = as.Date(holidayNYSE(as.numeric(substr(as.Date(emul.sdate),1,4)):as.numeric(substr(as.Date(emul.edate),1,4)))) ## candidates
  ehd = chd[chd >= emul.sdate & chd <= emul.edate]
  ehdidx = as.numeric(as.Date(ehd)-as.Date(emul.sdate))+1

  ## weekday idx
  ewdidx = (matrix(0:6,1,emul.len) + which(wd==weekdays(as.Date(emul.sdate))))%%7
  ewdidx[ewdidx==0]=7
  ## Summer: June 22 to September 21
  ## Fall: September 22 to December 21
  ## Winter: December 22 to March 21
  ## Spring: March 22 to June 21
  for (i in 1:emul.len){
    d = substr(as.Date(emul.sdate)+i-1,6,10)
    if (d > '06-21' & d < '09-22') {se.idx=1}
    else if (d > '09-21' & d < '12-22') {se.idx=2}
    else if (d > '12-21' | d < '03-22') {se.idx=3}
    else {se.idx=4} ## spr

    prob = cbind(p.ye[emul.encode[i],],p.wd[ewdidx[i],],p.se[se.idx,],p.sp[ifelse(i%in%ehdidx,2,1),])%*%w
    prob[prob<0]=0
    emul.encode[i+1] = sample(1:nrow(dic),1,prob=prob)
  }
  list(w = w, emul.res = emul.encode[-1] )
}

##test = kjs.lp.emul(dictionary,as.numeric(knn(dictionary,quick.norm(t(matrix(dmatrix_1hour[,2],nrow=24))),1:1000)))

## response model 0: hourly model without any breakpoint
fit.m0 = function(udata,temp.info){
  ## m0 is U(t, d) = A+(t) * T(t, d) + C(t) + e
  ## udata is n by p usage data
  ## temp.info is n by p temperature data

  n = nrow(udata);p=ncol(udata)
  coefs = matrix(0,p,2)
  tmp.rss = matrix(0,1,p)
  rss = matrix(10000,1,p)
  est = matrix(0,n,p)
  std = matrix(0,1,p)
  fit.data = c()
  for (i in 1:p){
    fit = lm(udata[,i]~temp.info[,i])
    est[,i] = fit$fit
    rss[i] = sum(fit$res^2)
    coefs[i,] = fit$coef
    std[i] = summary(fit)[[4]][2,2]
  }
  list(rss=sum(rss),est=est,coefs=coefs,std=std)
}

## response model 0 returning rss as they are
fit.m0.rss = function(udata,temp.info){
  ## m0 is U(t, d) = A+(t) * T(t, d) + C(t) + e
  ## udata is n by p usage data
  ## temp.info is n by p temperature data

  n = nrow(udata);p=ncol(udata)
  coefs = matrix(0,p,2)
  tmp.rss = matrix(0,1,p)
  rss = matrix(10000,1,p)
  est = matrix(0,n,p)
  std = matrix(0,1,p)
  fit.data = c()
  for (i in 1:p){
    fit = lm(udata[,i]~temp.info[,i])
    est[,i] = fit$fit
    rss[i] = sum(fit$res^2)
    coefs[i,] = fit$coef
    std[i] = summary(fit)[[4]][2,2]
  }
  list(rss=rss,est=est,coefs=coefs,std=std)
}

## response model 1: hourly model with one breakpoint which changes by given period (dur), but same breakpoint for every hour
fit.m1 = function(udata,temp.info,cand.tref=55:80,dur=5,tr.ch=2){
  ## cand.tref: candidate integer breakpoints (=reference temperature)
  ## dur: period of updating the reference temperature => if we exclude weekends dur=5 means one week
  ## udata is n by 24 usage data
  ## temp.info is n by 24 temperature data
  ## m1 is U(t, d) = A+(t) * (T(t, d) - Tr(w(d)))+  +  A-(t) * (Tr(w(d)) - T(t,d))+  + C(t) + e
  ## Tr(w(d)) is reference temperature of the week containing date d

  n = nrow(udata)
  ## at first, set Tr(w(d)) among cand.tref
  ## start from U(t, d) = A+(t) * (T(t, d) – Tr)+  +  A-(t) * (Tr - T(t,d))+  + C(t) + e
  tmp.rss = c()
  for (i in cand.tref){
    res = try.tref(udata,temp.info,i)
    tmp.rss = c(tmp.rss,res$rss)
  }
  tref = cand.tref[which.min(tmp.rss)] ## now set coefs starting point
  Trseq = rep(tref,ceiling(n/dur))
  estep = try.tref(udata,temp.info,tref)

  mstep = fit.tref.seq(estep$coefs,udata,temp.info,cand.tref,dur,tr.ch)
  estep = fit.coefs(mstep$Trseq,udata,temp.info,dur)
  prev.mstep = mstep
  prev.estep = estep

  ## start while loop for EM style algorithm
  while(1){
    ## with fixed coefs find Tr(w(d)) sequence
    mstep = fit.tref.seq(estep$coefs,udata,temp.info,cand.tref,dur,tr.ch)
    estep = fit.coefs(mstep$Trseq,udata,temp.info,dur)
    if (prev.estep$rss<=estep$rss) break
    prev.mstep = mstep
    prev.estep = estep
  }
  list(rss=prev.estep$rss,est=prev.estep$est,coefs=prev.estep$coefs, trseq=prev.mstep$Trseq)
}

## for each period, find Tref
fit.tref.seq=function(coefs,udata,temp.info,cand.tref,d,tr.ch){
  cand.ref.t = cand.tref
  coefs[is.na(coefs)]=0
  cnt=0;p=ncol(udata)
  datalen = nrow(udata)
  trseq = matrix(0,1,ceiling(datalen/d))
  rss = matrix(0,1,ceiling(datalen/d))
  for (j in 1:ceiling(datalen/d)){
    ls = min(d,datalen-d*j+d); err = c()
    for (i in cand.ref.t){
      tmp.err = 0
      for (k in 1:ls){
        tmp.err = tmp.err + sum((udata[d*cnt+k,]-coefs[seq(2,3*p,3)]*pmax(temp.info[d*cnt+k,]-i,0)-coefs[seq(3,3*p,3)]*pmax(i-temp.info[d*cnt+k,],0)-coefs[seq(1,3*p,3)])^2,na.rm=T)
      }
      err = c(err,tmp.err)
    }
    cnt = cnt+1
    trseq[cnt] = cand.ref.t[which.min(err)]
    rss[cnt] = min(err)
    cand.ref.t = (trseq[cnt]-tr.ch):(trseq[cnt]+tr.ch) ## plus minus 2 allowed
    cand.ref.t = cand.ref.t[cand.ref.t>=min(cand.tref) & cand.ref.t<=max(cand.tref)]
  }
  list(Trseq = trseq, rss = sum(rss))
}

## given Tref sequence, find coefs
fit.coefs=function(ref.temp,udata,temp.info,d){
  fit.data= c();cnt=0
  n = nrow(udata);p=ncol(udata)
  ref.t = matrix(rep(ref.temp,each=d),n,1)
  coefs=matrix(0,p,3);rss=0;est=matrix(0,n,p)

  for (k in 1:p){
    fit = lm(udata[,k]~pmax(temp.info[,k]-ref.t,0)+pmax(ref.t-temp.info[,k],0))
    est[,k] = fit$fit
    coefs[k,] = fit$coef
    rss = rss + sum(fit$res^2,na.rm=T)
  }
  coefs[is.na(coefs)]=0
  list(coefs=coefs, rss=rss, est=est)
}

## fit U(t, d) = A+(t) * (T(t, d) - Tr)+  +  A-(t) * (Tr - T(t,d))+  + C(t) + e with Tr = tref
try.tref=function(udata,temp.info,tref){
  rss = 0;coefs=c()
  for (k in 1:24){
    fit = lm(udata[,k]~pmax(temp.info[,k]-tref,0)+pmax(tref-temp.info[,k],0))
    coefs=c(coefs,fit$coef)
    rss = rss + sum(fit$res^2,na.rm=T)
  }
  coefs[is.na(coefs)]=0
  list(rss=rss,coefs=coefs)
}

## response model 2: hourly model with one breakpoint, and different breakpoints for every hour
fit.m2 = function(udata,temp.info,cand.tref=55:80){
  ## udata is n by 24 usage data
  ## temp.info is n by 24 temperature data
  ## m2 is U(t, d) = A+(t) * (T(t, d) - Tr(t))+  +  A-(t) * (Tr(t) - T(t,d))+  + C(t) + e
  ## m2 is time wise model. Just fit 24 times
  n = nrow(udata);p=ncol(udata)

  coefs = matrix(0,p,3)
  tmp.rss = matrix(0,1,p)
  rss = matrix(10000,1,p)
  est = matrix(0,n,p)
  tref = matrix(0,1,p)
  std = matrix(0,1,p)
  for (i in cand.tref){
    for (k in 1:p){
      fit = lm(udata[,k]~pmax(temp.info[,k]-i,0)+pmax(i-temp.info[,k],0))
      tmp.rss = sum(fit$res^2,na.rm=T)
      if (tmp.rss < rss[k]) {
        rss[k] = tmp.rss
        coefs[k,] = fit$coef
        est[,k] = fit$fit
        tref[k] = i
        std[k] = summary(fit)[[4]][2,2]
      }
    }
  }
  coefs[is.na(coefs)]=0
  list(rss=sum(rss),est=est,coefs=coefs,trseq=tref,std=std)
}

## response model 2 returning rss as they are
fit.m2.rss = function(udata,temp.info,cand.tref=55:80){
  ## udata is n by 24 usage data
  ## temp.info is n by 24 temperature data
  ## m2 is U(t, d) = A+(t) * (T(t, d) - Tr(t))+  +  A-(t) * (Tr(t) - T(t,d))+  + C(t) + e
  ## m2 is time wise model. Just fit 24 times
  n = nrow(udata);p=ncol(udata)

  coefs = matrix(0,p,3)
  tmp.rss = matrix(0,1,p)
  rss = matrix(10000,1,p)
  est = matrix(0,n,p)
  tref = matrix(0,1,p)
  std = matrix(0,1,p)
  for (i in cand.tref){
    for (k in 1:p){
      fit = lm(udata[,k]~pmax(temp.info[,k]-i,0)+pmax(i-temp.info[,k],0))
      tmp.rss = sum(fit$res^2,na.rm=T)
      if (tmp.rss < rss[k]) {
        rss[k] = tmp.rss
        coefs[k,] = fit$coef
        est[,k] = fit$fit
        tref[k] = i
        std[k] = summary(fit)[[4]][2,2]
      }
    }
  }
  coefs[is.na(coefs)]=0
  list(rss=rss,est=est,coefs=coefs,trseq=tref,std=std)
}

## response model 2 to return heating coefficient standard deviation
fit.m2.heat = function(udata,temp.info,cand.tref=55:80){
  ## udata is n by 24 usage data
  ## temp.info is n by 24 temperature data
  ## m2 is U(t, d) = A+(t) * (T(t, d) - Tr(t))+  +  A-(t) * (Tr(t) - T(t,d))+  + C(t) + e
  ## m2 is time wise model. Just fit 24 times
  n = nrow(udata);p=ncol(udata)

  coefs = matrix(0,p,3) ## A+(t), A-(t), C(t)
  tmp.rss = matrix(0,1,p)
  rss = matrix(10000,1,p)
  est = matrix(0,n,p)
  tref = matrix(0,1,p)
  std = matrix(0,1,p)
  for (i in cand.tref){
    for (k in 1:p){
      fit = lm(udata[,k]~pmax(temp.info[,k]-i,0)+pmax(i-temp.info[,k],0))
      tmp.rss = sum(fit$res^2,na.rm=T)
      if (tmp.rss < rss[k]) {
        rss[k] = tmp.rss
        coefs[k,] = fit$coef
        est[,k] = fit$fit
        tref[k] = i
        if (is.na(fit$coef[3])) {
          std[k] = Inf
        } else std[k] = summary(fit)[[4]][3,2]
      }
    }
  }
  coefs[is.na(coefs)]=0
  list(rss=rss,est=est,coefs=coefs,trseq=tref,std=std)
}

## response model 3: hourly model with one breakpoint which changes by given period (dur), and different breakpoints for every hour
fit.m3 = function(udata,temp.info,cand.tref=55:80,d=5,coefs,fit.t=1:24,tr.ch=2){
  ## udata is n by 24 usage data
  ## temp.info is n by 24 temperature data
  ## coefs is from m2 result. 24 by 3 matrix
  ## m3 is U(t, d) = A+(t) * (T(t, d) – Tr(w(d),t))+  +  A-(t) * (Tr(w(d),t) - T(t,d))+  + C(t) + e
  ## m3 starts from m2. fix coefs first for each t, and find Tr seq for each t
  ## so, run m2, m3 at same time

  n = nrow(udata);p=ncol(udata)
#   fit.t = 1:min(p,24)
  std = matrix(0,1,p)
  prev.rss = 10000000
  trseq = matrix(0,p,ceiling(n/d))
  while(1){
    rss=0
    for (k in fit.t){
      cand.ref.t = cand.tref
      cnt=0

      for (j in 1:ceiling(n/d)){
        ls = min(d,n-d*j+d); err = c()
        for (i in cand.ref.t){
          err = c(err,sum((udata[d*cnt+1:ls,k]-coefs[k,2]*pmax(temp.info[d*cnt+1:ls,k]-i,0)-coefs[k,3]*pmax(i-temp.info[d*cnt+1:ls,k],0)-coefs[k,1])^2,na.rm=T))
        }
        cnt = cnt+1
        trseq[k,j] = cand.ref.t[which.min(err)]
        rss = rss + min(err)
        cand.ref.t = (trseq[k,j]-tr.ch):(trseq[k,j]+tr.ch) ## plus minus 2 allowed
        cand.ref.t = cand.ref.t[cand.ref.t>=min(cand.tref) & cand.ref.t<=max(cand.tref)]
      }
    }
    ## trseq is set now

    ## need to fit coefs with fixed trseq
    coefs=matrix(0,24,3);rss=0;est=c()
    tmpstd=c()
    for (k in fit.t){
      fit = lm(udata[,k]~pmax(temp.info[,k]-matrix(rep(trseq[k,],each=d),n,1),0)+pmax(matrix(rep(trseq[k,],each=d),n,1)-temp.info[,k],0))
      est = rbind(est,fit$fit)
      tmpstd = c(tmpstd,summary(fit)[[4]][2,2]) ## std of A+ coef
      coefs[k,]=fit$coef
      rss = rss + sum(fit$res^2,na.rm=T)
    }
    if (prev.rss <= rss) break
    prev.rss = rss
    prev.est = est
    prev.coefs = coefs
    prev.trseq = trseq
    std[fit.t] = tmpstd
  }
  prev.coefs[is.na(prev.coefs)]=0
  list(rss=prev.rss,est=t(prev.est),coefs=prev.coefs,trseq=prev.trseq,std=std)
}

## solve QP problem => after giving sol, it will select largest N guys
## given lambda, solve arg min -lambda*t(mu)%*%x + t(x)%*%cov%*%x s.t sum of x <= N, x_i = 0 or 1
## => relax the integer constraint to 0~1
## solve.QP: min -db+1/2*bDb s.t. Ab >= b0
solve.subQP=function(lambda,mu,covar,N){
  require(quadprog)
  Dmat = 2*covar
  dvec = lambda*mu
  Amat = t(rbind(rep(-1,length(mu)),diag(length(mu)),-diag(length(mu))))
  #Amat = t(rbind(diag(length(mu)),-diag(length(mu))))
  bvec = c(-N,rep(0,length(mu)),rep(-1,length(mu)))
  #bvec = c(rep(0,length(mu)),rep(-1,length(mu)))
  sol = solve.QP(Dmat, dvec,Amat,bvec)
  sol
}

## solve integer lp problem
## p(sx>T)>eps = (mu+beta*sigma)x > target
## minimize cx, when c=const it's minimizing the number of people'
## it can solve independent case
solve.intLP=function(mu,sigma,target,eps){
  require(lpSolve)
  beta = qnorm(1-eps)
  obj = rep(1,length(mu))
  cons = matrix(mu+beta*sigma,nrow=1)## constraint: one row is one constraint
  cond = '>='
  sol = lp (direction = "min", objective.in=obj, const.mat=cons, const.dir=cond, const.rhs=target,
      all.int=F, all.bin=T)
  sol
}

## solve the customer targeting algorithm
solve.alg=function(mu,cov,Tes,N,mode=1,M=10){
  ## mode=1 :assume the covariance matrix is diagonal
  ## mode=2 :solve subQP if it's the covariance matrix is not diagonal
  ## mu: response mean
  ## cov: covariance matrix of responses
  ## Tes: Target energy saving
  ## N: limit of the number of customers to enroll
  ## M: design parameter to decide the step size of increasing the angle of the slope in the transferred domain

  tmpsol = matrix(0,M+1,length(mu))
  #a = tan((0:M)*pi/4/M+pi/4) ## slope to test
  a = tan((0:M)*pi/2/M) ## slope to test
  sigma = diag(cov)

  ## greedy solution
  grx = matrix(0,1,length(mu))
  Tes2=Tes
  sigs = matrix(0,1,M+1) ## sigma_s history

  if (Tes<=sum(sort(mu,decreasing=TRUE)[1:N])){ ## feasible solution
    for (i in 1:N){
      tmpmu = mu
      tmpmu[mu<Tes2/(N+1-i) | grx==1]=0
      idx = which.max(tmpmu/sqrt(sigma))
      grx[idx] = 1
      Tes2 = Tes2 - mu[idx]
    }
    grv = (Tes - sum(mu*grx))/sqrt(sum(sigma*grx))

    opti = 1e+8
    optx = 0; opta=0

    if (mode==1){ ## diagonal case
      for (i in 1:(M+1)){
        tmp = a[i]*mu - sigma
        n = sum(tmp>0)
        if (N<=n) tmpsol[i,order(tmp,decreasing=TRUE)[1:N]] = 1
        else tmpsol[i,which(tmp>0)] = 1

        tmp2 = (Tes - sum(mu*tmpsol[i,]))/sqrt(sum(sigma*tmpsol[i,])) ## objective function

        if (sum(mu*tmpsol[i,])>=Tes) sigs[i] = sqrt(sum(sigma*tmpsol[i,]))

        if (tmp2<opti){
          opti = tmp2
          optx = tmpsol[i,]
          opta = a[i]
        }

      }
    } else { ## mode==2
      for (i in 1:(M+1)){
        test = solve.subQP(a[i],mu,cov,N)
        tmpsol[i,order(test$solution,decreasing=TRUE)[1:N]] = 1

        tmp2 = (Tes - sum(mu*tmpsol[i,]))/sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))
        if (tmp2<opti){
          opti = tmp2
          optx = tmpsol[i,]
          opta = a[i]
        }
        if (sum(mu*tmpsol[i,])>=Tes) sigs[i] = sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))
      }
    }
  } else { ## infeasible case=> just solve with diagonal assumption
    grx[order(tmpmu/sqrt(sigma),decreasing=TRUE)[1:N]]=1
    grv = (Tes - sum(mu*grx))/sqrt(sum(sigma*grx))

    opti = 1e+8
    optx = 0; opta=0
    sigma = diag(cov)

    for (i in 1:(M+1)){
      tmp = a[i]*mu + sigma
      tmpsol[i,order(tmp,decreasing=TRUE)[1:N]] = 1
      tmp2 = (Tes - sum(mu*tmpsol[i,]))/sqrt(sum(sigma*tmpsol[i,])) ## objective function
      if (tmp2<opti){
        opti = tmp2
        optx = tmpsol[i,]
        opta = a[i]
      }
    }
  }
  list(opti = opti, optx = optx, opta = opta, grx = grx, grv=grv,minsr = min(sigs[1:M]/sigs[2:(M+1)],na.rm=T))
}

## solve the algorithm:mixed with second problem minimizing the penalty
solve.alg2=function(mu,cov,Tes,N,mode=1,M=10,obj=1,greedy=1){
  ## mode=1 :assume the covariance matrix is diagonal
  ## mode=2 :solve subQP
  ## obj=1 : solve SKP
  ## obj=2 : solve the second problem minimizing the penalty
  ## greedy=1: solve greedy algorithm
  tmpsol = matrix(0,M+1,length(mu))
  a = tan((0:M)*pi/2/M) ## slope to test
  #a = tan((0:M)*pi/4/M+pi/4) ## slope to test
  sigma = diag(cov)
  sigs = matrix(0,1,M+1)

  ## greedy solution
  if (greedy==1)  gr = solve.gre(mu,cov,Tes,N,mode=mode,obj=obj)

  if (Tes<=sum(sort(mu,decreasing=TRUE)[1:N])){ ## feasible solution
    opti = 1e+8
    optx = 0; opta=0

    if (mode==1){ ## diagonal case
      for (i in 1:(1+M)){
        tmp = a[i]*mu - sigma
        n = sum(tmp>0)
        if (N<=n) tmpsol[i,order(tmp,decreasing=TRUE)[1:N]] = 1
        else tmpsol[i,which(tmp>0)] = 1

        if (obj==1){
          tmp2 = (Tes - sum(mu*tmpsol[i,]))/sqrt(sum(sigma*tmpsol[i,])) ## objective function
        } else {
          tmp3 = -(Tes - sum(mu*tmpsol[i,]))/sqrt(sum(sigma*tmpsol[i,]))
          tmp2 = sqrt(sum(sigma*tmpsol[i,]))*(1/sqrt(2*pi)*exp(-1/2*tmp3^2)+tmp3*(pnorm(tmp3)-1))
        }
        if (sum(mu*tmpsol[i,])>=Tes) sigs[i] = sqrt(sum(sigma*tmpsol[i,]))

        if (tmp2<opti){
          opti = tmp2
          optx = tmpsol[i,]
          opta = a[i]
        }
      }
    } else { ## mode==2
      for (i in 1:(M+1)){
        test = solve.subQP(a[i],mu,cov,N)
        tmpsol[i,order(test$solution,decreasing=TRUE)[1:N]] = 1

        if (obj==1){
          tmp2 = (Tes - sum(mu*tmpsol[i,]))/sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))
        } else {
          tmp3 = -(Tes - sum(mu*tmpsol[i,]))/sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))
          tmp2 = sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))*(1/sqrt(2*pi)*exp(-1/2*tmp3^2)+tmp3*(pnorm(tmp3)-1))
        }
        if (sum(mu*tmpsol[i,])>=Tes) sigs[i] = sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))
        if (tmp2<opti){
          opti = tmp2
          optx = tmpsol[i,]
          opta = a[i]
        }
      }
    }
  } else { ## infeasible case=> just solve with diagonal assumption

    opti = 1e+8
    optx = 0; opta=0

    for (i in 1:(M+1)){
      tmp = a[i]*mu + sigma
      tmpsol[i,order(tmp,decreasing=TRUE)[1:N]] = 1
      if (mode==1){
        if (obj==1){
          tmp2 = (Tes - sum(mu*tmpsol[i,]))/sqrt(sum(sigma*tmpsol[i,])) ## objective function
        } else {
          tmp3 = -(Tes - sum(mu*tmpsol[i,]))/sqrt(sum(sigma*tmpsol[i,]))
          tmp2 = sqrt(sum(sigma*tmpsol[i,]))*(1/sqrt(2*pi)*exp(-1/2*tmp3^2)+tmp3*(pnorm(tmp3)-1))
        }
      } else {
        if (obj==1){
          tmp2 = (Tes - sum(mu*tmpsol[i,]))/sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))
        } else {
          tmp3 = -(Tes - sum(mu*tmpsol[i,]))/sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))
          tmp2 = sqrt(matrix(tmpsol[i,],nrow=1)%*%cov%*%matrix(tmpsol[i,]))*(1/sqrt(2*pi)*exp(-1/2*tmp3^2)+tmp3*(pnorm(tmp3)-1))
        }
      }

      if (tmp2<opti){
        opti = tmp2
        optx = tmpsol[i,]
        opta = a[i]
      }
    }
  }
  if (greedy==1)  list(opti = opti, optx = optx, opta = opta, grx = gr$grx, grv=gr$grv,minsr = min(sigs[1:M]/sigs[2:(M+1)],na.rm=T),sr=sigs)
  else list(opti = opti, optx = optx, opta = opta,minsr = min(sigs[1:M]/sigs[2:(M+1)],na.rm=T),sr=sigs)
}

## solve the algorithm:mixed with second problem minimizing the penalty
solve.gre=function(mu,cov,Tes,N,mode=1,obj=1){
  ## mode=1 :assume the covariance matrix is diagonal
  ## mode=2 :non diagonal covariance matrix
  ## obj=1 : solve SKP
  ## obj=2 : solve the second problem minimizing the penalty
  sigma = diag(cov)

  ## greedy solution
  grx = matrix(0,1,length(mu))
  Tes2=Tes
  tmpmu = mu
  if (Tes<=sum(sort(mu,decreasing=TRUE)[1:N])){ ## feasible solution
    for (i in 1:N){
      tmpmu = mu
      tmpmu[mu<Tes2/(N+1-i) | grx==1]=0
      idx = which.max(tmpmu/sqrt(sigma))
      grx[idx] = 1
      Tes2 = Tes2 - mu[idx]
    }

    if (mode==1) {
      if (obj==1) grv = (Tes - sum(mu*grx))/sqrt(sum(sigma*grx))
      else {
        tmp3 = -(Tes - sum(mu*grx))/sqrt(sum(sigma*grx))
        grv = sqrt(sum(sigma*grx))*(1/sqrt(2*pi)*exp(-1/2*tmp3^2)+tmp3*(pnorm(tmp3)-1))
      }
    } else { ## cov case
      if (obj==1) grv = (Tes - sum(mu*grx))/sqrt(matrix(grx,nrow=1)%*%cov%*%matrix(grx))
      else {
        tmp3 = -(Tes - sum(mu*grx))/sqrt(matrix(grx,nrow=1)%*%cov%*%matrix(grx))
        grv = sqrt(matrix(grx,nrow=1)%*%cov%*%matrix(grx))*(1/sqrt(2*pi)*exp(-1/2*tmp3^2)+tmp3*(pnorm(tmp3)-1))
      }
    }
  } else { ## infeasible case=> just solve with diagonal assumption
    grx[order(mu/sqrt(sigma),decreasing=TRUE)[1:N]]=1

    if (mode==1) {
      if (obj==1) grv = (Tes - sum(mu*grx))/sqrt(sum(sigma*grx))
      else {
        tmp3 = -(Tes - sum(mu*grx))/sqrt(sum(sigma*grx))
        grv = sqrt(sum(sigma*grx))*(1/sqrt(2*pi)*exp(-1/2*tmp3^2)+tmp3*(pnorm(tmp3)-1))
      }
    } else { ## cov case
      if (obj==1) grv = (Tes - sum(mu*grx))/sqrt(matrix(grx,nrow=1)%*%cov%*%matrix(grx))
      else {
        tmp3 = -(Tes - sum(mu*grx))/sqrt(matrix(grx,nrow=1)%*%cov%*%matrix(grx))
        grv = sqrt(matrix(grx,nrow=1)%*%cov%*%matrix(grx))*(1/sqrt(2*pi)*exp(-1/2*tmp3^2)+tmp3*(pnorm(tmp3)-1))
      }
    }
  }
  list(grx = grx, grv=grv)
}

## calculate the distance among groups: each group is represented by a lifestyle vector (probability distribution vector)
ls.gr.distance=function(dist.mat,lsv1,lsv2){
  ## dist.mat: distance matrix btw lifestyle components
  ## calculate the distance from lsv1 to lsv2 and vice versa. then average them
  tmp.lsv1 = lsv1;  tmp.lsv2 = lsv2;  d1 = 0;  d2 = 0

  for (i in order(lsv1,decreasing=T)){
    tmp = lsv1[i]
    while (tmp > 1e-5){
      for (j in order(dist.mat[i,])){
        if (tmp.lsv2[j]>0) {
          if (tmp.lsv2[j] >= tmp) {
            d1 = d1 + dist.mat[i,j]*tmp
            tmp.lsv2[i] = tmp.lsv2[i]-tmp
            tmp = 0
          } else {
            d1 = d1 + dist.mat[i,j]*tmp.lsv2[j]
            tmp = tmp - tmp.lsv2[j]
            tmp.lsv2[j] = 0
          }
        }
      }
    }
  }

  for (i in order(lsv2,decreasing=T)){
    tmp = lsv2[i]
    while (tmp > 1e-5){
      for (j in order(dist.mat[i,])){
        if (tmp.lsv1[j]>0) {
          if (tmp.lsv1[j] >= tmp) {
            d2 = d2 + dist.mat[i,j]*tmp
            tmp.lsv1[i] = tmp.lsv1[i]-tmp
            tmp = 0
          } else {
            d2 = d2 + dist.mat[i,j]*tmp.lsv1[j]
            tmp = tmp - tmp.lsv1[j]
            tmp.lsv1[j] = 0
          }
        }
      }
    }
  }
  (d1+d2)/2
}

## calculate the distance among groups: each group is represented by a lifestyle vector (probability distribution vector) and the distance is calculated by a heuristic to estimate EMD
ls.gr.distance2=function(dist.mat,lsv1,lsv2){
  ## dist.mat: distance matrix btw lifestyle components
  ## calculate the distance from lsv1 to lsv2 and vice versa. then average them
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

## Normal mixture for daily consumption distribution ##
normal.mixture = function(totale, n=2){
  ## n: number of mixtures
  require('mixtools')
  totaleln = log(as.vector(totale)+1)
  totmixmdl = normalmixEM(totaleln,maxit=100,k=n)
  datRange = seq(min(totaleln),max(totaleln),0.1)
  if (n!=2 & n!=3) print('n should be 2 or 3')
  else {
    if (n==2) {
      learnD = totmixmdl$lambda[1]*dnorm(datRange,mean= totmixmdl$mu[1], sd= totmixmdl$sigma[1])+
      totmixmdl$lambda[2]*dnorm(datRange,mean= totmixmdl$mu[2], sd= totmixmdl$sigma[2])
    } else {
      learnD = totmixmdl$lambda[1]*dnorm(datRange,mean= totmixmdl$mu[1], sd= totmixmdl$sigma[1])+
      totmixmdl$lambda[2]*dnorm(datRange,mean= totmixmdl$mu[2], sd= totmixmdl$sigma[2])+
      totmixmdl$lambda[3]*dnorm(datRange,mean= totmixmdl$mu[3], sd= totmixmdl$sigma[3])
    }

    par(mfrow=c(2,1))
    par(mar=c(4.2,4.8,4,1.1))
    plot(datRange,learnD,main='log(1+daily consumption) distribution fit with mixture gaussian',
         cex.lab=1.2,cex.axis=1.2,type='o',lwd=2,col=2, ylab='Probability density',xlab='log(1+daily consumption)')
    lines(density(totaleln),lwd=2)
    grid()
    legend('topleft',c('Actual distribution','Fitted distribution'),col=1:2,lwd=2,bty='n')

    plot(totmixmdl,which=2,breaks=20,xlab2='log(1+daily consumption)',cex.lab=1.2,cex.axis=1.2,main2='Density curves')
    if (n==2) legend('topleft',c('Actual distribution','Gaussian component 1','Gaussian component 2'),col=1:3,lwd=2,bty='n')
    else legend('topleft',c('Actual distribution','Gaussian component 1','Gaussian component 2','Gaussian component 3'),col=1:4,lwd=2,bty='n')
  }
}

## two sample t-test
two.t.test = function(a,b,d,conf=0.05){
## a,b : encoded vector / d: dimension of codes / conf: confidence level by default 5%
## run two sample t-test and find significantly frequent load shapes in each group
  res = matrix(0,2,d)
  ten1 = table(a)
  ten2 = table(b)
  res[1,as.numeric(names(ten1))] = as.numeric(ten1)
  res[2,as.numeric(names(ten2))] = as.numeric(ten2)

  n1 = length(as.vector(a)); n2 = length(as.vector(b))
  distA = res[1,]/n1;  distB = res[2,]/n2

  s1 = distA*(1-distA);    s2 = distB*(1-distB)
  tstats = (distA-distB)/sqrt(1/n1+1/n2)/sqrt(((n1-1)*s1+(n2-1)*s2)/(n1+n2-2))

  tmp = qt(1-conf/2,n1+n2-2)
  wh1 = which(tstats>=tmp)
  wh2 = which(tstats<=(-tmp))

  list(freq_in_a=wh1,freq_in_b=wh2,tstats=tstats)
}

## function to draw top n shapes from encoded result
top.n.shapes = function(en,en.s=0,dic, n=4, show.avg=0){
## en: encoded data / en.s: encoded daily sum / dic: dictionary
## n: number of top load shapes / show.avg: show avg consumption per load shape
  en = as.vector(en)
  en.s = as.vector(en.s)
  ten = table(en)
  top.n = sort(ten,decreasing=T)[1:min(n,length(ten))]
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
cal.emd1 = function(a,b){
  l = length(a)
  emd = matrix(0,l+1,1)
  for (i in 1:l){
    emd[i+1] = a[i]+emd[i]-b[i]
  }
  sum(abs(emd))
}

## calculate the earth mover's distance btw two dictionaries
dicdist2 = function(A,B,pa = 1,pb = 1, emd=1, same=0,tmpdist=NULL){
  ## dictionary A,B: n1,p and n2,p matrix
  ## pa: probability distribution vector of dictionary A
  ## pb: probability distribution vector of dictionary B
  ## emd: whether to estimate the distance btw two load shapes as EMD or Euclidean
  ## same: whether A and B are same dictionaries or not
  ## tmpdist: user can provide the distance matrix (n1 by n2) btw the dictionaries A and B if already calculated

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
      } else { ## Euclidean
        for (i in 1:(n1-1)){
          a = A[i,]
          for (k in (i+1):n1){
            b = A[k,]
            tmpdist[i,k] = sum(abs(a-b))/2
            tmpdist[k,i] = tmpdist[i,k]
          }
        }
      }
    } else {
      if (emd) {
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
      } else {
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
  d
}

## generate create SQL for given dataframe
createSQL = function(sdata,tname,cnames=NULL,fname=NULL,date.col=NULL,datetime.col=NULL){ ## sdata: source data (data frame)
  createtxt = paste('CREATE TABLE',tname,'(\n')
  if (is.null(cnames)) cnames = names(sdata)

  nr = nrow(sdata); nc = ncol(sdata)

  for (i in 1:nc){
    createtxt = paste(createtxt,'\t',cnames[i])

    if (i%in%date.col) {
      createtxt = paste(createtxt,'date,\n')
    }
    else if (i%in%datetime.col) {
      createtxt = paste(createtxt,'datetime,\n')
    }
    else if (is.numeric(sdata[,i])) { ## use int or float
      if (is.integer(sdata[,i])) { ## use int
        createtxt = paste(createtxt,'int,\n')
      }
      else { ## use float
        createtxt = paste(createtxt,'float,\n')
      }
    }
    else {## use char or varchar
      num = nchar(as.character(sdata[,i]))
      minlen = min(num)
      maxlen = max(num)
      if (maxlen==minlen) { ## use char
        createtxt = paste(createtxt,'char(',minlen,'),\n')
      }
      else { ## use varchar
        createtxt = paste(createtxt,'varchar(',maxlen,'),\n')
      }
    }
  }

  len = nchar(as.character(createtxt))
  ## delete the last ,\n
  createtxt = substr(createtxt,1,len-2)

  createtxt = paste(createtxt,'\n);')
  if (!is.null(fname)) write(createtxt,file=fname)

  createtxt
}


## calculate the distances among load shapes within a dictionary with shift flavor
dicdist3 = function(dic,shift.allow=0,d.metric=1){
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
      if (d.metric==3) dist.mat[i,j] = cal.emd1(dic[i,],dic[j,])
      else if (d.metric==2) dist.mat[i,j] = sum((dic[i,]-dic[j,])^2)
      else dist.mat[i,j] = sum(abs(dic[i,]-dic[j,]))
      dist.mat[j,i] = dist.mat[i,j]
      if (shift.allow) {
        if (d.metric==3) {
          #dist.mat2[i,j] = cal.emd1(dic[i,1:23],dic[j,2:24])
          #dist.mat3[i,j] = cal.emd1(dic[i,2:24],dic[j,1:23])
          dist.mat2[i,j] = cal.emd1(dic[i,],dic[j,c(2:24,1)])
          dist.mat3[i,j] = cal.emd1(dic[i,c(2:24,1)],dic[j,])
        }
        else if (d.metric==2) {
          #dist.mat2[i,j] = sum((dic[i,1:23]-dic[j,2:24])^2)
          #dist.mat3[i,j] = sum((dic[i,2:24]-dic[j,1:23])^2)
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
    dist.mat4
  } else dist.mat
}

raw2target = function(rawdata, is.clean = F, id.col = 1, date.col = 2, zip.col = 3, uidx=4:27,
                      s.date='2011-05-01', e.date='2011-07-31', use.weekends=0, cand.tref=68:86, tch=3, M = 20,
                      N = 1000, Tes = 1000, target.obj=1, target.hour = 18, target.mode=1, target.greedy=0,
                      find.N.for.Tes = 0, p = 0.95, get.prob.Curve = 0, ss=2, response.provide=F,mu=NULL, sigma=NULL){
  ## for the time being, let's use pre-processed temperature data from tempinfo_for_1houraligned.RData
  ## if the date or zipcode is not covered by this, we may need to get temperature info by another function

  ## rawdata: assume it's a data.frame which has sp_id, per_id, date and zip5 information
  ##          or some other id column and date and zip5
  ## is.clean: whether to do cleansing or not
  ## id.col: column idx to be used to identify the customer
  ## date.col: date column idx
  ## zip.col: zip code information column idx
  ## uidx: column idx for usage information (assume it's 24 columns = hourly data)
  ## s.date: start date of response modeling
  ## e.date: end date of response modeling
  ## use.weekends: whether to include weekend data in training response modeling
  ## cand.tref: reference temperature candidates to be used in response modeling
  ## tch: indoor temperature setpoint change
  ## M: number of iteration to run the heuristic algorithm to select customers
  ## N: number of customer enrollment limit
  ## Tes: Targeted energy saving
  ## target.obj: targeting objective, transferred to solve.alg2 as obj
  ## target.hour: targeting hour, 18 means 5~6PM
  ## target.mode: targeting mode, transferred to solve.alg2 as mode
  ## target.greedy: solve targeting by greedy algorithm, transferred to solve.alg2 as greedy
    ## mode=1 :assume the covariance matrix is diagonal
    ## mode=2 :solve subQP
    ## obj=1 : solve SKP
    ## obj=2 : solve the second problem minimizing the penalty
  ## find.N.for.Tes: find N to achieve Tes with probability > p
  ## p: probability to achieve Tes
  ## get.prob.Curve: find the probability curve increasing with increasing N (from prob 0.01 to 0.99)
  ## ss: step size

  if (!response.provide) { ## if the response parameters are provided, we don't have to calculate response parameters

    ## 1. cleanse the data first
    if (!is.clean) {
      print('Cleansing started')
      rawdata = kjs.impute(rawdata,uidx)
      print('Cleansing finished')
    }

    ## 2. get the temperature for the given period
    ## if this is not enough, need to get temperature info by another function
    print('Temperature info load...')
    load('tempinfo_for_1houraligned.RData') ## will load aligned_tempinfo (first column is zip5, next 8760 = 365*24)
    allzips = unique(rawdata[,zip.col])
    zips.to.search = c()
    if (as.Date(s.date)<as.Date('2010-08-01') | as.Date(e.date)>as.Date('2011-07-31')) {
      zips.to.search = allzips
    } else { ## date range is within the period
      zips.to.search = allzips[!(allzips%in%aligned_tempinfo[,1])]
    }

    if (length(zips.to.search)!=length(allzips)) {
      sidx = 24*as.numeric(as.Date(s.date)-as.Date('2010-08-01'))+1
      eidx = 24*as.numeric(as.Date(e.date)-as.Date('2010-08-01')+1)
    }
    periodlen = as.numeric(as.Date(e.date) - as.Date(s.date))+1
    print('Temperature info loaded')

    ## 3. group data per customer and order by date
    sort.data = rawdata[order(rawdata[,id.col], rawdata[,date.col]),]

    ## 4. estimate the parameters
    ## assume most of customers are valid candidates=> estimate the parameters directly
    candidates = c()
    temp.sense.coef = c()
    temp.sense.std = c()
    cnt = 0
    print('Response parameter estimation started')
    for (i in unique(sort.data[,id.col])){
      zip = unique(sort.data[sort.data[,id.col]==i,zip.col])
      udata = sort.data[sort.data[,id.col]==i &
        as.Date(sort.data[,date.col])>=as.Date(s.date) &
        as.Date(sort.data[,date.col])<=as.Date(e.date), uidx] ## n by 24 matrix
      if (nrow(udata)<periodlen) next ## data is not enough to cover the period
      cnt = cnt + 1
      candidates[cnt] = i

      ## get temp info
      temp.info = t(matrix(aligned_tempinfo[aligned_tempinfo[,1]==zip,sidx:eidx+1],nrow=24)) ## n by 24 matrix

      ## temperature response
      if (!use.weekends) {
        didx = get.day(as.Date(s.date)+0:(periodlen-1))
        tres = fit.m2(udata[didx<6,],temp.info[didx<6,],cand.tref)
      } else tres = fit.m2(udata,temp.info,cand.tref)

      temp.sense.coef[cnt] = tres$coefs[target.hour,2]
      temp.sense.std[cnt] = tres$std[target.hour]
    }
    mu = temp.sense.coef*tch
    sigma = diag(temp.sense.std*tch)
    print('Response parameter estimation finished')
  }

  ## 4. select customers and return the solved result
  ## maybe create another function only for this because it may be called seperately a lot
  ## find.N.for.Tes: find N to achieve Tes with probability > p
  ## p: probability to achieve Tes
  ## get.prob.Curve: find the probability curve increasing with increasing N

  print('Targeting started')
  target.res = solve.alg2(mu = mu, cov = sigma^2,
                      Tes=Tes, N=N, mode=target.mode, M=M, obj=target.obj, greedy = target.greedy)
  print('Targeting done')

  if (get.prob.Curve) { ## only for obj 1
    print('Get probability curve')
    Ns = c()
    alg.sol = c(); alg.sol.value = c()
    gre.sol = c(); gre.sol.value = c()

    NT = min(which(cumsum(sort(mu,decreasing=T))>Tes))
    i = NT
    while (i<=length(mu)){
      Ns = c(Ns, i)
      tmp = solve.alg2(mu,sigma^2,Tes=Tes, N=i, mode=target.mode, M=M, obj=target.obj, greedy = target.greedy)

      alg.sol = rbind(alg.sol,tmp$optx)
      alg.sol.value = c(alg.sol.value,1 - pnorm(tmp$opti))

      if (target.greedy==1) {
        gre.sol = rbind(gre.sol,tmp$grx)
        gre.sol.value = c(gre.sol.value,1 - pnorm(tmp$grv))
      }

      if (pnorm(tmp$opti)<0.01) break
      i = i+ss
    }
    i = NT - ss
    while (i>0){
      Ns = c(i, Ns)
      tmp = solve.alg2(mu,sigma^2,Tes=Tes, N=i, mode=target.mode, M=M, obj=target.obj, greedy = target.greedy)

      alg.sol = rbind(tmp$optx,alg.sol)
      alg.sol.value = c(1 - pnorm(tmp$opti),alg.sol.value)

      if (target.greedy==1) {
        gre.sol = rbind(gre.sol,tmp$grx)
        gre.sol.value = c(gre.sol.value,1 - pnorm(tmp$grv))
      }

      if (pnorm(tmp$opti)>0.99) break
      i = i-ss
    }
    print('Get probability curve done')
  }

  N.for.Tes = 0
  Cus.for.Tes.alg = c()
  Cus.for.Tes.gre = c()

  if (find.N.for.Tes) { ## assume p > 0.5 because it doesn't make sense if it's <0.5
    print('Start finding N for Tes ')
    if (get.prob.Curve) {
      idx = which(alg.sol.value>p)
      N.for.Tes = Ns[idx]
      Cus.for.Tes.alg = alg.sol[idx,]
      if (target.greedy==1) Cus.for.Tes.gre = gre.sol[idx,]
    } else {
      NT = min(which(cumsum(sort(mu,decreasing=T))>Tes))
      i = NT
      while (i<=length(mu)){
        tmp = solve.alg2(mu,sigma^2,Tes=Tes, N=i, mode=target.mode, M=M, obj=target.obj, greedy = target.greedy)

        if (pnorm(tmp$opti)<1-p) break
        i = i+ss
      }
      N.for.Tes = i
      Cus.for.Tes.alg = tmp$optx
      if (target.greedy==1) Cus.for.Tes.gre = tmp$grx
    }
    print('Finding N for Tes done')
  }

  if (!response.provide) {
    if (target.greedy) {
      list(N.for.Tes = N.for.Tes, Cus.for.Tes.alg= Cus.for.Tes.alg, Cus.for.Tes.gre = Cus.for.Tes.gre,
         Ns = Ns, alg.sol = alg.sol, alg.sol.value = alg.sol.value, gre.sol = gre.sol, gre.sol.value = gre.sol.value,
         target.prob = 1 - pnorm(target.res$opti), target.cus.alg = target.res$optx, target.cus.gre = target.res$grx,
           candidates = candidates, temp.sense.coef = temp.sense.coef, temp.sense.std = temp.sense.std)
    } else {
      list(N.for.Tes = N.for.Tes, Cus.for.Tes.alg= Cus.for.Tes.alg,
         Ns = Ns, alg.sol = alg.sol, alg.sol.value = alg.sol.value,
         target.prob = 1 - pnorm(target.res$opti), target.cus.alg = target.res$optx,
           candidates = candidates, temp.sense.coef = temp.sense.coef, temp.sense.std = temp.sense.std)
    }
  } else {
    if (target.greedy) {
      list(N.for.Tes = N.for.Tes, Cus.for.Tes.alg= Cus.for.Tes.alg, Cus.for.Tes.gre = Cus.for.Tes.gre,
         Ns = Ns, alg.sol = alg.sol, alg.sol.value = alg.sol.value, gre.sol = gre.sol, gre.sol.value = gre.sol.value,
         target.prob = 1 - pnorm(target.res$opti), target.cus.alg = target.res$optx, target.cus.gre = target.res$grx)
    } else {
      list(N.for.Tes = N.for.Tes, Cus.for.Tes.alg= Cus.for.Tes.alg,
         Ns = Ns, alg.sol = alg.sol, alg.sol.value = alg.sol.value,
         target.prob = 1 - pnorm(target.res$opti), target.cus.alg = target.res$optx )
    }
  }

}

plot.usage = function(toplot, normal=F, metric=1, add=F, lwd=1, col=1, type='o',main=NULL){
  ## toplot is matrix or data.frame of pure consumption data
  ## normal: whether to normalize or not
  ## metric: 1 : normalize by L1, 2 : normalize by L2
  if (normal) toplot = quick.norm(toplot, mod=ifelse(metric==1, 2, 1))

  if (add) {
    for (i in 1:nrow(toplot)) lines(toplot[i,],type=type,lwd=lwd,col=col)
  } else {
    if (normal) plot(toplot[1,], xlab=ifelse(ncol(toplot)==24,'Hour','Index'), ylab='Norm.Usage',ylim=c(min(toplot),max(toplot)),type=type,cex.lab=1.2,cex.axis=1.2,lwd=lwd,col=col,main=main)
    else plot(toplot[1,], xlab=ifelse(ncol(toplot)==24,'Hour','Index'), ylab='Usage (kWh)',ylim=c(min(toplot),max(toplot)),type=type,cex.lab=1.2,cex.axis=1.2,lwd=lwd,col=col,main=main)

    for (i in 2:nrow(toplot)) lines(toplot[i,],type=type,lwd=lwd,col=col)
  }

}
