# Copyright 2016 The Board of Trustees of the Leland Stanford Junior University.
# Direct inquiries to Sam Borgeson (sborgeson@stanford.edu) 
# or professor Ram Rajagopal (ramr@stanford.edu)

require('akmeans')  # adaptive K-means v1.1 library from CRAN
#vsource('akmeans.R') # bug-fix override required
#vsource('loadShape.R')


#' @export
cleanData = function(rawData,forceSameDuration=F, subtractMins=F, minPower=NULL) {
  tic('clean data')
  badRowIdx = which( rowSums(is.na(rawData))>0 )
  toc('clean data')
  print(paste('removing',length(badRowIdx),'out of',nrow(rawData),'due to NAs'))
  if(length(badRowIdx) != 0) { rawData = rawData[-badRowIdx,] } # removes rows.
  if( ! is.null(minPower)) {
    power = rowMeans(rawData[,-1:-4],na.rm=T)
    bads = power < minPower
    print(paste('Removing',sum(bads),'days for falling under a mean of',minPower,'kW'))
    rawData = rawData[which(! bads),]
  }
  if(subtractMins) {
    mins = apply(rawData[,-1:-4],1,min)
    rawData[,-1:-4] = rawData[,-1:-4] - mins
  }
  if(forceSameDuration) {
    tic('id counts')
    idCounts = table(rawData$id)
    toc('id counts')
    fullCount = Mode(idCounts) # full data for the period in question will have this many rows
    fullDataIds = names(idCounts)[idCounts == fullCount]
    print(paste('Preserving',
                sprintf("%.1f",length(fullDataIds)/length(idCounts)*100),
                '% of the ids with complete data for the time period.'))
    return(rawData[rawData$id %in% fullDataIds,])
  }
  return( rawData )
}

#' @export
runShapesSimple = function(rawData,nSamples=100000,dictSize=NULL) {
  print(dim(rawData))
  twoStep = ! is.null(dictSize)
  if( is.null(dictSize)) { dictSize = 1000 }
  shapeResults = raw2encoded(rawData[,-1:-4],
                             is.clean=T, use.all=F, s.size=nSamples,
                             target.size=dictSize, mode=1, d.metric=1,
                             ths=0.2, iter.max=100, nstart=1,
                             two.step.compress=twoStep,verbose=T)

  shapeResults$clean.data = NA # remove memory hog copy of the data

  #load shape code, daily sum, square err, relative square err
  encodings = cbind(rawData[,1:4],shapeResults$encoded.data)
  names(encodings)[-1:-4] = c('cluster','kWh','SE','RSE')
  dict = data.frame(shapeResults$dictionary)
  return(list(encodings=encodings,dict=dict))
}
