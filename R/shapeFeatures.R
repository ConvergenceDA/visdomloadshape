# Copyright 2016 The Board of Trustees of the Leland Stanford Junior University.
# Direct inquiries to Sam Borgeson (sborgeson@stanford.edu)
# or professor Ram Rajagopal (ramr@stanford.edu)

#library(class) # for knn
# encoding.dict is the dictionary that will be used to return the cluster, SE, and RSE columns of the
# encodings and can be used to calculate customer shape stats.If it is not provided, the category-mapped
# dict will be used.
#' @export
shapeCategoryEncoding = function(rawData,metaCols=1:4,encoding.dict=NULL) {
  load(file.path(VISDOM_PATH,'dictionary100.RData'))
  category.mapped.dict = merged.center # rename loaded dict variable
  rm(merged.center)

  shapeCategories = c('morning peak','day peak','evening peak','night peak',
                      'double peak (morn,eve)','double peak (eve,night)' ,'double peak (day,eve)',
                      'double peak (morn,day)','double peak (morn,night)','double peak (day,night)')
  shapeCatAbbrev = c('morn','day','eve','night','morn_eve','eve_night','day_eve','morn_day','morn_night','day_night')
  dict.category.mappings = c(  7, 4, 4, 1, 2, 7, 8, 6, 2, 6,
                               9, 2, 6, 7, 3, 3, 3, 5, 2, 4,
                               4, 3, 2, 1, 3, 2, 3, 7, 2, 2,
                               3, 6, 1, 7,10, 4, 3, 8, 2, 1,
                               2, 4, 7, 5, 3, 3,10, 3, 1,10,
                               3, 2, 2, 2, 6, 2, 2, 2, 4, 3,
                               1, 2, 7, 6, 1, 3, 7, 4, 2, 3,
                               1, 3, 2, 7, 4, 3, 5, 3, 3, 1,
                               4, 2, 1, 3, 1, 3, 3, 3, 4, 7,
                               4, 3, 4, 3, 2, 5, 2, 2, 4, 1)

  print('Encoding data using categorized shapes')
  categorical.encodings = data.frame(encode(rawData[,-metaCols],category.mapped.dict))
  names(categorical.encodings) = c('cluster','kWh','SE','RSE')

  if(is.null(encoding.dict)) {
    print('Using categorical dictionary as encoding dictionary')
    encoding.dict = category.mapped.dict
    shape.encodings = categorical.encodings
  } else {
    print('Encoding data using specified encoding dictionary')
    shape.encodings = data.frame(encode(rawData[,-metaCols],encoding.dict))
    names(shape.encodings) = c('cluster','kWh','SE','RSE')
  }
  shape.encodings = cbind(rawData[,metaCols],shape.encodings)

  # add the qualitative categories to the data
  shape.encodings$category = dict.category.mappings[categorical.encodings$cluster]

  # build a map between the encoding dict and the category mapped dict and use it to look up the matching categories
  encoding.dict.categories  = dict.category.mappings[
                                class::knn(category.mapped.dict,encoding.dict, 1:nrow(category.mapped.dict)) ]

  encoding.dict.category.info      = data.frame(cluster=as.numeric(row.names(encoding.dict)),category=encoding.dict.categories,name=shapeCategories[encoding.dict.categories])
  encoding.dict.category.info$name = paste(encoding.dict.category.info$name)


  return( list(
            encoding.dict               = encoding.dict,
            category.mapped.dict        = category.mapped.dict,        # the dictionary whose cluster centers have categorical mappings
            encoding.dict.category.info = encoding.dict.category.info, # mapping between encoding dict clusters and qualitative categories
            encodings                   = shape.encodings,             # dict encodings with category col added
            shapeCategories             = shapeCategories,             # names for the indices in the encodings$category column
            shapeCatAbbrev              = shapeCatAbbrev               # abberviated names for the indices in the encodings$categories column
    ))
}

#' @export
shapeFeatures = function(d,metaCols=1:4) {

  print('Dividing encodings by customer id')
  # encodings broken out per customer
  cust.encodings = dlply(d$encodings[,c('id','cluster','category','kWh')],
                         .(id),
                         .fun=function (x) { return(x[,c('cluster','category','kWh')] ) },
                         .progress='text' )

  shape.features = d$encodings[match(unique(d$encodings$id),d$encodings$id),metaCols] # get the first row of data for every customerto use their metadata
  print('Computing entropy')
  shape.features$entropy = laply(cust.encodings,.fun=function(x) { shannon.entropy2(x$cluster) },.progress='text')

  print('Computing category.counts')
  catCount = length(d$shapeCategories)
  clustCount = nrow(d$encoding.dict)
  shape.stats = list()
  emptyCategories = data.frame(cluster=0,category=1:catCount,kWh=0)
  emptyClusters   = data.frame(cluster=1:clustCount,category=0,kWh=0)
  shape.stats$category.counts = ldply(cust.encodings,
                               .fun=function(x) {
                                 counts = table(c(1:catCount,x$category)) - 1 # add 1:num of clusters so the table fills with all possible values, but then subtract 1 to recover the zeros
                                 return(counts) },
                               .progress='text' )
  colnames(shape.stats$category.counts) = c('id',paste(d$shapeCatAbbrev,'_count',sep=''))
  shape.stats$category.counts$id = paste(shape.stats$category.counts$id)

  print('Computing category.energy')
  shape.stats$category.energy = ldply(cust.encodings,
                               .fun=function(x) {
                                 sums = aggregate(kWh~category, data=rbind(emptyCategories,x), sum)
                                 catEnergy = data.frame(t(sums$kWh) )
                                 names(catEnergy) = sums$category
                                 return(catEnergy) },
                               .progress='text' )
  colnames(shape.stats$category.energy) = c('id',paste(d$shapeCatAbbrev,'_energy',sep=''))
  shape.stats$category.energy$id = paste(shape.stats$category.energy$id)

  print('Computing cluster.counts')
  shape.stats$cluster.counts = ldply(cust.encodings,
                                      .fun=function(x) {
                                        counts = table(c(1:clustCount,x$cluster)) - 1 # add 1:num of clusters so the table fills with all possible values, but then subtract 1 to recover the zeros
                                        return(counts) },
                                      .progress='text' )
  colnames(shape.stats$cluster.counts) = c('id', paste('clust_',colnames(shape.stats$cluster.counts)[-1],'_count',sep=''))
  shape.stats$cluster.counts$id = paste(shape.stats$cluster.counts$id)

  print('Computing cluster.energy')
  shape.stats$cluster.energy = ldply(cust.encodings,
                                      .fun=function(x) {
                                        sums = aggregate(kWh~cluster, data=rbind(emptyClusters,x), sum)
                                        clustEnergy = data.frame(t(sums$kWh) )
                                        names(clustEnergy) = sums$cluster
                                        return(clustEnergy) },
                                      .progress='text' )
  colnames(shape.stats$cluster.energy) = c('id',paste('clust_',colnames(shape.stats$cluster.energy)[-1],'_energy',sep=''))
  shape.stats$cluster.energy$id = paste(shape.stats$cluster.energy$id)

  # return the shape features by customer and the data segmented by customer id for further use
  return( list( shape.features       = shape.features,
                shape.stats          = shape.stats,
                encoding.dict        = d$encoding.dict,
                category.mapped.dict = d$category.mapped.dict,
                encoding.dict.category.info = d$encoding.dict.category.info,
                category.names       = d$shapeCategories,
                category.nms         = d$shapeCatAbbrev
              ) )
}

