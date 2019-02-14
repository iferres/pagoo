# PgR6M

#' @name PgR6M
#' @title PgR6 class with methods.
#' @description PgR6 with Methods.
#' @importFrom R6 R6Class
#' @importFrom micropan distJaccard distManhattan heaps fluidity binomixEstimate
#' @importFrom ggplot2 ggplot aes geom_bar geom_raster geom_tile theme element_blank scale_fill_grey xlab ylab
#' @importFrom reshape2 melt
#' @export
PgR6M <- R6Class('PgR6M',

                 inherit = PgR6,

                 public = list(

                   initialize = function(cluster_list,
                                         prefix = 'group',
                                         sep = '__'){

                     super$initialize(cluster_list = cluster_list,
                                      prefix = prefix,
                                      sep = sep)

                   },

                   #Analyses Methods
                   dist_jaccard = function(){
                     distJaccard(self$pan_matrix)
                   },

                   dist_manhattan = function(scale = 0,
                                             weights = rep(1, dim(self$pan_matrix)[2])){
                     distManhattan(self$pan_matrix,
                                   scale = scale,
                                   weights = weights)
                   },

                   heaps = function(n.perm=100){
                     #micropan::heaps()
                     heaps(self$pan_matrix ,n.perm = n.perm)
                   },

                   fluidity = function(n.sim = 10){
                     fluidity(self$pan_matrix, n.sim = n.sim)
                   },

                   binomix_estimate = function(K.range = 3:5,
                                               core.detect.prob = 1,
                                               verbose = TRUE){
                     binomixEstimate(pan.matrix = self$pan_matrix,
                                     K.range = K.range,
                                     core.detect.prob = core.detect.prob,
                                     verbose = verbose)
                   },


                   #Plot Methods
                   gg_barplot = function(){
                     dd <- as.data.frame(table(colSums(self$pan_matrix)))
                     ggplot(dd, aes(x=Var1, y=Freq)) +
                       geom_bar(stat='identity') +
                       ylab('Number of genes') +
                       xlab('Number of genomes')
                   },

                   gg_binmap = function(){
                     tpm <- t(self$pan_matrix)
                     tpm[which(tpm > 0, arr.ind = TRUE)] <- 1L
                     bm <- as.data.frame(tpm)
                     or <- order(rowSums(bm), decreasing = TRUE)
                     lvls <- rownames(bm)[or]
                     bm$gene <- factor(rownames(bm), levels = lvls)
                     bm <- melt(bm, 'gene')
                     bm$value <- factor(bm$value, levels = c(1, 0))
                     ggplot(bm, aes(gene, variable, fill=value)) +
                       geom_raster() +
                       theme(axis.ticks.x = element_blank(),
                             axis.text.x = element_blank()) +
                       scale_fill_grey(start = .2, end = .9)
                   },

                   gg_dist = function(dist = 'Jaccard', ...){
                     dist <- match.arg(dist, c('Jaccard', 'Manhattan'))
                     if (dist == 'Jaccard'){
                       m <- as.matrix(self$dist_jaccard())
                     }else{
                       m <- as.matrix(self$dist_manhattan(...))
                     }
                     me <- melt(m)
                     ggplot(me, aes(Var1, Var2, fill=value)) +
                       geom_tile()
                   }


                 )

)
