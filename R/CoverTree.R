#' Cover Tree
#'
#' Implementation of a Cover Tree algorithm for a specified set of data
#'
#' @details The implementation here is based on the rules and insertion algorithm defined in
#' \href{http://hunch.net/~jl/projects/cover_tree/paper/paper.pdf}{Alina Beygelzimer, Sham Kakade, and John Langford: Cover Trees for Nearest Neighbor. (Web)}
#' The current implementation, however, negates the levels outlined in that algorithm.
#' \itemize{
#' \item{All nodes in level i are implicitly in level i+1}
#' \item{All pairs of nodes in level i have a separation of more than 1/(2^i)}
#' \item{All nodes in level i+1 have exactly one parent in level i that is no farther away than 1/(2^i)}
#' }
#'
#' The provided data must be a data.frame with at least two rows of unique data.
#' For the purpose of the CoverTree algorithm, unique data means having a non-zero separation.
#' \strong{Redundant data will not be added to the cover tree}.
#'
#' The provided distance function \strong{must} take two rows from the input data as
#' arguments and must return a distance.   The returned distance must:
#' \itemize{
#' \item{be non-negative}
#' \item{be symmetric: \emph{dist(a,b) = dist(b,a)}}
#' \item{obey the triangle rule: \emph{dist(a,b) + dist(b,c) >= dist(a,c)}}
#' }
#'
#' The provided distance function, \strong{may} takeoptional parameters.
#' If it does, these parameter must also be provided to cover tree constructor in the format recognized by the distance function.
#' These parameter will be passed to every invocation of the distance function.
#'
#' The columns defined in the data.frame are irrelevant to the CoverTree algorithm, but must be
#' consistent with the distance function.
#'
#' @param data a data.frame containing the data to be covered
#' @param dist.func a function that computes the distance between 2 entry rows in the data
#' @return An initialized CoverTree as described above
#' @usage ct <- CoverTree(data,dist.func,...)
#' @export
CoverTree <- function(data,distfunc,...)
{
  self <- new.env()

  if(is.data.frame(data) == FALSE)
  {
    stop(paste("Data must be a data.frame, not a",class(data)))
  }
  if(is.function(distfunc) == FALSE)
  {
    stop(paste("dist.func must be a function, not a",class(dist.func)))
  }
  if(nrow(data)<2)
  {
    stop(paste("A CoverTree requires at least two data points, only",nrow(data),"provided"))
  }

  nodeCount <- nrow(data)
  data      <- data

  if( length(list(...)) > 0 )
  {
    dist.func <- function(a,b) { return(distfunc(a,b,...)) }
  }
  else
  {
    dist.func <- distfunc
  }

  root <- CoverTreeNode(1,data[1,])

  for( row in 2:nodeCount )
  {
    addRow(self,row,data[row,])
  }

  class(self) <- append( class(self), 'CoverTree' )

  return(self)
}

#' @export
add.data <- function(self,data)
{
  UseMethod('add.data',self)
}

#' Add Data
#'
#' Inserts additional data into the cover tree.
#'
#' @details
#' The column names in the input data frame must exactly match those in the data used to create
#' the CoverTree object.
#'
#' Redundant data will not be added to the cover tree.
#'
#' @param ct The CoverTree to be displayed
#' @param data a data.frame containing the data to be added to the cover tree
#' @return n/a
#' @usage  add.data(ct,data)
#' @export
add.data.CoverTree <- function(self,data)
{
  if(is.data.frame(data) == FALSE)
  {
    stop(paste("Data must be a data.frame, not a",class(data)))
  }
  if( all(names(data)==names(self$data)) == FALSE )
  {
    stop(paste("Data must match existing data in CoverTree: ",paste(names(self$root$data),collapse=", ")))
  }
  n.old <- nrow(self$data)
  n.new <- nrow(data)

  rows.new <- (n.old+1):(n.old+n.new)

  self$data <- rbind(self$data,data)

  rownames(data) <- rows.new
  for( row in rows.new )
  {
    addRow(self, row, self$data[row,])
  }
}

addRow <- function(self,row,data)
{
  root <- self$root

  if(is.null(self$param))
  {
    rootDist <- self$dist.func(data,root$data)
  }
  else
  {
    rootDist <- self$dist.func(data,root$data,self$param)
  }

  if(rootDist > 0)
  {
    r <- row
    if( is.null(root$level) )
    {
      c<-addChild(root,2,data,rootDist)
      self$min.level = root$level
      self$max.level = root$level
    }
    else if( insert(self,r,data,list(root),rootDist,root$level) == FALSE )
    {
      c<-addChild(root,r,data,rootDist)
      if( c$level    > self$max.level)  { self$max.level = c$level    }
      if( root$level < self$min.level ) { self$min.level = root$level }
    }
  }
}

insert <- function(self,row,data,Qi.nodes,Qi.dists,level)
{
  sep <- 1./(2^level)

  candQi.node <- NULL
  candQi.dist <- NULL

  Qj.nodes <- list()
  Qj.dists <- NULL

  nQi <- length(Qi.dists)
  for( qi in 1:nQi )
  {
    qi.node <- Qi.nodes[[qi]]
    qi.dist <- Qi.dists[qi]
    if( qi.dist <= sep )
    {
      if( is.null(candQi.dist) || (qi.dist < candQi.dist) )
      {
        candQi.node <- qi.node
        candQi.dist <- qi.dist
      }

      Qj.nodes <- append(Qj.nodes, qi.node)
      Qj.dists <- append(Qj.dists, qi.dist)
    }

    children <- getChildren(qi.node,level+1)
    if( is.null(children) == FALSE )
    {
      for(q in children)
      {
        if(is.null(self$param))
        {
          dist <- self$dist.func(data,q$data)
        }
        else
        {
          dist <- self$dist.func(data,q$data,self$param)
        }
        if( dist == 0.0 ) { return(TRUE) }

        if(dist <= sep)
        {
          Qj.nodes <- append(Qj.nodes, q)
          Qj.dists <- append(Qj.dists, dist)
        }
      }
    }
  }

  nQj <- length(Qj.dists)
  if( nQj == 0 ) { return(FALSE) }

  if( insert(self,row,data,Qj.nodes,Qj.dists,level+1) ) { return(TRUE) }
  if( is.null(candQi.node) ) { return(FALSE) }

  q <- addChild(candQi.node, row, data, candQi.dist)
  if( q$level > self$max.level) { self$max.level = q$level }
  if( self$root$level < self$min.level ) { self$min.level = self$root$level }

  return(TRUE)
}

#' Print a CoverTree
#' @param ct The CoverTree to be displayed
#' @return n/a
#' @usage
#' \code{> print(ct)}
#' \emph{or}
#' \code{> ct}
#' @export
print.CoverTree <- function(self)
{
  cat(sprintf('Min.Level = %3d (root)\nMax.Level = %3d\n',self$min.level,self$max.level))
  print(self$root)
}

#' @export
as.nodes <- function(x)
{
  UseMethod('as.nodes',x)
}

#' CoverTree node list
#'
#' @details Returns a data.frame with the following columns:
#' \itemize{
#'  \item{\strong{level}: level of the current node}
#'  \item{\strong{parent}: index of the row in the data which is the parent node of the current node}
#'  \item{\strong{distance}: distance from the current node to the parent node}
#'  }
#'  The parent node and distance to parent are both set to 0 for the root node
#'
#' @param ct The CoverTree containing the nodes
#' @return data.frame with columns: parent, level, distance
#' @export
as.nodes.CoverTree <- function(self)
{
  nodes <- nodeToDataframe(self$root)
  nodes <- nodes[ order(nodes$row), ]

  levels <- unique(nodes$level)
  levels <- levels[ order(levels) ]

  counts <- sapply(levels,function(t) { sum(nodes$level>=t) })
  nodes <- merge(nodes, data.frame(level=levels, n.thresh=counts))

  nodes$score <- (1- (nodes$n.thresh - nodes$n.cluster)/self$root$n.cluster ) * (nodes$n.cluster>1)

  return(nodes)
}

#' @export
as.dendrogram <- function(self,labels=NULL)
{
  UseMethod('as.dendrogram')
}

#' Dendrogram data
#'
#' @details Returns a data structure that can be plotted as hclust dendrograms
#' @param ct CoverTree to be converted to be converted
#' @param labels [OPTIONAL] Labels applied to each node (indexed by row number in covered data)
#' @return data structure that can be plotted as a dendrogram
#' @usage
#' \code{> dend <- as.dendrogram(ct[,labels])}
#' \code{> plot(dend) }
#' @export
as.dendrogram.CoverTree <- function(self,labels=NULL)
{
  dend <- new.env()
  dend$count  <- 0
  dend$merge  <- matrix(nrow=nrow(self$data)-1,ncol=2)
  dend$order  <- NULL
  dend$height <- NULL
  dend$labels <- NULL

  nodes    <- new.env()
  clusters <- new.env()

  addToDendrogram(self$root,dend,nodes,clusters,labels)

  dend$merge <- matrix( dend$merge[ ! is.na(dend$merge) ], ncol=2 )

  class(dend) <- c( class(dend), 'CoverTreeDendrogram' )
  return(dend)
}

#' Plot CoverTree Dendrogram
#'
#' @details Returns a data structure that can be plotted as hclust dendrograms
#'
#' Note that the CoverTreeDendrogram does \strong{not} provide data labels
#' If you wish to have the nodes labeled, you must append a labels dimension
#' to the CoverTreeDendrogram object as shown above
#'
#' This plot function takes all of the same optional arguments taken by
#' plot.hclust.  See the link below for additional info.
#'
#' @param dend A cover tree dendrogram object
#' @seealso \link{as.dendrogram.CoverTree}
#' @seealso \link{plot.hclust}
#' @return n/a
#' @usage
#' \code{> dend <- as.dendrogram(ct) }
#' \code{> dend$labels <- x }
#' \code{> plot(dend) }
#' @export
plot.CoverTreeDendrogram <- function(dend,
                                     main='CoverTree Dendrogram',
                                     ylab='Log2(distance)',
                                     ... )
{
  hc <- new.env()
  hc$merge  <- dend$merge
  hc$order  <- dend$order
  hc$height <- dend$height
  hc$labels <- dend$labels

  class(hc) <- 'hclust'

  plot( hc, main=main, ylab=ylab, ... )
}


#' @export
split <- function(self,level,prune=0)
{
  cat(sprintf("Generic1 split\n"))
  UseMethod('split',self)
}

#' Split CoverTree
#'
#' @details Returns a list of CoverTrees formed by segmenting the tree
#' such that the root node of each segmented tree is at ggthe specified
#' level.
#'
#' @param ct The cover tree to be split
#' @param level The level at which to cut the tree
#' @param prune How to handle nodes above the split level
#' \itemize{
#'    \item 0 = No nodes above the split level are ever included in the new subtrees
#'    \item 1 = A single node above the split level is included if there is no node in the subtree does not contain a node at the split level
#'    \item 2 = All nodes above the split level are replicated to each new subtree
#' }
#' @usage
#' \code{> subtrees <- split(ct,level,prune) }
#' @export
split.CoverTree <- function(self,level,prune=0)
{
  cat(sprintf("CoverTree split\n"))
  rval <- split(self$root,level,prune)
  return(rval)
}
