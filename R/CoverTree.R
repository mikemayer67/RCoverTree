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

  if( length(list(...)) > 0 )
  {
    self$dist.func <- function(a,b) { return(distfunc(a,b,...)) }
  }
  else
  {
    self$dist.func <- distfunc
  }

  self$root      <- CoverTreeNode(1,data[1,])
  self$nextID    <- 2
  self$n.nodes   <- 1

  for( row in 2:nrow(data) )
  {
    addRow(self,data[row,])
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

  n = nrow(data)
  for( i in 1:rows.new )
  {
    addRow(self, data[i,])
  }
}

addRow <- function(self,data)
{
  root <- self$root

  rootDist <- self$dist.func(data,root$data)

  if(rootDist > 0)
  {
    node <- NULL
    if( is.null(root$level) )
    {
      node<-addChild(root,self$nextID,data,rootDist)
    }
    else if( insert(self,data,list(root),rootDist,root$level) == FALSE )
    {
      node<-addChild(root,self$nextID,data,rootDist)
    }
    if( is.null(node) == FALSE )
    {
      self$nextID  <- 1 + self$nextID
      self$n.nodes <- 1 + self$n.nodes
    }
  }
}

insert <- function(self,data,Qi.nodes,Qi.dists,level)
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
        dist <- self$dist.func(data,q$data)
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

  if( insert(self,data,Qj.nodes,Qj.dists,level+1) ) { return(TRUE) }
  if( is.null(candQi.node) ) { return(FALSE) }

  q <- addChild(candQi.node, self$nextID, data, candQi.dist)
  if( is.null(q) == FALSE )
  {
    self$nextID  <- 1 + self$nextID
    self$n.nodes <- 1 + self$n.nodes
  }

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
  cat(sprintf('Min.Level = %3d (root)\nMax.Level = %3d\n',self$root$level,self$root$max.level))
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

  levels <- unique(nodes$level)
  levels <- levels[ order(levels) ]

  counts <- sapply(levels,function(t) { sum(nodes$level>=t) })
  nodes <- merge(nodes, data.frame(level=levels, n.thresh=counts))

  nodes <- nodes[ order(nodes$id), ]
  rownames(nodes) <- nodes$id

  nodes <- nodes[,c(1,3:ncol(nodes))]

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
  dend$merge  <- matrix(nrow=self$n.nodes-1,ncol=2)
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
  root.nodes <- split(self$root,level,prune)

  rval <- list()

  for( node in root.nodes )
  {
    subtree <- new.env()
    subtree$dist.func <- self$dist.func
    subtree$n.nodes   <- self$n.cluster
    subtree$root      <- node

    class(subtree) <- class(self)

    rval <- append(rval,subtree)
  }

  return(rval)
}
