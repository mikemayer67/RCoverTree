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
#' The provided distance function must take two rows from the input data as
#' arguments and must return a distance.   The returned distance must:
#' \itemize{
#' \item{be non-negative}
#' \item{be symmetric: \emph{dist(a,b) = dist(b,a)}}
#' \item{obey the triangle rule: \emph{dist(a,b) + dist(b,c) >= dist(a,c)}}
#' }
#' The columns defined in the data.frame are irrelevant to the CoverTree algorithm, but must be
#' consistent with the distance function.
#'
#' @param data a data.frame containing the data to be covered
#' @param dist.func a function that computes the distance between 2 entry rows in the data
#' @return An initialized CoverTree as described above
#' @usage ct <- CoverTree(data,dist.func)
#' @export
CoverTree <- function(data,dist.func)
{
  self  <- environment()

  if(is.data.frame(data) == FALSE)
  {
    stop(paste("Data must be a data.frame, not a",class(data)))
  }
  if(is.function(dist.func) == FALSE)
  {
    stop(paste("dist.func must be a function, not a",class(dist.func)))
  }
  if(nrow(data)<2)
  {
    stop(paste("A CoverTree requires at least two data points, only",nrow(data),"provided"))
  }

  nodeCount <- nrow(data)
  dist.func <- dist.func
  data      <- data

  root <- CoverTreeNode(1,data[1,])

  addChild(root,2,data[2,],dist.func(data[1,],data[2,]))

  for( row in 3:nodeCount )
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

  rootDist = self$dist.func(data,root$data)

  if(rootDist > 0)
  {
    r <- row
    if( insert(self,r,data,list(root),rootDist,root$level) == FALSE )
    {
      addChild(root,r,data,rootDist)
    }
  }
}

insert <- function(self,row,data,Qi.nodes,Qi.dists,level)
{
  sep <- 1./(2^level)

  candQi.node <- NULL
  candQi.dist <- NULL

  Qj.nodes = list()
  Qj.dists = NULL

  nQi <- length(Qi.dists)
  for( qi in 1:nQi )
  {
    qi.node <- Qi.nodes[[qi]]
    qi.dist <- Qi.dists[qi]
    if( qi.dist <= sep )
    {
      if( is.null(candQi.dist) || (qi.dist < candQi.dist) )
      {
        candQi.node = qi.node
        candQi.dist = qi.dist
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

  if( insert(self,row,data,Qj.nodes,Qj.dists,level+1) ) { return(TRUE) }
  if( is.null(candQi.node) ) { return(FALSE) }

  q <- addChild(candQi.node, row, data, candQi.dist)

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
  nodes <- node_to_dataframe(self$root)
  nodes <- nodes[ order(nodes$row), ]
  rownames(nodes) <- nodes$row
  nodes$row <- NULL

  return(nodes)
}

#' @export
as.hclust <- function(self)
{
  UseMethod('as.hclust',self)
}

#' Dendragram data
#'
#' @details Returns a data structure that can be plotted as hclust dendragrams
#' @param ct The CoverTree to be converted to be converted
#' @return data structure that can be plotted as a dendragram
#' @usage
#' \code{> hc <- as.hclust(ct)}
#' \code{> plot(hc) }
#' @export
as.hclust.CoverTree <- function(self)
{
  nodes <- as.nodes(self)
  nodes <- nodes[ order(-nodes$level,nodes$parent), ]
  n     <- nrow(nodes)
  nodes$node    <- as.numeric( rownames(nodes) )
  nodes$cluster <- rep(NA,n)

  for( i in 1:n-1 )
  {
    node   <- nodes$node[i]
    parent <- nodes$parent[i]


  }
}

