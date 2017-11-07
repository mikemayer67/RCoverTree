#' Cover Tree Node
#'
#' Constructs a cover tree node.  This constructor should proably never
#' be called directly.   Rather it should be called by CoverTree during
#' its construction and subsequent updates.
#' @param id A unique index identifying the node being constructed within the cover tree
#' @param data A valid R object containing the data represented by the node
#' @param parent The parent node for the one being constructed (n/a for root node)
#' @param distance Either a numeric value or a distance function [ f(data,data) -> numeric ]
#' @return An initialized CoverTreeNode, fully connected as a child of 'parent'
#' @seealso \link{CoverTree}
#' @export
CoverTreeNode <- function(id,data,parent=NULL,distance=NULL)
{
  self <- new.env()

  if(is.null(data))  { stop('The data argument must be provided') }

  self$id       <- id
  self$data     <- data
  self$isRoot   <- is.null(parent)
  self$isLeaf   <- TRUE
  self$children <- new.env()

  self$n.cluster <- 1

  if( self$isRoot )
  {
    if( is.numeric(distance) ) { warning('No need to specify distance value for root node') }
  }
  else  # non-root node
  {
    if( is(parent,'CoverTreeNode') == FALSE )
    {
      stop('Parent must be a CoverTreeNode, not: ',class(parent))
    }
    if( is.null(distance) )
    {
      stop("Must specify distance for non-root node")
    }
    if( is.numeric(distance) == FALSE )
    {
      stop("Distance must be a number, not: ",distance)
    }
    if( distance < 0.0 )
    {
      stop("Distance must be positive, not: ",distance)
    }

    self$parent    <- parent
    self$distance  <- distance
    self$level     <- floor(1-log2(distance))
    self$max.level <- self$level
  }

  class(self) <- append( class(self), 'CoverTreeNode' )

  return(self)
}

cloneNode <- function(node,root=FALSE)
{
  clone <- new.env()

  clone$id        <- node$id
  clone$data      <- node$data
  clone$isRoot    <- root
  clone$isLeaf    <- node$isLeaf

  if( root == FALSE )
  {
    clone$parent   <- node$parent
    clone$distance <- node$distance
  }

  if( is.null(node$level) == FALSE )
  {
    clone$level     <- node$level
    clone$max.level <- node$max.level
  }

  clone$children <- new.env()

  clone$n.cluster <- 1

  if( node$isLeaf == FALSE )
  {
    for( level in ls(node$children) )
    {
      levelChildren <- list()
      for( child in get(level,node$children) )
      {
        child.clone <- cloneNode(child)
        clone$n.cluster = clone$n.cluster + child.clone$n.cluster
        levelChildren <- append(levelChildren,child.clone)
      }
      assign(level,levelChildren,node$children)
    }
  }

  class(clone) <- class(node)
  return(clone)
}

#' @export
addChild <- function(node,id,data,distance)
{
  UseMethod('addChild',node)
}

#' @export
addChild.CoverTreeNode <- function(node,id,data,distance)
{
  if(distance == 0)
  {
    warning("Attempted to add new node ",id,", as a child node of row ",node$id," at distance 0", call.=FALSE)
    return(NULL)
  }
  child <- CoverTreeNode(id,data,node,distance)

  if( is.null(node$level) ) #this must be root node (first child being added)
  {
    node$level     <- child$level - 1
    node$max.level <- child$level
  }
  else if( child$level <= node$level )
  {
    if(node$isRoot) # all is well, this is root node... just move it up to "higher" level
    {
      node$level <- child$level - 1
    }
    else
    {
      stop("Internal error: attempted to add node at level ",child$level,
           " as child to non-root node at ",node$level)
    }
  }

  node$n.cluster <- 1 + node$n.cluster
  p = node$parent
  while( is.null(p) == FALSE )
  {
    p$n.cluster <- 1 + p$n.cluster
    if( child$level > p$max.level ) { p$max.level <- child$level }
    p = p$parent
  }

  level    <- as.character(child$level)
  children <- node$children

  if( exists(level,children) )
  {
    levelChildren <- append(get(level,children), child)
  }
  else
  {
    levelChildren <- list(child)
  }
  assign(level,levelChildren,children)

  node$isLeaf = FALSE

  return(child)
}

#' @export
getChildren <- function(node,level)
{
  UseMethod('getChildren',node)
}

#' @export
getChildren.CoverTreeNode <- function(node,level)
{
  rval <- NULL
  level <- as.character(level)
  if( exists(level,node$children) )
  {
    rval <- get(level,node$children)
  }
  return(rval)
}

#' @export
print.CoverTreeNode <- function(node)
{
  cat("\n")
  printNode(node,'  ')
  cat("\n")
}

#' @export
split.CoverTreeNode <- function(self,level,prune=0)
{
  rval <- list()
  cand.root.levels <- as.numeric(ls(self$children))
  new.root.levels <- cand.root.levels[cand.root.levels >= level]

  for( rlevel in new.root.levels )
  {
    new.roots <- get(as.character(rlevel),self$children)
    for( root in new.roots )
    {
      new.root = cloneNode(root,root=TRUE);
      rval <- append(rval, new.root)
    }
  }
  return(rval)
}

printNode <- function(node,prefix)
{
  if(is.null(node$distance))
  {
    msg <- sprintf("%s[%2d]: %s (ROOT)\n",prefix,node$level,node$id);
  }
  else
  {
    msg <- sprintf("%s[%2d]: %3d (d=%s)\n",prefix,node$level,node$id,node$distance);
  }
  cat(noquote(msg));
  for( level in ls(node$children ) )
  {
    children <- get(level,node$children)
    for( child in children )
    {
      printNode(child,paste(prefix,"  "));
    }
  }
}

#' @export
nodeToDataframe <- function(node)
{
  p <- ifelse(node$isRoot, NA, node$parent$id)
  d <- ifelse(node$isRoot, NA, node$distance)
  nc <- node$n.cluster
  rval <- data.frame( id=node$id, level=node$level, parent=p, distance=d, n.cluster=nc )

  for(level in ls(node$children))
  {
    children <- get(level,node$children)
    for(child in children)
    {
      st <- nodeToDataframe(child)
      rval <- rbind(rval,st)
    }
  }

  return(rval)
}

#' @export
addToDendrogram <- function(node,dend,nodes,clusters,labels)
{
  children <- NULL;
  levels <- ls(node$children)
  for( level in levels )
  {
    levelChildren <- get(level,node$children)
    children <- c(children, levelChildren)
  }

  dists <- as.numeric( lapply(children,function(x) { x$distance }) )
  children <- children[ order(dists) ]

  for( child in children )
  {
    if( child$isLeaf == FALSE )
    {
      addToDendrogram(child,dend,nodes,clusters,labels)
    }

    knr <- as.character(node$id)
    if(exists(knr,nodes)) { a <- get(knr,nodes) } else { a <- 1 + length(ls(nodes)); assign(knr,a,nodes) }

    kcr <- as.character(child$id)
    if(exists(kcr,nodes)) { b <- get(kcr,nodes) } else { b <- 1 + length(ls(nodes)); assign(kcr,b,nodes) }

    ka <- as.character(a)
    while(exists(ka,clusters))
    {
      hasA <- TRUE
      a <- get(ka,clusters)
      ka <- as.character(a)
    }

    kb <- as.character(b)
    while(exists(kb,clusters))
    {
      hasB <- TRUE
      b <- get(kb,clusters)
      kb <- as.character(b)
    }

    dend$count <- dend$count + 1
    dend$merge[ dend$count, 1 ] <- -a
    dend$merge[ dend$count, 2 ] <- -b
    dend$height <- c( dend$height, log2(child$distance))

    if( a > 0 )
    {
      dend$order <- c(dend$order, a)
      if(is.null(labels)) { dend$labels[a] <- as.character(node$id) }
      else                { dend$labels[a] <- labels[node$id] }
    }
    if( b > 0 )
    {
      dend$order <- c(dend$order, b)
      if(is.null(labels)) { dend$labels[b] <- as.character(child$id) }
      else                { dend$labels[b] <- labels[child$id] }
    }

    assign(ka,-dend$count, clusters)
    assign(kb,-dend$count, clusters)
  }
}
