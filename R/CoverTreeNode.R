#' Cover Tree Node
#'
#' Constructs a cover tree node.  This constructor should proably never
#' be called directly.   Rather it should be called by CoverTree during
#' its construction and subsequent updates.
#' @param row index associated with the data
#' @param data A valid R object containing the data represented by the node
#' @param parent The parent node for the one being constructed (n/a for root node)
#' @param distance Either a numeric value or a distance function [ f(data,data) -> numeric ]
#' @return An initialized CoverTreeNode, fully connected as a child of 'parent'
#' @seealso \link{CoverTree}
#' @export
CoverTreeNode <- function(row,data,parent=NULL,distance=NULL)
{
  self <- environment()

  if(is.null(data))  { stop('The data argument must be provided') }

  level    <- NULL
  isRoot   <- is.null(parent)
  children <- new.env()
  isLeaf   <- TRUE

  if( isRoot )
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

    level  <- floor(1-log2(distance))
  }

  class(self) <- append( class(self), 'CoverTreeNode' )

  return(self)
}

#' @export
addChild <- function(node,row,data,distance)
{
  UseMethod('addChild',node)
}

#' @export
addChild.CoverTreeNode <- function(node,row,data,distance)
{
  if(distance == 0)
  {
    warning("Attempted to add row ",row,", as a child node of row ",node$row," at distance 0", call.=FALSE)
    return(node)
  }
  child <- CoverTreeNode(row,data,node,distance)

  if( is.null(node$level) ) #this must be root node (first child being added)
  {
    node$level <- child$level - 1
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

  level <- as.character(child$level)
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
  node.print(node,'  ')
  cat("\n")
}

node.print <- function(node,prefix)
{
  if(is.null(node$distance))
  {
    msg <- sprintf("%s[%2d]: %s (ROOT)\n",prefix,node$level,node$row);
  }
  else
  {
    msg <- sprintf("%s[%2d]: %3d (d=%s)\n",prefix,node$level,node$row,node$distance);
  }
  cat(noquote(msg));
  for( level in ls(node$children ) )
  {
    children <- get(level,node$children)
    for( child in children )
    {
      node.print(child,paste(prefix,"  "));
    }
  }
}

#' @export
node_to_dataframe <- function(node)
{
  p <- ifelse(node$isRoot, NA, node$parent$row)
  d <- ifelse(node$isRoot, NA, node$distance)
  rval <- data.frame( row=node$row, level=node$level, parent=p, distance=d)

  for(level in ls(node$children))
  {
    children = get(level,node$children)
    for(child in children)
    {
      st <- node_to_dataframe(child)
      rval <- rbind(rval,st)
    }
  }

  return(rval)
}
