f.get.subtrees <-function(tree, k = NULL, h = NULL,
           add.element = NULL){
    ## tree is a tree, k and h are as in cutree.
    ## additional elements is the name of any additional info
    ## connected with tree nodes. The additional elements referred
    ## to should have the same length as tree$labels.
    ##
    groups <- cutree(tree, k, h)
    subtrees <- list()
    for(i in 1:length(unique(groups))){
      subtrees[[i]] <- f.make.subtree(tree, groups, i,
                                      add.element = add.element)
    }
    subtrees
  }


f.make.subtree <-function (tree, groups, n = 1, add.element = NULL)
  {
    which.nodes <- which(groups == n)
    names(which.nodes) <- NULL
    ## steal some code from plot.hclust
    use.labels <-
      if (is.null(tree$labels))
        paste(1:(nrow(tree$merge) + 1))
    else as.character(tree$labels)
    
    
    subtree.labels <- use.labels[which.nodes]
    
    
    subtree.order <- tree$order[use.labels[tree$order] %in%
                                  subtree.labels]
    which.merge.elements <- which(-tree$merge %in% which.nodes)
    which.merge.rows <- sort(unique(which.merge.elements%%nrow(tree$merge)))
    which.merge.rows[which.merge.rows == 0] <- nrow(tree$merge)
    old.length <- 0
    new.length <- length(which.merge.rows)
    if (length(which.merge.elements) == 1) {
      res <- list(merge = NULL, heights = tree$height[which.merge.rows],
                  order = 1,
                  labels = {if(is.null(tree$labels)) NULL else
                    subtree.labels},
                  method = tree$method,
                  call = tree$call, dist.method = tree$dist.method)
    }
    else {
      while (new.length - old.length > 0) {
        old.length <- new.length
        new.rows <- numeric(0)
        more.new.rows <- numeric(0)
        row.elements <- as.vector(tree$merge[which.merge.rows,
                                             ])
        pos.elements <- row.elements[row.elements > 0]
        if (length(pos.elements) > 0)
          new.rows <- pos.elements[apply(tree$merge[pos.elements,
                                                    , drop = F], 1,
                                         function(x) {
                                           any(x > 0) & all(-x[x <
                                                                 0] %in% which.nodes)
                                         })]
        more.new.rows <- which(apply(tree$merge, 1, function(x) {
          all(x %in% which.merge.rows)
        }))
        which.merge.rows <- sort(unique(c(which.merge.rows,
                                          new.rows, more.new.rows)))
        new.length <- length(which.merge.rows)
      }
      merge.list <- tree$merge[which.merge.rows, , drop = F]
      merge.height <- tree$height[which.merge.rows]
      pos.elements <- sort(merge.list[merge.list > 0])
      for (i in seq(along = pos.elements)) merge.list[merge.list ==
                                                        pos.elements[i]] <- which(pos.elements[i] ==
                                                                                    which.merge.rows)
      neg.elements <- merge.list[merge.list < 0]
      minus.neg.elements <- -neg.elements
      num.elements <- length(neg.elements)
      change.table <- cbind(sort(minus.neg.elements), 1:num.elements)
      for (i in 1:num.elements) {
        merge.list[merge.list == -change.table[i, 1]] <- -change.table[i,
                                                                       2]
      }
      old.order <- subtree.order
      subtree.order <- match(use.labels[old.order], subtree.labels)
      res <- list(merge = merge.list, height = merge.height,
                  order = subtree.order,
                  labels = {if(is.null(tree$labels)) NULL else
                    subtree.labels},
                  method = tree$method,
                  call = list(match.call(), tree$call), dist.method =
                    tree$dist.method)
      class(res) <- "hclust"
      for (i in seq(along = add.element)) {
        res[[add.element[i]]] <- tree[[add.element[i]]][which.nodes]
      }
    }
    res
  }