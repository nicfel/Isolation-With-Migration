## read.nexus.R (2014-10-20)

##   Read Tree File in Nexus Format

## Copyright 2003-2014 Emmanuel Paradis and 2010 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.treeBuildWithTokens <- function(x) {
  phy <- .Call(treeBuildWithTokens, x)
  dim(phy[[1]]) <- c(length(phy[[1]])/2, 2)
  nms <- c("edge", "edge.length", "Nnode", "node.label", "root.edge")
  if (length(phy) == 4) nms <- nms[-5]
  names(phy) <- nms
  if (all(phy$node.label == "")) phy$node.label <- NULL
  class(phy) <- "phylo"
  attr(phy, "order") <- "cladewise"
  phy
}

my.read.nexus <- function(file, tree.names = NULL) {
  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)

  endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
  if(length(endblock)==0) endblock = length(X)
  semico <- grep(";", X)
  i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
  if(length(i1)==0) i1 =0
  i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
  translation <- if (length(i2) == 1 && i2 > i1) TRUE else FALSE
  if (translation) {
    end <- semico[semico > i2][1]
    x <- X[(i2 + 1):end] # assumes there's a 'new line' after "TRANSLATE"
    ## x <- gsub("TRANSLATE", "", x, ignore.case = TRUE)
    x <- unlist(strsplit(x, "[,; \t]"))
    x <- x[nzchar(x)]
    TRANS <- matrix(x, ncol = 2, byrow = TRUE)
    TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    n <- dim(TRANS)[1]
  }
  start <-
    if (translation) semico[semico > i2][1] + 1
  else i1 + 1 # semico[semico > i1][1] ## fix done on 2014-08-25
  end <- endblock[endblock > i1][1] - 1
  if(is.na(end)) end = length(X)
  tree <- X[start:end]
  rm(X)
  
  ## check whether there are empty lines from the above manips:
  tree <- tree[tree != ""]
  semico <- grep(";", tree)
  Ntree <- length(semico) # provisional -- some ";" may actually mark end of commands
  ## are some trees on several lines?
  ## -- this actually 'packs' all characters ending with a ";" in a single string --
  if (Ntree == 1 && length(tree) > 1) STRING <- paste(tree, collapse = "") else {
    if (any(diff(semico) != 1)) {
      STRING <- character(Ntree)
      s <- c(1, semico[-Ntree] + 1)
      j <- mapply(":", s, semico)
      if (is.list(j)) {
        for (i in 1:Ntree)
          STRING[i] <- paste(tree[j[[i]]], collapse = "")
      } else {
        for (i in 1:Ntree)
          STRING[i] <- paste(tree[j[, i]], collapse = "")
      }
    } else STRING <- tree
  }
  rm(tree)
  ## exclude the possible command lines ending with ";":
  #STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
  Ntree <- length(STRING) # update Ntree
  ## get the tree names:
  nms.trees <- sub(" *= *.*", "", STRING) # only the first occurence of "="
  nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", nms.trees, ignore.case = TRUE) # fix by Graham Gower (2014-10-20)
  STRING <- sub("^.* = ", "", STRING) # delete title and 'TREE' command with 'sub'
  STRING <- gsub(" ", "", STRING) # delete all white spaces
  colon <- grep(":", STRING)
  if (length(colon) == Ntree) {
    trees <-
      if (translation) lapply(STRING, .treeBuildWithTokens)
    else lapply(STRING, my.tree.build)
  } else {
    trees <- vector("list", Ntree)
    trees[colon] <- lapply(STRING[colon], tree.build)
    nocolon <- (1:Ntree)[!1:Ntree %in% colon]
    trees[nocolon] <- lapply(STRING[nocolon], clado.build)
    if (translation) {
      for (i in 1:Ntree) {
        tr <- trees[[i]]
        for (j in 1:n) {
          ind <- which(tr$tip.label[j] == TRANS[, 1])
          tr$tip.label[j] <- TRANS[ind, 2]
        }
        if (!is.null(tr$node.label)) {
          for (j in 1:length(tr$node.label)) {
            ind <- which(tr$node.label[j] == TRANS[, 1])
            tr$node.label[j] <- TRANS[ind, 2]
          }
        }
        trees[[i]] <- tr
      }
      translation <- FALSE
    }
  }
  for (i in 1:Ntree) {
    tr <- trees[[i]]
    ## Check here that the root edge is not incorrectly represented
    ## in the object of class "phylo" by simply checking that there
    ## is a bifurcation at the root
    if (!translation) n <- length(tr$tip.label)
    ROOT <- n + 1
  }
  if (Ntree == 1) {
    trees <- trees[[1]]
    if (translation) {
      trees$tip.label <-
        if (length(colon)) TRANS[, 2] else
          TRANS[, 2][as.numeric(trees$tip.label)]
    }
  } else {
    if (!is.null(tree.names)) names(trees) <- tree.names
    if (translation) {
      if (length(colon) == Ntree) # .treeBuildWithTokens() was used
        attr(trees, "TipLabel") <- TRANS[, 2]
      else { # reassign the tip labels then compress
        for (i in 1:Ntree)
          trees[[i]]$tip.label <-
          TRANS[, 2][as.numeric(trees[[i]]$tip.label)]
        trees <- .compressTipLabel(trees)
      }
    }
    class(trees) <- "multiPhylo"
    if (!all(nms.trees == "")) names(trees) <- nms.trees
  }
  trees
}

my.tree.build <- function(tp) {
  add.internal <- function() {
    edge[j, 1] <<- current.node
    edge[j, 2] <<- current.node <<- node <<- node + 1L
    index[node] <<- j # set index
    j <<- j + 1L
  }
  add.terminal <- function() {
    edge[j, 1] <<- current.node
    edge[j, 2] <<- tip
    index[tip] <<- j # set index
    X <- unlist(strsplit(tpc[k], ":"))
    label = unlist(strsplit(X[1],"[]\\=\\[]"))
    if(label[1]=="404") browser()
    tip.label[tip] <<- label[1]
    tip.meta[tip] <<- label[3]
    edge.length[j] <<- as.numeric(X[2])
    k <<- k + 1L
    tip <<- tip + 1L
    j <<- j + 1L
  }
  go.down <- function() {
    l <- index[current.node]
    X <- unlist(strsplit(tpc[k], ":"))
    label = unlist(strsplit(X[1],"[]\\=\\[]"))
    node.label[current.node - nb.tip] <<- label[1]
    node.meta[current.node - nb.tip] <<- label[3]
    edge.length[l] <<- as.numeric(X[2])
    k <<- k + 1L
    current.node <<- edge[l, 1]
  }
  if (!length(grep(",", tp))) {
    obj <- list(edge = matrix(c(2L, 1L), 1, 2))
    tp <- unlist(strsplit(tp, "[\\(\\):;]"))
    obj$edge.length <- as.numeric(tp[3])
    obj$Nnode <- 1L
    obj$tip.label <- tp[2]
    if (tp[4] != "") obj$node.label <- tp[4]
    class(obj) <- "phylo"
    return(obj)
  }
  
  tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
  tpc <- tpc[nzchar(tpc)]
  ## the following 2 lines are (slightly) faster than using gsub()
  tsp <- unlist(strsplit(tp, NULL))
  skeleton <- tsp[tsp %in% c("(", ")", ",", ";")]
  nsk <- length(skeleton)
  nb.node <- sum(skeleton == ")")
  nb.tip <- sum(skeleton == ",") + 1
  ## We will assume there is an edge at the root;
  ## if so, it will be removed and put into a vector
  nb.edge <- nb.node + nb.tip
  node.label <- character(nb.node)
  tip.label <- character(nb.tip)
  node.meta <- character(nb.node)
  tip.meta <- character(nb.tip)
  
  edge.length <- numeric(nb.edge)
  edge <- matrix(0L, nb.edge, 2)
  current.node <- node <- as.integer(nb.tip + 1) # node number
  edge[nb.edge, 2] <- node
  index <- numeric(nb.edge + 1) # hash index to avoid which
  index[node] <- nb.edge
  
  ## j: index of the line number of edge
  ## k: index of the line number of tpc
  ## tip: tip number
  j <- k <- tip <- 1L
  
  for (i in 2:nsk) {
    if (skeleton[i] == "(") add.internal() # add an internal branch (on top)
    if (skeleton[i] == ",") {
      if (skeleton[i - 1] != ")") add.terminal() # add a terminal branch
    }
    if (skeleton[i] == ")") {
      if (skeleton[i - 1] == "," || skeleton[i - 1] == "(") { # add a terminal branch and go down one level
        add.terminal()
        go.down()
      }
      if (skeleton[i - 1] == ")") go.down() # go down one level
    }
  }
  
  edge <- edge[-nb.edge, ]
  obj <- list(edge = edge, Nnode = nb.node, tip.label = tip.label)
  root.edge <- edge.length[nb.edge]
  edge.length <- edge.length[-nb.edge]
  if (!all(is.na(edge.length))) # added 2005-08-18
    obj$edge.length <- edge.length
  if (is.na(node.label[1])) node.label[1] <- ""
  if (any(nzchar(node.label))) obj$node.label <- node.label
  if (is.na(node.meta[1])) node.meta[1] <- ""
  if (any(nzchar(node.meta))) obj$node.meta <- node.meta
  if (is.na(tip.meta[1])) tip.meta[1] <- ""
  if (any(nzchar(tip.meta))) obj$tip.meta <- tip.meta
  if (!is.na(root.edge)) obj$root.edge <- root.edge
  class(obj) <- "phylo"
  attr(obj, "order") <- "cladewise"
  obj
}
