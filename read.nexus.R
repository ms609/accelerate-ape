read.nexus.fast <- function (file, tree.names = NULL) {
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT)) {
        w <- LEFT == RIGHT
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        }
        w <- !w
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        }
    }
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- if (length(i2) == 1 && i2 > i1) TRUE else FALSE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- X[(i2 + 1):end]
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        n <- dim(TRANS)[1]
    }
    start <- if (translation) {
        semico[semico > i2][1] + 1
    } else {
      i1 + 1
    }
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]
    rm(X)
    tree <- tree[tree != ""]
    semico <- grep(";", tree)
    Ntree <- length(semico)
    if (Ntree == 1 && length(tree) > 1) {
        STRING <- paste(tree, collapse = "")
    } else {
        if (any(diff(semico) != 1)) {
            STRING <- character(Ntree)
            s <- c(1, semico[-Ntree] + 1)
            j <- mapply(":", s, semico)
            if (is.list(j)) {
                for (i in 1:Ntree) STRING[i] <- paste(tree[j[[i]]], 
                  collapse = "")
            }
            else {
                for (i in 1:Ntree) STRING[i] <- paste(tree[j[, 
                  i]], collapse = "")
            }
        }
        else STRING <- tree
    }
    rm(tree)
    STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
    Ntree <- length(STRING)
    nms.trees <- sub(" *= *.*", "", STRING)
    nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", 
        nms.trees, ignore.case = TRUE)
    STRING <- sub("^.*= *", "", STRING)
    STRING <- gsub(" ", "", STRING)
    colon <- grep(":", STRING)
    if (!length(colon)) {
        trees <- lapply(STRING, clado.build)
    }
    else if (length(colon) == Ntree) {
        trees <- if (translation) 
            lapply(STRING, .treeBuildWithTokens)
        else lapply(STRING, tree.build)
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
                    ind <- which(tr$node.label[j] == TRANS[, 
                      1])
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
        if (!translation) 
            n <- length(tr$tip.label)
        ROOT <- n + 1
        if (sum(tr$edge[, 1] == ROOT) == 1 && dim(tr$edge)[1] > 
            1) {
            stop(paste("The tree has apparently singleton node(s): cannot read tree file.\n  Reading NEXUS file aborted at tree no.", 
                i, sep = ""))
        }
    }
    if (Ntree == 1) {
        trees <- trees[[1]]
        if (translation) {
            trees$tip.label <- if (length(colon)) 
                TRANS[, 2]
            else TRANS[, 2][as.numeric(trees$tip.label)]
        }
    }
    else {
        if (!is.null(tree.names)) 
            names(trees) <- tree.names
        if (translation) {
            if (length(colon) == Ntree) 
                attr(trees, "TipLabel") <- TRANS[, 2]
            else {
                for (i in 1:Ntree) trees[[i]]$tip.label <- TRANS[, 
                  2][as.numeric(trees[[i]]$tip.label)]
                trees <- .compressTipLabel(trees)
            }
        }
        class(trees) <- "multiPhylo"
        if (!all(nms.trees == "")) 
            names(trees) <- nms.trees
    }
    trees
}

clado.build <- function (tp) {
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2L, 1L), 1L, 2L), Nnode = 1L)
        tp <- unlist(strsplit(tp, "[\\(\\);]"))
        obj$tip.label <- tp[2]
        if (tp[3] != "") 
            obj$node.label <- tp[3]
        class(obj) <- "phylo"
        return(obj)
    }
    skeleton <- unlist(strsplit(gsub('([^\\(\\),;])', '', tp), NULL))
    tp <- gsub(")", ")NA", tp)
    tp <- gsub(" ", "", tp)
    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
    tpc <- tpc[tpc != ""]
    nsk <- length(skeleton)
    nb.node <- length(skeleton[skeleton == ")"])
    nb.tip <- length(skeleton[skeleton == ","]) + 1L
    nb.edge <- nb.node + nb.tip
    node.label <- character(nb.node)
    tip.label <- character(nb.tip)
    ### edge <- matrix(NA_integer_, nb.edge, 2L) ##DEL
    parent <- rep(NA_integer_, nb.edge)
    child <- parent
    node <- nb.tip + 1L
    current.node <- node
    parent[nb.edge] <- 0L
    child [nb.edge] <- node
    index <- numeric(nb.edge + 1L)
    index[node] <- nb.edge
    j <- k <- tip <- 1L
    for (i in 2:nsk) {
        switch(skeleton[i],
        "(" = { # Add internal
            parent[j] <- current.node
            node <- node + 1L
            current.node <- node
            child[j] <- current.node
            index[node] <- j
            j <- j + 1L
        },
        "," = {
            if (skeleton[i - 1] != ")") {  
                # Add terminal
                parent[j] <- current.node
                child [j] <- tip
                index[tip] <- j
                tip.label[tip] <- tpc[k]
                k <- k + 1L
                tip <- tip + 1L
                j <- j + 1L
            }
        },
        ")" = switch(skeleton[i - 1],
            "," = {
                # Add terminal
                parent[j] <- current.node
                child [j] <- tip
                index[tip] <- j
                tip.label[tip] <- tpc[k]
                k <- k + 1L
                tip <- tip + 1L
                j <- j + 1L
                # Go down
                l <- index[current.node]
                node.label[current.node - nb.tip] <- tpc[k]
                k <- k + 1L
                current.node <- parent[l]
            },
            ")" = {
                # Go down
                l <- index[current.node]
                node.label[current.node - nb.tip] <- tpc[k]
                k <- k + 1L
                current.node <- parent[l]
            })
        )
    }
    edge <- cbind(parent, child)[-nb.edge, ]
    obj <- list(edge = edge, tip.label = tip.label, Nnode = nb.node, node.label = node.label)
    obj$node.label <- if (all(obj$node.label == "NA")) NULL else gsub("^NA", "", obj$node.label)
    class(obj) <- "phylo"
    attr(obj, "order") <- "cladewise"
    obj
}


























