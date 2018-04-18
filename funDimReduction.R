calculatePcaReduction <- function(
    data.use,
    nPcs = 20,
    reduction.name = "PCA",
    reduction.key = "PC",
    use.subst = FALSE,  # to use only not 0 substitutions for all species
    subs = NULL,  # or select substitutions directly
    center = TRUE,
    species = NULL,
    rev.pca = FALSE,  # cells X genes matrix or species X substitution
    fastpath = TRUE,
    maxit = 100,
    weight.by.var = TRUE,
    seed.use = 42,
    ...
) {
    if (!is.null(seed.use)) {
        set.seed(seed = seed.use)
    }
    if (use.subst && is.null(subs)) {
        wthZeros <- purrr::map_lgl(data.use, any(. == 0))
        subs <- !wthZeros
    }
    if (!is.null(species)) {
        # species subset is just for PC determination
        data.use <- 
            data.use %>% 
            filter(Species %in% species) %>% 
            as.data.frame() %>% 
            column_to_rownames("Species")
    } else {
        data.use <- column_to_rownames(as.data.frame(data.use), "Species")
    }
    if (!is.null(subs)) {
        data.use <- select(data.use, vars(subs))
    }
    data.use <- t(data.use)
    library(irlba)
    if (rev.pca) {
        nPcs <- min(nPcs, ncol(x = data.use) - 1)
        tdata.use <- t(data.use)
        if (center) {
            cm <- Matrix::colMeans(tdata.use)
            pcs <-
                irlba(
                    A          = tdata.use,
                    nv         = nPcs,
                    center     = cm,
                    right_only = FALSE,
                    fastpath   = fastpath,
                    maxit      = maxit,
                    reorth     = TRUE,
                    tol        = 1e-8,
                    ...
                )
        } else {
            pcs <-
                irlba(
                    A          = tdata.use,
                    nv         = nPcs,
                    right_only = FALSE,
                    fastpath   = fastpath,
                    maxit      = maxit,
                    reorth     = TRUE,
                    tol        = 1e-8,
                    ...
                )
        }
        sdev <<- pcs$d / sqrt(max(1, nrow(data.use) - 1))
        if(weight.by.var){
            subs.loadings <- pca.results$u %*% diag(pca.results$d)
        } else{
            subs.loadings <- pca.results$u
        }
        # adjust for centering!
        if (center) {
            pcs$center <- cm
            species.embeddings <- as.matrix(
                t(t(tdata.use %*% pcs$v) - as.vector(t(cm %*% pcs$v))))
        } else {
            species.embeddings <- as.matrix(tdata.use %*% pcs$v)
        }
    }
    else {
        nPcs <- min(nPcs, nrow(x = data.use) - 1)
        if (center) {
            cm <- Matrix::colMeans(data.use)
            pcs <-
                irlba(
                    A          = data.use,
                    nv         = nPcs,
                    center     = cm,
                    right_only = FALSE,
                    fastpath   = fastpath,
                    maxit      = maxit,
                    reorth     = TRUE,
                    tol        = 1e-8,
                    ...
                )
        } else {
            pcs <-
                irlba(
                    A          = data.use,
                    nv         = nPcs,
                    right_only = FALSE,
                    fastpath   = fastpath,
                    maxit      = maxit,
                    reorth     = TRUE,
                    tol        = 1e-8,
                    ...
                )
        }
        # adjust for centering!
        if (center) {
            pcs$center <- cm
            subs.loadings <- as.matrix(
                t(t(data.use %*% pcs$v) - as.vector(t(cm %*% pcs$v))))
        } else {
            subs.loadings <- as.matrix(data.use %*% pcs$v)
        }
        sdev <<- pcs$d / sqrt(max(1, ncol(data.use) - 1))
        if (weight.by.var) {
            species.embeddings <- pcs$v %*% diag(pcs$d)
        } else {
            species.embeddings <- pcs$v
        }
    }
    
    
    rownames(x = subs.loadings) <- rownames(x = data.use)
    colnames(x = subs.loadings) <- paste0(reduction.key, 1:nPcs)
    rownames(x = species.embeddings) <- colnames(x = data.use)
    colnames(x = species.embeddings) <- colnames(x = subs.loadings)
    subs.loadings <<- subs.loadings <- as.data.frame(subs.loadings) %>% rownames_to_column("Subs")
    species.embeddings <<- species.embeddings <- as.data.frame(species.embeddings) %>% rownames_to_column("Species")
    
    if (reduction.name == 'ICA' && rev.pca == FALSE) {
        pcas <- subs.loadings %>% select(-Species) %>% as.matrix()
        nIcs <- nPcs
        library(fastICA)
        a <-
            fastICA::ica.R.def(
                t(pcas),
                nIcs,
                tol = 1e-3,
                fun = 'logcosh',
                maxit = 200,
                verbose = TRUE,
                alpha = 1,
                w.init = matrix(rnorm(nIcs * nPcs), nIcs, nPcs)
            )
        ICA <- as.matrix(data.use %*% pcs$v %*% a)
        colnames(ICA) <- paste('IC', seq(ncol(ICA)), sep = '')
        return(ICA)
    } else {
        data.use <- as.data.frame(t(data.use))
        data.use <- 
            rownames_to_column(data.use, "Species") %>% 
            inner_join(species.embeddings, by = "Species")
        return(data.use)
    }
}

