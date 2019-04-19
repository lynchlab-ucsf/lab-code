#This is the code for pvclust.parallel, which is referenced in the pvclust function. I will need to source this file after loading the pvclust package, but before running the pvclust function.
#This code is located by typing: getAnywhere("pvclust.parallel") in the R environment, after loading the pvclust package.
# My edits entail removing the necessity for "method.dist", and instead allowing a distance matrix to be included as data instead. If method.dist is used, the algorithm will compute the dm with the distance method/function provided. However, if method.dist is not used (or NULL), the algorithm will consider the data as a distance matrix..
# Started on January 5th, 2017.
#  Katie McCauley

pvclust.parallel <- function (cl, data, method.hclust, method.dist, use.cor, nboot, 
    r, store, weight, init.rand = NULL, iseed, quiet, parallel.check) 
{
    if (parallel.check) {
        check.result <- check.parallel(cl = cl, nboot = nboot)
        if (!check.result) {
            msg <- paste(attr(check.result, "msg"), ". non-parallel version is executed", 
                sep = "")
            warning(msg)
            return(pvclust.nonparallel(data = data, method.hclust = method.hclust, 
                method.dist = method.dist, use.cor = use.cor, 
                nboot = nboot, r = r, store = store, weight = weight, 
                iseed = iseed, quiet = quiet))
        }
    }
    pkg.ver <- parallel::clusterCall(cl, packageVersion, pkg = "pvclust")
    r.ver <- parallel::clusterCall(cl, getRversion)
    if (length(unique(pkg.ver)) > 1 || length(unique(r.ver)) > 
        1) {
        node.name <- parallel::clusterEvalQ(cl, Sys.info()["nodename"])
        version.table <- data.frame(node = seq_len(length(node.name)), 
            name = unlist(node.name), R = unlist(lapply(r.ver, 
                as.character)), pvclust = unlist(lapply(pkg.ver, 
                as.character)))
        if (nrow(version.table) > 10) 
            table.out <- c(capture.output(print(head(version.table, 
                n = 10), row.names = FALSE)), "    ...")
        else table.out <- capture.output(print(version.table, 
            row.names = FALSE))
        warning("R/pvclust versions are not unique:\n", paste(table.out, 
            collapse = "\n"))
    }
    if (!is.null(init.rand)) 
        warning("\"init.rand\" option is deprecated. It is available for back compatibility but will be unavailable in the future.\nSpecify a non-NULL value of \"iseed\" to initialize random seed.")
    if (!is.null(iseed) && (is.null(init.rand) || init.rand)) 
        parallel::clusterSetRNGStream(cl = cl, iseed = iseed)
    n <- nrow(data)
    p <- ncol(data)
#If you have a function, this part allows the function to be used to calculate the distance matrix.
    if (is.function(method.dist)) {
        distance <- method.dist(data)
    }
#This part is what needs to be changed to allow for a specific distance matrix to be used (but I wanted to keep the old stuff, so the "method.dist" option could still be used.
    else {
         if(is.null(method.dist)) {
              distance <- as.dist(data)
         }
         else {
         distance <- dist.pvclust(data, method = method.dist, 
            use.cor = use.cor)
	 }
    }
    data.hclust <- hclust(distance, method = method.hclust)
    if (method.hclust == "ward" && getRversion() >= "3.1.0") {
        method.hclust <- "ward.D"
    }
    size <- floor(n * r)
    rl <- length(size)
    if (rl == 1) {
        if (r != 1) 
            warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
        r <- list(1)
    }
    else r <- as.list(size/n)
    ncl <- length(cl)
    nbl <- as.list(rep(nboot%/%ncl, times = ncl))
    if ((rem <- nboot%%ncl) > 0) 
        nbl[1:rem] <- lapply(nbl[1:rem], "+", 1)
    if (!quiet) 
        cat("Multiscale bootstrap... ")
    mlist <- parallel::parLapply(cl, nbl, pvclust.node, r = r, 
        data = data, object.hclust = data.hclust, method.dist = method.dist, 
        use.cor = use.cor, method.hclust = method.hclust, store = store, 
        weight = weight, quiet = quiet)
    if (!quiet) 
        cat("Done.\n")
    mboot <- mlist[[1]]
    for (i in 2:ncl) {
        for (j in 1:rl) {
            mboot[[j]]$edges.cnt <- mboot[[j]]$edges.cnt + mlist[[i]][[j]]$edges.cnt
            mboot[[j]]$nboot <- mboot[[j]]$nboot + mlist[[i]][[j]]$nboot
            mboot[[j]]$store <- c(mboot[[j]]$store, mlist[[i]][[j]]$store)
        }
    }
    result <- pvclust.merge(data = data, object.hclust = data.hclust, 
        mboot = mboot)
    return(result)
}

