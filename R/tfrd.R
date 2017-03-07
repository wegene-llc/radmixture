#' Transfer ped file to genotype matrix
#' @param referenceped ped file for your reference panel
#' @return genotype
#' @export

generateG <- function(referenceped) {
    ped <- referenceped[, -(1:6)]
    # calculate allele counts
    a1 <- ped[, seq(1, ncol(ped) - 1, 2)]
    a2 <- ped[, seq(2, ncol(ped), 2)]
    a1 <- as.matrix(a1)
    a2 <- as.matrix(a2)
    a <- rbind(a1, a2)
    # find major allele
    major <- rep(NA, ncol(a))
    minor <- rep(NA, ncol(a))
    del <- numeric()
    for(i in 1:ncol(a)) {
        freqco <- count(a[, i])
        if(nrow(freqco) != 2 || freqco[1, 2] == freqco[2, 2] || 0 %in% freqco[, 1]) {
            del[i] <- i
            major[i] <- NA
            minor[i] <- NA
        } else {
            major[i] <- freqco[which.max(freqco[, 2]), 1]
            minor[i] <- freqco[which.min(freqco[, 2]), 1]
        }
    }
    del <- which(!is.na(del))
    # generate two genotypes
    a1 <- a1[, -del]
    a2 <- a2[, -del]
    a <- a[, -del]
    major <- major[-del]
    minor <- minor[-del]
    genotype <- matrix(NA, 2, ncol(a))
    genotype[1, ] <- 10 * major + major
    genotype[2, ] <- 10 * minor + minor
    # generate G matrix
    genotype1 <- 10 * a1 + a2
    g <- matrix(NA, nrow(ped), ncol(a))
    for (i in 1:ncol(a)) {
        g[which(genotype1[, i] == genotype[1, i]), i] <- 0
        g[which(genotype1[, i] == genotype[2, i]), i] <- 2
    }
    g[is.na(g)] <- 1
    return(g = g)
}

#' Initialize Q and F
#' @param g genotype matrix
#' @param pop A data frame. If you intend to do supervised learning, 
#' you must specify the ancestries of the reference individuals.
#' @param alpha Parameter for dirichlet distribution.
#' Vector of shape parameters, or matrix of shape
#' parameters corresponding to the number of draw.
#' @param K If you intend to do unsupervised learning, 
#' set the number of populations you will use.
#' @param model Choose supervised or unsupervised learning
#' @export

initQF <- function(g, pop = NULL, alpha = NULL, K = NULL, model = c("supervised", "unsupervised")) {
    if (is.data.frame(g)) {
        g <- as.matrix(g)
    }
    if (model != "supervised" && model != "unsupervised") {
       stop("You must choose one model")
    }
    if (model == "supervised" && is.null(pop)) {
        stop("pop is needed under supervised learning")
    }
    if (model == "unsupervised" && is.null(K)) {
        stop("You must set up K for unsupervised learning")
    }
    if (model == "unsupervised" && is.null(alpha)) {
        stop("You must set up alpha for unsupervised learning")
    }
    if (model == "supervised") {
        if (nrow(pop) != nrow(g)) {
            stop("Check the number of individuals!")
        }
        pop1 <- as.character(pop[, 1]) %>%
            unique()
        num <- length(pop1) - 1
        q <- matrix(NA, nrow(pop), num)
        for(i in 1:num) {
            q[which(pop[, 1] == pop1[i]), i] <- 1 - 1e-5 * (num - 1)
        }
        q[nrow(q), ] <- rep(1 / num, num)
        q[is.na(q)] <- 1e-5
        f <- matrix(NA, num, ncol(g))
        for(i in 1:num) {
            f[i, ] <- colSums(g[which(pop[, 1] == pop1[i]), ]) / (2 * length(which(pop[, 1] == pop1[i])))
        }
    } else if (model == "unsupervised") {
        q <- rdirichlet(nrow(g), rep(alpha, K))
        f <- (t(q) %*% g) / (2 * nrow(g))
    }
    f[f == 0] <- lb
    f[f == 1] <- ub
    return(list(q = q, f = f))
}

#' transfer raw genotype file for f fixed admixture.
#' @param genotype A data frame contains your genotype information.
#' @param K the number of populations
#' @param map Index file, it should contain rsid, choromosome position,
#' major allele and minor allele information for both plus and minus strands.
#' @param f Frequency matrix
#' @return A list contains g, q, f which can be used for calculation.
#' @export

tfrdpub <- function(genotype, K, map, f) {
    if (K != ncol(f)) {
        stop("The number of populations does not match F matrix you use!")
    }
    overlap <- intersect(genotype[, 1], map[, 1])
    gindex <- match(overlap, genotype[, 1])
    mapindex <- match(overlap, map[, 1])
    genotype <- genotype[gindex, ]
    map <- map[mapindex, ]
    f <- f[mapindex, ]
    nocallindel <- which(genotype[, 4] == "--" |
                             genotype[, 4] == "__" | genotype[, 4] == "II" | genotype[, 4] == "DD")
    if (length(nocallindel) == 0) {
        g <- genotype
    } else {
        f <- f[-nocallindel, ]
        g <- genotype[-nocallindel, ]
        map <- map[-nocallindel, ]
    }
    f <- 1 - t(f)
    g1 <- rep(NA, nrow(g))
    g <- as.character(g[, 4])
    gt1 <- paste(map[, 2], map[, 2], sep = "")
    gt2 <- paste(map[, 4], map[, 4], sep = "")
    gt3 <- paste(map[, 3], map[, 3], sep = "")
    gt4 <- paste(map[, 5], map[, 5], sep = "")
    gt5 <- paste(map[, 2], map[, 3], sep = "")
    gt6 <- paste(map[, 3], map[, 2], sep = "")
    gt7 <- paste(map[, 4], map[, 5], sep = "")
    gt8 <- paste(map[, 5], map[, 4], sep = "")
    for (i in 1:length(g)) {
        if (g[i] == gt1[i] || g[i] == gt2[i]) {
            g1[i] <- 2
        } else if (g[i] == gt5[i] || g[i] == gt6[i] ||
                   g[i] == gt7[i] || g[i] == gt8[i]) {
            g1[i] <- 1
        }
        else if (g[i] == gt3[i] || g[i] == gt4[i]) {
            g1[i] <- 0
        }
    }
    xx <- which(is.na(g1))
    if (length(xx) == 0) {
        g1 <- g1
        f <- f
    } else {
        g1 <- g1[-xx]
        f <- f[, -xx]
    }
    g <- t(as.matrix(g1))
    q <- t(as.matrix(rep(1 / nrow(f), nrow(f))))
    if (nrow(map) < 40000) {
        warning("The number of SNPs is small, you may get unreasonable result!")
    }
    return(list(q = q, f = f, g = g))
}
