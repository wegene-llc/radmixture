#' Transfer raw genotype file for admixture
#' @param genotype A data frame which contains your genotype information. It should have 4 columns,
#' the first column is rsid, second is chromosome, third is position and the forth is genotype.
#' @param map Index file, it should contain rsid, choromosome and position.
#' @param referenceped ped file for your reference data.
#' @param f Initial frequency matrix matches your ped file.
#' @return A list which contains genotype matrix and frequency matrix.
#' @export
tfrd <- function(genotype, map, referenceped, f) {
        cha <- paste(genotype[, 2], genotype[, 3], sep = ",")
        cha1 <- paste(map[, 1], map[, 2], sep = ",")
        overlap <- intersect(cha, cha1)
        index_user <- match(overlap, cha)
        index_wg <- match(overlap, cha1)
        ped <- referenceped[, - (1:6)] %>%
            as.matrix()
        ped <- ped[, indp(index_wg)]
        f <- f[, index_wg]
        genotype <- genotype[index_user, 4]
        genotype <- as.character(genotype)
        newped <- paste(genotype, collapse = "") %>%
            strsplit(split = "") %>%
            unlist()
        newped[newped == "A"] <- 1
        newped[newped == "C"] <- 2
        newped[newped == "G"] <- 3
        newped[newped == "T"] <- 4
        newped[newped == "-"] <- 0
        newped[newped == "_"] <- 0
        newped[newped == "I"] <- 0
        newped <- as.numeric(newped)
        del <- which(is.na(newped))
        newped[del] <- 0
        ped <- rbind(ped, newped)
        gf <- gg(ped, f)
        if (ncol(gf$g) < 40000) {
            warning("The number of SNPs is
                    inappropriate for ancestry estimation!")
        } else if (ncol(gf$g) < 50000) {
            warning("The number of SNPs used for calculating is small,
                    so your result might be unreasonable!")
        }
        return(list(g = gf$g, f = gf$f))
}

#' transfer raw genotype file for f fixed admixture. Mostly for public dataset.
#' @param genotype A data frame contains your genotype information.
#' @param K the number of populations
#' @param map Index file, it should contain rsid, choromosome position,
#' major allele and minor allele information.
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
                             genotype[, 4] == "__" | genotype[, 4] == "II")
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
