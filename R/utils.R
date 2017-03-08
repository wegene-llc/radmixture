cacheEnv <- new.env()
assign("lb", 1e-5, envir = cacheEnv)
assign("ub", 0.99999, envir = cacheEnv)

lb <- get("lb", envir = cacheEnv)
ub <- get("ub", envir = cacheEnv)
# The function lfun() is used to calculate log-likelihood value.
lfun <- function(g, q, f) {
    l1 <- g * log(q %*% f)
    l2 <- (2 - g) * log(q %*% (1 - f))
    l <- l1 + l2
    return(sum(l))
}
# update q
em.q <- function(g, q0, q, f0, p, s) {
    for (k in 1:ncol(q0)) {
        c1 <- g[s, ] * q0[s, k] * f0[k, ]
        c2 <- q0[s, ] %*% f0
        c3 <- (2 - g[s, ]) * q0[s, k] * (1 - f0[k, ])
        c4 <- q0[s, ] %*% (1 - f0)
        c <- c1 / c2 + c3 / c4
        q[s, k] <- sum(c) / (2 * p)
    }
    q[q < lb] <- lb
    q[q > ub] <- ub
    q[s, ] <- q[s, ] / sum(q[s, ])
    return(q = q)
}

# Calculate Jacobian and Hessian for F when Q fixed

diff.f <- function(g, q, f, s) {
    qf <- q %*% f[, s]
    qf2 <- q %*% (1 - f[, s])
    gq <- crossprod(g[, s] / qf, q)
    gq2 <- crossprod( (2 - g[, s]) / qf2, q)
    df <- t(gq - gq2)
    qqf <- qf * qf
    qqf2 <- qf2 * qf2
    q_temp <- sweep(t(q), MARGIN = 2, g[, s] / qqf, "*")
    q_temp2 <- sweep(t(q), MARGIN = 2, (2 - g[, s]) / qqf2, "*")
    ddf <- -1 * (q_temp %*% q + q_temp2 %*% q)
    res <- list(Hess = ddf, Jaco = t(df))
    return(res)
}

# Calculate Jacobian and Hessian for Q when F fixed

diff.q <- function(g, q, f, s) {
    qf <- q[s, ] %*% f
    qf2 <- q[s, ] %*% (1 - f)
    gf <- tcrossprod(g[s, ] / qf, f)
    gf2 <- tcrossprod( (2 - g[s, ]) / qf2, 1 - f)
    dq <- gf + gf2
    qqf <- qf * qf
    qqf2 <- qf2 * qf2
    f_temp <- sweep(f, MARGIN = 2, g[s, ] / qqf, "*")
    f_temp2 <- sweep( (1 - f), MARGIN = 2, (2 - g[s, ]) / qqf2, "*")
    ddq <- -1 * (tcrossprod(f_temp, f) + tcrossprod(f_temp2, (1 - f)))
    res <- list(Hess = ddq, Jaco = dq)
    return(res)
}

# Quadratic Programming (QP) for F
brf <- function(g, q, f) {
    lb <- 1e-5
    ub <- 1 - 1e-5
    Amatf <- rbind(diag(ncol(q)), -1 * diag(ncol(q)))
    flongvec <- matrix(NA, nrow(f), ncol(f))
    for (j in 1:ncol(f)) {
        bvec <- c(-f[, j], f[, j] - 1)
        flongvec[, j] <- f[, j] + solve.QP(-diff.f(g, q, f, s = j)$Hess * 1e-2,
                                           diff.f(g, q, f, s = j)$Jaco * 1e-2,
                                           t(Amatf), bvec)$solution
    }
    flongvec[flongvec < lb] <- lb
    flongvec[flongvec > ub] <- ub
    return(f = flongvec)
}

# Quadratic Programming (QP) for Q
brq <- function(g, q, f, model) {
    lb <- 1e-5
    ub <- 1 - 1e-5
    Amatq <- rbind(rep(1, ncol(q)), diag(ncol(q)))
    if (model == "supervised") {
        bvec <- c(0, -q[nrow(q), ])
        q[nrow(q), ] <- q[nrow(q), ] + solve.QP(-diff.q(g, q, f, nrow(q))$Hess * 1e-2,
                                                diff.q(g, q, f, nrow(q))$Jaco * 1e-2,
                                                t(Amatq), bvec,
                                                meq = 1)$solution
        q[q < lb] <- lb
        q[q > ub] <- ub
        q[nrow(q), ] <- q[nrow(q), ] / sum(q[nrow(q), ])
    } else {
        qlongvec <- matrix(NA, nrow(q), ncol(q))
        for (i in 1:nrow(q)) {
            bvec <- c(0, -q[i, ])
            qlongvec[i, ] <- q[i, ] + solve.QP(-diff.q(g, q, f, i)$Hess * 1e-2,
                                               diff.q(g, q, f, i)$Jaco * 1e-2,
                                               t(Amatq), bvec, meq = 1)$solution
            qlongvec[qlongvec < lb] <- lb
            qlongvec[qlongvec > ub] <- ub
            qlongvec[i, ] <- qlongvec[i, ] / sum(qlongvec[i, ])
        }
        q <- qlongvec
    }
    return(q = q)
}

# quasi-newton (update Q)
qnq <- function(q, qvector, model) {
    if (model == "supervised") {
        u <- qvector[, 2] - qvector[, 1]
        v <- qvector[, 3] - qvector[, 2]
        if (u == 0 && v == 0) {
            c <- 1
        } else {
            c <- - (t(u) %*% u / t(u) %*% (v - u))
        }
        if (c == Inf || c == -Inf || is.nan(c)) {
            c <- 1
        }
        updateq <- (1 - c) * qvector[, 2] + c * qvector[, 3]
        updateq[updateq < lb] <- lb
        updateq[updateq > ub] <- ub
        updateq <- updateq / sum(updateq)
    } else {
        q1 <- matrix(qvector[, 1], nrow(q), ncol(q))
        q2 <- matrix(qvector[, 2], nrow(q), ncol(q))
        q3 <- matrix(qvector[, 3], nrow(q), ncol(q))
        updateq <- matrix(NA, nrow(q), ncol(q))
        for (i in 1:nrow(q)) {
            u <- q2[i, ] - q1[i, ]
            v <- q3[i, ] - q2[i, ]
            if (u == 0 && v == 0) {
                c <- 1
            } else {
                c <- - (t(u) %*% u / t(u) %*% (v - u))
            }
            if (c == Inf || c == -Inf || is.nan(c)) {
                c <- 1
            }
            updateq[i, ] <- (1 - c) * q2[i, ] + c * q3[i, ]
        }
        updateq[updateq < lb] <- lb
        updateq[updateq > ub] <- ub
        updateq[i, ] <- updateq[i, ] / sum(updateq[i, ])
    }
    return(updateq)
}
# quasi-newton (update F)
qnf <- function(f, fvector) {
    f1 <- matrix(fvector[, 1], nrow(f), ncol(f))
    f2 <- matrix(fvector[, 2], nrow(f), ncol(f))
    f3 <- matrix(fvector[, 3], nrow(f), ncol(f))
    updatef <- matrix(NA, nrow(f), ncol(f))
    for (i in 1:ncol(f)) {
        u <- f2[, i] - f1[, i]
        v <- f3[, i] - f2[, i]
        if (u == 0 && v == 0) {
            c <- 1
        } else {
            c <- - (t(u) %*% u / t(u) %*% (v - u))
        }
        if (c == Inf || c == -Inf || is.nan(c)) {
            c <- 1
        }
        updatef[, i] <- (1 - c) * f2[, i] + c * f3[, i]
    }
    updatef[updatef < lb] <- lb
    updatef[updatef > ub] <- ub
    return(updatef)
}

# index for ped

indp <- function(index) {
    res <- numeric(length = length(index) * 2)
    for(i in 1:length(index)) {
        res[2 * i - 1] <- index[i] * 2 - 1
        res[2 * i] <- index[i] * 2
    }
    res
}


gg <- function(ped, f) {
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
    if (length(del) == 0) {
        a1 <- a1
        a2 <- a2
        a <- a
        major <- major
        minor <- minor
    } else {
        a1 <- a1[, -del]
        a2 <- a2[, -del]
        a <- a[, -del]
        major <- major[-del]
        minor <- minor[-del]
    }
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
    f <- f[, -del]
    return(list(g = g, f = f, del = del))
}

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
    newped[newped == "D"] <- 0
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