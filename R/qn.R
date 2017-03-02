#' quasi-Newton acceleration
#' @description Use quasi-Newton algorithm to accelerate EM or block relaxation.
#' @param g Genotype matrix with dimensions \eqn{n × p}, where n is sample size
#' and p is the number of SNPs.
#' @param q Ancestry coefficient matrix with dimensions \eqn{n × K}, where n
#' is sample size and K is the number of populations.
#' @param f Minor allele frequency matrix with dimensions \eqn{K × p},
#' where K is the number of populations and p is the number of SNPs.
#' @param tol Tolerance, the default value is 1e-4.
#' @param method Choose which algorithm you want to use. EM or BR.
#' @param model Choose which model you want to use. Supervised learning or unsupervised learning.
#' @return Estimation results of q, f and the loglikelihood value of each iteration.
#' @export

qn <- function(g, q, f, tol = 1e-4, method = c("EM", "BR"),
               model = c("supervised", "unsupervised")) {
    em_res <- em(g, q, f, acc = T, max.iter = 5, model = model)
    if (method == "EM") {
        qvector <- em_res$qvector
        fvector <- em_res$fvector
    } else {
        br_res <- br(g, em_res$q, em_res$f, acc = T,
                     max.iter = 3, model = model)
        qvector <- br_res$qvector
        fvector <- br_res$fvector
    }
    qvector <- qvector[, (ncol(qvector) - 2):(ncol(qvector))]
    fvector <- fvector[, (ncol(fvector) - 2):(ncol(fvector))]
    iter <- 1
    loglike <- numeric()
    if (model == "supervised") {
        loglike[1] <- lfun(g, rbind(q[-nrow(q), ], qvector[, 1]),
                           matrix(fvector[, 1], nrow(f), ncol(f)))
    } else {
        loglike[1] <- lfun(g, matrix(qvector[, 1], nrow(q), ncol(q)),
                           matrix(fvector[, 1], nrow(f), ncol(f)))
    }
    q <- em_res$q
    f <- em_res$f
    repeat{
        iter <- iter + 1
        updateq <- qnq(q, qvector, model = model)
        updatef <- qnf(f, fvector)
        if (model == "supervised") {
            loglike[iter] <- lfun(g, rbind(q[-nrow(q), ], updateq), updatef)
        } else {
            loglike[iter] <- lfun(g, updateq, updatef)
        }
        if (loglike[iter] < loglike[iter - 1]) {
            if (model == "supervised") {
                q <- rbind(q[-nrow(q), ], qvector[, 3])
            } else {
                q <- matrix(qvector[, 3], nrow(q), ncol(q))
            }
            f <- matrix(fvector[, 3], nrow(f), ncol(f))
            loglike[iter] <- lfun(g, q, f)
        } else {
            if (model == "supervised") {
                q <- rbind(q[-nrow(q), ], updateq)
            } else {
                q <- updateq
            }
            f <- updatef
        }
        if (abs(lfun(g, q, f) - loglike[iter - 1]) < tol) {
            break
        } else {
            if (method == "EM") {
                fitmodel <- em(g, q, f, acc = T, max.iter = 3, model = model)
            } else {
                fitmodel <- br(g, q, f, acc = T, max.iter = 3, model = model)
            }
            qvector <- fitmodel$qvector
            fvector <- fitmodel$fvector
        }
    }
    if (model == "supervised") {
        qnew <- round(q[nrow(q), ], 5)
        res <- list(q = qnew, f = f, loglike = loglike[1:iter])
    } else {
        res <- list(q = q, f = f, loglike = loglike[1:iter])
    }
    return(res)
}


#' qn for fixed F
#' @param gnew Integer which length is the number of SNPs used in calculation.
#' @param qnew Initial q used in calculation. Numeric. Sum(q) must be 1.
#' @param f Allele frequencies learned from the reference panels.
#' @param tol Tolerance, the default value is 1e-4.
#' @param method Choose which algorithm you want to use. EM or BR.
#' @param pubdata You can choose a public dataset here, E11 or K13. You also can use other public
#' dataset which is not in this package.
#' @export

fFixQN <- function(gnew, qnew, f, tol = 1e-4,
                   method = c("EM", "BR"), pubdata = NULL) {
    if (method == "EM") {
        em_res <- fFixEm(gnew, qnew, f, acc = T, max.iter = 3)
        qvector <- em_res$qvector
    } else {
        br_res <- fFixBr(gnew, qnew, f, acc = T, max.iter = 3)
        qvector <- br_res$qvector
    }
    qvector <- qvector[, (ncol(qvector) - 2):(ncol(qvector))]
    iter <- 1
    loglike <- numeric()
    loglike[1] <- lfun(gnew, rbind(qnew[-nrow(qnew), ], qvector[, 1]), f)
    repeat{
        iter <- iter + 1
        updateq <- qnq(qnew, qvector, model = "supervised")
        loglike[iter] <- lfun(gnew, rbind(qnew[-nrow(qnew), ], updateq), f)
        if (loglike[iter] < loglike[iter - 1]) {
            qnew <- rbind(qnew[-nrow(qnew), ], qvector[, 3])
            loglike[iter] <- lfun(gnew, qnew, f)
        } else {
            qnew <- rbind(qnew[-nrow(qnew), ], updateq)
            loglike[iter] <- lfun(gnew, qnew, f)
        }
        if (abs(loglike[iter] - loglike[iter - 1]) < tol) {
            break
        } else {
            if (method == "EM") {
                fitmodel <- fFixEm(gnew, qnew, f, acc = T, max.iter = 3)
            } else {
                fitmodel <- fFixBr(gnew, qnew, f, acc = T, max.iter = 3)
            }
            qvector <- fitmodel$qvector
        }
    }
    qnew <- round(qnew, 5)
    qnew <- as.data.frame(qnew)
    rownames(qnew) <- "result"
    if(is.null(pubdata)) {
        colnames(qnew) <- NULL
    } else if (pubdata == "E11") {
        colnames(qnew) <- c("African", "European", "India", "Malay",
                            "SouthChineseDai", "SouthwestChineseYi",
                            "EastChinese", "Japanese", "NorthChineseOroqen",
                            "Yakut", "American")
    } else if (pubdata == "K13") {
        colnames(qnew) <- c("Siberian", "Amerindian", "West_African",
                            "Palaeo_African", "Southwest_Asian",
                            "East_Asian", "Mediterranean", "Australasian",
                            "Arctic", "West_Asian", "North_European",
                            "South_Asian", "East_African")
    }
    return(list(q = qnew, f = f, loglike = loglike[1:iter]))
}
