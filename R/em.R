#' @title Do ancestry analysis with EM algorithm
#' @description The EM algorithm could be used for estimating the Q and F matrix.
#' @usage em(g, q, f, acc, max.iter, tol, model)
#' @param g Genotype matrix with dimensions \eqn{n × p}, where n is sample size
#' and p is the number of SNPs.
#' @param q Ancestry coefficient matrix with dimensions \eqn{n × K}, where n
#' is sample size and K is the number of populations.
#' @param f Minor allele frequency matrix with dimensions \eqn{K × p},
#' where K is the number of populations and p is the number of SNPs.
#' @param acc a logical value indicating whether use accelerated EM or not.
#' @param max.iter an integer. If acc is TRUE, the number of iterations must be set.
#' @param tol Tolerance. If acc is FALSE, tol must be set. The default is 1e-4.
#' @param model Choose which model you want to use. Supervised learning or unsupervised learning.
#' @return Estimation results of q, f and the loglikelihood value of each iteration.
#' @export

# The EM algorithm
em <- function(g, q, f, acc, max.iter = 3, tol = 1e-4,
               model = c("supervised", "unsupervised")) {
    if (!is.matrix(g)) {
        stop("g should be a matrix")
    }
    if (!is.matrix(q)) {
        stop("q should be a matrix")
    }
    if (!is.matrix(f)) {
        stop("f should be a matrix")
    }
    p <- ncol(f)
    n <- nrow(q)
    K <- ncol(q)
    if (acc == TRUE) {
        loglikem <- numeric()
        loglikem[1] <- lfun(g, q, f)
        flongvec1 <- matrix(NA, nrow(f) * ncol(f), max.iter)
        flongvec1[, 1] <- matrix(f, nrow(f) * ncol(f), 1)
        if (model == "supervised") {
            qvec <- matrix(NA, ncol(q), max.iter)
            qvec[, 1] <- t(q[nrow(q), ])
        }
        if (model == "unsupervised") {
            qvec <- matrix(NA, nrow(q) * ncol(q), max.iter)
            qvec[, 1] <- matrix(q, nrow(q) * ncol(q), 1)
        }
        iter <- 1
        repeat{
            iter <- iter + 1
            # memory the current q and f
            q0 <- q
            f0 <- f
            if (model == "supervised") {
                q <- em.q(g, q0, q, f0, p, n)
                qvec[, iter] <- t(q[nrow(q), ])
            }
            if (model == "unsupervised") {
                for (i in 1:nrow(q0)) {
                    q[i, ] <- em.q(g, q0, q, f0, p, i)[i, ]
                }
                qvec[, iter] <- matrix(q, nrow(q) * ncol(q), 1)
            }
            # update F
            f <- .Call("updatef", n, p, K, q0, f0, g)
            f[f < lb] <- lb
            f[f > ub] <- ub
            flongvec1[, iter] <- matrix(f, nrow(f) * ncol(f), 1)
            loglikem[iter] <- lfun(g, q, f)
            if (iter == max.iter)
                break
        }
        if (model == "supervised") {
            res <- list(qvector = qvec, fvector = flongvec1,
                        loglike = loglikem, q = q, f = f)
        } else {
            res <- list(qvector = qvec, fvector = flongvec1,
                        loglike = loglikem, q = q, f = f)
        }
        return(res)
    }
    if (acc == FALSE) {
        if (is.null(tol)) {
            stop("you must set up the tolerance")
        }
        loglikem <- numeric()
        loglikem[1] <- lfun(g, q, f)
        iter <- 1
        repeat{
            iter <- iter + 1
            # memory the current q and f
            q0 <- q
            f0 <- f
            # update Q
            if (model == "supervised") {
                q <- em.q(g, q0, q, f0, p, n)
            }
            if (model == "unsupervised") {
                for (i in 1:nrow(q0)) {
                    q[i, ] <- em.q(g, q0, q, f0, p, i)[i, ]
                }
            }
            # update F
            f <- .Call("updatef", n, p, K, q0, f0, g)
            loglikem[iter] <- lfun(g, q, f)
            if (abs(loglikem[iter] - loglikem[iter - 1] < tol))
                break
        }
        return(list(f = f, loglike = loglikem[1:iter], q = q))
    }
}

#' @title EM when f is fixed
#' @description This function can be used for ancestry analysis when frequency matrix is fixed.
#' @usage fFixEm(gnew, qnew, f, acc, tol = 1e-4, pubdata)
#' @param gnew Genotype matrix. The number of row present in gnew is 1 and the number
#' of column is the number of SNPs.
#' @param qnew Initial q used in calculation. A vector. sum(q) must be 1.
#' @param f Allele frequencies learned from the reference panels.
#' @param acc a logical value indicating whether use quasi-Newton accelerated EM or not.
#' @param max.iter an integer. If acc is TRUE, the number of iterations must be set.
#' @param tol Tolerance. If acc is FALSE, tol must be set. The default is 1e-4.
#' @param pubdata You can choose a public dataset here, E11 or K13. You also can use other public
#' dataset which is not in this package.
#' @return Estimation results of q and the loglikelihood value of each iteration.
#' @export

fFixEm <- function(gnew, qnew, f, acc, max.iter = 3,
                   tol = 1e-4, pubdata = NULL) {
    if (is.null(acc)) {
        stop("You must set up a logical value for acc!")
    }
    p <- ncol(f)
    n <- nrow(qnew)
    if (acc == TRUE) {
        loglikem <- numeric()
        loglikem[1] <- lfun(gnew, qnew, f)
        qvec <- matrix(NA, ncol(qnew), max.iter)
        qvec[, 1] <- t(qnew[nrow(qnew), ])
        iter <- 1
        repeat {
            iter <- iter + 1
            # memory the current q and f
            q0 <- qnew
            qnew <- em.q(gnew, q0, qnew, f, p, n)
            qvec[, iter] <- t(qnew[nrow(qnew), ])
            loglikem[iter] <- lfun(gnew, qnew, f)
            if (iter == max.iter)
                break
        }
        return(list(qvector = qvec, loglike = loglikem, q = qnew))
    }
    if (acc == FALSE) {
        loglikem <- numeric()
        loglikem[1] <- lfun(gnew, qnew, f)
        iter <- 1
        repeat {
            iter <- iter + 1
            # memory the current q
            q0 <- qnew
            # update Q
            qnew <- em.q(gnew, q0, qnew, f, p, n)
            loglikem[iter] <- lfun(gnew, qnew, f)
            if (abs(loglikem[iter] - loglikem[iter - 1] < tol))
                break
        }
        qnew <- round(qnew, 5)
        qnew <- as.data.frame(qnew)
        rownames(qnew) <- "result"
        if (is.null(pubdata)) {
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
        return(list(loglike = loglikem[1:iter], q = qnew))
    }
}
