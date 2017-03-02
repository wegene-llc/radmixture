#' Block Relaxation for ancestry estimation
#' @description This function is also used for estimating Q and F but faster than EM.
#' @param g Genotype matrix with dimensions \eqn{n × p}, where n is sample size
#' and p is the number of SNPs.
#' @param q Ancestry coefficient matrix with dimensions \eqn{n × K}, where n
#' is sample size and K is the number of populations.
#' @param f Minor allele frequency matrix with dimensions \eqn{K × p},
#' where K is the number of populations and p is the number of SNPs.
#' @param acc a logical value indicating whether use quasi-Newton accelerated BR or not.
#' @param max.iter If acc = T, max.iter must be set, the default is 3.
#' @param tol Tolerance, if acc = F, tolerance must be set, the default is 1e-4.
#' @param model Choose which model you want to use. Supervised learning or unsupervised learning.
#' @return Estimation results of q, f and the loglikelihood value of each iteration.
#' @export

# block relaxation (sequential quadratic programming)
br <- function(g, q, f, acc, max.iter = 3, tol = NULL,
               model = c("supervised", "unsupervised")) {

   if (is.null(acc)) {
       stop("You must set a logical value for acc!")
   }
   if (acc == TRUE) {
        logbl <- numeric()
        logbl[1] <- lfun(g, q, f)
        flongvec1 <- matrix(NA, nrow(f) * ncol(f), max.iter)
        flongvec1[, 1] <- matrix(f, nrow(f) * ncol(f), 1)
        if (model == "supervised") {
            qlongvec1 <- matrix(NA, ncol(q), max.iter)
            qlongvec1[, 1] <- t(q[nrow(q), ])
        } else if (model == "unsupervised") {
            qlongvec1 <- matrix(NA, nrow(q) * ncol(q), max.iter)
            qlongvec1[, 1] <- matrix(q, nrow(q) * ncol(q), 1)
        }

        iter <- 1

        repeat {
            iter <- iter + 1
            f <- brf(g, q, f)
            flongvec1[, iter] <- matrix(f, nrow(f) * ncol(f), 1)
            q <- brq(g, q, f, model)
            if (model == "supervised") {
                qlongvec1[, iter] <- t(q[nrow(q), ])
            } else if (model == "unsupervised") {
                qlongvec1[, iter] <- matrix(q, nrow(q) * ncol(q), 1)
            }
            logbl[iter] <- lfun(g, q, f)
            if (iter == max.iter) {
                break
            }
        }
        res <- list(qvector = qlongvec1, fvector = flongvec1, loglike = logbl)
        return(res)
    }

    if (acc == FALSE) {
        if (is.null(tol)) {
            stop("You should set up tolerance!")
        }
        logbl <- numeric()
        logbl[1] <- lfun(g, q, f)
        iter <- 1
        repeat {
            iter <- iter + 1
            f <- brf(g, q, f)
            q <- brq(g, q, f, model)
            logbl[iter] <- lfun(g, q, f)
            if (abs(logbl[iter] - logbl[iter - 1]) < tol) {
                break
            }
        }
        if (model == "supervised") {
            res <- list(q = q[nrow(q), ], f = f, loglike = logbl[1:iter])
        } else {
            res <- list(q = q, f = f, loglike = logbl[1:iter])
        }
        return(res)
    }
}


#' block relaxation with fixed F
#' @param gnew Integer which length is the number of SNPs used in calculation.
#' @param qnew Initial q used in calculation. Numeric. Sum(q) must be 1.
#' @param f Allele frequencies learned from the reference panels.
#' @param acc a logical value indicating whether use quasi-Newton accelerated BR or not.
#' @param max.iter If acc = T, max.iter must be set, the default is 3.
#' @param tol If acc = F, tolerance must be set, the default is 1e-4.
#' @param pubdata You can choose a public dataset here, E11 or K13. You also can use other public
#' dataset which is not in this package.
#' @export

fFixBr <- function(gnew, qnew, f, acc, max.iter = 3,
                   tol = 1e-4, pubdata = NULL) {
    if (acc == TRUE) {
        logbl <- numeric()
        logbl[1] <- lfun(gnew, qnew, f)
        qvec <- matrix(NA, ncol(qnew), max.iter)
        qvec[, 1] <- t(qnew[nrow(qnew), ])
        iter <- 1
        repeat {
            iter <- iter + 1
            qnew <- brq(gnew, qnew, f, model = "supervised")
            qvec[, iter] <- t(qnew[nrow(qnew), ])
            logbl[iter] <- lfun(gnew, qnew, f)
            if (iter == max.iter) {
                break
            }
        }
        res <- list(qvector = qvec, q = qnew[nrow(qnew), ],
                    loglike = logbl[1:iter])
        return(res)
    } else {
        if (is.null(tol)) {
            stop("You should set up the tolerance!")
        }
        logbl <- numeric()
        logbl[1] <- lfun(gnew, qnew, f)
        iter <- 1
        repeat {
            iter <- iter + 1
            qnew <- brq(gnew, qnew, f, model = "supervised")
            logbl[iter] <- lfun(gnew, qnew, f)
            if (abs(logbl[iter] - logbl[iter - 1]) < tol) {
                break
            }
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
        res <- list(q = qnew, loglike = logbl[1:iter])
        return(res)
    }
}
