# radmixture

[![CRAN](https://www.r-pkg.org/badges/version/radmixture)](https://cran.r-project.org/package=radmixture)
![Downloads](https://cranlogs.r-pkg.org/badges/radmixture)
[![Build status](https://ci.appveyor.com/api/projects/status/wj5afxb1v42h8oui?svg=true)](https://ci.appveyor.com/project/wegene-llc/radmixture)

This is an introduction to radmixture which could help you with estimating individual ancestries from large SNP genotype data. 

This document introduces you how to use this package. 

### Installation

- Install the latest development version from GitHub with
```r
devtools::install_github("wegene-llc/radmixture")
```
- Install the latest version on CRAN with
```r
install.packages("radmixture")
```
First, you must prepare your raw data file as follow:

![](https://cloud.githubusercontent.com/assets/18478302/22725657/b4abd21a-ee09-11e6-9ef8-a4092be538e8.png)

and read it into R with `read.table` or `read.csv`.

```r
library(radmixture)
genotype <- read.table(file = '/path/to/file')
# genotype <- read.csv(file = 'path/to/file')
```

## Calculate your ancestry components with public dataset

<p><span style='color:red; font-size: 20px'><strong>Pay attention: </strong></span></p>
Since the public datasets are not on CRAN, if you download this package from CRAN, you should download the datasets that you need from GitHub by yourself. For example:

```r
download.file(url = 'https://github.com/wegene-llc/radmixture/raw/master/data/globe4.alleles.RData', destfile = '/path/to/globe4.alleles.RData')
download.file(url = 'https://github.com/wegene-llc/radmixture/raw/master/data/globe4.4.F.RData', destfile = '/path/to/globe4.4.F.RData')
```
And then, load them into R.
```r
load('/path/to/globe4.alleles.RData')
load('/path/to/globe4.4.F.RData')
```

### Example

Six public datasets could be used in this package.

- Use `tfrdpub()`to transfer your raw data to a format understood by radmixture
```r
# Use K4
res <- tfrdpub(genotype, 4, globe4.alleles, globe4.4.F)
# Use K7b
res <- tfrdpub(genotype, 7, K7b.alleles, K7b.7.F)
# Use world9
res <- tfrdpub(genotype, 9, world9.alleles, world9.9.F)
# Use E11
res <- tfrdpub(genotype, 11, e11.alleles, e11.11.F)
# Use K12b
res <- tfrdpub(genotype, 12, K12b.alleles, K12b.12.F)
# Use K13
res <- tfrdpub(genotype, 13, globe13.alleles, globe13.13.F)
```

- Use `fFixQN()` to calculate your result
```r
# Use K4
ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "K4")
# Use K7b
ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "K7b")
# Use world9
ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "World9")
# Use E11
ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "E11")
# Use K12b
ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "K12b")
# Use K13
ances <- fFixQN(res$g, res$q, res$f, tol = 1e-4, method = "BR", pubdata = "K13")
```
`ances$q` is your result.

### References

- D.H. Alexander, J. Novembre, and K. Lange. Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655–1664, 2009.
- H. Zhou, D. H. Alexander, and K.  Lange. A quasi-Newton method for accelerating the convergence of iterative optimization algorithms. Statistics and Computing, 2009.

### Contributions
We welcome contributions on radmixture. You can fork this repo and make your changes and submit a pull request.
Report bug issues on [issue page](https://github.com/wegene-llc/radmixture/issues).

### License

MIT + LICENSE

### TODO

- add tutorial
