library(knitr)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(parallel)
library(methods)
#library(GetoptLong)

#qq.options(code.pattern = "\\$\\{CODE\\}")
N <- function(n) prettyNum(n, big.mark = ",")

# top directory for the project
projectDir <- "/home/mgierlinski/projects/chromcomR/"

# Public HTML for file downloads
public_html <- "http://www.compbio.dundee.ac.uk/user/mgierlinski/chromcom/"

