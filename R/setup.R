library(knitr)
library(ggplot2)
library(gridExtra)
#library(latex2exp)
library(kableExtra)
library(parallel)
library(methods)

#qq.options(code.pattern = "\\$\\{CODE\\}")
N <- function(n) prettyNum(n, big.mark = ",")

# Instructions to run this document
#
# Most of the code in this document is calculated locally, however some
# computationally-intensive stuff is done on a computer cluster. The cluster is
# accessed via ssh, using host name and user name provided below. The ssh
# connection needs to be setup in a way that password is not required (that is
# with the public key copied to the remote server). Shell scripts are copied to
# the remote file sytem and sent to the cluster using qsub command with queue
# names defined below. These need to be changed according to local needs.
#
# Code chunks containing qsub calls are marked as "eval=FALSE". You need to run
# them manually, one by one, and wait for the code to finish before continuing
# to the next step. Once all of them are done, the entire document can be
# knitted quickly.
#
# Remote file system is accessed from the local machine via Samba, or simliar,
# mount. Hence, there are two sets of directories, local to the cluster and
# remote, from the point of view of the machine where this document is built.
# Again, these need to be specified according to the local needs.

# Project directory on the cluster filesystem
clusterProjectDir <- "/cluster/gjb_lab/mgierlinski/projects/chromcomR/"

# The same project directory as seen from the local machine, via Samba mount. If
# you run RStudio on the same filesystem, use the same dir as clusterProjectDir
remoteProjectDir <- "/home/mgierlinski/projects/chromcomR/"

# Public HTML for file downloads
public_html <- "http://www.compbio.dundee.ac.uk/user/mgierlinski/chromcom/"

# Options for qsub to cluster
my.username <- "mgierlinski"
my.host <- "login-gjb-1.compbio.dundee.ac.uk"

# Rscript location
rscript <- "/cluster/gjb_lab/mgierlinski/software/miniconda2/envs/work/bin/Rscript"
