#!/usr/bin/env Rscript

argv<-commandArgs(TRUE)
d <- read.table(argv[1])

fisher.test(d)
