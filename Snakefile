import os
import glob
import sys

configfile: "config.yaml"

rscript = config['rscript']
ncells = config['ncells']
ntry = config['ntry']
nbatch = config['nbatch']

SAMPLES = config['samples']
BATCHES = range(1, nbatch+1)
T0 = [0, 5, 10, 15]


rule all:
    input: 
      expand("fits/fitt_scramble_t.{t0}_tau1_k1_k2_tau2.rds", t0=T0),
      expand("fits/fits_{sample}_ref1_t0_tau1_k1_k2_tau2.rds", sample=SAMPLES),
      expand("fits/fits_{sample}_ref1_t0_tau1_k1_k2.rds", sample=SAMPLES),
      ["fits/fits_RAD21_tau2_15_ref1_t0_tau1_k1_k2.rds"],
      expand("bootstrap/boot_{sample}_{batch}.pars", sample=SAMPLES, batch=BATCHES)


####################################################################

rule fit_scramble_t0:
    input: "data/scramble.csv"
    params:
      t0 = "-{t0}"
    output: "fits/fitt_scramble_t.{t0}_tau1_k1_k2_tau2.rds"
    threads: 8
    log: "logs/fit_scramble_{t0}.log"
    shell:
        "{rscript} R/fitting.R {input} {output} 1 {ncells} {ntry} {params.t0} 4 8 >2 {log}" 


####################################################################

rule fit_all:
    input: "data/{sample}.csv"
    output: "fits/fits_{sample}_ref1_t0_tau1_k1_k2_tau2.rds"
    threads: 8
    log: "logs/fit_{sample}.log"
    shell:
        "{rscript} R/fitting.R {input} {output} 1 {ncells} {ntry} 0 4 8 >2 {log}" 

####################################################################

rule fit_all_3par:
    input: "data/{sample}.csv"
    output: "fits/fits_{sample}_ref1_t0_tau1_k1_k2.rds"
    threads: 8
    log: "logs/fit_3par_{sample}.log"
    shell:
        "{rscript} R/fitting.R {input} {output} 1 {ncells} {ntry} 0 3 8 >2 {log}" 

####################################################################

rule fit_RAD21_t15:
    input: "data/RAD21.csv"
    output: "fits/fits_RAD21_tau2_15_ref1_t0_tau1_k1_k2.rds"
    threads: 8
    log: "logs/fit_RAD21_t15.log"
    shell:
        "{rscript} R/fitting.R {input} {output} 1 {ncells} {ntry} 0 3 15 >2 {log}" 

####################################################################

rule boot_all:
    input: "data/{sample}.csv"
    output: "bootstrap/boot_{sample}_{batch}.pars"
    threads: 8
    log: "logs/boot_{sample}.log"
    shell:
        "{rscript} R/bootstrap.R {input} {output} 1 {ncells} {ntry} 0 4 >2 {log}" 


