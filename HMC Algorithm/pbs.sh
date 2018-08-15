#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=02:00:00

cd ~/Dissertation/"HMC Algorithm"

R CMD BATCH Optimal_Accept_Rate.R
