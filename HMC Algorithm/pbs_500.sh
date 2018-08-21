#!/bin/bash

#PBS -l nodes=1:ppn=4

#PBS -l walltime=48:00:00
#PBS -d /home/macdona8/Dissertation/"HMC Algorithm"



R CMD BATCH Optimal_Accept_Rate_500.R
