#Obl Bioinfo

library(seqinr)
setwd("/home/jovyan/Bioinformatica1_ORTee/Obl_Bioinfo1")

CglutamicumGenome <- read.fasta('Cglutamicum.fna')
CglutamicumCDS <- read.fasta('Cglutamicum.cds')
CglutamicumGB <- read.fasta('Cglutamicum.gbff')
SpneumoniaeGenome <- read.fasta('Spneumoniae.fna')
SpneumoniaeCDS <- read.fasta('Spneumoniae.cds')
SpneumoniaeGB <- read.fasta('Spneumoniae.gbff')
read.fasta('')
read.fasta('')

LargoGenomaCglutamicum <- length(unlist(CglutamicumGenome))
LargoGenomaSpneumoniae <- length(unlist(SpneumoniaeGenome))

Cant_CDS_Cglutamicum <- length(CglutamicumCDS)
Cant_CDS_Spneumoniae <- length(SpneumoniaeCDS)


