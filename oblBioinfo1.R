#Obl Bioinfo
install.packages('seqinr')
library(seqinr)

CglutamicumGenome <- read.fasta("GitHub/Obl_Bioinfo1/Cglutamicum.fna")
CglutamicumCDS <- read.fasta('GitHub/Obl_Bioinfo1/Cglutamicum.cds')
CglutamicumGB <- read.fasta("GitHub/Obl_Bioinfo1/Cglutamicum.gbff")
CglutamicumPEP <- read.fasta('GitHub/Obl_Bioinfo1/Cglutamicum.pep')
SpneumoniaeGenome <- read.fasta('GitHub/Obl_Bioinfo1/Spneumoniae.fna')
SpneumoniaeCDS <- read.fasta('GitHub/Obl_Bioinfo1/Spneumoniae.cds')
SpneumoniaeGB <- read.fasta('GitHub/Obl_Bioinfo1/Spneumoniae.gbff')
SpneumoniaePEP <- ('GitHub/Obl_Bioinfo1/Spneumoniae.pep')

LargoGenomaCglutamicum <- length(unlist(CglutamicumGenome))
LargoGenomaSpneumoniae <- length(unlist(SpneumoniaeGenome))

Cant_CDS_Cglutamicum <- length(CglutamicumCDS)
Cant_CDS_Spneumoniae <- length(SpneumoniaeCDS)

GC_Cglutamicum <- GC(unlist(CglutamicumGenome))
GC_Spneumonia <- GC(unlist(SpneumoniaGenome))

Codon1_Cglutamicum<- CglutamicumCDS[[1]]
Codon2_Cglutamicum<- CglutamicumCDS[[2]]
Codon3_Cglutamicum<- CglutamicumCDS[[3]]

Codon1_Spneumoniae<- SpneumoniaeCDS[[1]]
Codon2_Spneumoniae<- SpneumoniaeCDS[[2]]
Codon3_Spneumoniae<- SpneumoniaeCDS[[3]]
  
#Codon1_Spneumoniae_coronada <- splitseq(Codon1_Spneumoniae,frame=0,word=3)

GC_Cglutamicum_CD1 <- GC(unlist(Codon1_Cglutamicum))
GC_Cglutamicum_CD2 <- GC(unlist(Codon2_Cglutamicum))
GC_Cglutamicum_CD3 <- GC(unlist(Codon3_Cglutamicum))


GC_SpneumoniaGB_CD1 <- GC(unlist(Codon1_Spneumoniae))
GC_SpneumoniaGB_CD2 <- GC(unlist(Codon2_Spneumoniae))
GC_SpneumoniaGB_CD3 <- GC(unlist(Codon3_Spneumoniae))

boxplot(GC_Cglutamicum_CD1, GC_Cglutamicum_CD2, GC_Cglutamicum_CD3, names=c('Codon 1', 'Codon 2', 'Codon 3'))

Uco_Spneumoniae1 <- uco(SpneumoniaeCDS[[1]])
Uco_Spneumoniae2 <- uco(SpneumoniaeCDS[[2]])
Uco_Spneumoniae3 <- uco(SpneumoniaeCDS[[3]])

Uco_Cglutamicum1 <- uco(CglutamicumCDS[[1]])
Uco_Cglutamicum2 <- uco(CglutamicumCDS[[2]])
Uco_Cglutamicum3 <- uco(CglutamicumCDS[[3]])

Uco_CglutamicumGENOME <- uco(unlist(CglutamicumCDS), index="rscu")


Uco_SpneumoniaeGENOME <- uco(unlist(SpneumoniaeGenome))

Cglutamicum_AA <- translate(unlist(CglutamicumCDS))

Spneumoniae_AA <- translate(unlist(SpneumoniaeGenome))

nAAA <- aaa(Cglutamicum_AA)

names(Uco_CglutamicumGENOME) <- nAAA

vectArg <- grep("Arg",nAAA)

cArg <- Uco_CglutamicumGENOME[vectArg]

