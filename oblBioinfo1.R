#Obl Bioinfo
install.packages('seqinr')
library(seqinr)

CglutamicumGenome <- read.fasta("GitHub/Obl_Bioinfo1/Cglutamicum.fna")
CglutamicumCDS <- read.fasta('GitHub/Obl_Bioinfo1/Cglutamicum.cds')
CglutamicumGB <- read.fasta("GitHub/Obl_Bioinfo1/Cglutamicum.gbff")
CglutamicumPEP <- read.fasta('GitHub/Obl_Bioinfo1/Cglutamicum.pep')
SpneumoniaeGenome <- read.fasta('Spneumoniae.fna')
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

gc_genomas <- data.frame(Genoma = 0, GC = 0)

gc_genomas[1,1] <- "Cglutamicum"
gc_genomas[2,1] <- "Spneumonia"

gc_genomas[1,2] <- unlist(lapply(CglutamicumGenome, GC))
gc_genomas[2,2] <- unlist(lapply(SpneumoniaeGenome, GC))

GC.Cglutamicum <- sapply(CglutamicumCDS, GC)
GC.Spneumoniae <- sapply(SpneumoniaeCDS, GC)

boxplot(GC.Cglutamicum, GC.Spneumoniae, names=c('cglutamicum', 'Spneumoniae'),ylab='%GC')

#Hasta aca hice cada cds cuanto tiene de GC. Cada punto en la grafica es un cds.

#Ahora tengo que hacer lo mismo, pero solo para la primera posicion de cada AA.

#Osea, tengo que tener una lista que determine los AAA, osea ir separando de a 3 letras.








.























Uco_Cglutamicumcds <- uco(unlist(CglutamicumCDS,recursive = 'FALSE'))

lista_cds <- unlist(CglutamicumCDS)

Uco_SpneumoniaeGENOME <- uco(unlist(SpneumoniaeCDS))

Cglutamicum_AA <- translate(unlist(CglutamicumCDS))

Spneumoniae_AA <- translate(unlist(SpneumoniaeGenome))

nAAA <- aaa(Cglutamicum_AA)



names(Uco_CglutamicumGENOME) <- nAAA

n1 <- sapply(names(Uco_Cglutamicumcds),s2c)

nombres_CortosAA_CglutamicumGenome <- translate(n1)

aaa(nomrbres_CortosAA_CglutamicumGenome)
vectArg <- grep("Arg",nAAA)

cArg <- Uco_CglutamicumGENOME[vectArg]



uco_Cglutamicum_s <- sapply(CglutamicumCDS, function(x){uco(x, index = "rscu")})








Uco_Cglutamicumcds <- uco(unlist(CglutamicumCDS), index="rscu")

qweqwe <- translate(sapply(names(Uco_Cglutamicumcds),s2c))

names(Uco_Cglutamicumcds) <- aaa(qweqwe)

vecAla <-grep("Ala",aaa(qweqwe))

 













