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


#Esta linea de codigo esta para borrar xd
#Codon1_Spneumoniae_coronada <- splitseq(Codon1_Spneumoniae,frame=0,word=3)