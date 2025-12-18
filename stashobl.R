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

# Sliding window

## Ventana deslizante de largo 10K y step 10K

ln_ecoli_ch1 <- length(ecoli_ch1)
w_size <- 10000
step <- 10000

# Redondeo "hacia arriba"
?ceiling

num_frag <- ceiling(ln_ecoli_ch1/w_size)
wins <- data.frame(
  fragment_index = seq(1, num_frag),
  from = seq(from = 1, to = ln_ecoli_ch1, by = w_size),
  to = c(seq(from = step, to = ln_ecoli_ch1, by = w_size), ln_ecoli_ch1) # Le agregamos la ultima pos
)

fragmentos <- list()

i=2
for(i in 1:num_frag){
  fragmentos[[i]] <- getFrag(ecoli_ch1, begin = wins[i,2], end = wins[i,3])
}


length(fragmentos)
lengths(fragmentos)