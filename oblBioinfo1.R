#Obl Bioinfo
install.packages('seqinr')
library(seqinr)

#Obtener las secuencias necesarias
CglutamicumGenome <- read.fasta("Cglutamicum.fna")
CglutamicumCDS <- read.fasta('Cglutamicum.cds')
CglutamicumGB <- read.fasta("Cglutamicum.gbff")
CglutamicumPEP <- read.fasta('Cglutamicum.pep')
SpneumoniaeGenome <- read.fasta('Spneumoniae.fna')
SpneumoniaeCDS <- read.fasta('Spneumoniae.cds')
SpneumoniaeGB <- read.fasta('Spneumoniae.gbff')
SpneumoniaePEP <- ('Spneumoniae.pep')

#Calculo largos de genomas
LargoGenomaCglutamicum <- length(unlist(CglutamicumGenome))
LargoGenomaSpneumoniae <- length(unlist(SpneumoniaeGenome))

#Calculo largo de CDS
Cant_CDS_Cglutamicum <- length(CglutamicumCDS)
Cant_CDS_Spneumoniae <- length(SpneumoniaeCDS)

#Calculo porcentage CG de genomas
GC_Cglutamicum <- GC(unlist(CglutamicumGenome))
GC_Spneumonia <- GC(unlist(SpneumoniaGenome))

#Calculo densidades de secuencuencias codificantes
densidad_CDS_Cglutamicum <- Cant_CDS_Cglutamicum/LargoGenomaCglutamicum
densidad_CDS_Spneumoniae <- Cant_CDS_Spneumoniae/LargoGenomaSpneumoniae


#Creo dataframe que va a contener la informacion
gc_genomas <- data.frame(Genoma = 0, GC = 0)

gc_genomas[1,1] <- "Cglutamicum"
gc_genomas[2,1] <- "Spneumonia"

#Calculo contenido GC y lo agrego al DF
gc_genomas[1,2] <- unlist(lapply(CglutamicumGenome, GC))
gc_genomas[2,2] <- unlist(lapply(SpneumoniaeGenome, GC))

GC.Cglutamicum <- sapply(CglutamicumCDS, GC)
GC.Spneumoniae <- sapply(SpneumoniaeCDS, GC)

GC.Cglutamicum <- GC.Cglutamicum * 100
GC.Spneumoniae <- GC.Spneumoniae * 100

#Boxplot de Porcentage GC en genoma 

boxplot(GC.Cglutamicum, GC.Spneumoniae, names=c('cglutamicum', 'Spneumoniae'),ylab='%GC')

#Hasta aca hice cada cds cuanto tiene de GC. Cada punto en la grafica es un cds.

#Ahora tengo que hacer lo mismo, pero solo para la primera posicion de cada AA.

#Osea, tengo que tener una lista que determine los AAA, osea ir separando de a 3 letras.

GC_pos <- function(seq, pos) {
  bases <- seq[seq(pos, length(seq), by = 3)]
  sum(bases %in% c("g", "c")) / length(bases) * 100
}



GC1.Cglutamicum <- sapply(CglutamicumCDS, GC_pos, pos = 1)
GC2.Cglutamicum <- sapply(CglutamicumCDS, GC_pos, pos = 2)
GC3.Cglutamicum <- sapply(CglutamicumCDS, GC_pos, pos = 3)

GC1.Spneumoniae <- sapply(SpneumoniaeCDS, GC_pos, pos = 1)
GC2.Spneumoniae <- sapply(SpneumoniaeCDS, GC_pos, pos = 2)
GC3.Spneumoniae <- sapply(SpneumoniaeCDS, GC_pos, pos = 3)


boxplot(
  GC.Cglutamicum, 
  GC1.Cglutamicum, 
  GC2.Cglutamicum,
  GC3.Cglutamicum, 
  names = c("GC G1", 
            "GC1 G1", 
            "GC2 G1",
            "GC3 G1"),
  ylab = "%GC",
  main = "Contenido GC total y posicional en CDSs"
)

boxplot(
  GC.Spneumoniae, 
  GC1.Spneumoniae, 
  GC2.Spneumoniae,
  GC3.Spneumoniae, 
  names = c("GC G2", 
            "GC1 G2", 
            "GC2 G2",
            "GC3 G2"),
  ylab = "%GC",
  main = "Contenido GC total y posicional en CDSs"
)

boxplot(
  GC.Cglutamicum, GC.Spneumoniae,
  GC1.Cglutamicum, GC1.Spneumoniae,
  GC2.Cglutamicum, GC2.Spneumoniae,
  GC3.Cglutamicum, GC3.Spneumoniae,
  names = c("GC G1", "GC G2",
            "GC1 G1", "GC1 G2",
            "GC2 G1", "GC2 G2",
            "GC3 G1", "GC3 G2"),
  ylab = "%GC",
  main = "Contenido GC total y posicional en CDSs"
)


#USO DE CODONES

Cglu_all <- do.call(c, CglutamicumCDS)
Spneu_all <- do.call(c, SpneumoniaeCDS)


RSCU_Cglutamicum <- uco(Cglu_all, index = "rscu")
RSCU_Spneumoniae <- uco(Spneu_all, index = "rscu")

#Verifico que tengo para las 64 combinaciones
length(RSCU_Cglutamicum)  
length(RSCU_Spneumoniae)

#Obtengo tabla
tabla_rscu <- data.frame(
  codon  = names(RSCU_Cglutamicum),
  uco_G1 = as.numeric(RSCU_Cglutamicum),
  uco_G2 = as.numeric(RSCU_Spneumoniae),
  stringsAsFactors = FALSE
)

tabla_rscu$AA <- sapply(tabla_rscu$codon, function(cod) {
  translate(s2c(cod))
})

tabla_rscu <- tabla_rscu[order(tabla_rscu$AA, tabla_rscu$codon), ]

library(tidyr)
library(ggplot2)

par(mfrow = c(1, 2))  # 1 fila, 2 columnas

# Genoma 1
barplot(tabla_rscu$uco_G1,
        names.arg = tabla_rscu$codon,
        ylim = c(0, max(tabla_rscu$uco_G1, tabla_rscu$uco_G2)),
        ylab = "RSCU",
        main = "RSCU R Genoma 1",
        col = "grey",
        las = 1)

# Genoma 2
barplot(tabla_rscu$uco_G2,
        names.arg = tabla_rscu$codon,
        ylim = c(0, max(tabla_rscu$uco_G1, tabla_rscu$uco_G2)),
        ylab = "RSCU",
        main = "RSCU R Genoma 2",
        col = "grey",
        las = 1)

par(mfrow = c(1, 1))

#Limitados

tabla_LRS <- subset(tabla_rscu, AA %in% c("L", "R", "S"))
tabla_LRS <- tabla_LRS[order(tabla_LRS$AA, tabla_LRS$codon), ]

bp1 <- barplot(
  tabla_LRS$uco_G1,
  names.arg = tabla_LRS$codon,
  las = 2,
  col = "grey",
  ylab = "RSCU",
  main = "RSCU Leu / Arg / Ser – Genoma 1"
)

bp2 <- barplot(
  tabla_LRS$uco_G2,
  names.arg = tabla_LRS$codon,
  las = 2,
  col = "grey",
  ylab = "RSCU",
  main = "RSCU Leu / Arg / Ser – Genoma 2"
)

line_pos <- bp1[cumsum(table(tabla_LRS$AA))] + diff(bp)[1]/2
abline(v = line_pos, lty = 2)

line_pos <- bp2[cumsum(table(tabla_LRS$AA))] + diff(bp)[1]/2
abline(v = line_pos, lty = 2)
.

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

 













