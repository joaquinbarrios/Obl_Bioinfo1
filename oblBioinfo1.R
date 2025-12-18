#Obl Bioinfo
install.packages('seqinr')
library(seqinr)

#Obtener las secuencias necesarias
CglutamicumGenome <- read.fasta("Cglutamicum.fna")
CglutamicumCDS <- read.fasta('Cglutamicum.cds')
CglutamicumGB <- read.fasta("Cglutamicum.gbk")
CglutamicumPEP <- read.fasta('Cglutamicum.pep')
SpneumoniaeGenome <- read.fasta('Spneumoniae.fna')
SpneumoniaeCDS <- read.fasta('Spneumoniae.cds')
SpneumoniaeGB <- read.fasta('Spneumoniae.gbk')
SpneumoniaePEP <- ('Spneumoniae.pep')

#Calculo largos de genomas
LargoGenomaCglutamicum <- length(unlist(CglutamicumGenome))
LargoGenomaSpneumoniae <- length(unlist(SpneumoniaeGenome))

#Calculo largo de CDS
Cant_CDS_Cglutamicum <- length(CglutamicumCDS)
Cant_CDS_Spneumoniae <- length(SpneumoniaeCDS)

#Calculo porcentage CG de genomas
GC_Cglutamicum <- GC(unlist(CglutamicumGenome))
GC_Spneumonia <- GC(unlist(SpneumoniaeGenome))

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

boxplot(GC.Cglutamicum, GC.Spneumoniae, names=c('cglutamicum', 'Spneumoniae'),ylab='%GC',col = c('lightblue','lightyellow'))

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
  names = c("GC", 
            "GC1 ", 
            "GC2 ",
            "GC3"),
  ylab = "%GC",
  main = "Contenido GC total y posicional en CDSs de C.glutamicum",
  col = c("lightblue")
)

boxplot(
  GC.Spneumoniae, 
  GC1.Spneumoniae, 
  GC2.Spneumoniae,
  GC3.Spneumoniae, 
  names = c("GC ", 
            "GC1 ", 
            "GC2 ",
            "GC3 "),
  ylab = "%GC",
  main = "Contenido GC total y posicional en CDSs de S.pneumoniae",
  col = c("lightyellow")
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
  main = "Contenido GC total y posicional en CDSs",
  col = c('lightblue',"lightyellow")
)


#USO DE CODONES

Cglu_all <- do.call(c, CglutamicumCDS)
Spneu_all <- do.call(c, SpneumoniaeCDS)


RSCU_Cglutamicum <- uco(Cglu_all, index = "freq")
RSCU_Spneumoniae <- uco(Spneu_all, index = "freq")

#STOP
names(RSCU_Cglutamicum) <- tabla_rscu[,1]
names(RSCU_Spneumoniae) <- tabla_rscu[,1]

aa_gCg <- tapply(RSCU_Cglutamicum, names(RSCU_Cglutamicum), sum)
aa_gSp <- tapply(RSCU_Spneumoniae, names(RSCU_Spneumoniae), sum)


barplot(
  rbind(aa_gCg, aa_gSp),
  beside = TRUE,
  col=c("blue","yellow"),
  ylab = "Frecuencias de aminoacidos",
  xlab = "Aminoácido",
  legend.text = c("Cglu", "Spne"))

names(RSCU_Cglutamicum) <- sapply(
  names(RSCU_Cglutamicum),
  function(codon) {
    translate(s2c(toupper(codon)))
  }
)



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

#ACA ARRANCA EL BARPLOT
ORILOC
ol <- oriloc(seq.fasta = CglutamicumGenome, gbk = "Cglutamicum.gbk")

draw.oriloc(ol, main = "Grafico Oriloc")
 

library(dplyr)

blast12 <- read.table("sp1_vs_sp2.tsv", sep = "\t", stringsAsFactors = FALSE)
blast21 <- read.table("sp2_vs_sp1.tsv", sep = "\t", stringsAsFactors = FALSE)

colnames(blast12) <- c("query", "subject", "pident", "length", "evalue", "bitscore")
colnames(blast21) <- c("query", "subject", "pident", "length", "evalue", "bitscore")

best12 <- blast12 %>%
  group_by(query) %>%
  slice_max(bitscore, n = 1) %>%
  ungroup()

best21 <- blast21 %>%
  group_by(query) %>%
  slice_max(bitscore, n = 1) %>%
  ungroup()

brh <- inner_join(
  best12,
  best21,
  by = c("query" = "subject", "subject" = "query")
)

orthologs <- brh %>%
  select(gene_sp1 = query.x,
         gene_sp2 = subject.x)

write.table(orthologs,
            "BRH_orthologs.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


TABLAORTOLOGOSOK <- read.table('BRH_orthologs.tsv')

length(TABLAORTOLOGOSOK)

TABLAORTOLOGOSOK1 <-TABLAORTOLOGOSOK %>%
  rename(
    query    = V1,
    subject  = V2,
    pident   = V3,
    length   = V4,
    evalue   = V5,
    bitscore = V6
  )

#Bombini guzini

Top_Ten <- read.fasta('allfiles.fasta')

total_allfiles <- sum(sapply(top10alineadas, length))

nuevoall <- read.fasta('seqdump2.fasta')
total_allfiles <- sum(sapply(nuevoall, length))

top10alineadas <- read.fasta('top10alinedas.fasta')

length(top10alineadas[[2]])
CglutamicumCDS[[1]]

tabla_Cglutamicum <- data.frame(
  CDS = names(CglutamicumCDS),
  G1_gc  = GC.Cglutamicum,
  G1_gc1 = GC1.Cglutamicum,
  G1_gc2 = GC2.Cglutamicum,
  G1_gc3 = GC3.Cglutamicum
)

tabla_Spneumoniae <- data.frame(
  CDS = names(SpneumoniaeCDS),
  G1_gc  = GC.Spneumoniae,
  G1_gc1 = GC1.Spneumoniae,
  G1_gc2 = GC2.Spneumoniae,
  G1_gc3 = GC3.Spneumoniae
)
