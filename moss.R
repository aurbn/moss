####### Parameters #######
FC2_TH = 1
FC15_TH = 1
TH = FC15_TH
##########################

group <- function(mrna, prot)
{
    if (abs(mrna) > TH)
        if (mrna > 0)
            m <- 'M'
        else
            m <- 'm'
    else
        m <- '0'
    
    if (abs(prot) > TH)
        if (prot > 0)
            p <- 'P'
        else
            p <- 'p'
    else
        p <- '0'
    
    paste0(m, p)
}

mrna <- read.table("transcripts.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
mrna <- data.frame(Gene = mrna$Gene,# TAIR = mrna$TAIR,
                  PP = mrna$Protoplasts_FPKM,
                  PN = mrna$Protonema_FPKM,
                  stringsAsFactors = FALSE)

mrna <- mrna[grepl("^Pp.+", mrna$Gene),]
mrna$mrnaPNtoPP <- mrna$PN/mrna$PP
mrna <- na.omit(mrna)  # Fix it

prot <- read.table("proteins.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
prot <- data.frame(Gene = prot$ProteinName, protPNtoPP = prot$PNtoPP_1.fc, stringsAsFactors = FALSE)
prot <- prot[grepl("^Pp.+", prot$Gene),]
prot$Gene <- sapply(strsplit(prot$Gene, split = '\\.'), "[", 1)
prot <- aggregate(. ~ Gene, data = prot, FUN = sum)

total <- merge(mrna[,c("Gene", "mrnaPNtoPP")], prot, by = "Gene" )#, all.y = TRUE)

total$mrnafc <- log(total$mrnaPNtoPP, base = 1.5)
total$protfc <- log(total$protPNtoPP, base = 1.5)

total$group <- as.factor(apply(total[,c("mrnafc", "protfc")], 1,
                               function(x) group(x[1], x[2])))

total <- total[ abs(total$mrnafc) > TH |
                abs(total$protfc) > TH,]

plot(total$mrnafc, total$protfc, pch = 19, col = total$group)
