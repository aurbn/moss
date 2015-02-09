require(ggplot2)
####### Parameters #######
FC2_TH = 1
FC15_TH = 1
TH = FC15_TH
SWATH_SOURCE = "raw"   # "raw" or "processed"
PN_SAMPLES = c("F163", "F164", "F165")
PP_SAMPLES = c("F166", "F167", "F168", "F169")
PROT_FILTER_PV = TRUE
PROT_REQ_PV = 0.05

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

tst <- function(row, k, e)
{
    a = length(na.omit(row[c(k)]))
    b = length(na.omit(row[c(e)]))
    t = list(p.value = NA) # If test cannot be performed
    try(
{
    t <- wilcox.test(as.numeric(na.omit(row[c(k)])),
                     as.numeric(na.omit(row[c(e)])),
                     correct = F, exact = F)
}, TRUE)
if(is.finite(t$p.value))
{
    return(t$p.value)
}
else
{
    return(1)
}
}

mrna <- read.table("transcripts.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
mrna <- data.frame(Gene = mrna$Gene,# TAIR = mrna$TAIR,
                  PP = mrna$Protoplasts_FPKM,
                  PN = mrna$Protonema_FPKM,
                  stringsAsFactors = FALSE)

mrna <- mrna[grepl("^Pp.+", mrna$Gene),]
mrna$mrnaPNtoPP <- mrna$PN/mrna$PP
mrna <- na.omit(mrna)  # Fix it

if (SWATH_SOURCE == "processed")
{
    prot <- read.table("proteins.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    prot <- data.frame(Gene = prot$ProteinName, protPNtoPP = prot$PNtoPP_1.fc, 
                       stringsAsFactors = FALSE)
} else if (SWATH_SOURCE == "raw")
{
    prot <- read.table("swath_raw.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    prot$Gene <- prot$protein_id
    prot$protein_id <- NULL
    prot$PN <- apply(prot[,PN_SAMPLES], 1, FUN = mean)
    prot$PP <- apply(prot[,PP_SAMPLES], 1, FUN = mean)
    prot$protPNtoPP <- prot$PN/prot$PP
    prot$pv <- apply(prot, 1, tst, PN_SAMPLES, PP_SAMPLES)
    if (PROT_FILTER_PV)
    {
        prot <- prot[prot$pv < PROT_REQ_PV,]
    }
} else
{
    stop("Wrong SWATH source!")
}

prot <- prot[grepl("^Pp.+", prot$Gene),]
prot$Gene <- sapply(strsplit(prot$Gene, split = '\\.'), "[", 1)
prot <- aggregate(. ~ Gene, data = prot, FUN = sum)

total <- merge(mrna[,c("Gene", "mrnaPNtoPP")], 
               prot[, c("Gene", "protPNtoPP")], by = "Gene" )#, all.y = TRUE)

total$mrnafc <- log(total$mrnaPNtoPP, base = 1.5)
total$protfc <- log(total$protPNtoPP, base = 1.5)

total$group <- as.factor(apply(total[,c("mrnafc", "protfc")], 1,
                               function(x) group(x[1], x[2])))

total <- total[ abs(total$mrnafc) > TH |
                abs(total$protfc) > TH,]

#plot(total$mrnafc, total$protfc, pch = 19, col = total$group)
p <- ggplot(total, aes(x=protfc, y=mrnafc, colour = group))
p <- p + geom_point(size  = 2)
print(p)
