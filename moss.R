require(ggplot2)
####### Parameters #######
FC2_TH = 1
FC15_TH = 1
TH = FC15_TH
SWATH_SOURCE = "raw"   # "raw" or "processed"
PN_SAMPLES = c("F163", "F164", "F165")
PP_SAMPLES = c("F166", "F167", "F168", "F169")
PROT_FILTER_PV = FALSE
PROT_REQ_PV = 2

##########################

group <- function(mrna, prot, pvmrna, protpv)
{
    if (abs(mrna) > TH)
        if (mrna > 0)
            m <- 'M'
        else
            m <- 'm'
    else
        m <- '0'
    
    if (abs(prot) > TH & protpv < PROT_REQ_PV)
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

id_table <- data.frame(Gene = mrna$Gene, TAIR = mrna$TAIR, stringsAsFactors = FALSE)
id_table <- id_table[grepl("^AT.+", id_table$TAIR), ]

mrna <- data.frame(Gene = mrna$Gene,# TAIR = mrna$TAIR,
                  PP = mrna$Protoplasts_FPKM,
                  PN = mrna$Protonema_FPKM,
                  stringsAsFactors = FALSE)

mrna <- mrna[grepl("^Pp.+", mrna$Gene),]
mrna$mrnaPPtoPN <- mrna$PP/mrna$PN
mrna <- na.omit(mrna)  # Fix it

if (SWATH_SOURCE == "processed")
{
    prot <- read.table("proteins.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    prot <- data.frame(Gene = prot$ProteinName, protPPtoPN = 1/prot$PNtoPP_1.fc, 
                       stringsAsFactors = FALSE)
    prot$protpv <- 0
} else if (SWATH_SOURCE == "raw")
{
    prot <- read.table("swath_raw.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    prot$Gene <- prot$protein_id
    prot$protein_id <- NULL
    prot$PN <- apply(prot[,PN_SAMPLES], 1, FUN = mean)
    prot$PP <- apply(prot[,PP_SAMPLES], 1, FUN = mean)
    prot$protPPtoPN <- prot$PP/prot$PN
    prot$protpv <- apply(prot, 1, tst, PN_SAMPLES, PP_SAMPLES)
    if (PROT_FILTER_PV)
    {
        prot <- prot[prot$protpv < PROT_REQ_PV,]
    }
} else
{
    stop("Wrong SWATH source!")
}

prot <- prot[grepl("^Pp.+", prot$Gene),]
prot$Gene <- sapply(strsplit(prot$Gene, split = '\\.'), "[", 1)
prot <- aggregate(. ~ Gene, data = prot, FUN = sum)

total <- merge(mrna[,c("Gene", "mrnaPPtoPN")], 
               prot[, c("Gene", "protPPtoPN", "protpv")], by = "Gene" )#, all.y = TRUE)

total$mrnafc <- log(total$mrnaPPtoPN, base = 1.5)
total$protfc <- log(total$protPPtoPN, base = 1.5)

total$group <- as.factor(apply(total[,c("mrnafc", "protfc", "protpv")], 1,
                               function(x) group(x[1], x[2], 0, x[3])))

total <- merge(total, id_table, by = "Gene", all.x = TRUE)
dir.create("groups", showWarnings = FALSE)
for (g in levels(total$group))
{
    write(total[total$group==g, "TAIR"], file = paste0("groups/", g, ".txt"))
}

#total <- total[ abs(total$mrnafc) > TH |
#                abs(total$protfc) > TH,]

#plot(total$mrnafc, total$protfc, pch = 19, col = total$group)
p <- ggplot(total, aes(x=protfc, y=mrnafc, colour = group))
p <- p + geom_point(size  = 3)
print(p)

#backgrpond

#IDA
ida <- read.table("ida_by_spectra.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ida$Gene <- sapply(strsplit(ida$Gene, split = '\\.'), "[", 1)
ida <- aggregate(. ~ Gene, data = ida, FUN = sum)
ida$idaPPtoPN <- ida$PP/ida$PN

# Pathways
url <- "http://pmn.plantcyc.org/PLANT/pathway-genes?object=GLYCOLYSIS&ENZORG=TAX-3218"
fl <- file(url, open = "r")
path_name <- readLines(fl, n=1)
enz <- read.table(fl, skip = 1, , sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rm(fl)
enz <- enz[enz$Organism == "Physcomitrella patens", 
           c("Gene.name", "Enzymatic.activity", "Reaction.EC")]
enz <- data.frame(Gene = enz$Gene.name,
                  Reaction = enz$Enzymatic.activity,
                  EC = enz$Reaction.EC,
                  stringsAsFactors = FALSE)
enz$Gene <- sapply(strsplit(enz$Gene, split = '\\.'), "[", 1)
enz <- unique(enz)

pathway <- merge(ida[,c("Gene", 'idaPPtoPN')], 
                 prot[,c("Gene", "protPPtoPN")],
                 by = "Gene", all = TRUE)
pathway <- merge(pathway, enz, by = "Gene", all.y = TRUE)
pathway <- pathway[(!is.na(pathway$idaPPtoPN)) | (!is.na(pathway$protPPtoPN)),]
pathway <- pathway[order(pathway$EC),]
write.table(pathway, "pathways/glycolysis.txt", sep = "\t")
