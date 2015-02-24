require(ggplot2)
####### Parameters #######
FC2_TH <- 1
FC15_TH <- 1
TH <- FC15_TH
SWATH_SOURCE <- "processed"   # "raw" or "processed"
IDA_SOURCE <- "empai" # "empai" or "count"
PN_SAMPLES <- c("F163", "F164", "F165")
PP_SAMPLES <- c("F166", "F167", "F168", "F169")
PROT_FILTER_PV <- FALSE
PROT_REQ_PV <- 2
EMPAI_PV_REQ <- 0.05

##### ANNOTATION PARAMETERS #####
ANN_METHOD <- "david"  # "david"
DAVID_REQ_PV <- 0.05
KEGG_SP <- "ppp"
PATHWAYS <- c("00010", "00020", "00030", "00040", "00190", "00195" )
#http://www.kegg.jp/kegg-bin/search_pathway_text?map=ppp&keyword=&mode=1&viewImage=true

USE_PV_EMPAI <- FALSE
##########################

ggplotRegression <- function (fit) {
    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
        geom_point() +
        stat_smooth(method = "lm", col = "red") +
        ggtitle(paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                           "; Intercept =",signif(fit$coef[[1]],2 ),
                           ";\n Slope =",signif(fit$coef[[2]], 2),
                           "; P =",signif(summary(fit)$coef[2,4], 2)))
}

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

dir.create("plots", showWarnings = FALSE)

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
mrna$mrnafc <- log(mrna$mrnaPPtoPN, base = 1.5)


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
prot$protfc <- log(prot$protPPtoPN, base = 1.5)

total <- merge(mrna[,c("Gene", "mrnaPPtoPN", "mrnafc")], 
               prot[, c("Gene", "protPPtoPN", "protfc","protpv")], by = "Gene" )#, all.y = TRUE)

total$group <- as.factor(apply(total[,c("mrnafc", "protfc", "protpv")], 1,
                               function(x) group(x[1], x[2], 0, x[3])))

total <- merge(total, id_table, by = "Gene", all.x = TRUE)

##### FUNCTIONAL ANNOTATION #####
bk_genes <- scan("background.txt", what = character())

if (ANN_METHOD == "david")
{
    library("RDAVIDWebService")
    david <- DAVIDWebService$new(email='anatoly.urban@phystech.edu')
    BG <- addList(david, bk_genes, idType="TAIR_ID", listName="bkgrd", listType="Background")
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL"))
}

dir.create("groups", showWarnings = FALSE)
for (g in levels(total$group))
{
    print(paste0("Annotating ", g))
    
    list_go = total[total$group==g, "TAIR"]
    write(list_go, file = paste0("groups/", g, ".txt"))
    
    if (ANN_METHOD == "david")
    {
        FG <- addList(david, list_go, idType="TAIR_ID", listName=g, listType="Gene")
        setCurrentBackgroundPosition(david,
            grep("bkgrd", getBackgroundListNames(david)))
        setCurrentGeneListPosition(david,
                                     grep(g, getGeneListNames(david)))
        ann <- getFunctionalAnnotationChart(david, threshold=DAVID_REQ_PV)
        #ann <- ann[ann$PValue < DAVID_REQ_PV,]
        if (nrow(ann) > 0)
        {
            ann <- ann[order(ann$PValue), ]
            write.table(ann[,c("Category", "Term", "Count", "PValue")],
                file = paste0("groups/", g, ".david.txt"), quote = FALSE, sep = '\t')
        }else
        {
            write("No significant results!", file = paste0("groups/", g, ".david.txt"))
        }
    }    
}

######################################
#total <- total[ abs(total$mrnafc) > TH |
#                abs(total$protfc) > TH,]

#plot(total$mrnafc, total$protfc, pch = 19, col = total$group)
p <- ggplot(total, aes(x=protfc, y=mrnafc, colour = group))
p <- p + geom_point(size  = 3)
p <- p + geom_hline()
p <- p + geom_vline()
ggsave("plots/groups.png", p)
print(p)


tmp <- total[, c("Gene","mrnafc", "protfc", "group")]
tmp <- na.omit(tmp)
tmpm <- lm(mrnafc ~ protfc, data = tmp)
p <- ggplotRegression(tmpm)
ggsave("plots/mrna_prot.png", p)
print(p)

#backgrpond

#IDA

ida <- read.table("ida_by_spectra.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ida$Gene <- sapply(strsplit(ida$Gene, split = '\\.'), "[", 1)
ida <- aggregate(. ~ Gene, data = ida, FUN = sum)
ida$idaPPtoPN <- ida$PP/ida$PN
ida$idafc <- log(ida$idaPPtoPN, base = 1.5)

empai <- read.table("emPAI.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
empai$Gene <- sapply(strsplit(empai$Gene, split = '\\.'), "[", 1)
empai <- aggregate(. ~ Gene, data = empai, FUN = sum)
empai$PN <- apply(empai[,PN_SAMPLES], 1, FUN = mean)
empai$PP <- apply(empai[,PP_SAMPLES], 1, FUN = mean)
empai$empaiPPtoPN <- empai$PP/empai$PN
empai$empaifc <- log(empai$empaiPPtoPN, base = 1.5)
empai$empaipv <- apply(empai, 1, tst, PN_SAMPLES, PP_SAMPLES)
empai_all <- empai
empai <- empai[empai$empaipv < EMPAI_PV_REQ, ]   # CHECK IT !!!!


tmp <- prot[, c("Gene", "protfc")]
tmp <- merge(tmp, ida[, c("Gene", "idafc")], by = "Gene")
tmp <- na.omit(tmp)
tmpm <- lm(idafc ~ protfc, data = tmp)
p <- ggplotRegression(tmpm)
ggsave("plots/ida_prot.png", p)
print(p)

tmp <- prot[, c("Gene", "protfc")]
tmp <- merge(tmp, empai_all[, c("Gene", "empaifc")], by = "Gene")
tmp <- na.omit(tmp)
tmp <- tmp[is.finite(tmp$empaifc) & is.finite(tmp$protfc),]
tmpm <- lm(empaifc ~ protfc, data = tmp)
p <- ggplotRegression(tmpm)
ggsave("plots/empai_prot.png", p)
print(p)

# Pathways
#url <- "http://pmn.plantcyc.org/PLANT/pathway-genes?object=GLYCOLYSIS&ENZORG=TAX-3218"
#fl <- file(url, open = "r")
#path_name <- readLines(fl, n=1)
#enz <- read.table(fl, skip = 1, , sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#rm(fl)
### convert PP to PHYPADRAFT
pp2phypa <- read.table("cosmoss2phypa.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
pp2phypa$Gene <- sapply(strsplit(pp2phypa$Gene, split = '\\.'), "[", 1)

if(USE_PV_EMPAI)
{
    data <- merge(empai[, c("Gene", "empaifc")], prot[,c("Gene", "protfc")],
                  by ="Gene", all = TRUE)
}else
{
    data <- merge(empai_all[, c("Gene", "empaifc")], prot[,c("Gene", "protfc")],
                  by ="Gene", all = TRUE)
}

data <- merge(data, pp2phypa, by = "Gene" )#, all.x = TRUE)
tmp <- is.finite(data$protfc)
data[!tmp, "protfc"] <- data[!tmp, "empaifc"]
data$empaifc <- NULL
data$ida <- !tmp
data$fc <- data$protfc
data$protfc <- NULL
data <- data[is.finite(data$fc),]
#require(KEGG.db)

library(pathview)

d <- data$fc
names(d) <- data$pd_id

dir.create("pathways", showWarnings = FALSE)
dir.create("pathways_tmp", showWarnings = FALSE)
for (p in PATHWAYS)
{
    pv.out <- pathview(gene.data = d, pathway.id = p, gene.idtype = "KEGG",
                       kegg.dir = "pathways_tmp", species = KEGG_SP, out.suffix = "kegg")
    fname = paste0(KEGG_SP, p, ".kegg.png" )
    file.rename(from = fname, to = paste0("pathways/", fname))
}
unlink("pathways_tmp", recursive = TRUE)



# 
# enz <- enz[enz$Organism == "Physcomitrella patens", 
#            c("Gene.name", "Enzymatic.activity", "Reaction.EC")]
# enz <- data.frame(Gene = enz$Gene.name,
#                   Reaction = enz$Enzymatic.activity,
#                   EC = enz$Reaction.EC,
#                   stringsAsFactors = FALSE)
# enz$Gene <- sapply(strsplit(enz$Gene, split = '\\.'), "[", 1)
# enz <- unique(enz)
# 
# pathway <- merge(ida[,c("Gene", 'idaPPtoPN')], 
#                  prot[,c("Gene", "protPPtoPN")],
#                  by = "Gene", all = TRUE)
# pathway <- merge(pathway, enz, by = "Gene", all.y = TRUE)
# pathway <- pathway[(!is.na(pathway$idaPPtoPN)) | (!is.na(pathway$protPPtoPN)),]
# pathway <- pathway[order(pathway$EC),]
# write.table(pathway, "pathways/glycolysis.txt", sep = "\t")
