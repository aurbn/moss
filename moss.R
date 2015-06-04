require(ggplot2)
require(KEGGREST)
require(plyr)
require(GOSemSim)
require(topGO)
require(cluster)
require(GO.db)
require(gplots)

####### Parameters #######
MRNA_BASE <- 2
PROT_BASE <- 1.4
MRNA_TH <- 1
PROT_TH <- 1
MRNA_REQ_PV <- 0.05
LOGNAME <- "plots/info.txt"

SWATH_SOURCE <- "raw"   # "raw" or "processed" or "own" or "own2" or "split"
NORM_SWATH <- FALSE #only for raw
IDA_SOURCE <- "empai" # "empai" or "count"
PN_SAMPLES <- c("F163", "F164", "F165")
PN_SAMPLES_ <- c("X1_PN", "X2_PN")
PP_SAMPLES <- c("F166", "F167", "F168", "F169")
PP_SAMPLES_ <- c("X1_PP", "X2_PP")
PROT_FILTER_PV <- FALSE
PROT_REQ_PV <- 0.05
EMPAI_PV_REQ <- 0.05

##### ANNOTATION PARAMETERS #####
ANN_METHOD <- "topGO" #  # "david" or "topGO"
PLASTID_GO <- "GO:0009536"
USE_PV_EMPAI <- TRUE
DAVID_REQ_PV <- 0.05
KEGG_SP <- "ppp"
PATHWAYS <- c("00010", "00020", "00030", "00040", "00190", "00195", "00196")
ALL_PATHWAYS <- scan("all_kegg.txt", what = character())
PATHWAYS <- ALL_PATHWAYS
USE_PV_EMPAI <- FALSE
SCALE_EMPAI <- FALSE
TAIR_ANN <- "tair" # uniprot or orgdb or tair
#http://www.kegg.jp/kegg-bin/search_pathway_text?map=ppp&keyword=&mode=1&viewImage=true

#USE_PV_EMPAI <- FALSE
##########################
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
unfold_gos <- function(x) trim(strsplit(x, ";")[[1]])

ggplotRegression <- function (fit, xname = NULL, yname = NULL, meth = "pearson") {
    ggplot(fit$model, aes_string(x = ifelse(is.null(xname), names(fit$model)[2], xname), 
                                 y = ifelse(is.null(yname), names(fit$model)[1], yname))) + 
        geom_point() +
        stat_smooth(method = "lm", col = "red") +
        ggtitle(paste(#"Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                           "r = ", signif(cor((fit$model)[2], (fit$model)[1], 
                                                  method = meth), 2),
                           ";\n Intercept =",signif(fit$coef[[1]],2 ),
                           "; Slope =",signif(fit$coef[[2]], 2)
                           #"; P =",signif(summary(fit)$coef[2,4], 2)
                           ))
}

group <- function(mrna, prot, mrnapv, protpv)
{
    if (abs(mrna) > MRNA_TH & mrnapv < MRNA_REQ_PV)
        if (mrna > 0)
            m <- 'M'
        else
            m <- 'm'
    else
        m <- '0'
    
    if (abs(prot) > PROT_TH & protpv < PROT_REQ_PV)
        if (prot > 0)
            p <- 'P'
        else
            p <- 'p'
    else
        p <- '0'
    
    paste0(m, p)
}

gdesk <- function(s,what = NULL)
{
    cs <- c(substr(s, 1, 1), substr(s, 2, 2))
    if (what == "m")
    {
        if ("m" %in% cs)
        {
            return("mrnaPP < mrnaPN")
        } else if ("M" %in% cs)
        {
            return("mrnaPP > mrnaPN")
        } else
        {
            return("mrnaPP = mrnaPN")
        }
    } else if (what == "p")
    {
        if ("p" %in% cs)
        {
            return("protPP < protPN")
        } else if ("P" %in% cs)
        {
            return("protPP > protPN")
        } else
        {
            return("protPP = protPN")
        }
    } 
    stop("Wrong what!")
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


kegg_get <- function(id, key)
{
    inf <- function(i)
    {
        kg <- keggGet(paste0(key, ":", i))
        ort <- kg[[1]]$ORTHOLOGY
        def <- kg[[1]]$DEFINITION
        c(ID = i, ORT = ifelse(is.null(ort), "NA", ort),
                  DEF =  ifelse(is.null(def), "NA", def))
    }
    
    t <- lapply(id, inf)
    t <- data.frame(matrix(unlist(t), nrow=length(t), byrow=T), stringsAsFactors = FALSE)
    names(t) <- c("ID", "Orthology", "Definition")
    t
}

gosim <- function(golist, ont = "BP")
{
    #golist <- golist[sapply(golist, length) != 0]
    m <- matrix(nrow = length(golist), ncol = length(golist))
    rownames(m) <- names(golist)
    colnames(m) <- names(golist)
    for (i in 1:length(golist))
        for  (j in 1:(length(golist)))
        {
            if(length(golist[[i]]) == 0 && length(golist[[j]]) == 0)
                t_ <- 1
            else if (length(golist[[i]]) == 0 || length(golist[[j]]) == 0)
                t_ <- 0
            else
                t_ <- mgoSim(golist[[i]], golist[[j]], ont = ont,
                                 organism = "arabidopsis", combine = "avg")
            
            m[i,j] <- t_
            m[j,i] <- t_
        }
        
    return(m)
}

dir.create("plots", showWarnings = FALSE)

########## READ MRNA #################
mrna_old <- read.csv("transcripts.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
mrna <- read.csv("gene_exp.diff", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
mrna <- mrna[ mrna$sample_1 %in% c("Protonema", "Protoplasts") &
              mrna$sample_2 %in% c("Protonema", "Protoplasts"),]
if (unique(mrna$sample_1) != "Protonema" |
    unique(mrna$sample_2) != "Protoplasts")
{
    stop("Input data error!")
}
mrna <- mrna[, c("gene", "value_1", "value_2", "q_value")]
mrna <- rename(mrna, c("gene" = "Gene", "value_1" = "PN", 
                       "value_2" = "PP", "q_value" = "mrnapv"))
mrna <- mrna[grepl("^Pp.+", mrna$Gene) ,]

id_table <- data.frame(Gene = mrna_old$Gene, TAIR = mrna_old$TAIR, stringsAsFactors = FALSE)
id_table <- id_table[grepl("^AT.+", id_table$TAIR), ]

mrna_old <- data.frame(Gene = mrna_old$Gene,# TAIR = mrna$TAIR,
                  PP = mrna_old$Protoplasts_FPKM,
                  PN = mrna_old$Protonema_FPKM,
                  stringsAsFactors = FALSE)

#mrna <- mrna[grepl("^Pp.+", mrna$Gene),]
#mrna <- aggregate(. ~ Gene, data = mrna, FUN = sum)
mrna <- mrna[!(duplicated(mrna$Gene) | duplicated(mrna$Gene,fromLast = TRUE)),]

mrna$mrnaPPtoPN <- mrna$PP/mrna$PN
mrna <- na.omit(mrna)  # Fix it
mrna$mrnafc <- log(mrna$mrnaPPtoPN, base = MRNA_BASE)

######### READ SWATH #############
if (SWATH_SOURCE == "processed")
{
    prot <- read.table("proteins.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    prot <- data.frame(Gene = prot$ProteinName, protPPtoPN = 1/prot$PNtoPP_1.fc, 
                       stringsAsFactors = FALSE)
    prot$protpv <- 0
} else if (SWATH_SOURCE == "raw")
{
    prot <- read.table("swath_raw.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    if (NORM_SWATH)
    {
        p <- prot[,-1]
        avgs <- apply(p, 2, mean)
        p <- sweep(p, 2, avgs, "/") 
        p$Gene <- prot$protein_id
        prot <- p
        rm(p)
    } else
    {
        prot$Gene <- prot$protein_id
    }
     
    prot$protein_id <- NULL
    prot$PN <- apply(prot[,PN_SAMPLES], 1, FUN = mean)
    prot$PP <- apply(prot[,PP_SAMPLES], 1, FUN = mean)
    prot$protPPtoPN <- prot$PP/prot$PN
    prot$protpv <- apply(prot, 1, tst, PN_SAMPLES, PP_SAMPLES)
    if (PROT_FILTER_PV)
    {
        prot <- prot[prot$protpv < PROT_REQ_PV,]
    }
} else if(SWATH_SOURCE == "own")
{
    prot <- read.table("protein_score_wide.txt", sep = '\t', header = TRUE, 
                       stringsAsFactors = FALSE)
    prot <- rename(prot ,c("protein_id" = "Gene"))
    
    #Correletions between SWATH bio.replicates
    PN_c <- lm(X1_PN ~ X2_PN, prot)
    PN_c1 <- lm(X1_PN ~ X2_PN, prot[prot$X1_PN < 1e7 & prot$X2_PN < 1e7,])
    p <- ggplotRegression(PN_c)
    print(p)
    ggsave("plots/prot_replicates_PN_all.png", p)
    p <- ggplotRegression(PN_c1)
    print(p)
    ggsave("plots/prot_replicates_PN_less1e7.png", p)
    
    PP_c <- lm(X1_PP ~ X2_PP, prot)
    PP_c1 <- lm(X1_PP ~ X2_PP, prot[prot$X1_PP < 1e7 & prot$X2_PP < 1e7,])
    p <- ggplotRegression(PP_c)
    print(p)
    ggsave("plots/prot_replicates_PP_all.png", p)
    p <- ggplotRegression(PP_c1)
    print(p)
    ggsave("plots/prot_replicates_PP_less1e7.png", p)
    
    prot$PN <- apply(prot[,PN_SAMPLES_], 1, FUN = mean)
    prot$PP <- apply(prot[,PP_SAMPLES_], 1, FUN = mean)
    prot$protPPtoPN <- prot$PP/prot$PN
    prot$protpv <- apply(prot, 1, tst, PN_SAMPLES_, PP_SAMPLES_)   

}else if(SWATH_SOURCE == "own2")
{
    prot <- read.table("protein_score_wide2.txt", sep = '\t', header = TRUE, 
                       stringsAsFactors = FALSE)
    prot <- rename(prot ,c("protein_id" = "Gene", "X1_PP" = "PP", "X1_PN" = "PN"))
    
    #Correletions between SWATH bio.replicates
    
    prot$protPPtoPN <- prot$PP/prot$PN
    prot$protpv <- apply(prot, 1, tst, "PN", "PP")   
    
} else
{
    stop("Wrong SWATH source!")
}

prot$tmp <- NA
prot$tmp <- gsub("[^\\|]+\\|([^\\|]+)\\|.+", "\\1", prot$Gene)
prot$Gene <- prot$tmp
prot$tmp <- NULL
uniprot2kegg <- read.csv("./uniprot2kegg.txt", sep = "\t", header = TRUE, 
                         stringsAsFactors = FALSE)
uniprot2kegg_ <- subset(uniprot2kegg, Entry %in% prot$Gene, c(Entry, KEGG))
prot <- merge(prot, uniprot2kegg_, by.x = "Gene", by.y = "Entry", all.x = TRUE)
prot$Gene <- ifelse(is.na(prot$KEGG), yes = prot$Gene, no = prot$KEGG)
prot$KEGG <- NULL
#prot <- prot[grepl("^Pp.+", prot$Gene),]
prot$Gene <- sapply(strsplit(prot$Gene, split = '\\.'), "[", 1)
prot <- aggregate(. ~ Gene, data = prot, FUN = sum)
prot$protfc <- log(prot$protPPtoPN, base = PROT_BASE)
write(paste(nrow(prot), "Quntified (SWATH) proteins"), LOGNAME, append=T)

########## READ COSMOSS ######
cos_names <- read.table("./cosmoss.genonaut.protein_name.txt", sep="\t",
                      header = TRUE,stringsAsFactors = FALSE)
cos_names$Gene <- sapply(strsplit(cos_names$accession, split = '\\.'), "[", 1) 
cos_names <- cos_names[!duplicated(cos_names$Gene), c("Gene", "value")]
cos_names <- rename(cos_names, c("value" = "cosmoss_name"))


cos_desc <- read.csv2("./cosmoss.genonaut.description.txt", sep="\t",
                        header = TRUE,stringsAsFactors = FALSE)
cos_desc$Gene <- sapply(strsplit(cos_desc$accession, split = '\\.'), "[", 1) 
cos_desc <- cos_desc[!duplicated(cos_desc$Gene), c("Gene", "value")]
cos_desc <- rename(cos_desc, c("value" = "cosmoss_description"))

########## READ TAIR ######
tairs <- read.csv("./pp_tair.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(tairs) <- c("Gene", "TAIR", "ident", "ev")
rownames(tairs) <- NULL
tairs <- subset(tairs, ident > 50 & ev < 0.05)
tairs <- tairs[order(tairs$Gene, -tairs$ident, tairs$ev),]
tairs <- tairs[!duplicated(tairs),]
tairs_one <- tairs[!duplicated(tairs$Gene),]
tairs_one <- subset(tairs_one, select = c("Gene", "TAIR"))

if (TAIR_ANN == "orgdb")
{
    require(org.At.tair.db)
    tair_names_ <- org.At.tairSYMBOL
    keys_ <- mappedkeys(tair_names_)
    symbols_ <- sapply(as.list(tair_names_[keys_]),
                       function(q) paste(unlist(q), collapse = ","))
    tair_names <- data.frame(TAIR = keys_, T_SYM = symbols_, stringsAsFactors = FALSE)
    
    
    tair_desc <- org.At.tairGENENAME
    tair_desc_keys <- mappedkeys(tair_desc)
    tair_desc <- sapply(as.list(tair_desc[tair_desc_keys]), function(q) q[[1]][1])
    tair_desc <- tair_desc[!duplicated(names(tair_desc))]
    tair_desc <- data.frame(TAIR = names(tair_desc), T_DESC = tair_desc, 
                            stringsAsFactors = FALSE)
    
    tairs_one <- merge(tairs_one, tair_names, by = "TAIR", all.x = TRUE)
    tairs_one <- merge(tairs_one, tair_desc, by = "TAIR", all.x = TRUE)
    tairs_one$T_GO <- NA
} else if (TAIR_ANN == "uniprot")
{
    tairs_up <- read.csv2("./TAIR_uniprot.txt", sep ="\t", header = TRUE,
                          stringsAsFactors = FALSE)
    tairs_up <- tairs_up[!duplicated(tairs_up$TAIR), c("TAIR", "T_SYM", 
                                                       "T_DESC", "T_GO")]
    tairs_one <- merge(tairs_one, tairs_up, by = "TAIR", all.x  = TRUE)
    
} else if (TAIR_ANN == "tair")
{
    if (file.exists("./tairs_one.txt"))
    {
        tairs_ta <- read.csv("./tairs_one.txt", sep = "\t", header = TRUE, 
                              stringsAsFactors = FALSE)
    } else {
        tairs_ta <- read.csv("./tair_gene_aliases.txt", sep ="\t", header = TRUE,
                             stringsAsFactors = FALSE)
        names(tairs_ta) <- c("TAIR", "T_SYM", "T_DESC")
        tmp1 <- aggregate(tairs_ta$T_SYM, by = list(tairs_ta$TAIR), 
                          FUN = function(x) paste(unique(x[!grepl("^\\s+$|^$", x)]), 
                                                  collapse = ", "))
        names(tmp1) <- c("TAIR", "T_SYM")
        tmp2 <- aggregate(tairs_ta$T_DESC, by = list(tairs_ta$TAIR), 
                          FUN = function(x) paste(unique(x[!grepl("^\\s+$|^$", x)]), 
                                                  collapse = ", "))
        names(tmp2) <- c("TAIR", "T_DESC")
        
        tairs_ta <- merge(tmp1, tmp2, by = "TAIR", all = TRUE)
        write.table(tairs_one, "./tairs_one.txt", sep = "\t", 
                    quote = FALSE, row.names = FALSE)
    }
    tairs_one <-  merge(tairs_one, tairs_ta, by = "TAIR", all.x  = TRUE)
} else {
    stop("Wrong tair annotation!")
}

################


total <- merge(mrna[,c("Gene", "mrnaPPtoPN", "mrnafc", "mrnapv", "PP", "PN")], 
               prot[, c("Gene", "protPPtoPN", "protfc","protpv", "PP", "PN")], by = "Gene" )#, all.y = TRUE)
total <- rename(total,c("PP.x" = "PPmrna", "PN.x" = "PNmrna", "PP.y" = "PPswath", "PN.y" = "PNswath"))
total <- unique(total)

prot_dep <- prot
prot_dep <- rename(prot_dep, c("PN" = "protPN", "PP" = "protPP"))
prot_dep <- merge(prot_dep, mrna[,c("Gene", "mrnaPPtoPN",
                                    "mrnafc", "mrnapv", "PP", "PN")],
                                by = "Gene", all.x = TRUE )
prot_dep <- rename(prot_dep, c("PN" = "mrnaPN", "PP" = "mrnaPP", "mrnafc" = "logmrnafc"))
#prot_dep <- merge(prot_dep, tairs_one, by = "Gene", all.x = TRUE)




totalpp <- total[,c("Gene", "PPmrna", "PPswath")]
totalpp <- rename(totalpp, c("PPmrna" = "mrna", "PPswath" = "swath"))
totalpp$cells <- "PP"
totalpn <- total[,c("Gene", "PNmrna", "PNswath")]
totalpn <- rename(totalpn, c("PNmrna" = "mrna", "PNswath" = "swath"))
totalpn$cells <- "PN"
totalp <- rbind(totalpp, totalpn)
tmpm <- lm(mrna ~ swath, data = totalp)
p <- ggplotRegression(tmpm)
p <- p + xlab(expression(SWATH,~units))
p <- p + ylab(expression(mRNA,~FPKM))
ggsave("plots/mRNA_SWATH_raw.png", p)
print(p)

total$group <- as.factor(apply(total[,c("mrnafc", "protfc", "mrnapv","protpv")], 1,
                               function(x) group(x[1], x[2], x[3], x[4])))

total <- merge(total, id_table, by = "Gene", all.x = TRUE)
total <- unique(total)
write(paste(nrow(total), "Quantified mRNA and SWATH genes"), LOGNAME, append=T)


dir.create("groups", showWarnings = FALSE)
tbl <- with(warpbreaks, table(substr(total$group, 1, 1), substr(total$group, 2, 2),
                              dnn = c("MRNA", "SWATH")))
colnames(tbl) <- c("PP=PN", "PP<PN", "PP>PN")
rownames(tbl) <- c("PP=PN", "PP<PN", "PP>PN")
write.ftable(ftable(tbl),file = "plots/table.txt", sep = '\t', quote = FALSE)

total_ <- total
total_$mRNA <- sapply(total_$group, gdesk, what="m")
total_$SWATH <- sapply(total_$group, gdesk, what="p")
write.table(total_, "plots/mrna_swath_all.txt", sep="\t", quote = FALSE, row.names = FALSE)

p <- ggplot(total, aes(x=protfc, y=mrnafc, colour = group))
p <- p + geom_point(size  = 3)
p <- p + geom_hline()
p <- p + geom_vline()
p <- p + xlab(expression(log [1.4]~(protein~fold~change)))
p <- p + ylab(expression(log [2]~(transcript~fold~change)))
ggsave("plots/groups.png", p)
print(p)


tmp <- total[, c("Gene","mrnafc", "protfc", "group")]
tmp <- tmp[is.finite(tmp$mrnafc)&is.finite(tmp$protfc),]
tmpm <- lm(mrnafc ~ protfc, data = tmp)
p <- ggplotRegression(tmpm)
p <- p + xlab(expression(log [1.4]~(protein~fold~change)))
p <- p + ylab(expression(log [2]~(transcript~fold~change)))
ggsave("plots/mrna_prot.png", p)
print(p)

##### FUNCTIONAL ANNOTATION #####
bk_genes <- scan("background.txt", what = character())

if (ANN_METHOD == "david")
{
    library("RDAVIDWebService")
    david <- DAVIDWebService$new(email='anatoly.urban@phystech.edu')
    BG <- addList(david, bk_genes, idType="TAIR_ID", listName="bkgrd", listType="Background")
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL"))
} else if(ANN_METHOD == "topGO")
{
    require(topGO)
#     gosdb = read.table("cosmoss_go_melt.txt", sep = "\t", col.names = c("Gene", "GO", "Descr"), 
#                        stringsAsFactors = FALSE)
#     gosdb <- gosdb[grepl("^Pp.+", gosdb$Gene),]
#     gosdb$Gene <- sapply(strsplit(gosdb$Gene, split = '\\.'), "[", 1)
#     
    geneID2GO <- readMappings(file = "cosmoss_go_cast.txt")
    GO2geneID <- inverseList(geneID2GO)
    geneNames <- names(geneID2GO)

    getSigInTerm <- function(term, godata)
    {
        sg = sigGenes(GOdata)
        tg = genesInTerm(GOdata, term)[[1]]
        paste(intersect(sg,tg), collapse = ", ")
    }
}

for (g in levels(total$group))
{
    print(paste0("Annotating ", g))
    
    if (ANN_METHOD == "david")
    {
        list_go <- na.omit(total[total$group==g, "TAIR"])
        n_k <- length((list_go))
        t_ <- kegg_get(list_go, "ath")
        names(t_)[1] <- "TAIR"
        t_ <- merge(t_, total[,c("TAIR", "Gene", "mrnafc", "protfc")], by = "TAIR")
        t_ <- t_[,c("TAIR", "Gene", "mrnafc", "protfc", "Orthology", "Definition")]
        t_ <- unique(t_)
        write.table(t_, sep = "\t", quote = FALSE, row.names = FALSE,
                    file = paste0("groups/", g, sprintf("_%02i", n_k),  ".txt"))
        rm(t_)
        
        FG <- addList(david, list_go, idType="TAIR_ID", listName=g, listType="Gene")
        setCurrentBackgroundPosition(david,
            grep("bkgrd", getBackgroundListNames(david)))
        setCurrentGeneListPosition(david,
                                     grep(g, getGeneListNames(david)))
        ann <- getFunctionalAnnotationChart(david, threshold=DAVID_REQ_PV)
        
        setCurrentBackgroundPosition(david,
                                     grep("bkgrd", getBackgroundListNames(david)))
        setCurrentGeneListPosition(david,
                                   grep(g, getGeneListNames(david)))
        clann <- getClusterReport(david, type = "Term" )#, threshold=DAVID_REQ_PV)
        
        for (cl in clann@cluster)
        {
            fname <-  paste0("groups/", g, ".clusters.david.txt")
            write(paste("Enrichment",cl[[1]]), file = fname, append = TRUE)
            cltab <- cl[[2]][,c("Category", "Term", "Count",
                                "PValue", "Bonferroni", "Benjamini", "FDR")]
            cltab <- cltab[cltab$PValue < DAVID_REQ_PV,]
            if (nrow(cltab) > 0)
            {
                suppressWarnings(
                    write.table(cltab, file = fname, append = TRUE)
                    )
            } else
            {
                write("No significant results!", file = fname, append = TRUE)
            }
            write("\n", file = fname, append = TRUE)
         }
        
        #ann <- ann[ann$PValue < DAVID_REQ_PV,]
        if (nrow(ann) > 0)
        {
            ann <- ann[order(ann$PValue), ]
            write.table(ann[,c("Category", "Term", "Count", "PValue", "Genes")],
                file = paste0("groups/", g, ".david.txt"), quote = FALSE, sep = '\t')
        }else
        {
            write("No significant results!", file = paste0("groups/", g, ".david.txt"))
        }
    }else if (ANN_METHOD == "topGO")
    {
        geneList = total[total$group == g, "Gene"]
        genes = geneList
        geneList <- factor(as.integer(geneNames %in% geneList))
        if (length(levels(geneList)) != 2)
            next
        
        names(geneList) <- geneNames
        
        for (ont in c("CC", "BP", "MF"))
        {
            GOdata <- new("topGOdata",ontology=ont, allGenes = geneList, 
                          annot = annFUN.gene2GO, gene2GO = geneID2GO)
            
            test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
            resultFisher <- getSigGroups(GOdata, test.stat)
            resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
            
            allRes <- GenTable(GOdata, classic = resultFis,
                               ranksOf = "classic", topNodes = 100)#, orderBy = "weight")
            allRes <- allRes[allRes$classic < 0.05,]
            allRes$Genes <- sapply(allRes$GO.ID, FUN = getSigInTerm, godata = GOdata)
                                
            write.table(allRes, paste0("groups/", g,"_", ont, "_topGO.txt"),sep="\t")
            
        }
    }
}

require(topGO)
geneID2GO <- readMappings(file = "cosmoss_go_cast.txt")
GO2geneID <- inverseList(geneID2GO)


geneID2GO_TAIR <- readMappings(file = "tairGO_cast.txt")
GO2geneID_TAIR <- inverseList(geneID2GO_TAIR)
geneID2GO_CC_TAIR <- annFUN.GO2genes(whichOnto = "CC", GO2genes = GO2geneID_TAIR)
geneID2GO_CC_TAIR <- inverseList(geneID2GO_CC_TAIR)
geneID2GO_BP_TAIR <- annFUN.GO2genes(whichOnto = "BP", GO2genes = GO2geneID_TAIR)
geneID2GO_BP_TAIR <- inverseList(geneID2GO_CC_TAIR)
#### WRITE TABLES #####    

prot_dep <- merge(prot_dep, cos_names, by = "Gene", all.x = TRUE)
prot_dep <- merge(prot_dep, cos_desc, by = "Gene", all.x = TRUE)
prot_dep <- merge(prot_dep, tairs_one, by = "Gene", all.x = TRUE)
prot_dep <- rename(prot_dep, c("protfc" = "logprotfc"))
geneID2GO_CC <- annFUN.GO2genes(whichOnto = "CC", GO2genes = GO2geneID)
geneID2GO_CC <- inverseList(geneID2GO_CC)
geneID2GO_BP <- annFUN.GO2genes(whichOnto = "BP", GO2genes = GO2geneID)
geneID2GO_BP <- inverseList(geneID2GO_BP)
PLASTID_GOs  <- c(GOCCOFFSPRING[[PLASTID_GO]], PLASTID_GO)
prot_dep$isPlastid <- sapply(prot_dep$Gene, function(x) 
    any(geneID2GO[[x]] %in% PLASTID_GOs))
prot_dep$GO_CC <- sapply(prot_dep$Gene, function(x) paste(geneID2GO_CC[[x]], collapse = ", "))
prot_dep$GO_BP <- sapply(prot_dep$Gene, function(x) paste(geneID2GO_BP[[x]], collapse = ", "))


prot_dep$T_GO <-sapply(prot_dep$TAIR, function(x) 
    paste(geneID2GO_CC_TAIR[[x]], collapse = ", "))
prot_dep$isPlastid_tair <- sapply(prot_dep$TAIR, function(x) 
    any(geneID2GO_TAIR[[x]] %in% PLASTID_GOs))
prot_dep$isPlastid <- prot_dep$isPlastid | prot_dep$isPlastid_tair
prot_dep$isPlastid <- ifelse(prot_dep$GO_CC == "" & prot_dep$T_GO == "",
                             NA, prot_dep$isPlastid)
prot_dep[grepl("^PhpapaCp", prot_dep$Gene), "isPlastid"] <- TRUE
write(paste(sum(prot_dep$isPlastid, na.rm = TRUE), "SWATH proteins is Plastid"),
      LOGNAME, append=T)
write(paste(sum(!prot_dep$isPlastid, na.rm = TRUE), "SWATH proteins is not Plastid"),
      LOGNAME, append=T)
#<- ifelse(is.na(prot_dep$T_GO), NA, prot_dep$isPlastid)

prot4cls <- prot_dep

prot_dep <- prot_dep[, c("Gene",  "protPP", "protPN", "protPPtoPN", "logprotfc", "protpv", 
                         "mrnaPP", "mrnaPN", "mrnaPPtoPN", "logmrnafc",  "mrnapv",
                         "cosmoss_name", "cosmoss_description",
                         "isPlastid", "GO_CC", "GO_BP",
                         "TAIR", "T_SYM", "T_DESC", "T_GO" )]
prot_dep <- rename(prot_dep, c("cosmoss_name" = "Cosmoss name", 
                               "cosmoss_description" = "Cosmoss description",
                               "TAIR" = "TAIR homolog", "T_SYM" = "Homolog symbol", 
                               "T_DESC" = "Homolog description", "T_GO"= "TAIR GO CC"))
write.table(prot_dep, "plots/swathALL.txt", sep="\t", quote = FALSE, row.names = FALSE)
prot_dep_filt <- subset(prot_dep, abs(logprotfc) > 1)
prot_dep_filt <- subset(prot_dep_filt, protpv < 0.05)
write.table(prot_dep_filt, "plots/swathDEP.txt", sep="\t", quote = FALSE, row.names = FALSE)
write(paste(nrow(prot_dep), "differentialy expressed (SWATH) proteins"), LOGNAME, append=T)

prot_dep_false <- subset(prot_dep, isPlastid == FALSE)
write.table(prot_dep_false, "plots/swathNoPLastid.txt", sep="\t", 
            quote = FALSE, row.names = FALSE)
######## GO CLUSTERING AND HEATMAP ####################

color <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", 
                       "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100)
hmap_sep = which(!duplicated(prot_dep_hmap$cluster))
#heatmap(cbind(prot_dep_hmap$protPPtoPN, prot_dep_hmap$mrnaPPtoPN))
#x11()
my_palette <- colorRampPalette(c("green", "black", "red"))(n = 1000)
hmap <- as.matrix(cbind(prot_dep_hmap$logprotfc, prot_dep_hmap$logmrnafc))
hmap[!is.finite(hmap)] <- 0
colnames(hmap) <- c("Protein", "mRNA")
png("plots/hmap.png")
x11(width = 50)
heatmap.2(hmap,
          cexRow = 1,
          #           labRow = c(rep("", 100), "1",
          #                      rep("", 120), "2",
          #                      rep("", 60),  "3",
          #                      rep("", 100), "4",
          #                      rep("", 100), "5"),
          #labRow = c(rep("", 100), "1",
          #           rep("", 180), "2",
          #           rep("", 140),  "3"),
          lwid=c(.3,.3),
          #labRow = rn,
          labCol = c("SWATH", "FPKM"),
          Rowv= TRUE, #as.dendrogram(hc),
          Colv=FALSE,
          #cexRow=1,
          cexCol=1,
          dendrogram="row",
          scale="column",
          #labRow = FALSE,
          #rowsep = hmap_sep,
          sepwidth=c(1, 1),
          trace="none",
          density.info="none",
          key=FALSE,
          symkey = TRUE,
          col=color,
          margins = c(1,15))
dev.off()

###### EMPAI #############
#IDA

ida <- read.table("ida_by_spectra.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ida$Gene <- sapply(strsplit(ida$Gene, split = '\\.'), "[", 1)
ida <- aggregate(. ~ Gene, data = ida, FUN = sum)
ida$idaPPtoPN <- ida$PP/ida$PN
ida$idafc <- log(ida$idaPPtoPN, base = PROT_BASE)

#EMPAI
empai <- read.table("emPAI.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
empai$Gene <- sapply(strsplit(empai$Gene, split = '\\.'), "[", 1)
empai <- aggregate(. ~ Gene, data = empai, FUN = sum)
empai <- merge(empai, tairs_one, by = "Gene", all.x = TRUE)
empai$PN <- apply(empai[,PN_SAMPLES], 1, FUN = mean)
empai$PP <- apply(empai[,PP_SAMPLES], 1, FUN = mean)
empai$empaiPPtoPN <- empai$PP/empai$PN
empai$empaifc <- log(empai$empaiPPtoPN, base = PROT_BASE)
empai$empaipv <- apply(empai, 1, tst, PN_SAMPLES, PP_SAMPLES)
empai$isPlastid_pp <- sapply(empai$Gene, function(x) 
    any(geneID2GO[[x]] %in% PLASTID_GOs))
empai$isPlastid_at <- sapply(empai$TAIR, function(x) 
    any(geneID2GO_TAIR[[x]] %in% PLASTID_GOs))
empai$isPlastid <- empai$isPlastid_pp | empai$isPlastid_at

empai_all <- empai

write(paste(nrow(empai), "EMPAI proteins"), LOGNAME, append=T)
write(paste(sum(empai$isPlastid), "EMPAI proteins is Plastid"), LOGNAME, append=T)
write(paste(sum(!empai$isPlastid), "EMPAI proteins is not Plastid"), LOGNAME, append=T)
empai <- empai[empai$empaipv < EMPAI_PV_REQ, ]   # CHECK IT !!!!



if(SCALE_EMPAI == TRUE)
{
    tmp <- prot[, c("Gene", "protPPtoPN", "protfc")]
    tmp <- merge(tmp, empai_all[, c("Gene", "empaifc", "empaiPPtoPN")], by = "Gene")
    tmp <- na.omit(tmp)
    tmp <- tmp[is.finite(tmp$empaifc) & is.finite(tmp$protfc),]
    tmpm <- lm(empaiPPtoPN ~ protPPtoPN, data = tmp)
    empaiInt = tmpm$coefficients[1]
    empaiSlope = tmpm$coefficients[2]
    print(ggplotRegression(tmpm))
    
    empai_all$empaiPPtoPN <- (empai_all$empaiPPtoPN - empaiInt)/empaiSlope
    empai$empaifc <- log(empai$empaiPPtoPN, base = PROT_BASE)
    empai$empaipv <- apply(empai, 1, tst, PN_SAMPLES, PP_SAMPLES)
    empai_all <- empai
    empai <- empai[empai$empaipv < EMPAI_PV_REQ, ]
    tmp <- prot[, c("Gene", "protPPtoPN", "protfc")]
    tmp <- merge(tmp, empai_all[, c("Gene", "empaifc", "empaiPPtoPN")], by = "Gene")
    
    tmpm <- lm(empaiPPtoPN ~ protPPtoPN, data = tmp)
    p <- ggplotRegression(tmpm)
    p <- p + xlab(expression(protein~SWATH~fold~change))
    p <- p + ylab(expression(protein~EMPAI~fold~change))
    ggsave("plots/empai_prot_abs.png", p)
    print(p)
}



otab <- empai_all[,c("Gene", "PP", "PN", "empaiPPtoPN", "empaifc", "empaipv")]
otab <- otab[!(otab$PP == 0 & otab$PN == 0),]
otab <- rename(otab, c("PP" = "PPempai", "PN" = "PNempai"))
otab <- merge(otab, total, by = "Gene", all.x = TRUE)
otab <- rename(otab, c("empaifc" = "logempaifc", "mrnafc" = "logmrnafc", "protfc" = "logprotfc"))
write.table(otab, "plots/empai_mrna_swath.txt", sep="\t", quote = FALSE, row.names = FALSE)


tmp <- mrna[, c("Gene", "mrnafc")]
tmp <- merge(tmp, empai_all[, c("Gene", "empaifc")], by = "Gene")
tmp <- tmp[is.finite(tmp$mrnafc)&is.finite(tmp$empaifc),]
tmpm <- lm(mrnafc ~ empaifc, data = tmp)
p <- ggplotRegression(tmpm)
p <- p + xlab(expression(log [1.4]~(protein~EMPAI~fold~change)))
p <- p + ylab(expression(log [2]~(transcript~fold~change)))
ggsave("plots/mrna_empai.png", p)
print(p)


tmp <- prot[, c("Gene", "protPPtoPN", "protfc")]
tmp <- merge(tmp, empai_all[, c("Gene", "empaifc", "empaiPPtoPN")], by = "Gene")
tmp <- na.omit(tmp)
tmp <- tmp[is.finite(tmp$empaifc) & is.finite(tmp$protfc),]
tmpm <- lm(empaifc ~ protfc, data = tmp)
p <- ggplotRegression(tmpm)
p <- p + xlab(expression(log [1.4]~(protein~SWATH~fold~change)))
p <- p + ylab(expression(log [1.4]~(protein~EMPAI~fold~change)))
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
data$Source <- ifelse(data$ida, "emPAI", "SWATH")
#require(KEGG.db)

library(pathview)

d <- data$fc# / (3*sd(data$fc))
names(d) <- data$pd_id
catched <- c()

dir.create("pathways", showWarnings = FALSE)
dir.create("pathways_tmp", showWarnings = FALSE)
for (p in PATHWAYS)
{
    pv.out <- pathview(gene.data = d, pathway.id = p, gene.idtype = "KEGG",
                       kegg.dir = "pathways_tmp", species = KEGG_SP, out.suffix = "kegg",)
    if(length(pv.out) != 2)
        next

    ok_ <- !is.na(pv.out$plot.data.gene$mol.data)
    fname = paste0("pathways/", sprintf("%02i_", sum(ok_)), KEGG_SP, p, ".kegg" )
    keggs <- unique(pv.out$plot.data.gene[ok_,]$kegg.names)
    
    mm = max(abs(d[names(d) %in% keggs]))
    
    pv.out <- pathview(gene.data = d, pathway.id = p, gene.idtype = "KEGG",
                       kegg.dir = "pathways_tmp", species = KEGG_SP, out.suffix = "kegg",
                       limit=list(gene=c(-mm, mm), cpd=1))
    

    if(length(keggs) > 0)
    {
        catched <- c(catched, keggs)
        
        t_ <- kegg_get(keggs, "ppp")
        names(t_)[1] <- "pd_id"
        t_ <- merge(t_, data[,c("Gene", "pd_id", "fc", "Source")], by = "pd_id")
        t_ <- t_[,c("Gene", "pd_id", "fc", "Source", "Orthology", "Definition")]
        t_ <- unique(t_)
        t_ <- rename(t_, c("fc" = "logfc"))
        write.table(t_, sep = "\t", quote = FALSE, row.names = FALSE,
                    file = paste0(fname, ".txt"))
        file.rename(from = paste0(KEGG_SP, p, ".kegg.png"),
                    to = paste0(fname, ".png"))
    } else
    {
        unlink(paste0(KEGG_SP, p, ".kegg.png"))
    }
    
    rm(ok_)
}
data$catched <- data$pd_id %in% catched
print(paste(as.character(sum(data$catched)), "proteins in pathways"))
unlink("pathways_tmp", recursive = TRUE)
