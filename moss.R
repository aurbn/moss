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

total <- merge(mrna[,c("Gene", "mrnaPNtoPP")], prot, by = "Gene")

plot(total$mrnaPNtoPP, total$protPNtoPP)
