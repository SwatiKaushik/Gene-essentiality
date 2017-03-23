setwd("/Users/swati/Desktop/cnv-enrichment-EGFRm/EGFR-NA-high-samples")
rm(list=ls())

file <- read.table("LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.all_thresholded.by_genes.txt", header=T, sep="\t")
#file <- read.table("test.dat", header=T, sep="\t")
file = file[,-c(2,3)]
t =t(file)
write.table(t, file="LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.all_thresholded.by_genes.txt.t", sep="\t", col.names =FALSE)
#write.table(t, file="test.t", sep="\t", col.names =FALSE)

cnv.file <- read.table("LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.all_thresholded.by_genes.txt.t", header=T, sep="\t")
#cnv.file <- read.table("test.t", header=T, sep="\t")
gene.sym <- cnv.file$Gene.Symbol
sub.str <- substr(gene.sym,1,15)
new.sample.name <- gsub(".",'-',sub.str, fixed=T)
cnv.file$Gene.Symbol <- new.sample.name

mutant.sample.file <- read.table("Final.input.heatmap-EGFR-highNA-nomutant.txt", header=T, sep="\t")
mutant.sample.file <- mutant.sample.file[c(1),]
mutant.sample.names <- gsub(".",'-',colnames(mutant.sample.file), fixed=T)
colnames(mutant.sample.file) <- mutant.sample.names
sampl <- t(mutant.sample.file)
write.table(sampl, file="Final.input.heatmap-EGFR-highNA-nomutant.txt.T", sep="\t",col.names =FALSE)
mutant.sample.file <- read.table("Final.input.heatmap-EGFR-highNA-nomutant.txt.T", header=T, sep="\t",)

merged.data <- merge(mutant.sample.file, cnv.file, by.x = 'A', by.y = 'Gene.Symbol')
write.table(merged.data, file="merged.file", sep="\t")

output <- vector()
#output <- data.frame(colnames=character(), pvalue=integer(), ampification=integer(), deletion=integer())

for (i in 3:ncol(merged.data)){
	
	amp <- subset(merged.data[,2], merged.data[,i] >1)
	del <- subset(merged.data[,2], merged.data[,i] < -1)
	amp.len <- length(amp)
	del.len <- length(del)
	cat (colnames(merged.data[i]),"\n")
	
	if ((amp.len >2) & (del.len >2)){

		p.value <- wilcox.test(amp,del)$p.value
		cat <- paste( colnames(merged.data[i]), p.value,amp.len, del.len,collapse='	')
		output <- rbind(output, cat)
		#output <- rbind(output, colnames(merged.data[i]), p.value,amp.len, del.len)
	} else { 
		#cat (colnames(merged.data[i]),"NA\n")
		no.data <- paste(colnames(merged.data[i]), NA, amp.len, del.len,collapse='	')
		output <- rbind(output, no.data)
		#output <- rbind(output, colnames(merged.data[i]), NA,amp.len, del.len)
	}
			
}

write.table(output, file="output", sep="\t", row.names= FALSE, col.names=FALSE, quote = FALSE)

output.file <- read.table("output", header=F, sep =" ")
sorted.output.file <- output.file[order(output.file$V2),]
sorted.output.file$p.adjust <- p.adjust(sorted.output.file$V2, method="fdr")
write.table(sorted.output.file, file = "output.padjust.txt" , sep="\t")








