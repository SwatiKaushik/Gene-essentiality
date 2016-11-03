# compare essentiality of given genes in specified tumor vs other tumors
# t-statistics to calculate pvalue
# Swati Kaushik Nov 1 2016

#! /usr/bin/env Rscript 

rm(list=ls())
args <- commandArgs(TRUE)

# print usage
print.usage <- function(){
	cat("Rscript Map-Essentiality.R [RNAi file] [Tumorname][outputFile]\n")
	cat("Compares gene essentiality across tumor types\n")
}

# verify that at least required parameters are passed 
if(length(args)<3) {
	cat(print.usage)
	stop("Incorrect or missing required input!")
}

# verify if input file name is present in directory
data.file <- args[[1]] 
if (! file.exists(data.file)) {
	cat("RNAi file ", data.file,"does not exist\n")
	q(save="no",status=1)
}

#open input file
data.file <- read.table(data.file, header=T)
df <- rownames(data.file) #get the rows of the matrix

#take tumor type and convert it to uppercase
lp.tumor.type <- args[2]
tumor.type <- toupper(lp.tumor.type)

# check if tumortype is present in the input file
tumor.type.length <- unique(grepl(tumor.type, row.names(data.file)))

#if not stop
if(length(tumor.type.length) < 2)
 {  
 	stop("Tumor type not found in the input file\n") 
 }	 

#split the name of the cell lines
data.splitted <- stringr::str_split_fixed(df,"_",2)
#Repaste the cell line names 
names <- paste(data.splitted[,2], data.splitted[,1], sep=".")
#assign rownames again
rownames(data.file) <- names
#reorder the matrix
data.file <- data.file[order(rownames(data.file)),]

#grep the data of a tumor
#data <- grep(args[2],rownames(data.file), perl=TRUE, value=TRUE)

#generate two groups - matched vs non matched
matched.tumor <- data.file[grep(tumor.type, rownames(data.file)), ]
#count if essentiality value is more than 0 ....can be made more strict
matched.count <- apply(matched.tumor, 2, function(x) length(which(x>0)))
matched.length <- nrow(matched.tumor)
#calculate % cell lines with essential genes
matched.percentage <- matched.count*100/matched.length

#write output for matched data
output <- args[3]
output <- cbind(matched.count)
output <- cbind(output, matched.percentage)

#unmatched tumor types
unmatched.tumor <- data.file[grep(tumor.type, rownames(data.file), invert = TRUE), ]
#count if essentiality value is more than 0 ....can be made more strict
unmatched.count <- apply(unmatched.tumor, 2, function(x) length(which(x>0)))
unmatched.length <- nrow(unmatched.tumor)
#calculate % cell lines with essential genes
unmatched.percentage <- unmatched.count*100/unmatched.length

#write output for unmatched data
output <- cbind(output, unmatched.count)
output <- cbind(output, unmatched.percentage)
p.out <- vector()
  
#plot boxplot and significance of the genes
plot <- paste(args[3],".pdf", sep="")
pdf(file= plot)
par(mfrow=c(2,2))

#Generate boxplot and calculate pvalue of matched and unmatched tumors
for (i in 1:length(matched.tumor)) {
	
     pvalue <- (t.test(matched.tumor[,i],unmatched.tumor[,i], alternative = "greater")$p.value)
     boxplot(matched.tumor[,i],unmatched.tumor[,i],main=names(matched.tumor[i]), type="l", names=c(tumor.type, "others") )
     text(x=1.5, y=1, labels=round(pvalue, digits=4))
     p.out <- cbind(p.out, pvalue)
     
 }

output <- cbind(output, t(p.out))
write.table(output, file="out", sep="\t")
dev.off()
