rm(list=ls())
library(stringr)

data_raw <- read.csv("flower_isolates_sequences.txt",header=F,sep="\n")

n <- which(sapply(data_raw,FUN=function(x) grepl(pattern=">",x)))
data_split <- vector(mode="list")
for (i in 1:(length(n)-1)) {
  data_split[[i]] <- data_raw[n[i]:(n[i+1]-1),]
}

data_tab <- data.frame(species_id=sapply(data_split,function(x) str_extract(x[1],"(?<=\\>).*(?=\\.)")),
                       species_num=sapply(data_split,function(x) as.numeric(str_extract(x[1],"(?<=\\().*(?=\\))"))),
                       sequence=sapply(data_split,function(x) gsub("\\s.*$","",x[2])),
                       sequence_info=sapply(data_split,function(x) str_extract(x[2],"(?<=\\s).*")),
                       sequence_id=sapply(data_split,function(x) str_extract(x[3],"(?<=Sequence ID: ).*(?=Length)")),
                       sequence_length=sapply(data_split,function(x) as.numeric(str_extract(x[3],"(?<=Length: ).*(?=Number)"))),
                       sequence_number_of_matches=sapply(data_split,function(x) as.numeric(str_extract(x[3],"(?<=Matches: ).*"))))
data_tab$species_name <- sapply(data_tab$sequence_info,function(x) str_extract(x,".*(?=(\\s16S|\\schromosome))"))
data_tab$species_family <- gsub("([A-Za-z]+).*","\\1",data_tab$species_name)

data_tab$species_family <- gsub("Pseudomonas","Pseudomonadaceae",data_tab$species_family)
data_tab$species_family <- gsub("Erwinia","Erwiniaceae",data_tab$species_family)
data_tab$species_family <- gsub("Stenotrophomonas","Xanthomonadaceae",data_tab$species_family)
data_tab$species_family <- gsub("Pantoea","Enterobacteriaceae",data_tab$species_family)

data <- data_tab[,c("species_num","species_family","species_name","species_id","sequence","sequence_id","sequence_info","sequence_length","sequence_number_of_matches")]
data <- data[order(data$species_num),]
data <- data[order(data$species_family),]
row.names(data) <- NULL

# filter
data <- data[!is.na(data$species_name),]
unique_species_names <- unique(data$species_name[data$species_family!="Erwiniaceae"])
n <- rep(0,length(unique_species_names))
for (i in 1:length(unique_species_names)) {
  n[i] <- which(data$species_name==unique_species_names[i])[1]
}
data <- data[n,]
row.names(data) <- NULL

write.table(data,file="flower_isolates.csv",sep=",",quote=F,
            na="",row.names=F,col.names=T)

