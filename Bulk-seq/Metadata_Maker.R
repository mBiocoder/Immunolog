library(stringr)
#Generate metdata file
fastqDir= file.path("./raw_fastq/")
fastq <- list.files(fastqDir, pattern = "*.fq.gz")
fastq

#subsetting
SampleID <- str_sub(fastq, 1)
#SampleID
Group = c(rep("R1",2), rep("R2",2), rep("R3",2), rep("R4",2), rep("R5",2), rep("R6",2), rep("R7",2), rep("R8",2), rep("R9",2))
#Group
Replicate = c(rep("Rep1",1), rep("Rep2",1), rep("Rep1",1), rep("Rep2",1), rep("Rep1",1), rep("Rep2",1), rep("Rep1",1), rep("Rep2",1), rep("Rep1",1), rep("Rep2",1), rep("Rep1",1), rep("Rep2",1), rep("Rep1",1), rep("Rep2",1), rep("Rep1",1), rep("Rep2",1), rep("Rep1",1), rep("Rep2",1))
#Replicate

#make meta data
metadata <- data.frame(SampleID = SampleID,
                       Group = Group,
                       Replicate = Replicate)
metadata

#write to file called metadata
write.table(metadata, "metadata.txt", append = FALSE, sep = "\t", row.names = FALSE)