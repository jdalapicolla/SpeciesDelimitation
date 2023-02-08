# Jeronymo Dalapicolla, 2021


#1. Splitting my phylip SNPs alignment in two databases
library(tidyverse)
library(tidyr)


#2. Read the phylip as matrix. Regular alignment functions take long time do to this.
matrix_alig = scan('proechimys_MPE1.usnps', what = 'character', sep = '\n')
information = matrix_alig[1]
information
#"184 34324"


#3. Remove the 1st line, with information, covert to dataframe, and rename the colunm
df = as.data.frame(matrix_alig[-1])
names(df) = "all_matrix"
names(df)


#4. Define the total number of SNPs.
snps_number = 34324


#5. Selecting the number of characters composing sample name, reducing the number of SNPs from total of characters: 
nchar(as.character(df[,1]))
nchar(as.character(df[,1])) - snps_number
#26 in all samples


#6. Splitting database in two columns: sample name (with the 26 first characters) and all SNPs 
df2 = separate(df, all_matrix, into = c("Sample", "Alignment"), sep = 26, remove = T)
#check
names(df2)
nchar(as.character(df2[,1])) # 26
nchar(as.character(df2[,2])) # 34324
df2[1:10,1] #10 fist names to verify if there is any SNP


#7. Splitting the Alignment in two columns: Subset 1 and Subset 2 with half of all SNPs 
snps_number/2 #half of SNPs
df3 = separate(df2, Alignment, into = c("Subset1", "Subset2"), sep = 17162, remove = T)
#check results
names(df3) 
nchar(as.character(df3[,2])) #17162
nchar(as.character(df3[,3])) #17162
nchar(as.character(df3[,2])) + nchar(as.character(df3[,3])) #34324 (original number of SNPs)


#8. Setting the first Subset of SNPs 
sub1 = df3 %>% unite("FirstSet", Sample,Subset1, sep = "", remove = T)
#check
names(sub1)


#9. Editing the subset1 to put it in phylip format
sub1 = sub1[,-2] #remove the second colunm
#check
names(sub1) #NULL
nchar(as.character(sub1)) #17188
nchar(as.character(sub1)) - 26 # 17162
#convert sub1 to data.frame
sub1 = as.data.frame(sub1)
header_1 = as.data.frame("184 17162") #define number of samples and SNPs
sub1 = rbind(header_1[1,], sub1)
#check
names(sub1)
sub1[1,1]
sub1[2,1]


#10. Save as phylip
write.table(sub1, "proechimys_sub1.phylip", sep = '\n', row.names = F, col.names = F, quote = F)



#11. Setting the second Subset of SNPs 
sub2 = df3 %>% unite("SecSet", Sample,Subset2, sep = "", remove = T)
#check
names(sub2)

#8. Editing the subset2 to put it in phylip format
sub2 = sub2[,-2] #remove the second colunm
#check
names(sub2) #NULL
nchar(sub2) #16408
nchar(sub2) - 26 # 16382
#convert sub1 to data.frame
sub2 = as.data.frame(sub2)
header_2 = as.data.frame (c("184 17162")) #define number of samples and SNPs
sub2 = rbind(header_2[1,], sub2)
#check
names(sub2)
sub2[1,1]
sub2[2,1]

#9. Save as phylip
write.table(sub2, "proechimys_sub2.phylip", sep = '\n', row.names = F, col.names = F, quote = F)


##END;
