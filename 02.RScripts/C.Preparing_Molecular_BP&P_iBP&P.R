####################### VALE INSTITUTE OF TECHNOLOGY ##########################

####################### CONVERT .LOCI TO BPP or iBPP ##########################


### Script prepared by Jeronymo Dalapicolla ###




#### PRE-ANALYSIS #### 

##1. INPUTS FOR THIS TUTORIAL ----
#A. A .LOCI FILE, A OUTPUT FROM ipyrad OR stacks PIPELINE.

#B. THE FILE ".CSV" WITH INFORMATION ON SAMPLES ID, ON WHICH CLADE/POP/LOCATION EACH INDIVIDUAL IS LOCATED IN THE TREE, AND A PROVISIONAL NAME INCLUDING ID AND CLADE, SPLIT BY UNDERSCORES(_)


##2. GOALS FOR THIS STEP ----
#A. EDIT .LOCI FILES FROM GENOMICS DATA TO USE IN BP&P ANALYSES;
#B. SELECT LOCI BY A MINIMAL NUMBER OF SEGREGATING SITES

#### OPTIONAL OPTIONS:
#C. SELECT LOCI THAT HAVE LESS MISSING DATA;
#D. SELECT A SPECIFIC NUMBER OF LOCI;
#F. REDUCE THE NUMBER OF INDIVIDUALS BY POPULATION TO MAXIMUM OF TWO INDIVIDUALS;



##3. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! ----
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H

#B. MOVE TO YOUR WORKING DIRECTORY ALL NECESSARY FILES.


##4. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT ----
rm(list=ls())


##5. INSTALL AND LOAD THE PACKAGES ----
#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('tidyr'))        install.packages("tidyr");             library('tidyr')



#### ANALYSIS ---- 

###1. Read the file as matrix. Regular alignment functions take long time do to this ----
matrix_alig = scan('proechimys_MPE1.loci', what = 'character', sep = '\n')
class(matrix_alig)



###2. Define the number of characters to names ----
matrix_alig[1]
strsplit(matrix_alig[1], "")
names_char = 26 #26 in my case



###3. Define break lines for loci and other objects necessary to other steps ----
break.lines = grep('//', matrix_alig)
index.lines = c(0, break.lines)
loci_number = length(break.lines)
count_loci = c(0:loci_number)
count_ind = c(1, break.lines)



###4. Create and save a dataframe with information about the loci ----
len_loci = as.data.frame(matrix(NA, loci_number, 7))
colnames(len_loci) = c("Loci", "BreakLine", "ReadLength", "SampleSize", "InitialLine", "FinalLine", "S")
loop = 0

for(i in 1:length(break.lines)){
  loop = loop+1
  num.1 = index.lines[i] + 1
  num.2 = index.lines[i + 1] - 1
  temp.text = matrix_alig[c(num.1:num.2)]
  df1 = as.data.frame(temp.text)
  df2 = separate(df1,temp.text, into = c("Sample", "Alignment"), sep = 26, remove = T)
  temp_segregating = matrix(NA, length(temp.text), nchar(as.character(df2[1,2])))
  for (iii in 1:length(temp.text)){temp_segregating[iii,] = unlist(strsplit(df2[iii,2], ''))}
  df3 = as.data.frame(temp_segregating)
  
  seg_sites = c()
  for (p in 1:length(df3)){
    clean_sites = df3[,p][!df3[,p] %in% grep(paste0(c("N", "-", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V"), collapse = "|"), df3[,p], value = T)]
    if (length(unique(clean_sites)) != 1) {seg_sites = c(seg_sites, 1)}
  }
  
  theta_s = length(seg_sites)
  
  
  len_loci[loop,] = c(count_loci[i], break.lines[i], nchar(as.character(matrix_alig[i])) - names_char, length(temp.text), num.1, num.2, theta_s)
}

#check and save dataframe
head(len_loci)
tail(len_loci)
write.csv(len_loci, "summary_all_loci.csv")




###5. Remove loci out of limits of segregating sites ----
#you can specify whatever you want
lower.limit = 5; #at least how many var sites in a alignment #change this according to your need
upper.limit = 135; #at most how many var sites are allowed in an alignment
selected.loci = which(len_loci$S > lower.limit & len_loci$S < upper.limit)
selected.loci = len_loci[selected.loci,]
write.csv(selected.loci, "summary_selected_loci_S.csv")




####6. Transform your loci in a list ----
loci_subset = length(selected.loci[,1])
index_loci = selected.loci

loci = list()
loop = 0
for (j in 1:loci_subset){
  loop = loop+1
  initial_line = index_loci$InitialLine[j]
  final_line = index_loci$FinalLine[j]
  loci[[loop]] = matrix_alig[c(initial_line:final_line)]
}

#verify first and last locus
loci[[1]]
loci[[loci_subset]]






###7. Load file with information on populations or clades, putting a code for this in front of IDs ----
index_replace = read.csv("replace_names.csv")
head(index_replace)
tail(index_replace)

#ID: original ID from .loci file
#TEMP: Temporary names including the code for pop or clades + _ + number of individual by population. For example, population A01 has 5 individuals, so we have "A01_01" "A01_02" "A01_03" "A01_04" "A01_05". This columns is important if you wanna remove individuals in populations with large sample size.
#FINAL: Final name with ID + "_" + Individual order in original alignment that generate .loci. You have to kkep only one underscore!




###8. Rename individuals by population or group to create the header by loci ----
loci_subset = length(selected.loci[,1])
renamed_loci = list()

for (k in 1:loci_subset){
  temp = loci[[k]]
  for (l in 1:length(index_replace$ORIGINAL)){
    temp = gsub(paste0('\\b', index_replace$ORIGINAL[l], '\\b'), index_replace$TEMP[l], temp)
  }
  renamed_loci[[k]] = temp
}

renamed_loci[[1]]
length(renamed_loci) #34039
renamed_loci[[length(renamed_loci)]]





###9. Reducing the number of individuals to 2 to save time in analyses ----
temp = renamed_loci
ind = c("_01", "_02")
loci_subset = length(temp)
reduced_loci = list()

for (k in 1:loci_subset){
  input = temp[[k]]
  target = input[input %in% grep(paste0(ind, collapse = "|"), input, value = T)]
  if(!identical(target, character(0))) {reduced_loci[[k]] = target}
}

reduced_loci = reduced_loci[!sapply(reduced_loci,is.null)] 
length(reduced_loci)
reduced_loci[[1]]
reduced_loci[[length(reduced_loci)]]






###10. Testing have many loci are common to all populations ----
#a. Define which POPs need to be present in all loci

##BPP
#analysis CLADES:
clades = c("G01_", "E02_", "D02_", "F01_", "E01_", "A14_", "A13_", "D01_", "A0", "B0", "C0")
#analysis AN1:
POPA1 = c("A01_", "A02_", "A03_", "A04_", "A05_", "A06_")
#analysis AN2:
POPA2 = c("A07_", "A08_", "A09_", "A10_", "A11_", "A12_")
#analysis BN1:
POPB1 = c("B01_", "B02_", "B03_", "B04_", "B05_", "B06_", "B07_", "B08_")
#analysis BN2:
POPB2 = c("B09_", "B10_", "B11_", "B12_", "B13_")
#analysis CN1:
POPC = c("C01_", "C02_", "C03_", "C04_", "C05_", "C06_", "C07_", "C08_")


##iBPP
#analysis CLADES:
clades = c("G01_", "E02_", "F01_", "E01_", "A14_", "A13_", "D01_", "A0", "B0", "C0")
#analysis AN2:
POPA2 = c("A07_", "A08_", "A09_", "A10_", "A12_")
#analysis BN1:
POPB1 = c("B01_", "B02_", "B03_", "B04_", "B05_", "B07_", "B08_")



#b. Filtering loci present in ind in all pops:
#temp = renamed_loci #if you will use all individuals
temp = reduced_loci
pops = POPB1
nloci = c()


for (j in pops){
  loci_subset = length(temp)
  filtering_loci = list()
  
  for (k in 1:loci_subset){
    input = temp[[k]]
    target = input[input %in% grep(paste0(j, collapse = "|"), input, value = T)]
    if(!identical(target, character(0))) {filtering_loci[[k]] = input}
  }
  temp = filtering_loci[!sapply(filtering_loci,is.null)] 
  nloci = c(nloci,length(temp))
  }

filtered_loci = temp
nloci
length(filtered_loci)
filtered_loci[[1]]
filtered_loci[[length(filtered_loci)]]




###10. Randomly choose a number of loci  ----
#without unique loci
dataset1 = filtered_loci[sample(1:length(filtered_loci), 170, replace = FALSE)]
dataset2 = filtered_loci[sample(1:length(filtered_loci), 340, replace = FALSE)]


#If you have enough loci, for example unique loci per datasets:
ds1 = sample(1:length(filtered_loci), 250, replace = FALSE)
allow_loci = setdiff(1:length(filtered_loci), ds1)
dataset1 = filtered_loci[ds1]

ds2 = sample(allow_loci, 500, replace = FALSE)
allow_loci2 = setdiff(allow_loci, ds2)
dataset2 = filtered_loci[ds2]

ds3 = sample(allow_loci2, 1000, replace = FALSE)
dataset3 = filtered_loci[ds3]







###11.  Remove indviduals in each loci and Rename individuals by order in the tree files and add a header by loci ----
#loci_list = list(dataset1, dataset2, dataset3)
loci_list = list(dataset1, dataset2)


#remove_ind = c("OUT_") #clades in BPP

#remove_ind = c("OUT_", "G01_",  "F01_", "E01_", "E02_", "D02_", "D01_", "C01_", "C02_", "C03_", "C04_", "C05_", "C06_", "C07_", "C08_", "B01_", "B02_", "B03_", "B04_", "B05_", "B06_", "B07_", "B08_", "B09_", "B10_", "B11_", "B12_", "B13_", "A07_", "A08_", "A09_", "A10_", "A11_", "A12_", "A14_", "A13_") #AN1 in BPP

#remove_ind = c("OUT_", "G01_",  "F01_", "E01_", "E02_", "D02_", "D01_", "C01_", "C02_", "C03_", "C04_", "C05_", "C06_", "C07_", "C08_", "B01_", "B02_", "B03_", "B04_", "B05_", "B06_", "B07_", "B08_", "B09_", "B10_", "B11_", "B12_", "B13_", "A01_", "A02_", "A03_", "A04_", "A05_", "A06_", "A14_", "A13_") #AN2 in BPP

#remove_ind = c("OUT_", "G01_",  "F01_", "E01_", "E02_", "D02_", "D01_", "C01_", "C02_", "C03_", "C04_", "C05_", "C06_", "C07_", "C08_", "B09_", "B10_", "B11_", "B12_", "B13_", "A01_", "A02_", "A03_", "A04_", "A05_", "A06_", "A14_", "A13_","A07_", "A08_", "A09_", "A10_", "A11_", "A12_", "A14_", "A13_") #BN1 in BPP

#remove_ind = c("OUT_", "G01_",  "F01_", "E01_", "E02_", "D02_", "D01_", "C01_", "C02_", "C03_", "C04_", "C05_", "C06_", "C07_", "C08_", "B01_", "B02_", "B03_", "B04_", "B05_", "B06_", "B07_", "B08_", "A01_", "A02_", "A03_", "A04_", "A05_", "A06_", "A14_", "A13_","A07_", "A08_", "A09_", "A10_", "A11_", "A12_", "A14_", "A13_") #BN2 in BPP

#remove_ind = c("OUT_", "G01_",  "F01_", "E01_", "E02_", "D02_", "D01_", "B01_", "B02_", "B03_", "B04_", "B05_", "B06_", "B07_", "B08_", "B09_", "B10_", "B11_", "B12_", "B13_", "A01_", "A02_", "A03_", "A04_", "A05_", "A06_", "A14_", "A13_","A07_", "A08_", "A09_", "A10_", "A11_", "A12_", "A14_", "A13_") #CN1 in BPP

#remove_ind = c("OUT_", "D02_", "A11_", "B06_") #clades in iBPP

#remove_ind = c("OUT_", "G01_",  "F01_", "E01_", "E02_", "D02_", "D01_", "C01_", "C02_", "C03_", "C04_", "C05_", "C06_", "C07_", "C08_", "B01_", "B02_", "B03_", "B04_", "B05_", "B06_", "B07_", "B08_", "B09_", "B10_", "B11_", "B12_", "B13_", "A01_", "A02_", "A03_", "A04_", "A05_", "A06_", "A11_", "A14_", "A13_") #AN2 in iBPP

remove_ind = c("OUT_", "G01_",  "F01_", "E01_", "E02_", "D02_", "D01_", "C01_", "C02_", "C03_", "C04_", "C05_", "C06_", "C07_", "C08_", "B06_", "B09_", "B10_", "B11_", "B12_", "B13_", "A01_", "A02_", "A03_", "A04_", "A05_", "A06_", "A14_", "A13_","A07_", "A08_", "A09_", "A10_", "A11_", "A12_", "A14_", "A13_") #BN1 in iBPP


for (ll in 1:length(loci_list)){
  ind_reduced = loci_list[[ll]]
  final_renamed = list()
  loci_subset = length(ind_reduced)
  
  
  for (k in 1:loci_subset){
    input = ind_reduced[[k]]
    input = input[!input %in% grep(paste0(remove_ind, collapse = "|"), input, value = T)]
    
    for (l in 1:length(index_replace$TEMP)){
      input = gsub(paste0('\\b', index_replace$TEMP[l], '\\b'), index_replace$ORIGINAL[l], input)
    }
    for (l in 1:length(index_replace$TEMP)){
      input = gsub(paste0('\\b', index_replace$ORIGINAL[l], '\\b'), index_replace$FINAL[l], input)
    }
    find.p = '[[:alnum:]_]+[[:space:]]+'
    sub.pa = ''
    temp.dna.seq = sub(find.p,sub.pa,input)
    dna.seq = strsplit(temp.dna.seq, '')
    loci.length = length(dna.seq[[1]])
    header = paste(length(input), loci.length , collapse = " ")
    final_renamed[[k]] = c(header, input)
  }
  
  assign(paste0("loci_renamed_",ll),final_renamed)
}


#verify first and last locus
loci_renamed_1[[1]]
loci_renamed_2[[1]]
loci_renamed_1[[250]]
loci_renamed_2[[500]]






###12. Merge all loci into a single files replaces "_" and save it ----
#renamed = list(loci_renamed_1, loci_renamed_2, loci_renamed_3)
renamed = list(loci_renamed_1, loci_renamed_2)

for (q in 1:length(renamed)){
  final_renamed = renamed[[q]]
  loci_subset = length(final_renamed)
  bpp_file = c()
  
  for (p in 1:loci_subset){
    add_file = as.character(final_renamed[[p]])
    without_first = add_file[-1]
    without_first = paste('', without_first, sep="       ")
    without_first = stringr::str_replace(without_first, "_", "^")
    add_file = c(add_file[1],without_first)
    bpp_file = c(bpp_file, add_file)
  }
  assign(paste0("bpp_file_",q),bpp_file)
}



###13. Save the BPP input ----
write.table(bpp_file_1, "BN1i_GENETICS_250L.loci", append = FALSE, sep = "\n", quote = F, row.names = F, col.names = F)
write.table(bpp_file_2, "BN1i_GENETICS_500L.loci", append = FALSE, sep = "\n", quote = F, row.names = F, col.names = F)
#write.table(bpp_file_3, "SPECIES_TREE_1000L.loci", append = FALSE, sep = "\n", quote = F, row.names = F, col.names = F)


#END;
