## CONCATENATE BOOTSTRAP TREES SHOWING SUPPORT BELOW 50% (NON-MAJORITY RULE)
#Jeronymo Dalapicolla, 2021



#libraries
library(phangorn)
library(ape)

#Load the 1,000 bootstrap for species trees:
tree1 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_1")
tree2 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_2")
tree3 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_3")
tree4 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_4")
tree5 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_5")
tree6 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_6")
tree7 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_7")
tree8 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_8")
tree9 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_9")
tree10 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_10")
tree11 = read.nexus("./boots_svd_speciestree/Proechimys_spboot_11")

all_trees = c(tree1, tree2, tree3, tree4, tree5, tree6, tree7, tree8, tree9, tree10, tree11)
class(all_trees)


#calculate the Maximum Clade Credibility tree from trees
max_clade_cred = maxCladeCred(all_trees, tree = TRUE, part = NULL, rooted = TRUE)
class(max_clade_cred)
plot(max_clade_cred, main="Maximum clade credibility tree")


#reroot tree
max_clade_cred$tip.label #check tips and choose the outgroup
reroot_tree = root(max_clade_cred, outgroup = c("tri"))
plot(reroot_tree, main="Maximum clade credibility tree")

## get proportions of each clade:
clad = prop.clades(reroot_tree, all_trees, rooted = TRUE)
## get proportions of each bipartition:
boot = prop.clades(reroot_tree, all_trees)


#visual
plot(reroot_tree, main = "Bipartition vs. Clade Support Values")
drawSupportOnEdges(boot/10, cex = 1.4, frame= "none")
nodelabels(clad)
legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
       pt.bg = c("green", "lightblue"), pt.cex = 2.5)

reroot_tree$boots = boot/10

#save trees
ape::write.tree(reroot_tree, file='MCC_speciestree_proechimys.txt')
ape::write.nexus(reroot_tree, file='MCC_speciestree_proechimys.nex')



###############################################
#Load the 1,000 bootstrap for individual trees:
tree1 = read.nexus("./boots_svd_individualtree/proechimys_MPE1_all_unsps_spboot_1")
tree2 = read.nexus("./boots_svd_individualtree/proechimys_MPE1_all_unsps_spboot_2")
tree3 = read.nexus("./boots_svd_individualtree/proechimys_MPE1_all_unsps_spboot_3")
tree4 = read.nexus("./boots_svd_individualtree/proechimys_MPE1_all_unsps_spboot_4")
tree5 = read.nexus("./boots_svd_individualtree/proechimys_MPE1_all_unsps_spboot_5")
tree6 = read.nexus("./boots_svd_individualtree/proechimys_MPE1_all_unsps_spboot_6")
tree7 = read.nexus("./boots_svd_individualtree/proechimys_MPE1_all_unsps_spboot_7")
tree8 = read.nexus("./boots_svd_individualtree/proechimys_MPE1_all_unsps_spboot_8")


all_trees = c(tree1, tree2, tree3, tree4, tree5, tree6, tree7, tree8)
class(all_trees)


#calculate the Maximum Clade Credibility tree from trees
max_clade_cred = maxCladeCred(all_trees, tree = TRUE, part = NULL, rooted = TRUE)
class(max_clade_cred)
plot(max_clade_cred, main="Maximum clade credibility tree")


#reroot tree
max_clade_cred$tip.label #check tips and choose the outgroup
reroot_tree = root(max_clade_cred, outgroup = c("tri"))
plot(reroot_tree, main="Maximum clade credibility tree")

## get proportions of each clade:
clad = prop.clades(reroot_tree, all_trees, rooted = TRUE)
## get proportions of each bipartition:
boot = prop.clades(reroot_tree, all_trees)


#visual
plot(reroot_tree, main = "Bipartition vs. Clade Support Values")
drawSupportOnEdges(boot/10, cex = 1.4, frame= "none")
nodelabels(clad)
legend("bottomleft", legend = c("Bipartitions", "Clades"), pch = 22,
       pt.bg = c("green", "lightblue"), pt.cex = 2.5)

reroot_tree$boots = boot/10

#save trees
ape::write.tree(reroot_tree, file='MCC_individualtree_proechimys.txt')
ape::write.nexus(reroot_tree, file='MCC_individualtree_proechimys.nex')


#END;
