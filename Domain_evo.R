###Script for plotting a tree with associated domain maps
#note that you will need to fill in bits of code demarked with "<your code here>"

###Set the working directory (typically the directory where this script is stored)
#Add the full path to your desired working directory to the quotes
setwd("<the full path to the directory this scipt lives>")

###Install the R packages needed for our analysis

#The BiocManager package is needed in order to install some other packages
#This is an if-statement that asks whether the package is already installed and installs if not
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
library("BiocManager")

#Install and load treeio if it isn't already
if (!requireNamespace("treeio", quietly = TRUE)){
  BiocManager::install("treeio")

}
library(treeio)

#Install and load ggtree if it isn't already
if (!requireNamespace("treeio", quietly = TRUE)){
  BiocManager::install("ggtree")
}
library(ggtree)

#Install and load ggtree if it isn't already
if (!requireNamespace("msa", quietly = TRUE)){
  BiocManager::install("msa")
}
library(msa)

###Install and load several packages that are installed with the standard base install function
#Make a list of package names that we'll need
package_list<-c("ape", "ips", "Biostrings", "phytools", "seqinr", "dplyr", "ggplot2")

#Loop to check if package is installed and loaded. If not, install/load
#If you get a warning saying "there is no package called <XYZ>", run the loop again
for(k in 1:length(package_list)){
  
  if (!require(package_list[k], character.only = TRUE)) {
    install.packages(package_list[k], dependencies = TRUE)
    library(package_list[k], character.only=TRUE)
  }
}


### Perform multiple sequence alignment with the msa R-package
#First, take a look at the manual page for the msa
<R command to look at manual page for the msa() function>

#Save the path to the directory where sequnce files live
#Be sure to end this path with a "/"
seq_dir<-"<directory where the input alignments live>"

#List the files in this directory
seq_file_list<-list.files(path = seq_dir, pattern = "fasta")

#Randomly choose one of the alignments to work on
#random number
n<-sample(size = 1, x = 1:length(seq_file_list))

#Store the path to the unaligned sequence file
#paste0 concatenates two strings into one long string and stores this as an object
seqs_path<-paste0(seq_dir, seq_file_list[n])

#Read in seq file as AAstringSet object 
seqs<-readAAStringSet(filepath = <an object we created that points to the file we want to read>, format = "fasta")

#Perform alignment (use the Muscle algorithm)
#Note that there are many choices of algorithms for multiple sequence alignment. I typically use mafft (outside of R)
#aln<-msaMuscle(seqs, type = "protein", verbose = TRUE)
aln<-msa(<an object that points to the seqs we read in>, method = "Muscle", type = "protein", verbose = TRUE)

#Convert the alignment object to a AAbin object with a  different variable name
aln_bin<-as.AAbin(aln)

#Write a text file of the alignent (in AAbin forman)
ape::write.FASTA(<an object that points to the AAbin alignment>, file = paste0("output_files/aln_",seq_file_list[n]))

### Perform phylogenetic tree inference
#Note that there are many sophisticated statistical approaches to tree inference
#Two common approaches are maximum likelihood inference (e.g. raxml, IQ-TREE) and Bayesian inference (e.g. MrBayes)
#These are typically implemented outside of R. 
#For simplicity today we will use a simpler statistical approach called "neighbor joining"
#When running actual analyses, more sophisticated methods are recommended

#Create a distance matrix based on pairwise genetic distance between all seqs
dist_mat<-dist.aa(aln_bin)

#Use the neighbor joining method to infer tree
tree<-nj(dist_mat)

#Take a quick look at the tree
plot.phylo(tree)

#Check out the different facets of an R tree object
class(tree)

#The taxon labels
tree$tip.label

#A table of all the branches
tree$edge

#A list of all the branch lengths
tree$edge.length

#Just for fun, let's make a table of all the branches and their lengths
BL_df<-data.frame(Ancestor_node=tree$edge[,1],
           Descendant_node=tree$edge[,2],
           Branch_length=tree$edge.length)

#The nodes (internal and external) are stored as numbers now but we can also store the external (tips) by their seq name
#Remain external nodes (tip labels)
#Loop through all the tip labels and asign them
for(t in 1:length(tree$tip.label)){
  BL_df$Descendant_node[which(BL_df$Descendant_node==t)]<-tree$tip.label[t]
}

#Root the tree using the midpoint rooting method
#Note: usually you would have a pre-defined outgroup but in this case we'll make our best guess based on branch lengths
root_tree<-midpoint.root(tree)

#Write the rooted version of the tree to a text file (newick format)
write.tree(<an object that points to the rooted tree>, file = paste0("output_files/tree_",seq_file_list[n]))

###We now have the tree that we'll need for our figure! 

### Perform domain analysis (outside of R)
#Open a browser and go to https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi

#Edit the file in a text editor so that it looks like:
#Query   Hit_type        PSSM-ID From    To      E-Value Bitscore        Accession       Short_name      Incomplete      Superfamily
#A_ang_AANG005961        specific        438889  113     158     1.68696e-14     66.4916 cd22117 F-box_FBXL4      -      cl45894

#Remove the top lines
#Remove the "Q#1 >" from each line 
#Replace the spaces within column headers

#Read the tsv file into R
domain_df<-read.table(file = "output_files/<name of the domain file you created>", header = TRUE, sep = "\t")

#Clean up this dataframe a bit
names(domain_df)[1]<-"Newick_label"

###Add a column that gives the length of each sequence
#Read in seq file (in a different format)
seqs2<-seqinr::read.fasta(file = seqs_path, seqtype = "AA")

#Create a df of sequence lengths and join it to the domain data
domain_dat_full<-right_join(domain_df, data.frame(Newick_label=names(seqs2), Seq_ln=getLength(seqs2)), by = "Newick_label")

#Change the classes in the dataframe
domain_dat_full[,1]<-paste(domain_dat_full[,1])
domain_dat_full[,4]<-as.numeric(paste(domain_dat_full[,4]))
domain_dat_full[,5]<-as.numeric(paste(domain_dat_full[,5]))
domain_dat_full[,6]<-as.numeric(paste(domain_dat_full[,6]))
domain_dat_full[,7]<-as.numeric(paste(domain_dat_full[,7]))
domain_dat_full[,8]<-paste(domain_dat_full[,8])
domain_dat_full[,9]<-paste(domain_dat_full[,9])
domain_dat_full[,10]<-paste(domain_dat_full[,10])
domain_dat_full[,11]<-paste(domain_dat_full[,11])
domain_dat_full[,12]<-as.numeric(paste(domain_dat_full[,12]))
#Make a new column that's the same as newick labels
domain_dat_full[,13]<-paste(domain_dat_full[,1])
names(domain_dat_full)[13]<-"TipLabels"


### Begin creating the tree/domain plot using ggplot
#Make a ggtree object 
p1<-ggtree(root_tree, branch.length ='none', ladderize = TRUE)

#Add tip names in as a facet
p2<-facet_plot(p1, panel='tip_labels',
               data=domain_dat_full, geom=geom_text, 
               mapping=aes(x=0, label= TipLabels), size=3)

#add seq length line
p3<-facet_plot(p2, panel = "domains", data = domain_dat_full, geom= geom_segment, 
               mapping = aes(x=0, xend=Seq_ln, y=y, yend=y), size=0.5, color='black')

#Add domains
p4<-facet_plot(p3, panel = "domains", data = domain_dat_full, geom=geom_segment, 
               aes(x=From, xend=To, y=y, yend=y, col=Short_name), size=3) +
  theme(legend.position = "right")

#Plot the final plot
p4





