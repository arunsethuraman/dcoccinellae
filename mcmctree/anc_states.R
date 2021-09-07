#First I had to follow the tutorial here on the mcmctree page for amino acids, out.BV to obtain a fossil-time calibrated tree
#Thinds I had to do: Tutorial Link: http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf, look for Tutorial 4 Approximate likelihood with protein data
#1. Take a multiple sequence alignment file (say from pasta/onlyhym*/*.aln), then paste that FASTA/upload it to http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi to convert it to the PHYLIP (sequential) format
#2. Then I pasted that into a trial3.aa file, and added spaces after each species name (>2-5). Also had to edit all the species names in the tree file (tree1.tre)
#3. Edited the ndata line in the mcmctree.ctl file to the number of loci, usedata to 3, then run src/mcmctree mcmctree.ctl
#4. This creates numerous files, one for each locus - TODO: Check for convergence
#5. Then I followed the tutorial, deleted the out.BV file, rst files, then copied the wag.dat file from the paml*/dat/ folder, created a new folder for the tmp files, pasted the wag.dat file there, then set the PATH to the ../paml*/bin/ folder for CODEML, then run codeml tmp001.ctl after editing the tmp001.ctl file accordingly - see page 13 of the tutorial
#6. Then renamed rst2 file as in.BV, and re-ran mcmctree.ctl, following the tutorial, and setting usedata = 2
#7. The corresponding output file now should have the species tree in NEWICK format (look for figtree output).

#Now that I have that, I can then use the phytools tutorial here: http://www.phytools.org/Cordoba2017/ex/8/Anc-states-discrete.html to do the ancestral state reconstruction. Have to create a character CSV or TAB file. Thereon:

require(phytools)
X<-read.table("characters.txt",row.names=1,as.is=FALSE)
reproduction<-setNames(X[,1],rownames(X))
hym.tree<-read.tree("speciestree.tre")

#fit Ancestral State Model - discrete
fitER<-ace(reproduction,hym.tree,model="ER",type="discrete")
cols<-setNames(c("red","blue","yellow"),levels(reproduction))
pdf("reproduction.pdf")
plotTree(hym.tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:hym.tree$Nnode+Ntip(hym.tree),pie=fitER$lik.anc,piecol=cols,cex=0.4)
tiplabels(pie=to.matrix(reproduction[hym.tree$tip.label],levels(reproduction)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)
dev.off()

#Empirical Bayes - using simulated datasets
mtrees<-make.simmap(hym.tree,reproduction,model="ER",nsim=100)
pd<-summary(mtrees)
pdf("reproduction_empbayes.pdf")
plot(mtree,cols,type="fan",fsize=0.7,ftype="i")
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],y=0.8*par()$usr[3],fsize=0.8)
dev.off()

#based on tutorial: http://www.phytools.org/Cordoba2017/ex/8/Anc-states-discrete.html

