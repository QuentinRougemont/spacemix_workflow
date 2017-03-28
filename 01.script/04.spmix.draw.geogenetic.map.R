#!/usr/bin/env Rscript

#Author= "Quentin Rougemont"
#purpose = "analyse spacemix results" (see https://github.com/gbradburd/SpaceMix)
#how to run = "03.run.spacemix.R sample.covariance.matrix sample.size x_y.coordinates fast.model.opts long.model.opts nloci 
#last uptade = "06.09.2016"

argv <- commandArgs(TRUE)

spmix.data<-argv[1] 
mcn.freq.list<-argv[2]
mcmc.out<-argv[3]
geo.coord<-argv[4]
ind.name<-argv[5]

library(SpaceMix)
load(spmix.data)
load(mcn.freq)
load(mcmc.out)

#load("salmon_space_MCMC_output1.Robj")
#load("salmon_spacemix.data.Robj")

pop.coord <-unique(as.matrix(read.table(argv4)[,c(2:3)]))
sample.names<-as.character(as.matrix(read.table(argv4)[,2]))
#pop.coord <-(as.matrix(read.table("pop.coord")[,c(2:3)]))
#sample.names<-as.character(as.matrix(read.table("pop.coord")[,1]))

sample.colors <- rainbow(n=length(sample.names),start=4/6,end=6/6)[as.numeric(cut(pop.coord[,1],length(sample.names)))]

# Contstructing the geogenetic Map #########################################
salmon.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "salmon_space_MCMC_output1.Robj",
                                geographic.locations = pop.coord,
                                name.vector = sample.names,
                                color.vector = sample.colors,
                                quantile=0.95,
                                burnin=0)      

# Now we generate a map of the output showing sample names  at the locations of the maximum a posteriori (MAP) geogenetic location parameter estimates
pdf(file="salmon.geogenetic.map.pdf")
make.spacemix.map(spacemix.map.list = salmon.spacemix.map.list,
                text=TRUE,
                ellipses=FALSE,
                source.option=FALSE)
dev.off()

pdf(file="salmon.geogenetic.map.arrow.pdf")
make.spacemix.map(spacemix.map.list=salmon.spacemix.map.list,
                  text=TRUE,
                  source.option=FALSE,
                  ellipse=FALSE, ) #xlim=c(-23,40),ylim=c(42,71)
plot.admix.arrows(salmon.spacemix.map.list$geographic.locations,salmon.spacemix.map.list$MAPP.geogen.coords,admix.proportions) #, colors=salmon.spacemix.map.list$admix.source.color.vector,length=0.4) 
dev.off() 

admix.cols <- fade.admixture.source.points(rep("red",length(sample.names)),seq(1,0,length.out=length(sample.names)))
pdf(file="salmon.geogenetic.map.arrow.2.pdf")
make.spacemix.map(spacemix.map.list=salmon.spacemix.map.list,
                  text=TRUE,
                  source.option=FALSE,
                  ellipse=FALSE , xlim=c(-90,28),ylim=c(41,71))
plot.admix.arrows(salmon.spacemix.map.list$geographic.locations,salmon.spacemix.map.list$MAPP.geogen.coords, admix.proportions,length=0.2, colors=admix.cols) 
dev.off()                             

# Now, to visualize uncertainty in location parameter estimates,  we generate a map of the output showing 95% credible
#   ellipses for the geogenetic locations of all samples  and plotting sample names at the locations of the 
#   maximum a posteriori (MAP) geogenetic location parameter estimates
pdf(file="salmon.geogenetic.95ci.pdf",14,14)
make.spacemix.map(spacemix.map.list = salmon.spacemix.map.list,
                text=TRUE,
                ellipses=TRUE,
                source.option=FALSE)       
dev.off()               

# Now, to visualize the sources of admixture, we can plot  those as well
pdf(file="salmon95ci.source.admix.full.pdf", 20,12)
make.spacemix.map(salmon.spacemix.map.list, source.option=TRUE, text=TRUE) #,xlim=c(0.5,12),ylim=c(0,10))  
dev.off()
pdf(file="salmon95ci.source.admix.full.large.pdf", 20,12)
make.spacemix.map(salmon.spacemix.map.list, source.option=TRUE, text=TRUE ,xlim=c(-100,30),ylim=c(35,71))  
dev.off()

# This map looks like a bit of a mess, because even though most samples are drawing negligible amounts of admixture,
#   they're all drawing SOME admixture, so they all get plotted and the output is difficult to visually interpret.
# To do a better job, we can be selective about which admixutre  sources we highlight using the `query.spacemix.map` function
pdf(file="salmon.source_origin_pop2.pdf",14,14)
make.spacemix.map(salmon.spacemix.map.list,
                source.option=FALSE,
                text=FALSE, xlim=c(-23,43),ylim=c(42,72) )
#highlight tgeogenetic location for some sample  as well as the location of its source of admixture, which is plotted in italics with a dashed border around its 95% credible ellipse.                              
query.spacemix.map(focal.pops=c("Narcea","Loire","Emtsa","Olfusa","Kunda","Vindelaiven","Sela","Blackwater","Tuloma"),
                        spacemix.map.list = salmon.spacemix.map.list,
                        ellipses=F,source.option=TRUE)
dev.off()

a=c(10, 20, 30, 40 )
b=c(1,  11, 21, 31)

for(i in a) {
     for (j in b) { 
     pdf(file=paste("salmon.source",j,i,".pdf"),14,14)
        make.spacemix.map(salmon.spacemix.map.list, source.option=F,text=F,xlim=c(-125,30),ylim=c(35,122) )
        query.spacemix.map(focal.pops=c(sample.names[j:i]) , spacemix.map.list=salmon.spacemix.map.list, ellipses=F,source.option=TRUE)
dev.off()
     }
}

a <-(salmon.spacemix.map.list$MCMC.output$admix.proportions)
b <-t(apply(a,1, function(a) quantile(a, c(0.025,0.975))))
m <-(apply(a,1, mean))
x <-1:77

#faire un cbind de b,m,x,sample.names puis sort by regional areas
z<-cbind(as.matrix(sample.names),m,b)
write.table(z,"conf.intervals_mean",quote=F)
require(plotrix)
z<-read.table("confidence.reshape",T)

pdf(file="test.pdf",14,6)

par(mar=c(7.1,4.1,4.1,2.1))
plotCI(x,z[,2], ui=z[,4], li=z[,3], xlab="populations", ylab="admixture proportion", col=sample.colors, xaxt="n")
axis(1,1:77, z[,1], las=3, col="black" )

dev.off()
































pdf(file="salmon.source_origin_pop0.pdf",14,14)
make.spacemix.map(salmon.spacemix.map.list,
                source.option=FALSE,
                text=FALSE , xlim=c(-90,28),ylim=c(42,79) )
query.spacemix.map(focal.pops=c("gak","stewiacke","stpaul","stgenevieve","salmonier","auxfeuilles","emtsa","olfusa","kunda","loire","narcea","vindelaiven","narraguagus","cross","conne"), #, #,"emtsa","olfusa","kunda","chaloupe","stpaul","dugouffre","gaspereau","stewiacke"),
                        spacemix.map.list = salmon.spacemix.map.list,
                        ellipses=F,source.option=TRUE)
dev.off()























pdf(file="salmon.source_origin_pop0.pdf",14,14)
make.spacemix.map(salmon.spacemix.map.list,
                source.option=FALSE,
                text=FALSE , xlim=c(-90,28),ylim=c(42,79) )
query.spacemix.map(focal.pops=c("gak","stewiacke","stpaul","stgenevieve","salmonier","auxfeuilles","dugouffre","emtsa","olfusa","kunda","loire","narcea","vindelaiven","vieuxfort","narraguagus","auxsaumons","hurons"), #, #,"emtsa","olfusa","kunda","chaloupe","stpaul","dugouffre","gaspereau","stewiacke"),
                        spacemix.map.list = salmon.spacemix.map.list,
                        ellipses=F,source.option=TRUE)
dev.off()

pdf(file="salmon.source_origin_pop1.pdf",14,14)
make.spacemix.map(salmon.spacemix.map.list,
                source.option=FALSE,
                text=FALSE , xlim=c(-90,28),ylim=c(42,79) )
query.spacemix.map(focal.pops=c("gak","stewiacke","stpaul","stgenevieve","salmonier","auxfeuilles","dugouffre","emtsa","olfusa","kunda","loire","narcea","vindelaiven","vieuxfort","narraguagus","auxsaumons","hurons","mecatina","miramachi"), #,"emtsa","olfusa","kunda","chaloupe","stpaul","dugouffre","gaspereau","stewiacke"),
                        spacemix.map.list = salmon.spacemix.map.list,
                        ellipses=F,source.option=TRUE)
dev.off()

pdf(file="salmon.source_origin_pop2.pdf",14,14)
make.spacemix.map(salmon.spacemix.map.list,
                source.option=FALSE,
                text=FALSE , xlim=c(-90,28),ylim=c(42,79) )
query.spacemix.map(focal.pops=c("gak","stewiacke","stpaul","stgenevieve","salmonier","auxfeuilles","dugouffre","emtsa","olfusa","kunda","loire","narcea","vindelaiven","vieuxfort","narraguagus","auxsaumons","hurons","mecatina","miramachi","natashqua","matapedia","stmary"), #,"emtsa","olfusa","kunda","chaloupe","stpaul","dugouffre","gaspereau","stewiacke"),
                        spacemix.map.list = salmon.spacemix.map.list,
                        ellipses=F,source.option=TRUE)
dev.off()

pdf(file="salmon.source_origin_pop.pdf",14,14)
make.spacemix.map(salmon.spacemix.map.list,
                source.option=FALSE,
                text=FALSE , xlim=c(-90,28),ylim=c(42,79) )
query.spacemix.map(focal.pops=c("gak","stewiacke","stpaul","stgenevieve","salmonier","auxfeuilles","dugouffre","emtsa","olfusa","kunda","loire","narcea","vindelaiven","vieuxfort","narraguagus","auxsaumons","hurons","mecatina","miramachi","natashqua","matapedia","stmary","cross","rockyricver","malbaie","matane","stanne","jacquet"), #,"emtsa","olfusa","kunda","chaloupe","stpaul","dugouffre","gaspereau","stewiacke"),
                        spacemix.map.list = salmon.spacemix.map.list,
                        ellipses=F,source.option=TRUE)
dev.off()







# Procrustes Transformation 

