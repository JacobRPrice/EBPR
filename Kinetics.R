#----------------------------------------------------------
# Notes or References:
#
#----------------------------------------------------------
########
# Preliminary Items
########

# load utilities and functions
source("/PATH/TO/DIR/r-EBPR/R/utilities.R")
source("/PATH/TO/DIR/r-EBPR/R/functions.R")

# prep environment and get data
analysis_prep(c("phyloseq","ggplot2","gridExtra", "magrittr"))

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# 
########
cycleobj <- read.csv(file.path(data_path, "Cycle-Phos-SV30-SS-SVI.txt"), header=TRUE, sep="\t")
names(cycleobj)
str(cycleobj)
cycleobj$Cycle <- as.numeric(cycleobj$Cycle)
rownames(cycleobj) <- cycleobj$Cycle
cycleobj$Sample.o<-as.character(cycleobj$Sample.o)
cycleobj$Sample.s<-as.character(cycleobj$Sample.s)
cycleobj$Removal.PO4.o <- 100*((cycleobj$Init.PO4.o - cycleobj$Fin.PO4.o)/cycleobj$Init.PO4.o)
cycleobj$Removal.PO4.s <- 100*((cycleobj$Init.PO4.s - cycleobj$Fin.PO4.s)/cycleobj$Init.PO4.s)
cycleobj$SV30.o <- 100*cycleobj$SV30.o
cycleobj$SV30.s <- as.numeric(cycleobj$SV30.s)
summary(cycleobj)

samp <- cycleobj[,c(1,2,8)]
samp
samp <- samp[c(1,4,10,17,23,24,33,38),]
samp
str(samp)

dev.off()
postscript(file.path(figs_path, "Kinetics_Phosphate.eps"),
					 onefile=FALSE,
					 width=7.48,
					 height=4.72,
					 paper="special",
					 horizontal=FALSE
					 )

# par(mar=c(3.5, 10, 1, 4.5)+0.1)
par(mar=c(3, 3.25, 1, 3.25)+0.1)
with(cycleobj, {
	plot(Cycle, Init.PO4.o, xlim=c(0, 120), ylim=c(0, 14), 
			 axes=FALSE, xlab="", ylab="", type="l", col="red", lwd=2, lty=1)
	points(Cycle, Init.PO4.o, pch=15, col="red")
	# axis(2, ylim=c(0, 14), col="red", lwd = 2, col.axis="red")
	# mtext(2, text="Initial PO4 (mg/L)", line=2, col="red")
	axis(2, ylim=c(0, 14), col="black", lwd = 2, col.axis="black")
	mtext(2, text="PO4 (mg/L)", line=2, col="black")
	
	par(new=TRUE)
	plot(Cycle, Init.PO4.s, xlim=c(0, 120), ylim=c(0, 14), 
			 axes=FALSE, xlab="", ylab="", type="l", col="red", lwd=2, lty=3)
	points(Cycle, Init.PO4.s, pch=17, col="red")
	
	for (i in 1:7) {abline(v=samp$Cycle, col="slategray4")}
	samp %$% text(x=Cycle, y=rep(-0.15,7), label=Sample.o, cex=1, font=2)
	samp %$% text(x=Cycle, y=rep(1.15,2), label=Sample.s, cex=1, font=2)
	
	par(new=TRUE)
	plot(Cycle, Fin.PO4.o, xlim=c(0, 120), ylim=c(0, 14), 
			 axes=FALSE, xlab="", ylab="", type="l", col="blue4", lwd=2, lty=1)
	points(Cycle, Fin.PO4.o, pch=15, col="blue4")
	# axis(2, ylim=c(0, 14), col="blue4", lwd = 2, line=3.5, col.axis="blue4")
	# mtext(2, text="Final PO4 (mg/L)", line=5.5, col="blue4")
	
	par(new=TRUE)
	plot(Cycle, Fin.PO4.s, xlim=c(0, 120), ylim=c(0, 14), 
			 axes=FALSE, xlab="", ylab="", type="l", col="blue4", lwd=2, lty=3)
	points(Cycle, Fin.PO4.s, pch=17, col="blue4")
	
	par(new=TRUE)
	plot(Cycle, Removal.PO4.o, xlim=c(0, 120), ylim=c(0, 100), 
			 axes=FALSE, xlab="", ylab="", type="l", col="darkgreen", lwd=2, lty=1)
	points(Cycle, Removal.PO4.o, pch=15, col="darkgreen")
	# axis(2, ylim=c(0, 100), col="darkgreen", lwd = 2, line=7, col.axis="darkgreen")
	# mtext(2, text="Removal PO4 (%)", line=9, col="darkgreen")
	axis(4, ylim=c(0, 100), col="black", lwd = 2, col.axis="black")
	mtext(4, text="SV30 and Removal PO4 (%)", line=2, col="black")
	
	par(new=TRUE)
	plot(Cycle, Removal.PO4.s, xlim=c(0, 120), ylim=c(0, 100), 
			 axes=FALSE, xlab="", ylab="", type="l", col="darkgreen", lwd=2, lty=3)
	points(Cycle, Removal.PO4.s, pch=17, col="darkgreen")
	
	par(new=TRUE)
	plot(Cycle, SV30.o, xlim=c(0, 120), ylim=c(0, 100), 
			 axes=FALSE, xlab="", ylab="", type="l", col="black", lwd=2, lty=1)
	points(Cycle, SV30.o, pch=15, col="black")
	# axis(4, ylim=c(0, 100), col="black", lwd = 2, line=1.5)
	# mtext(4, text="SV30 (%)", line=3.5)
	
	par(new=TRUE)
	# plot(Cycle, SV30.s, xlim=c(0, 120), ylim=c(0, 100), 
	# 		 axes=FALSE, xlab="", ylab="", type="l", col="black", lwd=2, lty=3)
	points(Cycle, SV30.s, pch=17, col="black")
	cycleobj %$% 
		cycleobj[complete.cases(Cycle, SV30.s),] %$% 
			lines(Cycle, SV30.s, xlim=c(0, 120), ylim=c(0, 100), 
			 type="l", col="black", lwd=2, lty=3)
	
	axis(1, xlim=c(min(Cycle), max(Cycle)), lwd=2)
	mtext(1, text="Cycle Number", line=2)

})	
dev.off()

#----------------------------------------------------------
########
# 
########
p <- read.csv(file.path(data_path, "Phos-Aero-vs-Anaero.txt"), header=TRUE, sep="\t")
str(p)
p$Time<-as.numeric(p$Time)
summary(p)
p

phosplot<-
ggplot(data=p, aes(Time, X15, color="15")) +
	theme_bw() +
	theme(panel.grid.major = element_blank()) + 
	theme(panel.grid.minor = element_blank()) +
	scale_x_continuous(limits=c(0, 150)) + 
	scale_y_continuous(limits=c(0, 100)) +
	theme(legend.position = "bottom") + 
	#theme(legend.title = element_blank()) +
	labs(color="Cycle") +
	xlab("Time (minutes)") + 
	ylab("Phosphate (mg/L)") +
	geom_line() +
	facet_grid(Phase~.) +
	geom_line(data=p, aes(y=X40, color="40")) +
	geom_line(data=p, aes(y=X75, color="75")) + 
	geom_line(data=p, aes(y=X90, color="90")) 
ggsave(file.path(figs_path,"Kinetics_Aero.vs.Anaero.eps"),
			 phosplot,
			 height=2*90+30,
			 width=90,
			 units="mm")

#----------------------------------------------------------
########
# 
########
phb <- read.csv(file.path(data_path, "Phos-Glycogen-PHB.txt"), header=TRUE, sep="\t")
str(phb)
summary(phb)

dev.off()
postscript(file.path(figs_path, "Kinetics_Phos-PHB-Glycogen.eps"),
					 onefile=FALSE,
					 width=7.48,
					 height=4.72,
					 paper="special",
					 horizontal=FALSE
					 )
par(mar=c(3.25, 10, 0.5, 0.5)+0.1)

with(phb, {
	plot(Time, PO4, xlim=c(0, 7), ylim=c(0, 80), 
			 axes=FALSE, xlab="", ylab="", type="l", col="red", lwd=2, lty=1,
			 panel.first=rect(0,-1e6,2.5,80, col="lightblue"))
	points(Time, PO4, pch=20, col="red")
	axis(2, ylim=c(0, 80), col="red", lwd = 2, col.axis="red")
	mtext(2, text="Phosphate (mg/L)", line=2, col="red")
	
	par(new=TRUE)
	plot(Time, Glycogen, xlim=c(0, 7), ylim=c(0, 20), 
			 axes=FALSE, xlab="", ylab="", type="l", col="blue4", lwd=2, lty=2)
	points(Time, Glycogen, pch=20, col="blue4")
	axis(2, ylim=c(0, 20), col="blue4", lwd = 2, line=3.5, col.axis="blue4")
	mtext(2, text="Glycogen (mg/L)", line=5.5, col="blue4")
	
	par(new=TRUE)
	plot(Time, PHB, xlim=c(0, 7), ylim=c(0, 500), 
			 axes=FALSE, xlab="", ylab="", type="l", col="darkgreen", lwd=2, lty=3)
	points(Time, PHB, pch=20, col="darkgreen")
	axis(2, ylim=c(0, 500), col="darkgreen", lwd = 2, line=7, col.axis="darkgreen")
	mtext(2, text="PHB (mg/L)", line=9, col="darkgreen")
	
	axis(1, xlim=c(min(Time), max(Time)), lwd=2)
	mtext(1, text="Time (Hr)", line=2)
})	
dev.off()

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########