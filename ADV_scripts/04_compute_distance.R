#### This script computes L1-distance between original and subsampled data sets
#### Inputs and parameters: data - original data set; subsample - subsampled data set (may be a reference subsample or any "real-world" subsample); important - defined important length classes where distance is calculated; 
#### reference - if TRUE, than subsample = reference subsample and ADV is computed, if FALSE, than subsample = "real-world" subsample and distance is computed and can be than compared to ADV.
#### Ouputs: distance value + plot with empirical cumulative functions of original and subsampled data sets

cdf_adv <- function(data, subsample, important, delta, theta, epsilon, reference = TRUE, plotting = TRUE)

{

c_penalizing_1 <- 10;
c_penalizing_2 <- 2;
c_penalizing_3 <- 1;

M <- max(data$LengthClass_cm);

if (missing(important)) important <- c(0:M);

U <- data.frame(table(factor(data$LengthClass_cm, levels=seq(min(important), max(important), by=1))));

U_sub <- data.frame(table(factor(subsample$LengthClass_cm, levels=seq(min(important), max(important), by=1))));

names(U) <- c("LengthClass_cm","freq");
names(U_sub) <- c("LengthClass_cm","freq");

cp_ori <- sort(c(find_modes_antimodes(data, smoothed=TRUE, delta=delta, important=important)$modes_robust, find_modes_antimodes(data, smoothed=TRUE, delta=delta, important=important)$antimodes_robust));
cp_sub <- sort(c(find_modes_antimodes(subsample, smoothed=TRUE, delta=delta, important=important)$modes_robust, find_modes_antimodes(subsample, smoothed=TRUE, delta=delta, important=important)$antimodes_robust));



PF_full <- ecdf(subset(data, LengthClass_cm %in% important)$LengthClass_cm)
PF_subsample <- ecdf(subset(subsample, LengthClass_cm %in% important)$LengthClass_cm)

 
T <- data.frame(cbind(knots(PF_full), PF_full(knots(PF_full))))
names(T) <- c("x","P")

T_subsample <- data.frame(cbind(knots(PF_subsample), PF_subsample(knots(PF_subsample))))
names(T_subsample) <- c("x","P")

d_full <- setdiff(seq(0,max(data$LengthClass_cm), by=1), T$x);

if (length(d_full)>0) 
{
T_diff_full <- data.frame(cbind(d_full,rep(NA,length(d_full)))); 
names(T_diff_full) <- c("x","P");  
T <- rbind(T, T_diff_full); 
T <- T[order(T$x),];
m <- which(is.na(T$P))
for (i in 1:length(m))
{
if (m[i] > 1)
{
T[m[i],]$P <- T[m[i]-1,]$P
#T[m[i],]$P <- T[m[i]+1,]$P
} else T[m[i],]$P <- 0;
}
}

d_subsample <- setdiff(seq(0,max(data$LengthClass_cm),by=1), T_subsample$x);

if (length(d_subsample)>0) 
{
T_diff_subsample <- data.frame(cbind(d_subsample,rep(NA,length(d_subsample)))); 
names(T_diff_subsample) <- c("x","P");  
T_subsample <- rbind(T_subsample, T_diff_subsample); 
T_subsample <- T_subsample[order(T_subsample$x),];
m <- which(is.na(T_subsample$P))
for (i in 1:length(m))
{
if (m[i] > 1)
{
T_subsample[m[i],]$P <- T_subsample[m[i]-1,]$P
#T_subsample[m[i],]$P <- T_subsample[m[i]+1,]$P
} else T_subsample[m[i],]$P <- 0;
}
}


M_dist <- rbind(seq(0, M,by=1), T$P, T_subsample$P);
M_t_dist <- t(M_dist)
colnames(M_t_dist) <- c("Arg", "F", "G");
M_t_dist <- data.frame(M_t_dist);

penalty_1 <- ifelse((length(cp_ori)==length(cp_sub)),0, 1);

distance <- round(dist(M_dist[2:3,], method ="manhattan"),4) + c_penalizing_1*penalty_1;


if(length(cp_ori)==length(cp_sub) & distance > 0)
{

diff_1 <- abs(cp_ori-cp_sub) - epsilon;
penalty_2 <- 0;
for (i in 1:length(diff_1)) {penalty_2 <- penalty_2 + (max(0, diff_1[i]))};
distance <- distance + c_penalizing_2*penalty_2;


fifi <- subset(U, LengthClass_cm %in% cp_ori);
gigi <- subset(U_sub, LengthClass_cm %in% cp_sub);

if (nrow(gigi)>1) rat1 <- abs(diff(gigi$freq)) else rat1 <- abs(gigi$freq)
if (nrow(fifi)>1) rat2 <- abs(diff(fifi$freq)) else rat2 <- abs(fifi$freq)

rr <- rat1/rat2;
if (any(is.nan(rr))) rr[which(is.nan(rr))] <- 1;

diff <- theta - rr;
penalty_3 <- 0;
for (i in 1:length(diff)) {penalty_3 <- penalty_3 + max(0, diff[i])};
distance <- distance + c_penalizing_3*penalty_3;
}


if (plotting==TRUE)

{

dev.new()

op <- par(mfrow= c(1,1),oma = c(2,4,2,0) + 0.1, las=1, mar=c(2, 2, 3, 2) + 0.1, mgp = c(3, 2, 0))

ma <- mean(subsample$LengthClass_cm)

plot(T$x,T$P, xlim=c(0,M), col="plum3", cex.axis=1.7, cex.main=2, font.main=4, 
#verticals = TRUE, do.points = FALSE, 
lwd=4, type="s", lty=1)
par(new=TRUE)
plot(T_subsample$x,T_subsample$P, xlim=c(0,M), col=ifelse(reference==TRUE, "forestgreen", "tomato"), cex.axis=1.7, main=NA, 
#verticals = TRUE, do.points = FALSE, 
lwd=3, lty=1, type="s")
legend(ma, 0.6, c(paste(ifelse(reference==TRUE, "ADV = ","Distance = "), distance),"original sample", ifelse(reference==TRUE, "reference subsample", "real-world subsample")), 
lwd=3, col=c("black","plum3",ifelse(reference==TRUE, "forestgreen", "tomato")), text.font=c(2,1,1), lty=c(NA,1,1),title=NA, cex=2.2, bty="n", x.intersp=0.4)

#title(paste("ECDF: ",unique(data$Species_code), ", ", "EU State: ", paste(unique(data$Vessel_flag_country), collapse="/"), ", Area ",  paste(unique(data$Area), collapse="/"), ", Quarter ", paste(unique(data$Quarter), collapse="/"),  ", Year", paste(unique(data$Year))), outer=TRUE, cex.main=2.7)

if (reference==FALSE)
{
M.c <- intersect(data$FAC_EC_lvl6, subsample$FAC_EC_lvl6);
T.c <- intersect(data$Trip_number, subsample$Trip_number);

mtext("Subsample includes", side=3, line=-3, adj=0.1, cex=2.6, col="navy");
if (length(M.c)<length(data$FAC_EC_lvl6))
{
mtext(paste0("Fleets: ", paste(M.c,collapse=",")), side=3, line=-6, adj=0.1, cex=2.6, col="navy");
}
if (length(T.c)<length(data$Trip_number))
{
mtext(paste0("Trips: ", paste(T.c,collapse=",")), side=3, line=-9, adj=0.1, cex=2.6, col="navy");
}
}


par(op)

}


return(distance)

}

