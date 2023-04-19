

HL <- read.csv2(paste0(path,"/data/HL_4_years.csv"),sep=";");
SL <- read.csv2(paste0(path,"/data/SL_4_years.csv"),sep=";");
HH <- read.csv2(paste0(path,"/data/HH_4_years.csv"),sep=";");


source(paste0(path,"\\01_set_RDB_Data.R"))
source(paste0(path,"\\01_set_RDB_Data_unraised.R"))
source(paste0(path,"\\02_find_modes_antimodes.R"))
source(paste0(path,"\\03_minimal_reference_subsample_construction.R"))
source(paste0(path,"\\04_compute_distance.R"))




mydata <- setData(HL,SL,HH,selected.species=c("COD"), selected.country=c("DEU"), 
selected.area=c("27.4.a", "27.4.b","27.4.c"), selected.quarter=3, selected.year=2015:2018, delta=5)

#mydata_unraised <- setData_unraised(HL,SL,HH,selected.species=c("COD"), selected.country=c("DEU"), 
#selected.area=c("27.4.a", "27.4.b","27.4.c"), selected.quarter=3, selected.year=2015:2018, delta=5)

########################

mydata_1 <- setData_unraised(HL,SL,HH,selected.species=c("COD"), selected.country=c("DEU"), 
selected.area=c("27.4.a", "27.4.b","27.4.c"), selected.quarter=3, selected.year=2015:2018, delta=5)

U1 <- subset(mydata_1, Year==2016)

aggregate(Number_at_length_1 ~ Year + Trip_number, data=mydata, length)

#######################


#mydata <- setData(selected.species="Pleuronectes platessa", selected.country=c("DEU"), selected.area=c("27.4.a", "27.4.b","27.4.c"), selected.quarter=3, selected.year=2016, bin=5)
#mydata <- setData(selected.species="Pleuronectes platessa", selected.country=c("DEU"), selected.area=c("27.4.a", "27.4.b","27.4.c"), selected.quarter=3, selected.year=2016, bin=3)

subsample_annual <- subset(mydata, Year==2018);
subsample_annual_1 <- subset(mydata_1, Year==2018);

find_modes_antimodes(subsample_annual, smoothed=TRUE, delta=5)
find_modes_antimodes(subsample_annual, smoothed=TRUE, delta=5, important=c(40:80))

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=subsample_annual, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=subsample_annual, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";

M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

M_star <- MEAN_L;

#a
mysubsample_2018 <- subsampling(subsample_annual, delta=5, gamma=0.2);

ks.test(subsample_annual$LengthClass_cm, mysubsample_2018$LengthClass_cm)

cdf_adv(subsample_annual, mysubsample_2018)

#b
mysubsample_2018 <- subsampling(subsample_annual, delta=5, important = c(0:90), gamma=0.2);

ks.test(subset(subsample_annual, LengthClass_cm>=0 & LengthClass_cm<=90)$LengthClass_cm, subset(mysubsample_2018, LengthClass_cm>0 & LengthClass_cm<=90)$LengthClass_cm)

cdf_adv(subsample_annual, mysubsample_2018, important = c(0:90))

#c
mysubsample_2018 <- subsampling(subsample_annual, delta=5, important = c(40:80), epsilon=3, theta=0.7, gamma=0.2);

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=mysubsample_2018, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=mysubsample_2018, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";

M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));


ks.test(subset(subsample_annual, LengthClass_cm>=40 & LengthClass_cm<=80)$LengthClass_cm, subset(mysubsample_2018, LengthClass_cm>=40 & LengthClass_cm<=80)$LengthClass_cm)

cdf_adv(subsample_annual, mysubsample_2018, important = c(40:80))



#################### MEAN, SEM of original data

aggr_per_trip <- aggregate(NoAtLength_RaisedNormed ~ Trip_number, data=subsample_annual, FUN=length)



###############################




z <- aggregate(NoAtLength_RaisedNormed ~ Trip_number, data=subsample_annual, FUN=length)

z1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number, data=subset(subsample_annual, LengthClass_cm %in% c(40:80)), FUN=length)


######## exclude

Trips <- c("DEU20187334916","DEU20180681417","DEU20187334917")

reduced_subsample <- subset(subsample_annual, Trip_number!=Trips[3]);

find_modes_antimodes(reduced_subsample, smoothed=TRUE, delta=5, important=c(40:80))

cdf_adv(subsample_annual, reduced_subsample, reference=FALSE, important = c(40:80))

w <- aggregate(NoAtLength_RaisedNormed ~ Trip_number, data=reduced_subsample, FUN=length)
w1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number, data=subset(reduced_subsample, LengthClass_cm %in% c(40:80)), FUN=length)


#reduced_subsample <- subset(subsample_annual, !(Trip_number %in% Trips[c(1,3)]));


reduced_subsample_trips_1_3 <- subset(subsample_annual, Trip_number %in% Trips[c(1,3)]);
dim(reduced_subsample_trips_1_3)
dim(subset(reduced_subsample_trips_1_3, LengthClass_cm %in% c(40:80)))
cdf_adv(subsample_annual, reduced_subsample_trips_1_3, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_trips_1_3, smoothed=TRUE, delta=5, important=c(40:80))

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_1_3, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_1_3, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";

M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;


reduced_subsample_trips_1_2 <- subset(subsample_annual, Trip_number %in% Trips[c(1,2)]);
dim(reduced_subsample_trips_1_2)
dim(subset(reduced_subsample_trips_1_2, LengthClass_cm %in% c(40:80)))
cdf_adv(subsample_annual, reduced_subsample_trips_1_2, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_trips_1_2, smoothed=TRUE, delta=5, important=c(40:80))

reduced_subsample_trips_2 <- subset(subsample_annual, Trip_number %in% Trips[c(2)]);

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_2, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_2, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";

M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;

dim(reduced_subsample_trips_2)
dim(subset(reduced_subsample_trips_2, LengthClass_cm %in% c(40:80)))
cdf_adv(subsample_annual, reduced_subsample_trips_2, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_trips_2, smoothed=TRUE, delta=5, important=c(40:80))


number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_1_2, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_1_2, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";

M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;

################# Trips remained

reduced_subsample_trips_2_3 <- subset(subsample_annual, Trip_number %in% Trips[c(2,3)]);
cdf_adv(subsample_annual, reduced_subsample_trips_2_3, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_trips_2_3, smoothed=TRUE, delta=5, important=c(40:80))

#################### MEAN, SEM of data

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_2_3, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_2_3, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";

sum_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number + Station_number, data=reduced_subsample_trips_2_3, FUN=sum)
sum_length_haul <- merge(number_per_trip_sub1, sum_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(sum_length_haul)[which(names(sum_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(sum_length_haul)[which(names(sum_length_haul)=="LengthClass_cm")] <- "sum_length";

M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;

#MEAN_L <- sum(mean_length_haul$mean_length)/nrow(mean_length_haul);
#SE_MEAN_L <- sqrt(sum((mean_length_haul$mean_length - MEAN_L)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

#MEAN_L <- sum(sum_length_haul$sum_length)/sum(sum_length_haul$n_stat);
#SE_MEAN_L <- sqrt(sum((sum_length_haul$sum_length - sum_length_haul$n_stat*MEAN_L)^2)/(nrow(sum_length_haul)*(nrow(sum_length_haul)-1)));

library(lme4) 
real.world.sample <- reduced_subsample_trips_2_3;
real.world.sample$Trip_number <- as.factor(real.world.sample$Trip_number);
real.world.sample$Station_number <- as.factor(real.world.sample$Station_number);
real.world.sample$Mix_Trip_Station <- with(real.world.sample, factor(Trip_number:Station_number)); 
length.mixed <- lmer(LengthClass_cm ~ 1 + (1 | Trip_number) + (1 | Mix_Trip_Station), data=real.world.sample);




######################### Catch weight

#SS <- subset(SL, Trip_number %in% Trips[c(2,3)] & Year==2018 & Species==126436);

SS <- subset(SL, Trip_number %in% Trips[c(2,3)] & Year==2018);

SS.sum <- aggregate(Weight ~ Trip_number + Station_number + Year, data=SS, FUN=sum, na.rm=TRUE);

SS.sum$Weight <- (SS.sum$Weight/1000)/1000;

SSub.sum <- subset(SS.sum, Weight>1.5);

reduced_subsample_trips_2_3_hauls_more_n_t <- c();

for (j in Trips[c(2,3)])
{
st <- subset(SSub.sum, Trip_number==j)$Station_number;
reduced_subsample_trips_2_3_hauls_more_n_t <- rbind(reduced_subsample_trips_2_3_hauls_more_n_t, subset(reduced_subsample_trips_2_3, Trip_number==j & Station_number %in% st));
}

cdf_adv(subsample_annual, reduced_subsample_trips_2_3_hauls_more_n_t, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_trips_2_3_hauls_more_n_t, smoothed=TRUE, delta=5, important=c(40:80))

reduced_subsample_trips_2_3_hauls_more_n_t%>%distinct(Trip_number, Station_number)
reduced_subsample_trips_2_3%>%distinct(Trip_number, Station_number)

dim(subset(reduced_subsample_trips_2_3_hauls_more_n_t, LengthClass_cm %in% c(40:80)))
dim(reduced_subsample_trips_2_3_hauls_more_n_t)

#################### MEAN, SEM of data


library(lme4) 
real.world.sample <- reduced_subsample_trips_2_3_hauls_more_n_t;
real.world.sample$Trip_number <- as.factor(real.world.sample$Trip_number);
real.world.sample$Station_number <- as.factor(real.world.sample$Station_number);
real.world.sample$Mix_Trip_Station <- with(real.world.sample, factor(Trip_number:Station_number)); 
length.mixed <- lme(fixed=LengthClass_cm ~ 1, random = ~ 1|Trip_number / Station_number, data=real.world.sample);
#lmer(LengthClass_cm ~ 1 + (1 | Trip_number) + (1 | Mix_Trip_Station), data=real.world.sample);

################################

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_2_3_hauls_more_n_t, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_2_3_hauls_more_n_t, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";
M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);

SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;


######################### Day or night


HHS <- subset(HH, Trip_number %in% Trips[c(2,3)] & Year==2018);

HHS$DayNight <- ifelse(strftime(as.POSIXct(as.character(HHS$Time), format = "%H:%M", usetz = FALSE),"%H:%M") > "03:00" & 
strftime(as.POSIXct(as.character(HHS$Time), format = "%H:%M", usetz = FALSE),"%H:%M") < "21:00", "Day", "Night");

HHSub <- subset(HHS, DayNight=="Day");

reduced_subsample_trips_2_3_hauls_daytime <- c();

for (j in Trips[c(2,3)])
{
st <- subset(HHSub, Trip_number==j)$Station_number;
reduced_subsample_trips_2_3_hauls_daytime <- rbind(reduced_subsample_trips_2_3_hauls_daytime, subset(reduced_subsample_trips_2_3, Trip_number==j & Station_number %in% st));
}

cdf_adv(subsample_annual, reduced_subsample_trips_2_3_hauls_daytime, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_trips_2_3_hauls_daytime, smoothed=TRUE, delta=5, important=c(40:80))

reduced_subsample_trips_2_3_hauls_daytime%>%distinct(Trip_number, Station_number)

#################### MEAN, SEM of data

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_2_3_hauls_daytime, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_2_3_hauls_daytime, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";
M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;

library(lme4) 
real.world.sample <- reduced_subsample_trips_2_3_hauls_daytime;
real.world.sample$Trip_number <- as.factor(real.world.sample$Trip_number);
real.world.sample$Station_number <- as.factor(real.world.sample$Station_number);
real.world.sample$Mix_Trip_Station <- with(real.world.sample, factor(Trip_number:Station_number)); 
length.mixed <- lme(fixed=LengthClass_cm ~ 1, random = ~ 1|Trip_number / Station_number, data=real.world.sample);
#lmer(LengthClass_cm ~ 1 + (1 | Trip_number) + (1 | Mix_Trip_Station), data=real.world.sample);



########### Number of individuals

SS <- subset(SL, Trip_number %in% Trips[c(2,3)] & Year==2018);

SS.sum <- aggregate(Weight ~ Trip_number + Station_number + Year, data=SS, FUN=sum, na.rm=TRUE);

SS.sum$Weight <- (SS.sum$Weight/1000)/1000;

SSub.sum <- subset(SS.sum, Weight<=3);

######################### mode function

stat_mode <- function(x, return_multiple = TRUE, na.rm = FALSE) {
  if(na.rm){
    x <- na.omit(x)
  }
  ux <- unique(x)
  freq <- tabulate(match(x, ux))
  mode_loc <- if(return_multiple) which(freq==max(freq)) else which.max(freq)
  return(ux[mode_loc])
}
############################

modes_vector <- c();
antimodes_vector <- c();
dist_vector <- c();
dim_important <- c();
mean_l <- c()
se_l <- c()
mse_l <- c()

procent <- 0.9;


for (t in c(1:2000))
{
cat("Replicate number = ", t, fill=TRUE);

reduced_subsample_trips_2_3_hauls_daytime_ind <- c();

for (j in Trips[c(2,3)])
{
trip_rev <- subset(reduced_subsample_trips_2_3_hauls_daytime, Trip_number==j);
SSub.sum.st <- subset(SSub.sum, Trip_number==j);
stations_to_reduce <- unique(trip_rev$Station_number);

for(i in stations_to_reduce)
{
st_trip_rev <- subset(trip_rev, Station_number==i);
if(i %in% unique(SSub.sum.st$Station_number))
{
reduced_subsample_trips_2_3_hauls_daytime_ind <- rbind(reduced_subsample_trips_2_3_hauls_daytime_ind, st_trip_rev) 
} else 
{
I <- sample(1:nrow(st_trip_rev),round(procent*nrow(st_trip_rev)), replace=FALSE);
UI <- st_trip_rev[I,];
reduced_subsample_trips_2_3_hauls_daytime_ind <- rbind(reduced_subsample_trips_2_3_hauls_daytime_ind, UI);
}

}

}

#ml <- mean(reduced_subsample_trips_2_3_hauls_daytime_ind$LengthClass_cm)
#se <- sd(reduced_subsample_trips_2_3_hauls_daytime_ind$LengthClass_cm)/sqrt(length(reduced_subsample_trips_2_3_hauls_daytime_ind$LengthClass_cm))

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_2_3_hauls_daytime_ind, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_2_3_hauls_daytime_ind, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";
M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;

library(lme4) 
real.world.sample <- reduced_subsample_trips_2_3_hauls_daytime_ind;
real.world.sample$Trip_number <- as.factor(real.world.sample$Trip_number);
real.world.sample$Station_number <- as.factor(real.world.sample$Station_number);
length.mixed <- lmer(LengthClass_cm ~ 1 + (1 | Trip_number) + (1 | Station_number), data=reduced_subsample_trips_2_3_hauls_daytime_ind)

cat("DIM: ", dim(reduced_subsample_trips_2_3_hauls_daytime_ind), fill=TRUE)
cat("MODES: ", find_modes_antimodes(reduced_subsample_trips_2_3_hauls_daytime_ind, smoothed=TRUE, delta=5, important=c(40:80))$modes_robust, fill=TRUE)
cat("ANTIMODES: ", find_modes_antimodes(reduced_subsample_trips_2_3_hauls_daytime_ind, smoothed=TRUE, delta=5, important=c(40:80))$antimodes_robust, fill=TRUE)

 
dist_vector <- c(dist_vector, cdf_adv(subsample_annual, reduced_subsample_trips_2_3_hauls_daytime_ind, important = c(40:80), reference=FALSE, plotting = FALSE)[1]);

modes_vector <- rbind(modes_vector, find_modes_antimodes(reduced_subsample_trips_2_3_hauls_daytime_ind, smoothed=TRUE, delta=5, important=c(40:80))$modes_robust);
antimodes_vector <- rbind(antimodes_vector, find_modes_antimodes(reduced_subsample_trips_2_3_hauls_daytime_ind, smoothed=TRUE, delta=5, important=c(40:80))$antimodes_robust);

mean_l <- c(mean_l, MEAN_L);
se_l <- c(se_l, SE_MEAN_L)
mse_l <- c(mse_l, MSE)

dim_important <- c(dim_important, nrow(subset(reduced_subsample_trips_2_3_hauls_daytime_ind, LengthClass_cm %in% c(40:80))))

cat("",fill=TRUE);
cat("*****************************************",fill=TRUE);

}

print(apply(modes_vector,2,stat_mode))
print(apply(antimodes_vector,2,stat_mode))
print(mean(dist_vector))
print(mean(mean_l))
print(mean(se_l))
print(mean(mse_l))

reduced_subsample_trips_2_3_hauls_daytime_ind%>%distinct(Trip_number, Station_number)
#reduced_subsample_trips_2_3_hauls_daytime%>%distinct(Trip_number, Station_number)


########################################## 




SS <- subset(SL, Trip_number %in% Trips[c(1,2)] & Year==2018);

SS.sum <- aggregate(Weight ~ Trip_number + Station_number + Year, data=SS, FUN=sum, na.rm=TRUE);

SS.sum$Weight <- (SS.sum$Weight/1000)/1000;

SSub.sum <- subset(SS.sum, Weight>1.5);

reduced_subsample_trips_1_2_hauls_more_n_t <- c();

for (j in Trips[c(1,2)])
{
st <- subset(SSub.sum, Trip_number==j)$Station_number;
reduced_subsample_trips_1_2_hauls_more_n_t <- rbind(reduced_subsample_trips_1_2_hauls_more_n_t, subset(reduced_subsample_trips_1_2, Trip_number==j & Station_number %in% st));
}

cdf_adv(subsample_annual, reduced_subsample_trips_1_2_hauls_more_n_t, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_trips_1_2_hauls_more_n_t, smoothed=TRUE, delta=5, important=c(40:80))

reduced_subsample_trips_1_2_hauls_more_n_t%>%distinct(Trip_number, Station_number)
reduced_subsample_trips_1_2%>%distinct(Trip_number, Station_number)

dim(subset(reduced_subsample_trips_1_2_hauls_more_n_t, LengthClass_cm %in% c(40:80)))
dim(reduced_subsample_trips_1_2_hauls_more_n_t)

#################### MEAN, SEM of data

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_1_2_hauls_more_n_t, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_1_2_hauls_more_n_t, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";
M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);

SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;


######################### Day or night


HHS <- subset(HH, Trip_number %in% Trips[c(2,3)] & Year==2018);

HHS$DayNight <- ifelse(strftime(as.POSIXct(as.character(HHS$Time), format = "%H:%M", usetz = FALSE),"%H:%M") > "03:00" & 
strftime(as.POSIXct(as.character(HHS$Time), format = "%H:%M", usetz = FALSE),"%H:%M") < "21:00", "Day", "Night");

HHSub <- subset(HHS, DayNight=="Day");

reduced_subsample_trips_1_2_hauls_daytime <- c();

for (j in Trips[c(2,3)])
{
st <- subset(HHSub, Trip_number==j)$Station_number;
reduced_subsample_trips_1_2_hauls_daytime <- rbind(reduced_subsample_trips_1_2_hauls_daytime, subset(reduced_subsample_trips_1_2, Trip_number==j & Station_number %in% st));
}

cdf_adv(subsample_annual, reduced_subsample_trips_1_2_hauls_daytime, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_trips_1_2_hauls_daytime, smoothed=TRUE, delta=5, important=c(40:80))

reduced_subsample_trips_1_2_hauls_daytime%>%distinct(Trip_number, Station_number)

#################### MEAN, SEM of data

number_per_trip_sub1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number + Station_number, data=reduced_subsample_trips_1_2_hauls_daytime, FUN=length)
mean_per_trip_sub1 <- aggregate(LengthClass_cm ~ Trip_number  + Station_number, data=reduced_subsample_trips_1_2_hauls_daytime, FUN=mean)
mean_length_haul <- merge(number_per_trip_sub1, mean_per_trip_sub1, by=c("Trip_number","Station_number"), all.x=TRUE)
names(mean_length_haul)[which(names(mean_length_haul)=="NoAtLength_RaisedNormed")] <- "n_stat";
names(mean_length_haul)[which(names(mean_length_haul)=="LengthClass_cm")] <- "mean_length";
M_wave <- mean(mean_length_haul$n_stat);

MEAN_L <- sum(mean_length_haul$mean_length*mean_length_haul$n_stat)/sum(mean_length_haul$n_stat);
SE_MEAN_L <- sqrt(sum(((mean_length_haul$mean_length - MEAN_L)^2) * (mean_length_haul$n_stat/M_wave)^2)/(nrow(mean_length_haul)*(nrow(mean_length_haul)-1)));

MSE <- (MEAN_L - M_star)^2 + SE_MEAN_L^2;

library(lme4) 
real.world.sample <- reduced_subsample_trips_1_2_hauls_daytime;
real.world.sample$Trip_number <- as.factor(real.world.sample$Trip_number);
real.world.sample$Station_number <- as.factor(real.world.sample$Station_number);
real.world.sample$Mix_Trip_Station <- with(real.world.sample, factor(Trip_number:Station_number)); 
length.mixed <- lme(fixed=LengthClass_cm ~ 1, random = ~ 1|Trip_number / Station_number, data=real.world.sample);
#lmer(LengthClass_cm ~ 1 + (1 | Trip_number) + (1 | Mix_Trip_Station), data=real.world.sample);

########################################








w1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number, data=subset(reduced_subsample_trips_2_3, LengthClass_cm %in% c(40:80)), FUN=length)

w1 <- aggregate(NoAtLength_RaisedNormed ~ Trip_number, data=subset(reduced_subsample_trips_2_3, LengthClass_cm %in% c(40:80)), FUN=length)

U1 <- reduced_subsample_trips_2_3%>%distinct(Trip_number, Station_number, Weight)

SL_2_3 <- subset(SL, Trip_number %in% Trips[c(2,3)] & Species=="126436")
U2 <- SL_2_3%>%distinct(Trip_number, Station_number, Weight)







reduced_subsample_for_hauls_2 <- subset(reduced_subsample_for_hauls_1, !(Trip_number=="DEU20187334917" & Station_number %in% c(1)))
#cdf_adv(subsample_annual, reduced_subsample_for_hauls_2, important = c(40:80), reference=FALSE);
find_modes_antimodes(reduced_subsample_for_hauls_2, smoothed=TRUE, delta=5, important=c(40:80))


reduced_subsample_for_hauls_3 <- subset(reduced_subsample_for_hauls_2, !(Trip_number=="DEU20180681417" & Station_number %in% c(1:6)))




reduced_subsample_trips_2_3%>%distinct(Trip_number, Station_number, Weight)









###

#excl <- c(215653,195651, 205557, 325049);

reduced_subsample <- subset(subsample_annual, !(Trip_number %in% excl));

find_modes_antimodes(subsample_annual, smoothed=TRUE, delta=5, important=c(30:80))
find_modes_antimodes(reduced_subsample, smoothed=TRUE, delta=5, important=c(30:80))

cdf_adv(subsample_annual, reduced_subsample, important = c(30:80))

###

excl <- c(245051);

reduced_subsample <- subset(subsample_annual, !(Trip_number %in% excl));

find_modes_antimodes(subsample_annual, smoothed=TRUE, delta=5, important=c(30:80))
find_modes_antimodes(reduced_subsample, smoothed=TRUE, delta=5, important=c(30:80))

cdf_adv(subsample_annual, reduced_subsample, important = c(30:80))


#######################################




















min_per_length_class <- c();
gamma=10
important <- c(0:5000)
M <- max(subset(mydata, Year==2015)$LengthClass_cm);
cl <- seq(0, M, by=1);

for (j in cl)
{
if(j %in% important) n.cl <- gamma else n.cl <- 0;
min_per_length_class <- c(min_per_length_class, round(n.cl));
}

U <- data.frame(table(factor(subset(mydata, Year==2015)$LengthClass_cm, levels=seq(0, M, by=1))));
U <- cbind(U, min_per_length_class);

names(U) <- c("LengthClass_cm", "freq", "min_num_lc");

U$min_num_lc <- ifelse(U$min_num_lc>U$freq, U$freq, U$min_num_lc)

##************************


U1 <- data.frame(table(factor(mysubsample_2015$LengthClass_cm, levels=seq(0, M, by=1))));

names(U1) <- c("LengthClass_cm", "freq");







find_modes_antimodes(mysubsample_2015, smoothed=TRUE, delta=5, important=c(30:90))

cdf_adv(mydata,SUBSAMPLE)


excl <- c("OTB_DEF_>=120_0_0", "TBB_CRU_16-31_0_0")

new_subsample <- subset(mydata, !(MetierLvl6 %in% excl))

cdf_minimizing_samp_effort(mydata, new_subsample, excl)

######################################################################################

M <- max(subsample_annual$LengthClass_cm);
gamma=0.6
important = c(0:120);

cl <- seq(0, M, by=1);

if (gamma < 1)
{
min_per_length_class <- c();

for (j in cl)
{
n_in_j <- nrow(subset(subsample_annual, LengthClass_cm==j));
if(!(j %in% important)) n.cl <- n_in_j else n.cl <- gamma*n_in_j;
min_per_length_class <- c(min_per_length_class, round(n.cl));
}

U <- data.frame(table(factor(subsample_annual$LengthClass_cm, levels=seq(0, M, by=1))));
U <- cbind(U, min_per_length_class);

names(U) <- c("LengthClass_cm", "freq", "min_num_lc");

U$LengthClass_cm <- as.numeric(as.character(U$LengthClass_cm));
U$freq <- as.numeric(as.character(U$freq));
U$min_num_lc <- as.numeric(as.character(U$min_num_lc));

} else

{
min_per_length_class <- c();

for (j in cl)
{
if(!(j %in% important)) n.cl <- n_in_j else n.cl <- gamma;
min_per_length_class <- c(min_per_length_class, round(n.cl));
}

U <- data.frame(table(factor(subsample_annual$LengthClass_cm, levels=seq(0, M, by=1))));
U <- cbind(U, min_per_length_class);

names(U) <- c("LengthClass_cm", "freq", "min_num_lc");

U$min_num_lc <- ifelse(U$min_num_lc>U$freq, U$freq, U$min_num_lc)

U$LengthClass_cm <- as.numeric(as.character(U$LengthClass_cm));
U$freq <- as.numeric(as.character(U$freq));
U$min_num_lc <- as.numeric(as.character(U$min_num_lc));
}





#find_modes_antimodes(mysubsample_2015, delta=5, important = c(30:90));
UU <- data.frame(table(factor(mysubsample_2015$LengthClass_cm, levels=seq(0, M, by=1))));
names(UU) <- c("LengthClass_cm", "freq");

UU <- cbind(UU, min_per_length_class);

##############################################################################################

dev.new();

p1 <- ggplot(data, aes(x=LengthClass_cm)) + 
geom_histogram(data = data, aes(x = LengthClass_cm), binwidth = 1, boundary=0, closed="left", colour="black",fill="turquoise") + 
xlim(0,M)+
#ylim(0,r)+
#annotate("text", label = paste("n =", nrow(data)), x = M-10, y = 7/8*r, size = 9, colour = "black", fontface = 4) +
labs(title="ORIGINAL SAMPLE", x="LENGTH", y = "COUNT") +
theme(title = element_text(face="bold",size=16), axis.title = element_text(face="bold",size = 15), axis.text = element_text(face="bold", size = 14), strip.text = element_text(size=18, face="bold.italic")) 


p2 <- ggplot(subsample, aes(x=LengthClass_cm)) + 
geom_histogram(data = subsample, aes(x = LengthClass_cm), binwidth = 1, boundary=0, closed="left", colour="black",fill="lightskyblue1") + 
xlim(0,M)+
#ylim(0,r)+
#annotate("text", label = paste("n =", nrow(subsample)), x = M-10, y = 7/8*r, size = 9, colour = "black", fontface = 4) +
#annotate("text", label = paste("Excluded: ", paste(excl, collapse="/")), x = M-10, y = 2/5*r, size = 7, colour = "black", fontface = 4) +
labs(title="SUBSAMPLE", x="LENGTH", y = "COUNT") +
theme(title = element_text(face="bold",size=15), axis.title = element_text(face="bold",size = 15), axis.text = element_text(face="bold", size = 14), strip.text = element_text(size=16, face="bold.italic")) 

library(gridGraphics)

grid.arrange(p1, p2, ncol=2, top=textGrob(paste("LFD:", paste(selected.species), ", ", "EU STATE: ", paste(selected.country, collapse="/"), ", Area ",  paste(selected.area, collapse="/"), ", Quarter ", selected.quarter,  ", Year", paste(selected.year,collapse="/" ),"\n"), gp=gpar(fontsize=20, font=2)))




pp <- ggplot(mysubsample_2015, aes(x=LengthClass_cm)) + 
geom_histogram(data = mysubsample_2015, aes(x = LengthClass_cm), binwidth = 1, boundary=0, closed="left", colour="black",fill="yellow") + 
xlim(0,M)+
labs(title="SUBSAMPLE", x="LENGTH", y = "COUNT") +
theme(title = element_text(face="bold",size=15), axis.title = element_text(face="bold",size = 15), axis.text = element_text(face="bold", size = 14), strip.text = element_text(size=16, face="bold.italic")) 

plot(pp)

Jpp <- data.frame(table(factor(mysubsample_2015$LengthClass_cm, levels=seq(0, M, by=1))));
names(Jpp) <- c("LengthClass_cm", "freq")
