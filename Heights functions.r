#Merge with plot type, add palm etc
merge4heights<-function(mfp,PlotType){
	names(PlotType)<- c('PlotID','PlotCode_B','ClusterID','ForestMoistureID','ForestMoistureName',
                      'ForestEdaphicID','ForestEdaphicName', 'EdaphicHeight',	
                      'EdaphicHeightCode','ForestElevationID', 'ForestElevationName',
                      'ElevationHeight','ElevationHeightCode',	'BiogeographicalRegionID',
                       'BiogeographicalRegionName', 'ContinentID','CountryID')
	PlotType$ClusterID<- as.numeric(gsub("=","",PlotType$ClusterID))
	PlotType$ClusterID<- as.numeric(gsub("NULL","NA",PlotType$ClusterID))

	PlotType$ForestMoistureID<- as.numeric(gsub("=","",PlotType$ForestMoistureID))
	PlotType$ForestMoistureID<- as.numeric(gsub("NULL","NA",PlotType$ForestMoistureID))

	PlotType$ForestEdaphicID<- as.numeric(gsub("=","",PlotType$ForestEdaphicID))
	PlotType$ForestEdaphicID<- as.numeric(gsub("NULL","NA",PlotType$ForestEdaphicID))

	PlotType$EdaphicHeightCode<- as.numeric(gsub("=","",PlotType$EdaphicHeightCode))
	PlotType$EdaphicHeightCode<- as.numeric(gsub("NULL","NA",PlotType$EdaphicHeightCode))

	PlotType$ForestElevationID<- as.numeric(gsub("=","",PlotType$ForestElevationID))
	PlotType$ForestElevationID<- as.numeric(gsub("NULL","NA",PlotType$ForestElevationID))

	PlotType$ElevationHeightCode<- as.numeric(gsub("=","",PlotType$ElevationHeightCode))
	PlotType$ElevationHeightCode<- as.numeric(gsub("NULL","NA",PlotType$ElevationHeightCode))

	PlotType$ContinentID <- as.numeric(gsub("=","",PlotType$ContinentID))
	PlotType$ContinentID <- as.numeric(gsub("NULL","NA",PlotType$ContinentID))

	PlotType$CountryID <- as.numeric(gsub("=","",PlotType$CountryID))
	PlotType$CountryID <- as.numeric(gsub("NULL","NA",PlotType$CountryID))

	Heights1<-AGBChv05MH(mfp)
## Recategorize F5 Methods. This bit of code is useful when doing the comparisons of laser vs clinometer

	Heights1$Method<- ifelse(Heights1$F5==1,1,
                        ifelse(Heights1$F5==6,6,
                               ifelse( Heights1$F5==2 | Heights1$F5==3, 3,4
                               )
                        )
	)

## CHeck number of directly measured trees
#aggregate(TreeID ~ Method, data=Heights1, FUN= length)


# add useful columns. These columns should be incorporated in the next version of the Rpackage      

	Heights1$Palm<- ifelse ( grepl('Arecaceae',Heights1$Family)==TRUE |
                        grepl('Strelitziaceae',Heights1$Family)==TRUE |
                          grepl('Poaceae',Heights1$Family)==TRUE |
                          grepl('Cyatheaceae',Heights1$Family)==TRUE
                        ,1,0)

#to check number of records in each family
#aggregate(TreeID ~ Family +Palm, data = Heights1, FUN = length )


	Heights1$PomChange <- ifelse(grepl('6',Heights1$F4)==TRUE,1,0)


#Merge with PlotType
	Heights2<-merge(Heights1, PlotType, by ='PlotID')
	return(Heights2)
}


#############################################################################
###Fit Weibull models

fit.weib<-function(data,return.mods=FALSE){

# Select data with diameters and height 
# Remove trees without heights or with height=0.  Remove PAlms. To generate D:H relationships main plot views should be used.

TreesHt <- data[data$Height>0 & data$Alive==1 & data$DBH1>90 & data$DBH1<5000    & !is.na(data$Height) & data$Height<90 
                   & data$Palm==0 &  !is.na(data$F5), ]

# Exclude F3-3
TreesHt <- TreesHt[ grepl('3',TreesHt$F3)==FALSE,  ]
# Exclude F4 not like '60' or '0'
TreesHt <- TreesHt[ grepl('0',TreesHt$F4)==TRUE|grepl('06',TreesHt$F4)==TRUE ,  ]
# Exclude f1 flag1= b, c, d, f,g,j,k,m,o,p
TreesHt <- TreesHt[ grepl('b',TreesHt$F1)==FALSE & grepl('c',TreesHt$F1)==FALSE &
                            grepl('d',TreesHt$F1)==FALSE & grepl('f',TreesHt$F1)==FALSE&
                            grepl('g',TreesHt$F1)==FALSE & grepl('j',TreesHt$F1)==FALSE &
                            grepl('k',TreesHt$F1)==FALSE & grepl('m',TreesHt$F1)==FALSE &
                            grepl('o',TreesHt$F1)==FALSE & grepl('p',TreesHt$F1)==FALSE
                    ,  ]

#Exclude treeswith method 1
TreesHt <- TreesHt[!(TreesHt$Method==1),]

## Aggregate to check rules are applied correctly
#aggregate(TreeID ~ F1, data=TreesHt, FUN=length)
#aggregate(TreeID ~ F3, data=TreesHt, FUN=length)
#aggregate(TreeID ~ F4, data=TreesHt, FUN=length)
#aggregate(TreeID ~ Method, data=TreesHt, FUN=length)
#aggregate(TreeID ~ Palm, data=TreesHt, FUN=length)




###1,ContinentData## 
HtCont<- TreesHt[,c('ContinentID','PlotID','Height','DBH1')]

library (nlme)
weib1 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|ContinentID,
                  data=HtCont,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
#summary(weib1)
ContinentCoef<-data.frame(coef(weib1))
colnames(ContinentCoef) <- c('a_Continent','b_Continent', 'c_Continent')
ContinentCoef$ContinentID = rownames(ContinentCoef)
#ContinentCoef
#plot(weib1)

#2. ContinentData and Forest Type
HtCont_Type<- (TreesHt[,c('ContinentID','ForestMoistureID', 'EdaphicHeightCode', 'ElevationHeightCode','Height','DBH1')])
good<-complete.cases(HtCont_Type)
HtCont_Typea<-HtCont_Type[good,]
#head(HtCont_Typea)

weib2 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|ContinentID/ForestMoistureID/EdaphicHeightCode/ElevationHeightCode,
                  data=HtCont_Typea,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
#summary(weib2)
ContinentCoef_Type<-data.frame(coef(weib2))
colnames(ContinentCoef_Type) <- c('a_Continent_T','b_Continent_T', 'c_Continent_T')

ContinentCoef_Type$ContType <-rownames(ContinentCoef_Type)
#head(ContinentCoef_Type)

ContinentCoef_Type$ContinentID<-unlist(lapply(strsplit(ContinentCoef_Type$ContType, "/"),"[", 1))
ContinentCoef_Type$ForestMoistureID<-unlist(lapply(strsplit(ContinentCoef_Type$ContType, "/"),"[", 2))
ContinentCoef_Type$EdaphicHeightCode<-unlist(lapply(strsplit(ContinentCoef_Type$ContType, "/"),"[", 3))
ContinentCoef_Type$ElevationHeightCode<-unlist(lapply(strsplit(ContinentCoef_Type$ContType, "/"),"[", 4))
#ContinentCoef_Type
#plot (weib2)

#3. Biogeographic Region

HtBiogeo<- TreesHt[,c('ContinentID','CountryID','BiogeographicalRegionID', 'PlotID','Height','DBH1')]
weib3 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|'BiogeographicalRegionID',
                  data=HtBiogeo,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
#summary(weib3)
BioRCoef<-data.frame(coef(weib3))
colnames(BioRCoef) <- c('a_BioR','b_BioR', 'c_BioR')
BioRCoef$BiogeographicalRegionID = rownames(BioRCoef)
#BioRCoef
#plot(weib3)
#head(BioRCoef)
##write.csv(BioRCoef,'BioRCoef.csv')


#4. Biogeographic Region/ForestType
HtBiogeoFt<- TreesHt[,c('ContinentID','CountryID','BiogeographicalRegionID', 'PlotID','ForestMoistureID', 'EdaphicHeightCode', 'ElevationHeightCode','Height','DBH1')]
weib4 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|BiogeographicalRegionID/ForestMoistureID/EdaphicHeightCode/ElevationHeightCode,
                  data=HtBiogeoFt,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
#summary(weib4)
BioRCoefFt<-data.frame(coef(weib4))
colnames(BioRCoefFt) <- c('a_BioRF','b_BioRF', 'c_BioRF')
BioRCoefFt$BioRCoefFtype = rownames(BioRCoefFt)
head(BioRCoefFt)
BioRCoefFt$BiogeographicalRegionID<-unlist(lapply(strsplit(BioRCoefFt$BioRCoefFtype, "/"),"[", 1))
BioRCoefFt$ForestMoistureID<-unlist(lapply(strsplit(BioRCoefFt$BioRCoefFtype, "/"),"[", 2))
BioRCoefFt$EdaphicHeightCode<-unlist(lapply(strsplit(BioRCoefFt$BioRCoefFtype, "/"),"[", 3))
BioRCoefFt$ElevationHeightCode<-unlist(lapply(strsplit(BioRCoefFt$BioRCoefFtype, "/"),"[", 4))
#BioRCoefFt
#plot(weib4)


#write.csv(BioRCoefFt,'BioRCoefFt.csv')



##5 Analyze data by country
HtCountry<- TreesHt[,c('ContinentID','CountryID', 'PlotID','Height','DBH1')]
#library (nlme)
weib5 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|CountryID,
                  data=HtCountry,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
#summary(weib5)
CountryCoef<-data.frame(coef(weib5))
colnames(CountryCoef) <- c('a_Country','b_Country', 'c_Country')
CountryCoef$CountryID = rownames(CountryCoef)
#CountryCoef
#plot(weib5)
#write.csv(CountryCoef,'CountryCoef.csv')

##6 Analyze data by countryand ForestType
HtCountryF<- TreesHt[,c('ContinentID','CountryID', 'PlotID','ForestMoistureID', 'EdaphicHeightCode', 'ElevationHeightCode','Height','DBH1')]
#library (nlme)
weib6 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|CountryID/ForestMoistureID/EdaphicHeightCode/ElevationHeightCode,
                  data=HtCountryF,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
#summary(weib6)

HtCountryFt<-data.frame(coef(weib6))
colnames(HtCountryFt) <- c('a_CountryF','b_CountryF', 'c_CountryF')
HtCountryFt$HtCountryFtypw = rownames(HtCountryFt)
#head(HtCountryFt)

HtCountryFt$CountryID<-unlist(lapply(strsplit(HtCountryFt$HtCountryFtypw, "/"),"[", 1))
HtCountryFt$ForestMoistureID<-unlist(lapply(strsplit(HtCountryFt$HtCountryFtypw, "/"),"[", 2))
HtCountryFt$EdaphicHeightCode<-unlist(lapply(strsplit(HtCountryFt$HtCountryFtypw, "/"),"[", 3))
HtCountryFt$ElevationHeightCode<-unlist(lapply(strsplit(HtCountryFt$HtCountryFtypw, "/"),"[", 4))
#HtCountryFt

#plot(weib6)


#write.csv(HtCountryFt,'HtCountryFt.csv')



#7.Analyze data by clusterid
HtCluster<- (TreesHt[,c('ClusterID','Height','DBH1')])
goodCl<-complete.cases(HtCluster)
HtClusterA<-HtCluster[goodCl,]
#head(HtClusterA)

#library (nlme)
weib7 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))|ClusterID,
                  data=HtClusterA,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)

#summary(weib7)

ClusterCoef<-data.frame(coef(weib7))
colnames(ClusterCoef) <- c('a_Cluster','b_Cluster', 'c_Cluster')
ClusterCoef$ClusterID = rownames(ClusterCoef)
#ClusterCoef
#plot(weib7)
#write.csv(ClusterCoef,'ClusterCoef.csv')


#8. Cluster id and ForestType
HtClusterFt <- TreesHt[,c('ClusterID','ForestMoistureID', 'EdaphicHeightCode', 'ElevationHeightCode','Height','DBH1')]
goodClFt<-complete.cases(HtClusterFt)
HtClusterFtA<-HtClusterFt[goodClFt,]
#head(HtClusterFtA)
##write.csv (HtClusterFA, 'HtClusterFA.csv')

#aggregate (DBH1 ~ ClusterID +ForestMoistureID + EdaphicHeightCode +ElevationHeightCode, data = HtClusterFA, FUN = length)
##

#library (nlme)
weib8 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))| ClusterID/ForestMoistureID/EdaphicHeightCode/ElevationHeightCode,
                  data=HtClusterFtA,
                  na.action=na.omit,
                  start = c(a=25, b= 0.05, c= 0.7),
                  pool=FALSE)
#summary(weib8)
ClusterFCoef<-data.frame(coef(weib8))
colnames(ClusterFCoef) <- c('a_ClusterF','b_ClusterF', 'c_ClusterF')
ClusterFCoef$ClusterF = rownames(ClusterFCoef)
#ClusterFCoef

ClusterFCoef$ClusterID<-unlist(lapply(strsplit(ClusterFCoef$ClusterF, "/"),"[", 1))
ClusterFCoef$ForestMoistureID<-unlist(lapply(strsplit(ClusterFCoef$ClusterF, "/"),"[", 2))
ClusterFCoef$EdaphicHeightCode<-unlist(lapply(strsplit(ClusterFCoef$ClusterF, "/"),"[", 3))
ClusterFCoef$ElevationHeightCode<-unlist(lapply(strsplit(ClusterFCoef$ClusterF, "/"),"[", 4))
#head(ClusterFCoef)
#write.csv(ClusterFCoef,'ClusterFCoef.csv')
#plot (weib8)


#Analyze data by  plot
HtPlot<- TreesHt[,c('PlotID','Height','DBH1')]
#library (nlme)
weib9 <- nlsList (Height ~ a*(1-exp(-b*(DBH1/10)^c))| PlotID,
                  data=HtPlot,
                  na.action=na.omit,
                  start = c(a =25, b= 0.05, c= 0.7),
                  pool=FALSE)
#summary(weib9)
PlotCoef<-data.frame(coef(weib9))
colnames(PlotCoef) <- c('a_Plot','b_Plot', 'c_Plot')
PlotCoef$PlotID = rownames(PlotCoef)
#head(PlotCoef)
#nrow(PlotCoef)

#plot (weib9)

#Get unique list of plots
Heights3<-unique(data[,c('PlotID','ClusterID','ForestMoistureID','EdaphicHeightCode','ElevationHeightCode','CountryID','BiogeographicalRegionID','ContinentID')])

Ht1 <- merge(Heights3,PlotCoef, by='PlotID', all.x=TRUE)

#Ht1$HtPlot<-Ht1$a_Plot*(1-exp(-Ht1$b_Plot *(Ht1$DBH4/10)^Ht1$c_Plot))


#Cluster and Forest type
Ht2 <- merge(Ht1,ClusterFCoef, by=c('ClusterID','ForestMoistureID','EdaphicHeightCode','ElevationHeightCode'), all.x=TRUE)

#Ht2$HtClFt<-Ht2$a_ClusterF*(1-exp(-Ht2$b_ClusterF *(Ht2$DBH4/10)^Ht2$c_ClusterF))


# CLUSTER
Ht3<- merge(Ht2, ClusterCoef, by='ClusterID', all.x=TRUE)

#Ht3$HtCl<-Ht3$a_Cluster*(1-exp(-Ht3$b_Cluster *(Ht3$DBH4/10)^Ht3$c_Cluster))


#Country/ForestType
Ht4 <- merge(Ht3,HtCountryFt, by=c('CountryID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'), all.x=TRUE)

#Ht4$HtCtF<-Ht4$a_CountryF*(1-exp(-Ht4$b_CountryF *(Ht4$DBH4/10)^Ht4$c_CountryF))


#Country
Ht5 <- merge(Ht4,CountryCoef, by='CountryID', all.x=TRUE)

#Ht5$HtCt<-Ht5$a_Country*(1-exp(-Ht5$b_Country *(Ht5$DBH4/10)^Ht5$c_Country))


#Biogeographic region and Forest type
Ht6<- merge (Ht5,BioRCoefFt, by=c('BiogeographicalRegionID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'), all.x=TRUE)

#Ht6$HtBf<- Ht6$a_BioRF*(1-exp(-Ht6$b_BioRF*(Ht6$DBH4/10)^Ht6$c_BioRF))


#Biogeographic region 
Ht7<- merge (Ht6,BioRCoef, by='BiogeographicalRegionID', all.x=TRUE)

#Ht7$HtB<- Ht7$a_BioR*(1-exp(-Ht7$b_BioR*(Ht7$DBH4/10)^Ht7$c_BioR))

#Continent and ForestType

Ht8<- merge (Ht7,ContinentCoef_Type, by=c('ContinentID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'), all.x=TRUE)

#Ht8$HtCoF<- Ht8$a_Continent_T*(1-exp(-Ht8$b_Continent_T*(Ht8$DBH4/10)^Ht8$c_Continent_T))


#Continent
Ht9 <- merge(Ht8, ContinentCoef, by='ContinentID', all.x=TRUE)

#Ht9$HtCo<-Ht9$a_Continent*(1-exp(-Ht9$b_Continent*(Ht9$DBH4/10)^Ht9$c_Continent))
if(return.mods==T){
return(list(Ht9,weib9))
}else{
return(Ht9)
}
}

param.merge<-function(data,wparm){

data$MatchLevel<-ifelse(data$PlotID%in%wparm$PlotID,"Plot",
	ifelse(paste(data$ClusterID,data$ForestMoistureID,data$EdaphicHeightCode,data$ElevationHeightCode)%in%paste(wparm$ClusterID,wparm$ForestMoistureID,wparm$EdaphicHeightCode,wparm$ElevationHeightCode),"Cluster",
		ifelse(paste(data$BiogeographicalRegionID,data$ForestMoistureID,data$EdaphicHeightCode,data$ElevationHeightCode)%in%paste(wparm$BiogeographicalRegionID,wparm$ForestMoistureID,wparm$EdaphicHeightCode,wparm$ElevationHeightCode),"Biogeog",
			"Cont")))


d1<-data[data$MatchLevel=="Plot",]
dw1<-merge(data,wparm,by=c('PlotID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'))

d2<-data[data$MatchLevel=="Cluster",]
wparm2<-wparm[!duplicated(paste(wparm$ClusterID,wparm$ForestMoistureID,wparm$EdaphicHeightCode, wparm$ElevationHeightCode)),]
dw2<-merge(d2,wparm2,by=c('ClusterID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'))

d3<-data[data$MatchLevel=="Biogeog",]
wparm3<-wparm[!duplicated(paste(wparm$BiogeographicalRegionID,wparm$ForestMoistureID,wparm$EdaphicHeightCode, wparm$ElevationHeightCode)),]
dw3<-merge(d3,wparm3,by=c('BiogeographicalRegionID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'))

d4<-data[data$MatchLevel=="Cont",]
wparm4<-wparm[!duplicated(paste(wparm$ContinentID,wparm$ForestMoistureID,wparm$EdaphicHeightCode, wparm$ElevationHeightCode)),]
dw4<-merge(d4,wparm4,by=c('ContinentID','ForestMoistureID','EdaphicHeightCode', 'ElevationHeightCode'))

dw1<-dw1[,-c(grep(".x",names(dw1)),grep(".y",names(dw1)))]
dw2<-dw2[,-c(grep(".x",names(dw2)),grep(".y",names(dw2)))]
dw3<-dw3[,-c(grep(".x",names(dw3)),grep(".y",names(dw3)))]
dw4<-dw4[,-c(grep(".x",names(dw4)),grep(".y",names(dw4)))]

dw5<-rbind(dw1,dw2,dw3,dw4)

#Need to add in plots that do not match on forest type at any level
d6<-data[!data$PlotCode%in%dw5$PlotCode,]
d7<-cbind(d6,matrix(nrow=nrow(d6),ncol=(ncol(dw5)-ncol(d6))))
names(d7)[(ncol(d6)+1):ncol(d7)]<-names(dw5)[(ncol(d6)+1):ncol(d7)]
dw.all<-rbind(dw5,d7)

###
dw.all$a_Best<-with(dw.all,
	ifelse(!is.na(a_Plot),a_Plot,
		ifelse(!is.na(a_ClusterF),a_ClusterF,
			ifelse(!is.na(a_BioRF),a_BioRF,
	a_Continent_T))))

dw.all$b_Best<-with(dw.all,
	ifelse(!is.na(b_Plot),b_Plot,
		ifelse(!is.na(b_ClusterF),b_ClusterF,
			ifelse(!is.na(b_BioRF),b_BioRF,
	b_Continent_T))))

dw.all$c_Best<-with(dw.all,
	ifelse(!is.na(c_Plot),c_Plot,
		ifelse(!is.na(c_ClusterF),c_ClusterF,
			ifelse(!is.na(c_BioRF),c_BioRF,
	c_Continent_T))))

dw.all$Weib.parm.source<-with(dw.all,
	ifelse(!is.na(c_Plot),"Plot",
		ifelse(!is.na(c_ClusterF),"Cluster",
			ifelse(!is.na(c_BioRF),"BiogR",
	"Cont"))))
dw.all$HtEst<-dw.all$a_Best*(1-exp(-dw.all$b_Best*(dw.all$DBH4/10)^dw.all$c_Best))
dw.all$HtEst[is.na(dw.all$DBH4)]<-NA
dw.all$HtEst[dw.all$HtEst==0]<-NA

return(dw.all)
}


#####
ht.offset<-function(dw.all){
TreeswithHeights <- dw.all[dw.all$Height>0 & dw.all$Alive==1 & dw.all$DBH1>90 & dw.all$DBH1<5000    & !is.na(dw.all$Height) & dw.all$Height<90 
                    & dw.all$Palm==0 &  !is.na(dw.all$F5), ]

# Exclude F3-3
TreeswithHeights<- TreeswithHeights[ grepl('3',TreeswithHeights$F3)==FALSE,  ]
# Exclude F4 not like '60' or '0'
TreeswithHeights <-TreeswithHeights[ grepl('0',TreeswithHeights$F4)==TRUE|grepl('06',TreeswithHeights$F4)==TRUE ,  ]
# Exclude f1 flag1= b, c, d, f,g,j,k,m,o,p
TreeswithHeights <- TreeswithHeights[ grepl('b',TreeswithHeights$F1)==FALSE & grepl('c',TreeswithHeights$F1)==FALSE &
                      grepl('d',TreeswithHeights$F1)==FALSE & grepl('f',TreeswithHeights$F1)==FALSE&
                      grepl('g',TreeswithHeights$F1)==FALSE & grepl('j',TreeswithHeights$F1)==FALSE &
                      grepl('k',TreeswithHeights$F1)==FALSE & grepl('m',TreeswithHeights$F1)==FALSE &
                      grepl('o',TreeswithHeights$F1)==FALSE & grepl('p',TreeswithHeights$F1)==FALSE
                    ,  ]

#remove tree height estimated by eye
TreeswithHeights<- TreeswithHeights[!(TreeswithHeights$F5 ==1),]

nrow(TreeswithHeights)

#TreeswithHeights[TreeswithHeights$TreeID== 55856,] # test tree with two censuses with height measured

#nrow(TreeswithHeights)

HtM1<-aggregate (Height ~ TreeID+PlotID +PlotViewID, data= TreeswithHeights, FUN=length)
colnames(HtM1) <-c('TreeID','PlotID','PlotViewID','Records')
nrow(HtM1)
#Add number of tree height records information to the table
HtM2<-merge (HtM1, TreeswithHeights, by=c('TreeID','PlotID','PlotViewID'))
head(HtM2)
nrow(HtM2)

#HtM2[HtM2$TreeID==55856,]

#maximum censusnumber for all tree_census_heights
mxHt<- aggregate(Census.No ~TreeID+PlotID+PlotViewID, data=HtM2,FUN=max)

#mxHt[mxHt$TreeID == 55856,]

HtCns<-merge(mxHt,TreeswithHeights, by= c('TreeID', 'PlotID','PlotViewID','Census.No') )

nrow(mxHt)
nrow(HtCns)
head(HtCns)

#HtCns[HtCns$TreeID==55856,]

#Select information necessary for Offset

HtCns2<-HtCns[ , c('TreeID','PlotViewID','PlotID','Census.No','Height','HtEst','Weib.parm.source')]

#Offset table

colnames(HtCns2)<- c('TreeID','PlotViewID','PlotID','CensusNoHtLast','HtLast','HtRefEstLast','HtRefEstTypeLast')
Toffset<- HtCns2
Toffset$OffsetH<-Toffset$HtLast-Toffset$HtRefEstLast
head(Toffset)

Heightsb<-merge(dw.all,Toffset,by=c('TreeID','PlotID','PlotViewID'),all=TRUE )

Heightsb$Ht1<- ifelse (is.na(Heightsb$HtLast), Heightsb$HtEst,
                         ifelse(Heightsb$CensusNoHtLast== Heightsb$Census.No, Heightsb$Height, Heightsb$HtEst+Heightsb$OffsetH)
                              )

Heightsb$Ht2 <- ifelse (is.na(Heightsb$HtLast), Heightsb$HtEst,
                       ifelse(Heightsb$CensusNoHtLast== Heightsb$Census.No, Heightsb$Height, 
                              Heightsb$HtEst+(Heightsb$OffsetH*(Heightsb$HtRefEstLast/Heightsb$HtEst)) )
                        )

return(Heightsb)
}
