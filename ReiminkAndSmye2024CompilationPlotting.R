


rm( list = ls( ) )


computer       <- "JesseReimink"
if( computer == "jxr1350" ) {
  google.drive.dir <- paste( "/Users/", computer, 
                             "/My Drive/", sep = "" )
} else {
  google.drive.dir <- paste( "/Users/", computer, 
                             "/My Drive (jxr1350@psu.edu)/", sep = "" )
}
#############################  SET THE DIRECTORIES #############################
source.dir  <- paste( google.drive.dir, 
                      "Research/PostDoc/R scripts/", sep = "" )
lib.dir  <- paste( google.drive.dir, 
                   "Research/PostDoc/R scripts/", sep = "" )
data.dir		<- paste( google.drive.dir, 
                    "Research/PSU/Projects/Neoarchean Sediment Geotherms/", sep = "" )

# gcdkit.dir		<- paste( google.drive.dir, 
#                       "Research/PSU/Projects/Reworking Through Time", sep = "" )
export.dir		<- paste( google.drive.dir,
                      "Research/PSU/Projects/Neoarchean Sediment Geotherms/", sep = "" )

setwd( data.dir )


library( ggplot2) 


## HPE backcalculations

## asi calculation

### ASI calculation using molar proportions of the major oxides
asi.index.simple  <- function( input.data ) {
  # input.data <- data.acasta.wr
  ### ASI calculation using molar concentrations of the major oxides
  
  cations <- data.frame( 
    Si.cation = input.data$SiO2 / 60.08,
    Ti.cation	= input.data$TiO2 / 79.866,
    Al.cation	= input.data$Al2O3 / 101.96,
    Fe.cation	= input.data$FeO / 71.844,
    Mn.cation	= input.data$MnO / 70.9374,
    Mg.cation	= input.data$MgO / 40.3,
    Ca.cation	= input.data$CaO / 56.0774,
    Na.cation	= input.data$Na2O / 61.9789,
    K.cation	= input.data$K2O / 94.2,
    P.cation	= input.data$P2O5 / 141.943 )
  
  sumcations	<- apply( cations, 1, sum, na.rm = T )
  
  cations <- cbind( cations, sumcations )
  
  Si.proportion	<- cations$Si.cation / sumcations
  Ti.proportion	<- cations$Ti.cation / sumcations
  Al.proportion	<- cations$Al.cation / sumcations
  Fe.proportion	<- cations$Fe.cation / sumcations
  Mn.proportion	<- cations$Mn.cation / sumcations
  Mg.proportion	<- cations$Mg.cation / sumcations
  Ca.proportion	<- cations$Ca.cation / sumcations
  Na.proportion	<- cations$Na.cation / sumcations
  K.proportion	<- cations$K.cation / sumcations
  P.proportion	<- cations$P.cation / sumcations
  
  ASI = cations$Al.cation / ( cations$Ca.cation + cations$Na.cation + cations$K.cation )
  ASI
}


## constants required
heat.prod.232 = 2.636817e-05 # Watt/kg
heat.prod.238 = 9.494597e-05 # W/kg
heat.prod.235 = 5.684024e-04  # W/kg

lambda_40Kbd = 4.9587E-10 # beta decay lambda /yr
lambda.232Th = 4.948e-11
lambda.40K = 5.5492e-10
Lambda238 = 1.55125e-10
Lambda235 = 9.8485e-10
u.isotope.composition = 137.818
heat.prod.40bd = ( 9.350e-14 * lambda_40Kbd ) / 31556926 / ( 40 * 1.660539040e-27 ) # W/kg

lambda_40Kec = 5.9540e-11 # electron capture lambda /yr
heat.prod.40ec = ( 2.330e-13 * lambda_40Kec ) / 31556926 / ( 40 * 1.660539040e-27 ) # W/kg
br = 0.8928 #
amu235 = 235
amu238 = 238


hpe.backcalc <- function( U.conc, Th.conc, K2O.conc, Age, Density ) {
  ## age in Ma
  ## Density in kg/m3 ie 3200
  ## convert K2O to K ppm
  K.conc = K2O.conc* 1e4 * 78.2 / 94.2 # K in ppm
  
  ## 232 Th calcs
  c.232 = Th.conc
  m232.per.kg = c.232 * Density / 1e6
  a.232 = m232.per.kg * heat.prod.232 * 1e6     # uW/m3
  a.time.232 = a.232 * exp( Age * 1e6 * lambda.232Th )
  
  ## 238 U calcs
  n_235 = U.conc/( u.isotope.composition * amu235 + amu235 ) # mol
  n_238 = ( u.isotope.composition * n_235 * amu235 ) / amu238 # mol
  c_238 = n_238 * amu238
  c_235 = n_235 * amu235
  m238perkg = c_238 * Density / 1E6  # kg
  a.238 = m238perkg * heat.prod.238 * 1e6 # uW/m3
  a.time.238 = a.238 * exp( Age * 1e6 * Lambda238 ) 
  
  ## 235 U calcs
  m235perkg = c_235 * Density / 1E6 # kg
  a.235 = m235perkg * heat.prod.235 * 1e6 # uW/m3
  a.time.235 = a.235 * exp( Age * 1e6 * Lambda235 ) 
  
  # 40K calcs
  c_40 = 1.1668e-04*K.conc          # ppm
  # beta-decay 40K -> 40Ca
  m40bdperkg = c_40 * br * Density / 1E6 # kg
  a.40bd = m40bdperkg * heat.prod.40bd * 1e6 # uW/m3
  a.time.40bd = a.40bd * exp( Age * 1e6 * lambda_40Kbd ) 
  # electron-capture 40K -> 40Ar
  m40ecperkg = c_40 * ( 1 - br ) * Density / 1e6 # kg
  a.40ec = m40ecperkg * heat.prod.40ec * 1e6 # uW/m3
  a.time.40ec = a.40ec * exp( Age * 1e6 * lambda_40Kec )
  a.time.total.40 = a.time.40bd + a.time.40ec
  
  # Total heat production
  a.time.232 + a.time.238 + a.time.235 + a.time.total.40
}

hpe.backcalc.K <- function( U.conc, Th.conc, K.conc, Age, Density ) {
  ## age in Ma
  ## Density in kg/m3 ie 3200
  ## convert K ppm from wt% for stanford compilation
  K.conc = K.conc* 1e4  # K in ppm
  
  ## 232 Th calcs
  c.232 = Th.conc
  m232.per.kg = c.232 * Density / 1e6
  a.232 = m232.per.kg * heat.prod.232 * 1e6    # uW/m3
  a.time.232 = a.232 * exp( Age * 1e6 * lambda.232Th )
  
  ## 238 U calcs
  n_235 = U.conc/( u.isotope.composition * amu235 + amu235 ) # mol
  n_238 = ( u.isotope.composition * n_235 * amu235 ) / amu238 # mol
  c_238 = n_238 * amu238
  c_235 = n_235 * amu235
  m238perkg = c_238 * Density / 1E6  # kg
  a.238 = m238perkg * heat.prod.238 * 1e6 # uW/m3
  a.time.238 = a.238 * exp( Age * 1e6 * Lambda238 ) 
  
  ## 235 U calcs
  m235perkg = c_235 * Density / 1E6 # kg
  a.235 = m235perkg * heat.prod.235 * 1e6 # uW/m3
  a.time.235 = a.235 * exp( Age * 1e6 * Lambda235 ) 
  
  # 40K calcs
  c_40 = 1.1668e-04*K.conc          # ppm
  # beta-decay 40K -> 40Ca
  m40bdperkg = c_40 * br * Density / 1E6 # kg
  a.40bd = m40bdperkg * heat.prod.40bd * 1e6 # uW/m3
  a.time.40bd = a.40bd * exp( Age * 1e6 * lambda_40Kbd ) 
  # electron-capture 40K -> 40Ar
  m40ecperkg = c_40 * ( 1 - br ) * Density / 1E6 # kg
  a.40ec = m40ecperkg * heat.prod.40ec * 1e6 # uW/m3
  a.time.40ec = a.40ec * exp( Age * 1e6 * lambda_40Kec )
  a.time.total.40 = a.time.40bd + a.time.40ec
  
  # Total heat production
  a.time.232 + a.time.238 + a.time.235 + a.time.total.40
}

hpe.40.backcalc <- function( K2O, Age, Density ) {
  ## convert K2O to K ppm
  K.conc = K2O* 1e4 * 78.2 / 94.2 # K in ppm
  # 40K calcs
  # 40K calcs
  c_40 = 1.1668e-04*K.conc          # ppm
  # beta-decay 40K -> 40Ca
  m40bdperkg = c_40 * br * Density / 1E6 # kg
  a.40bd = m40bdperkg * heat.prod.40bd * 1e6 # uW/m3
  a.time.40bd = a.40bd * exp( Age * 1e6 * lambda_40Kbd ) 
  # electron-capture 40K -> 40Ar
  m40ecperkg = c_40 * ( 1 - br ) * Density / 1E6 # kg
  a.40ec = m40ecperkg * heat.prod.40ec * 1e6 # uW/m3
  a.time.40ec = a.40ec * exp( Age * 1e6 * lambda_40Kec )
  a.time.total.40 = a.time.40bd + a.time.40ec
}

hpe.232.backcalc <- function( Th.conc, Age, Density ) {
  c.232 = Th.conc
  m232.per.kg = c.232 * Density / 1e6
  a.232 = m232.per.kg * heat.prod.232 * 1e6 #    # uW/m3
  
  a.232 * exp( Age * 1e6 * lambda.232Th )
}
hpe.235.backcalc <- function( U.conc, Age, Density ) {
  
  ## 235 U calcs
  n_235 = U.conc/( u.isotope.composition * amu235 + amu235 ) # mol
  n_238 = ( u.isotope.composition * n_235 * amu235 ) / amu238 # mol
  c_238 = n_238 * amu238
  c_235 = n_235 * amu235
  m235perkg = c_235 * Density / 1E6 # kg
  a.235 = m235perkg * heat.prod.235 * 1e6 # uW/m3
  
  a.235 * exp( Age * 1e6 * Lambda235 ) 
}

hpe.238.backcalc <- function( U.conc, Age, Density ) {
  ## 238 U calcs
  n_235 = U.conc/( u.isotope.composition * amu235 + amu235 ) # mol
  n_238 = ( u.isotope.composition * n_235 * amu235 ) / amu238 # mol
  c_238 = n_238 * amu238
  c_235 = n_235 * amu235
  m238perkg = c_238 * Density / 1E6  # kg
  a.238 = m238perkg * heat.prod.238 * 1e6 # uW/m3
  
  a.238 * exp( Age * 1e6 * Lambda238 ) 
}


## test HPE backcalculations and plot
data.test.hpe <- data.frame( age = seq( from = 1, to = 4500, by = 1 ) )
data.test.hpe$U = 2.7
data.test.hpe$Th = 10.5
data.test.hpe$K2O = 2.80 
data.test.hpe$hpe <- hpe.backcalc( data.test.hpe$U, data.test.hpe$Th, data.test.hpe$K2O,
                                   data.test.hpe$age, 2600 )
data.test.hpe$rel.hpe <- data.test.hpe$hpe / data.test.hpe$hpe[1]
data.test.hpe$U.hpe <- hpe.238.backcalc( data.test.hpe$U, data.test.hpe$age, 2600 ) + hpe.235.backcalc(data.test.hpe$U, data.test.hpe$age, 2600)
data.test.hpe$Th.hpe <- hpe.232.backcalc( data.test.hpe$Th, data.test.hpe$age, 2600 )
data.test.hpe$K.hpe <- hpe.40.backcalc( data.test.hpe$K2O, data.test.hpe$age, 2600 )

data.test.hpe$rel.U.hpe <- data.test.hpe$U.hpe / data.test.hpe$U.hpe[1]
data.test.hpe$rel.Th.hpe <- data.test.hpe$Th.hpe / data.test.hpe$Th.hpe[1]
data.test.hpe$rel.K.hpe <- data.test.hpe$K.hpe / data.test.hpe$K.hpe[1]


ggplot( data.test.hpe, aes( x = age, y = hpe ) ) +
  geom_line( linewidth = 1 ) +
  geom_line( aes( y = U.hpe ), color = viridis(9)[3] ) +
  geom_line( aes( y = Th.hpe ), color = viridis(9)[6] ) +
  geom_line( aes( y = K.hpe ), color = viridis(9)[9] ) +
  labs( x = expression( paste( "Age (Ma)" )), 
        y = expression( paste( "Heat Production (",mu,"W/m"^3,")" )) ) 





data.ttg                <- read.csv( "Moyen_2011_Lithos compilation.csv" )
# data.seds               <- read.csv( "Shale Compilation Bindeman et al 2016.csv" )


data.stypes             <- read.csv( "Bucholz and Spencer SData.csv" ) ## read Bucholz and Spencer data
data.slave.granites     <- read.csv( "Slave Granite compilation.csv" ) ## Davis Slave granite compilation
data.kaapvaal.granites  <- read.csv( "Laurent 2014 Pietersburg Block.csv" ) ## Laurent 2014 data from Pietersburg block in Kaapvaal
data.seds.chicago       <- read.csv( "Ptacek et al Sedimentary Compilation.csv" )

## read in data as separate lists from the excel sheets
data.amphibolite.xeno   <- read.csv( "Amphibolite xenoliths.csv" )
data.granulite.xeno     <- read.csv( "Granulite Xenoliths.csv" )
data.granulite.archean  <- read.csv( "Archean Granulite Terrains.csv" )
data.granulite.young    <- read.csv( "Postarchean Granulite Terrains.csv" )
data.hacker.meta        <- read.csv( "Hacker_xenolith_data.csv" )

## read in melt compositions from Andy's modeling
data.melts.7kbar        <- read.table( "meltcomp_sat7kbar", header = T, sep = "" )
data.melts.sat          <- read.table( "meltcomp_saturated", header = T, sep = "" )
data.melts.v2           <- read.table( "meltcomp_11_30_22", header = T, sep = "" )


## Read in Smye's PCA and heat production for the sedimentary rock samples
data.seds.pca            <- read.csv( "PCAseds.csv", header = T ) 
data.ign.pca            <- read.csv("PCAscores_Igneous.csv", header = T )


# convert to numeric
data.kaapvaal.granites$MnO   <- as.numeric( as.character( data.kaapvaal.granites$MnO ))

data.seds.chicago$Age   <- as.numeric( as.character( data.seds.chicago$Age ))
data.seds.chicago$K2O <- as.numeric( as.character( data.seds.chicago$K2O ))
data.seds.chicago$U <- as.numeric( as.character( data.seds.chicago$U ))
data.seds.chicago$Th <- as.numeric( as.character( data.seds.chicago$Th ))

data.stypes$Th <- as.numeric( as.character( data.stypes$Th ))

## Calculate HP 

data.ttg$heat.prod.2800Ma <- hpe.backcalc( data.ttg$U, data.ttg$Th, data.ttg$K2O,
                                           2800, 2600 )
data.slave.granites$heat.prod.2800Ma <- hpe.backcalc( data.slave.granites$U, data.slave.granites$Th, data.slave.granites$K2O,
                                                      2800, 2600 )
data.kaapvaal.granites$heat.prod.2800Ma <- hpe.backcalc( data.kaapvaal.granites$U, data.kaapvaal.granites$Th, data.kaapvaal.granites$K2O,
                                                         2800, 2600 )
data.stypes$heat.prod.2800Ma <- hpe.backcalc( data.stypes$U, data.stypes$Th, data.stypes$K2O,
                                              2800, 2600 )



data.seds.chicago$heat.prod.2800Ma <- hpe.backcalc( data.seds.chicago$U, data.seds.chicago$Th, data.seds.chicago$K2O,
                                                    2800, 2600 )
data.amphibolite.xeno$heat.prod.2800Ma <- hpe.backcalc( data.amphibolite.xeno$U, data.amphibolite.xeno$TH, data.amphibolite.xeno$K2O,
                                                        2800, 2600 )
data.granulite.xeno$heat.prod.2800Ma <- hpe.backcalc( data.granulite.xeno$U, data.granulite.xeno$TH, data.granulite.xeno$K2O,
                                                      2800, 2600 )
data.granulite.archean$heat.prod.2800Ma <- hpe.backcalc( data.granulite.archean$U, data.granulite.archean$TH, data.granulite.archean$K2O,
                                                         2800, 2600 )
data.granulite.young$heat.prod.2800Ma <- hpe.backcalc( data.granulite.young$U, data.granulite.young$TH, data.granulite.young$K2O,
                                                       2800, 2600 )
data.hacker.meta$heat.prod.2800Ma <- hpe.backcalc( data.hacker.meta$U, data.hacker.meta$Th, data.hacker.meta$K2O,
                                                   2800, 2600 )


## first filter Laurent data
data.kaapvaal.granites$FeO <- data.kaapvaal.granites$Fe2O3 * 0.8998
data.kaapvaal.granites$asi <- asi.index.simple( data.kaapvaal.granites )
data.kaapvaal.potassic.granites <- subset( data.kaapvaal.granites, Classification %in% c("Bt-(Ms-)granite", "Hybrid" ) )

## filter slave data 
data.slave.granites$asi <- asi.index.simple( data.slave.granites )
data.slave.granites.peraluminous <- subset( data.slave.granites, suite %in% c("Contwoyto", "Prosperous", "Two-mica granite", "Two-mica Granodiorite", "Gt two-mica granite" ) )
data.slave.granites.potassic <- subset( data.slave.granites, suite %in% c( "Defeat", "Yamba", "Morose", "Concession", "Awry", "Siege", "Ghost", "Gondor", "Stagg", "Wishbone", "Olga"  ) )


## combine Kaapvaal granites and Slave potassic granites 
data.slave.kaapvaal.potassic.granites.1 <- data.slave.granites.potassic[ , c( 'suite', "SiO2", "K2O", "U", "Th", "heat.prod.2800Ma", "asi" ) ]
data.slave.kaapvaal.potassic.granites.2 <- data.kaapvaal.potassic.granites[ , c( "Classification", "SiO2", "K2O", "U", "Th", "heat.prod.2800Ma", "asi" ) ]
colnames(data.slave.kaapvaal.potassic.granites.2) <- c( 'suite', "SiO2", "K2O", "U", "Th", "heat.prod.2800Ma", "asi" )

data.slave.kaapvaal.potassic.granites <- rbind( data.slave.kaapvaal.potassic.granites.1, data.slave.kaapvaal.potassic.granites.2 )

## get archean samples
data.ttg.archean        <- subset( data.ttg, Age > 2800 )
data.seds.chicago.archean <- subset( data.seds.chicago, Age > 2500 )
data.stypes.archean         <- subset( data.stypes, Archean..Proterozoic == "A" )



## get sediments from Hacker xenolith metamorphic databases
data.hacker.seds.archean <- subset( data.hacker.meta, class == "metased" )
data.hacker.ign.archean  <- subset( data.hacker.meta, class != "metased" )


###### You can use the data above to replicate the values found in this PCA plotting below



############################# Use PCA space to categorize sediments
## calculate and assign proportions to the sample database
shale.loci = c( -10, 18 )
silic.loci = c( 38, -7)
mafic.loci = c( -30, -20 ) 
avg.loci = c( 1, 2 ) 

data.seds.pca$dist.shale <- sqrt( ( shale.loci[1] - data.seds.pca$PC1 ) ^2 + ( shale.loci[2] - data.seds.pca$PC2 ) ^2 )
data.seds.pca$dist.silic <- sqrt( ( silic.loci[1] - data.seds.pca$PC1 ) ^2 + ( silic.loci[2] - data.seds.pca$PC2 ) ^2 ) 
data.seds.pca$dist.mafic <- sqrt( ( mafic.loci[1] - data.seds.pca$PC1 ) ^2 + ( mafic.loci[2] - data.seds.pca$PC2 ) ^2 ) 
data.seds.pca$dist.avg <- sqrt( ( avg.loci[1] - data.seds.pca$PC1 ) ^2 + ( avg.loci[2] - data.seds.pca$PC2 ) ^2 ) 


## calculate the distance from a line going from the average out through the shale point. If the data point lies on the opposite side of a line
#   orthogonal to the 'shale line' and intersecting at the average point, we will use the distance from the average as the 'shale distance'. This prevents
#   data on the far end of the projected line from being calculated as close to the shale arm of the data
shale.m <- ( avg.loci[2] - shale.loci[2] ) / ( avg.loci[1] - shale.loci[1] )
shale.b <- shale.loci[2] - shale.m * shale.loci[1]
shale.a  <- - shale.m
shale.b2 <- 1
shale.c  <- - shale.b
data.seds.pca$dist.shale.line <- abs( ( shale.a * data.seds.pca$PC1 + shale.b2 * data.seds.pca$PC2 + shale.c ) / sqrt( shale.a ^ 2 + shale.b2 ^ 2 ) )
data.seds.pca$shale.line.diff <- data.seds.pca$PC2 - ( ( -1/shale.m ) * data.seds.pca$PC1 + ( avg.loci[ 2 ] - ( -1/shale.m ) * avg.loci[ 1 ] ) )
data.seds.pca$shale.line.diff.class <- ifelse( data.seds.pca$shale.line.diff > 0, 1,
                                               ifelse( data.seds.pca$shale.line.diff < 0, -1, 0 ) )
data.seds.pca$dist.shale.use <- ifelse( data.seds.pca$shale.line.diff > 0, data.seds.pca$dist.shale.line,
                                        ifelse( data.seds.pca$shale.line.diff < 0, data.seds.pca$dist.avg, data.seds.pca$dist.avg ) )

silic.m <- ( avg.loci[2] - silic.loci[2] ) / ( avg.loci[1] - silic.loci[1] )
silic.b <- silic.loci[2] - silic.m * silic.loci[1]
silic.a  <- - silic.m
silic.b2 <- 1
silic.c  <- - silic.b
data.seds.pca$dist.silic.line <- abs( ( silic.a * data.seds.pca$PC1 + silic.b2 * data.seds.pca$PC2 + silic.c ) / sqrt( silic.a ^ 2 + silic.b2 ^ 2 ) )
data.seds.pca$silic.line.diff <- data.seds.pca$PC2 - ( ( -1/silic.m ) * data.seds.pca$PC1 + ( avg.loci[ 2 ] - ( -1/silic.m ) * avg.loci[ 1 ] ) )
data.seds.pca$silic.line.diff.class <- ifelse( data.seds.pca$silic.line.diff > 0, 1,
                                               ifelse( data.seds.pca$silic.line.diff < 0, -1, 0 ) )
data.seds.pca$dist.silic.use <- ifelse( data.seds.pca$silic.line.diff < 0, data.seds.pca$dist.silic.line,
                                        ifelse( data.seds.pca$silic.line.diff > 0, data.seds.pca$dist.avg, data.seds.pca$dist.avg ) )

mafic.m <- ( avg.loci[2] - mafic.loci[2] ) / ( avg.loci[1] - mafic.loci[1] )
mafic.b <- mafic.loci[2] - mafic.m * mafic.loci[1]
mafic.a  <- - mafic.m
mafic.b2 <- 1
mafic.c  <- - mafic.b
data.seds.pca$dist.mafic.line <- abs( ( mafic.a * data.seds.pca$PC1 + mafic.b2 * data.seds.pca$PC2 + mafic.c ) / sqrt( mafic.a ^ 2 + mafic.b2 ^ 2 ) )
data.seds.pca$mafic.line.diff <- data.seds.pca$PC2 - ( ( -1/mafic.m ) * data.seds.pca$PC1 + ( avg.loci[ 2 ] - ( -1/mafic.m ) * avg.loci[ 1 ] ) )
data.seds.pca$mafic.line.diff.class <- ifelse( data.seds.pca$mafic.line.diff > 0, 1,
                                               ifelse( data.seds.pca$mafic.line.diff < 0, -1, 0 ) )
data.seds.pca$dist.mafic.use <- ifelse( data.seds.pca$mafic.line.diff < 0, data.seds.pca$dist.mafic.line,
                                        ifelse( data.seds.pca$mafic.line.diff > 0, data.seds.pca$dist.avg, data.seds.pca$dist.avg ) )




data.seds.pca$sum.dist <- data.seds.pca$dist.shale.use + data.seds.pca$dist.silic.use + data.seds.pca$dist.mafic.use
data.seds.pca$prop.shale <- 1 - ( data.seds.pca$dist.shale.use / data.seds.pca$sum.dist )
data.seds.pca$prop.silic <- 1 - ( data.seds.pca$dist.silic.use / data.seds.pca$sum.dist )
data.seds.pca$prop.mafic <- 1 - ( data.seds.pca$dist.mafic.use / data.seds.pca$sum.dist )
data.seds.pca$prop.cat <- ifelse( data.seds.pca$prop.shale > 0.65, "Shale",
                                  ifelse( data.seds.pca$prop.silic > 0.7, "Siliciclastic", "Mafic" ) )




############################ Make plots using the PCA fields
data.ign.pca$heat.prod.2800Ma <- data.ign.pca$A
data.seds.pca$heat.prod.2800Ma <- data.seds.pca$A
data.ign.pca.lowsio2 <- subset( data.ign.pca, SiO2.wt. <= 55 )
data.ign.pca.highsio2 <- subset( data.ign.pca, SiO2.wt. > 55 )
data.seds.pca.shales <- subset( data.seds.pca, prop.cat == "Shale" )
data.seds.pca.silic <- subset( data.seds.pca, prop.cat == "Siliciclastic" )
data.seds.pca.mafic <- subset( data.seds.pca, prop.cat == "Mafic" )


data.heatprod.summary <- data.ign.pca[ , c( 1, 6 ) ]
data.heatprod.summary <- rbind( data.heatprod.summary, data.ign.pca.highsio2[ , c( 1, 6 ) ] )
data.heatprod.summary <- rbind( data.heatprod.summary, data.ign.pca.lowsio2[ , c( 1, 6 ) ] )
data.heatprod.summary <- rbind( data.heatprod.summary, data.seds.pca.silic[ , c( 1, 7 ) ] )
data.heatprod.summary <- rbind( data.heatprod.summary, data.seds.pca.shales[ , c( 1, 7 ) ] )
data.heatprod.summary <- rbind( data.heatprod.summary, data.seds.pca.mafic[ , c( 1, 7 ) ] )
data.heatprod.summary <- bind_rows( data.heatprod.summary, data.stypes.archean[ , c( 3, 84 ) ] )
data.heatprod.summary <- bind_rows( data.heatprod.summary, data.slave.kaapvaal.potassic.granites[ , c( 1, 6 ) ] )
data.heatprod.summary <- bind_rows( data.heatprod.summary, data.hacker.seds.archean[ , c( 1, 159 ) ] )
data.heatprod.summary <- bind_rows( data.heatprod.summary, data.granulite.archean[ , c( 1, 108 ) ] )

data.heatprod.summary$category <- c( rep("Archean Igneous", nrow( data.ign.pca ) ),
                                     rep("Archean TTG", nrow( data.ign.pca.highsio2 ) ),
                                     rep("Archean Mafic", nrow( data.ign.pca.lowsio2 ) ),
                                     rep("Archean Siliciclastic", nrow( data.seds.pca.silic ) ),
                                     rep("Archean Shale", nrow( data.seds.pca.shales ) ),
                                     rep("Archean Mafic Seds", nrow( data.seds.pca.mafic ) ),
                                     rep("Archean Peraluminous", nrow( data.stypes.archean ) ),
                                     rep("Archean Metaluminous", nrow( data.slave.kaapvaal.potassic.granites ) ),
                                     rep("Archean Xenoliths", nrow( data.hacker.seds.archean ) ),
                                     rep("Archean Granulites", nrow( data.granulite.archean ) ) )


#### Figure 2 base plot
ggplot( data.heatprod.summary, aes( x = heat.prod.2800Ma, y = fct_rev( fct_inorder( category ) ), 
                                    fill = fct_rev( fct_inorder( category ) ) ) )  +
  xlim( 0, 10 ) +
  labs( x = expression( paste( "Heat Production (",mu,"W/m"^3,")" ) ) ) +
  theme( axis.text.y = element_text( size = 6 ) ) +
  scale_fill_viridis( discrete = T ) +
  # geom_point()
  stat_density_ridges( quantile_lines = TRUE, scale = 3) 







