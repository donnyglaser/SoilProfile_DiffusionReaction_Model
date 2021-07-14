## please reference readMe file prior to running this code ##
## Designed by Donald Glaser ##
## Please cite Glaser et al. 2021. Astrobiology. DOI: ##

library(ggplot2)
library(plyr)

## PARAMETERS ##
## Before you start, double check these parameters ##
inpFile <- ## input file name here ##
site <- ## Abbreviated site name ##
siteFull <- ## full site name ##
thRsq <- 0.9
thRXN <- 0.1
Por <- ## soil porosity ##
mVer <- 2.3 ## model version ##

inp <- read.csv2(file = inpFile, header=TRUE, sep=",", stringsAsFactors = FALSE)
inp$Date_Time <- as.POSIXct(inp$Date_Time, format = "%Y-%m-%d %H:%M:%S")
for(i in 2:12) {
  inp[,i] <- as.numeric(inp[,i])
}

maxOBS <- max(inp$OBS)
allData <- data.frame()

## Begin section 1 ##
## first step to identify profile shapes (H2O minima(adsorbing), H2O maxima(desorbing), or none) ##
## calculates the AH v depth slope between i) surface and 2.5cm (LSlp1_2) and ii) 2.5cm and 30cm (LSlp2_6)

for(i in 1:maxOBS) {
  sub <- data.frame()
  temp <- data.frame()
  sub <- subset(inp, OBS == i)
  time <- sub$Date_Time[1]
  sub1_2 <- sub[-3:-6,]
  sub2_6 <- sub[-1,]
  minAH <- min(sub$MeanAH)
  minDep <- sub$Depth[which.min(sub$MeanAH)]
  maxAH <- max(sub$MeanAH)
  maxDep <- sub$Depth[which.max(sub$MeanAH)]
  rng <- maxAH - minAH
  Lreg <- lm(MeanAH ~ Depth, data = sub1_2)
  LSlp1_2 <- Lreg$coefficients[2]
  Lreg <- lm(MeanAH ~ Depth, data = sub2_6)
  LSlp2_6 <- Lreg$coefficients[2]
  Lght <- sub$MeanLight[1]
  temp <- cbind(time, i, Lght, minAH, minDep, maxAH, maxDep, rng, LSlp1_2, LSlp2_6)
  allData <- rbind(allData, temp)
}
allData$time <- as.POSIXct(allData$time, origin = "1970-01-01")
row.names(allData) <- c(1:maxOBS)

## End section 1 ##

## Begin section 2 ##
## Finish identifying profiles with adsorbing or desorbing shape ##
## Identify which lower point to use as boundry for regression (no flux boundry)##
## Calculate two-part linear approximation of upper and lower profiles and depth of reaction point ##

AdsID <- data.frame()
AdsID[1,1] <- "NA"
for(i in 2:maxOBS) {
  temp <- data.frame()
  LSlp1_2 <- allData$LSlp1_2[i]
  maxDep <- allData$maxDep[i]
  minDep <- allData$minDep[i]
  LSlp2_6 <- allData$LSlp2_6[i]
  if((LSlp1_2 > 15) & (LSlp2_6 < 0)) { ##DESORPTION##
    Flg1 <- "Des"
  }
  else if((LSlp1_2 < -15) & (LSlp2_6 > 0)) {
    Flg1 <- "Ads"
  }
  else {
    Flg1 <- "NA"
  }
  AdsID <- rbind(AdsID, Flg1)
}
allData <- cbind(allData, AdsID)
names(allData)[names(allData) == "V1"] <- "Flg1"


profiles <- data.frame()
for(i in 1:maxOBS) { ## Calculate upper and lower AH & T regressions for each adsorbing and desorbing profile ##
  out <- data.frame()
  sub <- subset(inp, OBS == i)
  Flg1 <- allData$Flg1[i]
  RegT <- summary(lm(MeanT ~ Depth, data = sub))
  T.m <- as.numeric(RegT$coefficients[2,1])
  T.b <- as.numeric(RegT$coefficients[1,1])
  T.fun <- function(x) (T.m * x) + T.b
  if(Flg1 %in% "NA") {
    out[1,1:9] <- "NA"
    out[1,10] <- T.m
    out[1,11] <- T.b
    out[1,12] <- RegT$r.squared
    out[1,13] <- 6
    out[1,14:17] <- "NA"
  }

  else {
    upper <- subset(sub, Depth <= 20)
    Lrsq <- summary(lm(MeanAH ~ Depth, data = upper))
    Lrsq <- Lrsq$r.squared
    if(Lrsq < thRsq) {
      nUp <- nrow(upper)
      for(j in 1:(nUp - 1)) {
        upper <- upper[-nrow(upper),]
        Lrsq <- summary(lm(MeanAH ~ Depth, data = upper))
        Lrsq <- Lrsq$r.squared
        if(Lrsq > thRsq) {
          break
        }
      }
    }
    botUp <- max(upper$Depth)
    upRegAH <- summary(lm(MeanAH ~ Depth, data = upper))
    AH.m.up <- as.numeric(upRegAH$coefficients[2,1])
    out[1,1] <- AH.m.up
    AH.b.up <- as.numeric(upRegAH$coefficients[1,1])
    out[1,2] <- AH.b.up
    out[1,3] <- upRegAH$r.squared
    out[1,4] <- nrow(upper)

    if(Flg1 %in% "Des"){
      lower <- sub
      lower5_6 <- lower[-1:-4,]
      lower4_5 <- subset(sub, Depth >= 10 & Depth < 30)
      m5_6 <- summary(lm(MeanAH ~ Depth, data = lower5_6))
      m5_6 <- as.numeric(m5_6$coefficients[2,1])
      m4_5 <- summary(lm(MeanAH ~ Depth, data = lower4_5))
      m4_5 <- as.numeric(m4_5$coefficients[2,1])
      if(m5_6 > 0) {
        if(m4_5 > 0) {
          lower <- sub[-5:-6,]
          d.lo <- 4
          botLo <- 10
        }else{
          lower <- sub[-6,]
          d.lo <- 5
          botLo <- 20
        }
      }else {
        lower <- sub
      }
      nLo <- nrow(lower)
      t.lo <- data.frame()
      t.lo1 <- data.frame()
      lower1 <- lower
      lower2 <- lower[-nrow(lower),]

      for(j in 1:(nLo - 2)) {
        Lrsq <- summary(lm(MeanAH ~ Depth, data = lower1))
        Lrsq <- Lrsq$r.squared
        t.lo[j,1] <- Lrsq
        t.lo[j,2] <- min(lower1$Depth)
        t.lo[j,3] <- max(lower1$Depth)
        lower1 <- lower1[-1,]
      }
      nLo <- nrow(lower2)
      for(j in 1:(nLo - 2)) {
        Lrsq <- summary(lm(MeanAH ~ Depth, data = lower2))
        Lrsq <- Lrsq$r.squared
        t.lo1[j,1] <- Lrsq
        t.lo1[j,2] <- min(lower2$Depth)
        t.lo1[j,3] <- max(lower2$Depth)
        lower2 <- lower2[-1,]
      }
      t.lo <- rbind(t.lo, t.lo1)
      t.lo <- t.lo[order(-t.lo[,1]),]
      max.rsq <- t.lo[1,1]
      lower <- subset(sub, Depth >= t.lo[1,2] & Depth <= t.lo[1,3])
      AH.lo.reg <- summary(lm(MeanAH ~ Depth, data = lower))
      AH.m.lo <- as.numeric(AH.lo.reg$coefficients[2,1])
      AH.botUp <- upper$MeanAH[which(upper$Depth == botUp)]
      AH.botUp <- AH.botUp - (AH.botUp * .1)
      AH.lo.botUp <- (AH.m.lo * botUp) + AH.b.lo
      if(max.rsq < thRsq | AH.m.lo > 0 | AH.lo.botUp < AH.botUp) {
        lower <- subset(sub, Depth >= botUp)
        lower1 <- lower[-1,]
        if(nrow(lower1) > 2) { ## lower1 is the 'plateau' slope
          for(j in 1:(nrow(lower1) - 2)) {
            lower1 <- lower1[-nrow(lower1),]
          }
        }
        lower2 <- lower[-nrow(lower),]
        if(nrow(lower2) > 2) {
          for(j in 1:(nrow(lower2) - 2)) {
            lower2 <- lower2[-nrow(lower2),]
          }
        }
        lreg <- summary(lm(MeanAH ~ Depth, data = lower1))
        AH.m.lo <- lreg$coefficients[2,1]
        AH.b.lo <- lreg$coefficients[1,1]
        tempAH.lo <- (AH.m.lo * botUp) + AH.b.lo
        AH.botUp <- sub$MeanAH[which(sub$Depth == botUp)]
        AH.botUp <- AH.botUp - (AH.botUp * thRXN)
        if(tempAH.lo < AH.botUp) {
          lower <- lower2
        }else {
          lower <- lower1
        }
        lreg <- summary(lm(MeanAH ~ Depth, data = lower))
        AH.m.lo <- lreg$coefficients[2,1]
        if(AH.m.lo > 0) {
          lower <- subset(sub, Depth >= botUp)
          lower <- lower[-1:-2,]
          if(nrow(lower) > 2) {
            for(j in 1:(nrow(lower) - 2)) {
              lower <- lower[-nrow(lower),]
            }
          }
        }
      }
    }

    else { ## If Flg1 %in% "Ads"
      lower <- sub
      lower5_6 <- lower[-1:-4,]
      lower4_5 <- subset(sub, Depth >= 10 & Depth < 30)
      m5_6 <- summary(lm(MeanAH ~ Depth, data = lower5_6))
      m5_6 <- as.numeric(m5_6$coefficients[2,1])
      m4_5 <- summary(lm(MeanAH ~ Depth, data = lower4_5))
      m4_5 <- as.numeric(m4_5$coefficients[2,1])
      if(m5_6 < 4) {
        if(m4_5 < 4) {
          lower <- sub[-5:-6,]
          d.lo <- 4
          botLo <- 10
        }else{
          lower <- sub[-6,]
          d.lo <- 5
          botLo <- 20
        }
      }else {
        lower <- sub
      }
      nLo <- nrow(lower)
      t.lo <- data.frame()
      t.lo1 <- data.frame()
      lower1 <- lower
      lower2 <- lower[-nrow(lower),]

      for(j in 1:(nLo - 2)) {
        Lrsq <- summary(lm(MeanAH ~ Depth, data = lower1))
        Lrsq <- Lrsq$r.squared
        t.lo[j,1] <- Lrsq
        t.lo[j,2] <- min(lower1$Depth)
        t.lo[j,3] <- max(lower1$Depth)
        lower1 <- lower1[-1,]
      }
      nLo <- nrow(lower2)
      for(j in 1:(nLo - 2)) {
        Lrsq <- summary(lm(MeanAH ~ Depth, data = lower2))
        Lrsq <- Lrsq$r.squared
        t.lo1[j,1] <- Lrsq
        t.lo1[j,2] <- min(lower2$Depth)
        t.lo1[j,3] <- max(lower2$Depth)
        lower2 <- lower2[-1,]
      }
      t.lo <- rbind(t.lo, t.lo1)
      t.lo <- t.lo[order(-t.lo[,1]),]
      max.rsq <- t.lo[1,1]
      lower <- subset(sub, Depth >= t.lo[1,2] & Depth <= t.lo[1,3])
      AH.lo.reg <- summary(lm(MeanAH ~ Depth, data = lower))
      AH.m.lo <- as.numeric(AH.lo.reg$coefficients[2,1])
      AH.b.lo <- as.numeric(AH.lo.reg$coefficients[1,1])
      AH.botUp <- upper$MeanAH[which(upper$Depth == botUp)]
      AH.botUp <- AH.botUp + (AH.botUp * .1)
      AH.lo.botUp <- (AH.m.lo * botUp) + AH.b.lo
      if(max.rsq < thRsq | AH.m.lo < 0 | AH.lo.botUp > AH.botUp) {
        lower <- subset(sub, Depth >= botUp)
        lower1 <- lower[-1,]
        if(nrow(lower1) > 2) {
          for(j in 1:(nrow(lower1) - 2)) {
            lower1 <- lower1[-nrow(lower1),]
          }
        }
        lower2 <- lower[-nrow(lower),]
        if(nrow(lower2) > 2) {
          for(j in 1:(nrow(lower2) - 2)) {
            lower2 <- lower2[-nrow(lower2),]
          }
        }
        lreg <- summary(lm(MeanAH ~ Depth, data = lower1))
        AH.m.lo <- lreg$coefficients[2,1]
        AH.b.lo <- lreg$coefficients[1,1]
        tempAH.lo <- (AH.m.lo * botUp) + AH.b.lo
        AH.botUp <- sub$MeanAH[which(sub$Depth == botUp)]
        AH.botUp <- AH.botUp + (AH.botUp * thRXN)
        if(tempAH.lo > AH.botUp) {
          lower <- lower2
        }
        else {
          lower <- lower1
        }
        lreg <- summary(lm(MeanAH ~ Depth, data = lower))
        AH.m.lo <- lreg$coefficients[2,1]
        if(AH.m.lo < 0) {
          lower <- subset(sub, Depth >= botUp)
          lower <- lower[-1:-2,]
          if(nrow(lower) > 2) {
            for(j in 1:(nrow(lower) - 2)) {
              lower <- lower[-nrow(lower),]
            }
          }
        }
      }
    }

    loRegAH <- summary(lm(MeanAH ~ Depth, data = lower))
    AH.m.lo <- as.numeric(loRegAH$coefficients[2,1])
    out[1,5] <- AH.m.lo
    AH.b.lo <- as.numeric(loRegAH$coefficients[1,1])
    AH.fun.lo <- function(x) (AH.m.lo * x) + AH.b.lo
    out[1,6] <- AH.b.lo
    out[1,7] <- loRegAH$r.squared
    out[1,8] <- nrow(lower)
    bottom <- lower$Depth[nrow(lower)]
    out[1,9] <- bottom
    rxnDep <- as.numeric((AH.b.lo - AH.b.up) / (AH.m.up - AH.m.lo))
    rxnAH <- AH.fun.lo(rxnDep)
    rxnDep <- as.numeric((AH.b.lo - AH.b.up) / (AH.m.up - AH.m.lo))
    upMidDep <- rxnDep / 2
    loMidDep <- ((bottom - rxnDep) / 2) + rxnDep
    upIdx <- as.numeric(which.min(abs(sub$Depth-upMidDep)))
    loIdx <- as.numeric(which.min(abs(sub$Depth-loMidDep)))

    if(upMidDep > as.numeric(sub$Depth[upIdx])) {
      upIdx2 <- upIdx + 1
    }
    else {
      upIdx2 <- upIdx
      upIdx <- upIdx - 1
    }
    if(loMidDep > as.numeric(sub$Depth[loIdx])) {
      loIdx2 <- loIdx + 1
    }
    else {
      loIdx2 <- loIdx
      loIdx <- loIdx - 1
    }
    T.up.sub <- subset(sub, Depth >= sub$Depth[upIdx] & Depth <= sub$Depth[upIdx2])
    T.up.reg <- lm(MeanT ~ Depth, data = T.up.sub)
    T.up.fun <- function(x) (as.numeric(T.up.reg$coefficients[2]) * x) + as.numeric(T.up.reg$coefficients[1])
    T.up <- T.up.fun(upMidDep)
    T.lo.sub <- subset(sub, Depth >= sub$Depth[loIdx] & Depth <= sub$Depth[loIdx2])
    T.lo.reg <- lm(MeanT ~ Depth, data = T.lo.sub)
    T.lo.fun <- function(x) (as.numeric(T.lo.reg$coefficients[2]) * x) + as.numeric(T.lo.reg$coefficients[1])
    T.lo <- T.lo.fun(loMidDep)
    RegT <- summary(lm(MeanT ~ Depth, data = sub))
    out[1,10] <- RegT$coefficients[2,1]
    out[1,11] <- RegT$coefficients[1,1]
    out[1,12] <- RegT$r.squared
    out[1,13] <- 6
    out[1,14] <- rxnDep
    out[1,15] <- rxnAH
    out[1,16] <- T.up
    out[1,17] <- T.lo
  }
  profiles <- rbind(profiles, out)
}
colnames(profiles) <- c("UpperAHSlope", "UpperAHIntercept", "UpperAHRSq", "UpperAHN", "LowerAHSlope", "LowerAHIntercept", "LowerAHRSq", "LowerAHN", "LowerDepth", "TSlope", "TIntercept", "TRSq", "TN", "rxnDepth", "rxnAH", "T.up", "T.lo")

allData <- cbind(allData, profiles)
names(allData)[names(allData) == "i"] <- "OBS"

## End section 2 ##

## Begin section 3 ##
## this section i) labels diurnal cycle, ii) calculates upper and lower flux, iii) calculates net flux (at reaction point), ##
## iv) calculates reaction and WVA flux, and v) cumulative WVA over diurnal cycle ##

modOut <- matrix()
WVAFlg <- 0
for(i in 1:(maxOBS - 1)) { ## This loop labels the diurnal cycle ##
  if(allData$Flg1[i] %in% "Des" & allData$Flg1[i+1] %in% "NA") {
    modOut[i] <- WVAFlg
    WVAFlg <- WVAFlg + 1
  }
  else {
    modOut[i] <- WVAFlg
  }
}
modOut[maxOBS] <- WVAFlg
porTab <- matrix()
porTab[1:maxOBS] <- Por

DJTab <- data.frame()
SCoeff <- ((Por - 0.1)/0.9)^2
for(i in 1:maxOBS) {
  if(allData$Flg1[i] %in% "NA" | i == maxOBS) {
    temp <- cbind("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
    colnames(temp) <- c("D.up", "D.lo", "J.up", "J.lo", "NetJ", "RXN", "WVAFlux", "WVA")
    DJTab <- rbind(DJTab, temp)
  }
  else {
    temp <- data.frame()
    T.up <- as.numeric(allData$T.up[i]) + 273.15
    T.lo <- as.numeric(allData$T.lo[i]) + 273.15
    delAH.up <- as.numeric(allData$UpperAHSlope[i]) / 1000
    delAH.lo <- as.numeric(allData$LowerAHSlope[i]) / 1000
    delTime <- allData$time[i+1] - allData$time[i]
    delTime <- as.numeric(delTime) * 60
    D.up <- (T.up/273)^1.75
    D.up <- 0.073 * SCoeff * D.up
    D.lo <- (T.lo/273)^1.75
    D.lo <- 0.073 * SCoeff * D.lo
    J.up <- D.up * delAH.up
    J.lo <- D.lo * delAH.lo
    NetJ <- (-1 * J.lo) - J.up
    RXN <- -1 * NetJ
    WVAFlux <- NetJ
    WVA <- WVAFlux * delTime
    temp <- cbind(D.up, D.lo, J.up, J.lo, NetJ, RXN, WVAFlux, WVA)
    colnames(temp) <- c("D.up", "D.lo", "J.up", "J.lo", "NetJ", "RXN", "WVAFlux", "WVA")
    DJTab <- rbind(DJTab, temp)
  }
}
allData <- cbind(allData, porTab, modOut, DJTab)

cumlTab <- data.frame()
dayWVA <- 0
cumlWVA <- 0

for(i in 2:(maxOBS - 1)) {
  dayID <- as.numeric(allData$modOut[i])
  dayID.lo <- as.numeric(allData$modOut[i-1])
  dayID.up <- as.numeric(allData$modOut[i+1])
  Flg1 <- allData$Flg1[i]
  tWVA <- as.numeric(allData$WVA[i])
  tNetJ <- as.numeric(allData$NetJ[i])
  if(dayID == 0) {
    dayWVA <- dayWVA
    cumlWVA <- cumlWVA
    temp <- cbind(dayWVA, cumlWVA)
    cumlTab <- rbind(cumlTab, temp)
  } else if(dayID > dayID.lo){
    dayWVA <- 0
    if(Flg1 %in% "Ads" & tNetJ > 0) {
        dayWVA <- tWVA
        cumlWVA <- cumlWVA + tWVA
        temp <- cbind(dayWVA, cumlWVA)
        cumlTab <- rbind(cumlTab, temp)
    } else {
        cumlWVA <- cumlWVA
        temp <- cbind(dayWVA, cumlWVA)
        cumlTab <- rbind(cumlTab, temp)
    }
  } else {
      if((Flg1 %in% "Ads") & (tNetJ > 0)) {
        dayWVA <- (tWVA + dayWVA)
        cumlWVA <- (cumlWVA + tWVA)
        temp <- cbind(dayWVA, cumlWVA)
        cumlTab <- rbind(cumlTab, temp)
      } else if((Flg1 %in% "Des") & (tNetJ < 0)) {
        dayWVA <- (tWVA + dayWVA)
        cumlWVA <- (cumlWVA + tWVA)
        temp <- cbind(dayWVA, cumlWVA)
        cumlTab <- rbind(cumlTab, temp)
      } else {
        dayWVA <- dayWVA
        cumlWVA <- cumlWVA
        temp <- cbind(dayWVA, cumlWVA)
        cumlTab <- rbind(cumlTab, temp)
      }
  }
}

yr <- format(Sys.time(), "%y")
mn <- format(Sys.time(), "%m")
dy <- format(Sys.time(), "%d")
nmStr <- paste0(site, "_ProfileProperties_", yr, mn, dy, ".csv")
nmStr2 <- paste0(site, "_ProfileProperties_", yr, mn, dy, ".Rda")

cumlTab <- rbind(c(0,0), cumlTab)
cumlLen <- nrow(cumlTab)
cumlTab <- rbind(cumlTab, c(cumlTab[cumlLen,1], cumlTab[cumlLen,2]))
allData <- cbind(allData, cumlTab)
save(allData, nmStr2)

lab <- data.frame()
lab[1:maxOBS,1] <- ""
lab[1:6,1] <- c("Rsq Threshold", thRsq, "Porosity", Por, "Model Version", mVer)
colnames(lab) <- "Parameters"
allData <- cbind(lab, allData)

write.csv(allData, nmStr, row.names = FALSE)

## End section 3 ##

## Begin section 4 ##
## Plot the AH v depth profiles with and without the modeled linear approximations and T v depth profiles ##
## this section is optional ##

maxAH <- max(inp$MeanAH)
maxAH <- round_any(maxAH,100)
maxT <- max(inp$MeanT)
maxT <- round_any(maxT, 5)
maxOBS <- max(inp$OBS)
dir <- getwd()
dir.create("AH_Reg_Plots")
dir.create("AH_Plots")
dir.create("T_Reg_Plots")
dir.create("T_Plots")

for(i in 1:maxOBS) { 
  sub <- subset(inp, OBS == i)
  sub2 <- subset(allData, OBS == i)
  temp <- data.frame()
  yr <- format(sub$Date_Time[1], "%y")
  mon <- format(sub$Date_Time[1], "%m")
  dy <- format(sub$Date_Time[1], "%d")
  hr <- format(sub$Date_Time[1], "%H")
  min <- format(sub$Date_Time[1], "%M")
  timSt <- paste0(i, "-", yr, mon, dy, "_", hr, "-", min)
  if(sub2$Flg1 %in% "NA") {
    T.m <- as.numeric(sub2$TSlope)
    T.b <- as.numeric(sub2$TIntercept)
    T.fun <- function(x) (T.m * x) + T.b
  }
  else {
    T.m <- as.numeric(sub2$TSlope)
    T.b <- as.numeric(sub2$TIntercept)
    AH.m.up <- as.numeric(sub2$UpperAHSlope)
    AH.b.up <- as.numeric(sub2$UpperAHIntercept)
    AH.m.lo <- as.numeric(sub2$LowerAHSlope)
    AH.b.lo <- as.numeric(sub2$LowerAHIntercept)
    T.fun <- function(x) (T.m * x) + T.b
    AH.fun.up <- function(x) (AH.m.up * x) + AH.b.up
    AH.fun.lo <- function(x) (AH.m.lo * x) + AH.b.lo
  }

  # AH plot w/ regression lines##
  if(sub2$Flg1 %in% "Ads") {
    ggplot(sub, aes(y = MeanAH, x = Depth)) +
      ggtitle(paste0(siteFull, " \n", sub$Date_Time[1], "\n", "Adsorption")) +
      geom_line() +
      geom_point() +
      xlab("Depth (cm)") +
      ylab("Mean Absolute Humidity (µM)") +
      geom_errorbar(aes(ymin=MeanAH-StDevAH, ymax=MeanAH+StDevAH), width = .5) +
      coord_flip() +
      scale_x_reverse(breaks = c(0,10,20,30)) +
      scale_y_continuous(position = "right", limits = c(0,maxAH), breaks = c(0,300,600,900)) + ## Need to figure this out ##
      geom_vline(xintercept=0, linetype = "dashed") +
      stat_function(fun = AH.fun.up, color = "red", linetype = "dashed") +
      stat_function(fun = AH.fun.lo, color = "red", linetype = "dashed") +
      theme(plot.title = element_text(hjust = 0.5, size = 21), text = element_text(size = 18), axis.text.x = element_text(size = 16), aspect.ratio = 1.6, axis.line = element_line(color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
      ggsave(paste0(timSt, "_", "AH.png"), path = paste0(dir, "/AH_Reg_Plots"), dpi = 300, width = 4, height = 6.4, units = "in")

      ggplot(sub, aes(y = MeanAH, x = Depth)) +
        ggtitle(paste0(siteFull, " \n", sub$Date_Time[1], "\n", "Adsorption")) +
        geom_line() +
        geom_point() +
        xlab("Depth (cm)") +
        ylab("Mean Absolute Humidity (µM)") +
        geom_errorbar(aes(ymin=MeanAH-StDevAH, ymax=MeanAH+StDevAH), width = .5) +
        coord_flip() +
        scale_x_reverse(breaks = c(0,10,20,30)) +
        scale_y_continuous(position = "right", limits = c(0,maxAH), breaks = c(0,300,600,900)) +
        geom_vline(xintercept=0, linetype = "dashed") +
        theme(plot.title = element_text(hjust = 0.5, size = 21), text = element_text(size = 18), axis.text.x = element_text(size = 16), aspect.ratio = 1.6, axis.line = element_line(color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
        ggsave(paste0(timSt, "_", "AH.png"), path = paste0(dir, "/AH_Plots"), dpi = 300, width = 4, height = 6.4, units = "in")
  }
  else if(sub2$Flg1 %in% "Des") {
    ggplot(sub, aes(y = MeanAH, x = Depth)) +
      ggtitle(paste0(siteFull, " \n", sub$Date_Time[1], "\n", "Desorption")) +
      geom_line() +
      geom_point() +
      xlab("Depth (cm)") +
      ylab("Mean Absolute Humidity (µM)") +
      geom_errorbar(aes(ymin=MeanAH-StDevAH, ymax=MeanAH+StDevAH), width = .5) +
      coord_flip() +
      scale_x_reverse(breaks = c(0,10,20,30)) +
      scale_y_continuous(position = "right", limits = c(0,maxAH), breaks = c(0,300,600,900)) +
      geom_vline(xintercept=0, linetype = "dashed") +
      stat_function(fun = AH.fun.up, color = "red", linetype = "dashed") +
      stat_function(fun = AH.fun.lo, color = "red", linetype = "dashed") +
      theme(plot.title = element_text(hjust = 0.5, size = 21), text = element_text(size = 18), axis.text.x = element_text(size = 16), aspect.ratio = 1.6, axis.line = element_line(color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
      ggsave(paste0(timSt, "_", "AH.png"), path = paste0(dir, "/AH_Reg_Plots"), dpi = 300, width = 4, height = 6.4, units = "in")

    ggplot(sub, aes(y = MeanAH, x = Depth)) +
      ggtitle(paste0(siteFull, " \n", sub$Date_Time[1], "\n", "Desorption")) +
      geom_line() +
      geom_point() +
      xlab("Depth (cm)") +
      ylab("Mean Absolute Humidity (µM)") +
      geom_errorbar(aes(ymin=MeanAH-StDevAH, ymax=MeanAH+StDevAH), width = .5) +
      coord_flip() +
      scale_x_reverse(breaks = c(0,10,20,30)) +
      scale_y_continuous(position = "right", limits = c(0,maxAH), breaks = c(0,300,600,900)) +
      geom_vline(xintercept=0, linetype = "dashed") +
      theme(plot.title = element_text(hjust = 0.5, size = 21), text = element_text(size = 18), axis.text.x = element_text(size = 16), aspect.ratio = 1.6, axis.line = element_line(color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
      ggsave(paste0(timSt, "_", "AH.png"), path = paste0(dir, "/AH_Plots"), dpi = 300, width = 4, height = 6.4, units = "in")
  }
  else {
    ggplot(sub, aes(y = MeanAH, x = Depth)) +
      ggtitle(paste0(siteFull, " \n", sub$Date_Time[1], "\n", "No Reaction")) +
      geom_line() +
      geom_point() +
      xlab("Depth (cm)") +
      ylab("Mean Absolute Humidity (µM)") +
      geom_errorbar(aes(ymin=MeanAH-StDevAH, ymax=MeanAH+StDevAH), width = .5) +
      coord_flip() +
      scale_x_reverse(breaks = c(0,10,20,30)) +
      scale_y_continuous(position = "right", limits = c(0,maxAH), breaks = c(0,300,600,900)) +
      geom_vline(xintercept=0, linetype = "dashed") +
      theme(plot.title = element_text(hjust = 0.5, size = 21), text = element_text(size = 18), axis.text.x = element_text(size = 16), aspect.ratio = 1.6, axis.line = element_line(color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
      ggsave(paste0(timSt, "_", "AH.png"), path = paste0(dir, "/AH_Plots"), dpi = 300, width = 4, height = 6.4, units = "in")
      ggsave(paste0(timSt, "_", "AH.png"), path = paste0(dir, "/AH_Reg_Plots"), dpi = 300, width = 4, height = 6.4, units = "in")
  }

  ## T plot ##

  ggplot(sub, aes(y = MeanT, x = Depth)) +
    ggtitle(paste0(siteFull, " \n", sub$Date_Time[1])) +
    geom_line() +
    geom_point() +
    xlab("Depth (cm)") +
    ylab("Mean Temperature (ºC)") +
    geom_errorbar(aes(ymin=MeanT-StDevT, ymax=MeanT+StDevT), width = .5) +
    coord_flip() +
    scale_x_reverse(breaks = c(0,10,20,30)) +
    scale_y_continuous(position = "right", limits = c(0,maxT), breaks = c(0,25,50)) +
    geom_vline(xintercept=0, linetype = "dashed") +
    stat_function(fun = T.fun, color = "red", linetype = "dashed") +
    theme(plot.title = element_text(hjust = 0.5, size = 21), text = element_text(size = 18), axis.text.x = element_text(size = 16), aspect.ratio = 1.6, axis.line = element_line(color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    ggsave(paste0(timSt, "_", "T.png"), path = paste0(dir, "/T_Reg_Plots"), dpi = 300, width = 4, height = 6.4, units = "in")

  ggplot(sub, aes(y = MeanT, x = Depth)) +
    ggtitle(paste0(siteFull, " \n", sub$Date_Time[1])) +
    geom_line() +
    geom_point() +
    xlab("Depth (cm)") +
    ylab("Mean Temperature (ºC)") +
    geom_errorbar(aes(ymin=MeanT-StDevT, ymax=MeanT+StDevT), width = .5) +
    coord_flip() +
    scale_x_reverse(breaks = c(0,10,20,30)) +
    scale_y_continuous(position = "right", limits = c(0,maxT), breaks = c(0,25,50)) +
    geom_vline(xintercept=0, linetype = "dashed") +
    theme(plot.title = element_text(hjust = 0.5, size = 21), text = element_text(size = 18), axis.text.x = element_text(size = 16), aspect.ratio = 1.6, axis.line = element_line(color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    ggsave(paste0(timSt, "_", "T.png"), path = paste0(dir, "/T_Plots"), dpi = 300, width = 4, height = 6.4, units = "in")

}

## End section 4 ##
