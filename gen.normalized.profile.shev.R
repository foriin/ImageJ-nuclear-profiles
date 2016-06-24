
setwd('~/IMG/24.05 hp1 distribution in Kc and elav/csv/')
files <- dir()
nuclei <- list()
for (j in 1:length(files)) {
  file1 <- files[j]
  DATA <- read.csv(files[j], header = T, colClasses = c("character", "integer", "integer"))
  DATA[, 1] <- as.numeric(gsub(",", ".", DATA[, 1]))
  # HP1.mean <- mean(DATA$HP1.intensity)
  # Lam.mean <- mean(DATA$Lamin.intensity)
  DATA2 <- data.frame(DATA$Coordinates, DATA$Lamin.intensity, DATA$HP1.intensity)
  # diff(diff(DATA2$DATA.Lamin.intensity.Lam.mean))
  
  names(DATA2) <- c("coord", "lam", "hp1")

  
  
  norm.table <- function(tab, by=0.01){
    out.data <- data.frame()
    for (i in seq(-0.2, 1.2, by=by)){
      lam.sig <- mean(tab[, 2][tab[, 1] >= i & tab[, 1] < i + by])
      hp1.sig <- mean(tab[, 3][tab[, 1] >= i & tab[, 1] < i + by])
      out.data <- rbind(out.data, c(i, lam.sig, hp1.sig))
      
    }
    names(out.data) <- names(tab)
    out.data[is.na(out.data)] <- 0
    return(out.data)
  }
  # DATA3 <- norm.table(DATA2)

  
  # which(diff(sign(diff(DATA2$lam)))==-2)+1
  a <- max(which(DATA2$lam[1:(length(DATA2$lam)/2)] == max(DATA2$lam[1:(length(DATA2$lam)/2)])))
  b <- min(which(DATA2$lam[(length(DATA2$lam)/2):length(DATA2$lam)] == max(DATA2$lam[(length(DATA2$lam)/2):length(DATA2$lam)])) + length(DATA2$lam)/2 - 1)
  border <- b - a
  diameter <- DATA2$coord[b] - DATA2$coord[a]
  DATA2$coord <- DATA2$coord/diameter
  DATA2$coord <- DATA2$coord - DATA2$coord[a]
  bl <- max(round(a - (border/5)), 1)
  br <- min(round(b + border/5), length(DATA2$lam))
  lam.borders <- DATA2$lam[bl:br]
  HP1.in.nuc <- DATA2$hp1[bl:br]
  new.coord <- DATA2$coord[bl:br]
  DATA4 <- data.frame(new.coord, lam.borders, HP1.in.nuc)
  names(DATA4) <- names(DATA2)
  
  HP1.mean <- mean(DATA2$hp1[a:b])
  Lam.mean <- mean(DATA2$lam[a:b])
  DATA4$lam <- DATA4$lam/Lam.mean/2
  DATA4$hp1 <- DATA4$hp1/HP1.mean
  
  DATA4 <- norm.table(DATA4, by=0.05)
  # DATA4[DATA4$coord == 0,]$lam <- max(DATA4$lam[1:(length(DATA4$lam)/2)])
  nuclei[[j]] <- DATA4
  
  
}

FINAL.DATA <- Reduce("+", nuclei)/25

FINAL.DATA[FINAL.DATA$coord == 0,]$lam <- max(FINAL.DATA$lam[1:(length(FINAL.DATA$lam)/2)])
plot(FINAL.DATA$coord, FINAL.DATA$lam, t="l", xlab="x coordinate", ylab="fluorescence")
lines(FINAL.DATA$coord, FINAL.DATA$hp1, t="h", col="green")

 write.table(FINAL.DATA, file="../Kc_total_HP1_profile_25_nuclei_30.csv", row.names = F)


