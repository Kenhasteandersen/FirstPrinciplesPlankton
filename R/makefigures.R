require(latex2exp)
require(lattice)
source("basetools.R")
source("model.R")
source("modelChemostat.R")

#
# Conditions for the two examples
#
dOligotrophic = 0.001
LOligotrophic = 60
dEutrophic = 0.1
LEutrophic = 60

# ===================================================
# Plots for documentation
# ===================================================

plotAll = function() {
  pdfplot("../aL.pdf", plot_aL, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../aF.pdf", plot_aF, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../aN.pdf", plot_aN, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../Mumax.pdf", plotMumax, width = 1.5*singlewidth, height=1.5*height)
  pdfplot("../Rstar.pdf", plotRstar, width=doublewidth, height = 1.25*height)
  #pdfplot("../Mumax_corrected.pdf", plotMuAlphaCorrelation, width = 1.5*singlewidth, height=1.5*height)
  
  plotSimulationExamples()
  
  pdfplot("../Strategies.pdf", plotStrategies, 
          width=singlewidth, height=3*height)
  
  pdfplot("../SheldonComparison.pdf", plotSheldonComparison, width=doublewidth, height=doublewidth)
  
  #fontsize = trellis.par.get("fontsize")
  #fontsize$text = 10
  #trellis.par.set("fontsize",fontsize)
  #plt=plotVaryLightAndDiffusion(n=15)
  #pdf(file="../VaryLightAndDiffusion.pdf", width=doublewidth+0.5, height=height+.3)
  #plt
  #dev.off()
  
  pdfplot("../Functions.pdf", plotFunctions,n=40, width=doublewidth, height=2.5*height)
  
  pdfplot("../GridPreference.pdf", plotGridPreference, width=1.5*singlewidth, height=height)
  
  #pdfplot("../Gridtest.pdf", plotGridtest, width=doublewidth, height=height)
  
  pdfplot("../DOC.pdf", plotDOC, width=doublewidth, height=1.75*height)
  
  pdfplot("../BacterialGenerationTime.pdf", 
          plotBacteriaGenerationTime_vs_area, height=height)
  
  pdfplot("../HTL.pdf", plotHTL, height=2*height)
  
  pdfplot("../Temperature.pdf", plotTemperature, width=doublewidth, height=2*height)
  
  #system("cp ../*pdf ../../dropbox")
}

convertVolume2Mass = function(vol, taxon="other") {
  C = 1e-6 * exp(-0.665) * vol^0.939 # Menden-deyer and Lessard (2000), table 4, for protists. mugC
  ixSmall = (vol<3000) & (!is.na(vol))
  C[ixSmall] = 1e-6 * exp(-0.583) * vol[ixSmall]^0.860
  ixDiatom = (taxon=="diatom") & (vol>3000) & (!is.na(vol))
  C[ixDiatom] = 1e-6 * exp(-0.933) * vol[ixDiatom]^0.881 # As above for diatom>3000 mum^3
  ixDiatom = (taxon=="diatom") & (vol<=3000)  & (!is.na(vol))
  C[ixDiatom] = 1e-6 * exp(-0.541) * vol[ixDiatom]^0.811 # As above for diatom>3000 mum^3
  return(C)
}

# plotAL = function() {
#   data = data.frame(C=NA,taxon=NA,AL=NA)
#   #
#   # Taguchi
#   #
#   #Ata=read.csv("../data/Taguchi.dat",  header=FALSE, col.names=c("V ((mum)^3)", "err", "C (pgC", "CperChl", "alpha (mgC/(mg chlA) /h /W m^2)"), sep=" ")
#   #C = (Ata$C..pgC)*1e-6 # mugC
#   #chl = 1e-3 * C/Ata$CperChl # mg chl
#   #A = 24 * 1000* Ata$alpha * chl # mugC/d/(Wm2)
#   #data = data.frame(C=C, taxon="diatom", A=A, source="Taguchi")
#   #ALtaguchi = exp(mean(log(A/(C^(2/3)))))
#   #cat("AL = ", AL, "mugC/d/(wm2)\n")
#   #
#   # Edwards:
#   #
#   Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
#                  sep=";", skip=3, header=TRUE, na.strings = "na")
#   C = convertVolume2Mass(Aed$volume, Aed$taxon)
#   
#   A = Aed$alpha * C * (Aed$daylength/24) # convert to units of mu gC/(mu mol photons/m2/s),
#   # corrected for daylength.
#   data = data.frame(C=C, taxon=Aed$taxon, A=A, source="Edwards")#rbind(data, 
#   
#   #
#   # Sort our nans:
#   #
#   data = data[(!is.na(data$A) & !is.na(data$C)), ]
#   # 
#   # Fits to complex shading formula:
#   #
#   ixDiatom = data$taxon=="diatom"
#   
#   #form = formula(log(A) ~ log( a*C^(2/3) * Cmax*C / (a*C^(2/3) + Cmax*C)))
#   form = formula(log(A) ~ log( a*C^(2/3) * (1 - exp(-Cmax*C^(1/3)) ) ))
#   
#   fit = nls(form,
#             data = data,
#             start = list(Cmax = .001, a=0.001),
#             lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
#   fit_diatoms = nls(form,
#                     data = data[ixDiatom,],
#                     start = list(Cmax = .001, a=0.001),
#                     lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
#   fit_no_diatoms = nls(form,
#                        data = data[!ixDiatom,],
#                        start = list(Cmax = .001, a=0.001),
#                        lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
#   
#   print(summary(fit_no_diatoms))
#   #
#   # Fits to 2/3 scaling
#   #
#   AL = exp(mean(log(data$A/data$C^(2/3))))
#   AL_diatoms = exp(mean(log(data$A/data$C^(2/3))[ixDiatom]))
#   AL_no_diatoms = exp(mean(log(data$A/data$C^(2/3))[!ixDiatom]))
#   cat(AL_no_diatoms)
#   #
#   # Plot
#   #
#   defaultplot()
#   loglogpanel(xlim = data$C, 
#               ylim = c(0.9*min(data$A[!is.na(data$A)]), 1.1*max(data$A[!is.na(data$A)])),
#               xlab="Cell weight ($\\mu$gC)",
#               ylab="Affinity for light, $\\textit{A_L}$ ($\\mu$gC/($\\mu$ mol photons m$^{-2}s^{-1})")
#   points(data$C[!ixDiatom], data$A[!ixDiatom],pch=16, col="blue")
#   points(data$C[ixDiatom], data$A[ixDiatom],pch=16, col="red")
#   points(data$C[data$source=="Taguchi"], data$A[data$source=="Taguchi"], pch=17, col="red")
#   
#   C = 10^seq(log10(min(data$C)), log10(max(data$C)), length.out = 100)
#   lines(C, exp(predict(fit, list(C=C))), lwd=2)
#   lines(C, exp(predict(fit_diatoms, list(C=C))), lwd=2, col="red")
#   lines(C, exp(predict(fit_no_diatoms, list(C=C))), lwd=3, col="blue")
#   
#   lines(C, AL * C^(2/3), lty=dotted)
#   lines(C, AL_diatoms * C^(2/3), col="red", lty=dotted)
#   lines(C, AL_no_diatoms * C^(2/3), col="blue", lty=dotted)
#   
#   # Camila's parameters:
#   cL=0.08733
#   AL=0.004864  
#   lines(C, cL*C * AL*C^(2/3) / ( cL*C + AL*C^(2/3) ) , col="orange")
#   
#   legend(x="bottomright", bty="n",
#          legend=c("Diatoms", "Other phototrophs", "Camila"),
#          pch=c(16,16,NA),
#          lwd=c(0,0,2),
#          col=c("red","Blue","Orange"))
# }
# 
# plotALsimple = function() {
#   data = data.frame(C=NA,taxon=NA,AL=NA)
#   #
#   # Taguchi
#   #
#   #Ata=read.csv("../data/Taguchi.dat",  header=FALSE, col.names=c("V ((mum)^3)", "err", "C (pgC", "CperChl", "alpha (mgC/(mg chlA) /h /W m^2)"), sep=" ")
#   #C = (Ata$C..pgC)*1e-6 # mugC
#   #chl = 1e-3 * C/Ata$CperChl # mg chl
#   #A = 24 * 1000* Ata$alpha * chl # mugC/d/(Wm2)
#   #data = data.frame(C=C, taxon="diatom", A=A, source="Taguchi")
#   #ALtaguchi = exp(mean(log(A/(C^(2/3)))))
#   #cat("AL = ", AL, "mugC/d/(wm2)\n")
#   #
#   # Edwards:
#   #
#   Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
#                  sep=";", skip=3, header=TRUE, na.strings = "na")
#   C = convertVolume2Mass(Aed$volume, Aed$taxon)
#   r = (Aed$volume*3/(4*pi))^(1/3)
#   
#   A = Aed$alpha * C * (Aed$daylength/24) # convert to units of mu gC/day/(mu mol photons/m2/s),
#   # corrected for daylength.
#   data = data.frame(r=r, C=C, taxon=Aed$taxon, A=A, source="Edwards")#rbind(data, 
#   
#   #
#   # Sort our nans:
#   #
#   data = data[(!is.na(data$A) & !is.na(data$C)), ]
#   # 
#   # Fits to complex shading formula:
#   #
#   ixDiatom = data$taxon=="diatom"
#   
#   #form = formula(log(A) ~ log( a*C^(2/3) * Cmax*C / (a*C^(2/3) + Cmax*C)))
#   form = formula(log(A) ~ log( AL*r^2 * (1 - exp(-cL*C^(1/3)) ) ))
#   
#   fit = nls(form,
#             data = data,
#             start = list(cL = 1, AL=1e-7),
#             lower=list(cL=1e-20, AL=1e-20), algorithm="port")
#   
#   print(summary(fit))
#   cat(c("r* = ",1/(coef(fit)[[1]]*(0.19e-6*4/3*pi)^(1/3)) ), "mu m\n")
#   #
#   # Fits to 2/3 scaling
#   #
#   AL = exp(mean(log(data$A/data$C^(2/3))))
#   AL_diatoms = exp(mean(log(data$A/data$C^(2/3))[ixDiatom]))
#   AL_no_diatoms = exp(mean(log(data$A/data$C^(2/3))[!ixDiatom]))
#   cat(AL_no_diatoms)
#   #
#   # Plot
#   #
#   defaultplot()
#   loglogpanel(xlim = c(0.1, 300), 
#               ylim = c(0.5*min(data$A[!is.na(data$A)]), 2*max(data$A[!is.na(data$A)])),
#               xlab="Cell radius ($\\mu$m)",
#               ylab="Affinity for light, $\\textit{A_L}$ ($\\mu$gC/day/($\\mu$ mol photons m$^{-2}s^{-1})")
#   points(data$r[!ixDiatom], data$A[!ixDiatom],pch=15, col="darkgreen")
#   points(data$r[ixDiatom], data$A[ixDiatom],pch=16, col="darkgreen")
#   points(data$r[data$source=="Taguchi"], data$A[data$source=="Taguchi"], pch=17, col="darkgreen")
#   
#   r = 10^seq(-1, 4, length.out = 100)
#   C = 4/3*pi*r^3*0.2e-6
#   lines(r, exp(predict(fit, list(r=r,C=C))), lwd=2)
#   
#   lines(r, coef(fit)[[2]]*r^2, lty=dotted)
#   lines(r, coef(fit)[[2]]*r^2*coef(fit)[[1]]*C^(1/3), lty=dotted)
#   
#   legend(x="bottomright", bty="n",
#          legend=c("Diatoms", "Other phototrophs"),
#          pch=c(15,16),
#          lwd=0,
#          col="darkgreen")
# }

plot_aL = function() {
  data = data.frame(r=NA,C=NA,taxon=NA,aL=NA)
  #
  # Taguchi
  #
  #Ata=read.csv("../data/Taguchi.dat",  header=FALSE, col.names=c("V ((mum)^3)", "err", "C (pgC", "CperChl", "alpha (mgC/(mg chlA) /h /W m^2)"), sep=" ")
  #C = (Ata$C..pgC)*1e-6 # mugC
  #chl = 1e-3 * C/Ata$CperChl # mg chl
  #A = 24 * 1000* Ata$alpha * chl # mugC/d/(Wm2)
  #data = data.frame(C=C, taxon="diatom", A=A, source="Taguchi")
  #ALtaguchi = exp(mean(log(A/(C^(2/3)))))
  #cat("AL = ", AL, "mugC/d/(wm2)\n")
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  r = (Aed$volume*3/(4*pi))^(1/3)
  
  aL = Aed$alpha * (Aed$daylength/24) # convert to units of mu gC/day/(mu mol photons/m2/s),
  # corrected for daylength.
  data = data.frame(r=r, C=C, taxon=Aed$taxon, aL=aL, source="Edwards")#rbind(data, 
  
  #
  # Sort out nans:
  #
  data = data[(!is.na(data$aL) & !is.na(data$C)), ]
  # 
  # Fits to complex shading formula:
  #
  ixDiatom = data$taxon=="diatom"
  
  #form = formula(log(A) ~ log( a*C^(2/3) * Cmax*C / (a*C^(2/3) + Cmax*C)))
  delta = parameters()$delta
  form = formula(log(aL) ~ log( alphaL/r * (1 - exp(-r/rStar) ) * (1-3*delta/r)  ))
  
  fit = nls(form,
            data = data,
            start = list(rStar = 1.5, alphaL=1),
            lower=list(rStar=0, alphaL=1e-20), 
            upper=list(rStar=10, alphaL=1e10), algorithm="port")
  
  print(summary(fit))
  cat(c("alpha_L = ", coef(fit)[[2]], "1/day/mugC*mum/(light)"))
  cat(c("r* = ", coef(fit)[1], "mu m\n"))
  #
  # Fits to 2/3 scaling
  #
  #alphaL = exp(mean(log(data$A/data$C^(2/3))))
  #AL_diatoms = exp(mean(log(data$A/data$C^(2/3))[ixDiatom]))
  #AL_no_diatoms = exp(mean(log(data$A/data$C^(2/3))[!ixDiatom]))
  #cat(AL_no_diatoms)
  #
  # Plot
  #
  defaultplot()
  loglogpanel(xlim = c(0.1, 300), 
              ylim = c(0.5*min(data$aL[!is.na(data$aL)]), 2*max(data$aL[!is.na(data$aL)])),
              xlab="Cell radius ($\\mu$m)",
              ylab="Light affinity $\\textit{a_L}$ (day$\\cdot \\mu$ mol photons m$^{-2}s^{-1}\\cdot \\mu$g$_C)^{-1}$")
  points(data$r[!ixDiatom], data$aL[!ixDiatom],pch=15, col="darkgreen")
  points(data$r[ixDiatom], data$aL[ixDiatom],pch=16, col="darkgreen")
  points(data$r[data$source=="Taguchi"], data$aL[data$source=="Taguchi"], pch=17, col="darkgreen")
  
  r = 10^seq(-1, 4, length.out = 100)
  lines(r, exp(predict(fit, list(r=r,C=C))), lwd=2)
  
  lines(r, coef(fit)[[2]]/r, lty=dotted)
  lines(r, coef(fit)[[2]]/coef(fit)[[1]]*r/r, lty=dotted)
  lines(r, coef(fit)[[2]]/coef(fit)[[1]] * (1-3*delta/r), lty=dotted)
  
  legend(x="topright", bty="n",
         legend=c("Diatoms", "Other phototrophs"),
         pch=c(15,16),
         lwd=0,
         col="darkgreen")
}

# plotALvolume = function() {
#   data = data.frame(r=NA, V=NA,taxon=NA,AL=NA)
#   #
#   # Taguchi
#   #
#   #Ata=read.csv("../data/Taguchi.dat",  header=FALSE, col.names=c("V ((mum)^3)", "err", "C (pgC", "CperChl", "alpha (mgC/(mg chlA) /h /W m^2)"), sep=" ")
#   #C = (Ata$C..pgC)*1e-6 # mugC
#   #chl = 1e-3 * C/Ata$CperChl # mg chl
#   #A = 24 * 1000* Ata$alpha * chl # mugC/d/(Wm2)
#   #data = data.frame(C=C, taxon="diatom", A=A, source="Taguchi")
#   #ALtaguchi = exp(mean(log(A/(C^(2/3)))))
#   #cat("AL = ", AL, "mugC/d/(wm2)\n")
#   #
#   # Edwards:
#   #
#   Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
#                  sep=";", skip=3, header=TRUE, na.strings = "na")
#   C = convertVolume2Mass(Aed$volume, Aed$taxon)
#   
#   A = (Aed$alpha) * (Aed$daylength/24)
#   data = data.frame(V=Aed$volume, r = (3/(4*pi)*Aed$volume)^(1/3), 
#                     Aed$taxon, A=A, source="Edwards", taxon=Aed$taxon,
#                     mumax=Aed$mu_max)#rbind(data, 
#   
#   #
#   # Sort our nans:
#   #
#   data = data[(!is.na(data$A) & !is.na(data$r)), ]
#   # 
#   # Fits to complex shading formula:
#   #
#   ixDiatom = data$taxon=="diatom"
#   
#   form = formula(log(A) ~ log( cL*r^-1*(1-exp(-theta*r))))
#   
#   fit = nls(form,
#             data = data,
#             start = list(cL =1e-1, theta=1),
#             lower=list(cL=1e-10, theta=0.01), algorithm="port", trace=TRUE)
#   fit_diatoms = nls(form,
#                     data = data[ixDiatom,],
#                     start = list(cL =1e-1, theta=1),
#                     lower=list(cL=1e-10, theta=0.01), algorithm="port", trace=TRUE)
#   fit_no_diatoms = nls(form,
#                        data = data[!ixDiatom,],
#                        start = list(cL =1e-1, theta=1),
#                        lower=list(cL=1e-10, theta=0.01), algorithm="port", trace=TRUE)
#   
#   print(summary(fit_no_diatoms))
#   #
#   # Fits to -1 scaling
#   #
#   AL = exp(mean(log(data$A/data$r^(-1))))
#   AL_diatoms = exp(mean(log(data$A/data$r^(-1))[ixDiatom]))
#   AL_no_diatoms = exp(mean(log(data$A/data$r^(-1))[!ixDiatom]))
#   cat(AL_no_diatoms)
#   #
#   # Plot
#   #
#   defaultplotvertical(2)
#   loglogpanel(xlim = data$r, ylim = c(0.9*min(data$A[!is.na(data$A)]), 1.1*max(data$A[!is.na(data$A)])),
#               xaxis=FALSE,
#               ylab="Affinity for light, $\\textit{A_L}$ ($\\mu$gC/day/($\\mu$ mol photons m$^{-2}s^{-1})")
#   points(data$r[!ixDiatom], data$A[!ixDiatom],pch=16, col="blue")
#   points(data$r[ixDiatom], data$A[ixDiatom],pch=16, col="red")
#   
#   r = 10^seq(log10(min(data$r)), log10(max(data$r)), length.out = 100)
#   lines(r, exp(predict(fit, list(r=r))), lwd=2)
#   lines(r, exp(predict(fit_diatoms, list(r=r))), lwd=2, col="red")
#   lines(r, exp(predict(fit_no_diatoms, list(r=r))), lwd=3, col="blue")
#   
#   lines(r, AL * r^(-1), lty=dotted)
#   lines(r, AL_diatoms * r^(-1), col="red", lty=dotted)
#   lines(r, AL_no_diatoms * r^(-1), col="blue", lty=dotted)
#   
#   legend(x="topright", bty="n",
#          legend=c("Diatoms", "Other phototrophs"),
#          pch=16, col=c("red","Blue"))
#   #
#   # Imax
#   #
#   fit_mu_diatoms = lm(log(mumax) ~ log(r), data=data[ixDiatom,])
#   fit_mu_no_diatoms = lm(log(mumax) ~ log(r), data=data[!ixDiatom,])
#   fit_mu = lm(log(mumax) ~ log(r), data=data)
#   
#   loglogpanel(xlim = data$r, ylim = data$mumax,
#               xlab="radius ($\\mu$m)",
#               ylab="max growth rate ($d^{-1}$)")
#   points(data$r[!ixDiatom], data$mumax[!ixDiatom], pch=16, col="blue")
#   points(data$r[ixDiatom], data$mumax[ixDiatom], pch=16, col="red")
#   
#   lines(r, exp(predict(fit_mu_diatoms, list(r=r))), lwd=2, col="red")
#   lines(r, exp(predict(fit_mu, list(r=r))), lwd=2)
#   lines(r, exp(predict(fit_mu_no_diatoms, list(r=r))), lwd=2, col="blue")
#   #
#   # Imax and alpha
#   #
#   defaultplot()
#   x = data$A/exp(predict(fit, list(r=data$r)))
#   y = data$mumax/exp(predict(fit_mu, list(r=data$r)))
#   loglogpanel(xlim=x, ylim=x,
#               xlab="alpha corrected", ylab="mumax corrected")
#   points(x[ixDiatom], y[ixDiatom], pch=16, col="red")
#   points(x[!ixDiatom], y[!ixDiatom], pch=16, col="blue")
#   
#   lines(c(0.01,100), c(0.01,100),lty=dotted)       
# }

plotMumax = function() {
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  Aed$C =   C
  # Sort our nans:
  Aed = Aed[!is.na(Aed$C), ]
  # Correct to 10 degrees with a Q10 of 1.5
  Aed$mu_max = Aed$mu_max * 1.5^((10-Aed$temperature)/10)
  
  A = data.frame(species=Aed$species, taxon=as.factor(Aed$taxon), 
                 C=Aed$C, mu_max=Aed$mu_max)
  #
  # Kiørboe and hirst (2014):
  #
  Akh = read.csv('../Data/Kiorboe and Hirst (2014) maximum growth rates.csv', 
                 sep=",", header=TRUE, skip=1, as.is=TRUE)
  Akh$Group[Akh$Group=="Dinoflagellates"] = "dinoflagellate"
  Akh$Group[Akh$Group=="Nanoflagellates (except dinoflagellates)"] = "nanoflagellate"
  Akh$Group[Akh$Group=="Ciliates"] = "ciliate"
  # Correct from 15 to 10 degrees with a Q10 of 2.8
  Akh$Specific.growth.1 = Akh$Specific.growth.1 * 2.8^((10-15)/10)
  
  A = rbind(A, 
            data.frame(species=Akh$Species, taxon=as.factor(Akh$Group),
                       C=Akh$Body.Mass*1000, mu_max=Akh$Specific.growth.1*24))
  #
  # Kirchman (2016) - but not max growth rates
  #
  Ak = read.csv('../Data/ma08_kirchman_supappendix3.csv',
                sep=",", header=TRUE, skip=15)
  Ak = Ak[1:149,]
  # Sort out incomplete entries:
  Ak = Ak[Ak[,7]!="Used biovol*",]
  Ak = Ak[!is.na(Ak[,5]),]
  
  A = rbind(A,
            data.frame(species="NA",
                       taxon="bacteria", 
                       C=as.numeric(as.character(Ak$fgC.cell))*1e-9, 
                       mu_max=as.numeric(Ak$Rate...d.)))
  #
  # Rose and Caron (2007)
  #
  Arc = read.csv('../Data/Rose and Caron bacterivores.csv',
                 sep=",", header=TRUE)
  Arc$volume = as.numeric( gsub(",","", as.character(Arc$volume)) )
  # Correct to 10 degrees with a Q10 of 2.8
  Arc$Growth.Rate = Arc$Growth.Rate * 2.8^((10-Arc$Temperature)/10)
  
  A = rbind(A, data.frame(
    species="NA",
    taxon="bacterivore",
    C = convertVolume2Mass(Arc$volume),
    mu_max=as.numeric(Arc$Growth.Rate)))
  
  Arc = read.csv('../Data/Rose and Caron herbivores.csv',
                 sep=",", header=TRUE)
  Arc$volume = as.numeric( gsub(",","", as.character(Arc$volume)) )
  # Correct to 10 degrees with a Q10 of 2.8
  Arc$Growth.Rate = Arc$Growth.Rate * 2.8^((10-Arc$Temperature)/10)
  
  A = rbind(A,data.frame(
    species="NA",
    taxon="herbivore",
    C = convertVolume2Mass(Arc$volume),
    mu_max = as.numeric(Arc$Growth.Rate)))
  #
  # Define groups:
  #
  ixDiatom = A$taxon=="diatom"
  ixMixotroph = A$taxon=="nanoflagellate" | A$taxon=="dinoflagellate" | A$taxon=="ciliate" | A$taxon=="bacteriovore"
  ixHeterotroph = A$taxon=="herbivore" 
  ixBacteria = A$taxon=="bacteria"
  ixPhototroph = !ixDiatom & !ixMixotroph  &!ixBacteria &
    A$taxon!="bacterivore" & A$taxon!="herbivore"
  #
  # Fits
  #
  #fit_diatoms = lm(log(mu_max) ~ log(C), data=A[ixDiatom,])
  #fit_mixotrophs = lm(log(mu_max) ~ log(C), data=A[ixMixotroph,])
  #fit_heterotrophs = lm(log(mu_max) ~ log(C), data=A[ixHeterotroph,])
  #fit_phototrophs = lm(log(mu_max) ~ log(C), data=A[ixPhototroph,])
  
  #fit = nls(mu_max ~ (1-parameters()$c * C^(-1/3))*mu/(mu*c*C^0.3333+1), 
  #          data=A[!ixBacteria,],
  #          start=list(mu=1, c=1))
  
  defaultplot()
  semilogxpanel(xlim = c(1e-9,1), ylim = A$mu_max,
                xlab="Cell weight ($\\mu$gC)",
                ylab="max growth rate ($d^{-1}$)")
  #points(A$C[ixHeterotroph], A$mu_max[ixHeterotroph], pch=16, col="red")
  points(A$C[ixMixotroph], A$mu_max[ixMixotroph], pch=16, col="blue")
  points(A$C[A$taxon=="bacterivore"], A$mu_max[A$taxon=="bacterivore"], pch=1, col="blue", cex=.1)
  points(A$C[A$taxon=="herbivore"], A$mu_max[A$taxon=="herbivore"], pch=1, col="red", cex=0.1)
  points(A$C[ixPhototroph], A$mu_max[ixPhototroph], pch=16, col="green")
  points(A$C[ixDiatom], A$mu_max[ixDiatom], pch=16, col="darkgreen")
  points(A$C[ixBacteria], A$mu_max[ixBacteria], pch=1, col="brown", cex=0.1)
  
  
  C = 10^seq(-9, log10(max(A$C)), length.out = 100)
  #lines(C, exp(predict(fit_diatoms, list(C=C))), col="darkgreen")
  #lines(C, exp(predict(fit_mixotrophs, list(C=C))), col="blue")
  #lines(C, exp(predict(fit_heterotrophs, list(C=C))), col="red")
  #lines(C, exp(predict(fit_phototrophs, list(C=C))), col="green")
  
  m = 10^seq(-9,1,length.out = 100)
  lines(m, parameters()$alphaJ*(1-parameters()$c * m^(-1/3)), lwd=3)
  #lines(m, parameters()$alphaJ*(1-parameters()$c * m^(-1/3)) - parameters()$cR, lwd=2)
  #lines(C, predict(fit, list(C=C)))
  
  # Add the used curve:
  p = parameters()
  #lines(p$m, p$Jmax/p$m, col="black", lwd=2)
  
  # Maximum uptake of phagotrophy:
  lines(p$m, p$jFmaxm, col="red", lwd=2)
  
  # Camila
  #lines(C, 0.12*C^-0.25, col="orange", lwd=1)
  #lines(C, (0.12-0.03)*C^-0.25, col="orange", lwd=2)
  
  legend(x="topleft", bty="n",
         legend=c("Diatoms","Other phototrophs",
                  "Mixotrophs","Heterotrophs",
                  "Bacteria","Model","Max. phagotrophy"),
         lty=c(NA,NA,NA,NA,NA,solid,solid),
         pch=c(16,16,16,16,16,NA,NA),
         lwd=c(0,0,0,0,0,2,2),
         col=c("darkgreen","green","blue","red","brown","black","red"))
}

plotMuAlphaCorrelation = function() {
  #
  # Edwards:
  #
  Aed = read.csv("../data/Data from Edwards et al (2015).csv", 
                 sep=";", skip=3, header=TRUE, na.strings = "na")
  C = convertVolume2Mass(Aed$volume, Aed$taxon)
  Aed$C = C
  
  A = (Aed$alpha * 4.15) * C * (Aed$daylength/24)
  data = data.frame(C=C, taxon=Aed$taxon, A=A, mumax=Aed$mu_max)
  #
  # Sort our nans:
  #
  data = data[(!is.na(data$A) & !is.na(data$C)), ]
  ixDiatom = data$taxon=="diatom"
  data = data[ixDiatom,]
  #
  # Fits:
  #
  fit_mu = lm(log(mumax) ~ log(C), data=data)
  form = formula(log(A) ~ log( Cmax*C * a*C^(2/3) / (Cmax*C + a*C^(2/3))))
  fit_alpha = nls(form,
                  data = data,
                  start = list(Cmax = .001, a=0.001),
                  lower=list(Cmax=1e-20, a=1e-20), algorithm="port")
  #
  # Plot corrected values (using residuals)
  #
  defaultplot(mfcol=c(1,2))
  loglogpanel(xlim=c(0.1,10),ylim=c(0.1,10),
              xlab="Residual light-affinity",
              ylab="Residual $\\mu_{max}$")
  x = data$A/exp(predict(fit_alpha, list(C=data$C)))
  y = data$mumax/exp(predict(fit_mu, list(C=data$C)))
  points(x[ixDiatom], y[ixDiatom], pch=16, col="red")
  points(x[!ixDiatom], y[!ixDiatom], pch=16, col="blue")
  
  #lines(c(0.1,10), c(0.1,10),lty=dotted)     
  
  fit = lm(log(x) ~ log(y))
  xx = 10^seq(-1,2,length.out = 100)
  #lines(xx, exp(fit$coefficients[1])*xx^fit$coefficients[2])
  legend(x="topleft", pch=16, col=c("red","blue"),
         legend=c("Diatoms","Others"))
  #
  # Plot mass-specific values
  #
  loglogpanel(xlim=c(0.001,2),ylim=c(0.1,5),
              xlab="Specific light-affinity",
              ylab="Specific $\\mu_{max}$")
  x = data$A/data$C
  y = data$mumax
  points(x[ixDiatom], y[ixDiatom], pch=16, col="red")
  points(x[!ixDiatom], y[!ixDiatom], pch=16, col="blue")
  
  #lines(c(0.1,10), c(0.1,10),lty=dotted)     
  legend(x="topleft", pch=16, col=c("red","blue"),
         legend=c("Diatoms","Others"))
  
  fit = lm(log(x) ~ log(y))
  xx = 10^seq(-3,1,length.out = 100)
  #lines(xx, exp(fit$coefficients[1])*xx^fit$coefficients[2])
}


# plotAF = function() {
#   dat <- read.csv("../data/TK Appendix feeding rates - revised.csv",header=TRUE,sep=";")
#   data = data.frame(w=1e3*dat$Body.mass, beta=24*0.001*dat$Fmax.1, Group=dat$Group)  
#   
#   ixProtist = (data$Group=="Nanoflagellates") | 
#     (data$Group=="Dinoflagellates") | 
#     (data$Group=="Dinoflagellate") | 
#     (data$Group=="Ciliates") | 
#     (data$Group=="Ciliate")
#   
#   x = data$beta/data$w
#   x = x[ixProtist]
#   x = x[!is.na(x)]
#   SpecificBeta = exp(mean(log(x)))
#   cat(SpecificBeta)
#   
#   defaultplot()
#   loglogpanel(xlim=c(1e-6, 1), ylim=c(1e-8,1e-2),
#               xlab = "Mass ($\\mu$gC)", ylab="Clearance rate $\\textit{A_F}$ (l/d)")
#   points(data$w[ixProtist], data$beta[ixProtist])
#   lines(data$w[ixProtist], SpecificBeta*data$w[ixProtist], lwd=2)
#   
#   # Camila
#   C = range(data$w[ixProtist])
#   lines(C, 0.018*C, col="orange")
# }

plot_aF = function() {
  dat <- read.csv("../data/TK Appendix feeding rates - revised.csv",header=TRUE,sep=";")
  data = data.frame(w=1e3*dat$Body.mass, beta=24*0.001*dat$Fmax.1, Group=dat$Group)  
  
  ixProtist = (data$Group=="Nanoflagellates") | 
    (data$Group=="Dinoflagellates") | 
    (data$Group=="Dinoflagellate") | 
    (data$Group=="Ciliates") | 
    (data$Group=="Ciliate")
  
  x = data$beta/data$w
  x = x[ixProtist]
  x = x[!is.na(x)]
  SpecificBeta = exp(mean(log(x)))
  cat(c("aF = ", SpecificBeta, "L/day/mugC\n"))
  
  defaultplot()
  loglogpanel(xlim=c(1e-6, 1), ylim=c(1e-4,1),
              xlab = "Mass ($\\mu$gC)", ylab="Specfic clearance rate $\\textit{a_F}$ (L/d/{\\mu}gC)")
  points(data$w[ixProtist], data$beta[ixProtist]/data$w[ixProtist], pch=16)
  w = 10^seq(-7,1)
  lines(w, SpecificBeta*w/w, lwd=2)
}
#
# Plot specific affinity:
#
plot_aN = function() {
  dat = read.csv("../data/Nutrient data from Edwards et al (2015b).csv",
                 header=TRUE,sep=",")
  #
  # Convert carbon from mol to g:
  #
  dat$c_per_cell = dat$c_per_cell*12
  #
  # Convert to weights when only volume is given:
  # 
  ix = (!is.na(dat$volume)) & (is.na(dat$c_per_cell))
  dat$c_per_cell[ix] = convertVolume2Mass(dat$volume[ix], dat$taxon[ix])
  C = dat$c_per_cell
  #
  # Convert to ESD/2 (radius) assuming spherical cells
  #
  r = ( 3/(4*pi)*dat$volume )^(1/3)
  
  #
  # Calc affinities
  #
  a_nit = dat$vmax_nit / dat$k_nit / dat$c_per_cell
  a_amm = dat$vmax_amm / dat$k_amm / dat$c_per_cell 
  a_phos = dat$vmax_p / dat$k_p / dat$c_per_cell
  #
  # Plot
  #
  defaultplot()
  loglogpanel(xlim=c(0.1 ,200), ylim=c(1e-5, 20), #c(1e-9 ,1)
              xlab="Radius ($\\mathrm{\\mu}$m)",
              ylab="Nutrient affinity, $\\textit{a}_N$ (L/day/$\\mathrm{\\mu}$gC)")
  
  col = 1
  taxons = unique(dat$taxon[!is.na(C) & !(is.na(a_nit) & is.na(a_amm)  & is.na(a_phos))])
  aff = data.frame(r=NULL, A=NULL)
  for (i in taxons) {
    ix = dat$taxon == i
    points(r[ix], a_nit[ix], pch=16, col=col)
    points(r[ix], a_amm[ix], pch=17, col=col)
    points(r[ix], a_phos[ix], pch=18, col=col)
    
    #points(C[ix], a_nit[ix], pch=16, col=col)
    #points(C[ix], a_amm[ix], pch=17, col=col)
    #points(C[ix], a_phos[ix], pch=18, col=col)
    col = col + 1
    
    aff = rbind( aff, data.frame(r=r[ix], a= a_nit[ix]))
    aff = rbind( aff, data.frame(r=r[ix], a= a_amm[ix]))
    aff = rbind( aff, data.frame(r=r[ix], a= a_phos[ix]))
  }
  ix = !is.na(aff$r) & !is.na(aff$a)
  aff = aff[ix,]
  #points(C, Aphos,pch=18)
  legend(x="topright",bty="n",
         legend=taxons, pch=16, col=seq(1,length(taxons)))
  #
  # "Fit":
  # 
  r = 10^seq(-1,3,length.out = 100) # mum
  Diff = 1.5e-5*60*60*24 # cm^2/day, at 10 degrees (https://www.unisense.com/files/PDF/Diverse/Seawater%20&%20Gases%20table.pdf)
  rho = 0.40 #  g/cm3 Menden-Deuer (2000). If we assume that m propto V we get approximately rho = 0.4e-6 ugC/um^3. // 0.57 # g/cm3 (Andersen et al 2016; rho = m/V = 0.3*d^3/(4/3*pi*(d/2)^3) )
  aNmax = 3*Diff*1e8*1e-3/rho * 1e-6 # L/day/mugC
  rstar = 2; # mum
  corr = 1 - parameters()$c*calcMass(r)^(-1/3)
  #cat("alphaN_max = ", (ANmax/m^(1/3))[1], "\n")
  #lines(m, ANmax/m, lwd=1, lty=dotted)  
  #lines(m, 0.2*m/m, lwd=1, lty=dotted)  
  lines(r, aNmax*r^-2, lwd=1, lty=dotted)  
  #lines(m, 0.2*m/m, lwd=1, lty=dotted)  
  lines(r, aNmax*r^-2*(r/rstar)^2, lwd=1, lty=dotted, col="blue")
  lines(r, aNmax*r^-2/(1+((r)/rstar)^-2), lwd=2)
  
  #lines(m, p$AN/p$cN*m^(2/3)/m, lty=dotted, col="blue")
  #cat( mean(ANmax/(parameters()$AN*m^(1/3))) ,"\n" )
  cat("aN = ", aNmax, "r^-2 \n")
  cat("r_D^* = ", rstar," mum\n")
  
  # Make a bi-linear fit:
  #form = formula( log(C) ~ log(a/(1+b*C^-0.333))+0.3333*log(C) )
  #fit = nls(form, data=aff, 
  #          start = list(a=1e-4, b=0.05),
  #          lower = list(a=0,b=0), algorithm = "port", trace=TRUE)
  #lines(m, exp(predict(fit,list(C=m))))
  
  #V = 10^seq(-2,8,length.out = 100)
  #alphaN_Ward2018 = 1.1 * 1000 / 1000 / 12 # convert from m3->liter, from mmol->umol, from mol->g 
  #lines(convertVolume2Mass(V), alphaN_Ward2018*V^-0.35*convertVolume2Mass(V), col="blue", lty=dashed)
  
  # Camila's parameters:
  #lines(C, 3.75e-5*C^(1/3), col="orange")
  
  #alphaN_Banas2011 = 2.6/0.1 / 14 / 6 # convert from uM N-> ugN, from gN->gC
  #lines(m , alphaN_Banas2011*m^(1-0.45), col="blue", lty=dashdotted)
  
  legend(x="bottomleft", bty="n",
         legend=c("Model","Diffusion limit","Porter limit"),
         lty=c(solid,dotted,dotted),
         lwd=c(2,1,1,1),
         col=c("black","black","blue","orange"))
  
  #r = (3*dat$volume/(4*pi))^(1/3) # mu m
  #m = 0.3e-6*(2*r)^3
  #ANmax = 4*pi*Diff*(r*1e-4) * 1e-3
  #points(m, ANmax, col="red")
}



# plotAN = function() {
#   dat = read.csv("../data/Nutrient data from Edwards et al (2015b).csv",
#                  header=TRUE,sep=",")
#   #
#   # Convert carbon from mol to g:
#   #
#   dat$c_per_cell = dat$c_per_cell*12
#   #
#   # Convert to weights when only volume is given:
#   # 
#   ix = (!is.na(dat$volume)) & (is.na(dat$c_per_cell))
#   dat$c_per_cell[ix] = convertVolume2Mass(dat$volume[ix], dat$taxon[ix])
#   C = dat$c_per_cell
#   #
#   # Calc affinities
#   #
#   Anit = dat$vmax_nit / dat$k_nit
#   Aamm = dat$vmax_amm / dat$k_amm
#   Aphos = dat$vmax_p / dat$k_p
#   #
#   # Plot
#   #
#   defaultplot()
#   loglogpanel(xlim=C, ylim=c(1e-9,1e-3),
#               xlab="Mass ($\\mathrm{\\mu g_C}$)",
#               ylab="Nutrient affinity, $\\textit{A_N}$ (L/day)")
#   
#   col = 1
#   taxons = unique(dat$taxon[!is.na(C) & !(is.na(Anit) & is.na(Aamm)  & is.na(Aphos))])
#   aff = data.frame(C=NULL, A=NULL)
#   for (i in taxons) {
#     ix = dat$taxon == i
#     points(C[ix], Anit[ix], pch=16, col=col)
#     points(C[ix], Aamm[ix], pch=17, col=col)
#     points(C[ix], Aphos[ix], pch=18, col=col)
#     col = col + 1
#     
#     aff = rbind( aff, data.frame(C=C[ix], A= Anit[ix]))
#     aff = rbind( aff, data.frame(C=C[ix], A= Aamm[ix]))
#     aff = rbind( aff, data.frame(C=C[ix], A= Aphos[ix]))
#   }
#   ix = !is.na(aff$C) & !is.na(aff$A)
#   aff = aff[ix,]
#   #points(C, Aphos,pch=18)
#   legend(x="bottomright",bty="n", 
#          legend=taxons, pch=16, col=seq(1,length(taxons)))
#   
#   
#   m = 10^seq(-7,1) # mu gC
#   r = 1.5/2*(1e-6*m)^(1/3) # cm
#   Diff = 1.5e-5*60*60*24 # cm^2/day, at 10 degrees (https://www.unisense.com/files/PDF/Diverse/Seawater%20&%20Gases%20table.pdf)
#   ANmax = 4*pi*Diff*r*1e-3 # L/day
#   cat("alphaN_max = ", (ANmax/m^(1/3))[1], "\n")
#   lines(m, ANmax, lwd=1, lty=dotted)  
#   lines(parameters()$m, parameters()$ANm, lwd=2)
#   lines(m, p$AN/p$cN*m^(2/3), lty=dotted, col="blue")
#   #cat( mean(ANmax/(parameters()$AN*m^(1/3))) ,"\n" )
#   cat("alphaN = ", mean(parameters()$AN*m^(1/3)), "\n")
#   
#   # Make a bi-linear fit:
#   #form = formula( log(C) ~ log(a/(1+b*C^-0.333))+0.3333*log(C) )
#   #fit = nls(form, data=aff, 
#   #          start = list(a=1e-4, b=0.05),
#   #          lower = list(a=0,b=0), algorithm = "port", trace=TRUE)
#   #lines(m, exp(predict(fit,list(C=m))))
#   
#   #V = 10^seq(-2,8,length.out = 100)
#   #alphaN_Ward2018 = 1.1 * 1000 / 1000 / 12 # convert from m3->liter, from mmol->umol, from mol->g 
#   #lines(convertVolume2Mass(V), alphaN_Ward2018*V^-0.35*convertVolume2Mass(V), col="blue", lty=dashed)
#   
#   # Camila's parameters:
#   #lines(C, 3.75e-5*C^(1/3), col="orange")
#   
#   #alphaN_Banas2011 = 2.6/0.1 / 14 / 6 # convert from uM N-> ugN, from gN->gC
#   #lines(m , alphaN_Banas2011*m^(1-0.45), col="blue", lty=dashdotted)
#   
#   legend(x="topleft", bty="n",
#          legend=c("Model","Diffusion limit","Porter limit"),
#          lty=c(solid,dotted,dotted),
#          lwd=c(2,1,1,1),
#          col=c("black","black","blue","orange"))
#   
#   #r = (3*dat$volume/(4*pi))^(1/3) # mu m
#   #m = 0.3e-6*(2*r)^3
#   #ANmax = 4*pi*Diff*(r*1e-4) * 1e-3
#   #points(m, ANmax, col="red")
# }
# 
# plotANsimple = function() {
#   dat = read.csv("../data/Nutrient data from Edwards et al (2015b).csv",
#                  header=TRUE,sep=",")
#   #
#   # Convert carbon from mol to g:
#   #
#   dat$c_per_cell = dat$c_per_cell*12
#   #
#   # Convert to weights when only volume is given:
#   # 
#   ix = (!is.na(dat$volume)) & (is.na(dat$c_per_cell))
#   dat$c_per_cell[ix] = convertVolume2Mass(dat$volume[ix], dat$taxon[ix])
#   C = dat$c_per_cell
#   #
#   # Calc affinities
#   #
#   Anit = dat$vmax_nit / dat$k_nit
#   Aamm = dat$vmax_amm / dat$k_amm
#   Aphos = dat$vmax_p / dat$k_p
#   #
#   # Plot
#   #
#   defaultplot()
#   loglogpanel(xlim=C, ylim=c(1e-9,1e-3),
#               xlab="Mass ($\\mathrm{\\mu g_C}$)",
#               ylab="Nutrient affinity, $\\textit{A_N}$ (L/day)")
#   
#   col = 1
#   taxons = unique(dat$taxon[!is.na(C) & !(is.na(Anit) & is.na(Aamm)  & is.na(Aphos))])
#   for (i in taxons) {
#     ix = dat$taxon == i
#     points(C[ix], Anit[ix], pch=16, col="darkblue")
#     points(C[ix], Aamm[ix], pch=17, col="darkblue")
#     #points(C[ix], Aphos[ix], pch=18, col=col)
#     col = col + 1
#   }
#   #points(C, Aphos,pch=18)
#   #legend(x="bottomright",bty="n", 
#   #       legend=taxons, pch=16, col=seq(1,length(taxons)))
#   
#   
#   m = 10^seq(-7,1) # mu gC
#   r = 1.5/2*(1e-6*m)^(1/3) # cm
#   Diff = 1.5e-5*60*60*24 # cm^2/day, at 10 degrees (https://www.unisense.com/files/PDF/Diverse/Seawater%20&%20Gases%20table.pdf)
#   ANmax = 4*pi*Diff*r*1e-3 # L/day
#   cat("alphaN_max = ", (ANmax/m^(1/3))[1], "\n")
#   lines(m, ANmax, lwd=1, lty=dotted)  
#   lines(m, 1e-3*m^(2/3))
#   #lines(m, parameters()$AN*m^(1/3), lwd=2)
#   #cat( mean(ANmax/(parameters()$AN*m^(1/3))) ,"\n" )
#   cat("alphaN = ", mean(parameters()$AN*m^(1/3)), "\n")
#   
#   V = 10^seq(-2,8,length.out = 100)
#   alphaN_Ward2018 = 1.1 * 1000 / 1000 / 12 # convert from m3->liter, from mmol->umol, from mol->g 
#   #lines(convertVolume2Mass(V), alphaN_Ward2018*V^-0.35*convertVolume2Mass(V), col="blue", lty=dashed)
# }

# plotLimitation = function(p=parameters()) {
#   m = p$m
#   
#   DOC = 2
#   N = 1
#   L = 60
#   B = 10
#   r = calcRates(L, N, DOC, rep(B, p$n), p)
#   
#   defaultplot()
#   semilogxpanel(xlim=m, ylim=c(0, 1.5))
#   
#   lines(m, r$JN/m, col="blue")
#   lines(m, r$JDOC/m, col="brown")
#   lines(m, r$JL/m, col="darkgreen")
#   lines(m, r$JF/m, col="red")
#   lines(m, r$Jmax/m, col="black")
#   lines(m, r$Jtot/m, col="black")
# }

plotRstar = function() {
  p = parameters(n=100, mmin10=-9)
  
  m = p$m
  r = p$r
  #nu = pmin(1,p$c * m^(-1/3))
  mu = c(0,0.4)#c(0.2, 0.4, 0.6)
  
  tightaxes()
  defaultplot()
  loglogpanel(xlim=r, 
              ylim=c(0.0001,1000),
              xlab = "Cell radius ($\\mu$m)",
              ylab = "Limiting concentration ($\\mu$M) or ($\\mu$E/m$^2$/s)")
  
  
  rmin = -((p$cLeakage + 3*p$alphaJ*p$delta)/(p$jR - p$alphaJ))
  polygon(c(0.01, rmin, rmin,0.011), c(0.0001,0.00011,2000,2000), 
          col=rgb(0.5,0.5,0.5,alpha=0.5), border=NA)
  
  
  lty=1
  lwd=2
  for (i in 1:length(mu)) {
    #correct = mu[i]*(1-nu)/((1-nu)-mu[i])
    #Nstar = m^0.667 /(p$AN*p$rhoCN) * correct
    
    #correct = p$cLeakage + p$alphaJ*mu[i]*(1-nu)/(p$alphaJ*(1-nu)-mu[i])
    #Nstar = m^(2/3)/(p$AN*p$rhoCN) * correct
    Nstar = ((p$cLeakage + mu[i]*r)*(r^2 + p$rNstar^2))/(r*p$aN)
    Nstar[Nstar<0] = Inf
    lines(r, Nstar/14, lwd=lwd, lty=lty,col="blue")
    
    DOCstar = ((p$cLeakage + (p$jR+mu[i])*r)*(r^2 + p$rNstar^2))/(r*p$aN)
    DOCstar[DOCstar<0] = Inf
    lines(r, DOCstar/12, lwd=lwd, lty=lty,col="magenta")
    
    #Lstar = m^0.333 * p$alphaJ/(p$AL) * mu[i]/(p$alphaJ*(1-nu)-mu[i])
    #Lstar = p$alphaJ*(p$AL+p$cL*m^(1/3))*mu[i] / (p$AL*p$cL*(p$alphaJ-mu[i])*(1-nu))
    #Lstar = m^(1/3)/(p$AL*p$epsilonL) * 1/(1 -exp(-p$cL*m^(1/3))) * correct
    Lstar = exp(r/p$rL)*(p$cLeakage + (p$jR + mu[i])*r) / ((exp(r/p$rL)-1)*p$alphaL)
    Lstar[Lstar<0] = Inf
    lines(r , Lstar, lwd=lwd, lty=lty,col="darkgreen")
    
    
    Fstar = -((p$cF*(p$cLeakage + (p$jR + mu[i])*r))/(p$aF*r*(p$cLeakage + (p$jR + mu[i])*r - p$cF*p$epsilonF)))
    Fstar[Fstar<0] = Inf
    lines(r, Fstar/12, lwd=lwd, lty=lty, col="darkred")
    
    lty=dotted
    lwd=1
  }
  
  legend(x="bottomright", bty="n",
         legend=c(TeX("$\\textit{N}^*$"), 
                  TeX("$$DOC^*$"), 
                  TeX("$$\\textit{L}^*$"),
                  TeX("$$\\textit{F}^*$")),
         col=c("blue","magenta","darkgreen","darkred"),
         lwd=2)
}

calcMufactor = function(p) {
  f0 = 0.6
  delta = (p$m[2]-p$m[1]) / sqrt(p$m[2]*p$m[1])
  cat("Delta = ", delta, "\n")
  return( (1-f0) * p$AF *sqrt(2*pi)*p$sigma )
}

panelStrategies = function(p,N,L,Bsheldon,DOC,y ,
                           ylabel,xlabel='',bXaxis=FALSE) {
  loadNUMmodel()
  # Set parameters
  dummy = .Fortran("f_setupgeneralistsonly", as.integer(p$n))
  dummy = .Fortran("f_sethtl", 
                   as.double(p$mHTL),
                   as.double(p$mortHTL), 
                   as.logical(FALSE), as.logical(FALSE))
  p$tEnd = 1;
  
  col = matrix(0,nrow=length(L), ncol=p$n)
  for (i in 1:length(L)) {
    p$L = L[i]
    p$N0 = N[i]
    p$B0 = Bsheldon[,i]
    p$DOC0 = DOC[i]
    p$tEnd = 1
    
    derivativeF(t,c(p$N0,p$DOC0,p$B0),p)
    rates = getFrates(p)
    #sim = simulateChemostat(p)
    #rates = sim$rates
    m = p$m
    
    strategy = calcStrategy(p,rates)
    for (j in 1:p$n) {
      if (strategy[j]=="Heterotroph")
        col[i,j] = 1
      if (strategy[j]=="Mixotroph")
        col[i,j] = 2
      if (strategy[j]=="Light limited")
        col[i,j] = 3
      if (strategy[j]=="Nutrient limited")
        col[i,j] = 4
      if (strategy[j]=="Osmoheterotroph")
        col[i,j] = 5
      if (p$B0[j]/log(p$Delta) < 1e-4)
        col[i,j]=6
    }
  }
  # col = col-0.5
  alpha = 0.5
  color = list()
  color = c(rgb(1,0,0,alpha=alpha),
            rgb(1,0.5,0.5,alpha=alpha),
            rgb(0,1,0,alpha=alpha),
            rgb(0,0,1,alpha=alpha),
            rgb(165/256,42/256,42/256,alpha=alpha),
            rgb(1,0,1,alpha=alpha))
  
  #defaultplot()
  #semilogxpanel(xlim=p$m, ylim=L,xlab="x")
  #defaultpanel(xlim=p$m,ylim=c(y[1],y[n]))
  image(p$m, y, t(col), log="xy",col=color,
        ylim=c(y[1],y[length(y)]), zlim=c(0,6),
        xlab="",ylab="",xaxt="n", yaxt="n")
  loglogpanel(xlim=p$m, ylim=c(y[1],y[length(y)]),
              ylab=ylabel, xlab=xlabel,
              xaxis=bXaxis,
              new=TRUE)
  makepanellabel()
  #image(p$m, y, t(col), log="xy", col=color)
}

plotStrategies = function(n=50) {
  p = parametersChemostat(parameters(n=n))
  p$tEnd = 1
  sim = simulateChemostat(p) # ...just to get the library loaded
  #
  # Base parameters:
  #
  p$N0 = 1
  p$L = 30
  p$B0 = rep(1,p$n)/log(p$Delta)
  p$DOC0 = 1
  #
  # Variations:
  #
  L = 10^seq(1,2,length.out=n)
  Bsheldon = t((10^seq(0,2,length.out=n)) %*% t(rep(1,length.out=p$n)))
  N = 10^seq(-1,1,length.out=n)
  DOC = 10^seq(-1,1,length.out=n)
  ones = rep(1,n)
  #
  # Plot of gains:
  #
  defaultplot(mfcol=c(5,1))
  
  derivativeF(t,c(p$N0,p$DOC0,p$B0),p)
  rates = getFrates(p)
  
  mm = 10^seq(-8,2, length.out = 100)  
  r = rates
  
  ylim = c(0,1.6)
  semilogxpanel(xlim=p$m, ylim=ylim,
                xlab="",
                ylab="Uptakes (1/day)",
                xaxis=FALSE)
  calcStrategy(p,r,bPlot=TRUE)
  
  lines(p$m, r$jLreal, lwd=2, col="green")
  lines(p$m, r$jN, lwd=2, col="blue")
  lines(p$m, r$jDOC, lwd=2, col="brown")
  lines(p$m, r$jF,lwd=2,col="red")
  makepanellabel()
  
  # legend(x="topright", cex=cex,
  #         legend=c("Light harvesting","Nutrient uptake","DOC uptake","Food consumption"),
  #         col=c("green","blue","brown","red"),
  #         lwd=c(2,2,2,2),
  #         bty="n")
  #
  # The four strategy panels:
  #
  panelStrategies(p,N,p$L*ones,p$B0 %*% t(ones),p$DOC0*ones,y=N,
                  'Nutr. ($ \\mu g_N/l)$')
  hline(p$N0)
  panelStrategies(p,p$N0*ones,L,p$B0 %*% t(ones),p$DOC0*ones,y=L,
                  'Light ($\\mu E/m^2/s)$')
  hline(p$L)
  panelStrategies(p,p$N0*ones,p$L*ones,Bsheldon,p$DOC0*ones,y=Bsheldon[1,],'B_{Sheldon}')
  hline(p$B0[1])
  panelStrategies(p,p$N0*ones,p$L*ones,p$B0 %*% t(ones),DOC,y=DOC,
                  'DOC ($\\mu g_C/l)$','Cell mass ($\\mu g_C$)',
                  bXaxis=TRUE)
  hline(p$DOC0)
}



plotSimulationExamples = function() {
  
  plotSpectrum <- function(sim, t=max(sim$t), 
                           bXlabel=FALSE, bLegend=TRUE, bStatsTop=FALSE) {
    p = sim$p
    m = p$m
    r = sim$rates
    
    ylim = c(0.1,200)
    
    ixt = which(floor(sim$t)==t+365)[1]
    B = sim$B / log(p$Delta)
    N = sim$N
    DOC = sim$DOC
    r = sim$rates
    
    if (bXlabel)
      xlab = "Carbon mass ($\\mu$gC)"
    else
      xlab = ""
    loglogpanel(xlim=p$m, ylim=ylim, xaxis = bXlabel,
                xlab=xlab,
                ylab="Sheldon biomass ($\\mu$gC/l)")
    
    lines(m, B, lwd=4)
    if (p$n<15)
      points(m,B)
    #
    # Theoretical prediction:
    #
    kappa = calcSheldonKappa(sim$p)
    lines(m, rep(kappa, p$n), lty=dotted, col=grey(0.5))
    text(x = 8e-8, y = 0.8*kappa,
         labels='Theory',
         cex=0.75*cex, pos=3,col=grey(0.5))
    #     mar=c(4,5,8,2)+0.1)
    #
    # Add gray-scale variation
    #
    if (p$latitude==0)par()
    polygon(c(m, m[seq(p$n,1,by = -1)]), c(sim$Bmin, sim$Bmax[seq(p$n,1,by = -1)])/ log(p$Delta), 
            col=rgb(0.5,0.5,0.5,alpha=0.25), border=NA)
    
    # Determine limiting process:
    strategy = calcStrategy(p,r,bPlot=TRUE)
    #
    # Add extra size labels
    #
    d = 10^seq(-6,-1,by=1)
    axis(side=3, line=0,
         at=0.3e6*d^3,
         labels = c("0.01","0.1","1","10","100","1e3"))
    mtext(TeX("Diameter ($\\mu$m)"), side=3, line=0.75, at=1e-3, adj=1,cex=cex)
    #
    # Legend:
    #
    if (bLegend)
      legend(xpd=NA, x="topright", inset=c(-.75, 0), 
             bty="n", 
             legend=c("Osmoheterotrophs", "Light limited phototrophs","N limited phototrophs","Mixotrophs","Heterotrophs"),
             fill=c(colOsmo, colPhoto,colN,colMixo,colHetero,"transparent"),
             border=c("black","black","black","black","black","transparent"),
             lwd = c(0,0,0,0,0,3),
             col=c(NA,NA,NA,NA,NA,1))
    #
    # Summary state variables: 
    #
    sx = "bottomleft"
    if (bStatsTop)
      sx = "topleft"
    legend(x=sx, cex=0.8, bty="n", text.col=grey(0.5),
           y.intersp = 0.7,
           legend=c(TeX(sprintf("DIN: %2.2f $\\mu$M", N/14)) ,
                    TeX(sprintf("DOC: %2.2f $\\mu$M", DOC/12))))
    
    #func = calcFunctionsChemostat(sim$p, sim$rates, sim$p$L, sim$N, sim$B)
    func = calcFunctions(sim)
    sx = "bottomright"
    if (bStatsTop)
      sx = "topright"
    legend(x=sx, cex=0.8, bty="n", text.col=grey(0.5),
           y.intersp = 0.7,
           legend=c(
             TeX(sprintf("Chl-a: %1.3f $gC/m^2$", func$Chl_per_m2/1000)),
             TeX(sprintf("Pico: %1.2f $gC/m$^2$", func$Bpico)),
             TeX(sprintf("Nano: %1.2f $gC/m$^2$", func$Bnano)),
             TeX(sprintf("Micro: %1.2f $gC/m$^2$", func$Bmicro)))
    )
    
    makepanellabel()
    
    box()
  }
  
  plotRates = function(sim, p=sim$p, 
                       B=sim$B, N=sim$N, DOC=sim$DOC,
                       t=max(sim$t), 
                       bLosses=TRUE, bXaxis=TRUE) {
    mm = 10^seq(-8,2, length.out = 100)  
    
    L = p$L
    r = sim$rates
    
    ylim =c(-1.6,1.6)
    if (max(sim$rates$jTot) < 0.5)
      ylim = c(-0.5,0.5)
    
    semilogxpanel(xlim=p$m, ylim=ylim,
                  xlab="Carbon mass ($\\mu$gC)",
                  ylab="Rates (1/day)",
                  xaxis=bXaxis)
    calcStrategy(p,r,bPlot=TRUE, ylim=ylim)
    #
    # Gains
    #
    lines(p$m, r$jTot, lwd=4, type="l", col="black")# log="x", xlim=range(p$m),
    lines(p$m, r$jMax, lty=3)
    
    #lines(mm, p$AL*mm^(2/3)*input$L/mm, lty=3, lwd=1, col="green")
    #JLreal = r$Jtot - r$JF+p$Jresp-r$JDOC
    #lines(p$m, p$ALm*p$L/p$m, lty=dotted, lwd=1, col="green")
    lines(p$m, r$jL , lty=dotted, lwd=1, col="green")
    lines(p$m, r$jLreal, lwd=2, col="green")
    
    #lines(mm, p$AN*mm^(1/3)*p$rhoCN*sim$N/mm, lwd=1, lty=3, col="blue")
    #lines(p$m, p$Jmax * p$ANm*N / (p$Jmax/p$rhoCN + p$ANm*N)/p$m, lwd=1, lty=dotted, col="blue")
    #lines(p$m, r$JN/p$m*p$rhoCN, lwd=4, col="blue")
    #lines(p$m, p$JN/p$m, lwd=1, lty=dotted, col="blue")
    lines(p$m, r$jN, lwd=2, col="blue")
    
    #lines(mm, p$AN*mm^(1/3)*DOC/mm, lwd=1, lty=3, col="brown")
    lines(p$m, r$jDOC, lwd=2, col="brown")
    
    lines(p$m, r$jFreal,lwd=2,col="red")
    lines(p$m, p$epsilonF * r$jFmax ,col="red", lty=dotted)
    
    legend(x="topright", cex=cex, 
           xpd=NA, inset=c(-0.71,-.10),
           y.intersp = 0.7,
           legend=c("Light harvesting","Nutrient uptake","DOC uptake","Food consumption","Division rate"),
           col=c("green","blue","brown","red","black"),
           lwd=2,
           bty="n")
    #
    # Losses
    #
    if (bLosses) {
      #polygon(c(1e-9,10,10,1e-9), c(-1.5,-1.5,0,0), 
      #        col=rgb(1,0,0,alpha=0.25), border=NA)
      
      JNexude = r$JNloss 
      lines(p$m, -(r$mortpred + p$mortHTL + r$mort2 + p$mort), lwd=5)
      lines(p$m, -r$mortpred, col="red", lwd=2)
      lines(p$m, -r$mortHTL, col="magenta", lwd=2)
      lines(p$m, -r$mort2, col="orange", lwd=2)
      lines(p$m, -r$jR, col="grey", lwd=2)
      lines(p$m, -r$jLossPassive, col="darkgreen", lwd=2)
      
      #BSheldon =exp(mean(log(B)))
      #delta = (p$m[2]-p$m[1]) / sqrt(p$m[2]*p$m[1])
      #mortPredTheoretical = BSheldon * (1-0.6) * p$AF *sqrt(2*pi)*p$sigma / delta
      #lines(range(p$m), -mortPredTheoretical*c(1,1), lty=dotted, col="red")
      
      legend(x="bottomright", cex=cex,
             xpd=NA, inset=c(-0.75,-0.1),
             y.intersp = 0.7,
             legend=c("","", "Predation", "Virulysis", 
                      "Higher trophic levels","Respiration","Passive"),
             col=c(NA,NA,"red", "orange", "magenta","grey","darkgreen"),
             lwd=c(0,0,2,2,2,2,2), bty="n")
      
      lines(p$m, 0*p$m, col="white", lwd=2)
      lines(p$m, 0*p$m, lty=dashed, lwd=1)
      
      makepanellabel()
    }
    return(r)
  }
  
  
  plotSim = function() {
    defaultplot(c(2,1))
    par(cex.axis=cex,
        cex.lab=cex,
        oma=c(0, 0, 2.2, 10) + 0.1)
    
    if (sim$p$d < 0.01)
      bTop = TRUE
    else
      bTop = FALSE
    
    plotSpectrum(sim, bXlabel=FALSE, bStatsTop = bTop)
    plotRates(sim)
    #
    # Write functions:
    #
    func = calcFunctions(sim)
    cat('ProdGross' , func$ProdGross, 'gC/m2/yr\n')
    cat('ProdNet  ' , func$ProdNet, 'gC/m2/yr\n')
    cat('ProdHTL  ' , func$ProdHTL, 'gC/m2/yr\n')
    cat('ProdNew  ' , func$ProdNew, 'gC/m2/yr\n')
    cat('ProdBact ' , func$ProdBact, 'gC/m2/yr\n')
    
    cat('e_PP     ', func$ePP,'\n')
    cat('e_HTL    ', func$eHTL,'\n')
    cat('e_Bact   ', func$eBact,'\n')
    cat('-----------------------------------------------\n')
  }
  
  #
  # Oligotrophic
  #
  p = parametersChemostat()
  p$d = dOligotrophic
  p$L = LOligotrophic
  sim = simulateChemostat(p)
  pdfplot("../simOligo.pdf", plotSim,width = doublewidth, height=1.7*height)
  #
  # Eutrophic
  #
  p$d = dEutrophic
  p$L = LEutrophic
  sim = simulateChemostat(p)
  pdfplot("../simEutro.pdf", plotSim,width = doublewidth, height=1.7*height)
}


# plotVaryDiffusion = function(p=parametersChemostat(),
#                              L=100,
#                              d=10^seq(-6,-3,length.out=10)) {
#   
#   B = matrix(data=0, nrow=p$n, ncol=length(d))
#   for (i in 1:length(d)) {
#     p$d=d[i]
#     sim = simulateChemostat(p)
#     B[,i] = sim$B
#   }
#   image(log10(d),log10(p$m),t(B),
#         xlab="log10(mixing)", ylab="log10(cell mass)")
# }
# 
# plotVaryLightAndDiffusion = function(n=10) {
#   library(lattice)
#   library(latticeExtra)
#   library(reshape)
#   
#   d = seq(0.004,.5,length.out=n) #10^seq(-2,log10(2),length.out = 10)
#   L = seq(10,200,length.out=n)
#   #
#   # Simulations
#   #
#   df = data.frame(size=NULL, d=NULL, L=NULL, biomass=NULL)
#   df2 = data.frame(func=NULL, d=NULL, L=NULL, z=NULL)
#   #tmp = matrix(ncol=4, nrow=0)
#   #B = array(dim=c(3,length(d), length(L)))
#   for (i in 1:length(d))
#     for (j in 1:length(L)) {
#       p = parametersChemostat()
#       p$d = d[i]
#       p$L = L[j]
#       print(p$d)
#       print(p$L)
#       #p$mHTL = 0.008
#       sim = simulateChemostat(p, useF=TRUE)
#       func = calcFunctionsChemostat(sim$p, sim$rates, sim$N, sim$B)
#       
#       df = rbind(df, data.frame(size="Pico",d=d[i],L=L[j],biomass=func$Bpico))
#       df = rbind(df, data.frame(size="Nano",d=d[i],L=L[j],biomass=func$Bnano))
#       df = rbind(df, data.frame(size="Micro",d=d[i],L=L[j],biomass=func$Bmicro))
#       #df = rbind(df, data.frame(size="Chl",d=d[i],L=L[j],biomass=func$Chl_per_l))
#       
#       #tmp = rbind(tmp, c(func$Bpico, func$Bnano, func$Bmicro, func$Chl_per_l))
#       
#       df2 = rbind(df2, data.frame(func="Biomass", d=d[i],L=L[j],
#                                   z=func$Bpico+func$Bnano+func$Bmicro))
#       df2 = rbind(df2, data.frame(func="N", d=d[i],L=L[j], z=sim$N))
#       df2 = rbind(df2, data.frame(func="DOC", d=d[i],L=L[j], z=sim$DOC))
#       
#       #B[1,i,j] = func$Bpico+func$Bnano+func$Bmicro
#       #B[2,i,j] = sim$N
#       #B[3,i,j] = sim$DOC
#     }
#   #
#   # Plots
#   #
#   plt=levelplot(biomass ~ d*L | size, data=df,
#                 aspect=1, region=TRUE, 
#                 col.regions = terrain.colors(100),
#                 #scales = list(x = list(log = 10, at=c(0.01,0.1,1))),
#                 #xscale.components = xscale.components.logpower,#10ticks(n=5,n2=10),
#                 layout=c(3,1),
#                 ylab=TeX("Light (${\\mu}E m^{-2}s^{-1}$)"), xlab=TeX("Mixing rate (day$^{-1}$)"))
#   
#   plt2 = levelplot(z ~ d*L | func, data=df2,
#                    aspect=1, region=TRUE,
#                    col.regions = terrain.colors(100),
#                    layout=c(3,1),
#                    ylab=TeX("Light (${\\mu}E m^{-2}s^{-1}$)"), 
#                    xlab=TeX("Mixing rate (day$^{-1}$)"))
#   
#   return(plt)
# }

calcFunctions = function(sim) {
  ProdGross = 0
  ProdNet = 0
  ProdHTL = 0
  eHTL = 0
  Bpico = 0
  Bnano = 0
  Bmicro = 0
  x= 0
  func = .Fortran("f_getfunctions", 
                  u = as.numeric( c(sim$N, sim$DOC, sim$B) ),
                  ProdGross=as.numeric(1), 
                  ProdNet=as.numeric(1),
                  ProdHTL=as.numeric(1),
                  ProdBact=as.numeric(1),
                  eHTL=as.numeric(1),
                  Bpico=as.numeric(1),
                  Bnano=as.numeric(1),
                  Bmicro=as.numeric(1))
  f = function(x) x*sim$p$M
  for (j in c(2,3,4,5,7,8,9))
    func[j] = lapply(func[j],f)
  func$ProdNew = 365*sim$p$M*1e-6*1000*sim$p$d*(sim$p$N0-sim$N) * sim$p$rhoCN
  func$ePP = func$ProdNet/func$ProdGross
  func$eBact = func$ProdBact / func$ProdNet
  func$d = sim$p$d
  func$N = sim$N
  func$DOC = sim$DOC
  # Production rate in 1/day
  func$ProdRate = func$ProdNet / (func$Bpico+func$Bnan+func$Bmicro) / 365
  #
  # Chl-a: Use rough conversion from Edwards et al (2015) that Chl-a propto alpha
  #
  func$Chl_per_l = sum( sim$rates$jLreal/(sim$p$epsilonL*sim$p$L)*sim$B ) # ug/l = mg/m3
  func$Chl_per_m2 = sim$p$M * func$Chl_per_l*0.001 # Convert to gChl/m2
  
  
  return(func)
}

plotFunctions = function(L=c(20, 100), n=10) {
  
  panelsFunctions = function(L, n=10, yaxis=TRUE) {
    d = 10^seq(-3,0,length.out = n) #seq(0.02,2,length.out=n) #
    
    p = parametersChemostat()
    p$L = L
    
    F = data.frame()
    for (i in 1:length(d)) {
      print(d[i])
      p$d = d[i]
      sim = simulateChemostat(p, useF=TRUE)
      #func = calcFunctionsChemostat(sim$p, sim$rates, sim$N, sim$B)
      func = calcFunctions(sim)
      F = rbind(F, as.data.frame(func))
    }
    
    ylab = ""
    # Biomass
    B = F$Bpico+F$Bnano+F$Bmicro
    if (yaxis)
      ylab = "Biomass ($g_C/m^2$)"
    semilogxpanel(xlim=c(d), ylim=c(0,15), xaxis=FALSE, yaxis=yaxis,
                  ylab=ylab)
    lines(F$d, B, lwd=2)
    lines(F$d, F$Bpico, lwd=0.5, col="grey")
    lines(F$d, F$Bnano, lwd=1, col="grey")
    lines(F$d, F$Bmicro, lwd=1.5, col="grey")
    makepanellabel()
    
    
    # Concentrations:
    if (yaxis)
      ylab = "Concentrations"
    loglogpanel(xlim=d, ylim=c(0.01,1000), xaxis=FALSE, yaxis=yaxis,
                ylab=ylab)
    lines(F$d, F$N/14, lwd=2, col="blue")
    lines(F$d, F$DOC/12, lwd=2, col="brown")
    lines(F$d, F$Chl_per_l, lwd=2, col="darkgreen")
    if (yaxis) {
      legend("topleft",bty="n", cex=cex,
             legend=c(TeX("N ($\\mu$M)"),TeX("DOC ($\\mu$M)")),lwd=2,col=c("blue", "brown","darkgreen"))
      legend("topright",bty="n", cex=cex,
             legend=TeX("Chl ($\\mu$g/l)"),lwd=2,col="darkgreen")
    }
    makepanellabel()
    #DOC
    #if (yaxis)
    #  ylab = "DOC  ($\\mu$M)"
    #loglogpanel(xlim=d, ylim=c(1,5000), xaxis=FALSE, yaxis=yaxis, 
    #            ylab=ylab)
    #lines(F$d, 1000*F$DOC/12, lwd=2)
    
    #Production:
    if (yaxis)
      ylab = "Prod. ($g_C/m^2$/yr)"
    loglogpanel(xlim=d, ylim=c(10,5000), xaxis=FALSE, yaxis=yaxis, 
                ylab=ylab)
    lines(F$d, F$ProdGross, lwd=1.5)
    lines(F$d, F$ProdNet, col="darkgreen", lwd=1.5)
    lines(F$d, F$ProdNew, col="blue", lwd=1.5)
    lines(F$d, F$ProdHTL, col="red", lwd=1.5)
    lines(F$d, F$ProdBact, col="magenta", lwd=1.5)
    if (yaxis) {
      legend("topleft",bty="n", cex=cex,
             y.intersp = 0.8,
             legend=c("Gross","Net"),
             lwd=2,col=c("black", "darkgreen","blue","red"))
      legend("topright",bty="n", cex=cex,
             y.intersp = 0.8,
             legend=c("New","HTL","Bact"),
             lwd=2,col=c("blue","red",'magenta'))
    }
    makepanellabel()
    
    # Eff
    if (yaxis)
      ylab = "Efficiencies"
    semilogxpanel(xlim=d, ylim=c(0,1.5),
                  ylab=ylab,
                  xlab="Mixing rate (day$^{-1}$)",
                  yaxis=yaxis, xaxis=FALSE)
    lines(F$d, F$eHTL, lwd=2, col="red")
    lines(F$d, F$ePP, lwd=2, col="darkgreen")
    lines(F$d, F$ProdBact/F$ProdNet, lwd=2, col="magenta")
    hline(1)
    if (yaxis) {
      legend("topright", ,bty="n", cex=cex,
             y.intersp = 0.8,
             legend=c(TeX("$\\epsilon_{HTL}$"), TeX("$\\epsilon_{PP}$"), TeX("$\\epsilon_{bact}$")),
             lwd=2, col=c("red", "darkgreen","magenta"))
    }
    makepanellabel()
    
    # Production rate
    if (yaxis)
      ylab = "Prod. rate (day$^{-1}$)"
    semilogxpanel(xlim=d, ylim=c(0,0.5),
                  xaxis=TRUE, yaxis=yaxis,
                  ylab=ylab)
    lines(F$d, F$ProdRate, lwd=2)
    makepanellabel()
    
    return(F)
  }
  
  defaultplotvertical(mfcol=c(5,2))
  
  func= panelsFunctions(L=L[1],n=n)
  
  panelsFunctions(L=L[2], yaxis=FALSE,n=n)
}


# plotPreferences = function(p = parameters()) {
#   defaultplot()
#   semilogxpanel(xlim=p$m, xlab="Cell weight ($\\mu$gC)",
#                 ylim=c(0,1.1), ylab="Preference")
#   
#   m = 10^seq(-10,1,length.out = 1000)
#   for (i in 1:p$n) {
#     lines(m, phi(p$m[i]/m,p$beta,p$sigma), col=grey(0.8))
#     lines(p$m, phi(p$m[i]/p$m,p$beta,p$sigma), lty=dotted)
#     points(p$m, phi(p$m[i]/p$m,p$beta,p$sigma),pch=16)
#   }
#   
# }

plotGridPreference = function(p = parameters()) {
  m = 10^seq(-6,1,length.out = 100)
  defaultplot()
  
  semilogxpanel(xlim=c(1e-6,1), xlab="Prey mass:predator mass",
                ylim=c(0,1), ylab="Preference")
  
  y = lapply(m, function(m) Phi(1/m, p, Delta=11))
  lines(m, y, lwd=thick, col=stdgrey)
  y = lapply(m, function(m) Phi(1/m, p, Delta=80))
  lines(m, y, lwd=thick)
  lines(m, phi(1/m,beta=p$beta,sigma=p$sigma), lty=dotted)
  
  legend("topleft", bty="n",
         legend=TeX(c("n=$\\infty, \\Delta$=1",
                      "n=10, $\\Delta$=11",
                      "n=6, $\\Delta$=80")),
         lty=c(dotted,solid,solid),
         lwd=c(1,thick,thick),
         col=c(black,stdgrey,black)
  )
  # sigma = function(Delta) {
  #   y = as.numeric( lapply(m, function(m) fTot(1,m,Delta)) )
  #   n0 = trapz(log(m),y)
  #   n1 = trapz(log(m), log(m)*y)/n0
  #   n2 = trapz(log(m), (log(m)-n1)^2*y)/n0
  #   return( sqrt(n2) )
  # }
  # 
  # semilogxpanel( xlim=Delta, xlab="Grid expansion factor $\\Delta$",
  #                ylim=c(1,2.5), ylab="Width $\\sigma$")
  # Delta = 10^seq(0,2,length.out=20)
  # y = lapply( Delta, sigma )
  # lines(Delta, y, lwd=thick)
  # lines(Delta, Delta/Delta*sigma, lty=dotted)
  # lines(c(10,10), c(1,3), lty=dotted)
}

plotGridtest = function() {
  n = c(6, 8, 10, 15, 25, 50)
  
  defaultplot()
  loglogpanel( xlim=c(5e-9,5000), xlab="Cell weight ($\\mu$gC)",
               ylim=c(2,200), ylab="Normalized biomass ($\\mu$gC/l)")
  
  textLegend = NULL
  for (i in 1:length(n)) {
    p = parametersChemostat( parameters(n[i]) )
    sim = simulateChemostat(p)
    Delta = p$m[2]/p$m[1]
    textLegend = c(textLegend, TeX(sprintf("n = %2.0f, $\\Delta$ = %2.1f", n[i],Delta)))
    
    Bnorm = sim$B/(sqrt(Delta)-1/sqrt(Delta))
    Bnorm[ Bnorm<1e-4 ] = NA
    lines( p$m, Bnorm )
    points( p$m, Bnorm, pch=i)
  }
  
  legend(x="topright", bty="n",
         legend=textLegend, #c("n=6","n=8","n=10","n=15","n=25","n=50"),
         pch = seq(1,6),
  )
}

# plotFR = function() {
#   p = parameters()
#   B = seq(0,1000)
#   r = matrix(0,length(B),p$n)
#   for (i in 1:length(B)) {
#     rates = calcRates(0,0,0,B[i]*rep(1,p$n),p)
#     r[i,] = rates$Jtot/p$m
#   }
#   
#   defaultplot()
#   defaultpanel(xlim=B, ylim=c(0,1.5))
#   for (i in 1:p$n) {
#     lines(B, r[,i], lwd=i/3)  
#   }
# }

plotSheldonComparison = function(L = 100, n=20) {
  p = parametersChemostat(parameters())
  p$L = L
  d = 10^seq(log10(3.1e-4),0,length=n)
  Delta = log(p$m[2]/p$m[1])
  #
  # Run across a range of d:
  #
  sim = list()
  fit = data.frame()
  for (i in 1:length(d)) {
    p$d = d[i]
    p$tEnd = 1000
    
    sim[[i]] = simulateChemostatEuler(p)#, useF=TRUE)
    
    fit = rbind(fit, as.data.frame(calcSheldonFit(sim[[i]], FALSE)))
  }
  fit$d = d
  
  defaultplot(mfcol=c(3,2))
  #
  # Three examples:
  #
  iSamples = c(3,9,17)
  for (i in iSamples) {
    #loglogpanel(xlim=p$m, ylim=c(1e-3,10))
    calcSheldonFit(sim[[i]], TRUE)
    makepanellabel()
    # ylab=""
    # if (i==2)
    #   ylab="Sheldon spectrum (${$gC/l)"
    # xlab=""
    # if (i==3)
    #   xlab="Carbon mass (${\\mu}$gC)"
    # loglogpanel(xlim=p$m, ylim=c(0.1,100),xlab=xlab,ylab=ylab)
    # lines(p$m, sim[[iSamples[i]]]$B/fit$Delta[iSamples[i]], lwd=3)
    # lines(p$m, )
  }
  #
  # Kappa
  #
  vertlines = function() {
    for (i in iSamples)
      vline(d[i])
  }
  
  loglogpanel(xlim=d, ylim=c(0.01,1e3), 
              ylab="Height ($\\kappa$)")
  lines(d,fit$kappa,lwd=3)
  lines(d,fit$kappaTheo,lty=dashed)
  vertlines()
  makepanellabel()
  #
  # Lambda
  #
  semilogxpanel(xlim=d, ylim=c(-1,0.5),
                ylab="Slope ($\\lambda$)")
  lines(d, -fit$lambda, lwd=3)
  lines(d, d*0, lty=dashed)
  makepanellabel()
  vertlines()
  #
  # Size limits
  #
  rMin = (p$cLeakage+3*p$alphaJ*p$delta) / (p$alphaJ-p$cR*p$alphaJ)
  rMax = p$epsilonF*p$cF / ( p$cR*p$alphaJ  )
  
  ix = fit$kappa < 0.01
  fit$mMin[ix] = fit$mMax[ix]
  
  loglogpanel(xlim=d, ylim=c(1e-9, 50),
              xlab = "Mixing rate (day$^{-1}$)", ylab="Size limits (mug)")
  ribbon(x=fit$d, ymin=fit$mMin, ymax=fit$mMax)
  vertlines()
  hline( calcMass(rMin), lty=dashed)
  hline( calcMass(rMax), lty=dashed)
  makepanellabel()
}

plotDOC = function() {
  
  panelsExudation_vs_size = function(s = baserunChemostat(), sTitle, yaxis=TRUE) {
    m = s$p$m
    r = s$rates
    B = s$B
    p = s$p
    #
    # Calculate functions:
    #
    func = calcFunctionsChemostat(s$p, s$r, s$p$L, s$N, s$B)
    PP = func$prodCgross # Use gross PP
    
    #defaultplotvertical(2)
    
    
    #
    # Rates:
    #
    # ylab = ""
    # if (yaxis) 
    #   ylab = "Rate (1/day)"
    # semilogxpanel(xlim=m, xaxis = FALSE, 
    #               ylim=c(0,0.4), ylab=ylab, yaxis=yaxis)
    # 
    # lines(m, r$jDOCprodPhoto, col='green', lwd=thick)
    # lines(m, r$jDOCprodPassive, col='darkgreen', lwd=thick)
    # lines(m, r$jDOCprodFeeding, col='red', lwd=thick)
    # lines(m, r$jDOCprodVirulysis, col="orange", lwd=thick)
    # jDOCprodTot = r$jDOCprodPhoto+r$jDOCprodPassive+
    #   r$jDOCprodFeeding+r$jDOCprodVirulysis
    # lines(m, jDOCprodTot, col="black")
    # lines(m, r$jTot, col='black', lwd=thick)
    # 
    # mtext(sTitle, side=top, line=-1, padj=-2)
    
    #legend(x="topright", legend=c("Photouptake","Passive","Feeding","Lysis","Total"),
    #       lwd=thick,bty="n",
    #       col=c("green","darkgreen","darkred","orange","black"))
    
    #
    # Fractions of growth
    #
    ylab = ""
    if (yaxis) 
      ylab = "Fraction of growth"
    semilogxpanel(xlim=m, xaxis = FALSE, 
                  ylim=c(0,1), ylab=ylab, yaxis=yaxis)
    
    g = r$jTot
    jDOCprodTot = r$jDOCprodPhoto+r$jDOCprodPassive+
      r$jDOCprodFeeding+r$jDOCprodVirulysis
    # Find only sizes with positive growth:
    ix = (r$jTot-r$mortpred-r$mortHTL-r$mort2) > 0
    
    lines(m[ix], r$jDOCprodPhoto[ix]/g[ix], col='green', lwd=thick)
    lines(m[ix], r$jDOCprodPassive[ix]/g[ix], col='darkgreen', lwd=thick)
    lines(m[ix], r$jDOCprodFeeding[ix]/g[ix], col='red', lwd=thick)
    lines(m[ix], r$jDOCprodVirulysis[ix]/g[ix], col="orange", lwd=thick)
    lines(m[ix], jDOCprodTot[ix]/g[ix],
          col="black", lwd=thick)
    mtext(sTitle, side=top, line=-1, padj=-2)
    makepanellabel()
    #
    # Total losses
    #
    if (yaxis)
      ylab = "Total ($\\mu$gC/l/day)"
    loglogpanel(xlim=m, xlab="Carbon mass ($\\mu$gC)",
                ylim=c(0.001,10000000), 
                ylab=ylab, yaxis=yaxis)
    lines(m, r$jDOCprodPhoto*B, col="green", lwd=thick)
    lines(m, r$jDOCprodPassive*B, col="darkgreen", lwd=thick)
    lines(m, r$jDOCprodFeeding*B, col="red", lwd=thick)
    lines(m, r$jDOCprodVirulysis*B, col="orange", lwd=thick)
    lines(m, jDOCprodTot*B, col="black")
    
    totloss = func$lossPhotouptake + func$lossPassive + p$reminF*func$lossFeeding  + sum(r$jDOCprodVirulysis*B)
    legend(x="topright", bty="n", y.intersp = 0.7,
           legend=c("Fraction of PP:",
                    sprintf("Photo losses %1.2f", func$lossPhotouptake/PP),
                    sprintf("Passive %1.2f", func$lossPassive/PP),
                    sprintf("Feeding %1.2f", p$reminF*func$lossFeeding/PP),
                    sprintf("Virulysis %1.2f", sum(r$jDOCprodVirulysis*B)/PP),
                    sprintf("Total %1.2f", totloss/PP)),
           #sprintf("HTL feeding %1.2f", func$lossFeedingHTL)),
           lwd=c(solid,thick,thick,thick,thick,solid),
           col=c("white","green","darkgreen","red",'orange','black'))
    makepanellabel()
  }
  
  #    cLeakage = 0.00015
  
  pOlig = parametersChemostat()
  #    pOlig$cLeakage = cLeakage
  pOlig$L = LOligotrophic
  pOlig$d = dOligotrophic
  #pOlig$mortHTL = mHTLOligotrophic
  sOlig = simulateChemostatEuler(pOlig)
  
  pEut = parametersChemostat()
  #    pEut$cLeakage = cLeakage
  pEut$L = LEutrophic
  pEut$d = dEutrophic
  #pEut$mortHTL = mHTLEutrophic
  sEut = simulateChemostatEuler(pEut)
  
  defaultplot(c(2,2))
  panelsExudation_vs_size(sOlig, sTitle="Oligotrophic")
  panelsExudation_vs_size(sEut, sTitle="Eutrophic", yaxis=FALSE)
}

#
# Plot bacterial generation time vs the total surface area
#
# Generation time is calculated as the average effective growth rate 
# due to DOC uptake (minus exudation and respiration). The mean is taken over
# all plankton (leads to some overestimation)
#
# Surface area is calculated as the surface area of cell in the size range
# 2.5 to 60 mum and multiplied by their abundance.
#
# Corresponds to figure 7 in Kiørboe et al (1990)
#
plotBacteriaGenerationTime_vs_area = function(p=parametersChemostat(),
                                              n=10,
                                              L=seq(10,60,length.out=n), #seq(0.2*LOlig, 2*LEut,length.out=n),
                                              d = seq(0.001,0.1,length.out = n),
                                              bPlot=TRUE) {
  
  area = matrix(0,nrow=length(L),ncol=length(d))
  gentime = area
  for (j in 1:length(L)) {   
    for (i in 1:length(d)) {
      p$L = L[j]
      p$d = d[i]
      s = simulateChemostatEuler(p)
      func = calcFunctionsChemostat(s$p, s$r, p$L, s$N, s$B)
      
      ESD = calcESD(p$m)
      ix = ESD>5 & ESD<60
      
      area[j,i] = sum(pi*ESD[ix]^2 * s$B[ix]/p$m[ix]) # mum^2/liter
      
      jBact = func$jBact_m
      gentime[j,i] = 1/mean(jBact[jBact>0])
    }
  }
  #
  # Fit a power-law
  #
  fit = lm(gen ~ area, 
           data.frame(gen=as.vector(log(gentime[,])), 
                      area=as.vector(log(area[,]))))
  cat('Fit. Exponent: ', fit$coefficients[2],
      '. Factor: ', exp(fit$coefficients[1]))
  #
  # Fit with exponent -0.82:
  #
  y = gentime[,]/area[,]^-0.82
  factor = exp( mean( log( y[(!is.nan(y)) & (gentime[,]<40)]   ) ) )
  
  defaultplot()
  defaultpanel(xlim=area[,], ylim=c(0,50),
               xlab = "Particle surface area ($\\mu$m$^2$/l)",
               ylab = "Bacterial generation time (d)")
  points(area[,], gentime[,], pch=dots)
  
  #new = data.frame(area=log(seq(0,1e9,length.out=100)))
  #lines(exp(new$area), exp(predict(fit,new)))
  
  x = seq(5e7,4e8,length.out=100)
  lines(x, factor*x^-0.82)
}

plotHTL = function(d=dEutrophic, L=LEutrophic) {
  p = parametersChemostat(parameters(n=50))
  p$d = d
  p$L = L
  mHTL = seq(0,0.5, length.out=25)
  
  defaultplot(mfcol=c(2,1))
  loglogpanel(xlim=p$m, ylim=c(1,100), 
              xlab = 'Mass ($\\mu g_C$)',
              ylab = 'Sheldon biomass ($\\mu$gC/l)')
  B = NA
  NPP = NA
  prodHTL = NA
  ix = floor(seq(1, length(mHTL), length.out=5))
  for (i in 1:length(mHTL)) {
    p$mortHTL = mHTL[i]
    sim = simulateChemostat(p)
    if (sum(i == ix) != 0)
      lines(p$m, sim$B / log(p$Delta), lwd=1+i/length(mHTL))
    
    func = calcFunctionsChemostat(p,sim$rates,p$L,sim$N,sim$B)
    B[i] = sum(sim$B)
    NPP[i] = func$prodCnet
    prodHTL[i] = func$prodHTL
  }
  ribbon(x=p$m, ymin=0*p$m+0.01, ymax=10^(p$mortHTLm-2))
  
  defaultpanel(xlim=mHTL, ylim=c(0,100*floor(1.8*max(B)/100)+1),
               xlab="Higher trophic level mortality $\\mu_{htl}$ (day$^{-1}$)")
  
  lines(mHTL, B, lwd=2)
  lines(mHTL, NPP, col='green', lwd=2)
  lines(mHTL, prodHTL, col='red',lwd=2)
  
  legend(x="topright", bty="n", cex=cex,
         legend=c(TeX('Biomass (${\\mu}g_C$/l)'), 
                  TeX('NPP (g_C/yr/m^2)'), 
                  TeX('prod. HTL (g_C/yr/m^2)')),
         col=c('black','green','red'), lwd=rep(2,3))
  
}

plotTemperature = function(n=50) {
  
  panelTemperature = function(d=dEutrophic, L=LEutrophic, bLegend=TRUE) {
    p = parametersChemostat(parameters(n))
    p$d = d
    p$L = L
    T = seq(0,25,length.out=30)
    #
    # Simulate temperature increase:
    #
    B = NA
    Bsheldon = matrix(0, length(T), n)
    NPP = NA
    prodHTL = NA
    jTot = matrix(NA,length(T),p$n)
    ix = floor(seq(1, length(T), length.out=5))
    for (i in 1:length(T)) {
      p$T = T[i]
      sim = simulateChemostat(p)
      #    if (sum(i == ix) != 0)
      Bsheldon[i,] = sim$B / log(p$Delta)
      #          lines(p$m, sim$B / p$m[2]*p$m[1], lwd=1+i/length(T))
      
      func = calcFunctionsChemostat(p,sim$rates,p$L,sim$N,sim$B)
      B[i] = sum(sim$B)
      NPP[i] = func$prodCnet
      prodHTL[i] = func$prodHTL
      r = getFrates(p)
      jTot[i,] = r$jTot
    }
    
    
    #
    # Q10's:
    #
    ylab = ""
    if (bLegend)
      ylab = "Division rate  $\\textit{j}$_{net} (yr^{-1})"
    semilogypanel(xlim=T,ylim=c(0.05,2),
                  ylab=ylab, xaxis=FALSE, yaxis=bLegend)
    # Reference lines:
    lines(T, 0.3*fTemp(1.25,T,0),col="red")
    lines(T, 0.3*fTemp(1.5,T,0),col="red")
    lines(T, 0.3*fTemp(2,T,0),col="red")
    # Actual growth responses:
    ix = seq(7,45,by=10)
    for (i in 1:length(ix))
      lines(T, jTot[,ix[i]], lwd=i/2)
    makepanellabel()
    #
    # Functions:
    #
    semilogypanel(xlim=T, ylim=c(1,2500),
                  xlab="Temperature", yaxis=bLegend)
    
    lines(T, B, lwd=2)
    lines(T, NPP, col='green', lwd=2)
    lines(T, prodHTL, col='red',lwd=2)
    
    makepanellabel()
    
    if (bLegend) {
      legend(x="bottomright", bty="n", cex=cex,
             legend=c(TeX('Biomass (${\\mu}g_C$/l)'), 
                      TeX('NPP (g_C/yr/m^2)'), 
                      TeX('prod HTL (g_C/yr/m^2)')),
             col=c('black','green','red'), lwd=rep(2,3))
    }
    #
    # Spectra
    #
    if (bLegend)
      ylab = 'Sheldon biomass ($\\mu$gC/l)'
    loglogpanel(xlim=p$m, ylim=c(.1,100), 
                xlab = 'Cell mass ($\\mu g_C$)',
                ylab = ylab, yaxis=bLegend)
    for (i in seq(1,length(T),by=5))
      lines(p$m, Bsheldon[i,], lwd=1+i/length(T))
    makepanellabel()
  }
  
  defaultplot(mfcol=c(3,2))
  panelTemperature(dEutrophic)
  panelTemperature(dOligotrophic, bLegend=FALSE)
  
}


# plotLeakage = function() {
#   p = parameters()
#   
#   defaultplot()
#   semilogxpanel(xlim = p$r, ylim = c(0,0.25), 
#                 xlab = "Radius ({\\mu}m)",
#                 ylab = "Loss rate (day^{-1})")
#   lines(p$r, p$cLeakage/p$r, col="red", lwd=3)
# }





