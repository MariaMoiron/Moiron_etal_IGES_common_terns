# Code for "Understanding the Social Dynamics of Breeding Phenology: Indirect Genetic Effects and Assortative Mating in a Long-Distance Migrant"
# [doi: 10.1086/711045]
# Moiron M, Araya-Ajoy YG, Teplitsky C, Bouwhuis S, Charmantier A   

# The code provided here is sufficient to replicate the simulations presented in the above paper. 
# If you are interested in running the same or additional simulations, and need help
# feel free to get in touch with the following researchers:
#  - Dr Maria Moiron (CEFE-CNRS, Montpellier): mariamoironc<at>gmail.com
#  - Dr Yimen Araya-Ajoy (Norwegian University of Science and Technology, Trondheim): yimencr<at>gmail.com 


# We need to highlight that there is, at least, one limitation of our simulation approach: that two divorced individuals can not re-mate with each other.
# i.e divorced individuals con only mate with new individuals
# Also, please note, that assortative mating is caused by females "choosing a male", but also males "choosing" a female.

# Finally, we also want to highlight that there are several ways of running these analyses. This is just one of those, and it is not necessarily the fastest or most elegant.


######################################################
# DATA SIMULATIONS
######################################################

# Load R-packages
require(arm)
library(tidyr)
library(dplyr)
library (broom)
library(MCMCglmm)
library(beepr) 

# Simulate phenotypic data 
Vaf<-1 #female genetic variance
Vam<-1 #male genetic variance
Vrf<-1 #female residual variance
Vrm<-1 #residual male variance

n.ind<-729 #number of individuals per year; actual sample size = 729 IDs
n.years<-28 #number of sampling years; actual sample size = 28 years
n.sim<-1 #number of replicates per each simulated scenario

dat<-list() #Generate list to store data
simres<-list(data=list(), var=list(), cov=list(), asma=list(),psi=list(), batch=list()) #Generate list to store simulated data sets

asma  <-0 #among-pairs correlation, i.e., magnitude of assortative mating; observed estimate of assortative mating=0.77
psi <- 0.2 #regresion value of male trait on female trait, we simulated there are no male effects

par <- expand.grid(asma, psi) #Parameters that can be changed
par2<-par[rep(seq_len(nrow(par)), n.sim), ] #Parameters that can be changed
colnames(par2) <- c("asma", "psi")
batch <- nrow(par2)

#Simulation loop
for(i in 1:batch){
  asma<-par2$asma[i]
  psi<-par2$psi[i]
  
  asmaB1 <- (asma*sqrt(Vaf*Vam))/Vaf #we converted the correlation to a slope to reflect that assortative mating is caused by mate choice. 
  asmaB2 <- (asma*sqrt(Vaf*Vam))/Vam #we converted the correlation to a slope to reflect that assortative mating is caused by mate choice. 
  
  ##Initial year
  af<-rnorm(n.ind, 0, sqrt(Vaf)) #individual BLUP of female (equivalent of breeding value but at individual level)
  am<-af*asmaB1 + rnorm(n.ind, 0, sqrt(Vam-asmaB1^2)) #individual BLUP of male (equivalent of breeding value but at individual level)
  
  zm<- am + rnorm(n.ind, 0, sqrt(Vrf)) #male phenotype, we add observations-specific enviromental variance
  zf<- af + am*psi + rnorm(n.ind, 0, sqrt(Vrm)) #female phenotype,we add observations-specific enviromental variance
  
  df<-data.frame(maleId=1:n.ind, femaleId=1:n.ind, zm, am, zf, af, year=1) #data of initial population (at year 1)
  
  # two divorced individuals can not remate, but they can mate with a new individual 
  for(a in 2:n.years){
    
    ##################################################
  ##Next years
     s<-sample(1:n.ind, n.ind/2) #half of the population survive and re-mate with the same partner in the next year
    dftmp1<-df[s,] #select surviving pairs that mate with the same partner
    dftmp1$zm<-dftmp1$am + rnorm(nrow(dftmp1), 0, sqrt(Vrm)) #male phenotype include unexplained environmental variation for that year
    dftmp1$zf<- dftmp1$af + dftmp1$am*psi + rnorm(nrow(dftmp1), 0, sqrt(Vrf)) #female phenotype include unexplained environmental variation for that year
    
    dftmp2<-df[-s,] #rest of the population that either do not survive or that survive but do not re-mate with the same partner in the next year
    
    Maletmp<-sample(dftmp2$maleId, nrow(dftmp2)/2) #half of the surviving males mate with a different partner in the next year
    amtmp1<-df$am[match(Maletmp, df$maleId)] #individual BLUPS of divorced males in the next year
    zmtmp1<-amtmp1 + rnorm(length(amtmp1), 0, sqrt(Vrm)) #phenotypic values of divorced males in the next year
    aftmp2<-amtmp1*asmaB2 + rnorm(length(amtmp1), 0, sqrt(Vaf-asmaB2^2)) #individual BLUPS of their new female mates, based on mate choice generating assortative mating
    zftmp2<- aftmp2 + amtmp1*psi + rnorm(length(amtmp1), 0, sqrt(Vrf)) #phenotypic values of their new female mates
      
    Femaletmp<-sample(dftmp2$femaleId, nrow(dftmp2)/2) #half of the surviving females mate with a different partner in the next year
    aftmp1<-df$af[match(Femaletmp, df$femaleId)] #individual BLUPS of divorced females in the next year
    amtmp2<-aftmp1*asmaB1 + rnorm(length(amtmp1), 0, sqrt(Vam-asmaB1^2)) #individual BLUPS of their new male mates, based on mate choice rather than assortative mating
    zmtmp2<- amtmp2 + rnorm(length(amtmp1), 0, sqrt(Vrf)) #phenotypic values of their new male mates
    zftmp1<-aftmp1 + amtmp2*psi + rnorm(length(amtmp1), 0, sqrt(Vrf)) #phenotypic values of divorced females in the next year
    
    df2<-data.frame(maleId=c(Maletmp, paste(Maletmp, a)), am=c(amtmp1, amtmp2), zm=c(zmtmp1, zmtmp2), femaleId=c(paste(Femaletmp, a), Femaletmp), af=c(aftmp2, aftmp1), zf=c(zftmp2, zftmp1), year=a)
    
    df3<-rbind(dftmp1, df2)
    df3$year<-a
    dat[[a]]<-df3
  }
  
  data<-do.call(rbind.data.frame, dat)
  
  #We fit the simulated data to the statistical model
  
  #Set prior
  priorA <- list(G = list(G1 = list(V = 1, nu = 0.002, alpha.mu=0, alpha.V=diag(1)*1000),
                          G2 = list(V = 1, nu = 0.002, alpha.mu=0, alpha.V=diag(1)*1000),
                          G3 = list(V = 1, nu = 0.002, alpha.mu=0, alpha.V=diag(1)*1000)),
                 R=list(V=1, nu=0.002))
 
  #Set model and retrive results over loops
  mod<-MCMCglmm(zf~ 1,
                random=~maleId+femaleId+year,
                data=data,
                nitt=13000, thin=10, burnin=3000,
                prior=priorA,
                verbose=FALSE, family="gaussian")
  

  simres$data[[i]]<-data
  
  results<-data.frame(varM=round(posterior.mode(mod$VCV[,1]),3), varF=round(posterior.mode(mod$VCV[,2]),3), varY= round(posterior.mode(mod$VCV[,3]),3), 
                      covar= cov(data$zm, data$zf), asma=asma, psi=psi, batch=i)
  
  simres$results[[i]]<-results
  
}
