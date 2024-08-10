
library(tidyverse) # plotting and manipulation
library(grid) # combining plots
library(gridExtra) # combining plots
library(ggpubr) # combining plots
library(patchwork) # combining plots
library(ggfortify) # nice extension for ggplot
library(mgcv) #fitting gam models
library(GGally) # displaying pairs panel
library(caTools) # split dataset
library(readxl)
library(bestNormalize)
library(lme4)
library(AICcmodavg)
library(rsq)
library(lmerTest)
library(optimx)
library(glmmTMB)
library(MASS)
library(cluster)
library(ggstatsplot)						
library(FactoMineR)#PCANguyen
library(factoextra) #PCANguyen
library(viridis)


library(readxl)
Data_Penelitian <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "R")
View(Data_Penelitian)
str(Data_Penelitian)
		tibble [27 x 22] (S3: tbl_df/tbl/data.frame)
		 $ Zona_Pengelolaan: chr [1:27] "Piringan" "Piringan" "Piringan" "Piringan" ...
		 $ Kedalaman       : chr [1:27] "0-30" "30-60" "60-90" "0-30" ...
		 $ Kadar_air       : num [1:27] 352 723 1084 417 558 ...
		 $ pH              : num [1:27] 3.36 3.29 3.3 3.35 3.57 3.26 3.71 3.28 3.33 3.42 ...
		 $ C_organik       : num [1:27] 47.8 54.8 53.2 52.2 46.5 ...
		 $ Kadar_abu       : num [1:27] 17.65 5.46 8.23 9.92 19.79 ...
		 $ CN_rasio        : num [1:27] 57.1 87.9 142.2 123.7 115.7 ...
		 $ N_total         : num [1:27] 0.836 0.624 0.374 0.422 0.402 ...
		 $ P_total         : num [1:27] 56.5 52.9 46.3 42.1 59.2 ...
		 $ K_total         : num [1:27] 55.3 63.3 61.3 54.3 58.3 ...
		 $ Fe_total        : num [1:27] 1661 477 492 870 1366 ...
		 $ Mn_total        : num [1:27] 33.86 11.75 6.1 6.87 15.41 ...
		 $ Cu_total        : num [1:27] 5.78 7.18 7.91 3.37 2.97 ...
		 $ Zn_total        : num [1:27] 13.46 13.14 8.85 8.94 8.85 ...
		 $ PBL             : num [1:27] 6.60e+06 9.65e+06 1.48e+07 1.50e+08 3.44e+08 ...
		 $ PBS             : num [1:27] 8700000 12150000 10350000 6700000 23400000 ...
		 $ PFS             : num [1:27] 6818 4950 4265 5970 3167 ...
		 $ Lks             : num [1:27] 208 521 255 197 4873 ...
		 $ Mn_P            : num [1:27] 20.7 27.5 34.4 41.3 20.7 ...
		 $ Li_P            : num [1:27] 1245.5 17.9 80.6 152.3 116.5 ...
		 $ EBS             : num [1:27] 50.741 12.222 3.056 0.648 1.111 ...
		 $ EFS             : num [1:27] 9.63 4.54 5.93 9.44 2.96 ...

# omit missing data 
Data_Penelitian <- na.omit(Data_Penelitian)

## LME ReML berbasis data raw
Data_Penelitian <- as.data.frame(Data_Penelitian)
#selecting data

datafact<- dplyr::select(Data_Penelitian, c(1:2)) #quantitative selecting column
datanum <- dplyr::select(Data_Penelitian, c(3:22)) #quantitative selecting column

#assign IDs
datafact <- tibble::rowid_to_column(datafact, "ID")



# box cox transform
baru <- sapply(datanum,\(x)boxcox(x)$x.t)
baru

df <- data.frame(baru)
df <- tibble::rowid_to_column(df, "ID")

datafix_hadi <- merge(datafact,df,by="ID")

datafix_hadi$Zona_Pengelolaan <- as.factor(datafix_hadi$Zona_Pengelolaan)
datafix_hadi$Kedalaman        <- as.factor(datafix_hadi$Kedalaman)
## Populasi bakteri lignolitik
PBL.lm<-  lm(PBL ~ Zona_Pengelolaan * Kedalaman, data = datafix_hadi)

PBL.models <- list( )
PBL.models[[1]]  <- lmer(PBL ~ (1 | pH), data = datafix_hadi, REML=FALSE)
PBL.models[[2]]  <- lmer(PBL ~ (1 | K_total), data = datafix_hadi, REML=FALSE)
PBL.models[[3]]  <- lmer(PBL ~ (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PBL.models[[4]]  <- lmer(PBL ~ (1 | Zn_total), data = datafix_hadi, REML=FALSE)
PBL.models[[5]]  <- lmer(PBL ~ (1 | pH) + (1 | K_total), data = datafix_hadi, REML=FALSE)
PBL.models[[6]]  <- lmer(PBL ~ (1 | pH) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PBL.models[[7]]  <- lmer(PBL ~ (1 | pH) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)
PBL.models[[8]]  <- lmer(PBL ~ (1 | K_total) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PBL.models[[9]]  <- lmer(PBL ~ (1 | K_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
PBL.models[[10]]  <- lmer(PBL ~ (1 | Mn_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
PBL.models[[11]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman +  (1 | pH) , 
								data = datafix_hadi, REML=FALSE)									
PBL.models[[12]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)
PBL.models[[13]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman +  (1 | Mn_total) , 
								data = datafix_hadi, REML=FALSE)	
PBL.models[[14]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman +  (1 | Zn_total) , 
								data = datafix_hadi, REML=FALSE)	
PBL.models[[15]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)								
PBL.models[[16]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
PBL.models[[17]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)
PBL.models[[18]]  <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
PBL.models[[19]]  <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
PBL.models[[20]]  <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | Mn_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)									

								
								
##create a vector of names to trace back models in set
PBLnames <- paste("mod", 1:length(PBL.models), sep = " ")

##generate AICc table
aictab(cand.set = PBL.models, modnames = PBLnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = PBL.models, modnames = PBLnames, sort = TRUE),
      digits = 4, LL = TRUE)

Model selection based on AICc:
Model selection based on AICc:

				K     AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 12 11   5.6312     0.0000      1      1  17.6129
		mod 18 12  29.4768    23.8455      0      1   9.2616
		mod 19 12  48.4768    42.8455      0      1  -0.2384
		mod 2   3  77.0533    71.4220      0      1 -34.9812
		mod 1   3  79.3040    73.6728      0      1 -36.1066
		mod 8   4  79.6594    74.0282      0      1 -34.8773
		mod 9   4  79.7258    74.0946      0      1 -34.9105
		mod 3   3  79.7642    74.1330      0      1 -36.3367
		mod 4   3  79.8433    74.2120      0      1 -36.3762
		mod 5   4  79.8671    74.2359      0      1 -34.9812
		mod 6   4  80.7894    75.1582      0      1 -35.4423
		mod 7   4  81.4319    75.8006      0      1 -35.7636
		mod 10  4  82.5781    76.9468      0      1 -36.3367
		mod 14 11  84.0611    78.4298      0      1 -21.6020
		mod 11 11  97.7581    92.1268      0      1 -28.4505
		mod 13 11 101.7651    96.1338      0      1 -30.4540
		mod 16 12 103.7242    98.0930      0      1 -27.8621
		mod 15 12 103.7716    98.1403      0      1 -27.8858
		mod 17 12 104.9009    99.2697      0      1 -28.4505
		mod 20 12 108.9080   103.2767      0      1 -30.4540

AICc(PBL.lm)
[1] 95.70571

PBL.models.fix <- list( )									
PBL.models.fix[[1]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi , REML=FALSE)	
PBL.models.fix[[2]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
PBL.models.fix[[3]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
PBL.models.fix[[4]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
PBL.models.fix[[5]]   <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	
##create a vector of names to trace back models in set
PBLnames.fix <- paste("mod", 1:length(PBL.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = PBL.models.fix, modnames = PBLnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = PBL.models.fix, modnames = PBLnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)

		Model selection based on AICc:

			   K    AICc Delta_AICc AICcWt Cum.Wt      LL
		mod 1 11  5.6312     0.0000 0.9991 0.9991 17.6129 ## Final Model
		mod 3 11 19.5858    13.9545 0.0009 1.0000 10.6357
		mod 2 11 28.4653    22.8340 0.0000 1.0000  6.1959
		mod 4 11 33.7434    28.1122 0.0000 1.0000  3.5568	  
	  
PBLmodelfit.all <- lme4::allFit(PBL.models.fix[[1]]  )
		bobyqa : [OK]
		Nelder_Mead : [failed] ## <<- fail
		nlminbwrap : [OK]
		optimx.L-BFGS-B : [OK]
		nloptwrap.NLOPT_LN_NELDERMEAD : [OK]
		nloptwrap.NLOPT_LN_BOBYQA : [OK]

PBL_ss <- summary(PBLmodelfit.all)
 
## Final model 
PBL.model.fix  <- lmer(PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi )
summary(PBL.model.fix)
		Linear mixed model fit by REML. t-tests use Satterthwaite method ['lmerModLmerTest']
		Formula: PBL ~ Zona_Pengelolaan * Kedalaman + (1 | K_total)
		   Data: datafix_hadi

		REML criterion at convergence: 55

		Scaled residuals: 
			   Min         1Q     Median         3Q        Max 
		-2.116e-03 -7.671e-04 -2.463e-05  8.035e-04  1.835e-03 

		Random effects:
		 Groups   Name        Variance  Std.Dev.
		 K_total  (Intercept) 9.325e-01 0.965667
		 Residual             1.429e-06 0.001195
		Number of obs: 26, groups:  K_total, 23

		Fixed effects:
													  Estimate Std. Error        df  t value Pr(>|t|)
		(Intercept)                                   0.802736   0.431861 17.000000    1.859  0.08046
		Zona_PengelolaanGawangan Mati                -0.855547   0.610743 17.000000   -1.401  0.17925
		Zona_PengelolaanPiringan                     -1.045725   0.001691 17.000000 -618.523  < 2e-16
		Kedalaman30-60                               -1.370220   0.705225 17.000000   -1.943  0.06877
		Kedalaman60-90                               -1.111439   0.610743 17.000000   -1.820  0.08644
		Zona_PengelolaanGawangan Mati:Kedalaman30-60  2.087243   0.997338 17.000000    2.093  0.05167
		Zona_PengelolaanPiringan:Kedalaman30-60       1.880606   0.737543 17.000000    2.550  0.02072
		Zona_PengelolaanGawangan Mati:Kedalaman60-90  0.081608   1.115058 17.000000    0.073  0.94251
		Zona_PengelolaanPiringan:Kedalaman60-90       1.907193   0.647793 17.000000    2.944  0.00907
														
		(Intercept)                                  .  
		Zona_PengelolaanGawangan Mati                   
		Zona_PengelolaanPiringan                     ***
		Kedalaman30-60                               .  
		Kedalaman60-90                               .  
		Zona_PengelolaanGawangan Mati:Kedalaman30-60 .  
		Zona_PengelolaanPiringan:Kedalaman30-60      *  
		Zona_PengelolaanGawangan Mati:Kedalaman60-90    
		Zona_PengelolaanPiringan:Kedalaman60-90      ** 
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

		Correlation of Fixed Effects:
					(Intr) Zn_PGM Zn_PnP K30-60 K60-90 Z_PGM:K3 Z_PP:K3 Z_PGM:K6
		Zn_PngllnGM -0.707                                                      
		Zn_PngllnPr -0.002  0.001                                               
		Kedlmn30-60 -0.612  0.433  0.001                                        
		Kedlmn60-90 -0.707  1.000  0.001  0.433                                 
		Z_PGM:K30-6  0.433 -0.612 -0.001 -0.707 -0.612                          
		Z_PP:K30-60  0.000  0.000 -0.002 -0.598  0.000  0.423                   
		Z_PGM:K60-9  0.387 -0.822 -0.001 -0.237 -0.822  0.503    0.000          
		Z_PP:K60-90  0.000 -0.471 -0.003  0.000 -0.471  0.289    0.488   0.516  
		optimizer (nloptwrap) convergence code: 0 (OK)
		unable to evaluate scaled gradient
		Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
			
anova(PBL.model.fix, test = "F")
		Type III Analysis of Variance Table with Satterthwaite method
									   Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
		Zona_Pengelolaan           4.9770e-07 2.4890e-07     2    17  0.2863 0.75459  
		Kedalaman                  1.8623e-06 9.3120e-07     2    17  1.0713 0.36460  
		Zona_Pengelolaan:Kedalaman 1.4769e-05 3.6923e-06     4    17  4.2478 0.01453 *
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

rsq(PBL.model.fix, adj='TRUE')
		$model
		[1] 0.9999978

		$fixed
		[1] -0.00856059

		$random
		[1] 1.008558

AICc(PBL.model.fix)	
95.81427	
	
	
## populasi bakteri selulolitik
PBS.lm <- lm(PBS ~ Zona_Pengelolaan * Kedalaman, data = datafix_hadi)	

PBS.models <- list( )
PBS.models[[1]]  <- lmer(PBS ~ (1 | pH), data = datafix_hadi, REML=FALSE)
PBS.models[[2]]  <- lmer(PBS ~ (1 | K_total), data = datafix_hadi, REML=FALSE)
PBS.models[[3]]  <- lmer(PBS ~ (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PBS.models[[4]]  <- lmer(PBS ~ (1 | Zn_total), data = datafix_hadi, REML=FALSE)
PBS.models[[5]]  <- lmer(PBS ~ (1 | pH) + (1 | K_total), data = datafix_hadi, REML=FALSE)
PBS.models[[6]]  <- lmer(PBS ~ (1 | pH) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PBS.models[[7]]  <- lmer(PBS ~ (1 | pH) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)
PBS.models[[8]]  <- lmer(PBS ~ (1 | K_total) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PBS.models[[9]]  <- lmer(PBS ~ (1 | K_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
PBS.models[[10]]  <- lmer(PBS ~ (1 | Mn_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
PBS.models[[11]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman +  (1 | pH) , 
								data = datafix_hadi, REML=FALSE)									
PBS.models[[12]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)
PBS.models[[13]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman +  (1 | Mn_total) , 
								data = datafix_hadi, REML=FALSE)	
PBS.models[[14]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman +  (1 | Zn_total) , 
								data = datafix_hadi, REML=FALSE)	
PBS.models[[15]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)								
PBS.models[[16]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
PBS.models[[17]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)
PBS.models[[18]]  <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
PBS.models[[19]]  <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
PBS.models[[20]]  <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman +(1 | Mn_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
								

##create a vector of names to trace back models in set
PBSnames <- paste("mod", 1:length(PBS.models), sep = " ")

##generate AICc table
aictab(cand.set = PBS.models, modnames = PBSnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = PBS.models, modnames = PBSnames, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

				K     AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 12 11   5.0464     0.0000 0.9871 0.9871  17.9054
		mod 18 12  13.7283     8.6819 0.0129 1.0000  17.1359
		mod 15 12  41.4728    36.4263 0.0000 1.0000   3.2636
		mod 14 11  78.5403    73.4939 0.0000 1.0000 -18.8416
		mod 17 12  79.0298    73.9834 0.0000 1.0000 -15.5149
		mod 4   3  79.1057    74.0593 0.0000 1.0000 -36.0074
		mod 2   3  79.4077    74.3612 0.0000 1.0000 -36.1584
		mod 1   3  79.8560    74.8095 0.0000 1.0000 -36.3825
		mod 3   3  79.8560    74.8095 0.0000 1.0000 -36.3825
		mod 9   4  81.7783    76.7319 0.0000 1.0000 -35.9368
		mod 10  4  81.9196    76.8732 0.0000 1.0000 -36.0074
		mod 7   4  81.9196    76.8732 0.0000 1.0000 -36.0074
		mod 5   4  82.2215    77.1751 0.0000 1.0000 -36.1584
		mod 8   4  82.2215    77.1751 0.0000 1.0000 -36.1584
		mod 6   4  82.6698    77.6234 0.0000 1.0000 -36.3825
		mod 11 11 102.0412    96.9948 0.0000 1.0000 -30.5920
		mod 13 11 102.0413    96.9948 0.0000 1.0000 -30.5921
		mod 16 12 109.1841   104.1377 0.0000 1.0000 -30.5920  

AICc(PBS.lm)
[1] 95.85078

PBS.models.fix <- list( )									
PBS.models.fix[[1]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi , REML=FALSE)	
PBS.models.fix[[2]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
PBS.models.fix[[3]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
PBS.models.fix[[4]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
PBS.models.fix[[5]]   <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	
										
##create a vector of names to trace back models in set
PBSnames.fix <- paste("mod", 1:length(PBS.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = PBS.models.fix, modnames = PBSnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = PBS.models.fix, modnames = PBSnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

			   K    AICc Delta_AICc AICcWt Cum.Wt      LL
		mod 1 11  5.0464     0.0000 0.6650 0.6650 17.9054  ## Model Final
		mod 4 11  6.4999     1.4535 0.3215 0.9866 17.1786
		mod 3 11 12.8495     7.8031 0.0134 1.0000 14.0038
		mod 2 11 28.8033    23.7569 0.0000 1.0000  6.0269
		
PBSmodelfit.all <- lme4::allFit(PBS.models.fix[[1]] )
		bobyqa : [OK]
		Nelder_Mead : [failed] ## fail
		nlminbwrap : [OK]
		optimx.L-BFGS-B : [OK]
		nloptwrap.NLOPT_LN_NELDERMEAD : [OK]
		nloptwrap.NLOPT_LN_BOBYQA : [OK]

PBS_ss <- summary(PBSmodelfit.all)


##Final Model
PBS.model.fix  <- lmer(PBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi)	
summary(PBS.model.fix)		
		Linear mixed model fit by REML. t-tests use Satterthwaite method ['lmerModLmerTest']
		Formula: PBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total)
		   Data: datafix_hadi

		REML criterion at convergence: 54.6

		Scaled residuals: 
			   Min         1Q     Median         3Q        Max 
		-2.042e-03 -6.620e-04 -3.039e-05  5.329e-04  2.162e-03 

		Random effects:
		 Groups   Name        Variance  Std.Dev.
		 K_total  (Intercept) 9.157e-01 0.956896
		 Residual             1.275e-06 0.001129
		Number of obs: 26, groups:  K_total, 23

		Fixed effects:
													  Estimate Std. Error        df  t value Pr(>|t|)
		(Intercept)                                   0.588271   0.427938 17.000000    1.375    0.187
		Zona_PengelolaanGawangan Mati                -0.597343   0.605196 17.000000   -0.987    0.337
		Zona_PengelolaanPiringan                     -1.558117   0.001597 17.000000 -975.855   <2e-16
		Kedalaman30-60                               -0.245249   0.698819 17.000000   -0.351    0.730
		Kedalaman60-90                               -0.380545   0.605196 17.000000   -0.629    0.538
		Zona_PengelolaanGawangan Mati:Kedalaman30-60  0.447491   0.988279 17.000000    0.453    0.656
		Zona_PengelolaanPiringan:Kedalaman30-60       0.453029   0.730844 17.000000    0.620    0.544
		Zona_PengelolaanGawangan Mati:Kedalaman60-90  1.305957   1.104929 17.000000    1.182    0.254
		Zona_PengelolaanPiringan:Kedalaman60-90       0.525902   0.641909 17.000000    0.819    0.424
														
		(Intercept)                                     
		Zona_PengelolaanGawangan Mati                   
		Zona_PengelolaanPiringan                     ***
		Kedalaman30-60                                  
		Kedalaman60-90                                  
		Zona_PengelolaanGawangan Mati:Kedalaman30-60    
		Zona_PengelolaanPiringan:Kedalaman30-60         
		Zona_PengelolaanGawangan Mati:Kedalaman60-90    
		Zona_PengelolaanPiringan:Kedalaman60-90         
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

		Correlation of Fixed Effects:
					(Intr) Zn_PGM Zn_PnP K30-60 K60-90 Z_PGM:K3 Z_PP:K3 Z_PGM:K6
		Zn_PngllnGM -0.707                                                      
		Zn_PngllnPr -0.002  0.001                                               
		Kedlmn30-60 -0.612  0.433  0.001                                        
		Kedlmn60-90 -0.707  1.000  0.001  0.433                                 
		Z_PGM:K30-6  0.433 -0.612 -0.001 -0.707 -0.612                          
		Z_PP:K30-60  0.000  0.000 -0.002 -0.598  0.000  0.423                   
		Z_PGM:K60-9  0.387 -0.822 -0.001 -0.237 -0.822  0.503    0.000          
		Z_PP:K60-90  0.000 -0.471 -0.002  0.000 -0.471  0.289    0.488   0.516 
		


anova(PBS.model.fix, test = "F")
		Type III Analysis of Variance Table with Satterthwaite method
									   Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
		Zona_Pengelolaan           1.3120e-05 6.5602e-06     2    17  5.1465 0.01788 *
		Kedalaman                  6.3980e-07 3.1990e-07     2    17  0.2510 0.78089  
		Zona_Pengelolaan:Kedalaman 2.7758e-06 6.9400e-07     4    17  0.5444 0.70543  
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

rsq(PBS.model.fix, adj='TRUE')
		$model
		[1] 0.9999981

		$fixed
		[1] -0.001786371

		$random
		[1] 1.001784


AICc(PBS.model.fix)	
[1] 95.50403


## populasi fungi selulolitik
PFS.lm <- lm(PFS ~ Zona_Pengelolaan * Kedalaman, data = datafix_hadi)	

PFS.models <- list( )
PFS.models[[1]]  <- lmer(PFS ~ (1 | pH), data = datafix_hadi, REML=FALSE)
PFS.models[[2]]  <- lmer(PFS ~ (1 | K_total), data = datafix_hadi, REML=FALSE)
PFS.models[[3]]  <- lmer(PFS ~ (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PFS.models[[4]]  <- lmer(PFS ~ (1 | Zn_total), data = datafix_hadi, REML=FALSE)
PFS.models[[5]]  <- lmer(PFS ~ (1 | pH) + (1 | K_total), data = datafix_hadi, REML=FALSE)
PFS.models[[6]]  <- lmer(PFS ~ (1 | pH) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PFS.models[[7]]  <- lmer(PFS ~ (1 | pH) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)
PFS.models[[8]]  <- lmer(PFS ~ (1 | K_total) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
PFS.models[[9]]  <- lmer(PFS ~ (1 | K_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
PFS.models[[10]]  <- lmer(PFS ~ (1 | Mn_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
PFS.models[[11]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman +  (1 | pH) , 
								data = datafix_hadi, REML=FALSE)									
PFS.models[[12]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)
PFS.models[[13]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman +  (1 | Mn_total) , 
								data = datafix_hadi, REML=FALSE)	
PFS.models[[14]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman +  (1 | Zn_total) , 
								data = datafix_hadi, REML=FALSE)	
PFS.models[[15]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)								
PFS.models[[16]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
PFS.models[[17]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)
PFS.models[[18]]  <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
PFS.models[[19]]  <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
PFS.models[[20]]  <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman +(1 | Mn_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
								
##create a vector of names to trace back models in set
PFSnames <- paste("mod", 1:length(PFS.models), sep = " ")

##generate AICc table
aictab(cand.set = PFS.models, modnames = PFSnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = PFS.models, modnames = PFSnames, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

				K     AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 12 11  14.8223     0.0000      1      1  13.0174 ## payah di akhir
		mod 4   3  78.8081    63.9858      0      1 -35.8586 ## model surrogate
		mod 3   3  79.7183    64.8960      0      1 -36.3137
		mod 1   3  79.8426    65.0202      0      1 -36.3758
		mod 2   3  79.8560    65.0337      0      1 -36.3825
		mod 7   4  81.3834    66.5610      0      1 -35.7393
		mod 9   4  81.6219    66.7996      0      1 -35.8586
		mod 10  4  81.6219    66.7996      0      1 -35.8586
		mod 6   4  81.8360    67.0137      0      1 -35.9656
		mod 8   4  82.5322    67.7098      0      1 -36.3137
		mod 5   4  82.6564    67.8341      0      1 -36.3758
		mod 11 11 107.3562    92.5339      0      1 -33.2495	

AICc(PFS.lm)
[1] 101.1657


PFS.models.fix <- list( )									
PFS.models.fix[[1]]   <- lmer(PFS ~  (1 | Zn_total), 
										data = datafix_hadi , REML=FALSE)	
PFS.models.fix[[2]]   <- lmer(PFS ~  (1 | Zn_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
PFS.models.fix[[3]]   <- lmer(PFS ~ (1 | Zn_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
PFS.models.fix[[4]]   <- lmer(PFS ~  (1 | Zn_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
PFS.models.fix[[5]]   <- lmer(PFS ~  (1 | Zn_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	


##create a vector of names to trace back models in set
PFSnames.fix <- paste("mod", 1:length(PFS.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = PFS.models.fix, modnames = PFSnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = PFS.models.fix, modnames = PFSnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)

		Model selection based on AICc:

			  K    AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 4 3 78.8081     0.0000 0.2086 0.2086 -35.8586
		mod 1 3 78.8081     0.0000 0.2086 0.4173 -35.8586
		mod 3 3 78.8081     0.0000 0.2086 0.6259 -35.8586
		mod 2 3 78.8081     0.0000 0.2086 0.8345 -35.8586
		mod 5 3 79.2719     0.4638 0.1655 1.0000 -36.0905


PFS.models.fix <- list( )									
PFS.models.fix[[1]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi , REML=FALSE)	
PFS.models.fix[[2]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
PFS.models.fix[[3]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
PFS.models.fix[[4]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
PFS.models.fix[[5]]   <- lmer(PFS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	
										
##create a vector of names to trace back models in set
PFSnames.fix <- paste("mod", 1:length(PFS.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = PFS.models.fix, modnames = PFSnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = PFS.models.fix, modnames = PFSnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

			   K    AICc Delta_AICc AICcWt Cum.Wt      LL
		mod 1 11 14.8223     0.0000 0.5564 0.5564 13.0174  ## Final Model
		mod 4 11 15.2776     0.4553 0.4432 0.9996 12.7898
		mod 3 11 29.3744    14.5521 0.0004 1.0000  5.7414
		mod 2 11 35.2152    20.3929 0.0000 1.0000  2.8210

PFSmodelfit.all <- lme4::allFit(PFS.models.fix[[1]] )
		bobyqa : [OK]
		Nelder_Mead : [failed]
		nlminbwrap : [OK]
		optimx.L-BFGS-B : [OK]
		nloptwrap.NLOPT_LN_NELDERMEAD : [failed]
		nloptwrap.NLOPT_LN_BOBYQA : [OK]

PFS_ss <- summary(PFSmodelfit.all)


##Final Model
PFS.model.fix  <- lmer(PFS ~ (1 | Zn_total), 
										control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)),
										data = datafix_hadi)	
summary(PFS.model.fix)		
		Linear mixed model fit by REML. t-tests use Satterthwaite method [
		lmerModLmerTest]
		Formula: PFS ~ (1 | Zn_total)
		   Data: datafix_hadi
		Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))

		REML criterion at convergence: 73.1

		Scaled residuals: 
			 Min       1Q   Median       3Q      Max 
		-0.81856 -0.13465  0.01233  0.15517  0.69825 

		Random effects:
		 Groups   Name        Variance Std.Dev.
		 Zn_total (Intercept) 0.93602  0.9675  
		 Residual             0.09202  0.3033  
		Number of obs: 26, groups:  Zn_total, 25

		Fixed effects:
					Estimate Std. Error       df t value Pr(>|t|)
		(Intercept)  0.01723    0.20260 23.92069   0.085    0.933

## ga dipakai
anova(PFS.model.fix, test = "F")
		Type III Analysis of Variance Table with Satterthwaite method
								   Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
		Zona_Pengelolaan           1.3275 0.66374     2    17  0.5743 0.5736
		Kedalaman                  0.9759 0.48796     2    17  0.4222 0.6623
		Zona_Pengelolaan:Kedalaman 3.5368 0.88419     4    17  0.7651 0.5624

rsq(PFS.model.fix, adj='TRUE')

		$model
		[1] 0.9048526

		$fixed
		[1] -0.0003087346

		$random
		[1] 0.9051613
		
AICc(PFS.model.fix)	
		80.18392



## Aktivitas enzim lakase
LKS.lm <- lm(LKS ~ Zona_Pengelolaan * Kedalaman, data = datafix_hadi)	

LKS.models <- list( )
LKS.models[[1]]  <- lmer(LKS ~ (1 | pH), data = datafix_hadi, REML=FALSE)
LKS.models[[2]]  <- lmer(LKS ~ (1 | K_total), data = datafix_hadi, REML=FALSE)
LKS.models[[3]]  <- lmer(LKS ~ (1 | Mn_total), data = datafix_hadi, REML=FALSE)
LKS.models[[4]]  <- lmer(LKS ~ (1 | Zn_total), data = datafix_hadi, REML=FALSE)
LKS.models[[5]]  <- lmer(LKS ~ (1 | pH) + (1 | K_total), data = datafix_hadi, REML=FALSE)
LKS.models[[6]]  <- lmer(LKS ~ (1 | pH) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
LKS.models[[7]]  <- lmer(LKS ~ (1 | pH) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)
LKS.models[[8]]  <- lmer(LKS ~ (1 | K_total) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
LKS.models[[9]]  <- lmer(LKS ~ (1 | K_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
LKS.models[[10]]  <- lmer(LKS ~ (1 | Mn_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
LKS.models[[11]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman +  (1 | pH) , 
								data = datafix_hadi, REML=FALSE)									
LKS.models[[12]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)
LKS.models[[13]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman +  (1 | Mn_total) , 
								data = datafix_hadi, REML=FALSE)	
LKS.models[[14]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman +  (1 | Zn_total) , 
								data = datafix_hadi, REML=FALSE)	
LKS.models[[15]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)								
LKS.models[[16]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
LKS.models[[17]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)
LKS.models[[18]]  <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
LKS.models[[19]]  <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
LKS.models[[20]]  <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman +(1 | Mn_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
								

##create a vector of names to trace back models in set
LKSnames <- paste("mod", 1:length(LKS.models), sep = " ")

##generate AICc table
aictab(cand.set = LKS.models, modnames = LKSnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = LKS.models, modnames = LKSnames, sort = TRUE),
      digits = 4, LL = TRUE)
			Model selection based on AICc:

				K     AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 13 11  12.1672     0.0000      1      1  14.3450 
		mod 14 11  70.5392    58.3721      0      1 -14.8410
		mod 1   3  75.5141    63.3469      0      1 -34.2116
		mod 7   4  77.2564    65.0892      0      1 -33.6758
		mod 6   4  78.3279    66.1608      0      1 -34.2116
		mod 5   4  78.3279    66.1608      0      1 -34.2116
		mod 3   3  78.9853    66.8182      0      1 -35.9472
		mod 4   3  79.8294    67.6622      0      1 -36.3692
		mod 2   3  79.8560    67.6888      0      1 -36.3825
		mod 10  4  81.7565    69.5894      0      1 -35.9259
		mod 8   4  81.7992    69.6320      0      1 -35.9472
		mod 9   4  82.6432    70.4761      0      1 -36.3692
		mod 11 11 100.8481    88.6809      0      1 -29.9955
		mod 12 11 104.8529    92.6858      0      1 -31.9979
		mod 15 12 107.9909    95.8238      0      1 -29.9955
	

AICc(LKS.lm)
[1] 98.66244

LKS.models.fix <- list( )									
LKS.models.fix[[1]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | Mn_total), 
										data = datafix_hadi , REML=FALSE)	
LKS.models.fix[[2]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | Mn_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
LKS.models.fix[[3]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | Mn_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
LKS.models.fix[[4]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | Mn_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
LKS.models.fix[[5]]   <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | Mn_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	

										
##create a vector of names to trace back models in set
LKSnames.fix <- paste("mod", 1:length(LKS.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = LKS.models.fix, modnames = LKSnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = LKS.models.fix, modnames = LKSnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

			   K    AICc Delta_AICc AICcWt Cum.Wt      LL
		mod 4 11 11.9658     0.0000 0.3219 0.3219 14.4457  ## best model
		mod 1 11 12.1672     0.2014 0.2911 0.6129 14.3450
		mod 2 11 12.7999     0.8342 0.2121 0.8251 14.0286
		mod 3 11 13.1854     1.2196 0.1749 1.0000 13.8359
		
LKSmodelfit.all <- lme4::allFit(LKS.models.fix[[1]] )
		bobyqa : [OK]
		Nelder_Mead : [failed]
		nlminbwrap : [OK]
		optimx.L-BFGS-B : [OK]
		nloptwrap.NLOPT_LN_NELDERMEAD : [failed]
		nloptwrap.NLOPT_LN_BOBYQA : [OK]

LKS_ss <- summary(LKSmodelfit.all)


##Final Model
LKS.model.fix  <- lmer(LKS ~ Zona_Pengelolaan * Kedalaman + (1 | Mn_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)))	
summary(LKS.model.fix)		
		Linear mixed model fit by REML. t-tests use Satterthwaite method ['lmerModLmerTest']
		Formula: LKS ~ Zona_Pengelolaan * Kedalaman + (1 | Mn_total)
		   Data: datafix_hadi
		Control: lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))

		REML criterion at convergence: 57.3

		Scaled residuals: 
			 Min       1Q   Median       3Q      Max 
		-0.77710 -0.27479 -0.01398  0.18309  0.77710 

		Random effects:
		 Groups   Name        Variance Std.Dev.
		 Mn_total (Intercept) 0.8110   0.9006  
		 Residual             0.2155   0.4642  
		Number of obs: 26, groups:  Mn_total, 23

		Fixed effects:
													 Estimate Std. Error      df t value Pr(>|t|)
		(Intercept)                                   -0.1259     0.5849 16.3593  -0.215    0.832
		Zona_PengelolaanGawangan Mati                  0.5554     0.8272 16.3593   0.671    0.511
		Zona_PengelolaanPiringan                      -0.7643     0.7812 15.0959  -0.978    0.343
		Kedalaman30-60                                 0.3557     0.8272 16.3593   0.430    0.673
		Kedalaman60-90                                 0.6768     0.7882 15.4332   0.859    0.404
		Zona_PengelolaanGawangan Mati:Kedalaman30-60  -1.3397     1.1233 16.7999  -1.193    0.250
		Zona_PengelolaanPiringan:Kedalaman30-60        1.3690     1.1378 16.9332   1.203    0.245
		Zona_PengelolaanGawangan Mati:Kedalaman60-90  -1.0423     1.1426 16.9611  -0.912    0.374
		Zona_PengelolaanPiringan:Kedalaman60-90       -0.6667     0.9282  4.1632  -0.718    0.511

		Correlation of Fixed Effects:
					(Intr) Zn_PGM Zn_PnP K30-60 K60-90 Z_PGM:K3 Z_PP:K3 Z_PGM:K6
		Zn_PngllnGM -0.707                                                      
		Zn_PngllnPr -0.749  0.529                                               
		Kedlmn30-60 -0.707  0.500  0.529                                        
		Kedlmn60-90 -0.742  0.525  0.593  0.525                                 
		Z_PGM:K30-6  0.521 -0.736 -0.269 -0.736 -0.329                          
		Z_PP:K30-60  0.514 -0.364 -0.687 -0.727 -0.407  0.452                   
		Z_PGM:K60-9  0.512 -0.724 -0.409 -0.362 -0.690  0.493    0.281          
		Z_PP:K60-90  0.630 -0.446 -0.803 -0.446 -0.702  0.286    0.551   0.484

anova(LKS.model.fix, test = "F")
		Type III Analysis of Variance Table with Satterthwaite method
									Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
		Zona_Pengelolaan           0.33763 0.16882     2 11.4043  0.7835 0.4799
		Kedalaman                  0.16911 0.08456     2 11.9459  0.3925 0.6838
		Zona_Pengelolaan:Kedalaman 1.66177 0.41544     4  9.0656  1.9282 0.1893

rsq(LKS.model.fix, adj='TRUE')
		$model
		[1] 0.6890083

		$fixed
		[1] -0.1493072

		$random
		[1] 0.8383155

AICc(LKS.model.fix)	
		98.12945


# aktivitas enzim Mn-PB

Mn_P.lm <- lm(Mn_P ~ Zona_Pengelolaan * Kedalaman, data = datafix_hadi)	

Mn_P.models <- list( )
Mn_P.models[[1]]  <- lmer(Mn_P ~ (1 | pH), data = datafix_hadi, REML=FALSE)
Mn_P.models[[2]]  <- lmer(Mn_P ~ (1 | K_total), data = datafix_hadi, REML=FALSE)
Mn_P.models[[3]]  <- lmer(Mn_P ~ (1 | Mn_total), data = datafix_hadi, REML=FALSE)
Mn_P.models[[4]]  <- lmer(Mn_P ~ (1 | Zn_total), data = datafix_hadi, REML=FALSE)
Mn_P.models[[5]]  <- lmer(Mn_P ~ (1 | pH) + (1 | K_total), data = datafix_hadi, REML=FALSE)
Mn_P.models[[6]]  <- lmer(Mn_P ~ (1 | pH) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
Mn_P.models[[7]]  <- lmer(Mn_P ~ (1 | pH) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)
Mn_P.models[[8]]  <- lmer(Mn_P ~ (1 | K_total) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
Mn_P.models[[9]]  <- lmer(Mn_P ~ (1 | K_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
Mn_P.models[[10]]  <- lmer(Mn_P ~ (1 | Mn_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
Mn_P.models[[11]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman +  (1 | pH) , 
								data = datafix_hadi, REML=FALSE)									
Mn_P.models[[12]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)
Mn_P.models[[13]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman +  (1 | Mn_total) , 
								data = datafix_hadi, REML=FALSE)	
Mn_P.models[[14]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman +  (1 | Zn_total) , 
								data = datafix_hadi, REML=FALSE)	
Mn_P.models[[15]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)								
Mn_P.models[[16]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
Mn_P.models[[17]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)
Mn_P.models[[18]]  <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
Mn_P.models[[19]]  <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
Mn_P.models[[20]]  <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman +(1 | Mn_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
								

##create a vector of names to trace back models in set
Mn_Pnames <- paste("mod", 1:length(Mn_P.models), sep = " ")

##generate AICc table
aictab(cand.set = Mn_P.models, modnames = Mn_Pnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Mn_P.models, modnames = Mn_Pnames, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

				K     AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 12 11  12.0158     0.0000      1      1  14.4207
		mod 3   3  77.5940    65.5782      0      1 -35.2515
		mod 2   3  79.3686    67.3528      0      1 -36.1389
		mod 4   3  79.7018    67.6860      0      1 -36.3054
		mod 1   3  79.8560    67.8402      0      1 -36.3825
		mod 8   4  79.9057    67.8899      0      1 -35.0005
		mod 6   4  80.3764    68.3606      0      1 -35.2358
		mod 10  4  80.4078    68.3920      0      1 -35.2515
		mod 9   4  81.3089    69.2931      0      1 -35.7021
		mod 5   4  82.1825    70.1667      0      1 -36.1389
		mod 7   4  82.4766    70.4608      0      1 -36.2859
		mod 11 11 106.3170    94.3012      0      1 -32.7299
	

AICc(Mn_P.lm)
101.4093

Mn_P.models.fix <- list( )									
Mn_P.models.fix[[1]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi , REML=FALSE)	
Mn_P.models.fix[[2]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
Mn_P.models.fix[[3]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
Mn_P.models.fix[[4]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
Mn_P.models.fix[[5]]   <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	

										
##create a vector of names to trace back models in set
Mn_Pnames.fix <- paste("mod", 1:length(Mn_P.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = Mn_P.models.fix, modnames = Mn_Pnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Mn_P.models.fix, modnames = Mn_Pnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

			   K    AICc Delta_AICc AICcWt Cum.Wt      LL
		mod 1 11 12.0158          0      1      1 14.4207
		
Mn_Pmodelfit.all <- lme4::allFit(Mn_P.models.fix[[1]] )
		bobyqa : [OK]
		Nelder_Mead : [failed]
		nlminbwrap : [OK]
		optimx.L-BFGS-B : [failed]
		nloptwrap.NLOPT_LN_NELDERMEAD : [failed]
		nloptwrap.NLOPT_LN_BOBYQA : [OK]

Mn_P_ss <- summary(Mn_Pmodelfit.all)


##Final Model
Mn_P.model.fix  <- lmer(Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi)	
summary(Mn_P.model.fix)		
		Linear mixed model fit by REML. t-tests use Satterthwaite method ['lmerModLmerTest']
		Formula: Mn_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total)
		   Data: datafix_hadi

		REML criterion at convergence: 59.7

		Scaled residuals: 
			   Min         1Q     Median         3Q        Max 
		-0.0027633 -0.0012276 -0.0003171  0.0012153  0.0032006 

		Random effects:
		 Groups   Name        Variance  Std.Dev.
		 K_total  (Intercept) 1.230e+00 1.108837
		 Residual             3.694e-06 0.001922
		Number of obs: 26, groups:  K_total, 23

		Fixed effects:
													  Estimate Std. Error        df t value Pr(>|t|)
		(Intercept)                                   0.108131   0.495889 17.000000   0.218   0.8300
		Zona_PengelolaanGawangan Mati                 0.438404   0.701294 17.000000   0.625   0.5402
		Zona_PengelolaanPiringan                     -0.216224   0.002718 17.000000 -79.554   <2e-16
		Kedalaman30-60                               -0.372499   0.809782 17.000000  -0.460   0.6513
		Kedalaman60-90                               -0.058252   0.701294 17.000000  -0.083   0.9348
		Zona_PengelolaanGawangan Mati:Kedalaman30-60  0.511191   1.145205 17.000000   0.446   0.6610
		Zona_PengelolaanPiringan:Kedalaman30-60       0.733032   0.846894 17.000000   0.866   0.3988
		Zona_PengelolaanGawangan Mati:Kedalaman60-90 -0.536326   1.280377 17.000000  -0.419   0.6805
		Zona_PengelolaanPiringan:Kedalaman60-90      -1.398550   0.743840 17.000000  -1.880   0.0773
														
		(Intercept)                                     
		Zona_PengelolaanGawangan Mati                   
		Zona_PengelolaanPiringan                     ***
		Kedalaman30-60                                  
		Kedalaman60-90                                  
		Zona_PengelolaanGawangan Mati:Kedalaman30-60    
		Zona_PengelolaanPiringan:Kedalaman30-60         
		Zona_PengelolaanGawangan Mati:Kedalaman60-90    
		Zona_PengelolaanPiringan:Kedalaman60-90      .  
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

		Correlation of Fixed Effects:
					(Intr) Zn_PGM Zn_PnP K30-60 K60-90 Z_PGM:K3 Z_PP:K3 Z_PGM:K6
		Zn_PngllnGM -0.707                                                      
		Zn_PngllnPr -0.003  0.002                                               
		Kedlmn30-60 -0.612  0.433  0.002                                        
		Kedlmn60-90 -0.707  1.000  0.002  0.433                                 
		Z_PGM:K30-6  0.433 -0.612 -0.001 -0.707 -0.612                          
		Z_PP:K30-60  0.000  0.000 -0.003 -0.598  0.000  0.423                   
		Z_PGM:K60-9  0.387 -0.822 -0.001 -0.237 -0.822  0.503    0.000          
		Z_PP:K60-90  0.000 -0.471 -0.004  0.000 -0.471  0.289    0.488   0.516 

anova(Mn_P.model.fix, test = "F")
		Type III Analysis of Variance Table with Satterthwaite method
									   Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
		Zona_Pengelolaan           9.9990e-06 4.9997e-06     2    17  1.3536 0.28478  
		Kedalaman                  1.7677e-05 8.8383e-06     2    17  2.3928 0.12144  
		Zona_Pengelolaan:Kedalaman 3.7395e-05 9.3487e-06     4    17  2.5310 0.07865 .
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

rsq(Mn_P.model.fix, adj='TRUE')
		$model
		[1] 0.9999944

		$fixed
		[1] -0.2665684

		$random
		[1] 1.266563

AICc(Mn_P.model.fix)	
100.5147


# Aktivitas enzim Li-P

Li_P.lm <- lm(Li_P ~ Zona_Pengelolaan * Kedalaman, data = datafix_hadi)	

Li_P.models <- list( )
Li_P.models[[1]]  <- lmer(Li_P ~ (1 | pH), data = datafix_hadi, REML=FALSE)
Li_P.models[[2]]  <- lmer(Li_P ~ (1 | K_total), data = datafix_hadi, REML=FALSE)
Li_P.models[[3]]  <- lmer(Li_P ~ (1 | Mn_total), data = datafix_hadi, REML=FALSE)
Li_P.models[[4]]  <- lmer(Li_P ~ (1 | Zn_total), data = datafix_hadi, REML=FALSE)
Li_P.models[[5]]  <- lmer(Li_P ~ (1 | pH) + (1 | K_total), data = datafix_hadi, REML=FALSE)
Li_P.models[[6]]  <- lmer(Li_P ~ (1 | pH) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
Li_P.models[[7]]  <- lmer(Li_P ~ (1 | pH) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)
Li_P.models[[8]]  <- lmer(Li_P ~ (1 | K_total) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
Li_P.models[[9]]  <- lmer(Li_P ~ (1 | K_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
Li_P.models[[10]]  <- lmer(Li_P ~ (1 | Mn_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
Li_P.models[[11]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman +  (1 | pH) , 
								data = datafix_hadi, REML=FALSE)									
Li_P.models[[12]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)
Li_P.models[[13]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman +  (1 | Mn_total) , 
								data = datafix_hadi, REML=FALSE)	
Li_P.models[[14]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman +  (1 | Zn_total) , 
								data = datafix_hadi, REML=FALSE)	
Li_P.models[[15]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)								
Li_P.models[[16]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
Li_P.models[[17]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)
Li_P.models[[18]]  <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
Li_P.models[[19]]  <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
Li_P.models[[20]]  <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman +(1 | Mn_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
								

##create a vector of names to trace back models in set
Li_Pnames <- paste("mod", 1:length(Li_P.models), sep = " ")

##generate AICc table
aictab(cand.set = Li_P.models, modnames = Li_Pnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Li_P.models, modnames = Li_Pnames, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

				K     AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 12 11   9.0941     0.0000 0.9556 0.9556  15.8815
		mod 19 12  16.3020     7.2079 0.0260 0.9816  15.8490
		mod 16 12  16.9972     7.9030 0.0184 1.0000  15.5014
		mod 13 11  32.4521    23.3580 0.0000 1.0000   4.2025
		mod 18 12  38.9361    29.8419 0.0000 1.0000   4.5320
		mod 14 11  77.3504    68.2562 0.0000 1.0000 -18.2466
		mod 4   3  77.4898    68.3956 0.0000 1.0000 -35.1994
		mod 3   3  79.3749    70.2808 0.0000 1.0000 -36.1420
		mod 1   3  79.7051    70.6110 0.0000 1.0000 -36.3071
		mod 2   3  79.8560    70.7618 0.0000 1.0000 -36.3825
		mod 9   4  80.0890    70.9948 0.0000 1.0000 -35.0921
		mod 7   4  80.1063    71.0122 0.0000 1.0000 -35.1008
		mod 10  4  80.2527    71.1585 0.0000 1.0000 -35.1740
		mod 6   4  82.1888    73.0946 0.0000 1.0000 -36.1420
		mod 8   4  82.1888    73.0946 0.0000 1.0000 -36.1420
		mod 5   4  82.5190    73.4248 0.0000 1.0000 -36.3071
		mod 11 11  99.7166    90.6225 0.0000 1.0000 -29.4297
		mod 17 12 106.8383    97.7441 0.0000 1.0000 -29.4191
		mod 15 12 106.8595    97.7653 0.0000 1.0000 -29.4297

AICc(Li_P.lm)
97.52085

Li_P.models.fix <- list( )									
Li_P.models.fix[[1]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi , REML=FALSE)	
Li_P.models.fix[[2]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
Li_P.models.fix[[3]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
Li_P.models.fix[[4]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
Li_P.models.fix[[5]]   <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	

										
##create a vector of names to trace back models in set
Li_Pnames.fix <- paste("mod", 1:length(Li_P.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = Li_P.models.fix, modnames = Li_Pnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Li_P.models.fix, modnames = Li_Pnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

			   K    AICc Delta_AICc AICcWt Cum.Wt      LL
		mod 1 11  9.0941     0.0000 0.4906 0.4906 15.8815
		mod 4 11  9.1939     0.0998 0.4668 0.9574 15.8316
		mod 3 11 14.0174     4.9232 0.0418 0.9992 13.4199
		mod 2 11 22.0569    12.9627 0.0008 1.0000  9.4001
		
Li_Pmodelfit.all <- lme4::allFit(Li_P.models.fix[[1]] )
		bobyqa : [OK]
		Nelder_Mead : [failed]
		nlminbwrap : [OK]
		optimx.L-BFGS-B : [OK]
		nloptwrap.NLOPT_LN_NELDERMEAD : [failed]
		nloptwrap.NLOPT_LN_BOBYQA : [OK]

Li_P_ss <- summary(Li_Pmodelfit.all)


##Final Model
Li_P.model.fix  <- lmer(Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
										data = datafix_hadi)	
summary(Li_P.model.fix)		
		Linear mixed model fit by REML. t-tests use Satterthwaite method ['lmerModLmerTest']
		Formula: Li_P ~ Zona_Pengelolaan * Kedalaman + (1 | K_total)
		   Data: datafix_hadi

		REML criterion at convergence: 56.5

		Scaled residuals: 
			   Min         1Q     Median         3Q        Max 
		-0.0026108 -0.0013100 -0.0000794  0.0009669  0.0033299 

		Random effects:
		 Groups   Name        Variance  Std.Dev.
		 K_total  (Intercept) 1.020e+00 1.00979 
		 Residual             3.687e-06 0.00192 
		Number of obs: 26, groups:  K_total, 23

		Fixed effects:
													  Estimate Std. Error        df t value Pr(>|t|)
		(Intercept)                                  -0.610083   0.451593 17.000000  -1.351    0.194
		Zona_PengelolaanGawangan Mati                 1.534096   0.638649 17.000000   2.402    0.028
		Zona_PengelolaanPiringan                      0.958154   0.002715 17.000000 352.854   <2e-16
		Kedalaman30-60                               -0.430139   0.737446 17.000000  -0.583    0.567
		Kedalaman60-90                                0.704998   0.638649 17.000000   1.104    0.285
		Zona_PengelolaanGawangan Mati:Kedalaman30-60 -0.147535   1.042906 17.000000  -0.141    0.889
		Zona_PengelolaanPiringan:Kedalaman30-60      -0.704699   0.771244 17.000000  -0.914    0.374
		Zona_PengelolaanGawangan Mati:Kedalaman60-90 -1.443964   1.166003 17.000000  -1.238    0.232
		Zona_PengelolaanPiringan:Kedalaman60-90      -0.340108   0.677395 17.000000  -0.502    0.622
														
		(Intercept)                                     
		Zona_PengelolaanGawangan Mati                *  
		Zona_PengelolaanPiringan                     ***
		Kedalaman30-60                                  
		Kedalaman60-90                                  
		Zona_PengelolaanGawangan Mati:Kedalaman30-60    
		Zona_PengelolaanPiringan:Kedalaman30-60         
		Zona_PengelolaanGawangan Mati:Kedalaman60-90    
		Zona_PengelolaanPiringan:Kedalaman60-90         
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

		Correlation of Fixed Effects:
					(Intr) Zn_PGM Zn_PnP K30-60 K60-90 Z_PGM:K3 Z_PP:K3 Z_PGM:K6
		Zn_PngllnGM -0.707                                                      
		Zn_PngllnPr -0.003  0.002                                               
		Kedlmn30-60 -0.612  0.433  0.002                                        
		Kedlmn60-90 -0.707  1.000  0.002  0.433                                 
		Z_PGM:K30-6  0.433 -0.612 -0.001 -0.707 -0.612                          
		Z_PP:K30-60  0.000  0.000 -0.004 -0.598  0.000  0.423                   
		Z_PGM:K60-9  0.387 -0.822 -0.001 -0.237 -0.822  0.503    0.000          
		Z_PP:K60-90  0.000 -0.471 -0.004  0.000 -0.471  0.289    0.488   0.516  
		optimizer (nloptwrap) convergence code: 0 (OK)
		unable to evaluate scaled gradient
		Model failed to converge: degenerate  Hessian with 1 negative eigenvalues


anova(Li_P.model.fix, test = "F")
		Type III Analysis of Variance Table with Satterthwaite method
									   Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
		Zona_Pengelolaan           2.9789e-05 1.4895e-05     2    17  4.0399 0.03669 *
		Kedalaman                  1.8791e-05 9.3956e-06     2    17  2.5484 0.10765  
		Zona_Pengelolaan:Kedalaman 1.9729e-05 4.9323e-06     4    17  1.3378 0.29643  
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

rsq(Li_P.model.fix, adj='TRUE')
$model
[1] 0.9999944

$fixed
[1] -0.1104761

$random
[1] 1.11047

AICc(Li_P.model.fix)	
97.33325

## Aktivitas enzim (bakteri) selulase

EBS.lm <- lm(EBS ~ Zona_Pengelolaan * Kedalaman, data = datafix_hadi)	

EBS.models <- list( )
EBS.models[[1]]  <- lmer(EBS ~ (1 | pH), data = datafix_hadi, REML=FALSE)
EBS.models[[2]]  <- lmer(EBS ~ (1 | K_total), data = datafix_hadi, REML=FALSE)
EBS.models[[3]]  <- lmer(EBS ~ (1 | Mn_total), data = datafix_hadi, REML=FALSE)
EBS.models[[4]]  <- lmer(EBS ~ (1 | Zn_total), data = datafix_hadi, REML=FALSE)
EBS.models[[5]]  <- lmer(EBS ~ (1 | pH) + (1 | K_total), data = datafix_hadi, REML=FALSE)
EBS.models[[6]]  <- lmer(EBS ~ (1 | pH) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
EBS.models[[7]]  <- lmer(EBS ~ (1 | pH) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)
EBS.models[[8]]  <- lmer(EBS ~ (1 | K_total) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
EBS.models[[9]]  <- lmer(EBS ~ (1 | K_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
EBS.models[[10]]  <- lmer(EBS ~ (1 | Mn_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
EBS.models[[11]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman +  (1 | pH) , 
								data = datafix_hadi, REML=FALSE)									
EBS.models[[12]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)
EBS.models[[13]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman +  (1 | Mn_total) , 
								data = datafix_hadi, REML=FALSE)	
EBS.models[[14]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman +  (1 | Zn_total) , 
								data = datafix_hadi, REML=FALSE)	
EBS.models[[15]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)								
EBS.models[[16]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
EBS.models[[17]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)
EBS.models[[18]]  <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
EBS.models[[19]]  <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
EBS.models[[20]]  <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman +(1 | Mn_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
								

##create a vector of names to trace back models in set
EBSnames <- paste("mod", 1:length(EBS.models), sep = " ")

##generate AICc table
aictab(cand.set = EBS.models, modnames = EBSnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = EBS.models, modnames = EBSnames, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

				K    AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 4   3 78.6791     0.0000 0.2572 0.2572 -35.7941    ## bermasalah semua modelnya
		mod 2   3 79.7615     1.0824 0.1497 0.4068 -36.3353
		mod 1   3 79.8560     1.1769 0.1428 0.5496 -36.3825
		mod 3   3 79.8560     1.1769 0.1428 0.6924 -36.3825
		mod 10  4 81.2482     2.5691 0.0712 0.7636 -35.6717
		mod 9   4 81.4253     2.7462 0.0651 0.8287 -35.7603
		mod 7   4 81.4930     2.8139 0.0630 0.8917 -35.7941
		mod 8   4 82.5753     3.8962 0.0367 0.9284 -36.3353
		mod 5   4 82.5753     3.8962 0.0367 0.9650 -36.3353
		mod 6   4 82.6698     3.9907 0.0350 1.0000 -36.3825
		mod 11 11 97.3680    18.6888 0.0000 1.0000 -28.2554

AICc(EBS.lm)
103.0256

# simulasi dijalankan sesuai AICc terendah 

EBS.models.fix <- list( )									
EBS.models.fix[[1]]   <- lmer(EBS ~  (1 | Zn_total), 
										data = datafix_hadi , REML=FALSE)	
EBS.models.fix[[2]]   <- lmer(EBS ~  (1 | Zn_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
EBS.models.fix[[3]]   <- lmer(EBS ~  (1 | Zn_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
EBS.models.fix[[4]]   <- lmer(EBS ~  (1 | Zn_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
EBS.models.fix[[5]]   <- lmer(EBS ~  (1 | Zn_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	
##create a vector of names to trace back models in set
EBSnames.fix <- paste("mod", 1:length(EBS.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = EBS.models.fix, modnames = EBSnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = EBS.models.fix, modnames = EBSnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:
Model selection based on AICc:

      K    AICc Delta_AICc AICcWt Cum.Wt       LL
mod 1 3 78.6791     0.0000 0.2009 0.2009 -35.7941
mod 3 3 78.6791     0.0000 0.2009 0.4018 -35.7941
mod 4 3 78.6791     0.0000 0.2009 0.6028 -35.7941
mod 2 3 78.6791     0.0000 0.2009 0.8037 -35.7941
mod 5 3 78.7253     0.0462 0.1963 1.0000 -35.8172



EBS.models.fix <- list( )									
EBS.models.fix[[1]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi , REML=FALSE)	
EBS.models.fix[[2]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
EBS.models.fix[[3]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
EBS.models.fix[[4]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
EBS.models.fix[[5]]   <- lmer(EBS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	

										
##create a vector of names to trace back models in set
EBSnames.fix <- paste("mod", 1:length(EBS.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = EBS.models.fix, modnames = EBSnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = EBS.models.fix, modnames = EBSnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

			   K   AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 4 11 97.368          0    0.2    0.2 -28.2554
		mod 1 11 97.368          0    0.2    0.4 -28.2554
		mod 2 11 97.368          0    0.2    0.6 -28.2554
		mod 5 11 97.368          0    0.2    0.8 -28.2554
		mod 3 11 97.368          0    0.2    1.0 -28.2554
		
EBSmodelfit.all <- lme4::allFit(EBS.models.fix[[1]] )


EBS_ss <- summary(EBSmodelfit.all)


##Final Model
EBS.model.fix  <- lmer(EBS ~ (1 | Zn_total), 
										data = datafix_hadi)	
summary(EBS.model.fix)		
		Linear mixed model fit by REML. t-tests use Satterthwaite method [
		lmerModLmerTest]
		Formula: EBS ~ (1 | Zn_total)
		   Data: datafix_hadi

		REML criterion at convergence: 73

		Scaled residuals: 
			 Min       1Q   Median       3Q      Max 
		-1.06943 -0.22954 -0.01755  0.36154  0.82775 

		Random effects:
		 Groups   Name        Variance Std.Dev.
		 Zn_total (Intercept) 0.7439   0.8625  
		 Residual             0.2433   0.4932  
		Number of obs: 26, groups:  Zn_total, 25

		Fixed effects:
					Estimate Std. Error       df t value Pr(>|t|)
		(Intercept)  0.03855    0.19815 24.14317   0.195    0.847

## tidak dipakai karena tidak ada Fixed variabel --> ini pakai random var pH
anova(EBS.model.fix, test = "F")
		Type III Analysis of Variance Table with Satterthwaite method
								   Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
		Zona_Pengelolaan           1.5504 0.77521     2     2  4.8342 0.1714
		Kedalaman                  0.8924 0.44620     2     2  2.7825 0.2644
		Zona_Pengelolaan:Kedalaman 2.8239 0.70597     4     2  4.4024 0.1936

rsq(EBS.model.fix, adj='TRUE')
		[1] 0.7485406

		$fixed
		[1] -0.001545824

		$random
		[1] 0.7500864

AICc(EBS.model.fix)	
80.09912



## Aktivitas enzim (fungi) selulase

EFS.lm <- lm(EFS ~ Zona_Pengelolaan * Kedalaman, data = datafix_hadi)	

EFS.models <- list( )
EFS.models[[1]]  <- lmer(EFS ~ (1 | pH), data = datafix_hadi, REML=FALSE)
EFS.models[[2]]  <- lmer(EFS ~ (1 | K_total), data = datafix_hadi, REML=FALSE)
EFS.models[[3]]  <- lmer(EFS ~ (1 | Mn_total), data = datafix_hadi, REML=FALSE)
EFS.models[[4]]  <- lmer(EFS ~ (1 | Zn_total), data = datafix_hadi, REML=FALSE)
EFS.models[[5]]  <- lmer(EFS ~ (1 | pH) + (1 | K_total), data = datafix_hadi, REML=FALSE)
EFS.models[[6]]  <- lmer(EFS ~ (1 | pH) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
EFS.models[[7]]  <- lmer(EFS ~ (1 | pH) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)
EFS.models[[8]]  <- lmer(EFS ~ (1 | K_total) + (1 | Mn_total), data = datafix_hadi, REML=FALSE)
EFS.models[[9]]  <- lmer(EFS ~ (1 | K_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
EFS.models[[10]]  <- lmer(EFS ~ (1 | Mn_total) + (1 | Zn_total), data = datafix_hadi, REML=FALSE)	
EFS.models[[11]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman +  (1 | pH) , 
								data = datafix_hadi, REML=FALSE)									
EFS.models[[12]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)
EFS.models[[13]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman +  (1 | Mn_total) , 
								data = datafix_hadi, REML=FALSE)	
EFS.models[[14]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman +  (1 | Zn_total) , 
								data = datafix_hadi, REML=FALSE)	
EFS.models[[15]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | K_total), 
								data = datafix_hadi, REML=FALSE)								
EFS.models[[16]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
EFS.models[[17]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)
EFS.models[[18]]  <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Mn_total), 
								data = datafix_hadi, REML=FALSE)
EFS.models[[19]]  <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman +(1 | K_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
EFS.models[[20]]  <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman +(1 | Mn_total) + (1 | Zn_total), 
								data = datafix_hadi, REML=FALSE)	
								

##create a vector of names to trace back models in set
EFSnames <- paste("mod", 1:length(EFS.models), sep = " ")

##generate AICc table
aictab(cand.set = EFS.models, modnames = EFSnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = EFS.models, modnames = EFSnames, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

				K     AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 3   3  78.8439     0.0000 0.2082 0.2082 -35.8765 
		mod 2   3  78.8598     0.0159 0.2065 0.4147 -35.8844
		mod 4   3  79.5299     0.6859 0.1477 0.5624 -36.2195
		mod 1   3  79.8560     1.0120 0.1255 0.6879 -36.3825
		mod 8   4  81.1074     2.2635 0.0671 0.7550 -35.6013
		mod 10  4  81.4458     2.6018 0.0567 0.8117 -35.7705
		mod 6   4  81.6578     2.8139 0.0510 0.8627 -35.8765
		mod 5   4  81.6737     2.8297 0.0506 0.9132 -35.8844
		mod 9   4  81.6737     2.8297 0.0506 0.9638 -35.8844
		mod 7   4  82.3437     3.4998 0.0362 1.0000 -36.2195
		mod 11 11 107.1665    28.3226 0.0000 1.0000 -33.1547


AICc(EFS.lm)
100.976

EFS.models.fix <- list( )									
EFS.models.fix[[1]]   <- lmer(EFS ~  (1 | Mn_total), 
										data = datafix_hadi , REML=FALSE)	
EFS.models.fix[[2]]   <- lmer(EFS ~  (1 | Mn_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
EFS.models.fix[[3]]   <- lmer(EFS ~ (1 | Mn_total), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
EFS.models.fix[[4]]   <- lmer(EFS ~ (1 | Mn_total), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
EFS.models.fix[[5]]   <- lmer(EFS ~ (1 | Mn_total), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	
										
EFSnames.fix <- paste("mod", 1:length(EFS.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = EFS.models.fix, modnames = EFSnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = EFS.models.fix, modnames = EFSnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)

			Model selection based on AICc:

				  K    AICc Delta_AICc AICcWt Cum.Wt       LL
			mod 2 3 78.8439          0    0.2    0.2 -35.8765
			mod 4 3 78.8439          0    0.2    0.4 -35.8765
			mod 1 3 78.8439          0    0.2    0.6 -35.8765
			mod 3 3 78.8439          0    0.2    0.8 -35.8765
			mod 5 3 78.8439          0    0.2    1.0 -35.8765


EFS.models.fix <- list( )									
EFS.models.fix[[1]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi , REML=FALSE)	
EFS.models.fix[[2]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), REML=FALSE)								
EFS.models.fix[[3]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi, control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='nlminb')), REML=FALSE)	
EFS.models.fix[[4]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi, control=lmerControl(optimizer="bobyqa",
										optCtrl=list(maxfun=2e5)), REML=FALSE)
EFS.models.fix[[5]]   <- lmer(EFS ~ Zona_Pengelolaan * Kedalaman + (1 | pH), 
										data = datafix_hadi, control = lmerControl(optimizer ="Nelder_Mead"),
										REML=FALSE)	

										
##create a vector of names to trace back models in set
EFSnames.fix <- paste("mod", 1:length(EFS.models.fix), sep = " ")

##generate AICc table
aictab(cand.set = EFS.models.fix, modnames = EFSnames.fix, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = EFS.models.fix, modnames = EFSnames.fix, sort = TRUE),
      digits = 4, LL = TRUE)
		Model selection based on AICc:

			   K     AICc Delta_AICc AICcWt Cum.Wt       LL
		mod 1 11 107.1665          0    0.2    0.2 -33.1547
		mod 2 11 107.1665          0    0.2    0.4 -33.1547
		mod 3 11 107.1665          0    0.2    0.6 -33.1547
		mod 4 11 107.1665          0    0.2    0.8 -33.1547
		mod 5 11 107.1665          0    0.2    1.0 -33.1547

EFSmodelfit.all <- lme4::allFit(EFS.models.fix[[1]] )


EFS_ss <- summary(EFSmodelfit.all)


##Final Model
EFS.model.fix  <- lmer(EFS ~ (1 | Mn_total),control = lmerControl(
										optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), 
										data = datafix_hadi)	
summary(EFS.model.fix)		
		Linear mixed model fit by REML. t-tests use Satterthwaite method [
		lmerModLmerTest]
		Formula: EFS ~ (1 | Mn_total)
		   Data: datafix_hadi
		Control: lmerControl(optimizer = "optimx", optCtrl = list(method = "L-BFGS-B"))

		REML criterion at convergence: 73.1

		Scaled residuals: 
			 Min       1Q   Median       3Q      Max 
		-1.26664 -0.37953 -0.07467  0.42834  1.35849 

		Random effects:
		 Groups   Name        Variance Std.Dev.
		 Mn_total (Intercept) 0.7324   0.8558  
		 Residual             0.3085   0.5554  
		Number of obs: 26, groups:  Mn_total, 23

		Fixed effects:
					Estimate Std. Error       df t value Pr(>|t|)
		(Intercept) -0.05764    0.21036 21.26940  -0.274    0.787


anova(EFS.model.fix, test = "F")
		Type III Analysis of Variance Table with Satterthwaite method
								   Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
		Zona_Pengelolaan           3.4119 1.70593     2    17  1.4870 0.2540
		Kedalaman                  0.7211 0.36057     2    17  0.3143 0.7345
		Zona_Pengelolaan:Kedalaman 1.1506 0.28766     4    17  0.2507 0.9052

rsq(EFS.model.fix, adj='TRUE')
		$model
		[1] 0.6861151

		$fixed
		[1] -0.003455559

		$random
		[1] 0.6895707

AICc(EFS.model.fix)	
[1] 80.14712



## Model Bermasalah: PFS, LKS, EBS, EFS --> coba pakai pendekatan lain lain 

Data_Penelitian <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "R")
View(Data_Penelitian)
str(Data_Penelitian)
	'data.frame':	26 obs. of  22 variables:
	 $ Zona_Pengelolaan: chr  "Piringan" "Piringan" "Piringan" "Piringan" ...
	 $ Kedalaman       : chr  "0-30" "30-60" "60-90" "0-30" ...
	 $ Kadar_air       : num  352 723 1084 417 558 ...
	 $ pH              : num  3.36 3.29 3.3 3.35 3.57 3.26 3.71 3.28 3.42 3.32 ...
	 $ C_organik       : num  47.8 54.8 53.2 52.2 46.5 ...
	 $ Kadar_abu       : num  17.65 5.46 8.23 9.92 19.79 ...
	 $ CN_rasio        : num  57.1 87.9 142.2 123.7 115.7 ...
	 $ N_total         : num  0.836 0.624 0.374 0.422 0.402 ...
	 $ P_total         : num  56.5 52.9 46.3 42.1 59.2 ...
	 $ K_total         : num  55.3 63.3 61.3 54.3 58.3 ...
	 $ Fe_total        : num  1661 477 492 870 1366 ...
	 $ Mn_total        : num  33.86 11.75 6.1 6.87 15.41 ...
	 $ Cu_total        : num  5.78 7.18 7.91 3.37 2.97 5.44 6.78 5.33 3.59 2.36 ...
	 $ Zn_total        : num  13.46 13.14 8.85 8.94 8.85 ...
	 $ PBL             : num  6.60e+06 9.65e+06 1.48e+07 1.50e+08 3.44e+08 ...
	 $ PBS             : num  8700000 12150000 10350000 6700000 23400000 ...
	 $ PFS             : num  6818 4950 4265 5970 3167 ...
	 $ LKS             : num  208 521 255 197 4873 ...
	 $ Mn_P            : num  20.7 27.5 34.4 41.3 20.7 ...
	 $ Li_P            : num  1245.5 17.9 80.6 152.3 116.5 ...
	 $ EBS             : num  50.741 12.222 3.056 0.648 1.111 ...
	 $ EFS             : num  9.63 4.54 5.93 9.44 2.96 ...
	 - attr(*, "na.action")= 'omit' Named int 9
	  ..- attr(*, "names")= chr "9"


## LME ReML
Data_Penelitian <- as.data.frame(Data_Penelitian)
#selecting data

datafact<- dplyr::select(Data_Penelitian, c(1:2)) #quantitative selecting column
datadeptv <- dplyr::select(Data_Penelitian, c(15:22)) #quantitative selecting column
dataindeptv <- dplyr::select(Data_Penelitian, c(3:14)) #quantitative selecting column

#assign IDs
datafact <- tibble::rowid_to_column(datafact, "ID")
datadeptv <- tibble::rowid_to_column(datadeptv, "ID")
data1 <- merge(datafact,datadeptv,by="ID")

# box cox transform
baru <- sapply(dataindeptv,\(x)boxcox(x)$x.t)
baru

df <- data.frame(baru)
df <- tibble::rowid_to_column(df, "ID")


datafix_hadi1 <- merge(data1,df,by="ID")

datafix_hadi1$Zona_Pengelolaan <- as.factor(datafix_hadi1$Zona_Pengelolaan)
datafix_hadi1$Kedalaman        <- as.factor(datafix_hadi1$Kedalaman)



## Tes binomial logit untuk PFS (mungkin bisa berhasil)
PFS.test.models.glm.pois <- glm(PFS ~ Zona_Pengelolaan * Kedalaman ,
						data = datafix_hadi1, family=poisson(link="log"))
PFS.test.models.glm.qpois <- glm(PFS ~ Zona_Pengelolaan * Kedalaman ,
						data = datafix_hadi1, family=quasipoisson(link="log"))
PFS.test.models.glm.nb <- glm.nb(PFS ~ Zona_Pengelolaan * Kedalaman ,
						data = datafix_hadi1)


PFS.test.models <- list( )	
PFS.test.models[[1]] <- glmer.nb(PFS ~ Zona_Pengelolaan * Kedalaman  + (1 | pH),
						data = Data_Penelitian)
PFS.test.models[[2]] <- glmer.nb(PFS ~ Zona_Pengelolaan * Kedalaman  + (1 | K_total),
						data = Data_Penelitian)
PFS.test.models[[3]] <- glmer.nb(PFS ~ Zona_Pengelolaan * Kedalaman  + (1 | Mn_total),
						data = Data_Penelitian)
PFS.test.models[[4]] <- glmer.nb(PFS ~ Zona_Pengelolaan * Kedalaman  + (1 | Cu_total),
						data = Data_Penelitian)
PFS.test.models[[5]] <- glmer.nb(PFS ~ Zona_Pengelolaan * Kedalaman  + (1 | Zn_total),
						data = Data_Penelitian)


## Final model --> No[[4]]
summary (PFS.test.models[[4]])
Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: Negative Binomial(43.3037)  ( log )
Formula: PFS ~ Zona_Pengelolaan * Kedalaman + (1 | Cu_total)
   Data: Data_Penelitian

     AIC      BIC   logLik deviance df.resid 
   527.9    541.7   -252.9    505.9       15 

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-1.22894 -0.09670 -0.00319  0.10775  1.10236 

Random effects:
 Groups   Name        Variance Std.Dev.
 Cu_total (Intercept) 0.6466   0.8041  
Number of obs: 26, groups:  Cu_total, 23

Fixed effects:
                                             Estimate Std. Error z value Pr(>|z|)    
(Intercept)                                   8.61758    0.46175  18.663   <2e-16 ***
Zona_PengelolaanGawangan Mati                 0.43835    0.73904   0.593   0.5531    
Zona_PengelolaanPiringan                     -0.07306    0.64822  -0.113   0.9103    
Kedalaman30-60                                0.90162    0.59881   1.506   0.1321    
Kedalaman60-90                                0.25764    0.59089   0.436   0.6628    
Zona_PengelolaanGawangan Mati:Kedalaman30-60 -1.85760    0.91664  -2.027   0.0427 *  
Zona_PengelolaanPiringan:Kedalaman30-60      -1.15892    0.88836  -1.305   0.1920    
Zona_PengelolaanGawangan Mati:Kedalaman60-90 -1.08017    0.91011  -1.187   0.2353    
Zona_PengelolaanPiringan:Kedalaman60-90       1.13513    0.92677   1.225   0.2206    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) Zn_PGM Zn_PnP K30-60 K60-90 Z_PGM:K3 Z_PP:K3 Z_PGM:K6
Zn_PngllnGM -0.626                                                      
Zn_PngllnPr -0.694  0.435                                               
Kedlmn30-60 -0.772  0.483  0.536                                        
Kedlmn60-90 -0.770  0.482  0.530  0.594                                 
Z_PGM:K30-6  0.511 -0.810 -0.358 -0.659 -0.171                          
Z_PP:K30-60  0.507 -0.317 -0.730 -0.664 -0.387  0.440                   
Z_PGM:K60-9  0.500 -0.807 -0.344 -0.162 -0.650  0.363    0.100          
Z_PP:K60-90  0.469 -0.294 -0.678 -0.362 -0.613  0.103    0.495   0.398 


AICc(PFS.test.models[[4]])
[1] 546.7352

anova(PFS.test.models[[4]], test = "F")
		Analysis of Variance Table
								   npar Sum Sq Mean Sq F value
		Zona_Pengelolaan              2 44.833 22.4167 22.4167
		Kedalaman                     2  1.743  0.8717  0.8717
		Zona_Pengelolaan:Kedalaman    4 10.848  2.7120  2.7120						
						
				
