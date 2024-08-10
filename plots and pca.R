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


R_PCA_en <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "R_PCA_en")
View(R_PCA_en)




R_PCA_en <- na.omit(R_PCA_en)						
R_PCA_en$Management_Zone <- as.factor(R_PCA_en$Management_Zone)
R_PCA_en$Sampling_Depth <- as.factor(R_PCA_en$Sampling_Depth)			

str(R_PCA_en)
		tibble [21 × 15] (S3: tbl_df/tbl/data.frame)
		 $ Management_Zone: Factor w/ 3 levels "FRS","FTC","HVP": 2 2 2 2 2 2 2 1 1 1 ...
		 $ Sampling_Depth : Factor w/ 3 levels "0-30","30-60",..: 1 2 3 1 3 2 3 1 2 3 ...
		 $ PWC            : num [1:21] 352 723 1084 417 962 ...
		 $ pH             : num [1:21] 3.36 3.29 3.3 3.35 3.26 3.28 3.33 3.42 3.32 3.33 ...
		 $ Org_C          : num [1:21] 47.8 54.8 53.2 52.2 55.6 ...
		 $ PAC            : num [1:21] 17.65 5.46 8.23 9.92 4.09 ...
		 $ CNR            : num [1:21] 57.1 87.9 142.2 123.7 78 ...
		 $ TN             : num [1:21] 0.836 0.624 0.374 0.422 0.713 ...
		 $ TP             : num [1:21] 56.5 52.9 46.3 42.1 62.2 ...
		 $ TK             : num [1:21] 55.3 63.3 61.3 54.3 58.3 ...
		 $ TFe            : num [1:21] 1661 477 492 870 702 ...
		 $ TMn            : num [1:21] 33.86 11.75 6.1 6.87 3.36 ...
		 $ TCu            : num [1:21] 5.78 7.18 7.91 3.37 5.44 ...
		 $ TZn            : num [1:21] 13.46 13.14 8.85 8.94 11.19 ...
		 $ RH             : num [1:21] 14.4 11.66 6.17 15.77 8.23 ...
		 - attr(*, "na.action")= 'omit' Named int [1:6] 5 7 13 16 18 22
		  ..- attr(*, "names")= chr [1:6] "5" "7" "13" "16" ..

soil.pca <- PCA(R_PCA_en[,-(1:2)], graph = FALSE)

soil.pca$eig #2
				  eigenvalue percentage of variance cumulative percentage of variance
		comp 1  5.315525e+00           4.088866e+01                          40.88866
		comp 2  2.651957e+00           2.039967e+01                          61.28833
		comp 3  1.585835e+00           1.219873e+01                          73.48706
		comp 4  9.708829e-01           7.468330e+00                          80.95539
		comp 5  7.389635e-01           5.684334e+00                          86.63972
		comp 6  5.586245e-01           4.297112e+00                          90.93684
		comp 7  3.984972e-01           3.065363e+00                          94.00220
		comp 8  2.777085e-01           2.136219e+00                          96.13842
		comp 9  2.362316e-01           1.817166e+00                          97.95558
		comp 10 1.654167e-01           1.272436e+00                          99.22802
		comp 11 6.788233e-02           5.221718e-01                          99.75019
		comp 12 3.247505e-02           2.498081e-01                         100.00000
		comp 13 8.833779e-31           6.795214e-30                         100.00000


bar_dim_soil <- fviz_eig(soil.pca, ncp = 10, addlabels = TRUE, barfill="grey", barcolor="grey") #3
bar_dim_soil1 <- bar_dim_soil + theme1_pca

ctrb_dim1_soil1 <- fviz_contrib(soil.pca,choice='var',top=10, fill="grey", color="grey", title = NULL)
dimsoil1_1<-ctrb_dim1_soil1+theme2_pca
ctrb_dimsoil2_1  <- fviz_contrib(soil.pca,choice='var',top=10,axes=2, fill="grey", color="grey", title = NULL)
dimsoil2_1 <-ctrb_dimsoil2_1+theme2_pca

dimsoil1_1+dimsoil2_1

windowsFonts(Times=windowsFont("Times New Roman"))
soil.pcavar <- fviz_pca_var(soil.pca,
             col.var = "contrib", # Color by contributions to the PC
            # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
			 
plotsoil_pcavar <- soil.pcavar+
								theme(text=element_text(family="Times"))



plotsoil_pcavar_Fix <- ggpubr::ggpar(plotsoil_pcavar,
                                  xlab = "PC1 (40.9%)", ylab = "PC2 (20.4%)")



biplot_soil_Zonasi <- fviz_pca_biplot(soil.pca, 
                                col.ind = R_PCA_en$Management_Zone, 
                                palette = "NULL",
                                addEllipses = TRUE,
								ellipse.level = 0.95, 
								ellipse.type = c("confidence"),
                                label = "var",
                                col.var = "black", 
                                repel = TRUE,
                                legend.title = "Management Zone")+ ggtitle ("") + 
								theme(text=element_text(family="Times"))

biplot_soil_Zonasi_fix <- ggpubr::ggpar(biplot_soil_Zonasi,
                                 xlab = "PC1 (40.9%)", ylab = "")+ 
								 theme(text=element_text(family="Times"))
biplot_soil_Zonasi_fix 

biplot_soil_Kdl <- fviz_pca_biplot(soil.pca, 
                                col.ind = R_PCA_en$Sampling_Depth, 
                                palette = "NULL",
                                addEllipses = TRUE,
								ellipse.level = 0.95, 
								ellipse.type = c("confidence"),
                                label = "var",
                                col.var = "black", 
                                repel = TRUE,
                                legend.title = "Sampling Depth (cm)")+ ggtitle ("")

biplot_soil_Kdl_fix <- ggpubr::ggpar(biplot_soil_Kdl,
                                 xlab = "PC1 (40.9%)", ylab = "PC2 (20.4%)")+ 
								 theme(text=element_text(family="Times"))
biplot_soil_Kdl_fix

biplot_soil_Kdl_fix+  theme(legend.position="bottom")+biplot_soil_Zonasi_fix + theme(legend.position="bottom")


R_PCA_en2 <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "R_PCA_en2")
View(R_PCA_en2)

R_PCA_en2 <- na.omit(R_PCA_en2)						
R_PCA_en2$Management_Zone <- as.factor(R_PCA_en2$Management_Zone)
R_PCA_en2$Sampling_Depth <- as.factor(R_PCA_en2$Sampling_Depth)			

str(R_PCA_en2)
		tibble [11 × 10] (S3: tbl_df/tbl/data.frame)
		 $ Management_Zone: Factor w/ 3 levels "FRS","FTC","HVP": 2 2 2 2 1 1 1 3 3 3 ...
		 $ Sampling_Depth : Factor w/ 3 levels "0-30","30-60",..: 1 3 1 2 1 2 3 2 1 2 ...
		 $ PBS            : num [1:11] 8700000 10350000 6700000 200000 23800000 ...
		 $ PFS            : num [1:11] 6818 4265 5970 4000 5856 ...
		 $ LKS            : num [1:11] 208 255 197 185 197 ...
		 $ Mn_P           : num [1:11] 20.7 34.4 41.3 20.7 110.2 ...
		 $ Li_P           : num [1:11] 1245.5 80.6 152.3 976.7 779.6 ...
		 $ EBS            : num [1:11] 50.741 3.056 0.648 5.093 36.019 ...
		 $ EFS            : num [1:11] 9.63 5.93 9.44 9.07 2.5 ...
		 $ RH             : num [1:11] 14.4 6.17 15.77 7.54 5.49 ...
		 - attr(*, "na.action")= 'omit' Named int [1:16] 2 5 6 7 9 11 12 13 16 17 ...
		  ..- attr(*, "names")= chr [1:16] "2" "5" "6" "7" ...
		

soil.pca2 <- PCA(R_PCA_en2[,-(1:2)], graph = FALSE)

soil.pca2$eig #2

bar_dim_soil2 <- fviz_eig(soil.pca2, ncp = 10, addlabels = TRUE, barfill="grey", barcolor="grey") #3
bar_dim_soil12 <- bar_dim_soil2 + theme1_pca
windowsFonts(Times=windowsFont("Times New Roman"))
soil.pcavar2 <- fviz_pca_var(soil.pca2,
             col.var = "contrib", # Color by contributions to the PC
            # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
			 
plotsoil_pcavar2 <- soil.pcavar2+
								theme(text=element_text(family="Times"))



plotsoil_pcavar_Fix2 <- ggpubr::ggpar(plotsoil_pcavar2,
                                  xlab = "PC1 (31.5%)", ylab = "PC2 (24.7%)")



biplot_soil_Zonasi2 <- fviz_pca_biplot(soil.pca2, 
                                col.ind = R_PCA_en2$Management_Zone, 
                                palette = "NULL",
                                addEllipses = TRUE,
								ellipse.level = 0.95, 
								ellipse.type = c("confidence"),
                                label = "var",
                                col.var = "black", 
                                repel = TRUE,
                                legend.title = "Management Zone")+ ggtitle ("") + 
								theme(text=element_text(family="Times"))

biplot_soil_Zonasi_fix2 <- ggpubr::ggpar(biplot_soil_Zonasi2,
                                 xlab = "PC1 (31.5%)", ylab = "")+ 
								 theme(text=element_text(family="Times"))
biplot_soil_Zonasi_fix2 

biplot_soil_Kdl2 <- fviz_pca_biplot(soil.pca2, 
                                col.ind = R_PCA_en2$Sampling_Depth, 
                                palette = "NULL",
                                addEllipses = TRUE,
								ellipse.level = 0.95, 
								ellipse.type = c("confidence"),
                                label = "var",
                                col.var = "black", 
                                repel = TRUE,
                                legend.title = "Sampling Depth (cm)")+ ggtitle ("")

biplot_soil_Kdl_fix2 <- ggpubr::ggpar(biplot_soil_Kdl2,
                                 xlab = "PC1 (31.5%)", ylab = "PC2 (24.7%)")+ 
								 theme(text=element_text(family="Times"))
biplot_soil_Kdl_fix2

biplot_soil_Kdl_fix2+  theme(legend.position="bottom")+biplot_soil_Zonasi_fix2 + theme(legend.position="bottom")



library(readxl)
R_boxplot1 <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "R_boxplot1")
View(R_boxplot1)


R_boxplot1 <- na.omit(R_boxplot1)

R_boxplot1$ID <- factor(R_boxplot1$ID,levels=unique(R_boxplot1$ID))

str(R_boxplot1)
		tibble [318 × 5] (S3: tbl_df/tbl/data.frame)
		 $ Management_Zone: chr [1:318] "FTC" "FTC" "FTC" "FTC" ...
		 $ Sampling_Depth : chr [1:318] "0-30" "30-60" "60-90" "0-30" ...
		 $ Variables      : chr [1:318] "WTC (%)" "WTC (%)" "WTC (%)" "WTC (%)" ...
		 $ Value          : num [1:318] 352 723 1084 417 558 ...
		 $ ID             : Factor w/ 12 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
		 - attr(*, "na.action")= 'omit' Named int [1:6] 59 67 86 94 126 310
		  ..- attr(*, "names")= chr [1:6] "59" "67" "86" "94" ...


# bikin list label

IndvarNames = list(
'1'="PWC (%)",
'2'="pH",
'3'="Org_C (%)",
'4'="PAC (%)",
'5'="CNR",
'6'="TN (%)",
'7'="TP (ppm)",
'8'="TK (ppm)",
'9'="TFe (ppm)",
'10'="TMn (ppm)",
'11'="TCu(ppm)",
'12'="TZn (ppm)")

#bikin fungsi label, passing ke facet_wrap ggplot2

indvar_labeller <- function(variable,value){
  return(IndvarNames[value])
}


# fungsi eliminasi outlier, passing ke stat_summary ggplot2

calc_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(stats)
}

windowsFonts(Times=windowsFont("Times New Roman"))

## plot faktor lingkungan 
indvar_p1 <- ggplot(R_boxplot1, aes(x=Management_Zone, y=Value, fill = Management_Zone)) + 
	#geom_violin(trim=FALSE) + 
    stat_summary(fun.data = calc_stat, geom="boxplot" ,alpha = 0.1 ) + 	
	geom_jitter( alpha = 0.6, aes(color = Management_Zone))+
    facet_wrap(ID~., labeller=indvar_labeller, scales="free", ncol=4) +
    #scale_fill_viridis(discrete = TRUE,  option = "turbo", alpha=0.7)+
	#scale_fill_brewer(palette="BuPu")+
	#scale_fill_brewer(palette="OrRd")+
	theme_bw()+ 
	theme(text=element_text(family="Times")) +
    theme(legend.position = "none")+
	theme(panel.grid = element_blank())+
	theme(
				strip.text.x = element_text(size = 12, face='bold'),
				axis.text.y=element_text(family="Times", size=12),
				axis.text.x=element_text(family="Times", size=12),
				axis.title.y=element_blank(),
				axis.title.x=element_text(family="Times", size=13))+
	labs(x = "Management Zone")
	

## plot faktor lingkungan 
indvar_p2 <- ggplot(R_boxplot1, aes(x=Sampling_Depth, y=Value, fill = Sampling_Depth)) + 
	#geom_violin(trim=FALSE) + 
    stat_summary(fun.data = calc_stat, geom="boxplot" ,alpha = 0.1 ) + 	
	geom_jitter( alpha = 0.6, aes(color = Sampling_Depth))+
    facet_wrap(ID~., labeller=indvar_labeller, scales="free", ncol=4) +
    #scale_fill_viridis(discrete = TRUE,  option = "turbo", alpha=0.7)+
	#scale_fill_brewer(palette="BuPu")+
	#scale_fill_brewer(palette="OrRd")+
	theme_bw()+ 
	theme(text=element_text(family="Times")) +
    theme(legend.position = "none")+
	theme(panel.grid = element_blank())+
	theme(
				strip.text.x = element_text(size = 12, face='bold'),
				axis.text.y=element_text(family="Times", size=12),
				axis.text.x=element_text(family="Times", size=12),
				axis.title.y=element_blank(),
				axis.title.x=element_text(family="Times", size=13))+
	labs(x = "Sampling Depth (cm)")


indvar_p1
indvar_p2


library(readxl)
R_boxplot2 <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "R_boxplot2")
View(R_boxplot2)


R_boxplot2 <- na.omit(R_boxplot2)

R_boxplot2$ID <- factor(R_boxplot2$ID,levels=unique(R_boxplot2$ID))

str(R_boxplot2)
		tibble [185 × 5] (S3: tbl_df/tbl/data.frame)
		 $ Management_Zone: chr [1:185] "FTC" "FTC" "FTC" "FTC" ...
		 $ Sampling_Depth : chr [1:185] "0-30" "30-60" "60-90" "0-30" ...
		 $ Variables      : chr [1:185] "PBL (CFU/ml)" "PBL (CFU/ml)" "PBL (CFU/ml)" "PBL (CFU/ml)" ...
		 $ Value          : num [1:185] 6600000 9650000 14750000 4750000 15250000 ...
		 $ ID             : Factor w/ 8 levels "13","14","15",..: 1 1 1 1 1 1 1 1 1 1 ...
		 - attr(*, "na.action")= 'omit' Named int [1:31] 4 5 6 14 16 17 20 22 25 39 ...
		  ..- attr(*, "names")= chr [1:31] "4" "5" "6" "14" ...
	
	
DepvarNames = list(
'13'="PBL (CFU/ml)",
'14'="PBS (CFU/ml)",
'15'="PFS (CFU/ml)",
'16'="LKS (U/ml, min)",
'17'="Mn_P (U/ml, min)",
'18'="Li_P (U/ml, min)",
'19'="EBS (U/ml, min)",
'20'="EFS (U/ml, min)"
)


depvar_labeller <- function(variable,value){
  return(DepvarNames[value])
}


depvar_p1 <- ggplot(R_boxplot2, aes(x=Management_Zone, y=Value, fill = Management_Zone)) + 
	#geom_violin(trim=FALSE) + 
    stat_summary(fun.data = calc_stat, geom="boxplot" ,alpha = 0.1 ) + 	
	geom_jitter( alpha = 0.6, aes(color = Management_Zone))+
    facet_wrap(ID~., labeller=depvar_labeller, scales="free", ncol=4) +
    #scale_fill_viridis(discrete = TRUE,  option = "turbo", alpha=0.7)+
	#scale_fill_brewer(palette="BuPu")+
	#scale_fill_brewer(palette="OrRd")+
	theme_bw()+ 
	theme(text=element_text(family="Times")) +
    theme(legend.position = "none")+
	theme(panel.grid = element_blank())+
	scale_y_continuous(limits = c(0, NA))+ theme(
				strip.text.x = element_text(size = 12, face='bold'),
				axis.text.y=element_text(family="Times", size=12),
				axis.text.x=element_text(family="Times", size=12),
				axis.title.y=element_blank(),
				axis.title.x=element_text(family="Times", size=13))+
	labs(x = "Management Zone")

depvar_p2 <- ggplot(R_boxplot2, aes(x=Sampling_Depth, y=Value, fill = Sampling_Depth)) + 
	#geom_violin(trim=FALSE) + 
    stat_summary(fun.data = calc_stat, geom="boxplot" ,alpha = 0.1 ) + 	
	geom_jitter( alpha = 0.6, aes(color = Sampling_Depth))+
    facet_wrap(ID~., labeller=depvar_labeller, scales="free", ncol=4) +
    #scale_fill_viridis(discrete = TRUE,  option = "turbo", alpha=0.7)+
	#scale_fill_brewer(palette="BuPu")+
	#scale_fill_brewer(palette="OrRd")+
	theme_bw()+ 
	theme(text=element_text(family="Times")) +
    theme(legend.position = "none")+
	theme(panel.grid = element_blank())+
	scale_y_continuous(limits = c(0, NA))+ theme(
				strip.text.x = element_text(size = 12, face='bold'),
				axis.text.y=element_text(family="Times", size=12),
				axis.text.x=element_text(family="Times", size=12),
				axis.title.y=element_blank(),
				axis.title.x=element_text(family="Times", size=13))+
	labs(x = "Sampling Depth (cm)")

depvar_p2+depvar_p1+plot_layout(ncol=1)



Rekap_no_outlier <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "Rekap_R_1")
View(Rekap_no_outlier)

Rekap_no_outlier <- na.omit(Rekap_no_outlier)
					
Rekap_no_outlier$Management_Zone <- as.factor(Rekap_no_outlier$Management_Zone)
Rekap_no_outlier$Sampling_Depth <- as.factor(Rekap_no_outlier$Sampling_Depth)	

#PCA 
Data_hadi.pca <- PCA(Rekap_no_outlier[,-(1:2)], graph = FALSE)

eig.val <- get_eigenvalue(Data_hadi.pca )
eig.val  


res.pca$Data_hadi.pca #2



head(Data_hadi.pca$var$coord)


windowsFonts(Times=windowsFont("Times New Roman"))

theme1_pca <- theme(title =element_blank(),
				axis.text.y=element_text(family="Times", size=14),
				axis.text.x=element_text(family="Times", size=14),
				axis.title.y=element_text(family="Times", size=14),
				axis.title.x=element_blank())

theme2_pca <- theme(title =element_blank(),
				axis.text.y=element_text(family="Times", size=14),
				axis.text.x=element_text(family="Times", size=14),
				axis.title.y=element_blank(),
				axis.title.x=element_text(family="Times", size=14))

theme3_pca <- theme(title =element_blank(),
				axis.text.y=element_text(family="Times", size=14),
				axis.text.x=element_text(family="Times", size=14),
				axis.title.y=element_blank(),
				axis.title.x=element_blank())

bar_dim <- fviz_eig(Data_hadi.pca, ncp = 10, addlabels = TRUE, barfill="grey", barcolor="grey") #3
bar_dim1 <- bar_dim + theme1_pca
ctrb_dim1 <- fviz_contrib(Data_hadi.pca,choice='var',top=10, fill="grey", color="grey", title = NULL)
dim1<-ctrb_dim1+theme1_pca
ctrb_dim2 <- fviz_contrib(Data_hadi.pca,choice='var',top=10,axes=2, fill="grey", color="grey", title = NULL)
dim2<-ctrb_dim2+theme3_pca
ctrb_dim3 <- fviz_contrib(Data_hadi.pca,choice='var',top=10,axes=3, fill="grey", color="grey", title = NULL)
dim3<-ctrb_dim3+theme3_pca
ctrb_dim4 <- fviz_contrib(Data_hadi.pca,choice='var',top=10,axes=4, fill="grey", color="grey", title = NULL)
dim4<-ctrb_dim4+theme3_pca
ctrb_dim5 <- fviz_contrib(Data_hadi.pca,choice='var',top=10,axes=5, fill="grey", color="grey", title = NULL)
dim5<-ctrb_dim5+theme3_pca

dim1+dim2+dim3 + dim4+ plot_layout(ncol=4)



pcavar <- fviz_pca_var(Data_hadi.pca,
             col.var = "contrib", # Color by contributions to the PC
            # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )+ ggtitle ("")
			 
plot_pcavar <- pcavar+
								theme(text=element_text(family="Times"))

plot_pcavar_Fix <- ggpubr::ggpar(plot_pcavar,
                                  xlab = "PC1 (39.4%)", ylab = "PC2 (27.1%)")
bar_dim1+plot_pcavar_Fix

# grouping based on Kedalaman
pcaKD <- fviz_pca_ind(Data_hadi.pca , #11
                      geom.ind = "point", # show points only (nbut not "text")
                      col.ind = Rekap_no_outlier$Sampling_Depth       , # color by groups
                      #palette = viridis,
                      addEllipses = TRUE,
						ellipse.level = 0.95, 
						ellipse.type = c("confidence"),
                      legend.title = "Sampling Depth (cm)")+ ggtitle ("PCA - Sampling Depth")
plot_pcaKD <- pcaKD+
								theme(text=element_text(family="Times"))+  theme(legend.position="bottom")

plot_pcaKD_Fix <- ggpubr::ggpar(plot_pcaKD,
                                  xlab = "PC1 (34.8%)", ylab = "PC2 (19.7%)")
								  
# grouping based on Management Zone
pcaZP <- fviz_pca_ind(Data_hadi.pca , #11
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Rekap_no_outlier$Management_Zone, # color by groups
                     #    palette = c("seagreen4", "red4", "plum4","slategray3"),
                         addEllipses = TRUE,
						ellipse.level = 0.95, 
						ellipse.type = c("confidence"), # Concentration ellipses
                         legend.title = "Management Zone") + ggtitle ("PCA - Management Zone")



plot_pcaZP <- pcaZP+
								theme(text=element_text(family="Times"))+  theme(legend.position="bottom")

plot_pcaZP_Fix <- ggpubr::ggpar(plot_pcaZP,
                                   xlab = "PC1 (34.8%)", ylab = "")

plot_pcavar_Fix+plot_pcaKD_Fix+plot_pcaZP_Fix

					
# Supplying Fuzzy Cmeans to PCA	
Data_hadi.pca <- PCA(Rekap_no_outlier[,-(1:2)], graph = FALSE)

All_clust2 <- fanny(x = Rekap_no_outlier[,-(1:2)], k = 2, metric = "SqEuclidean")

All_clust2$membership %>% head()

		
All_clust2$coeff

		
All_clust2$clustering %>% unique()


All_clust3 <- fanny(x = Rekap_no_outlier[,-(1:2)], k = 3, metric = "SqEuclidean")
All_clust3$membership %>% head()

All_clust3$coeff

All_clust3$clustering %>% unique()



fuzzy_all <- Rekap_no_outlier 
fuzzy_all$cluster2 <- All_clust2$clustering
fuzzy_all$cluster3 <- All_clust3$clustering

fuzzy_all$cluster2 <- as.factor(fuzzy_all$cluster2)
fuzzy_all$cluster3 <- as.factor(fuzzy_all$cluster3)
								
biplot_All_fuzzy_cl2 <- fviz_pca_ind(Data_hadi.pca, 
								geom.ind = "point", 
                                col.ind = fuzzy_all$cluster2, 
                             #   palette = "NULL",
                              #  addEllipses = FALSE,
                               # label = "var",
                               # col.var = "black", 
                               # repel = TRUE,
							   addEllipses = TRUE,
						ellipse.level = 0.95, 
						ellipse.type = c("confidence"),
                                legend.title = "Cluster")+ ggtitle ("PCA-Fuzzy C-means Cluster=2")+
								theme(text=element_text(family="Times"))

biplot_Fuzzy_All_fix_cl2 <- ggpubr::ggpar(biplot_All_fuzzy_cl2,
                                 xlab = "PC1 (34.8%)", ylab = "")+  theme(legend.position="bottom")

biplot_All_fuzzy_cl3 <- fviz_pca_ind(Data_hadi.pca, 
								geom.ind = "point", 
                                col.ind = fuzzy_all$cluster3, 
                                #palette = "NULL",
                               # addEllipses = FALSE,
                               # label = "var",
                               # col.var = "black", 
                               # repel = TRUE,
							   addEllipses = TRUE,
						ellipse.level = 0.95, 
						ellipse.type = c("confidence"),
                                legend.title = "Cluster")+ ggtitle ("PCA-Fuzzy C-means Cluster=3")+
								theme(text=element_text(family="Times"))
	
biplot_Fuzzy_All_fix_cl3 <- ggpubr::ggpar(biplot_All_fuzzy_cl3,
                                   xlab = "PC1 (34.8%)", ylab = "")+  theme(legend.position="bottom")

plot_pcavar_Fix+bar_dim1

plot_pcaKD_Fix + plot_pcaZP_Fix+ biplot_Fuzzy_All_fix_cl2+biplot_Fuzzy_All_fix_cl3+plot_layout(ncol=4)

# Visualize kmeans clustering
# use repel = TRUE to avoid overplotting
PlotKmeans <- fviz_cluster(km.res, Data_Penelitian[, -(1:2)], 
				addEllipses = FALSE, #Concentration ellipses
				#ellipse.type = "norm")
			ggtheme = theme_minimal()
					)


PlotKmeans1 <- PlotKmeans+
								theme(text=element_text(family="Times"))

PlotKmeans_Fix <- ggpubr::ggpar(PlotKmeans1,
                                  xlab = "", ylab = "")

plot_pcavar_Fix+plot_pcaKD_Fix+plot_pcaZP_Fix+PlotKmeans_Fix+plot_layout(ncol=4)






R_PCA_en3 <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "R_PCA_en3")
View(R_PCA_en3)

R_PCA_en3 <- na.omit(R_PCA_en3)						
R_PCA_en3$Management_Zone <- as.factor(R_PCA_en3$Management_Zone)
R_PCA_en3$Sampling_Depth <- as.factor(R_PCA_en3$Sampling_Depth)			

str(R_PCA_en3)
		tibble [12 × 5] (S3: tbl_df/tbl/data.frame)
		 $ Management_Zone: Factor w/ 3 levels "FRS","FTC","HVP": 2 2 2 2 2 1 1 1 3 3 ...
		 $ Sampling_Depth : Factor w/ 3 levels "0-30","30-60",..: 1 2 3 1 2 1 2 3 3 2 ...
		 $ PBL            : num [1:12] 6600000 9650000 14750000 4750000 15250000 ...
		 $ PBS            : num [1:12] 8700000 12150000 10350000 650000 200000 ...
		 $ PFS            : num [1:12] 6818 4950 4265 3333 4000 ...
		 - attr(*, "na.action")= 'omit' Named int [1:15] 4 5 6 9 12 13 14 16 17 18 ...
		  ..- attr(*, "names")= chr [1:15] "4" "5" "6" "9" ...
		

soil.pca3 <- PCA(R_PCA_en3[,-(1:2)], graph = FALSE)

soil.pca3$eig #2

bar_dim_soil3 <- fviz_eig(soil.pca3, ncp = 10, addlabels = TRUE, barfill="grey", barcolor="grey") #3
bar_dim_soil13 <- bar_dim_soil3 + theme1_pca
windowsFonts(Times=windowsFont("Times New Roman"))
soil.pcavar3 <- fviz_pca_var(soil.pca3,
             col.var = "contrib", # Color by contributions to the PC
            # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
			 
plotsoil_pcavar3 <- soil.pcavar3+
								theme(text=element_text(family="Times"))



plotsoil_pcavar_Fix3 <- ggpubr::ggpar(plotsoil_pcavar3,
                                  xlab = "PC1 (37.1%)", ylab = "PC2 (25.8%)")



biplot_soil_Zonasi3 <- fviz_pca_biplot(soil.pca3, 
                                col.ind = R_PCA_en3$Management_Zone, 
                                palette = "NULL",
                                addEllipses = TRUE,
								ellipse.level = 0.95, 
								ellipse.type = c("confidence"),
                                label = "var",
                                col.var = "black", 
                                repel = TRUE,
                                legend.title = "Management Zone")+ ggtitle ("") + 
								theme(text=element_text(family="Times"))

biplot_soil_Zonasi_fix3 <- ggpubr::ggpar(biplot_soil_Zonasi3,
                                 xlab = "PC1 (37.1%)", ylab = "")+ 
								 theme(text=element_text(family="Times"))
biplot_soil_Zonasi_fix3 

biplot_soil_Kdl3 <- fviz_pca_biplot(soil.pca3, 
                                col.ind = R_PCA_en3$Sampling_Depth, 
                                palette = "NULL",
                                addEllipses = TRUE,
								ellipse.level = 0.95, 
								ellipse.type = c("confidence"),
                                label = "var",
                                col.var = "black", 
                                repel = TRUE,
                                legend.title = "Sampling Depth (cm)")+ ggtitle ("")

biplot_soil_Kdl_fix3 <- ggpubr::ggpar(biplot_soil_Kdl3,
                                 xlab = "PC1 (37.1%)", ylab = "PC2 (25.8%)")+ 
								 theme(text=element_text(family="Times"))
biplot_soil_Kdl_fix3

biplot_soil_Kdl_fix3+  theme(legend.position="bottom")+biplot_soil_Zonasi_fix3 + theme(legend.position="bottom")





R_PCA_en4 <- read_excel("D:/Publikasi/Pak Heru B Pulunggono/20 microbial quantities and ensime activity in tropical peat/Data Penelitian.xlsx", 
    sheet = "R_PCA_en4")
View(R_PCA_en4)

R_PCA_en4 <- na.omit(R_PCA_en4)						
R_PCA_en4$Management_Zone <- as.factor(R_PCA_en4$Management_Zone)
R_PCA_en4$Sampling_Depth <- as.factor(R_PCA_en4$Sampling_Depth)			

str(R_PCA_en4)
		tibble [17 × 7] (S3: tbl_df/tbl/data.frame)
		 $ Management_Zone: Factor w/ 3 levels "FRS","FTC","HVP": 2 2 2 2 2 2 1 1 1 1 ...
		 $ Sampling_Depth : Factor w/ 3 levels "0-30","30-60",..: 1 3 1 3 1 2 1 1 2 3 ...
		 $ LKS            : num [1:17] 208.3 254.6 196.8 23.1 34.7 ...
		 $ Mn_P           : num [1:17] 20.66 34.44 41.32 6.89 34.44 ...
		 $ Li_P           : num [1:17] 1245.5 80.6 152.3 1066.3 1173.8 ...
		 $ EBS            : num [1:17] 50.741 3.056 0.648 14.167 1.574 ...
		 $ EFS            : num [1:17] 9.63 5.93 9.44 6.48 3.52 ...
		 - attr(*, "na.action")= 'omit' Named int [1:10] 2 5 9 11 12 16 17 19 20 21
		  ..- attr(*, "names")= chr [1:10] "2" "5" "9" "11" ...
		

soil.pca4 <- PCA(R_PCA_en4[,-(1:2)], graph = FALSE)

soil.pca4$eig #2

bar_dim_soil4 <- fviz_eig(soil.pca4, ncp = 10, addlabels = TRUE, barfill="grey", barcolor="grey") #3
bar_dim_soil14 <- bar_dim_soil4 + theme1_pca
windowsFonts(Times=windowsFont("Times New Roman"))
soil.pcavar4 <- fviz_pca_var(soil.pca4,
             col.var = "contrib", # Color by contributions to the PC
            # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
			 
plotsoil_pcavar4 <- soil.pcavar4+
								theme(text=element_text(family="Times"))



plotsoil_pcavar_Fix4 <- ggpubr::ggpar(plotsoil_pcavar4,
                                  xlab = "PC1 (37.1%)", ylab = "PC2 (25.8%)")



biplot_soil_Zonasi4 <- fviz_pca_biplot(soil.pca4, 
                                col.ind = R_PCA_en4$Management_Zone, 
                                palette = "NULL",
                                addEllipses = TRUE,
								ellipse.level = 0.95, 
								ellipse.type = c("confidence"),
                                label = "var",
                                col.var = "black", 
                                repel = TRUE,
                                legend.title = "Management Zone")+ ggtitle ("") + 
								theme(text=element_text(family="Times"))

biplot_soil_Zonasi_fix4 <- ggpubr::ggpar(biplot_soil_Zonasi4,
                                 xlab = "PC1 (37.1%)", ylab = "")+ 
								 theme(text=element_text(family="Times"))
biplot_soil_Zonasi_fix4 

biplot_soil_Kdl4 <- fviz_pca_biplot(soil.pca4, 
                                col.ind = R_PCA_en4$Sampling_Depth, 
                                palette = "NULL",
                                addEllipses = TRUE,
								ellipse.level = 0.95, 
								ellipse.type = c("confidence"),
                                label = "var",
                                col.var = "black", 
                                repel = TRUE,
                                legend.title = "Sampling Depth (cm)")+ ggtitle ("")
								
biplot_soil_Kdl_fix4 <- ggpubr::ggpar(biplot_soil_Kdl4,
                                 xlab = "PC1 (37.1%)", ylab = "PC2 (25.8%)")+ 
								 theme(text=element_text(family="Times"))
biplot_soil_Kdl_fix4

biplot_soil_Kdl_fix4+  theme(legend.position="bottom")+biplot_soil_Zonasi_fix4 + theme(legend.position="bottom")

biplot_soil_Kdl_fix3+  theme(legend.position="bottom")+biplot_soil_Zonasi_fix3 + theme(legend.position="bottom")+
biplot_soil_Kdl_fix4+  theme(legend.position="bottom")+biplot_soil_Zonasi_fix4 + theme(legend.position="bottom")+
plot_layout(ncol=4)
