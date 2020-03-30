#DATABASE NAME: Oxylipins_CGAs_FoamCells 
  
#LIBRARIES: (execute all, twice)
    library(factoextra) 
    library(utils) 
    library(stats) 
    library(FactoMineR) 
    library(ggfortify) 
    library(ggplot2) 
    library("corrplot") 
    library(magrittr) 
    library(cluster) 
    library(ggpubr) 
    library("NbClust") 
    library(REdaS)
    library(DiscriMiner) 
    library(corrr)
    library(tidyverse)
    library(ggraph)
    library(devtools)
    library(dplyr)
    library(igraph)

# 1. DATABASE

    attach(Oxylipins_CGAs_FoamCells)

# 2. DATA FRAMES GENERATED FROM THE ORIGINAL DATABASE 
     
    # 1. Database for analysis of all treatments vs positive control
    All_Treatments <- Oxylipins_CGAs_FoamCells [c(-1)] # The categoric variable is exclude. 
    write.table(All_Treatments,file = "All_Treatments.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
    
    # 2. Database for analysis for correlation, distance measures, and cluster ANALYSIS. Without the MIX treatment and positive control
    Treatments <- Oxylipins_CGAs_FoamCells [c(-1)] [-(66:70),] [-(1:5),] #Exclusion of the MIX treatment.
    write.table(Treatments,file = "Treatments.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
    
    #Scaled and normalized data 
    
    All_Treatments_Scale = scale(All_Treatments) #For PCA
    write.table(All_Treatments_Scale,file = "All_Treatments_Scale.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
    
    Treatments_Scale = scale(Treatments) #For measure distance analysis
    write.table(Treatments_Scale,file = "Treatments_Scale.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
    

# 3. CORRELATIONS BETWEEN VARIABLES 
      
    #Correlation between all oxylipins produced in foam cells after treatments with phenolic compounds and CGAs
      correlations <- (cor(Treatments)) 
      correlations
      corrplot(correlations, is.corr=FALSE, method="circle", type = "upper", tl.cex = 0.4, tl.srt = 90)
      write.table(correlations,file = "correlations.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
      
      #P_value
      cor.mtest <- function(Treatments, ...) {
        Treatments <- as.matrix(Treatments)
        n <- ncol(Treatments)
        p.Treatments<- matrix(NA, n, n)
        diag(p.Treatments) <- 0
        for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
            tmp <- cor.test(Treatments[, i], Treatments[, j], ...)
            p.Treatments[i, j] <- p.Treatments[j, i] <- tmp$p.value
          }
        }
        colnames(p.Treatments) <- rownames(p.Treatments) <- colnames(Treatments)
        p.Treatments
      }
      
      # matrix of the p-value of the correlation
      p.Treatments <- cor.mtest(Treatments)

      #Corrplot
      corrplot(correlations, is.corr=FALSE, method="circle", type = "upper", p.mat = p.Treatments, sig.level = 0.05, tl.cex = 0.6, tl.srt = 90, order="hclust")

# 4. DISTANCE MEASURE OR (DIS)SIMILARITY BETWEEN OBSERVATIONS
     
      HM <- get_dist(Treatments_Scale, stand = TRUE, method = "pearson")
      fviz_dist(HM, gradient = list(low = "#02364f", mid = "white", high = "#b86005"), order = FALSE, lab_size = 10)    
      
     
# 5. KMO PROBE (Kaiser-Meyer-Olkin Test is a measure of how suited your data is for Factor Analysis)
      
      KMOS_VALIDACION =  KMOS(All_Treatments_Scale, use = "everything")
      write.table(KMOS_VALIDACION[["KMO"]],file = "KMOS_VALIDACION.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
      KMOS_VALIDACION

# 6.  Clusters analysis
      
      fviz_nbclust(All_Treatments, kmeans, method = "gap_stat")
      set.seed(123)
      kmeansGraph <- kmeans(All_Treatments_Scale, 4, nstart = 3)
      # Visualize
      fviz_cluster(kmeansGraph, data = All_Treatments_Scale, ellipse.type = "convex", palette = "jco",ggtheme = theme_minimal())
      
      
# 7. PCA
      
      #Eigenvalues
      res.pca <- PCA(All_Treatments_Scale, scale.unit = TRUE, ncp = 10, graph = TRUE)
      eig.val <- get_eigenvalue(res.pca)
      eig.val
      write.table(eig.val,file = "Autovalores.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
      fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 60))
      
      #PC1 vs PC2
      Scoreplot = fviz_pca_ind(res.pca, geom.ind = c("text","point"), fill.ind = Oxylipins_CGAs_FoamCells$Clusters, col.ind = "black", pointshape = 21, pointsize = "cos2", labelsize = 6, palette = c("#810569",  "#052c81", "#2a7405", "#c70e04"), addEllipses = TRUE, legend.title = list(size = "Quality"),  repel = TRUE, mean.point = FALSE, ellipse.level = 0.85, ellipse.type = "convex")
      ggpubr::ggpar(Scoreplot, title = "Oxylipins in Macrophage foam cells", subtitle = "dataset", caption = FALSE, xlab = "PC1(22.5%)", ylab = "PC2(54.4%)", legend.position = "top", legendsize = 8)
      
      LoadingPlot = fviz_pca_var(res.pca, geom.var = c("text","point"), pointshape = 20, pointsize = 10, labelsize = 6, addEllipses = FALSE, fill.var = "black", col.var = "cos2", gradient.cols = c("#d06a15", "#90460a", "#4d2605"), legend.title = list(color = "Quality"), repel = TRUE, mean.point = FALSE, ellipse.level = 0.85, ellipse.type = "norm", circle = FALSE)
      ggpubr::ggpar(LoadingPlot, title = "Oxylipins in Macrophage foam cells", subtitle = "dataset", caption = FALSE, xlab = "PC1(22.5%)", ylab = "PC2(54.4%)", legend.position = "top", legendsize = 8)
      
      Biplot = fviz_pca_biplot(res.pca, geom.var = c("text","point"), geom.ind = c("point"), fill.ind = Oxylipins_CGAs_FoamCells$Clusters, palette = c( "#810569",  "#052c81", "#2a7405", "#c70e04"), pointshape = 21, pointsize = 5, labelsize = 5, center = c(0, 1), fill.var = "black", col.var = "cos2", gradient.cols = c("#d06a15", "#90460a", "#4d2605"), legend.title = list(fill = "Groups", color = "QUALITY"), repel = TRUE, mean.point = TRUE, ellipse.type = "convex", addEllipses = FALSE)
      ggpubr::ggpar(Biplot, title = "Oxylipins in Macrophage foam cells", subtitle = "dataset", caption = FALSE, xlab = "PC1(22.5%)", ylab = "PC2(54.4%)", legend.position = "top", legendsize = 8)
    

# 8  Correlations and P_value for each variable in PC1 and PC2
        
        res.desc <- dimdesc(res.pca, axes = c(1,2,3), proba = 0.01)
        # Description of dimension 
        res.desc$Dim.1
        res.desc$Dim.2
        write.table(res.desc[["Dim.1"]][["quanti"]],file = "P_Value_PC1.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
        write.table(res.desc[["Dim.2"]][["quanti"]],file = "P_Value_PC2.xl", sep = "\t", eol = "\n", dec = ".", row.names = TRUE, col.names = TRUE) 
   
