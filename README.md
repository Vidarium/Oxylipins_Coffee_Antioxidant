# Oxylipins Coffee Antioxidant

This repository includes a code and the database for statistical analyses of oxylipins in the foam cell model as supplementary material for the manuscript:  

Oscar J. Lara-Guzmán, Sonia Medina, Rafael Álvarez, Camille Oger, Thierry Durand, Jean-Marie Galano, Natalia Zuluaga, Ángel Gil-Izquierdo and Katalina Muñoz-Durango. 2020. Oxylipin regulation by phenolic compounds from coffee beverage:  positive outcomes from a randomized controlled trial in healthy adults and macrophage derived foam cells. 
The manuscript is currently submitted to Free Radical Biology & Medicine journal as a research article.

The file Oxylipins_CGAs_FoamCells.R is a commented R script that requires the provided database file Oxylipins_CGAs_FoamCells.csv to properly run.  

RStudio and R software used: R i386 3.6.1


### Summary step by step:
1.	Select the database and set up:
-	Encoding: automatic
-	Heading: YES
-	Row names: Use first column
-	Separator: Semicolon
-	Decimal: Period
-	Quote: Double quote (")
-	na.strings: NA
2.	Load the required libraries: factoextra, utils, stats, FactoMineR, ggfortify, ggplot2, "corrplot", magrittr, cluster, ggpubr, NbClust, REdaS, DiscriMiner, corrr, tidyverse, ggraph, devtools, dplyr, igraph.
3.	Attach the database.
4.	Two data frames derived from the database are required for analysis:
-	All_treatments, which contain all data, except categorical variables.
-	Treatments, which contain all data, except observations for MIX treatment and positive control. This database is used for the correlation analysis.
5.	Scale and normalize the databases All_treatments and Treatments to:
-	All_Treatments_Scale.  For KMO (Kaiser-Meyer-Olkin) test and PCAs.
-	Treatments_Scale. For clustering and distance analysis.
6.	Create correlations from treatments using the command “cor”.
-	Use the command “corrplot” to generate graphics.
-	Create the function cor.mtest and matrix of the p-value of the correlation. Execute the command “corrplot” for visualization.
7.	Create a heatmap of distances measure or (dis)similarity between observations using the commands “get_dist” and  “fviz_dist”.
8.	Estimate the KMO value (Kaiser-Meyer-Olkin), this test measures how suited your dataset is for Factor Analysis.
9.	For the PCA analysis is required a cluster analysis. “fviz_nbclust” and “fviz_cluster” functions are used to estimate the number of the clusters and visualization.
10.	Calculate the eigenvalues with “get_eigenvalue” command and visualize this with “fviz_eig”.
11.	Use “fviz_pca_ind” and “fviz_pca_var” to visualize the PCA score and loading plots.
12.	Use  “dimdesc”  to estimate the correlations and P-value for each variable in PC1 and PC2.

