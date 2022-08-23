setwd("C:/Users/Vicent/Documents/Uni/Article Iris/Estadistic") 
#Aixo serveix per dir-li on ha de buscar el arxiu amb les dades. MoTissue important que els
#separadors "/" estiguin amb aquesta orientació
getwd()
list.files()
install.packages("readxl") 
library(readxl)

mydata <- read_excel("C:/Users/Vicent/Documents/Uni/Article Iris/Estadistic/Arxius R.xlsx", 
                     sheet = "ABA", col_names=T)
View(mydata)


factors <- c("Tissue", "Time")
mydata[factors] <- lapply(mydata[factors], factor)

newdata <- mydata[ which((mydata$Tissue=="FM" | mydata$Tissue=="FS")), ]
View(newdata)


library("ggpubr")
ggdensity(mydata$Fvfm, 
          main = "Density plot of RWC",
          xlab = "RWC")
library(ggpubr)
ggqqplot(mydata$Z)

shapiro.test(newdata$NxVx)

#Correlation tests

#Pearson si les variables segueixen una distribució normal
res <- cor.test(newdata$NxVx, newdata$ChTissueot, method = "pearson")
res

#Spearman si no segueixen una distribució normal
res2 <-cor.test(newdata$ABA, newdata$CK,  method = "spearman")
res2

library("ggpubr")
ggscatter(newdata, x = "C", y = "IPA", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ABA", ylab = "CK")

str(mydata) #Comprovar factors de les meves dades

library(lattice) 
library(plyr) 
library(dplyr)
library(Rmisc)
summarySE(newdata, measurevar="ABA", 
          groupvars=c("Tissue", "Time"))

#Per comprovar homocedasticitat de variancies
bartlett.test(DPS ~ interaction(Tissue:Time), data=newdata)

#One-way ANOVA
plant.lm <- lm(Nx ~ Tissue + Time + Tissue:Time, data=newdata)
summary(plant.lm)
plant.av <- aov(plant.lm)
summary(plant.av)

#Fer ANOVA amb muTissueiples factors:
plant.lm <- lm(DPS ~ Tissue + Time + Tissue:Time
               , data = newdata)
summary(plant.lm)
plant.av <- aov(plant.lm)
summary(plant.av)

#Per checkear normalitat de residus!
shapiro.test(residuals(plant.av))

#post-hoc
newdata2 <- filter(newdata, Time == "Dec")
View (newdata2)

plant.lm <- lm(DPS ~ Tissue
               , data = newdata2)
summary(plant.lm)
plant.av <- aov(plant.lm)
summary(plant.av)
tukey.test <- TukeyHSD(plant.av)
tukey.test

#Duncan test
library(agricolae)
duncan.test(plant.av, "Tissue",console=TRUE)


#Per fer ART
library(ARTool)
m = art(ABA ~ Tissue + Time + Tissue:Time, data = newdata)
anova(m)
install.packages(lsmeans)
library(lsmeans)

newdata2 <- filter(newdata, Time == "Dec")
View (newdata2)
install.packages(lsmeans)
library(lsmeans)
m = art(ABA ~ Tissue, data = newdata2)
anova(m)
lsmeans(artlm(m, "Tissue"), pairwise ~ Tissue)

#Kruskal Wallis
newdata2 <- filter(newdata, Time == "Nov")
newdata <- mydata[ which((mydata$Tissue=="FV" | mydata$Tissue=="FM" | mydata$Tissue=="FS") & (mydata$Time=="Nov")), ]
View(newdata2)
kruskal.test(Z ~ Tissue, data = newdata2)

#post-hoc KW
install.packages("FSA")
library(FSA)

dunnTest(Z ~ Tissue,
         data=newdata2,
         method="bonferroni")

#T-test
# Compute t-test
res <- t.test(H ~ Tissue, data = newdata2, var.equal = TRUE)
res

summarySE(newdata, measurevar="ABA", 
          groupvars=c("Tissue", "Time"))
newdata <- mydata[ which((mydata$Tissue=="R" | mydata$Tissue=="A")), ]
View(newdata)

library(ggplot2)
# Basic box plot
ggplot(newdata, aes(x=Time, y=Z, fill=Tissue)) +
  geom_boxplot()
# Rotate the box plot
p + coord_flip()
# Notched box plot
ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot(notch=TRUE)
# Change outlier, color, shape and size
ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)

