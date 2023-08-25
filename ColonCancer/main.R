#In the name of God
#Colon Cancer prediction
#### Load data
library(readr)
name <- c(paste0("Gene",1:2000))
colon <- read_table2("~/colon-cancer", col_names = name)

#### Preprocessing
for (i in 2:2001) {
  if(i<=10){
    colon[i] <- as.double(substr(colon[[i]],3,20)) 
  }else if(i<=100){
    colon[i] <- as.double(substr(colon[[i]],4,20))
  }else if(i<=1000){
    colon[i] <- as.double(substr(colon[[i]],5,20))
  }else{
    colon[i] <- as.double(substr(colon[[i]],6,20))
  }
}
colon$cancer <- as.factor(colon$cancer)

#### Boxplot for genes
pdf("results/boxplot.pdf")
boxplot(colon)
dev.off()

#### Heatmap for gene correlation
library(pheatmap)
pdf("results/CorHeatmap.pdf")
pheatmap(cor(colon[,2:2001]),width=15,height=15)
dev.off()
#### Heatmap for samples correaltion
pdf("results/CorHeatmap2.pdf")
pheatmap(cor(t(colon)[2:2001,]))
dev.off()
#### PCs of Genes
pc <- prcomp(colon[,2:2001])
pdf("results/PC.pdf")
plot(pc)
dev.off()

#### PCs of Samples
pct <- prcomp(t(colon[,2:2001]))
pdf("results/PCT.pdf")
plot(pct)
dev.off()
#### 1st & 2nd PCA
pdf("results/PCA.pdf")
plot(pc$x[,1:2])
plot(pct$x[,1:2])
dev.off()

#### 1st & 2nd PCR
pdf("results/PCR.pdf")
plot(pc$r[,1:2])
plot(pct$r[,1:2])
dev.off()
#### ggplot
library(ggplot2)
pcr <- data.frame(pc$x[,1:3], Group=colon$cancer)
pdf("results/PCA_samples.pdf")
ggplot(pcr, aes(PC1,PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()

library(limma)
design <- model.matrix(~cancer+0,colon)
colnames(design) <- c("neg","pos")
fit <- lmFit(t(colon[,2:2001]) ,design)
cont.matrix <- makeContrasts(pos-neg, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2,0.01)
tT <- topTable(fit2, adjust="fdr", sort.by = "B", number = Inf)
write.table(tT, "results/Colon.txt", row.names = F, sep = "\t", quote = F)
#markGenes <- subset(tT, logFC>1 & adj.P.Val<0.05)
#aml <- sub("///.*",""aml)
#library(caTools)
#set.seed(123)
#split <- sample.split(colon$cancer,SplitRatio=0.8)
#train <- subset(colon,split==1)
#test <- subset(colon, split==-1)

#train[,2:2001] = scale(train[,2:2001])
#test[,2:2001] = scale(train[,2:2001])

## Logistic Regression
ind <- sample(1:62,size = 50)
train <- rect[ind,]
test <- rect[-ind,]
classifier <- glm(formula = cancer ~ . ,
                  family = binomial,
                  data = train)
prob_pred <- predict(classifier, type = 'response', newdata = test[-1])
y_pred <- ifelse(prob_pred>0.5, 1, -1)
tp <- sum(y_pred==test[[1]] & y_pred==1)
tn <- sum(y_pred==test[[1]] & y_pred==-1)
fn <- sum(y_pred!=test[[1]] & y_pred==-1)
fp <- sum(y_pred!=test[[1]] & y_pred==1)
acc <- (tp+tn)/(tp+tn+fp+fn)
sens <- tp/(tp+fn)
spec <- tn/(tn+fp)
f.score <- 2*(acc*sens)/(acc+sens)
#table(test[[1]],y_pred)
print("        |   Real+  |   Real-   ")
print(paste0("predict+|   ",tp,"   |   ", fp,"   "))
print(paste0("predict-|   ",fn,"   |   ", tn,"   "))
print("******************************************")
print(paste0("Accuracy = ", acc))
print(paste0("Sensitivity = ", sens))
print(paste0("Specificity = ", spec))
print(paste0("f-score = ", f.score))
print("******************************************")

