
##Load Data
library(readr)
name <- c(paste0("Gene",1:2000))
colon <- read_table2("colon-cancer", col_names = c("cancer",name))

##Preprocessing
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




#val <- train[i*(1:5)]
#ktrain <- train[-i*(1:5)]
##train_test_split
library(randomForest)
set.seed(51)
ind <- sample(1:62,size = 50)
train <- colon[ind,]
test <- colon[-ind,]
#tree <- data.frame()
krf <- randomForest(cancer~. ,data = train,
                      xtest = test[-1] , ytest = test[[1]],
                      mtry=2000, importance =TRUE)  
#conf <- krf$confusion[1:2,1:2]
#error <- 1-sum(diag(tabl))/50
#tree <- rbind(tree,c(i,error))

tp <- krf$confusion[2,2]
tn <- krf$confusion[1,1]
fn <- krf$confusion[1,2]
fp <- krf$confusion[2,1]
acc <- (tp+tn)/(tp+tn+fp+fn)
sens <- tp/(tp+fn)
spec <- tn/(tn+fp)
f.score <- 2*(acc*sens)/(acc+sens)
cat(paste0("        |   Real+  |   Real-   ","\n",
           "predict+|   ",tp,"   |   ", fp,"   ","\n",
           "predict-|   ",fn,"   |   ", tn,"   ","\n",
           "******************************************","\n",
           "Accuracy = ", acc,"\n",
           "Sensitivity = ", sens,"\n",
           "Specificity = ", spec,"\n",
           "f-score = ", f.score,"\n",
           "******************************************"))