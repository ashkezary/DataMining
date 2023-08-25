library(dplyr)
library(readr)
colon <- read_table2("colon-cancer", col_names = c('cancer',paste0("Gene",1:2000)))

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


rank_func <- function(d){
  scores <- data.frame()
  for (x in colnames(d)[-1]) {
    pos <- d %>% filter(cancer==1) %>%select(x)
    neg <- d %>% filter(cancer==-1) %>%select(x)
    pos <- pos[[1]]
    neg <- neg[[1]]
    pos_mean = mean(pos)
    neg_mean = mean(neg)
    pos_var = var(pos)
    neg_var = var(neg)
    s = abs(pos_mean - neg_mean)/sqrt(pos_var/length(pos) + neg_var/length(neg))
    scores <- rbind(scores,c(x,s))
  }
  return(scores)
}

att_score1 <- rank_func(colon)
colnames(att_score1) <- c("Gene","Score")
tops1 <- att_score1[order(-as.numeric(att_score1$Score)),]

entropy <- function(d){
  m <- dim(d)[2]
  for (x in 2:m) {
    d[,x] = d[,x][[1]]/(max(d[,x][[1]])-min(d[,x][[1]]))
  }
  distances <- as.matrix(dist(t(d[,2:2001]),upper = T,diag = T),nrow = 2000)
  s <- exp(-0.5*distances)
  ss <- s[lower.tri(s)]
  sss <- ss[ss<1]
  total = - sum(sss*log2(sss)) - sum((1-sss)*log2(1-sss))
  return (total)
}  

entirent <- entropy(colon)

IG <- function(d){
  n = dim(d)[2]
  ranks = data.frame()
  for (v in 1:n) {
    e = entropy(d)
    minimum = e
    k = -1
    for (h in 2:2001) {
      dif = e-entropy(d[,-h])
      if(dif<minimum){
        minimum = dif
        k = h
      }
    }
    ranks <- rbind(k,dif)
  }
  return(ranks)
}
    
att_score2 <- IG(colon[,-1])

library(e1071)
ind <- sample(1:62,size = 50)
rectu <- colon[,c('cancer',tops1[1:10,]$Gene)]
train <- rect[ind,]
test <- rect[-ind,]
SVMclassifier <- svm(formula = cancer ~ . ,
                      data = train,
                     type = 'C-classification',
                     kernel = 'linear')
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
cat(paste0("        |   Real+  |   Real-   ","\n",
           "predict+|   ",tp,"   |   ", fp,"   ","\n",
           "predict-|   ",fn,"   |   ", tn,"   ","\n",
           "******************************************","\n",
           "Accuracy = ", acc,"\n",
           "Sensitivity = ", sens,"\n",
           "Specificity = ", spec,"\n",
           "f-score = ", f.score,"\n",
           "******************************************"))