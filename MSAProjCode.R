data=read.table("Adjusted_data.txt",sep=",",header=TRUE)

data$advance=0
data$G1
data<-data[which(data$G1!=0),]
data<-data[which(data$G2!=0),]
data<-data[which(data$G3!=0),]
data$advance<-100*((data$G2/data$G1+data$G3/data$G2)/2-1)

dim(data)
data<-data[,-2]
data

write.table(data,'Adjusted_data.txt',sep=',')

for (i in 2:(length(colnames(data))-4)) {
  data[,colnames(data)[i]]<-as.numeric(data[,colnames(data)[i]])
  data[,colnames(data)[i]]<-(data[,colnames(data)[i]]-min(data[,colnames(data)[i]]))/(max(data[,colnames(data)[i]])-min(data[,colnames(data)[i]]))
}

trainset<-data[data$school=='GP',-c(1,27,26)]
testset<-data[data$school=='MS',-c(1,27,26)]
trainset
trainset$G3<-(trainset$G3-min(trainset$G3))/(max(trainset$G3)-min(trainset$G3))
testset$G3<-(testset$G3-min(testset$G3))/(max(testset$G3)-min(testset$G3))
trainset.x<-trainset[,-c(3,10,17,26,25)]
trainset.y<-as.data.frame(trainset[,c(26,25)])
testset.x<-testset[,-c(3,10,17,26,25)]
testset.y<-testset[,c(26,25)]

allset.x<-rbind(trainset.x,testset.x)
allset.y<-rbind(trainset.y,testset.y)

colnames(trainset)
data
# fa
trainset.x
KMO(trainset.x)
cor(trainset.x)

library(psych)
cortest.bartlett(cor(trainset.x),dim(trainset.x)[1])

fap<-principal(trainset.x,9,rotate="varimax")
fap

matt<-fap$loadings
scores=t(t(matt)%*%solve(cor(trainset.x))%*%t(scale(trainset.x,scale = TRUE)))
testscores=t(t(matt)%*%solve(cor(trainset.x))%*%t(scale(testset.x,scale = TRUE)))
scores<-cbind(scores,as.data.frame(trainset.y))
testscores<-cbind(testscores,as.data.frame(testset.y))
matt[abs(matt)<0.35]<-0
matt
fap

# Reg
scores
fit<-lm(G3~  RC1 +RC3   +RC2  +RC4+  RC6+ RC5+RC9+RC7 +   RC8 ,data=scores)
summary(fit)
modelscore<-fit$coefficients[-1]
fittedscore<-as.matrix(scores[,-c(9,10)])%*%modelscore
plot(fittedscore[order(trainset.y[,2])],ylim = c(-0.5,0.5))
lines(trainset.y[order(trainset.y[,2]),2]-0.5,type = 'l')
cor(fittedscore,trainset.y[,2])
trainset.y[fittedscore[fittedscore<0]]
box1=as.vector(unlist(trainset.y[fittedscore<(-0.1),'G3']))
box2=as.vector(unlist(trainset.y[fittedscore>(-0.1)&&fittedscore<0,'G3']))
box3=as.vector(unlist(trainset.y[fittedscore>0,'G3']))
boxplot(box1,box2,box3,names = c('Gscore<-3','-3<Gscore<3','Gscore>3'))
trainset.y

fit<-lm(advance~ RC1 +RC3 +RC2 +RC4+ RC6+ RC5+RC9+RC7 +RC8,data=scores)
summary(fit)
modelscore<-fit$coefficients[-1]
fittedscore<-as.matrix(scores[,-c(9,10)])%*%modelscore
trainset.y[,1]=trainset.y[,1]
plot(fittedscore[order(trainset.y[,1])])
lines(trainset.y[order(trainset.y[,1]),1]-0.5,type = 'l')
cor(fittedscore,trainset.y[,1])
hist(fittedscore)
box1=as.vector(unlist(trainset.y[fittedscore<(-3),'G3']))
box2=as.vector(unlist(trainset.y[fittedscore>(-3)&&fittedscore<(3),'G3']))
box3=as.vector(unlist(trainset.y[fittedscore>3,'G3']))
boxplot(box1,box2,box3,names = c('Ascore<-3','-3<Ascore<3','Ascore>3'))

kernesti.regr(as.matrix(testscores[6,-c(10,11)]),as.matrix(scores[,-c(10,11)]),scores[,10], h=1.06*506^(-1/6), kernel="gauss", g=NULL, gernel="gauss", vect=FALSE)
testscores[6,10]

# Kmeans
par(c(1,1))
cltsc<-scale(trainset.x)
par(mfrow=c(1,1))
clt<-stats::kmeans(cltsc,3)
clt$cluster
plot(x=trainset.y[,2],y=trainset.y[,1],col=clt$cluster,pch=16,xlab = 'Final Score',ylab = 'Advance')


# cansor
matt

person=trainset.x[,c(1,7,8,10,11,13,16,17,18,19,20,21)]
trainset.x
fam=trainset.x[,-c(1,7,8,10,11,13,16,17,18,19,20,21)]
cancor(person,fam)

