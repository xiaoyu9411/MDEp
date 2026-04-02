source("functionsmdep.R")

ds<-read.csv("realdata.csv") #Data for PHQ-9 scores

colnames(ds)<-c("Study","a","q1","m","q3","b","n")


dstrue<-read.csv("truevalues.csv")#Include true mean and 95% CI
dstrue$SE<-(dstrue$U-dstrue$L)/(2*1.96)

model<-c()
modelqe<-c()
modelwpe<-c()
modelmdep<-c()

for (i in 1:nrow(ds)){
  dat<-ds[i,2:7]
  n<-as.numeric(dat[6])
  
  
  a1<-as.numeric(dat[1])/27
  q11<-as.numeric(dat[2])/27
  m1<-as.numeric(dat[3])/27
  q31<-as.numeric(dat[4])/27
  b1<-as.numeric(dat[5])/27
  
  
  ds$luo[i]<- (((0.7+0.39/n)*((q11+q31)/2))+((0.3-0.39/n)*m1))*27
  
  theta2<-2*qnorm((0.75*n-0.125)/(n+0.25))
  ds$shi[i] <-((q31-q11)/theta2)*27  
  
  ds$mqe[i]<-(as.numeric(qe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[1]))*27
  ds$sdqe[i]<-(as.numeric(qe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[2]))*27
  ss<-summary(qe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n))
  modelqe[i]<-rownames(ss)[1]
  
  ds$mwpe[i]<-(as.numeric(wpe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[1]))*27
  ds$sdwpe[i]<-(as.numeric(wpe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[2]))*27
  swpe<-wpe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)
  modelwpe[i]<-swpe$selected.dist
  
  
  ds$mmdep[i]<-(as.numeric(MDEp.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[1]))*27
  ds$sdmdep[i]<-(as.numeric(MDEp.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[2]))*27
  smdep<-MDEp.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)
  modelmdep[i]<-smdep$selected.dist
}


mymeta1<-metagen(TE=luo,seTE=shi/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta1)

mymeta2<-metagen(TE=mqe,seTE=sdqe/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta2)

mymeta3<-metagen(TE=mwpe,seTE=sdwpe/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta3)

mymeta4<-metagen(TE=mmdep,seTE=sdmdep/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta4)

mymeta5<-metagen(TE=Mean,seTE=SE,studlab=Study,data=dstrue,method.tau="REML")
summary(mymeta5)

#table making
extract_meta_row <- function(meta_obj, method_name) {
  s <- summary(meta_obj)
  
  data.frame(
    Methods = method_name,
    `Pooled Estimate (95% CI)` = sprintf(
      "%.2f\n(%.2f, %.2f)",
      s$random[1],
      s$lower.random[1],
      s$upper.random[1]
    ),
    tau = sprintf("%.2f", s$tau),
    I2 = sprintf("%.1f%%", s$I2 * 100),
    stringsAsFactors = FALSE
  )
}

table_S2 <- rbind(
  extract_meta_row(mymeta5, "True sample mean/SD"),
  extract_meta_row(mymeta1, "Luo/Shi (Wan)"),
  extract_meta_row(mymeta2, "QE"),
  extract_meta_row(mymeta3, "wPE"),
  extract_meta_row(mymeta4, "MDEp")
)

write.csv(table_S2, "table_S2.csv", row.names = FALSE)

