# This version of EM algorithm for SPCD includes dependence on baseline values for Stage I and on Y1 for Stage II;
# EM function is a single iteration of EM algorithm;
# doEM function iterates EM function till convergence;

#INPUT:
# d - dataset with following columns: g1 (treatment Stage I: 1-drug, 0-placebo), g2 (treatment Stage II: 1-drug, 0-placebo),
# y01 (baseline Stage I), y02 (baseline Stage II), y1 (change Stage I), y2 (change Stage II);
# d HAS TO BE SORTED by g1 so that np Stage 1 palcebo subjects are in the first np rows; 
# s - parameter vector initial values:
# (b01,b11,s1,b02,b12,s21,b03,b13,s01,b04,b14,b24,s201,b05,b15,s02,b06,b16,b26,s202,pr,
# m01,m02,v01,cov01,v11,v02,cov02,v12,v03,cov03,v13,v04,v14,v24,v05,cov05,v15,v06,v16,v26,
# E1_em,SE1_em,P1_em,E2_em,SE2_em,P2_em,E_em,SE_em,P_em);
# Values of parameters np (number of Stage 1 placebo subjects) and spcd_w (SPCD weight for Stage I) have to be specified
# before running the script.

#OUTPUT:
# $data - dataset;
# $parameters - 50 by 1 vector of estimated parameters and treatment effects;
# $response - N by 1 vector of estimated response probabilities (use first np values).

#TO RUN:
# np<-dim(my.data[my.data$g1==0,])[1]
# my.data<-my.data[order(my.data$g1,decreasing=F),]
# spcd_w<-0.7
# my.s<-c(rep(0.5,20),0.5,0.5,0.5,rep(0.5,27))
# myEM<-doEMy(my.data,my.s)

EMy = function(d,s) {

#Center Y01;
d$y0c<-d$y01-mean(d$y01)

#Center Y1-Y3;
d$y1c<-d$y1-mean(d$y1)
	
#Response Probability;
# Stage I data only;

A=dnorm(d$y1,s[7]+s[8]*d$y0c,sqrt(s[9]))
B=dnorm(d$y1,s[14]+s[15]*d$y0c,sqrt(s[16]))

pR=s[21]*A/(s[21]*A+(1-s[21])*B)
s[21]<-mean(pR[1:np],na.rm=T)
pR[which(is.na(pR))]<-0

######################################## estimates #############################################

#Stage 1 Drug;
s[1]<-sum(d$g1*(d$y1-s[2]*d$y0c))/sum(d$g1)
s[2]<-sum(d$g1*d$y0c*(d$y1-s[1]))/sum(d$g1*(d$y0c**2))
s[3]<-sum(d$g1*((d$y1-s[1]-s[2]*d$y0c)**2))/sum(d$g1)

#Stage 2 Drug;
s[4]<-sum(d$g1*(d$y2-s[5]*d$y1c))/sum(d$g1)
s[5]<-sum(d$g1*d$y1c*(d$y2-s[4]))/sum(d$g1*(d$y1c**2))
s[6]<-sum(d$g1*((d$y2-s[4]-s[5]*d$y1c)**2))/sum(d$g1)

#Stage 1 Placebo Responders (A);
s[7]<-sum((1-d$g1)*pR*(d$y1-s[8]*d$y0c))/sum((1-d$g1)*pR)
s[8]<-sum((1-d$g1)*pR*d$y0c*(d$y1-s[7]))/sum((1-d$g1)*pR*(d$y0c**2))
s[9]<-sum((1-d$g1)*pR*((d$y1-s[7]-s[8]*d$y0c)**2))/sum((1-d$g1)*pR)

#Stage 2 Placebo Responders (A);
s[10]<-sum((1-d$g1)*pR*(d$y2-s[11]*d$y1c-s[12]*d$g2))/sum((1-d$g1)*pR)
s[11]<-sum((1-d$g1)*pR*d$y1c*(d$y2-s[10]-s[12]*d$g2))/sum((1-d$g1)*pR*(d$y1c**2))
s[12]<-sum((1-d$g1)*pR*d$g2*(d$y2-s[10]-s[11]*d$y1c))/sum((1-d$g1)*pR*(d$g2**2))
s[13]<-sum((1-d$g1)*pR*((d$y2-s[10]-s[11]*d$y1c-s[12]*d$g2)**2))/sum((1-d$g1)*pR)

#Stage 1 Placebo Non-Responders (B);
s[14]<-sum((1-d$g1)*(1-pR)*(d$y1-s[15]*d$y0c))/sum((1-d$g1)*(1-pR))
s[15]<-sum((1-d$g1)*(1-pR)*d$y0c*(d$y1-s[14]))/sum((1-d$g1)*(1-pR)*(d$y0c**2))
s[16]<-sum((1-d$g1)*(1-pR)*((d$y1-s[14]-s[15]*d$y0c)**2))/sum((1-d$g1)*(1-pR))

#Stage 2 Placebo Non-Responders (B);
s[17]<-sum((1-d$g1)*(1-pR)*(d$y2-s[18]*d$y1c-s[19]*d$g2))/sum((1-d$g1)*(1-pR))
s[18]<-sum((1-d$g1)*(1-pR)*d$y1c*(d$y2-s[17]-s[19]*d$g2))/sum((1-d$g1)*(1-pR)*(d$y1c**2))
s[19]<-sum((1-d$g1)*(1-pR)*d$g2*(d$y2-s[17]-s[18]*d$y1c))/sum((1-d$g1)*(1-pR)*(d$g2**2))
s[20]<-sum((1-d$g1)*(1-pR)*((d$y2-s[17]-s[18]*d$y1c-s[19]*d$g2)**2))/sum((1-d$g1)*(1-pR))

s[22]<-s[7]+s[8]*sum((1-d$g1)*pR*d$y0c)/sum((1-d$g1)*pR)
s[23]<-s[14]+s[15]*sum((1-d$g1)*(1-pR)*d$y0c)/sum((1-d$g1)*(1-pR))

######################################## variance ##############################################

#Stage 1 Drug;
s11<-sum(d$g1)/s[3]
s21<-sum(d$g1*d$y0c)/s[3]
s12<-sum(d$g1*d$y0c)/s[3]
s22<-sum(d$g1*(d$y0c**2))/s[3]
v<-solve(matrix(c(s11,s21,s12,s22),2,2))
s[24]<-v[1,1]
s[25]<-v[1,2]
s[26]<-v[2,2]

#Stage 2 Drug;
s11<-sum(d$g1)/s[6]
s21<-sum(d$g1*d$y1c)/s[6]
s12<-sum(d$g1*d$y1c)/s[6]
s22<-sum(d$g1*(d$y1c**2))/s[6]
v<-solve(matrix(c(s11,s21,s12,s22),2,2))
s[27]<-v[1,1]
s[28]<-v[1,2]
s[29]<-v[2,2]

#Stage 1 Placebo Responders;
s11<-sum((1-d$g1)*pR)/s[9]
s21<-sum((1-d$g1)*pR*d$y0c)/s[9]
s12<-sum((1-d$g1)*pR*d$y0c)/s[9]
s22<-sum((1-d$g1)*pR*(d$y0c**2))/s[9]
v<-solve(matrix(c(s11,s21,s12,s22),2,2))
s[30]<-v[1,1]
s[31]<-v[1,2]
s[32]<-v[2,2]

#Stage 2 Placebo Responders;
s11<-sum((1-d$g1)*pR)/s[13]
s21<-sum((1-d$g1)*pR*d$y1c)/s[13]
s31<-sum((1-d$g1)*pR*d$g2)/s[13]
s12<-sum((1-d$g1)*pR*d$y1c)/s[13]
s22<-sum((1-d$g1)*pR*(d$y1c**2))/s[13]
s32<-sum((1-d$g1)*pR*d$y1c*d$g2)/s[13]
s13<-sum((1-d$g1)*pR*d$g2)/s[13]
s23<-sum((1-d$g1)*pR*d$y1c*d$g2)/s[13]
s33<-sum((1-d$g1)*pR*(d$g2**2))/s[13]
v<-solve(matrix(c(s11,s21,s31,s12,s22,s32,s13,s23,s33),3,3))
s[33]<-v[1,1]
s[34]<-v[2,2]
s[35]<-v[3,3]

#Stage 1 Placebo Non-Responders;
s11<-sum((1-d$g1)*(1-pR))/s[16]
s21<-sum((1-d$g1)*(1-pR)*d$y0c)/s[16]
s12<-sum((1-d$g1)*(1-pR)*d$y0c)/s[16]
s22<-sum((1-d$g1)*(1-pR)*(d$y0c**2))/s[16]
v<-solve(matrix(c(s11,s21,s12,s22),2,2))
s[36]<-v[1,1]
s[37]<-v[1,2]
s[38]<-v[2,2]

#Stage 2 Placebo Non-Responders;
s11<-sum((1-d$g1)*(1-pR))/s[20]
s21<-sum((1-d$g1)*(1-pR)*d$y1c)/s[20]
s31<-sum((1-d$g1)*(1-pR)*d$g2)/s[20]
s12<-sum((1-d$g1)*(1-pR)*d$y1c)/s[20]
s22<-sum((1-d$g1)*(1-pR)*(d$y1c**2))/s[20]
s32<-sum((1-d$g1)*(1-pR)*d$y1c*d$g2)/s[20]
s13<-sum((1-d$g1)*(1-pR)*d$g2)/s[20]
s23<-sum((1-d$g1)*(1-pR)*d$y1c*d$g2)/s[20]
s33<-sum((1-d$g1)*(1-pR)*(d$g2**2))/s[20]
v<-solve(matrix(c(s11,s21,s31,s12,s22,s32,s13,s23,s33),3,3))
s[39]<-v[1,1]
s[40]<-v[2,2]
s[41]<-v[3,3]

###################################### treatment effects ########################################

#Stage 1;
#drug;
del1<-s[1]+s[2]*sum(d$g1*d$y0c)/sum(d$g1)
var1<-s[24]+s[26]+2*s[25]
#placebo;
del0<-s[21]*s[22]+(1-s[21])*s[23]
var01<-s[30]+s[32]+2*s[31]
var02<-s[36]+s[38]+2*s[37]
var0<-(s[21]**2)*var01+((1-s[21])**2)*var02
#effect;
del_s1<-del1-del0
var_s1<-var1+var0

s[42]<-del_s1
s[43]<-sqrt(var_s1)
s[44]<-2*pnorm(abs(del_s1/sqrt(var_s1)),lower.tail=F)

#Stage 2;
#beta s[12] or s[19];
#var s[35] or s[41];
del_s2<-ifelse(s[23]>s[22],s[19],s[12]) 
var_s2<-ifelse(s[23]>s[22],s[41],s[35])

s[45]<-del_s2
s[46]<-sqrt(var_s2)
s[47]<-2*pnorm(abs(del_s2/sqrt(var_s2)),lower.tail=F)

#SPCD effect;
#E_em;
s[48]<-spcd_w*del_s1+(1-spcd_w)*del_s2
#SE_em;
s[49]<-sqrt((spcd_w**2)*var_s1+((1-spcd_w)**2)*var_s2)
#P_em;
s[50]<-2*pnorm(abs(s[48]/s[49]),lower.tail=F)

z<-list(d,s,pR)
names(z)<-c("data","parameters","response")
return(z)
}

doEMy = function(d,s){
  s1 <- s
  emy <- EMy(d,s);

  s2 <- emy$parameters
  r <- emy$response

  while (max(abs(s1-s2))>0.0001) {
    s1 <- s2
    emy <- EMy(d,s1);
    s2 <- emy$parameters
    r <- emy$response
  }
  
  z<-list(d,s2,r);
  names(z)<-c("data","parameters","response");
  return(z);
}
