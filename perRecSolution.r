#Script to do lab 7 modified by JRB 3-31-20
{
#********************Read in data
mat<-read.delim("matsched.txt");
Sdat<-read.delim("Sel.txt");

# pull out or calculate values needed for YPR and SSBR
M <- 0.2; #from Gabriel et al.
mass <- mat[,3]; #returns mass at age near start of year
#calculate mass at harvest
massh<-c((mass[1:14]+mass[2:15])/2, mass[15]);  #updated some calcs to avoid loops

#*****************/Critical age stuff
#Instantaneous growth rates (to evaluate critical age)
G<-log(mass[2:15]/mass[1:14]);

# Based on this G<M sometime during age-5 but given use of constant massh
#   it would be slightly after harvest, so ypr would increase with F on ages over 6 but 
#   decrease with F on ages under 6 
# Plot G vs age with M as dashed ref line
age<-1:14;
plot(G~age,type="l", lty=1); 
abline(h=0.2,col=4,lty=2);
legend(x= 10.0 , y= 1.0, c("G", "M"), col=c("black","blue"), lty=c(1,2),lwd=2);


#********************YPR function
#arguments are Fval = fully selected F, Sv = selectivity at age, M is natural morality rate, and mh=mass at age for harvest
YPRf <-function(Fval,Sv,M,mh){
  maxa <-length(Sv); #set maxa to equal the number of ages based on length of selectivity vectore
  Nv <- numeric(maxa); #create vector with length or maxa
  Fv <-Fval*Sv;  #calculate F at age 
  Zv <- Fv+M;   #calculate Z at age 
  Nv[1]=1;  #per recruit calculation so start with 1 recruit
  for (i in 2:maxa) Nv[i] <- Nv[i-1]*exp(-Zv[i-1]);   #Analytical solution for exponential mortality diffeq
  Cv <- (Fv/Zv)*Nv*(1-exp(-Zv));   #Baranov catch equation-catch by age per recruit
  Yv <- Cv*mh;
  sum(Yv); #return YPR
}

#********************Calculate YPR for each mesh for range of fully selected F
YPR_res<-matrix(nrow=20,ncol=6);
out<-0;
#This is currently set to go to F=2.0 
# to determine shape and useful range to plot
for (Fval in seq(from=0.01, to=2.0, length.out=20)){
  out<-out+1; 
  YPR_res[out,1]<-Fval;
  for (seli in 2:6){
       YPR_res[out,seli] <-  YPRf(Fval,Sdat[,seli],M, massh);
  }
}


colnames(YPR_res)<-(c("F","m1.5in","m2.5in","m3.5in", "m4.5in", "m5.5in"));
YPR_res<-as.data.frame(YPR_res);
YPR_res;

#plot YPR versus Fully selected F for the selected selectivity vector
plot(YPR_res$m1.5in~YPR_res$F,col="black", xlab="F",ylab="YPR",ylim=c(0,1.3),type="l", lwd=3) ;
lines(YPR_res$m2.5in~YPR_res$F, col="navy", lwd=3) ;
lines(YPR_res$m3.5in~YPR_res$F, col="red", lwd=3) ;
lines(YPR_res$m4.5in~YPR_res$F, col="green",lwd=3);
lines(YPR_res$m5.5in~YPR_res$F, col="orange",lwd=3) ; 
legend(x= 0.3 , y= .6 , c("1.5 in.", "2.5 in.", "3.5 in.", "4.5 in.", "5.5 in."), 
       col=c("black","navy","red","green","orange"), lty=c(1,1,1,1,1),lwd=3);


#******************************SSBR function 
SSBRf<-function(Fval,Sv,M,mass,Pmat){
  maxa<-length(Sv);
  SSBv <-numeric(maxa);
  Nv <- numeric(maxa);  #See comment below for why this is duplicated
  
  Fv <-Fval*Sv; # We already did Fv, Zv, Nv but I am duplicating since you would need these in function if you pass it selectivity and fully selected F
  Zv <- Fv+M;
  Nv[1]=1;
  for (i in 2:maxa) Nv[i] <- Nv[i-1]*exp(-Zv[i-1]);   #Analytical solution for exponential mortality diffeq
  Nsp <- Nv*exp(-0.2*Fv-0.167*M); 
  SSBv <- Nsp*mass*Pmat;
  sum(SSBv);
}


#********************Calculate SSBR for each mesh for range of fully selected F
SSBR_res<-matrix(nrow=20,ncol=6);
out<-0;
for (Fval in seq(from=0.01, to=2, length.out=20)){
  out<-out+1 ;
  SSBR_res[out,1]<-Fval;
  for (seli in 2:6){
  SSBR_res[out,seli] <-  SSBRf(Fval,Sdat[,seli],M, mass, Pmat);
}
}

colnames(SSBR_res)<-(c("F","m1.5in","m2.5in","m3.5in", "m4.5in", "m5.5in"));
SSBR_res<-as.data.frame(SSBR_res);
SSBR_res;

#plot SSBR versus Fully selected F for the selected selectivity vector
plot(SSBR_res$m1.5in~SSBR_res$F,col="black", xlab="F",ylab="SSBR",ylim=c(0,10.0),type="l", lwd=3) ;
lines(SSBR_res$m2.5in~SSBR_res$F, col="navy", lwd=3); 
lines(SSBR_res$m3.5in~SSBR_res$F, col="red", lwd=3); 
lines(SSBR_res$m4.5in~SSBR_res$F, col="green",lwd=3);
lines(SSBR_res$m5.5in~SSBR_res$F, col="orange",lwd=3)   
legend(x= 1 , y= 8 , c("1.5 in.", "2.5 in.", "3.5 in.", "4.5 in.", "5.5 in."), 
       col=c("black","navy","red","green","orange"), lty=c(1,1,1,1,1),lwd=3);

#Reference point results --- Fmax ***************************************
Fmax_res<-matrix(nrow=5,ncol=4);
#Find Fmax for 5 mesh sizes
#We use the optimize function search within the interval 0,2.5 for the value of Fval that maximizes the YPRf function given the other parameters for YPRf that we pass
#The first argument to YPRf has to be the parameter that we are seeking to optimize for
for (i in 2:6){
Fmax_res[i-1,1]<-i-0.5; #mesh size
Fmax_res[i-1,2];
tmp<-optimize(YPRf, interval=c(0,2.5),Sv=Sdat[,i],M=0.2,mh=massh, maximum=T);  #Fmax for 1.5 mesh size
Fmax_res[i-1,2]<-tmp$max; #Fmax
Fmax_res[i-1,3]<-tmp$objective; #ypr at Fmax
Fmax_res[i-1,4]<-SSBRf(Fmax_res[i-1,2],Sdat[,i],M, mass, Pmat);
}

#label and clean up Fmax results
Fmax_res<-as.data.frame(Fmax_res);
colnames(Fmax_res)<-c("Mesh", "Fmax", "YPR", "SSBR");
Fmax_res<-format(Fmax_res,digits=3);
write.table(Fmax_res,row.names=F,quote=F);


#Reference point results --- F0.1 ***************************************
#Finite difference function to approx deriv of YPR w/ respect to F
DerYPR <- function(Fval,Sv,M,mh,delF){
  YPR1 <- YPRf(Fval,Sv,M,mh); #calculate YPR for Fval
  YPR2 <- YPRf(Fval+delF,Sv,M,mh); #Calculate YPR for Fval+delF
  (YPR2-YPR1)/delF; #return der YPR with respect to Fval using finite forward difference approx
}
#Modified version of deriv function that calulates squared dev from from target
devDer <- function(Fval,Sv,M,mh,delF,targder){
  YPR1 <- YPRf(Fval,Sv,M,mh); #calculate YPR for Fval
  YPR2 <- YPRf(Fval+delF,Sv,M,mh); #Calculate YPR for Fval+delF
  ((YPR2-YPR1)/delF - targder)^2; #return sq dev between der and target der
}

F01_res<-matrix(nrow=5,ncol=4);
for (i in 2:6){
  F01_res[i-1,1]<-i-0.5; #mesh size
  targd <- DerYPR(0,Sv=Sdat[,i],M,massh,0.0001)*0.1; 
  # note in next line the search goes up to Fmax but not beyond
  F01_res[i-1,2]<-optimize(devDer, interval=c(0,Fmax_res[i-1,"Fmax"]),Sv=Sdat[,i],
                           M=0.2,mh=massh,delF=0.0001,targder=targd)$minimum; #F0.1
  F01_res[i-1,3]<-YPRf(F01_res[i-1,2],Sv=Sdat[,i],M=M,mh=massh); #YPR at F0.1
  F01_res[i-1,4]<-SSBRf(Fval=F01_res[i-1,2],Sv=Sdat[,i],M=M,mass=mass,Pmat=Pmat); #SSBR at F0.1
}

#label and clean up F0.1 results;
Fmax_res<-as.data.frame(Fmax_res);
colnames(F01_res)<-c("Mesh", "F_0.1", "YPR", "SSBR");
F01_res<-as.data.frame(F01_res);
F01_res<-format(F01_res,digits=3);
write.table(F01_res,row.names=F,quote=F);

#Reference point results --- FX% ***************************************

#devSSBR modified from SSBRf to return the squared deviation in SSBR (for Fval) versus a specified target value
devSSBR<-function(Fval,Sv,M,mass,Pmat,tarSSBR){
  (SSBRf(Fval,Sv,M,mass,Pmat)-tarSSBR)^2;
}

#Calculate SSBR at F=0 targets for 35% and 50%
SSBR0<-SSBRf(0,Sdat[,2],M,mass,Pmat); #Spawning stock biomass/R with no fishing
SSBR35<-SSBR0*.35; #35% of unfished spawning stock biomass/R 
SSBR50<-SSBR0*.5; #50% of unfished spawning stock biomass/R

F35_res<-matrix(nrow=5,ncol=4);
F50_res<-matrix(nrow=5,ncol=4);
for (i in 2:6){
F35_res[i-1,1]<-i-0.5;
F50_res[i-1,1]<-i-0.5;
F35_res[i-1,2]<-optimize(devSSBR, interval=c(0,2),Sv=Sdat[,i],M=0.2,mass=mass,Pmat=Pmat,tarSSBR=SSBR35)$minimum ;
F50_res[i-1,2]<-optimize(devSSBR, interval=c(0,2),Sv=Sdat[,i],M=0.2,mass=mass,Pmat=Pmat,tarSSBR=SSBR50)$minimum ;
F35_res[i-1,3]<-YPRf(Fval=F35_res[i-1,2],Sv=Sdat[,i],M=0.2,mh=massh);
F50_res[i-1,3]<-YPRf(Fval=F50_res[i-1,2],Sv=Sdat[,i],M=0.2,mh=massh);
F35_res[i-1,4]<-SSBRf(Fval=F35_res[i-1,2],Sv=Sdat[,i],M=0.2,mass=mass,Pmat=Pmat);
F50_res[i-1,4]<-SSBRf(Fval=F50_res[i-1,2],Sv=Sdat[,i],M=0.2,mass=mass,Pmat=Pmat);
}
 
#label and clean up F35% results
colnames(F35_res)<-c("Mesh", "F35%", "YPR", "SSBR");
F35_res<-as.data.frame(F35_res);
F35_res<-format(F35_res,digits=3);
write.table(F35_res,row.names=F,quote=F);

#label and clean up F50% results
colnames(F50_res)<-c("Mesh", "F50%", "YPR", "SSBR");
F50_res<-as.data.frame(F50_res);
F50_res<-format(F50_res,digits=3);
write.table(F50_res,row.names=F,quote=F);

}










