####book (Ron Wehrens)

#first part Preliminaries
install.packages("ChemometricsWithR")
library("ChemometricsWithRData")
data(gasoline, package = "pls")
wavelengths <- seq(900, 1700, by = 2)
matplot(wavelengths, t(gasoline$NIR), type = "l",
        lty = 1, xlab = "Wav", ylab = "1/R")

install.packages("kohonen")
library("kohonen")
data(wines, package = "kohonen")
colnames(wines)
table(vintages)
pairs(wines[, 1:3], col = wine.classes)
 
data(Prostate2000Raw, package = "msProstate")
plot(Prostate2000Raw$mz, Prostate2000Raw$intensity[,1],
     type = "h", xlab = "m/z", ylab = "Int", main = "Prostate data")


##列转换成矩阵
install.packages("kohonen")
library(kohonen)
classes <- c(rep(1, 5), rep(2, 7), rep(3, 9))
classes
classmat <- classvec2classmat(classes)
classmat
classmat2classvec(classmat) 


##
install.packages(c('RcppArmadillo', 'bigmemory'))
install.packages('hiblup_1.2.0.zip', repos = NULL)
library(hiblup)
suppressMessages(library("hiblup"))
data("hidata")
X <- model.matrix(~as.factor(Sex), data = pheno) # fixed effects
R <- as.matrix(pheno$Sire) # random effects
# R can be either character or numeric. For interaction between two or
# more random effects, it can be fitted by pasting them together, for
# example, there are two random effects R1 and R2, we could fit their
# interaction in the model as: R=cbind(R1,R2,paste(R1,R2,sep='_')).
gebv <- hiblup(pheno = pheno[, c(1, 5)], geno = geno, map = map, geno.id = geno.id,
               pedigree = pedigree, vc.method = c("HI"), mode = "A", CV = X, R = R,
               snp.solution = TRUE)
install.packages("rMVP")
library(rMVP)

##########################################################
###############L R Schaeffer 2019 fall 课程  #####
A = matrix(data=c(3,-1,-2,4),byrow=TRUE,ncol=2)
B = matrix(data=c(1,-1,3,1,1,-1),byrow=TRUE,ncol=3)
M = A %*% B  #矩阵相乘
M
M = A %x% B # 矩阵直积
M
###########直和功能
block= function( ... ) {
  argv = list( ... )
  i = 0
  for( a in argv ){
    m = as.matrix(a)
    if(i == 0)
      rmat = m
    else
    {
      nr = dim(m)[1]
      nc = dim(m)[2]
      aa = cbind(matrix(0,nr,dim(rmat)[2]),m)
      rmat = cbind(rmat,matrix(0,dim(rmat)[1],nc))
      rmat = rbind(rmat,aa)
    }
    i = i+1
  }
  rmat
}
M = block(A,B) ###直和
M 
###连接
SA = c(23, 14, 38, 54, 17)
SB = c(1,-1,1,-1,1)
M1 = cbind(SA,SB) # order 5 x 2
M2 = rbind(SA,SB) # order 2 x 5


A = matrix(data=(1:10000), byrow=T,  ncol=100)
# keep only rows where first element
# is greater than 10
B = A[A[ ,1]>1500, ]
# keep rows 4,5, and 9, and columns 21 to 30
kr = c(4, 5, 9)
kc = c(21:30)
C = A[kr,kc]



S = c( 3, 6, -1, 2, 11, 4, 5)
ka = order(S) # ascending
kd = order(-S) # descending
ka
kd
S[ka]  ##升序
S[kd]  ##降序


##A的逆矩阵
AINV = function(sid,did,bi){
  # IDs 假定为连续编号的ID
  rules=matrix(data=c(1,-0.5,-0.5,-0.5,
                      0.25,0.25,-0.5,0.25,0.25),
               byrow=TRUE,nrow=3)
  nam = length(sid) : np = nam + 1
  ss = sid + 1 : dd = did + 1
  LAI = matrix(data=c(0),nrow=np,ncol=np)
  for(i in 1:nam){
    ip = i + 1 : X = 1/bi[i]
    k = cbind(ip,ss[i],dd[i])
    LAI[k,k] = LAI[k,k] + rules*X }
  k = c(2:np) : C = LAI[k,k]
  return(C) }

######y = Xb + Zu + e  ##
##  MME 方程
MME = function(X,Z,GI,RI,y){
  XX = t(X) %*% RI %*% X ; XZ = t(X) %*% RI %*% Z
  ZZ = t(Z) %*% RI %*% Z ; Xy = t(X) %*% RI %*% y
  Zy = t(Z) %*% RI %*% y
  N = length(y); R1 = cbind(XX,XZ)
  R2 = cbind(t(XZ),(ZZ+GI))
  LHS = rbind(R1,R2)
  RHS = rbind(Xy,Zy)
  C = ginv(LHS)
  bhat = C %*% RHS
  SSR = t(bhat) %*% RHS; VPE = diag(C)
  sep = matrix(data=VPE,ncol=1)
  return(list(LHS=LHS,RHS=RHS,SSR=SSR,C=C,
              VPE=sep,SOLNS=bhat)) }

########一个实例########
# 29 animals           #
# 15 observations      #
# 4 contemporary groups#
# 2 years              #
## X 矩阵的构建（LHS） #
########################
desgn=function(v,nc){
  if(is.numeric(v)){
    va = v; mrow = length(va)
    mcol = max(va)
    if(nc > mcol)mcol = nc }
  if(is.character(v)){
    vf = factor(v); va = as.numeric(vf)
    mrow = length(v)
    mcol = length(levels(vf))
    if(nc > mcol)mcol = nc }
  X = matrix(data=c(0),nrow=mrow,ncol=mcol)
  for(i in 1:mrow){
    ic = va[i]; X[i,ic] = 1 }
  return(X) }
##years
years=c(1,1,1,1,2,2,2,1,1,1,1,2,2,2,2)
ny=0
X = desgn(years,ny)
##### classvec2classmat(years) 与作者的自己功能一样
## contemporary groups
cg=c(1,1,1,1,2,2,2,3,3,3,3,4,4,4,4)
nc=0
W = desgn(cg,nc) # 15 x 4 matrix
#y
y = matrix(data=c(94, 89, 72, 100, 73, 70, 84,
                  88,102,82,130,93,105,118,69),ncol=1)
##Animal Additive Effects
awr=c(15:29)
anim=c(1:29)
na=29
Zwr = desgn(awr,na)
Zo = matrix(data=c(0),nrow=15, ncol=14)
        ##Za = cbind(Zo,Zwr) # 15 x 29
Z = cbind(W,Zwr) # CG and Anim. Add.

sire=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,
       1,1,3,3,13,5,5,1,1,7,7,5,5,9,11)
dam=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      2,2,4,6,2,4,4,8,10,12,14,8,10,12,14)
# None of the animals is inbred
bi=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     .5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5)
AI = AINV(sire, dam, bi)
AAI = AI *(144/64)
HI = diag(c(1,1,1,1))*(144/36)
GI = block(HI,AAI)
RI = id(15) # identity 15 x 15


EX31 = MME(X,Z,RI,GI,y)
EX31$SOLNS


### LR.Schaeffer RRM 2016 book
## LP 的公式计算

LPOLY = function(no) {
  if(no > 9 ) no = 9
  nom = no - 1
  phi = matrix(data=c(0),nrow=9,ncol=9)
  phi[1,1]=1
  phi[2,2]=1
  for(i in 2:nom){
    ia = i+1
    ib = ia - 1
    ic = ia - 2
    c = 2*(i-1) + 1
    f = i - 1
    c = c/i
    f = f/i
    for(j in 1:ia){
      if(j == 1){ z = 0 }
      else {z = phi[ib,j-1]}
      phi[ia,j] = c*z - f*phi[ic,j]
    }
  }
  for( m in 1:no){
    f = sqrt((2*(m-1)+1)/2)
    phi[m, ] = phi[m, ]*f
  }
  return(phi[1:no,1:no])
}
LP4 <- LPOLY(4)

####R中的一种将上半个矩阵更改为一个向量的方法
hsmat <- function(vcvfull) {
  mord = nrow(vcvfull)
  np = (mord *(mord + 1))/2
  desg = rep(0,np)
  k = 0
  for(i in 1:mord){
    for(j in i:mord){
      k = k + 1
      desg[k] = vcvfull[i,j] } }
  return(desg) }
hsmat(LP4)

##降低后相似比较ss值，越小越好
#计算ss
N=10000
can=c(1:N)*0
VR=V
nocov = 2 # order of fit + 1
phr = PH[ ,c(1:nocov)]
PVP = t(phr)%*%VR%*%phr
PP = t(phr)%*%phr
PPI=ginv(PP)
Kr = PPI%*%PVP%*%PPI
ndf=199
for(ko in 1:N){
  Ka = rWishart(1,ndf,Kr)/ndf
  Kb = Ka[, ,1]
  Vr = phr%*%Kb%*%t(phr)
  DEL = Vr - VR
  er = hsmat(DEL)
  vh = hsmat(VR)
  vv = sum(vh*vh)
  ssr = sum(er*er)/vv
  can[ko] = ssr
} #end of samples
hist(can,breaks=50)
##进行比较
kb = order(-can)
ncan = can[kb]
ncan[1:10]
kc=which(ncan < 0.4592)
prob = 0
if(length(kc)>0)prob = 1 - (kc[1]/length(ncan))
prob


###使矩阵变为positive matrix
# Let A be the matrix to be made p.d.
# no be the order of the matrix
E = eigen(A)
ev = E$values
U = E$vectors
no = dim(A)[1]
nev = which(ev < 0)
wr = 0
k=length(nev)
if(k > 0){
  p = ev[no - k]
  B = sum(ev[nev])*2.0
  wr = (B*B*100.0)+1
  val = ev[nev]
  ev[nev] = p*(B-val)*(B-val)/wr
  A = U%*%diag(ev)%*%t(U)
}

##p 3.1  48-49页
y = c(38,37,35,27,40,40,28,42,37,22,39,36,33,
      25,41,38,23,37,35,30,24)
days=c(4,7,16,25,2,11,21,6,17,25,5,15,19,24,
       3,14,22,4,9,17,23)
dgrp=c(1,2,4,5, 1,3,5, 2, 4, 5,1, 3, 4, 5,1,
       3, 5,1,2,4,5)
# Animals with records
anw = c(7,7,7,7,8,8,8,9,9,9,10,10,10,10,11,11,
        11,12,12,12,12)
# Gender codes for each observation
gend = c(1,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,1,1,
         1,1,1)
# Contemporary group levels
cg = c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
       2,2)
# Sires and dams of all animals (first 6 are unknown)
sirs=c(0,0,0,0,0,0,1,2,3,1,2,3)
dams=c(0,0,0,0,0,0,4,4,5,5,6,6)
# bii values = 0.5 - 0.25*(Fsire + Fdam)
bii=c(1,1,1,1,1,1,.5,.5,.5,.5,.5,.5)
# Initial residual variances
vare=c(1,1,1,1,1,.97,.97,.97,.97,.97,.95,.95,.95,.95,.95,
       .93,.93,.93,.93,.93,.90,.90,.90,.90,.90)
length(y)
length(days) # check that lengths are the same

# divide days by 100 to reduce magnitude
alld=c(1:25)/100
all2=alld*alld
alle=alld*0 + 1
fT=cbind(alle,alld,all2)
# Legendre polynomials
LAM=LPOLY(3) ## oder 3
LAM
ti=c(1:25)
tmin=1
tmax=25# you could also use 25, in this case
qi = 2*(ti - tmin)/(tmax - tmin) - 1
qi

x1=qi
x0=x1*0 + 1
x2=x1*x1
M=cbind(x0,x1,x2)
PH = M %*% t(LAM)
PH

## Design Matrices for Model
rdesgn = function(v,tim,fT,nc,no){
  if(is.numeric(v)){
    va = v
    mrow = length(va)
    mcol = max(va)
    if(nc > mcol)mcol = nc }
  if(is.character(v)){
    vf = factor(v)
    va = as.numeric(vf)
    mrow = length(v)
    mcol = length(levels(vf))
    if(nc > mcol)mcol = nc }
  mcc = mcol*no
  X = matrix(data=c(0),nrow=mrow,ncol=mcc)
  for(i in 1:mrow){
    ic = (va[i]-1)*no
    jc = c((ic+1):(ic+no))
    X[i,jc] = fT[tim[i], ] }
  return(X) }
# design matrix for gender effects
Xg=rdesgn(gend,days,fT,2,3)
Xg
# design matrix for random factors using
# Legendre polynomials in PH
Zc=rdesgn(cg,days,PH,2,3)
Za=rdesgn(anw,days,PH,12,3)
Zp=Za[ ,c(19:36)] # PE design matrix is a
# subset of design matrix for genetic
# Check the dimensions of the matrices


##inve - A
# Routine to set up inverse of additive relationship
# matrix from list of sires and dams and bii-values
AINV = function(sid,did,bi){
  # IDs assumed to be consecutively numbered, and
  # parents come before progeny
  rules=matrix(data=c(1,-0.5,-0.5,
                      -0.5,0.25,0.25,
                      -0.5,0.25,0.25),
               byrow=TRUE,nrow=3)
  nam = length(sid)
  np = nam + 1
  ss = sid + 1
  dd = did + 1
  LAI = matrix(data=c(0),nrow=np,ncol=np)
  for(i in 1:nam){
    ip = i + 1
    X = 1/bi[i]
    k = cbind(ip,ss[i],dd[i])
    LAI[k,k] = LAI[k,k] + rules*X
  }
  k = c(2:np)
  C = LAI[k,k]
  return(C) }
AI = AINV(sirs,dams,bii)
# AI is 12 by 12 for the example

Kc = matrix(data=c(.51696, -.1623, -.0895,
                   -.1623, .3504, .1135,
                   -.0895, .1135, 1.20888),byrow=TRUE,ncol=3)

Ka = matrix(data=c(.0522, -.00170,-.00142,
                   -.00170,.0350,-.00149,
                   -.00142,-.00149,.121),byrow=TRUE,ncol=3)

Kp = matrix(data=c(.06,-.00635,.003753,
                   -.00635,.04,-.00106,
                   .003753,-.00106,.15),byrow=TRUE,ncol=3)

# Invert the covariance matrices
library(MASS)

Kci=ginv(Kc)
Kai=ginv(Ka)
Kpi=ginv(Kp)
# Residual variances
R = diag(vare[days])
RI=ginv(R)
# Set up covariance matrices for each factor
# Contemporary groups (2 of them
#install.packages("qdap")
library(qdap)
C=as.numeric(id(2))
CI = C %x% Kci # direct product, order 6 x 6

# Additive genetic (12 animals
dim(AI)
GI = AI %x% Kai # order 36 x 36
# Permanent Environmental
P=id(6)
P=as.numeric(id(6))
PI = P %x% Kpi # order 18 x 18
W = cbind(Zc,Za,Zp)
X = Xg
HI=block(CI,GI,PI)
# Uses block function, or direct sum

# Function to form MME
MME = function(X,Z,GI,RI,y){
  XX = t(X) %*% RI %*% X
  XZ = t(X) %*% RI %*% Z
  ZZ = t(Z) %*% RI %*% Z
  Xy = t(X) %*% RI %*% y
  Zy = t(Z) %*% RI %*% y
  N = length(y)
  R1 = cbind(XX,XZ)
  R2 = cbind(t(XZ),(ZZ+GI))
  LHS = rbind(R1,R2)
  RHS = rbind(Xy,Zy)
  # now solve
  C = ginv(LHS)
  bhat = C %*% RHS
  # estimate residual variance
  SSR = t(bhat) %*% RHS
  VPE = diag(C)
  sep = matrix(data=VPE,ncol=1)
return(list(LHS=LHS,RHS=RHS,SSR=SSR,C=C,
              VPE=sep,SOLNS=bhat)) }
##出错。。。。
SA = MME(Xg,W,HI,RI,y)


bh = SA$SOLNS
ghat = bh[c(1:6),]
chat = bh[c(7:12),]
ahat = bh[c(13:48),]
phat = bh[c(49:66),]
# must reformat the solutions
gh = matrix(data=ghat,byrow=TRUE,ncol=3)
gh # GENDER EFFECTS

ch = matrix(data=chat,byrow=TRUE,ncol=3)
ch # CONTEMPORARY GROUP EFFECTS

ah = matrix(data=ahat,byrow=TRUE,ncol=3)
ah # ANIMAL ADDITIVE GENETIC

ph = matrix(data=phat,byrow=TRUE,ncol=3)
ph # ANIMAL PERMANENT ENVIRONMENTAL
###plot
gsoln = gh%*%t(fT) # points along 25 days
gs1 = gsoln[1,]
gs2 = gsoln[2,]
par(bg="cornsilk")
plot(gs2,col="blue",lwd=5,type="l",xlab="Days on Test",
     ylab="Resistance Level")
title(main="Gender Trajectories")
lines(gs1,col="red",lwd=5)

### 3.6 Estimation of Covariance Matrices
#CG
pchi=rchisq(1,4)
# pchi=5.14983
Kc =t(ch)%*%ch/pchi
#AG
pchi=rchisq(1,14)
# pchi=10.42511
Ka =t(ah)%*%AI%*%ah/pchi
Ka
#PE
pchi=rchisq(1,8)
# pchi = 8.717388
Kp =t(ph)%*%ph/pchi
Kp
#R
T = cbind(Xg,W)
res = y - T%*%bh # residuals
sse = sum(res*res)
sse/19 # overall variance
#r分为5组
Q = desgn(dgrp,5)
idual=Q*cbind(res,res,res,res,res)
D = t(Q)%*%Q
DI = ginv(D)
ee = t(idual)%*%idual
RR=ee*DI
RR


##评估每天的遗传力
# CONTEMPORARY GROUPS
varc = PH%*%Kc%*%t(PH)
vc = diag(varc) # 25 variances, each day
# ADDITIVE GENETIC
vara = PH%*%Ka%*%t(PH)
va = diag(vara) # 25 variances, each day
# PERMANENT ENVIRONMENTAL
varp = PH%*%Kp%*%t(PH)
vp = diag(varp) # 25 variances, each day
# combine into one table
vtab = cbind(vc,va,vp)
# RESIDUAL VARIANCES
vres = diag(RR) # from earlier section
v1=vres[1]
v2=vres[2]
v3=vres[3]
v4=vres[4]
v5=vres[5]
R = c(v1,v1,v1,v1,v1, v2,v2,v2,v2,v2,
      v3,v3,v3,v3,v3,
      v4,v4,v4,v4,v4, v5,v5,v5,v5,v5)
vtt = cbind(vtab,R)
# total sum of individual variances by day
pvar = vtt[ ,1]+vtt[ ,2]+vtt[ ,3]+vtt[ ,4]
D=diag(pvar)
DI=ginv(D)
# Convert absolute values to percentages of pvar
HH=DI%*%vtt
Hc = HH[ ,1] # contemporary group
Ha = HH[ ,2]
Hp = HH[ ,3]
Hr = HH[ ,4]
par(bg="oldlace")
plot(Hc,type="l",lwd=3,col="red",xlab="Days on Test",
     ylab="Percentage of Variance",ylim=c(0,1))
title(main="Percentage of Variance over Days")
lines(Ha,lwd=3,col="blue")
lines(Hp,lwd=3,col="cyan")
lines(Hr,lwd=3,col="magenta")


###4.4 75 页 三个胎次的遗传方差随DIM变化
# Legendre polynomials
LAM=LPOLY(5)
ti=c(5:365)
tmin=5
tmax=365
qi = 2*(ti - tmin)/(tmax - tmin) - 1
x=qi
x0=x*0 + 1
x2=x*x
x3=x2*x
x4=x3*x
M=cbind(x0,x,x2,x3,x4)
PH = M %*% t(LAM)
Ka1 = matrix(data=c(8.1910, 0.2880,-0.6694,0.2360, -0.1407,
                    0.2880, 1.4534, -0.1327, 0.4590, 0.4926,
                    -0.6694, -0.1327, 0.5108, -0.1512, 0.0713,
                    0.2360, 0.4590, -0.1512, 0.1855, -0.0524,
                    -0.1407, 0.4926, 0.0713, -0.0524, 0.0766),byrow=TRUE,ncol=5)
Va1 = PH%*%Ka1%*%t(PH) # order 361 x 361
vg1 = diag(Va1)
# similar arrays for vg2, vg3, Ka2, Ka3 (not shown)
par(bg="cornsilk")
plot(vg1,col="blue",lwd=5,type="l",axes=FALSE,xlab="Days on Test",
     ylab="Genetic Variance",ylim=c(4,16))
axis(1,days)
axis(2)
title(main="Genetic Variances")
lines(vg2,col="red",lwd=5)
lines(vg3,col="darkgreen",lwd=5)
points(55,15,pch=0,col="blue",lwd=3)
text(55,15,"First Parity",col="blue",pos=4)
points(55,14,pch=0,col="red",lwd=3)
text(55,14,"Second Parity",col="red",pos=4)
points(55,13,pch=0,col="darkgreen",lwd=3)
text(55,13,"Third Parity",col="darkgreen",pos=4)


##算EBVs
install.packages("moonsun")
library(moonsun)
PH
ka=c(1:301)
P305 = PH[ka, ]
dim(P305)
C305 = t(P305)%*%jd(301,1)
C305
