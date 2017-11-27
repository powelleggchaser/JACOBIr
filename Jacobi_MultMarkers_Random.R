#Gaussian Seidel Residual Updater - Multiple Markers
#Rversion

###########################################################

##  FUNCTIONS

###########################################################

t_mult<-function(matrix){
  t_matrix<-t(matrix)
  prod_matrix<-t_matrix%*%matrix
  return(prod_matrix)
}

t_mult_dupl<-function(matrix1,matrix2){
  t_matrix<-t(matrix1)
  prod_matrix<-t_matrix%*%matrix2
  return(prod_matrix)
}

t_mult_diag<-function(matrix){
  t_matrix<-t(matrix)
  prod_matrix<-t_matrix%*%matrix
  diag_prod_matrix<-diag(prod_matrix)
  return(diag_prod_matrix)
}

dupl_mult_diag<-function(matrix1,matrix2){
  prod_matrix<-matrix2%*%matrix1
  diag_matrix<-diag(prod_matrix)
  return(diag_matrix)
}

lamda<-function(vare,vara){
  vare<-as.numeric(vare)
  vara<-as.numeric(vara)
  lamda=vare/vara
  return(lamda)
}

normalise <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

ssqd <- function(x,y) {
  ssqd <-sum((x-y)^2)
  return(ssqd)
}

###########################################################
###########################################################


### DATA ENTRY ###


#read in individuals (Data/No Data Score)
x<-read.csv("X_Jacobi.csv",sep=",",header=F)
x<-as.matrix(x)
class(x)<-"numeric"
ndata=nrow(x)

#read in genotypes
ped<-read.table("ped.csv",sep=",",header=T)
#install.packages('pedigree')
library('pedigree')
makeA(ped,which=ped$ID!=0)
a <- read.csv("A.csv",header=F)
a<-as.matrix(a)
class(a)<-"numeric"


# phenotypes
y<-read.csv('y_Jacobi.csv',sep=',',header=F)
y<-as.matrix(y)
class(y)<-"numeric"

#calculate allele frequencies & variances
#p=sum(z)/(2*ndata)
#q=1-p
#pq2=2*p*q
#var_a=pq2*(allele_effect^2)
#var_p=var(phenotype)
#var_e=var_p-var_a
#lamda=lamda(var_e,var_a)


### MATRIX MANIPULATION ###

#xpy<-as.matrix(t_mult_dupl(x,y))
#zpy<-as.matrix(t_mult_dupl(z,y))
X<-cbind(x,a)
#X<-z
neq=ncol(X)
nmeans=ncol(x)

e=y

#generate vector of starting values (sub-class means) for fixed and random effects 
sum<-rep(0,nmeans)
mean<-rep(0,nmeans)
ncols<-ncol(x)

for (i in 1:nmeans){
  sum<-x[,i]*y
  count<-length((which(sum!=0)))
  mean[i]<-colSums((sum)/count)
}


sol_random=rep(0,neq)
old_sol_random=rep(0,neq)
j=rep(0,neq)
start<-min(which(y!=0))
stop<-max(which(y!=0))

for (i in start:stop){
  j[i]=mean[which(x[i,]==1)]
  sol_random[i]=y[i]-j[i]
}

sol_random[1:nmeans]<-mean
CONV=1
CONV_TOL=0.00000000000001

ssqd_on=as.matrix(rep(0,neq))
ssq_n=as.matrix(rep(0,neq))
xpx=rep(0,neq)


for (j in 1:neq){
  xpx[j]<-as.matrix(t_mult_dupl(X[,j],X[,j]))  
}
count=1


##################################

### FIT SNPS AS RANDOM EFFECTS ###

##################################

repeat{
  for (j in start:neq){
    if (j==nmeans){
      #form lhs
      lhs=xpx[j]
    }
    else{
      #form lhs
      lhs=xpx[j]
    }
    #form rhs with y corrected by other effects (formula 1)
    rhs=(X[j]%*%e[j]) + (xpx[j]*sol_random[j])
    #do Gauss Seidel
    val=rhs/lhs
    #MCMC sample solution from its conditional (commented out here)
    #val_num=(rhs/lhs)
    #val_den=(1/lhs)
    #val=rnorm(1,(rhs/lhs),(1/lhs))
    #update e with current estimate (formula 2)
    e[j] = e[j] - (X[j]*(val-sol_random[j]))
    #update sol
    sol_random[j]=val
    ssqd_on[j]<-ssqd(old_sol_random[j],sol_random[j])
    ssq_n[j]<-sum(sol_random[j]^2)
    CONV<-ssqd_on[j]/ssq_n[j]
  }
  if (CONV<=CONV_TOL) break
  count=round(count+1/neq,1)
  print(count)
}

print(paste0("SNP EFFECTS (RANDOM) took ",count," iterations to solve",sep=" "))

sim_effects=c(0.5,-1,1,-0.5)
plot(sol_random[2:5],sim_effects,main="SNP Effect Estimates vs True Values")

######################################################################################

