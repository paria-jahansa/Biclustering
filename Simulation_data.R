#########################################
#Protein Allocation
rm(list = ls())
G=20
N=120

w=rep(c(1,2,0),c(8,4,8))
W=rep(c(1,2,0),c(8,4,8))

#n_0=1
#n_1=1
#n_2=1
#set.seed(2)
#for(i in 1:20){
#  p_0=8/20
#  p_1=(1-(8/20))*(n_1/(.2+(n_0+n_1+n_2)))
#  p_2=(1-(8/20))*(n_2/(.2+(n_0+n_1+n_2)))
#  w[i]<-sample(0:2,size=1,prob=c(p_0,p_1,p_2))
#  if(w[i]==1){
#    n_1=n_1+1
#  }else if(w[i]==2){
#    n_2=n_2+1
#  }
#}

###########################################
#Sample Allocation

set.seed(1)
c_1=numeric(120)
c_2=numeric(120)
c_1=sample(0:3,size=120,replace = T)
c_2=sample(0:2,size=120,replace = T)


table(c_1)
table(c_2)

###########################################
#theta_stars 

#for active protein and active sample

library(MASS)
mu_mv=c(-.8,-.4,0,.4,.8)
sigma_mv=0.01*diag(5)
set.seed(3)

#Protein set 1

v=mvrnorm(n = 1, mu=mu_mv, Sigma=sigma_mv)
theta_star_Protein_set_1=sample(v,3)
for(i in 2:8){
  v=mvrnorm(n = 1, mu=mu_mv, Sigma=sigma_mv)
  theta_star_Protein_set_1=cbind(theta_star_Protein_set_1,sample(v,3))
}

#Protein set 2

v=mvrnorm(n = 1, mu=mu_mv, Sigma=sigma_mv)
theta_star_Protein_set_2=sample(v,2)
for(i in 2:4){
  v=mvrnorm(n = 1, mu=mu_mv, Sigma=sigma_mv)
  theta_star_Protein_set_2=cbind(theta_star_Protein_set_2,sample(v,2))
}

###########################################
#let's simulate data


active_sample=function(v){
  sigma=.3
  c=numeric()
  for(j in 1:length(v)){
    c[j]=rnorm(1,v[j],sigma)
  }
  return(c)
}

#
sigma=.3
Y_1=numeric()
theta_1=numeric()
for(i in 1:N){
  if(c_1[i]==1){
    Y_1=rbind(Y_1,active_sample(theta_star_Protein_set_1[1,]))
    theta_1=rbind(theta_1,theta_star_Protein_set_1[1,])
  }else if(c_1[i]==2){
    Y_1=rbind(Y_1,active_sample(theta_star_Protein_set_1[2,]))
    theta_1=rbind(theta_1,theta_star_Protein_set_1[2,])
  }else if(c_1[i]==3){
    Y_1=rbind(Y_1,active_sample(theta_star_Protein_set_1[3,]))
    theta_1=rbind(theta_1,theta_star_Protein_set_1[3,])
  }else{
    y=numeric()
    mu=numeric()
    for(j in 1:8){
      mu[j]=runif(1,-.8,.8)
      y[j]=rnorm(1,mu[j],sigma)
    }
    Y_1=rbind(Y_1,y)
    theta_1=rbind(theta_1,mu)
  }
}

dim(Y_1)
dim(theta_1)
#
theta_2=numeric()
Y_2=numeric()
for(i in 1:N){
  if(c_2[i]==1){
    Y_2=rbind(Y_2,active_sample(theta_star_Protein_set_2[1,]))
    theta_2=rbind(theta_2,theta_star_Protein_set_2[1,])
  }else if(c_2[i]==2){
    Y_2=rbind(Y_2,active_sample(theta_star_Protein_set_2[2,]))
    theta_2=rbind(theta_2,theta_star_Protein_set_2[2,])
  }else{
    y=numeric()
    mu=numeric()
    for(j in 1:4){
      mu[j]=runif(1,-.8,.8)
      y[j]=rnorm(1,mu[j],sigma)
    }
    Y_2=rbind(Y_2,y)
    theta_2=rbind(theta_2,mu)
  }
}

dim(Y_2)
dim(theta_2)



##
Y_3=numeric()
theta_3=numeric()

for(i in 1:8){
  mu=numeric()
  mu=runif(N,-.8,.8)
  Y_3=cbind(Y_3,active_sample(mu))
  theta_3=cbind(theta_3,mu)
}
dim(Y_3)
dim(theta_3)
##

simulated_data=cbind(Y_1,Y_2,Y_3)
theta_truth=cbind(theta_1,theta_2,theta_3)
dim(simulated_data)
head(simulated_data)
class(simulated_data)
colnames(simulated_data)=1:20
rownames(simulated_data)=1:120

colnames(theta_truth)=1:20
rownames(theta_truth)=1:120
#####################################################
#heatmap of simulated data


#install.packages("devtools")
#library(usethis)
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(grid)
library(ComplexHeatmap)

###
library(gplots)
col=redgreen(64)


#protein set 0
h_1=Heatmap(simulated_data[,13:20], cluster_rows = FALSE,cluster_columns = FALSE, 
            name="protein set 0",col=col,
            column_title = "protein set 0",show_row_names = FALSE)

#protein set 1
data_2=cbind(Y_1,c_1)
head(data_2)
data_2=data_2[order(data_2[,9]),]
colnames(data_2)=1:9
rownames(data_2)=1:120

h_2=Heatmap(data_2[,1:8], cluster_rows = FALSE,cluster_columns = FALSE
            ,name="protein set 1",col=col,
            column_title = "protein set 1",show_row_names = FALSE)



#Protein set 2
data_3=cbind(Y_2,c_2)
head(data_3)
data_3=data_3[order(data_3[,5]),]
colnames(data_3)=9:13
rownames(data_3)=1:120

h_3=Heatmap(data_3[,1:4], cluster_rows = FALSE,cluster_columns = FALSE, 
            name="protein set 2",col=col,
            column_title = "protein set 2",show_row_names = FALSE)


h=h_2+h_3+h_1
draw(h)

#############################################
#Hierarchical clustering

heatmap(simulated_data,col=col)

hist(simulated_data)
