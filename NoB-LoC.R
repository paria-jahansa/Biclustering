library(invgamma)
C_full=list()
w_full=list()

save(C_full ,file= "c50.rda")
save(w_full ,file= "w50.rda")
dddd=load(file = "c50.rda")
attach(dddd)
for(l in 49:50){


Y=simulated_data
dim(Y)
#####################################
coloumn.clusters = hclust(dist(t(Y)))
w_initial=as.vector(cutree(coloumn.clusters,k=11))
singleton_clusters_w=as.vector(which(table(w_initial)==1))
for(i in 1:length(w_initial)){
  for(j in 1:length(singleton_clusters_w)){
    if(w_initial[i]==singleton_clusters_w[j]){
      w_initial[i]=0
    }
  }
}

table(w_initial)
#w_initial=rep(c(1,2,0),c(8,4,8))#fixed
##########################
#initial value of c

C_initial=numeric()
for(i in unique(w_initial)){
  if(i!=0){
    y=Y[,which(w_initial==i)]
    row.clusters = hclust(dist(y))
    c_initial=as.vector(cutree(row.clusters,k=40))
    singleton_clusters_c=as.vector(which(table(c_initial)==1))
    for(z in 1:length(c_initial)){
      for(j in 1:length(singleton_clusters_c)){
        if(c_initial[z]==singleton_clusters_c[j]){
          c_initial[z]=0
        }
      }
    }
    C_initial=rbind(C_initial,c(i,c_initial))
  }
}
dim(C_initial)

C_initial[,1]
#####################################
#hyper parameters for mu

G=20
N=120
mu=numeric()
tau2_mean=numeric()
for(g in 1:G){
  mu[g]=median(Y[,g])
  tau2_mean[g]=var(Y[,g])
}
sigma2_mean=1

tau2_var=.4
sigma2_var=.4
tau2_hy_alpha=(tau2_var/((tau2_mean)^2))+2
tau2_hy_beta=.3*(tau2_hy_alpha-1)
sigma2_hy_alpha=(sigma2_var/((sigma2_mean)^2))+2
sigma2_hy_beta=.3*(sigma2_hy_alpha-1)
#pi values
pi_0=.6
pi_1=.8



###########################
#initial value of theta

c=numeric()
theta_initial=matrix(0,N,G)
g=1


for(g in 1:G){
  if(w_initial[g]>0){
    c=C_initial[max(which(C_initial[,1]==w_initial[g])),2:length(C_initial[1,])]
    for(i in 1:N){
      if(c[i]>0){
        mu_normal=((1/tau2_mean[g]+(sum(c==c[i]))/sigma2_mean)^(-1))*
          (mu[g]/tau2_mean[g]+(sum(Y[,g][c==c[i]]))/sigma2_mean)
        sigma2_normal=((1/tau2_mean[g]+(sum(c==c[i]))/sigma2_mean)^(-1))
        ##print(sigma2_normal)
        theta_initial[which(c==c[i]),g]=rnorm(1,mu_normal,sigma2_normal)
      }else{
        mu_normal=((1/tau2_mean[g]+1/sigma2_mean)^(-1))*
          ( mu[g]/tau2_mean[g]+Y[i,g]/sigma2_mean)
        sigma2_normal=((1/tau2_mean[g]+1/sigma2_mean)^(-1))
        theta_initial[i,g]=rnorm(1,mu_normal,sigma2_normal)
      }
    }
  }else{
    for(i in 1:N){
      mu_normal=((1/tau2_mean[g]+1/sigma2_mean)^(-1))*
        (mu[g]/tau2_mean[g]+Y[i,g]/sigma2_mean)
      sigma2_normal=((1/tau2_mean[g]+1/sigma2_mean)^(-1))
      theta_initial[i,g]=rnorm(1,mu_normal,sigma2_normal)
    }
  }
}


dim(theta_initial)
###############################################
MM=20000
M=0.1;alpha=10
thetas=list()
thetas[[1]]=theta_initial


sigma=matrix(0,MM,G)
sigma[1,]=rep(1,G)

tau0=tau1=tau2=matrix(0,MM,G)
tau0[1,]=tau1[1,]=tau2[1,]=rep(1,G)



c=list()
c[[1]]=C_initial


w=list()
w[[1]]=w_initial
###############################################

for(r in 2:MM){ 
  
  #################### thetas ###########################
  thetas[[r]]=matrix(0,N,G)
  for(g in 1:G){ 
    #print(g)
    if(w[[r-1]][g]==0){
      for(i in 1:N){
        v1=(1/tau2[r-1,g]+1/sigma[r-1,g])^(-1)
        mu1=v1*(mu[g]/tau2[r-1,g]+Y[i,g]/sigma[r-1,g])
        thetas[[r]][i,g]=rnorm(1,mu1,sqrt(v1))
        ##print(thetas[[r]][i,g])
      }}else{ 
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        index=which(w[[r-1]][g]==c[[r-1]][,1])
        for(i in 1:N){
          #print(c[[r-1]][index,-1][i])
          if(c[[r-1]][index,-1][i]>0){
            index1=c(1:N)[c[[r-1]][index,-1]==c[[r-1]][index,-1][i]]
            n2=length(index1)
            v2=(1/tau0[r-1,g]+n2/sigma[r-1,g])^(-1)
            mu2=v2*(mu[g]/tau0[r-1,g]+sum(Y[index1,g])/sigma[r-1,g])
            thetas[[r]][index1,g]=rnorm(1,mu2,sqrt(v2))
            
          }else{
            v3=(1/tau1[r-1,g]+1/sigma[r-1,g])^(-1)
            mu3=v3*(mu[g]/tau1[r-1,g]+Y[i,g]/sigma[r-1,g])
            thetas[[r]][i,g]=rnorm(1,mu3,sqrt(v3))
          }
        }     
      }}
  
  ########################## sigma ###########################
  for(g in 1:G){
    sigma[r,g]= rinvgamma(1, shape=sigma2_hy_alpha+N/2, 
                          rate =sigma2_hy_beta+sum((Y[,g]- thetas[[r]][,g])^2)/2   )
  }
  
  ########################## tau ###########################
  
  for(g in 1:G){
    #print(g)
    #active
    index=which(w[[r-1]][g]==c[[r-1]][,1])# th
    
    indexcluster=unique(c[[r-1]][index,-1])#labels of samples' cluster
    if(0 %in% indexcluster){Ks=length(indexcluster)-1}else{Ks=length(indexcluster)}
    sample=c[[r-1]][index,-1]
    ms=sum(sample!=0)#*sum(w[[r-1]]==w[[r-1]][g])
    
    indexKs=which(c[[r-1]][index,-1]!=0)
    tau0a= rinvgamma(1, shape=tau2_hy_alpha+Ks/2, rate =tau2_hy_beta+sum((unique(thetas[[r]][indexKs,g])-mu[g])^2)/2)#?
    
    tau1a= rinvgamma(1, shape=tau2_hy_alpha+(N-ms)/2, 
                     rate =tau2_hy_beta+sum((thetas[[r]][c(1:N)[c[[r-1]][index,-1]==0],g]-mu[g])^2)/2  )
    
    tau2a= rinvgamma(1, shape=tau2_hy_alpha, rate =tau2_hy_beta)
    
    ########## inactive
    
    tau0i= rinvgamma(1, shape=tau2_hy_alpha, rate =tau2_hy_beta)
    
    tau1i= rinvgamma(1, shape=tau2_hy_alpha, rate =tau2_hy_beta)
    
    
    
    tau2i= rinvgamma(1, shape=tau2_hy_alpha+N/2, rate =tau2_hy_beta+sum((thetas[[r]][,g]-mu[g])^2)/2)
    
    #####################
    
    tau0[r,g]= tau0a*sum(w[[r-1]][g]>0)+ tau0i*sum(w[[r-1]][g]==0)
    tau1[r,g]= tau1a*sum(w[[r-1]][g]>0)+ tau1i*sum(w[[r-1]][g]==0)
    tau2[r,g]= tau2a*sum(w[[r-1]][g]>0)+ tau2i*sum(w[[r-1]][g]==0)
    
    #print(g)
  }
  
  ################################ C 
  Cnew2=c()
  
  for(s in 1:length(c[[r-1]][,1])){
    #print(s)
    Cnew=c[[r-1]][s,-1]
    for(i in 1:N){ 
      #print(i)
      
      indexcluster=unique(Cnew)#labels of samples' cluster
      if(0 %in% indexcluster){Ksnew=length(indexcluster)}else{Ksnew=length(indexcluster)+1}
      
      prob=c()
      indexg=1:G
      Vs=indexg[w[[r-1]]==c[[r-1]][s,1]]#cluster sth of protein 
      
      #prob[e]=(1-pi_1)*prod(dnorm(Y[i,Vs],mu[Vs],sqrt(tau1[r,Vs]+sigma[r,Vs])))  #=0
      ############################################################################################
      indexN=1:N
      Rs=list()
      e=0
      for(k in 1:Ksnew){
        e=e+1
        Cnew1=Cnew
        Cnew1[i]=NA
        
        Rs[[k]]=c(indexcluster[k],length(indexN[Cnew1==indexcluster[k]][is.na(indexN[Cnew1==indexcluster[k]])==FALSE]),
                  indexN[Cnew1==indexcluster[k]][is.na(indexN[Cnew1==indexcluster[k]])==FALSE])# the first one is the label of C, the second is the length of Rsk-=nsk-, the others are Rsk-
        
        
        ms=sum(Cnew[-i]!=0)#ms-
        
        if(indexcluster[k]==0){prob[e]=(1-pi_1)*prod(dnorm(Y[i,Vs],mu[Vs],sqrt(tau1[r,Vs]+sigma[r,Vs])))  #=0
        }else{ 
          if((length(Rs[[k]][-c(1,2)])!=1) & (length(Rs[[k]][-c(1,2)])!=0) & (length(Vs)!=1)){prob[e]=pi_1*Rs[[k]][2]/(M+ms)*
            prod(sqrt((1/tau0[r,Vs]+Rs[[k]][2]/sigma[r,Vs])/
                        ((1/tau0[r,Vs]+(Rs[[k]][2]+1)/sigma[r,Vs])*(2*pi*sigma[r,Vs])))*
                   exp(-Y[i,Vs]^2/(2*sigma[r,Vs])-.5*((1/tau0[r,Vs]+Rs[[k]][2]/sigma[r,Vs]))^(-1)*
                         (mu[Vs]/tau0[r,Vs]+apply(Y[Rs[[k]][-c(1,2)],Vs],2,sum)/sigma[r,Vs])^2)*
                   exp(0.5*(1/tau0[r,Vs]+(Rs[[k]][2]+1)/sigma[r,Vs])^(-1)*
                         (mu[Vs]/tau0[r,Vs]+apply(Y[union(i,Rs[[k]][-c(1:2)]),Vs],2,sum)/sigma[r,Vs])^2))
          
          }
          if((length(Rs[[k]][-c(1,2)])==1) | (length(Vs)==1) ){prob[e]=pi_1*Rs[[k]][2]/(M+ms)*
            prod(sqrt((1/tau0[r,Vs]+Rs[[k]][2]/sigma[r,Vs])/
                        ((1/tau0[r,Vs]+(Rs[[k]][2]+1)/sigma[r,Vs])*(2*pi*sigma[r,Vs])))*
                   exp(-Y[i,Vs]^2/(2*sigma[r,Vs])-.5*((1/tau0[r,Vs]+Rs[[k]][2]/sigma[r,Vs]))^(-1)*
                         (mu[Vs]/tau0[r,Vs]+sum(Y[Rs[[k]][-c(1,2)],Vs])/sigma[r,Vs])^2)*
                   exp(0.5*(1/tau0[r,Vs]+(Rs[[k]][2]+1)/sigma[r,Vs])^(-1)*
                         (mu[Vs]/tau0[r,Vs]+sum(Y[union(i,Rs[[k]][-c(1:2)]),Vs])/sigma[r,Vs])^2))
          
          }
          if(length(Rs[[k]][-c(1,2)])==0){
            prob[e]=0
          }
        }
        ##print(k)
      }
      
      
      prob[Ksnew+1]=pi_1*M/(M+ms)*
        prod(dnorm(Y[i,Vs],mu[Vs],sqrt(tau0[r,Vs]+sigma[r,Vs]))) 
      
      prob[which(is.infinite(prob))]=1
      prob[which(is.na(prob))]=0
      Csample=cbind(c(indexcluster,max(indexcluster)+1),prob/sum(prob))
      Cnew[i]=sample(Csample[,1],size=1,prob=Csample[,2])
      #print(i)
    }
    
    Cnew2=rbind(Cnew2,c(c[[r-1]][s,1],Cnew))
    
  }
  c[[r]]=Cnew2
  #c[[r-1]][s,-1]
  
  
  ############################# w ################ 
  ################################################
  Cupdate=c()
  Wnew=c()
  
  wnew=w[[r-1]]
  for(g in 1:G){ 
    #print(g)
    cnew2=c()
  
    wnew[g]=NA
    
    indexG=1:G
    ps=length(indexG[wnew==w[[r-1]][g]][is.na(indexG[wnew==w[[r-1]][g]])==FALSE])
    Gm=length(indexG[wnew!=0][is.na(indexG[wnew!=0])==FALSE])
    
    prob1=c()
    e1=1
    indexclusterprotein=unique(wnew[is.na(wnew)==FALSE])#labels of samples' cluster
    if(0 %in% indexclusterprotein){S=length(indexclusterprotein)}else{S=length(indexclusterprotein)+1}
    prob1[e1]=(1-pi_0)*prod(dnorm(Y[,g],mu[g],sqrt(tau2[r,g]+sigma[r,g])))  #=0
    for(z in 2:S){ 
      e1=e1+1
      cnew2=c[[r]][which(c[[r]][,1]==sort(unique(w[[r-1]]))[z]),-1]
      ms1=length(cnew2[cnew2!=0])
      
      
      indexcluster2=unique(cnew2)  
      if(0 %in% indexcluster2){Ksnew2=length(indexcluster2)}else{Ksnew2=length(indexcluster2)+1}
      indexcluster2=sort(indexcluster2)
      Rs2=list()
      tt=nsk=c()
      for(k in 1:Ksnew2){ 
        Rs2[[k]]=c(indexcluster2[k],length(indexN[cnew2==indexcluster2[k]][is.na(indexN[cnew2==indexcluster2[k]])==FALSE]),
                   indexN[cnew2==indexcluster2[k]][is.na(indexN[cnew2==indexcluster2[k]])==FALSE])# the first one is the label of C, the second is the length of Rsk-=nsk-, the others are Rsk-
        
        if(k!=1){
          tt[k]= sqrt(2*pi*(1/tau0[r,g]+Rs2[[k]][2]/sigma[r,g])^(-1))/((sqrt(2*pi*sigma[r,g])^Rs2[[k]][2])*sqrt(2*pi*tau0[r,g]))*
            exp(-mu[g]^2/(2*tau0[r,g])-sum(Y[Rs2[[k]][-c(1,2)],g]^2)/(2*sigma[r,g])+
                  0.5*(1/tau0[r,g]+Rs2[[k]][2]/sigma[r,g])^(-1)*(mu[g]/tau0[r,g]+sum(Y[Rs2[[k]][-c(1,2)],g])/sigma[r,g])^2)   
        }else{
          tt[k]=prod(dnorm(Y[Rs2[[1]][-c(1,2)],g],mu[g],sqrt(tau1[r,g]+sigma[r,g])))
        }
        
        
        nsk[k]=Rs2[[k]][2]
        
        
      }

        prob1[e1]=pi_0*ps/(alpha+Gm)*prod(tt[2:Ksnew2])*tt[1]
      } 
    
    cnew2=c()
    coloumn.clusters = hclust(dist(Y[,g]))
    Csingle=as.vector(cutree(coloumn.clusters,k=40))
    singleton_clusters_w=as.vector(which(table(Csingle)==1))
    for(i in 1:length(Csingle)){
      for(j in 1:length(singleton_clusters_w)){
        if(Csingle[i]==singleton_clusters_w[j]){
          Csingle[i]=0
        }
      }
    }
    cnew2=Csingle
    ms1=length(cnew2[cnew2!=0])
    
    
    indexcluster2=unique(cnew2)  
    if(0 %in% indexcluster2){Ksnew2=length(indexcluster2)}else{Ksnew2=length(indexcluster2)+1}
    indexcluster2=sort(indexcluster2)
    Rs2=list()
    tt=nsk=c()
    for(k in 1:Ksnew2){ 
      Rs2[[k]]=c(indexcluster2[k],length(indexN[cnew2==indexcluster2[k]][is.na(indexN[cnew2==indexcluster2[k]])==FALSE]),
                 indexN[cnew2==indexcluster2[k]][is.na(indexN[cnew2==indexcluster2[k]])==FALSE])# the first one is the label of C, the second is the length of Rsk-=nsk-, the others are Rsk-
      
      if(k!=1){
        tt[k]= sqrt(2*pi*(1/tau0[r,g]+Rs2[[k]][2]/sigma[r,g])^(-1))/((sqrt(2*pi*sigma[r,g])^Rs2[[k]][2])*sqrt(2*pi*tau0[r,g]))*
          exp(-mu[g]^2/(2*tau0[r,g])-sum(Y[Rs2[[k]][-c(1,2)],g]^2)/(2*sigma[r,g])+
                0.5*(1/tau0[r,g]+Rs2[[k]][2]/sigma[r,g])^(-1)*(mu[g]/tau0[r,g]+sum(Y[Rs2[[k]][-c(1,2)],g])/sigma[r,g])^2)   
      }else{
        tt[k]=prod(dnorm(Y[Rs2[[1]][-c(1,2)],g],mu[g],sqrt(tau1[r,g]+sigma[r,g])))
      }
      
      nsk[k]=Rs2[[k]][2]
    }
    if(0 %in% indexcluster2){Ksnew3=length(indexcluster2)-1}else{Ksnew3=length(indexcluster2)}
    
    prob1[length(prob1)+1]=pi_0*alpha/(alpha+Gm)*pi_1^ms1*(1-pi_1)^(N-ms1)*(M^Ksnew3*prod(gamma(nsk[-1]))/prod(M+(1:ms1)-1))*
      prod(tt[2:Ksnew2])*tt[1]
    
    prob1[which(is.infinite(prob1))]=1
    prob1[which(is.na(prob1))]=0
    indexw=unique(c(0,unique(wnew[is.na(wnew)==FALSE]),max(unique(wnew[is.na(wnew)==FALSE]))+1 ))
    Wsample=cbind(indexw,prob1/sum(prob1))
    wnew[g]<-sample(Wsample[,1],size=1,prob=Wsample[,2])
    #print(wnew[g])
    if(wnew[g]!=0){
      #print(wnew[g])
      Cnew3=c()
      
      if((wnew[g] %in% c[[r]][,1]))(Cnew3=c[[r]][which(c[[r]][,1]==wnew[g]),])
      if(((wnew[g] %in% c[[r]][,1])==FALSE))(Cnew3=c(wnew[g],Cnew2))
      
      
      #print(Cnew3)
      Cupdate=rbind(Cupdate,Cnew3)
      #print(dim(Cupdate))
    }
    
    
    #print(Cupdate)
  }
  #edit(Cupdate)
  #print(w[[r-1]])
  #print(wnew)
  w[[r]]=wnew
  
  #AAA= setdiff(intersect(unique(wnew),unique(Cupdate[,1])),{0})
  
  c[[r]]=unique(Cupdate)
  
  
  print(table(w[[r]]))
  #print(c[[r]])
  print(r)
}



C_full[[l]]=c
w_full[[l]]=w
print("#########################################")
print(l)
}


########################### Similarity for C ############ 


burnin=floor(MM/2)
#burnin=7500
indexn=c(burnin:20000)

BB=list()
j=0
cc=C_full[[1]]
for(k in indexn){
  j=j+1
  BB[[j]]=c[[k]][,-1]
}
M=length(indexn)
DIS_1=array(0,c(N,N,M))
DIS_2=array(0,c(N,N,M))

k=1
i=1
j=1
for(k in 1:M){
  print(k)
  for(i in 1:N){
    for(j in 1:N){
      DIS_1[i,j,k]<-as.numeric(identical(BB[[k]][1,i], BB[[k]][1,j]))
      DIS_2[i,j,k]<-as.numeric(identical(BB[[k]][2,i], BB[[k]][2,j]))
    }
  }
}
PSM_1=apply(DIS_1,c(1,2),mean)
PSM_2=apply(DIS_2,c(1,2),mean)

m1=minbinder(PSM_1)
m2=minbinder(PSM_2)
table(c_1,m1$cl)
table(c_2,m2$cl)
table(c_1)

c1_sim=numeric(120)
c2_sim=numeric(120)
c2_sim[which(m2$cl==30)]=28


c1_sim[which(c1_sim==24)]=30
c2_sim[which(c2_sim==10)]=34
table(c_1,c1_sim)
table(c_2,c2_sim)


adj.rand.index(c_1,c1_sim)
adj.rand.index(c_2,c2_sim)

col=redgreen(64)


#protein set 0
hh1=Heatmap(simulated_data[,13:20], cluster_rows = FALSE,cluster_columns = FALSE, 
            name="protein set 0",col=col,
            column_title = "protein set 0",show_row_names = FALSE)

#protein set 1
data2=cbind(simulated_data[,1:8],c1_sim)
head(data2)
data2=data2[order(data2[,9]),]
colnames(data2)=1:9
rownames(data2)=1:120

hh2=Heatmap(data2[,1:8], cluster_rows = FALSE,cluster_columns = FALSE
            ,name="protein set 1",col=col,
            column_title = "protein set 1",show_row_names = FALSE)



#Protein set 2
data3=cbind(simulated_data[,9:12],c2_sim)
head(data3)
data3=data3[order(data3[,5]),]
colnames(data3)=9:13
rownames(data3)=1:120

hh3=Heatmap(data3[,1:4], cluster_rows = FALSE,cluster_columns = FALSE, 
            name="protein set 2",col=col,
            column_title = "protein set 2",show_row_names = FALSE)


hh=hh2+hh3+hh1
draw(hh)
s