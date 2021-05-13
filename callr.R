library(glmnet)
library(Matrix)

"dist2" = function( x, c = NA ) {
  
  # set the parameters for x
  if(is.na(c)) {
    c = x
  }
  
  # compute the dimension
  n1 = nrow(x)
  d1 = ncol(x)
  n2 = nrow(c)
  d2 = ncol(c)
  if(d1!=d2) {
    stop("Data dimension does not match dimension of centres.")
  }
  
  # compute the distance
  dist = t(rep(1,n2) %*% t(apply(t(x^2),MARGIN=2,FUN=sum))) + 
    (rep(1,n1) %*% t(apply(t(c^2),MARGIN=2,FUN=sum))) - 
    2 * (x%*%t(c))
  
  return(dist)
  
}

#the function to project a vector y to simplex
"ps" = function(y){
  
  
  z=sort(y)
  n=length(y)
  t=seq(n)
  for(i in 1:n-1){
    t[i]=(sum(z[(i+1):n])-1)/(n-i)
  }
  t[n]=z[n]-1
  tt=n-length(which(t<z))
  
  if(tt==0){
    tt=(sum(y)-1)/n
  }else{
    tt=t[tt]
  }
  
  x=y-tt
  x[x<0]<-0
  return(x)
}

#this is the function to get a training set matrix T from a known true label vector x
#r is the ratio of training set, can be set around 5%

"rand" = function(x,r){
  
  aaa=c()
  for (i in 1:length(unique(x))){
    aaa=c(aaa,sample(which(x==sort(unique(x))[i]),
                     max(ceiling(length(which(x==sort(unique(x))[i]))*r),2),replace=FALSE))
  }
  F=matrix(0,nrow=2,ncol=length(aaa))
  F[1,]=aaa
  F[2,]=x[aaa]
  return(F)
}



#input: expression matrix X, parameter u, training set matrix T

"callr"=function(X,u,T){
  
  N<-dim(X)[1]
  K=max(T[2,])
  
  Xt<-X[T[1,],]
  
  #do logistic regression on all samples to obtain the initial label matrix U
  logistic<-glmnet::glmnet(Xt,T[2,],family=c("multinomial"))
  
  b=matrix(0,nrow=K,ncol=dim(X)[2])
  for(i in 1:K){
    b[i,]=logistic[["beta"]][[i]][,length(logistic[["df"]])]
  }
  a=logistic[["a0"]][,length(logistic[["df"]])]
  
  U_hat=exp(t(t(X%*%t(b))+a))
  U_hat=U_hat/apply(U_hat,1,sum)
  U_hat=log(U_hat)
  
  for( i in 1:dim(T)[2]){
    U_hat[T[1,i],] = t(rep(-10000,K))
    U_hat[T[1,i],T[2,i]] = 0
  }
  U_hat[is.nan(U_hat)]=-10000
  U_hat[is.infinite(U_hat)]=-10000
  
  
  #Compute Laplacian Matrixx L
  
  
  
  D=dist2(X)
  
  D_sort = t(apply(D,MARGIN=2,FUN=sort));
  n = 17;
  
  
  mu=apply(D_sort[,2:n+1],MARGIN=1,FUN=mean)
  e=(matrix(rep(mu,N),N,N)+t(matrix(rep(mu,N),N,N)))/2
  
  W=exp(-D/(2*e^2))/(sqrt(2*pi)*e);
  
  
  for(i in 1:N){
    for(j in 1:N){
      if(D[i,j]>D_sort[i,n+1]&&D[i,j]>D_sort[j,n+1]){
        W[i,j] = 0;
      }else{
      }
    }
  }
  
  Degree = diag((rowSums(W))^(-0.5));
  L=diag(seq(1,1,length=N))-(Degree)%*%W%*%(Degree)
  
  
  #Some defaults and starting value
  
  dt = 0.005
  NS = 3
  
  
  U = matrix(data = runif(N*K,min=0,max=1),nrow=N,ncol=K)
  for(i in 1:N){
    U[i,]=ps(U[i,])
  }
  
  
  U[T[1,],] = U_hat[T[1,],]
  
  
  
  U0 = matrix(data = 0,nrow=N,ncol=K)
  U_old = matrix(data = 0,nrow=N,ncol=K)
  
  #The main iteration
  
  while (max(abs(U_old-U))>0.00001){
    U_old=U
    
    #MBO algorithm
    while(max(abs(U0-U))>0.00001){
      U0=U
      
      for(i in 1:NS){
        U = U+(dt/NS)*(-L%*%U+u*U_hat)
      }
      
      for(i in 1:N){
        U[i,] = ps(U[i,])
        
      }
    }
    #end of MBO algorithm
    
    for(i in 1:N){
      U[i,which.max(U[i,])]=1
      U[i,][U[i,]<1]=0
    }
    
    for( i in 1:dim(T)[2]){
      U[T[1,i],] = t(seq(0,0,length=K))
      U[T[1,i],T[2,i]] = 1
    }
    
    
    lab=seq(N)
    for(i in 1:N){
      lab[i]=which(U[i,]==1)
    }
    
    logistic<-glmnet::glmnet(X,lab,family=c("multinomial"))
    b=matrix(0,nrow=K,ncol=dim(X)[2])
    for(i in 1:K){
      b[i,]=logistic[["beta"]][[i]][,length(logistic[["df"]])]
    }
    a=logistic[["a0"]][,length(logistic[["df"]])]
    
    U_hat=exp(t(t(X%*%t(b))+a))
    U_hat=U_hat/apply(U_hat,1,sum)
    
    U_hat=log(U_hat)
    
    for( i in 1:dim(T)[2]){
      U_hat[T[1,i],] = t(seq(-10000,-10000,length=K))
      U_hat[T[1,i],T[2,i]] = 0
    }
    U_hat[is.nan(U_hat)]=-10000
    U_hat[is.infinite(U_hat)]=-10000
    
  }
  #end of the main iteration
  
  
  return(lab)
  #return the estimated labels
}