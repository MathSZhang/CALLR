#This is the example marker file in which we don't distinguish CD14+ Monocytes and FCGR3A+ Monocytes,
#we also don't distinguish CD4 T cells and CD8 T cells in this example.
#If you want to input more cell types and their marker genes, just write them in this format.



marker_file=t(matrix(c( 'CD34', 'THY1', 'ENG', 'KIT', 'PROM1', 
                        'NCAM1', 'CD56','CD16','CD3','NKG7',
                        'CD14', 'FCGR1A', 'CD68', 'S100A12','',
                        'CD19', 'MS4A1', 'CD79A','','',
                        'CD8A', 'CD8B','CD3D','CD3E','CD4',
                        'IL3RA',  'MHC class II', 'CLEC4C', 'NRP1', 'CD83',
                        'PF4','','','',''),5,7))
#in this case the maximum number of marker genes in one cell type is 5,
#and the total number of labels is 7, so the matrix is 5*7

rownames(marker_file)=c("CD34+","NK cells","Monocytes","B cells",
                        "T cells","Dendritic cells","Megakaryocytes")
View(marker_file)


#This is the function for data preprocessing 
preprocess=function(X){
  zs=which(apply(X,1,sum)==0)
  X=X[-zs,]
  
  m=dim(X)[2]
  for(i in 1:m){
    a=X[,i]
    a=as.matrix(a)[-which(a==0),]
    X[,i]=X[,i]/exp(mean(log(a)))
  }
  X=as.matrix(X)
  return(X)
}


#This is the function to compute marker gene score and choose representative cells
representative=function(marker_file,X,cutoff){
  n=dim(X)[1]
  m=dim(X)[2]
  #calculate TF-IDF transformation matrix
  rs=as.matrix(apply(X,1,sum))
  cs=as.matrix(apply(X,2,sum))
  T=X*log(1+m/rs)%*%t(1/cs)
  T[T<t(apply(T, 1, sort))[,ceiling(m*0.95)]]=0
  
  #calculate aggregated marker score
  S=matrix(0,dim(marker_file)[1],m)
  for(k in 1:dim(marker_file)[1]){
    for(i in 1:dim(marker_file)[2]){
      a=which(rownames(X)==marker_file[k,i])
      if (length(a)==0){
      } else {
        for(j in 1:m){
          S[k,j]=S[k,j]+T[a,j]
        }
      }
    }
  }

  st=sort(matrix(S,nrow=1))
  score=st[length(which(st==0))+
             ceiling(cutoff*(dim(marker_file)[1]*m-length(which(st==0))))]
  
  #choose cells in the cutoff percentile and above  in only one cell type 
  chosen_cells=rep(0,m)
  for(i in 1:m){
    if(length(S[,i][S[,i]<score])==dim(marker_file)[1]-1){
      chosen_cells[i]=which(S[,i]>=score)
    }else{
    }
  }
  
  #transform the label value to apply to CALLR function
  #the numeric label in CALLR must be continuous natural number from 1 to K, where K is the total number of labels
  a=as.numeric( names(table(chosen_cells)[-1]))
  correspondence=rep(0,max(a))
  correspondence[a]=seq(length(a))
  names(correspondence)=rownames(marker_file)
  

  F=matrix(0,2,length(chosen_cells[chosen_cells>0]))
  F[1,]=which(chosen_cells>0)
  F[2,]=correspondence[chosen_cells[which(chosen_cells>0)]]
  
  

  return(list(training_set=F,label=rownames(marker_file)[a]))
}



callr_marker=function(marker_file, X, u, cutoff){
  X=preprocess(X)
  rep=representative(marker_file,X,cutoff)
  F=rep[["training_set"]]
  output=callr(t(X),u,F)
  
  return(rep[["label"]][output])
}




#here is the command to load data and to do cluster analysis
load("sample_data.RData")
X=preprocess(X)
#use marker gene to choosed representative cells
rep=representative(marker_file,X,0.85)
#get the training set
T=rep[["training_set"]]
#semi-supervised learning
output=callr(t(X),0.3,T)

#predicted labels
View(rep[["label"]][output])
#True labels
View(type[,2])












