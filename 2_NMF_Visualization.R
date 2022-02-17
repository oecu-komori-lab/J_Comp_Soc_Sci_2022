#Performing NMF and visualization
#Masashi Komori

#Load the formatted data
load(file="sub_adj_vec.RData")

#7.NMF ---------------------------------------------------------------------
library(NMF)
library(dplyr)

sub_adj_vec <- sub_adj_vec + 0.000001

#Examination of the num of factors
#res<-nmf(sub_adj_vec,rank=c(6,7,8,9,10),seed=12345)
#plot(res)
#save(res, file = "res2-12.rda")

num_factor <-5 
nmf_res <- nmf(sub_adj_vec, num_factor, seed=12345)  #NMF


# 8.Transforming coef into a social network -------------------------------

gender_age <- read.csv("対応表/所持者データ.csv")# File containing data of device num, gender, and age etc.
Cname <- colnames(adj_array_temp)
Rname <- rownames(adj_array_temp)


par(mfrow=c(1,num_factor))
library(igraph)
for(i in 1 :num_factor){
  coef_adj_mat <- matrix(0, nrow=length(Rname), ncol=length(Cname))# 
  dimnames(coef_adj_mat) <- list(Cname,Rname)
  for(j in 1: ncol(coef(nmf_res))){
    coef_adj_mat[fromname[j],toname[j]] <-  coef(nmf_res)[i,j]
  }

  
  ##Reorder by age
  age <- rep(NA, nrow(coef_adj_mat))
  for(j in 1 : nrow(coef_adj_mat)){
    age[j] <- gender_age[gender_age[,1] == rownames(coef_adj_mat)[j],5]
  }
  coef_adj_mat <- coef_adj_mat[order(age,decreasing=T),order(age,decreasing=T)]

  #coef into adjacent matrix
  g.in<-graph.adjacency(coef_adj_mat*(1000),weighted=TRUE, mode="undirected")
  
  #centrality
  g.deg<-degree(g.in)
  g.bet<-betweenness(g.in)
  g.eigen <- evcent(g.in)$vector
  
  #clustering
  fc.in <- fastgreedy.community(g.in)
  mo.in <- multilevel.community(g.in)
  ed.in <- edge.betweenness.community(g.in)
  
  age <- gender <- gendermark <- rep(NA,  length(V(g.in)$name))
  for(j in 1 : length(V(g.in)$name)){
    age[j] <- gender_age[gender_age[,1] == V(g.in)$name[j],5]
    if(gender_age[gender_age[,1] == V(g.in)$name[j],4] == "女"){
      gendermark[j] <- "circle"
      gender[j] <- 0
    }else{
      gendermark[j] <- "square"
      gender[j] <- 1
    }
  }

  V(g.in)$name <-  1:length(rownames(coef_adj_mat))
  
  par(mfrow=c(1,1))
  agef <- 1-(floor(age/10) /10)
  plot(g.in,vertex.size=10, vertex.label.font=1,vertex.label.cex=0.5,edge.width=E(g.in)$weight/10,
       layout=layout.circle,main=paste("Factor",i),vertex.color=fc.in$membership)
  
  options(X11fonts = c("-*-gothic-%s-%s-normal--%d-*-*-*-*-*-*-*","-adobe-symbol-*-*-*-*-%d-*-*-*-*-*-*-*"))
  dev.copy2eps(file=paste("log_factor",i,"of",num_factor,".eps",sep = ""),family = "serif",width = 6)
  indivchara <- data.frame(ID=rownames(coef_adj_mat), DEGREE=g.deg, BETWEEN=g.bet, EIGEN=g.eigen, clust_fc=fc.in$membership, clust_mo=mo.in$membership, clust_ed=ed.in$membership, AGE =age, GENDER = gender) 
  write.csv(indivchara,  paste("log_factor",i,"of",num_factor,"centrality.csv", sep=""))
}

# 9.Figure ----------------------------------------------------------------
#Making a diagram of the occurrence of social contact
basis_restore <- matrix(0, nrow=length(adj_array[1,1,]), ncol=ncol(basis(nmf_res)))
rownames(basis_restore) <- as.POSIXct(as.numeric(dimnames(adj_array)[[3]]),origin="1970-01-01","Japan")
for(i in 1 : nrow(basis(nmf_res))){
  basis_restore[time[i],] <- basis(nmf_res)[i,]
}

day <- format(as.POSIXct(as.numeric(dimnames(adj_array)[[3]]),origin="1970-01-01","Japan"), "%m/%d")
#par(mfrow=c(num_factor,1))

par(mfrow=c(1,1))
for(i in 1: num_factor){
  options(X11fonts = c("-*-gothic-%s-%s-normal--%d-*-*-*-*-*-*-*","-adobe-symbol-*-*-*-*-%d-*-*-*-*-*-*-*"))
  plot(basis_restore[,i], type="l",xaxt= "n", xlab="DAY",ylab="SOCIAL CONTACT",main=paste("Factor", i), xaxp  = c(1, length(day), 48), cex.main = 1.8, cex.lab  = 1.5)
  par(xaxt = "s")
  axis(side=1, at=seq(1,length(day),by =48),labels =day[seq(1,length(day),by=48)],tcl=-.5)
  dev.copy2eps(file=paste("Factor_",i,"of",num_factor,".eps",sep = ""),family = "serif",width = 25)
}

