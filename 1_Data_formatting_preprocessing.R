#Preprocessing for converting raw data to adjacency matrix
#Masashi Komori

# 1.Joining Raw Data files  -----------------------------------------------
print(csvflist <- list.files(pattern = "*[[:digit:]].csv")) ##List the csv files in the folder

for(i in 1 : 100){#Combine csv files starting with "BL0" by number into Matome_BL00XX.csv
  fname <- paste("BL0", formatC(i, width = 3, flag = "0"), "_[0-9]*.csv", sep="")
  tmpflist <- grep(fname , csvflist)
  
  if(length(tmpflist) != 0){
    tmp <- NA
    for(j in 1 :  length(tmpflist)){
      tmp2 <- try(read.csv(csvflist[tmpflist[j]],header=F))
      if(length(tmp2)!=1){
        tmp <- rbind (tmp, try(read.csv(csvflist[tmpflist[j]],header=F)))
      }
    }
    tmp<-na.omit(tmp)
    colnames(tmp) <- c("Device","MAC","Year", "Mon", "Day", "Hour", "Min", "Sec" )
    tmp$Min <- tmp$Min - tmp$Min %% 5
    tmp$date <- paste(tmp$Year,"/",tmp$Mon,"/",tmp$Day," ",tmp$Hour,":",tmp$Min,sep="")
    write.csv(tmp,paste("Matome_BL0", formatC(i, width = 3, flag = "0"), ".csv", sep=""),row.names=F)
  }
}

# 2.Obtaining a combined data file ----------------------------------------
#Obtain Matome_BL00XX.csv
print(csvflist <- list.files(pattern = "Matome_BL[[:digit:]]*.csv")) 
n_csv <- length(csvflist) #Count the number of csv files
d.list <- lapply(csvflist,read.csv, header=T, sep=",", colClasses=c("character", "character")) 

start_node <- rep(NA, n_csv) #The nodes of starting point of the edges
end_node <- NA #The nodes of end point of edges

mintime <- min(strptime(d.list[[1]]$date, "%Y/%m/%d %H:%M")) #Day measurement start
maxtime <- max(strptime(d.list[[1]]$date, "%Y/%m/%d %H:%M")) #Day measurement end

for(i in 1 : n_csv){
  d.list[[i]]$time <- strptime(d.list[[i]]$date, "%Y/%m/%d %H:%M")   
  if(length(d.list[[i]]$time) != 0){ 
    if(mintime > min(na.exclude(d.list[[i]]$time))){ 
      mintime <- min(d.list[[i]]$time)
    }
    if(maxtime < max(d.list[[i]]$time)){
      maxtime <- max(d.list[[i]]$time)
    }
    
    end_node <- c(end_node, d.list[[i]]$Device) #List the detected devices.
  }
  start_node[i] <- gsub("Matome_","",gsub(".csv","",csvflist[i])) #List the detecting devices.
}

node <- unique(na.exclude(c(start_node, end_node)))
init <- floor(as.numeric(mintime)/(60*60*24))*60*60*24 
end <- ceiling(as.numeric(maxtime)/(60*60*24))*60*60*24 #transrate to UNIX time


# 3.Creating adjacent arrays ----------------------------------------------
subnode <- sort(node[grep("BL0",node)])

adj_array2 <- array(0, dim=c(length(subnode), length(subnode), (end-init)/(5*60)+1))
dimnames(adj_array2) <- list(subnode,subnode,seq(init, end, length=(end-init)/(5*60)+1)) 

for (i in 1:n_csv){
  if(nrow(d.list[[i]]) != 0){
    for(j in 1:nrow(d.list[[i]])) {
      if(length(grep("BL0",d.list[[i]]$Device[j])) != 0){
        adj_array2[as.character(start_node[i]),
                   as.character(d.list[[i]]$Device[j]),
                   as.character(as.numeric(d.list[[i]]$time[j]))] <-1
      }
    }
  }
}

for (i in 1 :length(adj_array2[1,1,])){
  adj_array2[,,i] <- adj_array2[,,i] + t(adj_array2[,,i] ) - adj_array2[,,i] * t(adj_array2[,,i] ) 
}


# 4.Apply a rolling filter ------------------------------------------------
ratio <- 6 #every 30 minutes
adj_array3 <- array(0, dim=c(nrow(adj_array2), ncol(adj_array2), ((end-init)/(5*60*ratio) +1))) 
dimnames(adj_array3) <- list(dimnames(adj_array2)[[1]],dimnames(adj_array2)[[2]],seq(init, end, length=((end-init)/(5*60*ratio)+1) )) 
for(i in 1 : (length(adj_array2[1,1,])/ratio))
{
  adj_array3[,,i] <- apply(adj_array2[,,(((i-1)*ratio+1):(i*ratio))], c(1,2) , sum)
}

adj_array <- adj_array3
#adj_array <- adj_array2


# 5.Set the period to be analyzed. ----------------------------------------
anlystart <- as.character(as.numeric(strptime("2017/12/01 00:00", "%Y/%m/%d %H:%M", "Japan"))) 
anlyend <- as.character(as.numeric(strptime("2018/6/30 00:00", "%Y/%m/%d %H:%M", "Japan"))) 
adj_array <- adj_array[,, (which(dimnames(adj_array)[[3]] == anlystart):(which(dimnames(adj_array)[[3]] == anlyend)))]


# 6.Converting a 3D array to a matrix -------------------------------------
adj_vec <- matrix(NA,length(adj_array[1,1,]), nrow(adj_array) * (ncol(adj_array)-1)/2)
rownum <-  matrix(1:ncol(adj_array), nrow=nrow(adj_array), ncol=ncol(adj_array))
colnum <- t(rownum)
fromname <- rownames(adj_array[,,1])[rownum[row(rownum)>col(rownum)]]
toname <- rownames(adj_array[,,1])[colnum[row(colnum)>col(colnum)]]

for(i in 1 : length(adj_array[1,1,]))
{
  tmp <- adj_array[,,i]
  adj_vec[i,] <- tmp[row(tmp)>col(tmp)] 
}

sub_adj_vec <- adj_vec[which(rowSums(adj_vec) != 0),which(colSums(adj_vec) != 0)]
fromname <- fromname[which(colSums(adj_vec) != 0)]
toname <- toname[which(colSums(adj_vec) != 0)]
time <- dimnames(adj_array) [[3]][which(rowSums(adj_vec) != 0)]
sub_adj_vec <- as.matrix(sub_adj_vec);

save(sub_adj_vec,file="sub_adj_vec.RData")
