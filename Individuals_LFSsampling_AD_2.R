
# Output: data.frame "data.miss"
# data.frame "data.miss" has 500 rows and 2 columns
# columns represent the information from ADMIN and LFS
# rows represent observations indiv*time
# 20 time points per each individuals

# n.occ - indicator for the number of (consecutive) observations
# n.ind - incicator for the number of individuals
# n.totocc - indicator for total time points

# res is the matrix with the LFS missing structure (stuck)
# individuals are randomly sampled for the matrix res

## False = missing

# LFS sampling design

lfs_str_miss <- function(N, Ntimes, prop_com){
  
  n.totocc <- Ntimes
  prop.com <- prop_com
  n.occ <- round(n.totocc*prop.com)
  n.ind <- N
  
  LFS.sampling<-rep(TRUE,n.occ)
  
  n.start <- n.totocc-n.occ
  a<-NULL
  res<-NULL
  for (i in 0:n.start){
    #n.false<-(n.totocc-n.occ)
    a<-c(rep(FALSE,i),LFS.sampling,rep(FALSE,n.totocc-length(LFS.sampling)-i ) )
    #print(a)
    #print(length(a))
    res<-rbind(res,a)
  }
  
  #dim(res)
  #head(res)
  
  indici <- sample(1:nrow(res), n.ind, replace = TRUE)
  LFS.sim<-res[indici,]
  LFS.sim<-as.vector(t(LFS.sim))
  AD<-c(rep(TRUE,length(LFS.sim)))
  data.miss<-cbind(LFS.sim,AD)
  
  ## TRUE is missing
  ## FALSE is observed
  return(!data.miss)
}



dd <- lfs_str_miss(N=25, Ntimes=20, prop_com = .2)
dd
