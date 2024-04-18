
# Output: data.frame "data.miss"
# data.frame "data.miss" has 500 rows and 2 columns
# columns represent the information from ADMIN and LFS
# rows represent observations indiv*time
# 20 time points per each individuals

# n.occ - indicator for the number of (consecutive) observations
# n.ind - incicator for the number of individuals

# res is the matrix with the LFS missing structure (stuck)
# individuals are randomly sampled for the matrix res

# LFS sampling design
n.occ<-10
LFS.sampling<-rep(TRUE,n.occ)

n.start<-20-n.occ
a<-NULL
res<-NULL
for (i in 0:n.start){
  n.false<-(20-n.occ)
  a<-c(rep(FALSE,i),LFS.sampling,rep(FALSE,20-length(LFS.sampling)-i))
  print(a)
  print(length(a))
  res<-rbind(res,a)
}

dim(res)
head(res)

n.ind<-25
indici <- sample(1:nrow(res), n.ind, replace = TRUE)
LFS.sim<-res[indici,]
LFS.sim<-as.vector(t(LFS.sim))
AD<-c(rep(TRUE,length(LFS.sim)))
data.miss<-cbind(LFS.sim,AD)





