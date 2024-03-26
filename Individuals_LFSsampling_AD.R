
# LFS sampling design
LFS.sampling<-c(rep(TRUE,2),rep(FALSE,2),rep(TRUE,2))

# Output: data.frame "data.miss"
# data.frame "res" has 500 rows and 2 columns
# columns represent the information from ADMIN and LFS
# rows represent observations indiv*time
# 25 individuals (each individual belongs to an LFS group)
# 20 time points per each individuals

a<-NULL
res<-NULL
for (i in 0:24){
  a<-c(rep(FALSE,i),LFS.sampling,rep(FALSE,24-i))
  print(a)
  print(length(a))
  res<-rbind(res,a)
}

dim(res)
head(res)
sel.res<-res[,6:25]
dim(sel.res)
head(sel.res)
LFS<-as.vector(t(sel.res))
AD<-c(rep(TRUE,length(col.LFS)))
data.miss<-cbind(LFS,AD)

# to increase the number of indivuals 
# we require ("times" individuals) in each LFS group
times<-2
n.ind.tot<-25*times
data.miss.tot<-data.miss
for (i in 1:(times-1)){
  data.miss.tot<-rbind(d,data.miss)
}
dim(data.miss.tot)

