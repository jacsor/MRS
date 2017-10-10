library(MRS)

file = "dnase_chr7_100807500_100809500_ENCODE_6_cells_tagcounts.csv"
begin = 100807500
end = 100809500

reads.mat = read.csv(file,header=TRUE)
colnames(reads.mat) = c("SampleID",1:(ncol(reads.mat)-1))

X = NULL
G = NULL
H = NULL

reads.new.all.long = data.frame(NULL)

## Figure 5 in Ma and Soriano (2017)
# pdf("k562vsMedullo721&341_VGF_long_count_histograms.pdf",width=7,height=7)
par(mfcol=c(3,3))
for (ID in reads.mat[,"SampleID"]) {
  print(ID)
  ID.temp = unlist(strsplit(ID,split="_Rep"))
  GroupID = ID.temp[1]
  RepID = ID.temp[2]
  
  reads.original = reads.mat[reads.mat$SampleID == ID,-1]
  reads.new = rep(1:length(reads.original),reads.original)
  
  
  if (length(grep("Medullo",GroupID))>0 | length(grep("K562",GroupID)>0))  {
    if (GroupID == "Medullo") GroupID = "Medullo_D721"
    reads.new.all.long = rbind(reads.new.all.long,data.frame(GroupID, SampleID = paste(GroupID, "Rep", RepID),Reads=reads.new+begin-1))
  }
  
  X = c(X,reads.new)
  G = c(G,rep(GroupID,length(reads.new)))
  H = c(H,rep(RepID,length(reads.new)))
}
# dev.off()


hist_K562 = ggplot(data = subset(reads.new.all.long,GroupID=="K562"), aes(x=Reads)) 
hist_K562 + geom_histogram(bins=1024,size=3) + ylim(0,40) + 
  scale_x_discrete(name ="Genomic location", breaks=c(begin,(begin+end)/2,end), limits=c(begin,end)) +
  facet_wrap(~SampleID,ncol=1,nrow=3) + theme_bw() + 
  geom_vline(xintercept = 100808844, color="red",linetype="dashed",size=0.5) + 
  geom_vline(xintercept = 100808876, color="red",linetype="dashed") + theme(plot.margin=unit(c(0,1,0,0),"cm"))

hist_Medullo_D721 = ggplot(data = subset(reads.new.all.long,GroupID=="Medullo_D721"), aes(x=Reads)) 
hist_Medullo_D721 + geom_histogram(bins=1024,size=3) + ylim(0,10) + 
  scale_x_discrete(name ="Genomic location", breaks=c(begin,(begin+end)/2,end), limits=c(begin,end)) +
  facet_wrap(~SampleID,ncol=1,nrow=3) + theme_bw() + 
  geom_vline(xintercept = 100808844, color="red",linetype="dashed",size=0.5) + 
  geom_vline(xintercept = 100808876, color="red",linetype="dashed") + theme(plot.margin=unit(c(0,1,0,0),"cm"))

hist_Medullo_D341 = ggplot(data = subset(reads.new.all.long,GroupID=="Medullo_D341"), aes(x=Reads)) 
hist_Medullo_D341 + geom_histogram(bins=1024,size=3) + ylim(0,22) + 
  scale_x_discrete(name ="Genomic location", breaks=c(begin,(begin+end)/2,end), limits=c(begin,end)) +
  facet_wrap(~SampleID,ncol=1) + theme_bw() + 
  geom_vline(xintercept = 100808844, color="red",linetype="dashed",size=0.5) + 
  geom_vline(xintercept = 100808876, color="red",linetype="dashed") + theme(plot.margin=unit(c(0,1,0,0),"cm"))




### Carry out the cross-group comparison

subset = (G=="Medullo_D721" | G=="Medullo_D341" | G=="K562")
G = G[subset]
X = X[subset]
H = H[subset]

G = as.numeric(factor(G))
X = as.numeric(X)
H = as.numeric(factor(H))

X = X + begin - 1
Omega = t(range(X));Omega[,2] = Omega[,2]+1

nu.vec = 10^seq(-1,4,length=10)
K=11
system.time({ans_andova = andova(X,G,H,K=K,nu=nu.vec,delta=0.4,gamma=0.07,beta=1,Omega = Omega)})
ans_andova$PostGlobNull


## Figure 6(a) in Ma and Soriano (2017)
# pdf(file="k562vsMedullo721vsMedullo341_VGF_ANDOVA_long.pdf",width=6)
plot1D(ans_andova,legend=TRUE,main="ANDOVA")
# dev.off()

# pdf(file="k562vsMedullo721vsMedullo341_VGF_ANDOVA_long_bfdr10.pdf",width=6)
plot1D_sig_windows(ans_andova,legend=TRUE,main="ANDOVA",fdr=0.1)
# dev.off()


rect(100808000,1,100809000,2)


## Figure 6(b) in Ma and Soriano (2017)
pdf(file="Effect_sizes_K562.pdf",width=6)
plot1D(ans_andova,type="eff",abs=FALSE,legend=TRUE,group=1,main="Effect sizes Leukemia K562")
dev.off()

pdf(file="Effect_sizes_MedulloD721.pdf",width=6)
plot1D(ans_andova,type="eff",abs=FALSE,legend=TRUE,group=2,main="Effect sizes Medullo D721")
dev.off()

pdf(file="Effect_sizes_MedulloD341.pdf",width=6)
plot1D(ans_andova,type="eff",abs=FALSE,legend=TRUE,group=3, main="Effect sizes Medullo D341")
dev.off()

system.time({ans_andova2 = andova(X,G,H,K=K,nu=nu.vec,delta=0.4,gamma=0.07,beta=1,Omega = Omega,base=1)})

pdf(file="Effect_sizes_K562_base1.pdf",width=6)
# plot1D(ans_andova2,type="eff",abs=FALSE,legend=TRUE,group=1,main="Effect sizes Leukemia K562")
plot1D(ans_andova2,type="eff",eff_scale="or",abs=FALSE,legend=TRUE,group=1,main="Effect sizes Leukemia K562")
dev.off()

pdf(file="Effect_sizes_MedulloD721_base1.pdf",width=6)
# plot1D(ans_andova2,type="eff",abs=FALSE,legend=TRUE,group=2,main="Effect sizes Medullo D721")
plot1D(ans_andova2,type="eff",eff_scale="or",abs=FALSE,legend=TRUE,group=2,main="Effect sizes Medullo D721")

dev.off()

pdf(file="Effect_sizes_MedulloD341_base1.pdf",width=6)
# plot1D(ans_andova2,type="eff",abs=FALSE,legend=TRUE,group=3, main="Effect sizes Medullo D341")
plot1D(ans_andova2,type="eff",eff_scale="or",abs=FALSE,legend=TRUE,group=3, main="Effect sizes Medullo D341")
dev.off()

## Plot the PMAP threshold versus Bayesian FDR
par(mfrow=c(1,1))
c_vec = seq(0,1,by=0.00001)
bFDR_vec = double(length(c_vec))
for (i in 1:length(c_vec)) {
  bFDR_vec[i] = bFDR(ans_andova2,c_vec[i])
}
pdf(file="bFDR_vs_c.pdf",width=6)
plot(c_vec,bFDR_vec,type='l',main="Bayesian FDR vs PMAP threshold",xlab = "PMAP threshold c", ylab="Bayesian FDR")
abline(a=0.1,b=0,lty="dashed")
dev.off()

## Draw and plot posterior samples of effect sizes and the states
n_post_samples = 25
ans_andova = andova(X,G,H, K=11,delta=0.4,gamma=0.07,beta=1,Omega = Omega,n_post_samples = n_post_samples,nu_vec = 10^seq(-1,4,length=10),baseline=1)
for (sample_id in 1:n_post_samples) {
  
  ans_andova$RepresentativeTree = ans_andova$PostSamples[[sample_id]]
  ans_andova$RepresentativeTree$EffectSizes[is.nan(ans_andova$RepresentativeTree$EffectSizes)]=0
  pdf(file=paste("posterior_sample_draw",sample_id,".pdf",sep=""),width=7,height=7)
  plot1D(ans_andova, abs=FALSE, type = "eff", eff_scale = 'or', legend = T, group = 2, main =paste("Effect size of Medullo D721 - Sample",sample_id))
  dev.off()
    # plot1D(ans_andova, legend = T, group = 1, main =paste("Sample",sample_id))
  
}

