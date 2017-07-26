library(GenomicRanges)
library(Biostrings)
library(ggplot2)
library(plyr)
library(doParallel)
registerDoParallel(cores=4)
library(qpcR)
library(reshape2)
library(flux)
library(dplyr)
library(data.table)
library(Gviz)

setwd("~/LBC_projects/TFBSAnalyzer/LargeScale/Promoters_example/new_scan_MC/")

###cria lista dos fastas experimentais - essas sequencias estão no cluster 10.0.2.146:/home/tfbs/work/dna_seqs
dir_fasta = "./dna_seqs/"
lst_fasta = list.files(dir_fasta,pattern = ".fasta")
fasta_list = sapply(lst_fasta, function(x) readDNAStringSet(paste(dir_fasta,x,sep="")))
###carrega distribuiçao ordenada - a distribuição está em 10.0.2.146:/home/tfbs/work/amber99-dist
NR_dist = read.table("../../distributions/P24468.dist",  header = F)
NR_name = "P24468"
names(NR_dist)=c("seq","energy")

#PPARG.dist=read.table("../../distributions/P37231.dist",sep="\t",header = F)
RXRA.dist=read.table("../../distributions/P19793.dist",sep="\t",header = F)



# Binding Simulation
COUPT.homo.dist=join.NR.dist(NR_dist,NR_dist,Energy.column=2)
tab1=as.data.table(COUPT.homo.dist)
mc.d=data.frame()
for(dr in 1:6){
  d.d = lapply(1:length(fasta_list), function(x,y) MMC.TFBS(y[[x]],dr,tab1,10000,8),y=fasta_list)

  c=0
  for(x in 1:length(fasta_list)) {
    mc.d=rbind(mc.d,d.d[[x]])
    c=c+dim(d.d[[x]])[1]
  }
  
}

mc.d=cbind(mc.d,end1=mc.d$start1+5,end2=mc.d$start2+5)
seq.names=sapply(1:length(fasta_list),function(x)names(fasta_list[[x]][1]))
partititions=sapply(1:length(fasta_list),function(x)partition(seq.names[x],d.d[[x]]))


PlotSeq(mc.d,"chr8",1)






###################################################

MMC.TFBS=function(seq.1,dr=3,lookup.table,nsteps,exclusion_cutoff,temp=298){
  #Monte Carlo Metropolis simulation
  
  #Monte Carlo parameters
  R      = 1.9872156e-3    # Gas constant  kcal/mol/degree
  beta   = 1 / ( R * temp)
  
  mc.d=data.frame()
  seq.name=names(seq.1)
  len.1=width(as.character(seq.1[[1]]))
  E_old=0
  timestamp()
  for( i in 1:nsteps){
    #dr=sample(0:6,1)
    start1=sample(1:(len.1-12-dr),1)
    end1=start1+5
    start2=start1+5+dr
    end2=start2+5
    bs1=as.character(seq.1[[1]][start1:end1])
    bs2=as.character(seq.1[[1]][start2:end2])
    
    dt1=lookup.table[X1==bs1, ]
    E_new=dt1[dt1$Y1==bs2,]$score
    
    if(E_new> E_old){
      d.tmp=data.frame(i=i,bs1=bs1,bs2=bs2,start1=start1,start2=start2,dr=dr,Energy=E_new,prob=1,rndn=0,delta=0,seq.name=seq.name)
      mc.d=rbind(mc.d,d.tmp)  
      E_old=E_new
    } else { # Try Metropolis criterium
      delta = E_old - E_new 
      prob.1 = exp( - beta * delta * exclusion_cutoff )
      rndn=runif(1,0:1)
      if (prob.1 > rndn) {
        d.tmp=data.frame(i=i,bs1=bs1,bs2=bs2,start1=start1,start2=start2,dr=dr,Energy=E_new,prob=prob.1,rndn=rndn,delta=delta,seq.name=seq.name)
        mc.d=rbind(mc.d,d.tmp)
        E_old=E_new
      }
      
    }
  }
 mc.d
}

partition=function(seqname,mc.d,temp=298){
  R      = 1.9872156e-3    # Gas constant  kcal/mol/degree
  beta   = 1 / ( R * temp)
  
  #Z partition funcion
  Energy=subset(mc.d,seq.name==seqname)$Energy
  Z=sum(exp(beta*Energy))
  # The probability to find the system in some energy eigenstate r is given by:
  Prnum=exp(beta*Energy)/Z
  
  #The average internal energy of the system is the expectation value of the energy and can be expressed in terms of Z as follows:
  U=sum(Prnum*Energy)
  S.1=-1*sum(Prnum*log(Prnum))
  # Free energy
  As=-1*R*temp*log(Z)
  
  cat(seqname,"\t",U,"\t",S.1,"\t",As,"\n")
  
}

# Create distribution for dimers
join.NR.dist=function(df1,df2,Energy.column,no.use.solv=T){
  ##combinaçao de todas as sequencias da distribuiçao all vs all
  x1 = df1$seq
  x2 = df1[,Energy.column]
  y1 = df2$seq
  y2 = df2[,Energy.column]
  d1 = expand.grid(X1 = x1, Y1 = y1)
  d2 = expand.grid(X2 = x2, Y2 = y2)
  
  ##soma do score de cada combinaçao e criaçao de um dataframe ordenado pela soma dos scores
  df_sum = d2$X2 + d2$Y2
  df = cbind(d1,d2,sum = df_sum)
  df_order = df[order(df$sum, decreasing=no.use.solv),]
  rm(list=c("d1","d2","df","df_sum"))
  
  ##nomeia as linhas
  rownames(df_order) = seq(1,dim(df_order)[1])
  ##normaliza a coluna de soma e cria uma coluna de score normalizado
  #min_sum = min(df_order$sum)
  #max_sum = max(df_order$sum)
  #Very slow
  #score = sapply(df_order$sum, function(x) (x-min_sum)/(max_sum-min_sum))
  score=1-(as.numeric(rownames(df_order))/dim(df_order)[1]) 
  df_order = cbind(df_order, score = score)
  df_order
}


PlotSeq=function(mc.d,seqname,seqnum,coord.1=1,coord.2=width(fasta_list[[seqnum]]),color.level=50){
  #subset dataset & sort lines by start1
  d.tmp=subset(mc.d,seq.name==seqname)
  d.tmp=d.tmp[order(d.tmp$start1),]
  gr.tfbs=makeGRangesFromDataFrame(d.tmp,start.field = "start1",end.field="end2",seqnames.field ="seq.name",keep.extra.columns = T)
  
  colgrad=c(colors()[grep("blue",colors())][1:color.level],colors()[grep("red",colors())][color.level:1])
  sTrack <- SequenceTrack(fasta_list[[seqnum]])
  axisTrack = GenomeAxisTrack()
  scores.track=DataTrack(gr.tfbs[,c("Energy")],name="Score",type=c("b"))
  tfbs.track= AnnotationTrack(gr.tfbs,name=seqname,fill=c(colgrad[round(gr.tfbs$Energy*100)]))
  plotTracks(list(sTrack,scores.track,tfbs.track),from=coord.1,to=coord.2)
}
