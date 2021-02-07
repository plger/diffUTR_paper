
rm(list=ls())
suppressMessages(library(tximport))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(polyester))
suppressMessages(library(ballgown))


set.seed(8675)

#import APA sites
extraColvect<-c( percentage="numeric", numberofprot="integer", tpm2="numeric", encod="character", addinfo="character")
APA<- rtracklayer::import(file.path("../data/atlas.clusters.2.0.GRCm38.96.bed.gz"),format="bed",extraCols=extraColvect)

#Ignore Antisense site for now
APA<-APA[(APA$encod!="AE"&APA$encod!="AI"&APA$encod!="AU"),]
#Extract only main chromosomes
APA <- keepStandardChromosomes(APA, pruning.mode="coarse")
#Collapse unused factors in seqlevels
seqlevels(APA)<-seqlevelsInUse(APA)
#Change seqnames from for example 1 to chr1
seqlevelsStyle(APA) <- "UCSC"

#load salmon estimates
quantnames<-"SRR39289"
quantnames<-paste0(quantnames,c("06","07","08","09","10","11"),"salmonquant")
quant_files = file.path("simulation", quantnames, "quant.sf")
txi = tximport(files = quant_files, type = "salmon", txOut = TRUE)

#delete version numbers at the end of transcript_id
row.names(txi$counts)<-substring(row.names(txi$counts), 1,18)

#get parameters
params<-polyester::get_params(txi$counts)
params<-as.data.frame(params[1:3])
params$transcript_id<-row.names(params)

#import gtf
gtf <- sort(rtracklayer::import("../data/large/gencode.vM25.annotation.gtf"))
gtf$gene_id<-substring(gtf$gene_id, 1,18)
gtf$transcript_id<-substring(gtf$transcript_id, 1,18)


#assign APA signs to genes
gtf1<-gtf[gtf$type=="CDS"|gtf$type=="stop_codon",]
#resize gene to 1 as only start of a gene would be a boundry
gtf1[gtf1$type=="gene",]<-resize(gtf1[gtf1$type=="gene",],1)



#find boundry region before each APA site...some have non, so ignore these and do it again
out<-follow(APA, gtf1)
APA<-APA[is.na(out)==FALSE]
out<-follow(APA, gtf1)

APA$gene_id<-gtf1$gene_id[out]

#only keep the ones salmon used aswell
gtf<-gtf[gtf$transcript_id %in% params$transcript_id,]


names(mcols(gtf))[names(mcols(gtf)) == "gene_type"] <- "gene_biotype"
#filter for only protein coding genes
gtf<-gtf[gtf$gene_biotype=="protein_coding"]

# Only keep relevant metadata
meta<-data.frame(gene_id=gtf$gene_id, type=gtf$type, transcript_id=gtf$transcript_id, gene_biotype=gtf$gene_biotype, stringsAsFactors=FALSE)
mcols(gtf)=meta

#merge with parameters
df<-as.data.frame(gtf)

df<-merge(df,params,by="transcript_id")
gtf<-makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)

#filter low expressed genes
gtf<-gtf[gtf$mu>0.1,]

#Extract only main chromosomes

gtf<-keepStandardChromosomes(gtf,pruning.mode="coarse")

#Collapse unused factors in seqlevels
seqlevels(gtf)<-seqlevelsInUse(gtf)





#get last exon of each transcrip
utrs <- gtf[gtf$type=="UTR"|gtf$type=="CDS",]


utrsminus<-utrs[as.character(strand(utrs))=="-",]
utrsplus<-utrs[as.character(strand(utrs))=="+",]

#split by transcripts
utrsplus<- split(utrsplus, utrsplus$transcript_id)
utrsminus<- split(utrsminus, utrsminus$transcript_id)

#find direct index of the last utr of each transcript and get them
l<-elementNROWS(utrsplus)
l<-cumsum(l)
lastUtrsplus<-unlist(sort(utrsplus))[l,]

l<-elementNROWS(utrsminus)
l<-cumsum(l)-l+1
lastUtrsminus<-unlist(sort(utrsminus))[l,]
lastUTR<-c(lastUtrsminus,lastUtrsplus)


#function to get ranges to the next utr that belong to the same gene

getnextUTR<- function(utrs,APA){

lastUtrsminus<-utrs[as.character(strand(utrs))=="-",]
lastUtrsplus<-utrs[as.character(strand(utrs))=="+",]
  
outmin<-precede(lastUtrsminus, APA)
lastUtrsminus<-lastUtrsminus[is.na(outmin)==FALSE]
outmin<-precede(lastUtrsminus, APA)

outplus<-precede(lastUtrsplus, APA)
lastUtrsplus<-lastUtrsplus[is.na(outplus)==FALSE]
outplus<-precede(lastUtrsplus, APA)


minusranges<-APA[outmin,]
plusranges<-APA[outplus,]

acc<-minusranges$gene_id==lastUtrsminus$gene_id
minusranges<-minusranges[acc,]
lastUtrsminus<-lastUtrsminus[acc,]
acc<-plusranges$gene_id==lastUtrsplus$gene_id
plusranges<-plusranges[acc,]
lastUtrsplus<-lastUtrsplus[acc,]

#create regions from last UTR to APA

beginplusbin<-end(lastUtrsplus)+1
endplusbin<-end(plusranges)

beginminusbin<-start(minusranges)
endminusbin<-start(lastUtrsminus)-1



plusbin<-GRanges(seqnames=seqnames(lastUtrsplus), IRanges(start=beginplusbin, end=endplusbin), strand=as.character(strand(lastUtrsplus)), gene_id=lastUtrsplus$gene_id,  transcript_id=lastUtrsplus$transcript_id,gene_biotype=lastUtrsplus$gene_biotype, type="exon", mu=lastUtrsplus$mu, p0=lastUtrsplus$p0, size=lastUtrsplus$size)

minusbin<-GRanges(seqnames=seqnames(lastUtrsminus), IRanges(start=beginminusbin, end=endminusbin), strand=as.character(strand(lastUtrsminus)), gene_id=lastUtrsminus$gene_id,  transcript_id=lastUtrsminus$transcript_id,gene_biotype=lastUtrsminus$gene_biotype, type="exon", mu=lastUtrsminus$mu, p0=lastUtrsminus$p0, size=lastUtrsminus$size)

#merge the strands to get all new defined UTRS
longUTR<-c(plusbin,minusbin)
  
return(longUTR)
}


#get next UTRs
longUTR<-getnextUTR(lastUTR,APA)

#filter UTR by length
longUTR<-longUTR[width(longUTR)>100]
longUTR<-longUTR[width(longUTR)<8000]

longUTR<-longUTR[!duplicated(longUTR$gene_id)]

#randome sample 1000 of them (transcripts to be regulated)
idx<-sample(length(longUTR),1000)
chosenUTR<-longUTR[idx,]

#find next APA for them 
UTR2<-getnextUTR(chosenUTR,APA)
#filter UTR by length
UTR2<-UTR2[width(UTR2)>100]
UTR2<-UTR2[width(UTR2)<8000]

#find next APA for them
UTR3<-getnextUTR(UTR2,APA)
#filter UTR by length
UTR3<-UTR3[width(UTR3)>100]
UTR3<-UTR3[width(UTR3)<8000]




#sample around 300 genes to be regulated (overlapping with the transcripts to be regulated)
idx<-sample(length(chosenUTR),as.integer(0.3*length(chosenUTR)))
diffgenes<-unique(chosenUTR$gene_id[idx])

#sample 400 genes to be regulated that don't overlapp with regulated transcripts
rest<-setdiff(unique(gtf$gene_id),unique(chosenUTR$gene_id))
idx<-sample(length(rest),300)
diffgenesrest<-rest[idx]

#all genes to be regulated
diffgenes<-c(diffgenes, diffgenesrest)





#duplicate all exons of the transcripts with new UTR and add new UTR to one of the duplicate for parsing by polyester (call the longer ones "transcript_id".x)

exons<-gtf[gtf$type=="exon",]
 #assign transcripts that are differentially expressed to categorys. comb03->all 3 additional bins differentialy expressed, comb13: after the first one until and with bin 3 differential and so on...
idx<-sample(length(UTR3$transcript_id),83)
comb03<-UTR3$transcript_id[idx]
restutr3<-setdiff(UTR3$transcript_id,comb03)
idx<-sample(length(restutr3),83)
comb13<-restutr3[idx]
comb23<-setdiff(restutr3,comb13)
used<-c(comb03,comb13,comb23)
restutr2<-setdiff(UTR2$transcript_id,used)
idx<-sample(length(restutr2),as.integer(0.5*length(restutr2)))
comb02<-restutr2[idx]
comb12<-setdiff(restutr2,comb02)
used<-c(used,comb02,comb12)
comb01<-setdiff(chosenUTR$transcript_id,used)


#create dataframe to be saved to see which transcripts are differential..set start end end to know after in which category the transcript belong
difftranscripts<-data.frame(gene_id=chosenUTR$gene_id,
                            transcript_id=chosenUTR$transcript_id, start=-1,end=-1,startpos=0,endpos=0, stringsAsFactors = FALSE)
difftranscripts$start[difftranscripts$transcript_id %in% c(comb01,comb02,comb03)]<-0
difftranscripts$start[difftranscripts$transcript_id %in% c(comb12,comb13)]<-1
difftranscripts$start[difftranscripts$transcript_id %in% comb23]<-2

difftranscripts$end[difftranscripts$transcript_id %in% c(comb23,comb13,comb03)]<-3
difftranscripts$end[difftranscripts$transcript_id %in% c(comb02,comb12)]<-2
difftranscripts$end[difftranscripts$transcript_id %in% comb01]<-1

#add additional UTR that is not differentially expressed to data
start1<-chosenUTR[chosenUTR$transcript_id %in% difftranscripts$transcript_id[difftranscripts$start>=1],]
start2<-UTR2[UTR2$transcript_id %in% difftranscripts$transcript_id[difftranscripts$start==2],]
exons<-c(exons,start1,start2)

#add the differentially expressed parts to the annotation
duplicates<-exons[exons$transcript_id %in% chosenUTR$transcript_id,]

end3<-UTR3[UTR3$transcript_id %in% difftranscripts$transcript_id[difftranscripts$end==3],]
end2<-UTR2[UTR2$transcript_id %in% difftranscripts$transcript_id[difftranscripts$end>=2] & UTR2$transcript_id %in% difftranscripts$transcript_id[difftranscripts$start<2],]
end1<-chosenUTR[chosenUTR$transcript_id %in% difftranscripts$transcript_id[difftranscripts$end>=1] & chosenUTR$transcript_id %in% difftranscripts$transcript_id[difftranscripts$start<1] ,]


duplicates<-c(duplicates,end3,end2,end1)
duplicates$transcript_id<-paste0(duplicates$transcript_id,".x")
exons<-c(exons,duplicates)


#sort all the data by transcript_id so the next step works
difftranscripts<-difftranscripts[order(difftranscripts$transcript_id),]
chosenUTR<-sort(chosenUTR, by=~transcript_id)
UTR2<-sort(UTR2, by=~transcript_id)
UTR3<-sort(UTR3, by=~transcript_id)


#get start and end positions of the differential expressed parts (use resize to 1 so it does take the actual start and end positions of the transcript regardless of strand)
rm(difftranscripts)
difftranscripts<-c(end3,end2,end1)




#duplicate the transcripts again this time to a data set of transcripts, which will be used to define fold change, #mean and size of each transcript

transcripts<-gtf[gtf$type=="transcript",]

tduplicates<-transcripts[transcripts$transcript_id %in% chosenUTR$transcript_id]
tduplicates$transcript_id<-paste0(tduplicates$transcript_id,".x")
transcripts<-c(tduplicates,transcripts)
transcripts<-data.frame(transcript_id=transcripts$transcript_id, gene_id=transcripts$gene_id,mu=transcripts$mu,size=transcripts$size)

#define foldchange
transcripts$fchg<-1

#create set of foldchanges
downregulation<-c(0.2,0.3,0.5,0.66,0.75)
upregulation<-1/downregulation
regulation<-c(upregulation,downregulation)

#sample fold changes for each gene that is to be regulated
generegulation<-sample(regulation,length(diffgenes),replace=TRUE)


#multiply the foldchanges for every transcript of that gene
for(i in 1:length(diffgenes)){
  transcripts$fchg[transcripts$gene_id==diffgenes[i]]=transcripts$fchg[transcripts$gene_id==diffgenes[i]]*generegulation[i]
}

diffgenes<-cbind(diffgenes, generegulation)

difftranscripts$foldchange<-1

#sample fold changes for each transcript that is regulated
transcriptregulation<-sample(regulation, length(chosenUTR),replace=TRUE)

#multiply the foldchanges for the transcripts, the counter part (.x) will be devided 
for(i in 1:length(chosenUTR)){
  tr=chosenUTR$transcript_id[i]
  trx=paste0(tr,".x")
  difftranscripts$foldchange[difftranscripts$transcript_id==tr]<-transcriptregulation[i]
  transcripts$fchg[transcripts$transcript_id==tr]=transcripts$fchg[transcripts$transcript_id==tr]*transcriptregulation[i]
  
  transcripts$fchg[transcripts$transcript_id==trx]=transcripts$fchg[transcripts$transcript_id==trx]/transcriptregulation[i]
  
}
diffgenes<-data.frame(gene_id = diffgenes[,1], foldchange=diffgenes[,2],stringsAsFactors =FALSE)
saveRDS(difftranscripts,file=file.path("difftranscripts.rds"))
saveRDS(diffgenes,file=file.path("diffgenes.rds"))

#safe the exon file for polyester
rtracklayer::export(exons, file.path("transcripts.gtf"))

gtfpolyester<-gffRead(file.path("transcripts.gtf"))

#sort transcripts with foldchanges, mu, size etc. so they have the right order

transcripts<-transcripts[order(transcripts$transcript_id),]

#create fold changes, first row all ones, second row what was calculated before
fold_changes<-rep(1,length(transcripts$fchg))
fold_changes<-cbind(fold_changes,transcripts$fchg)


#simulate experiment
simulate_experiment(seqpath=file.path("genome.by.chr"), gtf=file.path("transcripts.gtf"), reads_per_transcript=transcripts$mu, size=transcripts$size, num_reps=c(3,3), fold_changes=fold_changes, outdir=file.path("simulation","simulated.reads"), paired=FALSE,error_model="illumina5", bias="cdnaf", strand_specific=TRUE)





