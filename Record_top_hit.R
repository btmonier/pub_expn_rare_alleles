#Karl Kremling
#Record the top marker in a set window
#December 12 2015

gff_file=read.table("/home/kak268/work_data/maize_genome/Zea_mays_transcript_only.AGPv3.29.gff3", header=F, skip=6, sep="\t")
colnames(gff_file)=c("Chr", "Annot_source", "Feature_type", "Start", "End", "blankcol1", "Strand", "blankcol2", "Feature_Name")

transcript_file=gff_file[gff_file$Feature_type=='transcript',] # keep only the entries in the GFF file which are transcripts
transcript_file$Feature_Name=sub(pattern= "ID=transcript:", replacement='', x=transcript_file$Feature_Name)
transcript_file$Feature_Name=sub(pattern= ";.*", replacement='', x=transcript_file$Feature_Name)



Record.top.hit.in.window <- function(x, minlog10PVal, zoomradius) {
  x <- x[-log10(x[,"p"])>=minlog10PVal,] # This keeps only the pvalues above(below) yhe cutoff
  genename=toString(unique(x$Trait))
  cat('\n genename is ', genename, " ")
  transcriptFileLine=transcript_file[transcript_file$Feature_Name==genename,] #"GRMZM2G093344_T01",]
  #cat(toString(transcriptFileLine))
  genechr=toString(transcriptFileLine$Chr)
  cat(' genechr is ', genechr, '\n')
  if (genechr!=""){ # If there is a gff line for the trait that was mapped
    leftbound=max(transcriptFileLine$Start-zoomradius, 0) 
    rightbound=transcriptFileLine$End+zoomradius
    curr <- x[x$Chr==genechr,] #subset of only hits on the chr where the annotation says the gene is
    curr <- curr[curr$Pos>=leftbound,]
    curr <- curr[curr$Pos<=rightbound,]
    #cat("transcriptFileLine dim is ", dim(transcriptFileLine), "\n")
    #cat("min(curr$p)dim is", min(curr$p), "\n")
    #cat("curr dim is", dim(curr), "\n")
    colnames(curr)[colnames(curr) == "Chr"] <- "HitChr" # change name so it doesn't clash with the Chr in the gff entry
    if (dim(curr)[1]!=0){
      recorded.hit=cbind(transcriptFileLine, curr[curr$p==min(curr$p),], zoomradius)
      return(recorded.hit)
    }
    else{ # necessary to record the gff line of current gene, but note there was no cis hit
      #print("curr is ", curr) 
      curr=rbind(curr, data.frame(Trait='noCisHit', Marker='noCisHit', HitChr='noCisHit', Pos='noCisHit', df='noCisHit', r2='noCisHit', p='noCisHit'))
      #print("made it to line 31")
      #print("curr with blanks dim is ", dim(curr)) 
      recorded.hit=cbind(transcriptFileLine, curr, zoomradius)
      return(recorded.hit)
    }
    #line=cbind(transcriptFileLine, curr[max(curr$p),])
    #line # this returns the gff entry of the current gene, and the top SNP hit within the
  }
}
