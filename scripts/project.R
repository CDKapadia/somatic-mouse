##  MD
SPECIES="Mmusculus" # Hsapiens
SPECIES_BUILD="mm10"
PATIENTS=c("MD7180","MD7181","MD7182","MD7183","MD7184","MD7185")
CV=c("missense","nonsense","ess_splice","frameshift","inframe","loh","start_lost","cna","stop_lost")
GENES=readLines("GENES.txt")

get_driver_scheme2=function(){
  driver.scheme=read.table("driver_scheme_simple.txt",head=T,stringsAsFactors = FALSE,sep="\t")
  n=max(driver.scheme$number)
  pallete=c(RColorBrewer::brewer.pal(9,"Set1")[-6],RColorBrewer::brewer.pal(8,"Dark2"))
  driver.scheme$colour=pallete[driver.scheme$number]
  driver.scheme
}

##Import age of sample data
##Add all relevant meta
##  Could also add columns to CFG here for individul specific info

add_agedf=function(PD){
  tree=PD$pdx$tree_ml
  agefile="~/Desktop/config_files/sample_ages.txt"
  if(!file.exists(agefile)){
    cat("can't find",agefile,"\n")
    cat("need to creat tab delimited file ",agefile," with column headings sample and age_at_sample.  Note that sample column should contain sample long names")
    stop("unable to find age file ")
  }
  sample=PD$pdx$cfg$LABEL[match(PD$pdx$tree_ml$tip.label,PD$pdx$cfg$SHORT_LABEL)]
  if(length(which(is.na(sample)==1))){
    sample[is.na(sample)]="zeros"
  }else{
    stop("unexpected lookup error!")
  }
  
  agedf=data.frame(tip.label=tree$tip.label,sample=sample)
  af=read.table(agefile,stringsAsFactors=FALSE,sep="\t",head=TRUE)
  agedf=agedf %>% left_join(af) %>% dplyr::rename(age_at_sample_exact=age_at_sample)
  agedf$age_at_sample_pcy=agedf$age_at_sample_exact+(MEAN_AGE_AT_DELIVERY_DAYS)/365.25
  agedf$age_at_sample_pcy[which(agedf$sample=="zeros")]=1e-6
  PD$pdx$agedf=agedf
  PD
}
