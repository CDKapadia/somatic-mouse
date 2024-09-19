require("ape")
require("parallel")
require("treemut")
require("methods")
#require("phytools")




get_dat=function(prefix,qc.dir=NULL){
  if(is.null(qc.dir)){
    load(sprintf("../post_qc/%s/%s_qc.Rdat",prefix,prefix))
  }else{
    cat("loading data from",sprintf("%s/%s/%s_qc.Rdat",qc.dir,prefix,prefix),"\n")
    load(sprintf("%s/%s/%s_qc.Rdat",qc.dir,prefix,prefix))
  }
  #pdx$summary$profile=pdx$df$df$profile[pdx$summary$ml]
  list(pdx=pdx,inf=inf)
}


setup_inference=function(inf,stub,##label for simulation (will create files that start with this stub)
                         niter=1e5,runs=runs,type="scite"){
  n=dim(inf$gt)[1]+length(inf$cfg$LOH)
  m=dim(inf$gt)[2]
  #browser()
  out.file.stub=sprintf("%s_%smuts_%ssamples",stub,n,m)
  ##Export the "noisy" genotypes in SCITE
  export(inf,out.file.stub)
  out.stub=out.file.stub
  if(type=="scite"){
    cmd=scite_wrapper(list(stub=out.file.stub,n=n,m=m),out.stub,niter = niter,runs = runs)
  }else if(type=="scite_t"){
    cmd=scite_wrapper(list(stub=out.file.stub,n=n,m=m),out.stub,niter = niter,runs = runs,b_transpose=TRUE)
  }else if(type=="mpboot"){
    cmd=mpboot_wrapper(list(stub=out.file.stub,n=n,m=m),out.stub,niter=niter)
  }
  else
  {
    cmd=sifit_wrapper(list(stub=out.file.stub,n=n,m=m),out.stub,niter = niter,runs = runs)
  }
  ##TODO  Can parallelise externally or use mclappy for a quick solution.
  writeLines(cmd$cmd,con = sprintf("%s.cmd",out.file.stub))
  out.file.stub
}

infer_tree=function(prefix,label,type="scite_t",niter=1e5,runs=5,mem.gb=2,queue="normal",qc.dir=NULL,allowed.type="SNV",output.dir="../results"){
  pd=get_dat(prefix,qc.dir=qc.dir)
  pdx=pd$pdx
  inf=pd$inf
  geno=cbind(pdx$dat$geno_colony,zeros=0)
  snames=colnames(geno)
  out=list(gt=geno,details=remap_details(pdx$dat$details),mut=cbind(pdx$dat$mtr,zeros=0)[,snames],ref=cbind(pdx$dat$wtr,zeros=0)[,snames],
           depth=cbind(pdx$dat$dep,zeros=30)[,snames],
           gpe=cbind(pdx$dat$gp,zeros=1e-330)[,snames])
  idx.keep.shared=which(rowSums(out$gt>0,na.rm=T)>1 & pdx$dat$details$TYPE %in% allowed.type)
  for(x in c("mut","ref","depth","gpe","gt","details")){
    #colnames(out[[x]])=c(gsub("^_","",snames),"zeros")
    out[[x]]=out[[x]][idx.keep.shared,]
  }
  for(event in inf$meta$LOH){
    out=add_loh(out,event$chr,event$start,event$end,inf$cfg$SHORT_LABEL[match(event$samples,inf$cfg$LABEL)])
  }
  #browser()
  stub=setup_inference(inf=out,sprintf("%s/%s_%s_%s",output.dir,prefix,type,label),niter=niter,runs=runs,type=type)
  done.file=sprintf("%s.cmd.DONE",stub)
  if(file.exists(done.file)){
    #stop(sprintf("%s already exists!",done.file))
  }else{
    if(type=="mpboot"){
      ## Super quick - no need to submit to farm
      system(readLines(sprintf("%s.cmd",stub))[1])
    }else{
      done.file=submit_lsf_jobs(sprintf("%s.cmd",stub),mem.gb = mem.gb,queue=queue)
      while(!file.exists(done.file)){
        Sys.sleep(30)
      }
    }
  }
  if(type=="scite_t"){
    best_tree=get_best_tree_scite(stub,niter=niter,runs = runs)
  }else if(type=="mpboot"){

    best_tree=read.tree(sprintf("%s.fa.parstree",stub))
    best_tree=root(best_tree,"zeros",resolve.root = TRUE)
    ##
    best_tree$edge.length=rep(1,dim(best_tree$edge)[1])
  }else{
    best_tree=get_best_tree(stub,niter=niter,runs = runs)
  }
  best_tree=read.tree(text = write.tree(best_tree,file = ""))
  list(tree=best_tree,stub=stub,niter=niter,runs=runs)
}

scite_wrapper=function(cfg,out.stub,runs=1,niter=1e5,ad=0.001,fd=0.01,b_transpose=FALSE){
  if(system("uname",intern = TRUE)=="Darwin"){
    SCITE="/Users/nw14/software/SCITE/scite"
    if(!file.exists(SCITE)){
      stop("Need to install SCITE locally on laptop or run this program on the server (e.g. cgpbar)")
    }
  }else{
    SCITE="/nfs/team78pc/nw14/software/SCITE/scite"
  }
  out.stub=sprintf("%s_iter%d",out.stub,niter)
  seed=sapply(sprintf("%s_run%d",out.stub,1:runs),function(s) sample(sum(utf8ToInt(s))+999999,size = 1))
  if(b_transpose){
    scite_cmd="%s -i %s_geno.txt -names %s_gene_names.txt -n %s -m %s  -r 1 -l %d -fd %3.1g -ad %3.1g -a -seed %d -max_treelist_size 1 -transpose -o %s_run%d"
  }else{
    scite_cmd="%s -i %s_geno.txt -names %s_gene_names.txt -n %s -m %s  -r 1 -l %d -fd %3.1g -ad %3.1g -a -seed %d -max_treelist_size 1 -o %s_run%d"
  }
  cmd=sapply(1:runs,function(r){sprintf(scite_cmd,#"%s -i %s_geno.txt -names %s_gene_names.txt -n %s -m %s  -r 1 -l %d -fd %3.1g -ad %3.1g -a -seed %d -o %s_run%d",
                                        SCITE,cfg$stub,cfg$stub,cfg$n,cfg$m,niter,fd,ad,seed[r],out.stub,r )})
  list(cmd=cmd,out.stub=sapply(1:runs,function(r){sprintf("%s_run%d",out.stub,r)}))
}

mpboot_wrapper=function(cfg,out.stub,runs=1,niter=1e5,ad=0.001,fd=0.01,b_transpose=FALSE){
  if(system("uname",intern = TRUE)=="Darwin"){
    stop("MPNboot not currently installed")
  }else{
    #MPNBOOT="/nfs/casm/team273jn/nw14/software/mpboot-1.1.0-install/mpboot"
    MPNBOOT="mpboot"
  }
  ##seed=sapply(sprintf("%s_run%d",out.stub,1:runs),function(s) sample(sum(utf8ToInt(s))+999999,size = 1))
  if(niter<1000){
    niter=1000
    cat("Number of bootstraps must be at least 1000.  Boosting to 1000\n")
  }
  if(niter>10000){
    niter=10000
    cat("Number of bootstraps can be at most 10000.  Capping at 10000\n")
  }
  cmd=sprintf("%s -s %s.fa -bb %d",MPNBOOT,out.stub,niter)
  list(cmd=cmd,out.stub=out.stub)
}

sifit_wrapper=function(cfg,out.stub,runs=1,niter=1e5,ad=0.001,fd=0.001){
  #if(system("uname",intern = TRUE)=="Darwin"){
  #  SCITE="/Users/nw14/software/SCITE/scite"
  #  if(!file.exists(SCITE)){
  #    stop("Need to install SCITE locally on laptop or run this program on the server (e.g. cgpbar)")
  #  }
  #}else{
  #  SCITE="/nfs/team78pc/nw14/software/SCITE/scite"
  #}
  out.stub=sprintf("%s_iter%d",out.stub,niter)
  seed=sapply(sprintf("%s_run%d",out.stub,1:runs),function(s) sample(sum(utf8ToInt(s))+999999,size = 1))
  cmdsifit=sprintf("awk '{print NR\"\\t\"$0}' %s_geno.txt | sed 's/\\t/ /g' > %s_sifit_geno.txt",cfg$stub,cfg$stub)
  #samples=writeLines()
  cat("Executing:\n",cmdsifit,"\n")
  system(cmdsifit)
  cmd=sapply(1:runs,function(r){sprintf(
    "java -Xmx20000m -Xms20000m -XX:+UseSerialGC -jar /nfs/team78pc/nw14/software/sifit/hamimzafar-sifit-e466245520f3/SiFit.jar -ipMat %s_sifit_geno.txt -n %s -m %s  -r 1 -iter %d -fp %3.1g -fn %3.1g -df 0 -cellNames %s_sample_names_sifit.txt",
    cfg$stub,cfg$n,cfg$m,niter,fd,ad,cfg$stub )})
  list(cmd=cmd,out.stub=sapply(1:runs,function(r){sprintf("%s_run%d",out.stub,r)}))
}

submit_lsf_jobs=function(job.file,mem.gb,queue="normal"){
  runAsLsfJobArray="/nfs/casm/team273jn/nw14/projects/runAsLsfJobArray.pl"
  N=length(readLines(job.file))
  #cmd=
  job.id=gsub("\\.","",format(Sys.time(), "simc%d%m%Y%H%M%OS3"))
  lsf.log=sprintf("%s%%I",job.file)
  mem_flag=sprintf("-R\"select[mem>%d] rusage[mem=%d]\" -M%d",1000*mem.gb,1000*mem.gb,1000*mem.gb);
  mem_flag2=sprintf("-R\"select[mem>%d] rusage[mem=%d]\" -M%d",100,100,100);

  cmd1=sprintf("bsub -J\"%s[1-%s]\" -o %s.log -q %s %s %s %s",
               job.id,N,lsf.log,queue,mem_flag,runAsLsfJobArray,job.file
  )
  done.file=sprintf("%s.DONE",job.file)
  cmd2=sprintf("bsub -w \"ended(%s)\" -q %s %s -o %s.done.log touch %s",job.id,queue,mem_flag2,job.file,done.file)

  cat(cmd1,"\n",cmd2,"\n")
  system(cmd1)
  system(cmd2)
  return(done.file)
  #callsys("bsub -w \"ended($job_id)\" -q $queue $mem_flag -o $out_dir/combo.log $combineChunksScript $out_dir $chr");
}

wait_for_done=function(done.file){
  while(TRUE){
    Sys.sleep(5)
    #cat("testing for ",done.file,"\n")
    if(file.exists(done.file)){
      return(0)
    }
  }
}


get_best_tree_scite=function(scite_stub,niter=niter,runs = runs){
  best_tree=get_best_tree(scite_stub,niter=niter,runs = runs)
  names=readLines(sprintf("%s_sample_names.txt",scite_stub))
  best_tree$tip.label=names[as.numeric(best_tree$tip.label)]
  best_tree$edge.length=rep(1,dim(best_tree$edge)[1])
  best_tree
}

get_best_tree_mpboot_bs=function(stub){
  best_tree=read.tree(sprintf("%s.fa.contree",stub))
  best_tree=root(best_tree,"zeros",resolve.root = TRUE)
  best_tree$boot_strap=as.numeric(best_tree$node.label)
  best_tree$edge.length=rep(1,dim(best_tree$edge)[1])
  best_tree
}

get_best_tree_mpboot=function(stub){
  best_tree=read.tree(sprintf("%s.fa.parstree",stub))
  best_tree=root(best_tree,"zeros",resolve.root = TRUE)
  best_tree$edge.length=rep(1,dim(best_tree$edge)[1])
  best_tree
}


get_best_tree=function(cmd,niter,runs){
  #browser()
  out=get_output(sprintf("%s.cmd",cmd),runs)
  scores=sapply(out,function(x) x$score)
  hist(scores)
  i=which.max(scores)
  cat("best score=",scores[i],"\n")
  if(is.null(out[[i]]$tree)){
    get_scite_tree(cmd,niter,i)
  }else{
    read.tree(text=out[[i]]$tree)
  }
}
get_output=function(cmd,runs){
  lapply(1:runs,function(i){read_lsf_out(sprintf("%s%d.log",cmd,i))})
}

read_lsf_out=function(lsf.log){
  #browser()
  l=readLines(lsf.log)
  if(length(grep("Success",l))==0){
    stop(sprintf("Tree Inference LSF job FAILED:%s",lsf.log))
  }
  if(length(grep("best LOH parameter",l))>0){
    idx.score=max(grep("best log-likelihood score = ",l))# Just in case its a re-run take the most recent LSF output
    score=as.numeric(gsub("best log-likelihood score = ","",l[idx.score]))
    tree=gsub("best tree = ","",l[max(grep("best tree = ",l))])
    nopt=1
    ##Example sifit
    #best false negative rate = 0.000004
    #best deletion parameter = 0.239989
    #best LOH parameter = 0.000000
    #best log-likelihood score = -2269.842648
    #best tree = ((((((s2
  }else{
    idx.score=max(grep("best log score for tree:",l))# Just in case its a re-run take the most recent LSF output
    score=as.numeric(gsub("best log score for tree:\t","",l[idx.score]))
    nopt=as.numeric(gsub(" opt trees","",l[max(grep("opt trees",l))]))
    if(length(grep("-transpose",l))>0){
      idx.tree=max(grep("samples from posterior written to: ",l))
      tree_path=gsub(".samples$","_ml0.newick",gsub("samples from posterior written to: ","",l[idx.tree]))
      tree=sprintf("%s;",readLines(tree_path))
    }else{
      tree=NULL
    }
  }
  if(length(grep("Successfully completed",l))>0){
    status=0
  }else{
    status=1
  }

  list(score=score,nopt=nopt,status=status,tree=tree)
}

get_scite_tree=function(inf){
  cat("kicking off ",label," \n")
  scite_stub=setup_inference(inf=inf,sprintf("output/%s_scite",label),niter=NITER,runs=runs,type="scite_t")
  done.file=sprintf("%s.cmd.DONE",scite_stub)
  if(file.exists(done.file)){
    #stop(sprintf("%s already exists!",done.file))
  }else{
    done.file=submit_lsf_jobs(sprintf("%s.cmd",scite_stub),mem.gb = 2,queue="normal")

    #sifit_stub=
    while(!file.exists(done.file)){
      Sys.sleep(30)
    }
  }
  best_tree=get_best_tree(scite_stub,niter=NITER,runs = runs)
  names=readLines(sprintf("%s_sample_names.txt",scite_stub))
  best_tree$tip.label=names[as.numeric(best_tree$tip.label)]
  best_tree$edge.length=1
}
assign_genotypes_mtr=function(tree,pdx,maxits=5,germline_trunk_length=3e4){
  #df=reconstruct_genotype_summary(best_tree)
  mtr=pdx$dat$mtr
  dep=pdx$dat$dep
  nm=dim(mtr)[2]
  nv=dim(mtr)[1]
  if(germline_trunk_length>0){
    #Add in additional germline variants to give appropriate prior to trunk
    mtr=rbind(mtr,matrix(10,ncol=nm,nrow=germline_trunk_length))
    dep=rbind(dep,matrix(20,ncol=nm,nrow=germline_trunk_length))
  }
  mtr=cbind(mtr,zeros=0)
  dep=cbind(dep,zeros=20)
  if(is.null(pdx$meta$error_rate)){
    cat("Assuming error rate of 0.01 - can configure this on per colony basis in inf$meta$error_rate parallel to inf$meta$clones_short\n")
    pdx$meta$error_rate=rep(0.01,length(pdx$meta$clones_short))
  }
  
  p.err=c(pdx$meta$error_rate,1e-6)
  p.err=p.err[match(tree$tip.label,c(pdx$meta$clones_short,"zeros"))]
  info=assign_to_tree(tree,mtr,dep,error_rate=p.err,maxits=maxits)
  
  pdx$tree_ml=info$tree
  pdx$summary=info$summary[1:nv,]
  pdx$tree_ml$edge.length=sapply(1:length(pdx$tree_ml$edge.length),function(x) length(which(pdx$summary$edge_ml==x)))
  pdx$df=info$df
  pdx
}

glue_infered_tree=function(output_rdat,pdat,output.dir="../results"){
  load(sprintf("%s/%s",output.dir,output_rdat))
  if(length(setdiff(pdat$inf$meta$clones_short,res$tree$tip.label))>0){
    stop("Mismatch")
  }
  ##Assign genotypes...
  #For this iteration we just assign snvs - correct polytomies so that only branches supported by at least 1 snv (rather than just indels)
  #are retained.
  snvpdx=pdat$pdx
  idx=which(pdat$pdx$dat$details$TYPE=="SNV")
  snvpdx$dat$mtr=snvpdx$dat$mtr[idx,]
  snvpdx$dat$dep=snvpdx$dat$dep[idx,]
  snvpdx=assign_genotypes_mtr(res$tree,snvpdx)
  idx.germline=which(snvpdx$tree_ml$edge[,1]==length(snvpdx$tree_ml$tip.label)+1   & snvpdx$tree_ml$edge[,2]>length(snvpdx$tree_ml$tip.label)+1)
  if(length(idx.germline)!=1){
    stop("Germline branch missing!")
  }
  snvpdx$tree_ml$edge.length[idx.germline]=1
  if(TRUE){
    while(length(which(snvpdx$tree_ml$edge.length==0 & snvpdx$tree_ml$edge[,2]>length(snvpdx$tree_ml$tip.label)))>0){
      cat("Collapsing 0 length branches and reassigning...\n")
      #Read write to ensure readable ordering of tips..  plot_tree isn't clever enough to do this..
      res$tree=read.tree(text=write.tree(di2multi(snvpdx$tree_ml),file=""))
      snvpdx=assign_genotypes_mtr(res$tree,snvpdx)
      idx.germline=which(snvpdx$tree_ml$edge[,1]==length(snvpdx$tree_ml$tip.label)+1   & snvpdx$tree_ml$edge[,2]>length(snvpdx$tree_ml$tip.label)+1)
      if(length(idx.germline)!=1){
        stop("Germline branch missing!")
      }
      snvpdx$tree_ml$edge.length[idx.germline]=1
    }
  }
  message("Now assigning SNVs and Indels")
  pdat$pdx=assign_genotypes_mtr(res$tree,pdat$pdx)
  pdat$pdx$summary$profile=pdat$pdx$df$df$profile[pdat$pdx$summary$edge_ml]
  pdat$pdx$tree_ml$label=pdat$pdx$df$df$profile
  pdat$pdx$dat$details$node=pdat$pdx$tree_ml$edge[pdat$pdx$summary$edge_ml,2]

  pdat
}
