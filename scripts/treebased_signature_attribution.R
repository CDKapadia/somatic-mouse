get_clade_sigs_mouse=function(PD,prior,sigmethod="SIGFIT",extracats=c("Top50","Top100")){
  stats=NULL
  CTX=NULL
  mp=list()
  for(categ in c("Shared","Private",extracats)){
    cat("Processing",categ,"..\n")
    details=get_mouse_details(PD,categ)
    mutc=with(details,get_mut_matrix(Chrom,Pos,Ref,Alt))
    ctx=as.data.frame(mutc) %>% mutate(count=gr,trin=rownames(.)) %>% mutate(ctx=substr(trin,3,5)) %>% 
      mutate(ctx=ifelse(ctx=="C>T","C>T other",ctx)) %>%
      mutate(ctx=ifelse(grepl("\\[C>T\\]G",trin),"C>T at CpG",ctx)) %>% group_by(ctx) %>% summarise(N=n(),count=sum(count)) %>%
      mutate(donor=PD$patient,category=categ)
    CTX=rbind(CTX,ctx)
    sig=matrix(as.vector(mutc),ncol=1)
    if(sigmethod=="SIGFIT"){
      sf=fit_signatures(t(sig),t(prior),iter=10000,show_messages=F)
      exposures=retrieve_pars(sf,par="exposures")
      decomp=list(contribution=t(as.matrix(exposures$mean)))
      colnames(decomp$contribution)=categ
      decomp$csm=cos_sim(as.vector(retrieve_pars(sf,par="reconstructions")$mean),as.vector(sig))
      ## Following is a bit ugly - but emulates the structure emitted by MP method below.
      qq=rbind(do.call("rbind",exposures) %>% (function(x){rownames(x)=c("mean","2.5%","97.5%");x}),
               total=rep(sum(mutc),dim(exposures$mean)[2]))
    }else{
      ## Mutational patterns fit
      decomp=MutationalPatterns:::fit_to_signatures(sig,prior)
      decomp$contribution[,1]=decomp$contribution[,1]/sum(decomp$contribution[,1])
      decomp$csm=cos_sim(decomp$reconstructed[,1],as.vector(mutc))
      colnames(decomp$contribution)=categ
      ## Next do bootstrapping to get CIs.
      bs=fit_to_signatures_bootstrapped(sig,prior,verbose = FALSE,method = "regular") %>% (function(probs) sweep(probs,MARGIN=1,rowSums(probs),"/"))
      qq=sapply(colnames(bs),function(x) {qq=quantile(bs[,x],c(0.025,0.5,0.975));c(qq,mean=mean(bs[,x],),sd=sd(bs[,x]),total=sum(mutc))})
      
    }
    stats=rbind(stats,as.data.frame(qq) %>% mutate(field=rownames(.),category=categ,donor=PD$patient,rec_csm=decomp$csm) )
    MP=as.matrix(get_mutation_probability(prior,t(decomp$contribution))[,-(1:2)])
    mp[[categ]]=MP
    cat("\n")
  }
  list(sigstats=stats %>% pivot_longer(cols=colnames(prior),names_to = "signature",values_to = "proportion") %>%
         pivot_wider(names_from="field",values_from = "proportion"),
       ctx=CTX,MP=mp)
  
}

plot_mouse_sig_tree_posterior=function(PD,prior,b.ultra=FALSE,txtitle=PD$patient,b.add.sigs=TRUE,focal.signature="SBS1",
                                       b.ggbar_separate=F,lwd=1.3,default_edge_color="grey30"){
  info=PD$pdx$dat$info
  if(is.null(info)){
    stop('prior to running please do: PD$pdx$dat$info=get_clade_sigs_mouse(PD,prior=pcawg[, c("SBS1","SBSblood","SBS18")])')
    info=get_clade_sigs_mouse(PD,prior)
  }
  if(is.null(PD$pdx$dat$sigbynode)){
    PD2=add_mouse_sig_to_pdx(PD,prior=prior,info = info)
  }else{
    PD2=PD
  }
  colorscheme=get_sig_color_scheme() %>% filter(SIG %in% colnames(prior))
  if(b.ultra){
    PD2$pdx$tree_ml=PD2$fit$poisson_tree$altmodel$ultratree
  }
  # browser()
  tree=plot_tree(PD2$pdx$tree_ml,cex.label = 0,left.margin.prop = 0.3,ymax=250,lwd=lwd,default_edge_color=default_edge_color)#,mar = c(1, 2, 1, 1) + 0.1)
  title(txtitle);
  tree=plot_signature_annotated_tree_mouse(tree,PD2$pdx$dat,colorscheme,maxlen=1e6,b.include.csm.text = FALSE)
  ## Convert data to form ## For 
  if(!b.add.sigs){
    agg1=info$ctx %>% group_by(category) %>% summarise(total=sum(count))
    df=info$ctx %>% filter(ctx=="C>T at CpG") %>% left_join(agg1)
    bounds=t(sapply(1:dim(df)[1],function(i) binom.test(df$count[i],df$total[i])$conf.int))
    colnames(bounds)=c("lb","ub")
    df=cbind(df %>% mutate(value=count/total),as.data.frame(bounds))
    df=df %>% left_join(data.frame(category=c("Shared","Private"),label=c("Shared","Private")))
    add_bar_chart(df %>% mutate(type=category),col2 =c("darkorange","white"),txttitle="C>T @ CpG")
  }else{
    # browser()
    if(!b.ggbar_separate){
      df=info$sigstats %>% mutate(lb=`2.5%`,ub=`97.5%`, value=mean) %>% filter(signature==focal.signature) %>% mutate(label=category)
      add_bar_chart(df %>% mutate(type=category),c(colorscheme %>% filter(SIG==focal.signature) %>% pull(COL),"pink"))
    }else{
      df=info$sigstats %>% mutate(lb=`2.5%`,ub=`97.5%`, value=mean) %>% mutate(label=category)
      #  left_join(data.frame(category=c("Shared","Private"),label=c("Shared","Private")))
      ordering = c('SBS1','SBS5', 'SBS18')
      df$signature = factor(df$signature, levels=ordering)
      #Stupid ggplot..
      df[which(df$signature==ordering[2]),]$lb=df[which(df$signature==ordering[2]),]$lb+ df[which(df$signature==ordering[1]),]$mean
      df[which(df$signature==ordering[2]),]$ub=df[which(df$signature==ordering[2]),]$ub+ df[which(df$signature==ordering[1]),]$mean
      df[which(df$signature==ordering[3]),]$lb=df[which(df$signature==ordering[3]),]$lb+ df[which(df$signature==ordering[1]),]$mean +df[which(df$signature==ordering[2]),]$mean
      df[which(df$signature==ordering[3]),]$ub=df[which(df$signature==ordering[3]),]$ub+ df[which(df$signature==ordering[1]),]$mean + df[which(df$signature==ordering[2]),]$mean
      
      p1<-
        ggplot(df, aes(fill=signature, y=value, x=label)) + 
        geom_bar(position=position_stack(reverse = TRUE), stat="identity")+
        geom_linerange( aes(x=label, ymin=lb, ymax=ub),position = position_dodge(width=0.15), colour="black", alpha=0.9, size=1)+
        # geom_errorbar( aes(x=label, ymin=lb, ymax=ub),position = "dodge", colour="black", alpha=0.9, size=1)+
        theme_cowplot()+
        ylab("Exposure (%)")+xlab(NULL)+theme(legend.position = "none")+
        scale_y_continuous(breaks= seq(0,1,by=.2))+
        scale_fill_manual(breaks=colorscheme$SIG, values=colorscheme$COL)
      df$sample=PD2$patient
      # sigdf_list <<- append(sigdf_list,list(df))
      return(p1)
    }   
    
  }
  tree
}

## Adds grouped category level signature to the tree 
add_mouse_sig_to_pdx=function(PD,prior,info=get_clade_sigs(PD,prior)){
  sigbynode=list()
  sigs=colnames(prior)
  empty=matrix(NA,ncol=1,nrow=length(sigs))
  rownames(empty)=sigs
  #stats=NULL
  for(categ in c("Shared","Private")){
    cat("Processing",categ,"..\n")
    details=get_mouse_details(PD,categ)
    MP=info$MP[[categ]]
    for(node in unique(details$node)){
      i=which(PD$pdx$tree_ml$edge[,2]==node)
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc)){
        contr=t(t(mutc) %*% MP)
        contr=contr/sum(contr)
        csm=cos_sim(as.vector(mutc),as.vector(prior %*% contr))
        sigbynode[[i]]=list(node=node,contr=contr,csm=csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
      cat("#")
    }
    cat("\n")
  }
  ## Fill in any empty nodes
  for(i in 1:length(PD$pdx$tree_ml$edge.length)){
    if(i>length(sigbynode) || is.null(sigbynode[[i]])){
      sigbynode[[i]]=list(node=PD$pdx$tree_ml$edge[i,2],contr=empty,csm=NA,count=0)
    }
  }
  cat("\n")
  PD$pdx$dat$sigbynode=sigbynode
  PD$pdx$dat$siglist=sigs
  return(PD)
}


get_mouse_details=function(PD,categ){
  ntip=length(PD$pdx$tree_ml$tip.label)
  if(categ=="Shared"){
    PD$pdx$dat$details %>% filter(node > ntip & TYPE=="SNV" & is_localx_excluded==0)
  }else if(categ=="Private"){
    PD$pdx$dat$details %>% filter(node <= ntip & TYPE=="SNV" & is_localx_excluded==0)
  }else if(grepl("^Top",categ)){
    threshold=as.numeric(gsub("Top","",categ))
    nh=nodeHeights(PD$pdx$tree_ml)
    idx=which(nh[,2]<=threshold)
    nodes=unique(PD$pdx$tree_ml$edge[idx,2])
    PD$pdx$dat$details %>% filter(node %in% nodes & TYPE=="SNV" & is_localx_excluded==0)
  }
}

get_sig_color_scheme=function(){
  ## Save this to a file
  if(file.exists("../data/signature_scheme.txt")){
    scheme=read.table("../data/signature_scheme.txt",header = TRUE,stringsAsFactors = FALSE,sep="\t",comment.char = "")
    if(any(is.na(match(c("COL","SIG"),names(scheme))))){
      stop("Invalid signature scheme - need COL and SIG columns")
    }
  }else{
    cols=RColorBrewer::brewer.pal(8,"Set1")
    scheme=data.frame(SIG=c("SBS1","SBSblood","SBS5","SBS18"),COL=cols[1:4])
  }
  scheme
  
}

add_bar_chart=function(df,col2,col.err="black",txttitle=""){
  ## df has the form value,lb,ub,type,label
  # browser()
  xm=abs(par("usr")[1])
  ym=par("usr")[c(3,4)]
  ymm=mean(ym)
  h=diff(ym)*0.5
  w=0.15
  gap=0.1
  inc=0.05
  #col.err="black"
  col.mut=col2[1]
  col.other=col2[2]
  segments(x0 =-xm+0.09*xm,y0=ymm-(h/2),y1=ymm+(h/2),lwd=1)
  for(x in seq(0,1,0.1)){
    segments(x0=-xm+0.05*xm,x1=-xm+0.09*xm,y0=ymm-(h/2)+(h*x),lwd=1)
    text(sprintf("%3.2f",x),x=-xm+0.05*xm,y=ymm-(h/2)+(h*x),pos=2,cex=0.8,xpd=NA,offset=0.1)
  }
  xd=0
  for(TYPE in unique(df$type)){
    tmp=df %>% filter(type==TYPE) 
    rect(xleft=-xm+(xd+0.1)*xm,xright=-xm+(xd+0.1+w)*xm,ybottom = ymm-(h/2),ytop=ymm-(h/2)+(h*tmp$value),col=col.mut,border=NA)
    rect(xleft=-xm+(xd+0.1)*xm,xright=-xm+(xd+0.1+w)*xm,ybottom =ymm-(h/2)+(h*tmp$value),ytop=ymm+(h/2),col=col.other,border=NA)
    segments(x0 = -xm+(0.1+xd+0.5*w)*xm,y0=ymm-(h/2)+h*tmp$lb,
             y1 = ymm-(h/2)+h*tmp$ub,lwd=2,lend=2,col=col.err)
    text(tmp$label,x=-xm+(xd+0.1+0.5*w)*xm,y=ymm-(h/2),pos=1,cex=1)
    xd=xd+inc+w
  }
  text(txttitle,x=-xm+(0.1*w)*xm,y=ymm+1.1*(h/2),pos=4,cex=1.2)
}

plot_signature_annotated_tree_mouse=function(tree,pdx,sig.col.scheme,maxlen=NULL,b.include.csm.text=F){
  if(is.null(pdx$sigbynode)){
    stop("Need to build sigbynode first: see add_sig_to_pdx")
  }
  ##need df of SIG,col
  browser()
  control=list(col.scheme=sig.col.scheme)
  if(!is.null(maxlen)){
    control$maxlen=maxlen
    control$b.include.csm.text=b.include.csm.text
  }
  res=add_annotation(pdx,
                     tree,
                     add_signature_decomp,control=control)
  leg=legend("topleft",sig.col.scheme$SIG,col = sig.col.scheme$COL,pch=15,pt.cex=2)$rect
  tree
}

get_mutation_probability=function(signatures,### Each column is a signature, each row a context
                                  contr ## Each column is a signature, row as sample
){
  signames=colnames(signatures)
  
  ##Check the mutation context ordering
  idx=match(MutationalPatterns:::TRIPLETS_96,rownames(signatures))
  if(length(which(is.na(idx)))>0){
    stop("ERROR:Unexpected rownames")
  }
  signatures=signatures[idx,]
  mutcontext=rownames(signatures)
  contr=sweep(contr,MARGIN = 1,rowSums(contr),"/")## Make sure contr normalised
  signatures==sweep(signatures,MARGIN = 2,colSums(signatures),"/")
  #samples=rownames(contr)
  ###browser()
  do.call("rbind",
          lapply(rownames(contr),function(x) {
            
            probs=signatures %*% diag(contr[x,])
            probs=sweep(probs,MARGIN=1,rowSums(probs),"/")
            colnames(probs)=colnames(signatures)
            cbind(data.frame(
              sample=rep(x,dim(signatures)[1]),
              mutcontext=mutcontext),
              as.data.frame(probs))})
  )
  
}

get_clade_sigs=function(PD,prior,sigmethod="SIGFIT",corecats=c("Mutant","WildType","Trunk"),extracats=c("MutantI","Top50","Top100")){
  stats=NULL
  CTX=NULL
  mp=list()
  for(categ in c(corecats,extracats)){
    cat("Processing",categ,"..\n")
    details=get_details_for_category(PD,categ)
    if(dim(details)[1]==0){
      next
    }
    mutc=with(details,get_mut_matrix(Chrom,Pos,Ref,Alt))
    ctx=as.data.frame(mutc) %>% mutate(count=gr,trin=rownames(.)) %>% mutate(ctx=substr(trin,3,5)) %>% 
      mutate(ctx=ifelse(ctx=="C>T","C>T other",ctx)) %>%
      mutate(ctx=ifelse(grepl("\\[C>T\\]G",trin),"C>T at CpG",ctx)) %>% group_by(ctx) %>% summarise(N=n(),count=sum(count)) %>%
      mutate(donor=PD$patient,category=categ)
    CTX=rbind(CTX,ctx)
    sig=matrix(as.vector(mutc),ncol=1)
    if(sigmethod=="SIGFIT"){
      sf=fit_signatures(t(sig),t(prior))
      exposures=retrieve_pars(sf,par="exposures")
      decomp=list(contribution=t(as.matrix(exposures$mean)))
      colnames(decomp$contribution)=categ
      decomp$csm=cos_sim(as.vector(retrieve_pars(sf,par="reconstructions")$mean),as.vector(sig))
      ## Following is a bit ugly - but emulates the structure emitted by MP method below.
      qq=rbind(do.call("rbind",exposures) %>% (function(x){rownames(x)=c("mean","2.5%","97.5%");x}),
               total=rep(sum(mutc),dim(exposures$mean)[2]))
    }else{
      ## Mutational patterns fit
      decomp=MutationalPatterns:::fit_to_signatures(sig,prior)
      decomp$contribution[,1]=decomp$contribution[,1]/sum(decomp$contribution[,1])
      decomp$csm=cos_sim(decomp$reconstructed[,1],as.vector(mutc))
      colnames(decomp$contribution)=categ
      ## Next do bootstrapping to get CIs.
      bs=fit_to_signatures_bootstrapped(sig,prior,verbose = FALSE,method = "regular") %>% (function(probs) sweep(probs,MARGIN=1,rowSums(probs),"/"))
      qq=sapply(colnames(bs),function(x) {qq=quantile(bs[,x],c(0.025,0.5,0.975));c(qq,mean=mean(bs[,x],),sd=sd(bs[,x]),total=sum(mutc))})
      
    }
    stats=rbind(stats,as.data.frame(qq) %>% mutate(field=rownames(.),category=categ,donor=PD$patient,rec_csm=decomp$csm) )
    MP=as.matrix(get_mutation_probability(prior,t(decomp$contribution))[,-(1:2)])
    mp[[categ]]=MP
    cat("\n")
  }
  list(corecats=corecats,extracats=extracats,sigstats=stats %>% pivot_longer(cols=colnames(prior),names_to = "signature",values_to = "proportion") %>%
         pivot_wider(names_from="field",values_from = "proportion"),
       ctx=CTX,MP=mp)
  
}

plot_sig_tree_posterior=function(PD,prior,b.ultra=FALSE,txtitle=PD$patient,b.add.sigs=TRUE,focal.signature="SBS1"){
  info=PD$pdx$dat$info
  if(is.null(info)){
    stop('prior to running please do: PD$pdx$dat$info=get_clade_sigs(PD,prior=pcawg[, c("SBS1","SBSblood","SBS18")])')
    ###info=get_clade_sigs_mouse(PD,prior)
  }
  if(is.null(PD$pdx$dat$sigbynode)){
    PD2=add_v2_sig_to_pdx(PD,prior=prior,info = info)
  }else{
    PD2=PD
  }
  colorscheme=get_sig_color_scheme() %>% filter(SIG %in% colnames(prior))
  if(b.ultra){
    PD2$pdx$tree_ml=PD2$fit$poisson_tree$altmodel$ultratree
  }
  tree=plot_tree(PD2$pdx$tree_ml,cex.label = 0,left.margin.prop = 0.3)#,mar = c(1, 2, 1, 1) + 0.1)
  title(txtitle);
  tree=plot_signature_annotated_tree_v2(tree,PD2$pdx$dat,colorscheme,maxlen=1e6,b.include.csm.text = FALSE)
  ## Convert data to form ## For 
  if(!b.add.sigs){
    agg1=info$ctx %>% group_by(category) %>% summarise(total=sum(count))
    df=info$ctx %>% filter(ctx=="C>T at CpG") %>% left_join(agg1)
    bounds=t(sapply(1:dim(df)[1],function(i) binom.test(df$count[i],df$total[i])$conf.int))
    colnames(bounds)=c("lb","ub")
    df=cbind(df %>% mutate(value=count/total),as.data.frame(bounds))
    df=df %>% left_join(data.frame(category=c("Shared","Private"),label=c("Shared","Private")))
    add_bar_chart(df %>% mutate(type=category),col2 =c("darkorange","white"),txttitle="C>T @ CpG")
  }else{
    df=info$sigstats %>% mutate(lb=`2.5%`,ub=`97.5%`, value=mean) %>% filter(signature==focal.signature) %>% mutate(label=category)
    #  left_join(data.frame(category=c("Shared","Private"),label=c("Shared","Private")))
   
    add_bar_chart(df %>% mutate(type=category),c(colorscheme %>% filter(SIG==focal.signature) %>% pull(COL),"white"))
  }
  tree
}

## Adds grouped category level signature to the tree 
add_v2_sig_to_pdx=function(PD,prior,info){
  sigbynode=list()
  sigs=colnames(prior)
  empty=matrix(NA,ncol=1,nrow=length(sigs))
  rownames(empty)=sigs
  for(categ in info$corecats){
    cat("Processing",categ,"..\n")
    details=get_details_for_category(PD,categ)
    MP=info$MP[[categ]]
    for(node in unique(details$node)){
      i=which(PD$pdx$tree_ml$edge[,2]==node)
      mutc=get_mut_matrix_for_node(PD$pdx,node)
      if(!is.null(mutc)){
        contr=t(t(mutc) %*% MP)
        contr=contr/sum(contr)
        csm=cos_sim(as.vector(mutc),as.vector(prior %*% contr))
        sigbynode[[i]]=list(node=node,contr=contr,csm=csm,count=sum(mutc))
      }else{
        sigbynode[[i]]=list(node=node,contr=empty,csm=NA,count=sum(mutc))
      }
      cat("#")
    }
    cat("\n")
  }
  ## Fill in any empty nodes
  for(i in 1:length(PD$pdx$tree_ml$edge.length)){
    if(i>length(sigbynode) || is.null(sigbynode[[i]])){
      sigbynode[[i]]=list(node=PD$pdx$tree_ml$edge[i,2],contr=empty,csm=NA,count=0)
    }
  }
  cat("\n")
  PD$pdx$dat$sigbynode=sigbynode
  PD$pdx$dat$siglist=sigs
  return(PD)
}

get_details_for_category=function(PD,categ){
  ## 
  tnode=PD$nodes %>% filter(driver=="t(9;22)(q34;q11)") %>% pull(node)
  if(categ %in% c("Mutant","Trunk","MutantI") && length(tnode)==0){
    cat("no mutant node defined - returning empty matrix\n")
    return(PD$pdx$dat$details %>% head(0))
  }
  ntip=length(PD$pdx$tree_ml$tip.label)
  if(categ=="All"){
    PD$pdx$dat$details %>% filter(TYPE=="SNV" & is_localx_excluded==0)
  }else if(categ=="Trunk"){
    PD$pdx$dat$details %>% filter(node==tnode & TYPE=="SNV" & is_localx_excluded==0)
  }else if(categ=="WildType"){
    PD$pdx$dat$details %>% filter((!(node %in% c(tnode,get_all_node_children(tnode,PD$pdx$tree_ml))) & TYPE=="SNV" & is_localx_excluded==0))
  }else if(categ=="Mutant"){
    PD$pdx$dat$details %>% filter((node %in% get_all_node_children(tnode,PD$pdx$tree_ml) )& TYPE=="SNV" & is_localx_excluded==0 )
  }else if(categ=="MutantI"){
    PD$pdx$dat$details %>% filter((node %in% setdiff(get_all_node_children(tnode,PD$pdx$tree_ml),1:length(PD$pdx$tree_ml$tip.label) ) ) & TYPE=="SNV" & is_localx_excluded==0)
  }else if(categ=="Shared"){
    PD$pdx$dat$details %>% filter(node > ntip & TYPE=="SNV" & is_localx_excluded==0)
  }else if(categ=="Private"){
    PD$pdx$dat$details %>% filter(node <= ntip & TYPE=="SNV" & is_localx_excluded==0)
  }else if(grepl("^Top",categ)){
    threshold=as.numeric(gsub("Top","",categ))
    nh=nodeHeights(PD$pdx$tree_ml)
    idx=which(nh[,2]<=threshold)
    nodes=unique(PD$pdx$tree_ml$edge[idx,2])
    PD$pdx$dat$details %>% filter(node %in% nodes & TYPE=="SNV" & is_localx_excluded==0)
  }else{
    stop("invalid tree category")
  }
}

plot_signature_annotated_tree_v2=function(tree,pdx,sig.col.scheme,maxlen=NULL,b.include.csm.text=TRUE){
  if(is.null(pdx$sigbynode)){
    stop("Need to build sigbynode first: see add_sig_to_pdx")
  }
  ##need df of SIG,col
  control=list(col.scheme=sig.col.scheme)
  if(!is.null(maxlen)){
    control$maxlen=maxlen
    control$b.include.csm.text=b.include.csm.text
  }
  res=add_annotation(pdx,
                     tree,
                     add_signature_decomp,control=control)
  leg=legend("topleft",sig.col.scheme$SIG,col = sig.col.scheme$COL,pch=15,pt.cex=2)$rect
  tree
}

get_sig_burdens=function(PD,priors=pcawg[,c("SBS1","SBSblood")],ultramodel="nullmodel",corecats=c("Shared","Private")){
  if(is.null(PD$pdx$dat$sigbynode)){
    PD$pdx$dat$info=get_clade_sigs(PD,prior=priors,extracats = c(),corecats = corecats)
    PD=add_v2_sig_to_pdx(PD,prior=priors,info = PD$pdx$dat$info)
  }
  ## Get burden as tree tip height (tree is adjusted and corresponds to SNV burden)
  nh=nodeHeights(PD$pdx$tree_ml)
  tmp=data.frame(tip.label=PD$pdx$tree_ml$tip.label,
                 nsub_adj=nh[match(1:length(PD$pdx$tree_ml$tip.label),PD$pdx$tree_ml$edge[,2]),2]
                 )
  # browser()
  # nht=nodeHeights(PD$fit$poisson_tree[[ultramodel]]$ultratree)
  nht=nh
  traj=data.frame(start=nht[,1],mid=0.5*(nht[,1]+nht[,2]),end=nht[,2],nsub_adj_start=nh[,1],nsub_adj_mid=0.5*(nh[,1]+nh[,2]),nsub_adj_end=nh[,2],
                  parent=PD$pdx$tree_ml$edge[,1],child=PD$pdx$tree_ml$edge[,2])
  for(SIG in colnames(priors)){
    frac=sapply(1:length(PD$pdx$tree_ml$edge.length),function(i) if(i<=length(PD$pdx$dat$sigbynode)){PD$pdx$dat$sigbynode[[i]]$contr[SIG,]}else{0})
    frac=ifelse(is.na(frac),0,frac)
    tree=PD$pdx$tree_ml
    tree$edge.length=tree$edge.length*frac
    nh=nodeHeights(tree)
    ht=nh[match(1:length(tree$tip.label),tree$edge[,2]),2]
    #tmp2=data.frame(tip.label=tree$tip.label)
    tmp[[paste0("nsub_adj_",SIG)]]=ht
    traj[[paste0(SIG,"_start")]]=nh[,1]
    traj[[paste0(SIG,"_mid")]]=0.5*(nh[,1]+nh[,2])
    traj[[paste0(SIG,"_end")]]=nh[,2]
  }
  
  df=PD$pdx$agedf %>% filter(tip.label!="zeros")  %>% left_join(tmp) %>%
    left_join(PD$pdx$cfg %>% mutate(tip.label=SHORT_LABEL,colony=LABEL) %>% dplyr::select(tip.label,colony)) %>% mutate(donor=PD$patient)
  list(traj=traj,df=df)
}

make_ultrametric_PD_chiraag=function(md,generation_time=1){ #Current implementation leaves in the outlier zeros branch..
  tree<-md$pdx$tree_ml
  s<-md$patient
  tree$edge.length<-(tree$edge.length+.000001)
  germline_profile=paste(ifelse(md$pdx$df$samples=="zeros",0,1),collapse="")
  
  tree.ut<-make.ultrametric.tree(tree) 
  
  tree.ut$edge.length=tree.ut$edge.length/max(nodeHeights(tree.ut)) 
  generation_time=generation_time;

    if(grepl("MD7180|MD7181|MD7182",s)){
    age=132+3;  
    age=age/generation_time
  }else{
    age=12+3;
    age=age/generation_time
  }
  
  
  tree.ut$edge.length=tree.ut$edge.length*age
  
  md$pdx$tree_ml <- tree.ut
  md
}


