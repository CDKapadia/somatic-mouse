## Code Supporting analysis of MPP/HSC ontogeny
# See supplementary_note2.Rmd for example usage
# Author(s) : Nick Williams and Chiraag Kapadia.  
#             Nangalia Group CASM Wellcome Sanger Institute
require("Matrix")
require("matrixcalc")
require("matrixStats")
require("parallel")
require("bbmle")
require("rstan")

### Utility function  for annotating tree tip phenotypes based on info in a tab delimited file "hsc_mpp.txt".
annotate_celltype=function(tree,mouse){
  df=read.table("hsc_mpp.txt",head=TRUE,sep="\t")
  dff=data.frame(tip.label=tree$tip.label) %>% left_join(df %>% filter(grepl(mouse,LABEL)),by=c("tip.label"="SHORT_LABEL"))
  cols=c(HSC="blue",MPP="red")
  tree$CELLTYPE=dff$CELLTYPE
  tree$tip.color=cols[tree$CELLTYPE]
  tree
}

## Get 2 epoch transition probability data frame
get_prior_probs_full_2epoch=function(ph, #emb_to_hsc
                                     pm, #emb_to_mpp
                                     p_m_h, #mpp_to_hsc
                                     p_h_m, #hsc_to_mpp
                                     mt.start=10){
  data.frame(
    start=c(0,mt.start),
    emb_to_mpp=c(0,pm),
    emb_to_hsc=c(0,ph),
    hsc_to_mpp=c(0,p_h_m),
    mpp_to_hsc=c(0,p_m_h)
  )
}
## Converts data.frame of per epoch rates into a list of transition matrices. Use get_prior_probs_full_2epoch
## to get example format for prior.probs
getTransitionMatrices=function(prior.probs,epsilon=1e-12){
  P=lapply(1:length(prior.probs$start),function(i){
    z=with(prior.probs[i,],matrix(c(1-(hsc_to_mpp+epsilon),hsc_to_mpp,epsilon,
                                    mpp_to_hsc,1-(mpp_to_hsc+epsilon),epsilon,
                                    emb_to_hsc+epsilon,emb_to_mpp+epsilon,1-(emb_to_hsc+emb_to_mpp+2*epsilon)),
                                  byrow = TRUE,ncol=3))
    colnames(z)=c("HSC","MPP","EMB")
    rownames(z)=colnames(z)
    z
  })
  
  out=list(epochs=c(prior.probs$start,Inf),P=P)
  class(out)="TransitionMatrices"
  out
}

## Gets the discrete time transition matrix between start and end based on
## the epoch timing and transition matrices in transitionMatrices (typically obtained using getTransitionMatrices)
get_Discrete_P=function(start,end,transitionMatrices){
  P=diag(1,3)
  tm=transitionMatrices
  epochs=tm$epochs
  nepoch=length(tm$P)
  for(k in 1:nepoch){
    if(epochs[k]>end){
      ## starts after the branch ends nothing required
    }else{
      ##epoch starts before or in branch life
      if(epochs[k]<=start){
        if(epochs[k+1]>start){
          ## need to include
          P=P %*% matrix.power(tm$P[[k]], (min(end,epochs[k+1])-start))
        }## otherwise no action
      }else{
        ## epoch starts in range..
        P=P %*% matrix.power(tm$P[[k]],(min(end,epochs[k+1])-epochs[k]))
      }
    }
  }
  colnames(P)=c("HSC","MPP","EMB")
  rownames(P)=colnames(P)
  P
}
## Plots transition probabilities encoded in the provided transtionMatrices (see getTransitionMatrices)
plotTransitionMatrices=function(tm,label="",xmax=150){
  PS=t(sapply(0:xmax,function(i) get_Discrete_P(0,i,tm)[3,])) %>% as.data.frame() %>% mutate(T=0:xmax)
  ggplot(PS %>% pivot_longer(cols=colnames(tm$P[[1]]),names_to = "TYPE",values_to = "P"),
         aes(x=T,y=P,col=TYPE))+geom_line()+theme_bw()+xlim(c(0,xmax))+xlab("Molecular Time (#Mutations)")+ylab("Probability")+ggtitle(paste0(label," : Cell Type Probability vs Molecular Time"))+theme(legend.title = element_blank())+geom_vline(xintercept = tm$epochs,size=0.1)
}

#' @title 
#' Infers the most likely set of discrete states at the nodes of a phylogeny.
#' 
#' For each of the tips we have a cell type designation
#' 
#' @description 
#' This is the top level function  for inferring "hidden"  states for all nodes in a phylogenetic tree where the only
#' "emitted" data is at the tips (in our case phenotype of the tips).
#' @param tree APE phylo class.   
#' @param pheno phenotype vector aligned with tree$tip.label
#' @param numStates number of hidden states (k in Durand et al)
#' @param emissionFunction should take a pheno object and sample number and return an n x numState matrix (n is number of tips)
#' @param prior.rates data frame encoding the multi-epoch model with per epoch transition probabilities
#' @param root.prior vector of prior probability for starting state in the order HSC,MPP,EMB. Always c(0,0,1)  in our model.
#' @return list of results (tree - marked up tree, S.root - inferred state of root (zygote),changes - data.frame detailing changes in state)
viterbi.tree=function(tree,
                      pheno,
                      numStates=3,
                      emissionFunction=emissionFunction3state,
                      prior.rates,
                      root.prior=rep(1/numStates,numStates),
                      type="CONTINUOUS"
){
  if(is.null(tree$P)){
    if(type=="CONTINUOUS"){
      tree=addTransitionProbabilities(tree,prior.rates)
    }else{
      tree=addDiscreteTransitionProbabilities(tree,getTransitionMatrices(prior.rates))
    }
  }
  if(is.null(tree$children)){
    tree$children=lapply(1:max(tree$edge[,2]),get_node_children,tree=tree)
  }
  if(dim(tree$P[[1]])[1]!=numStates){
    stop("mismatch between specified number of states and provied transitionfunction")
  }
  
  if(dim(emissionFunction(pheno,1))[2]!=numStates){
    stop("mismatch between specified number of states and provied emissionfunction")
  }
  tree$delta=list()###matrix(NA,ncol=numStates,nrow=m)
  tree$node.lookup=match(1:max(tree$edge[,2]),tree$edge[,2])
  root.delta=rep(NA,numStates)
  # ensure pheno is in the correct order
  idx=match(tree$tip.label,pheno$tip.label)
  n=length(tree$tip.label)
  m=1##dim(pheno$MF)[1]
  tree$nsites=m
  tree$nstates=numStates
  idx.tip=tree$node.lookup[1:n]#match(1:n,tree$edge[,2])
  
  tree$delta[idx.tip]=lapply(1:length(idx.tip),
                             function(i) t(emissionFunction(pheno,i )))
  ## Check for missing values
  if(any(sapply(tree$delta[idx.tip],function(x) any(is.na(x))))){
    qstop("Deltas are NA")
  }
  if(any(sapply(tree$delta[idx.tip],function(x) any(colSums(x==0)>0)))){
    warning("Zeros in emmision function")
  }
  tree$psi=lapply(1:length(tree$edge.length),function(x){list()})
  root=length(tree$tip.label)+1
  children=tree$children[[root]]
  delta.root=matrix(1,nrow=numStates,ncol=m)
  psi.root=list()
  k=1
  tree$stabilityfactorlog=0
  for(cnode in children){
    psi.root[[k]]=list(child=cnode,most_likely_state=matrix(NA,nrow=numStates,ncol=m))
    tree=update.deltas(tree,cnode)
    idx=tree$node.lookup[cnode]
    P=tree$P[[idx]]
    for(j in 1:numStates){
      z=tree$delta[[idx]]*P[j,]##  This works because recycling works by column
      kmax=max.col(t(z),ties.method = "first")##  The ties condition should apply infrequently
      zmax=colMaxs(z)
      if(any(is.na(zmax)) || any(zmax==0)){
        qstop("zmax=0...")
      }
      delta.root[j,]=delta.root[j,]*zmax
      psi.root[[k]]$most_likely_state[j,]=kmax
    }
    k=k+1
  }
  ##set most likely states
  #cat("Root Probs",delta.root*root.prior,"\n")
  p=max.col(t(delta.root*root.prior),ties.method = "first")
  #cat("log P: raw:",log(max(delta.root*root.prior))," + stability=",log(max(delta.root*root.prior))+tree$stabilityfactorlog,"\n")
  LL=log(max(delta.root*root.prior))+tree$stabilityfactorlog
  
  tree$S=matrix(NA,nrow=length(tree$edge.length),ncol=m)
  for(entry in psi.root){
    tree=update.states(tree,entry$child,entry$most_likely_state[tree$nstates*((1:length(p))-1)+p])
  }
  S.root=p ##  Zygote genotype
  time=tree$edge.length
  durationByParentState=matrix(0,nrow=tree$nstates,ncol=tree$nsites)
  
  idx=match(tree$edge[,1],c(tree$edge[,2],root))
  S=rbind(tree$S,S.root)
  Sparent=S[idx,,drop=FALSE]
  Schild=tree$S
  idx2=which(Sparent!=Schild)
  N=length(idx)
  
  time=tree$edge.length
  durationByParentState=do.call("rbind",lapply(1:numStates,function(i) colSums((Sparent==i)*time)))
  genotype.summary=do.call("rbind",lapply(1:numStates,function(i) colSums(tree$S[idx.tip,,drop=FALSE]==i)))
  changes=data.frame(i=((idx2-1) %/% N)+1,
                     s.parent=Sparent[idx2]-1,
                     s.child=Schild[idx2]-1,
                     node=tree$edge[((idx2-1) %% N)+1,2])
  GENO=t(tree$S[idx.tip,])
  list(changes=changes,
       allbranches=data.frame(node=tree$edge[,2],s.parent=Sparent-1,s.child=Schild-1),
       S.root=S.root,
       durationByParentState=t(durationByParentState),
       genotype.summary=t(genotype.summary),
       geno=GENO,
       tree=tree,
       summary=data.frame(s.root=S.root-1,d0=durationByParentState[1],d1=durationByParentState[2]),
       prior.rates=prior.rates,
       LL=LL
  )
}
update.states=function(tree,node,state){
  idx=tree$node.lookup[node]#match(node,tree$edge[,2])
  tree$S[idx,]=state
  if(length(tree$psi[[idx]])>0){
    for(entry in tree$psi[[idx]]){
      tree=update.states(tree,entry$child,entry$most_likely_state[tree$nstates*((1:length(state))-1)+state])
    }
  }
  tree
}

update.deltas=function(tree,node){
  children=tree$children[[node]]#get_node_children(node,tree);
  idx1=tree$node.lookup[node]#match(node,tree$edge[,2])
  if(!is.null(tree$delta[[idx1]])){
    return(tree)
  }else{
    tree$delta[[idx1]]=matrix(1,nrow=tree$nstates,ncol=tree$nsites)
    tree$psi[[idx1]]=list()
    k=1
    for(cnode in children){
      tree$psi[[idx1]][[k]]=list(child=cnode,most_likely_state=matrix(NA,nrow=tree$nstates,ncol=tree$nsites))
      tree=update.deltas(tree,cnode)
      idx=tree$node.lookup[cnode]##match(cnode,tree$edge[,2])
      P=tree$P[[idx]]#transitionFunction(tree$edge.length[idx])
      ##  Do some rescaling to stop going to zero
      zmm=rep(0,tree$nsites)
      for(j in 1:tree$nstates){
        z=tree$delta[[idx]]*P[j,]##  This works because recycling works by column
        kmax=max.col(t(z),ties.method = "first")
        zmax=colMaxs(z)
        if(any(is.na(zmax)) || any(zmax==0)){
          qstop("zmax=0...")
        }
        zmm=ifelse(zmax>zmm,zmax,zmm)
        tree$delta[[idx1]][j,]=tree$delta[[idx1]][j,]*zmax
        
        tree$psi[[idx1]][[k]]$most_likely_state[j,]=kmax
      }
      ## Numerical Stability ### START
      
      tree$delta[[idx1]]=tree$delta[[idx1]]/rep(zmm,each=tree$nstates)
      tree$stabilityfactorlog=tree$stabilityfactorlog+log(zmm)
      ## Numerical Stability ### END
      k=k+1
    }
  }
  tree
}


addDiscreteTransitionProbabilities=function(tree,transitionMatrices){
  nh=nodeHeights(tree)
  ##
  tree$P=lapply(1:length(tree$edge.length),function(i) get_Discrete_P(nh[i,1],nh[i,2],transitionMatrices))
  tree
  
}



HMT_STAN_3TREE="functions{
  real hmt_LL(int N, int M,int C,matrix transitionMatrix,  matrix deltaIn,int [,] children, int [] visit_order,int [] bl,vector root_prior){
    //matrix [N+1,3] delta;
    matrix [N+1,3] delta=deltaIn;
    matrix [3,3] P;
    int k=0;
    int child=0;
    vector [3] tots;
    real stability_factor;
    real stability_factor_log;
    stability_factor_log=0.0;
    for(i in 1:((N-M)+1)){
      k=visit_order[i];
      delta[k]=rep_row_vector(1,3);
      for( j in 1:C){
          child=children[i,j];
          if(child<0) break;
          P=matrix_power(transitionMatrix,bl[child]);
          tots=P*delta[child]';//(row(delta,child)');
          stability_factor=max(tots);
          delta[k]=(delta[k].*tots')/stability_factor;
          stability_factor_log=stability_factor_log+log(stability_factor);
      }
    }
    //Finally do the root. 
    return log(dot_product(delta[N+1]',root_prior))+stability_factor_log;
    
  }
}

data{

  vector [3] root_prior; // Start in EMB state (0,0,1)

  int N1; //num branches
  int C1; // max number of children (largest polytomy)
  int M1; // number of tips
  int children1[N1-M1+1,C1]; //IDX of children. On each row children after and including first negative.
  int visit_order1[N1-M1+1];
  int bl1[N1]; // Length of each branch.  
  matrix [N1+1,3] delta1; // probability of observed state given true underlying state (3 state)
  
  int N2; //num branches
  int C2; // max number of children (largest polytomy)
  int M2; // number of tips
  int children2[N2-M2+1,C2]; //IDX of children. On each row children after and including first negative.
  int visit_order2[N2-M2+1];
  int bl2[N2]; // Length of each branch.  
  matrix [N2+1,3] delta2; // probability of observed state given true underlying state (3 state)
  
  int N3; //num branches
  int C3; // max number of children (largest polytomy)
  int M3; // number of tips
  int children3[N3-M3+1,C3]; //IDX of children. On each row children after and including first negative.
  int visit_order3[N3-M3+1];
  int bl3[N3]; // Length of each branch.  
  matrix [N3+1,3] delta3; // probability of observed state given true underlying state (3 state)
}

parameters {
  real<lower=0.000001,upper=0.999999> pa;
  real<lower=0.000001,upper=0.999999> po;
  real<lower=0.000001,upper=0.999999> phscmpp;
  real<lower=0.000001,upper=0.999999> pmpphsc;
}

model {
  matrix [3,3] TM;
  real PA;
  real PO;
  real PHM;
  real PMH;
  pa~uniform(0,1);
  po~uniform(0,1);
  phscmpp~uniform(0,0.5);
  pmpphsc~uniform(0,0.5);
  PA=pa;
  PO=po;
  PHM=phscmpp;
  PMH=pmpphsc;
  TM[1,3]=0;
  TM[2,3]=0;
  TM[1,1]=1-PHM;
  TM[1,2]=PHM;
  TM[2,1]=PMH;
  TM[2,2]=1-PMH;
  TM[3,1]=PA*PO;
  TM[3,2]=PO*(1-PA);
  TM[3,3]=1-PO;
  ///print(\"target:\",target());
  target+=hmt_LL(N1, M1,C1,TM, delta1,children1, visit_order1,bl1,root_prior);
  target+=hmt_LL(N2, M2,C2,TM, delta2,children2, visit_order2,bl2,root_prior);
  target+=hmt_LL(N3, M3,C3,TM, delta3,children3, visit_order3,bl3,root_prior);
}

generated quantities {
  real P_EMB_HSC=pa*po;
  real P_EMB_MPP=po*(1-pa);
  
}
"



HMT_STAN_LL="functions{
  real hmt_LL(int N, int M,int C,matrix transitionMatrix,  matrix deltaIn,int [,] children, int [] visit_order,int [] bl,vector root_prior){
    //matrix [N+1,3] delta;
    matrix [N+1,3] delta=deltaIn;
    matrix [3,3] P;
    int k=0;
    int child=0;
    vector [3] tots;
    real stability_factor;
    real stability_factor_log;
    stability_factor_log=0.0;
    for(i in 1:((N-M)+1)){
      k=visit_order[i];
      delta[k]=rep_row_vector(1,3);
      for( j in 1:C){
        child=children[i,j];
        if(child<0) break;
        P=matrix_power(transitionMatrix,bl[child]);
        //print(transitionMatrix,P,bl[child],delta[child]);
        tots=P*delta[child]';//(row(delta,child)');
        stability_factor=max(tots);
        delta[k]=(delta[k].*tots')/stability_factor;
        stability_factor_log=stability_factor_log+log(stability_factor);
      }
    }
    //Finally do the root. 
    //print(delta[N+1],root_prior,stability_factor_log);
    return log(dot_product(delta[N+1]',root_prior))+stability_factor_log;

      }
    }
    
    data{
      
      vector [3] root_prior; // Start in EMB state (0,0,1)
      
      int N1; //num branches
      int C1; // max number of children (largest polytomy)
      int M1; // number of tips
      int children1[N1-M1+1,C1]; //IDX of children. On each row children after and including first negative.
      int visit_order1[N1-M1+1];
      int bl1[N1]; // Length of each branch.  
      matrix [N1+1,3] delta1; // probability of observed state given true underlying state (3 state)
      matrix [3,3] TM; // Transition Matrix
    }
    
    generated quantities {
      real LL;
      LL=hmt_LL(N1, M1,C1,TM, delta1,children1, visit_order1,bl1,root_prior);
    }

    
"
get_model=function(modelString){
  fp=sprintf("./%s.RDS",modelString)
  if(file.exists(fp)){
    readRDS(fp)
  }else{
    out=stan_model(model_code=get(modelString),model_name = modelString)
    saveRDS(out,fp)
    out
  }
}
HMT_STAN_3TREE_MODEL=get_model("HMT_STAN_3TREE")
##HMT_STAN_SINGLE_TREE_MODEL=get_model("HMT_STAN_SINGLE_TREE")
HMT_STAN_LL_MODEL=get_model("HMT_STAN_LL")

### This defines the order of in which nodes should visited so that 
## parent nodes are visited after all their child nodes. This simplifies
## the STAN code as no recursion in STAN is then required.
add_visit_vector=function(tree){
  root=length(tree$tip.label)+1
  tree$control=1
  tree$visits=c()
  tree=get_visits(tree,root)
  tree$visits=setdiff(tree$visits,1:length(tree$tip.label))
  tree$visits=match(tree$visits,tree$edge[,2])
  ## 
  children=matrix(-1,nrow=length(tree$visits)+1,ncol=max(table(tree$edge[,1])))
  for(i in 1:length(tree$visits)){
    idx=which(tree$edge[,1]==tree$edge[tree$visits[i],2])
    children[i,1:length(idx)]=idx
  }
  ## Add in the root
  tree$visits=c(tree$visits,length(tree$edge.length)+1)
  idx=which(tree$edge[,1]==root)
  children[length(tree$visits),1:length(idx)]=idx
  tree$children=children
  tree
}

get_visits=function(tree,node){
  if(tree$control>1000){
    stop("too deep!")
  }
  tree$control=tree$control+1
  #cat(node,":\n")
  #cat(paste(tree$visits,collapse=","),"\n")
  children=get_node_children(node,tree)
  for(node in children){
    tree$visits=c(node,tree$visits)
    tree=get_visits(tree,node)
  }
  
  tree
}

### Run the STAN model across three trees.  
## Note this is hard coded for three trees. To generalise this to N trees will need to code up an additional array dimension for STAN 
## data elements.
run.upward.discrete.3mice.hmt=function(trees,mt.start=10,epsilon=1e-12,niter=5000,chains=4,cores=4){
  if(length(trees)!=3){
    stop("Only use for groups of 3 mice!")
  }
  dats=lapply(trees,function(tree){
  dff=data.frame(tip.label=tree$tip.label,CELLTYPE=tree$CELLTYPE)
  tree$edge.length=round(tree$edge.length)
  ## First trim the first mt.start mutations of the tree.. 
  nh=nodeHeights(tree)
  tree$edge.length=ifelse(nh[,1]<=mt.start,
                          ifelse(nh[,2]>=mt.start,tree$edge.length+nh[,1]-mt.start,0),
                          tree$edge.length)
  N=length(tree$edge.length)
  M=length(tree$tip.label)
  tree=add_visit_vector(tree)
  
  tree$delta=matrix(1,ncol=3,nrow=N+1)
  idx=match(1:M,tree$edge[,2])
  tree$delta[idx,]=t(sapply(1:M,function(i) {
    if(dff$CELLTYPE[i]=="HSC"){ 
      c(1-epsilon,epsilon*0.5,epsilon*0.5)
    }else{
      c(epsilon*0.5,1-epsilon,epsilon*0.5)
    }
  }))
  ## prepare the data for STAN
  list(N=N,
           C=ncol(tree$children),## Max number of children. Non-sparse representation!
           M=M,
           children=tree$children,#IDX of children in edge matrix . On each row children after and including first negative.
           visit_order=tree$visits, ## Proceeding in this way goes up the tree.
           bl=tree$edge.length, ##Length of each branch.  
           delta=tree$delta
  )
  })
  DAT=do.call("c",lapply(1:3,function(i) {
    fields=names(dats[[i]])
    for(field in fields){
      dats[[i]][[sprintf("%s%d",field,i)]]=dats[[i]][[field]]
    }
    dats[[i]][sprintf("%s%d",fields,i)]
    })
  )
  DAT$root_prior=c(0,0,1)
  LL=sampling(HMT_STAN_3TREE_MODEL,data=DAT,iter=niter,chains=chains,cores=cores,control=list(adapt_delta=0.95,max_treedepth=14))
}

###bbmle based MLE implementation
invlogit=plogis
logit=qlogis

LOGLIK_VIT_FULL=function(a,po,hm,mh,dat){
  trees=dat$trees
  a=invlogit(a)
  po=invlogit(po)
  pp=get_prior_probs_full_2epoch(ph=a*po,pm=(1-a)*po,p_m_h=invlogit(mh),p_h_m=invlogit(hm),
                                 mt.start = dat$mt.start)
  tm=getTransitionMatrices(pp)
  LLL=sapply(names(trees),function(id){
    tree=trees[[id]]
    tree$edge.length=round(tree$edge.length)
    ##LL=run.upward.discrete.mle(tree,mouse = id,prior.probs = pp)  Pure R implementation.
    LL=get.upward.discrete.LL(tree,prior.probs = pp,mt.start = dat$mt.start)
  }
  )
  -sum(LLL)
}


##Use STAN to calculate the log likelihood ( A pure R version of this function is also implemented)
get.upward.discrete.LL=function(tree,prior.probs,mt.start=10,epsilon=1e-12,root.prior=c(0,0,1)){
  dff=data.frame(tip.label=tree$tip.label,CELLTYPE=tree$CELLTYPE)
  tree$edge.length=round(tree$edge.length)
  ## First trim the first mt.start mutations of the tree.. 
  nh=nodeHeights(tree)
  #browser()
  tree$edge.length=ifelse(nh[,1]<=mt.start,
                          ifelse(nh[,2]>=mt.start,tree$edge.length+nh[,1]-mt.start,0),
                          tree$edge.length)
  
  N=length(tree$edge.length)
  M=length(tree$tip.label)
  tree=add_visit_vector(tree)
  
  tree$delta=matrix(1,ncol=3,nrow=N+1)
  idx=match(1:M,tree$edge[,2])
  tree$delta[idx,]=t(sapply(1:M,function(i) {
    if(dff$CELLTYPE[i]=="HSC"){ 
      c(1-epsilon,epsilon*0.5,epsilon*0.5)
    }else{
      c(epsilon*0.5,1-epsilon,epsilon*0.5)
    }
  }))
  tm=getTransitionMatrices(prior.probs)
  ## prepare the data for STAN
  dat=list(N1=N,
           C1=ncol(tree$children),## Max number of children. Non-sparse representation!
           M1=M,
           children1=tree$children,#IDX of children in edge matrix . On each row children after and including first negative.
           visit_order1=tree$visits, ## Proceeding in this way goes up the tree.
           bl1=tree$edge.length, ##Length of each branch.  
           delta1=tree$delta,  ##Delta initialised with tips. probability of observed state conditional on underlying state (3 state)
           root_prior=c(0,0,1), #Start in EMB state (0,0,1)
           TM=tm$P[[2]]
  )
  ## Use a the STAN based likelihood. Have equivalent in pure R if necessary.
  LL=sampling(HMT_STAN_LL_MODEL,data=dat,iter=1,chains=1,cores=1,algorithm="Fixed_param",verbose=FALSE,show_message=FALSE,refresh=0)
  LL=summary(LL)$summary["LL","mean"]
  if(is.infinite(LL)){
    LL=0.1*.Machine$double.xmax*sign(LL)
  }
  LL
}

## Only consider tip states of HSC or MPP 
emissionFunction3state=function(pheno,i){
  epsilon=1e-12
  out=matrix(NA,ncol=3,nrow=1)
  if(pheno$CELLTYPE[i]=="HSC"){
    out[1,]=c(1-epsilon,0.5*epsilon,0.5*epsilon)
  }else{
    out[1,]=c(0.5*epsilon,1-epsilon,0.5*epsilon)
  }
  out
}
###MLE Reporting Function.
report_results=function(result.list,b.verbose=FALSE,a=1){
  allres=sapply(result.list,function(x) c(invlogit(coef(x)),LL=logLik(x)))
  idx=which.max(allres["LL",])
  max.ll=max(allres["LL",])
  if(b.verbose){
    print(allres)
    cat("close matches...\n")
    close.idx=which(allres["LL",]>max.ll-0.1)
    print(allres[,close.idx] %>% (function(x){colnames(x)=sprintf("i%d",close.idx);x}) )
    print(allres[,idx])
    
  }
  ##cat("IDX=",idx,":Num close matches=",length(which(allres["LL",]>max.ll-0.1)),"\n")
  res=allres[,idx]
  if(is.na(res["a"])){
    res=c(a=a,res) 
  }
  res[1:2]=c(res["a"]*res["po"],(1-res["a"])*res["po"])
  names(res)[1:4]=c("EMB->HSC","EMB->MPP","HSC->MPP","MPP->HSC")
  res
}

report_results_FULL=function(result.list,b.verbose=FALSE,a=1){
  allres=sapply(result.list,function(x) {
    se=stdEr(x)
    vals=coef(x)
    out=c(invlogit(vals),invlogit(vals[1:2]),invlogit(vals-1.96*se),invlogit(vals+1.96*se),LL=logLik(x))
    names(out)=c(names(vals),paste0(names(vals)[1:2],"_est"),paste0(names(vals),"_lb"),paste0(names(vals),"_ub"),"LL")
    
    out})
  idx=which.max(allres["LL",])
  max.ll=max(allres["LL",])
  if(FALSE){
    ## disabled profile based CI code
    se=stdEr(result.list[[idx]])
    ci=invlogit(confint(result.list[[idx]],std.err=se))
    allres[paste0(rownames(ci),"_lb"),idx]=ci[,1]
    allres[paste0(rownames(ci),"_ub"),idx]=ci[,2]
  }
  
  if(b.verbose){
    print(allres)
    cat("close matches...")
    print(allres[,which(allres["LL",]>max.ll-0.1)])
    print(allres[,idx])
    
  }
  #cat("Num close matches=",length(which(allres["LL",]>max.ll-0.1)),"\n")
  res=allres[,idx]
  if(is.na(res["a"])){
    res=c(a=a,res) 
  }
  res[1:2]=c(res["a"]*res["po"],(1-res["a"])*res["po"])
  names(res)[1:4]=c("EMB->HSC","EMB->MPP","HSC->MPP","MPP->HSC")
  res
}

get_full_results=function(res){
  Q=sapply(res,function(x) suppressWarnings(report_results_FULL(x)))
  QQ=as.data.frame(Q) %>% mutate(TYPE=c(rep("est",6),rep("lb",4),rep("ub",4),"LL"),PARAM=c("EMB->HSC","EMB->MPP","HSC->MPP","MPP->HSC","EMB->HSC/(EMB->HSC+EMB->MPP)","EMB->HSC+EMB->MPP","EMB->HSC/(EMB->HSC+EMB->MPP)","EMB->HSC+EMB->MPP","HSC->MPP","MPP->HSC","EMB->HSC/(EMB->HSC+EMB->MPP)","EMB->HSC+EMB->MPP","HSC->MPP","MPP->HSC","LL"))
  QQQ=QQ %>% tidyr::pivot_longer(cols=names(QQ)[-((ncol(QQ)-1):ncol(QQ))],names_to = "AGE",values_to = "value") %>% filter(TYPE!="LL")
  list(longres=QQQ,wideres=QQ)
}


run.viterbi.discrete=function(tree,mouse,prior.probs,root.prior=c(0,0,1)){
  dff=data.frame(tip.label=tree$tip.label,CELLTYPE=tree$CELLTYPE)
  tree$edge.length=round(tree$edge.length)
  res=viterbi.tree(tree,dff,numStates = 3,
                   emissionFunction = emissionFunction3state,prior.rates = prior.probs ,
                   root.prior = root.prior,type="DISCRETE")
  res$summary$ID=mouse
  if(dim(res$changes)[1]>0){
    res$changes$ID=mouse
  }
  nhh=cbind(tree$edge[,2],nodeHeights(tree)) %>% as.data.frame %>% (function(x) {colnames(x)=c("node","start","end");x})
  nhh$child.count=sapply(nhh$node,function(node) length(get_samples_in_clade(node,tree)))
  res$changes=res$changes %>% inner_join(nhh)
  list(vit=res,pheno=dff)
}



get_last_transition_heights=function(vit){
  heights=sapply(1:length(vit$tree$tip.label),
                 function(i){
                   nodes=get_parents(i,vit$tree$edge)
                   idx=match(nodes,vit$changes$node)
                   as.numeric(vit$changes[idx[which(!is.na(idx))[1]],c("start","end")])
                 }
  )
  data.frame(start=heights[1,],end=heights[2,],celltype=vit$tree$CELLTYPE)
}



### New faster versions
infer_all_trees_FULL_RND=function(trees,mt.start=10,mt.end=50,method="L-BFGS-B",control=list(),epsilon=1e-12){
  START=list(a=runif(1),po=0.2*runif(1),hm=0.01*runif(1),mh=0.01*runif(1))
  LB=list(a=1e-8,po=1e-8,hm=1e-8,mh=1e-8)
  UB=list(a=1-1e-8,po=1-1e-8,hm=0.5,mh=0.5)
  START=lapply(START,logit)
  LB=lapply(LB,logit)
  UB=lapply(UB,logit)
  trees=init_trees(trees,epsilon = epsilon,mt.start=mt.start)
  mm=mle2(LOGLIK_VIT_FULL_CACHED,
          start=START,
          data=list(dat=list(trees=trees,
                             mt.start=mt.start,
                             mt.end=mt.end)
          ),
          method=method,
          lower=LB,
          upper=UB,control=control)
  ###cat("#")
  mm
}
infer_all_trees_FULL_RND_FIX_A=function(trees,a=1,mt.start=10,mt.end=50,method="L-BFGS-B",control=list(),epsilon=1e-12){
  START=list(po=0.2*runif(1),hm=0.01*runif(1),mh=0.01*runif(1))
  LB=list(po=1e-8,hm=1e-8,mh=1e-8)
  UB=list(po=1-1e-8,hm=0.5,mh=0.5)
  START=lapply(START,logit)
  LB=lapply(LB,logit)
  UB=lapply(UB,logit)
  trees=init_trees(trees,epsilon = epsilon,mt.start=mt.start)
  mm=mle2(LOGLIK_VIT_FULL_CACHED_FIX_A,
          start=START,
          data=list(dat=list(trees=trees,
                             a=a,
                             mt.start=mt.start,
                             mt.end=mt.end)
          ),
          method=method,
          lower=LB,
          upper=UB,control=control)
  mm
}
LOGLIK_VIT_FULL_CACHED_FIX_A=function(po,hm,mh,dat){
  LOGLIK_VIT_FULL_CACHED(logit(dat$a),po=po,hm=hm,mh=mh,dat=dat)
}
LOGLIK_VIT_FULL_CACHED=function(a,po,hm,mh,dat){
  trees=dat$trees
  a=invlogit(a)
  po=invlogit(po)
  pp=get_prior_probs_full_2epoch(ph=a*po,
                                 pm=(1-a)*po,
                                 p_m_h=invlogit(mh),
                                 p_h_m=invlogit(hm),
                                 mt.start = dat$mt.start)
  tm=getTransitionMatrices(pp)
  LLL=sapply(names(trees),function(id){
    tree=trees[[id]]
    dat=trees[[id]]$dat
    dat$TM=tm$P[[2]]
    ## Here we use a the STAN based likelihood calculation as it is faster. Have equivalent in pure R if necessary.
    LL=sampling(HMT_STAN_LL_MODEL,data=dat,iter=1,chains=1,cores=1,algorithm="Fixed_param",verbose=FALSE,show_message=FALSE,refresh=0)
    LL=summary(LL)$summary["LL","mean"]
    if(is.infinite(LL)){
      LL=0.1*.Machine$double.xmax*sign(LL)
    }
    LL
  }
  )
  ll=-sum(LLL)
  if(is.infinite(ll)){
    browser()
    #return(-.Machine$double.xmax)
  }
  ll
}
## Initialised parameter independent aspects of the data structures passed to STAN_HMT_LL 
init_trees=function(trees,epsilon=1e-12,mt.start=10){
  for(id in names(trees)){
    tree=trees[[id]]
    tree$edge.length=round(tree$edge.length)
    ##LL=run.upward.discrete.mle(tree,mouse = id,prior.probs = pp)  Pure R implementation.
    dff=data.frame(tip.label=tree$tip.label,CELLTYPE=tree$CELLTYPE)
    tree$edge.length=round(tree$edge.length)
    ## First trim the first mt.start mutations of the tree.. 
    nh=nodeHeights(tree)
    tree$edge.length=ifelse(nh[,1]<=mt.start,
                            ifelse(nh[,2]>=mt.start,tree$edge.length+nh[,1]-mt.start,0),
                            tree$edge.length)
    N=length(tree$edge.length)
    M=length(tree$tip.label)
    tree=add_visit_vector(tree)
    tree$delta=matrix(1,ncol=3,nrow=N+1)
    idx=match(1:M,tree$edge[,2])
    tree$delta[idx,]=t(sapply(1:M,function(i) {
      if(dff$CELLTYPE[i]=="HSC"){ 
        c(1-epsilon,epsilon*0.5,epsilon*0.5)
      }else{
        if(dff$CELLTYPE[i]=="MPP"){
          c(epsilon*0.5,1-epsilon,epsilon*0.5)
        }else{
          c(epsilon*0.5,epsilon*0.5,1-epsilon)
        }
      }
    })
    )
    EPS=1e-308
    ## prepare the data for STAN
    dat=list(N1=N,
             C1=ncol(tree$children),## Max number of children. Non-sparse representation!
             M1=M,
             children1=tree$children,#IDX of children in ee matrix . On each row children after and including first negative.
             visit_order1=tree$visits, ## Proceeding in this way goes up the tree.
             bl1=tree$edge.length, ##Length of each branch.  
             delta1=tree$delta,  ##Delta initialised with tips. probability of observed state conditional on underlying state (3 state)
             root_prior=c(EPS,EPS,1-(2*EPS)), #Start in EMB state (0,0,1)
             TM=NULL
    )
    tree$dat=dat
    trees[[id]]=tree
  }
  trees
}
