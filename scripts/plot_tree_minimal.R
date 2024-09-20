require("RColorBrewer")
require("scales")
require("tidyverse")

#' @title 
#' plots a horizontal aspect tree 
#' 
#' This is a bespoke function for plotting phylo trees. Many alternatives are available!   
#' @description 
#' Plots a phylo tree
#' @param tree APE phylo class.   
#' @param direction Specifies the direction ("down" or "across"). Currently on "down" is supported. 
#' @param cex.label Specifies the tip label size.
#' @param offset Specifies the offset of text labels relative to annotated feature.
#' 
#' @return phylo tree augmented with coords.
#' 
plot_tree=function(tree,
                   direction="down",
                   cex.label=1,
                   offset=0,
                   b_do_not_plot=FALSE,
                   lwd=1,
                   bars=NULL,
                   default_edge_color="darkgrey",
                   ymax=NULL,
                   cex.terminal.dots=0,
                   mar=NULL,
                   b.add.scale=TRUE,
                   left.margin.prop=0.1,
                   vspace.reserve=0,
                   bar.label="",
                   ytick.gap=NA){
  if(is.null(mar)){
    par(mar=c(1, 1, 1, 3) + 0.1)
  }else{
    par(mar=mar)
  }
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ytop=max(coords$a1)
      ymax=max(coords$a1)*1.05
    }else{  #add 24jun
      ytop=ymax
      ymax=1.05*ymax
    }  #end add
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  
  if(is.null(bars)){
    YMIN=0-ymax*0.05-vspace.reserve*ymax
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*left.margin.prop),xmax),ylim=c(YMIN,ymax),xlab="",ylab="")
    #plot(NULL,axes=FALSE,xlim=c(0-(xmax*left.margin.prop),xmax),ylim=c(0-(ymax*0.05),ymax),xlab="",ylab="")
  }else{
    ##ymax not supported..
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*left.margin.prop),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=ytop
    #top=max(coords$a1)
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
      text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=tipcol,cex=cex.terminal.dots,pch=19)
    }
  }
  tree$direction=direction
  tree$top=top
  scales=c(0,1,10,50,100,200,500,1000,2000,5000)
  #scale=scales[max(which(ymax/2>scales))]
  scale=scales[max(which(ymax/5>=scales))]
  if(scale==0){
    scale=10**floor(log10(ymax))
  }
  if(!is.na(ytick.gap)){
    scale=ytick.gap
  }
  #browser()
  #cat("scale=",scale,"\n")
  if(b.add.scale){
    #browser()
    axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2)
  }
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("%s Muts",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,col = "grey")
    axis(side = 2,at=c(-ymax*0.15,ytop=-ymax*0.15+ymax*0.1),label=trimws(sprintf("%3.2g",c(0,maxbar))),las=2,line=-2)
    mtext(side=2,at=-ymax*0.15+ymax*0.05,text = bar.label,line = 2)
  }
  tree$ymax=ymax
  tree$vspace.reserve=vspace.reserve
  tree
}

set_cedge=function(parent,tree){
  for(child in get_node_children(parent,tree)){
    #cat("\nsetting parent - child ",parent,child,"\n")
    child.idx=which(tree$edge[,2]==child)
    parent.idx=which(tree$edge[,2]==parent)
    if(length(parent.idx)==0){
      pedge=0
    }else{
      pedge=tree$cedge[parent.idx]
    }
    #cat("before:",length(which(tree$cedge>0)),"\n")
    tree$cedge[child.idx]=pedge+tree$edge.length[child.idx];
    #cat("after:",length(which(tree$cedge>0)),"\n")
    #cat(parent,child,":",tree$cedge)
    tree=set_cedge(child,tree)
  }
  tree
}

##Not very efficient -- recursively calculates height as average of children's height.
get_height=function(tree,node){
  #if(is.null(tree$height)){
  #  tree$height=rep(NA,1:length(tree$edge.length))
  #}
  N=length(tree$tip.label)
  if(node<=N){
    return(node)#tree$height[which(tree$edge[,2]==node)]=node
  }else{
    children=get_node_children(node,tree)
    return(mean(sapply(children,function(child) get_height(tree,child))))
  }
}

set_height=function(tree){
  tree$height_end=sapply(tree$edge[,2],function(i) get_height(tree,i))
  tree$height_start=tree$height[match(tree$edge[,1],tree$edge[,2])]
  N=length(tree$tip.label)
  root=N+1
  idx=which(tree$edge[,1]==root)
  tree$height_start[idx]=mean(tree$height_end[idx])
  tree
}

##horizontal elbow
elbow=function(x0,x1,y0,y1,...){
  arrows(x0=x0,y0=y0,y1=y1,length=0,...)
  arrows(x0=x0,x1=x1,y0=y1,length=0,...)
}
##vertical elbow
elbowv=function(x0,x1,y0,y1,...){
  #browser()
  arrows(x0=x0,x1=x1,y0=y0,length=0,...)
  arrows(x0=x1,y0=y0,y1=y1,length=0,...)

}

set_tree_coords=function(atree){
  ##get the cumulative distance from the root.
  tree=atree
  tree$cedge=rep(0,length(tree$edge.length))
  N=length(tree$tip.label)
  root=N+1
  tree=set_cedge(root,tree)
  tt=set_height(tree)
  atree$coords=data.frame(a0=tt$cedge-tt$edge.length,a1=tt$cedge,
                          b0=tt$height_start,b1=tt$height_end,stringsAsFactors = FALSE)
  if(!is.null(atree$color)){
    atree$coords$color=atree$color
  }
  atree
}

get_all_node_children=function(node,tree){
  children=tree$edge[which(tree$edge[,1]==node),2]
  offspring=children
  for(child in children){
    offspring=c(offspring,get_all_node_children(child,tree))
  }
  offspring
}
get_node_children=function(node,tree){
  tree$edge[which(tree$edge[,1]==node),2]
}

get_samples_in_clade=function(node,tree){
  if(node<=length(tree$tip.label)){
    return(tree$tip.label[node])
  }
  tree$tip.label[intersect(get_all_node_children(node,tree),1:length(tree$tip.label))]
}

get_y_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$top-tree$coords[idx,c("a1","a0")])
}

get_x_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$coords[idx,c("b1","b0")])
}

get_parents=function(node,mut_table,exclude_root=TRUE){
  idx=which(mut_table[,2]==node)
  parents=node##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=mut_table[idx,1]
    parents=c(parents,parent)
    idx=which(mut_table[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)]
  }else{
    parents
  }
}
##The following ties the PDX datastructure to the

get_idx_for_node=function(pdx,node){
  which(pdx$details$node==node)
}
###


get_edge_info=function(pdx,tree,node){
  y=get_y_range(tree,node)
  x=get_x_range(tree,node)
  idx=get_idx_for_node(pdx,node)
  samples=get_samples_in_clade(node,tree)
  list(yb=y[1],yt=y[2],x=x[1],xm=x[2],idx.in.details=idx,samples=samples)
}

add_annotation=function(pdx,tree,annot_function,control=NULL){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(pdx,tree,tree$edge[i,2],control=control))
}



add_defaults=function(control,defaults,mandatory_fields=c()){
  if(is.null(control)){
    control=list()
  }
  for(field in mandatory_fields){
    if(is.null(control[[field]])){
      stop(sprintf("Required parameter %s not supplied in control list",field))
    }
  }

  for(field in names(defaults)){
    if(is.null(control[[field]])){
      control[[field]]=defaults[[field]]
    }
  }
  control
}


node_labels=function(tree,col="blue",cex=1,b_include_tips=FALSE){
  idx=which(tree$edge[,2]>length(tree$tip.label))
  if(b_include_tips){
    idx=1:dim(tree$edge)[1]
  }else{
    idx=which(tree$edge[,2]>length(tree$tip.label))
  }
  text(tree$edge[idx,2],x=tree$coords$b1[idx],y=tree$top-tree$coords$a1[idx],col=col,cex=cex)
}

add_heatmap=function(tree,heatmap,heatvals=NULL,border="black",cex.label=2,pos=4,force.count=-1,txtcols=rep("black",dim(heatmap)[1]),b.add.lines=FALSE){

  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.05*ymax
  if(force.count<0){
    force.count=dim(heatmap)[1]
  }
  gap=ymax*tree$vspace.reserve/force.count
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-gap
    if(!is.na(border) && border=="none"){
      ##collapse 
      cc=unique(heatmap[i,])
      for(ccc in cc){
        idxx=which(heatmap[i,]==ccc)
        start=idxx[1]
        end=start
        for(j in idxx){
          if(j>(end+1)){
            rect(xleft=start-0.5,xright=end+0.5,ybottom = bot,ytop=top,col = heatmap[i,end],border=NA)
            start=j
            end=j
          }else{
            end=j
          }
        }
        rect(xleft=start-0.5,xright=end+0.5,ybottom = bot,ytop=top,col = heatmap[i,end],border=NA)
      }
    }else{
      rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd=0,ljoin=2)
      if(TRUE){
        offsets=0.05
      ##draw lines for boundaries...
      ch=heatmap[i,1]
      k=1
      LWD=1
      LTY="dotted"
      if(b.add.lines && ch!="white"){
        segments(x0 = idx[k]-0.5,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
      }
      for(k in 2:dim(heatmap)[2]){
        
        if(heatmap[i,k]!=ch){
          
          if(b.add.lines && ch!="white"){
            segments(x0 = idx[k]-0.5-offsets,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
          }
          ch=heatmap[i,k]
          if (b.add.lines && ch!="white"){
            segments(x0 = idx[k]-0.5+offsets,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
          }
        }
      }
      if(b.add.lines && ch!="white"){
        segments(x0 = idx[k]+0.5,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
      }
      }
    }
    
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      text(labels[i],x=par("usr")[1]+0.1,y=0.5*(top+bot),pos = pos,cex = cex.label,col=txtcols[i])
    }
    top=bot
  }
  tree
}

add_mini_axis=function(x,ybot,ytop,breaks=c(0,0.5,1),labels=breaks,cex=0.8){
  segments(x0=x,y0=ybot,y1=ytop,lwd=1)
  segments(x0=x-0.2,x1=x,y0=ybot+breaks*(ytop-ybot),lwd=1)
  text(x=x-0.2,y=ybot+breaks*(ytop-ybot),labels = labels,pos = 2,cex=cex)
}


add_unit_heat=function(tree,heatmap,heatvals,border="black",cex.label=2,pos=4,force.count=-1,txtcols=rep("black",dim(heatmap)[1]),
                       b.add.lines=FALSE,
                       b.add.halfway=FALSE,ranges=rep(20,dim(heatmap)[1])){
  
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.05*ymax
  if(force.count<0){
    force.count=dim(heatmap)[1]
  }
  gap=0.9*ymax*tree$vspace.reserve/force.count
  gap2=0.1*ymax*tree$vspace.reserve/force.count
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-gap
    if(!is.na(border) && border=="none"){
      ##collapse 
      cc=unique(heatmap[i,])
      for(ccc in cc){
        idxx=which(heatmap[i,]==ccc)
        start=idxx[1]
        end=start
        for(j in idxx){
          if(j>(end+1)){
            rect(xleft=start-0.5,xright=end+0.5,ybottom = bot,ytop=bot+(top-bot)*heatvals[i,end],col = heatmap[i,end],border=NA)
            
            start=j
            end=j
          }else{
            end=j
          }
        }
        rect(xleft=start-0.5,xright=end+0.5,ybottom = bot,ytop=bot+(top-bot)*heatvals[i,end],col = heatmap[i,end],border=NA)
        add_mini_axis(0,ybot = bot,ytop=top,breaks = c(0,10,20))
      }
    }else{
      if(b.add.halfway && !is.na(ranges[i])){
        segments(x0=min(idx)-0.5,x1=max(idx)+0.5,y0=bot+0.5*(top-bot),lwd=0.5,lty="dotted")
      }
      #browser()
      rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=bot+(top-bot)*heatvals[i,],col = heatmap[i,],border=border,lwd=0,ljoin=2)
      if(!is.na(ranges[i])) {
        add_mini_axis(-1,ybot = bot,ytop = top,labels=c("",ranges[i]/2,ranges[i]))
        segments(x0=min(idx)-0.5,x1=max(idx)+0.5,y0=bot,lwd=0.5)
      }
      
      if(TRUE){
        offsets=0.05
        ##draw lines for boundaries...
        ch=heatmap[i,1]
        k=1
        LWD=1
        LTY="dotted"
        if(b.add.lines && ch!="white"){
          segments(x0 = idx[k]-0.5,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
        }
        for(k in 2:dim(heatmap)[2]){
          
          if(heatmap[i,k]!=ch){
            
            if(b.add.lines && ch!="white"){
              segments(x0 = idx[k]-0.5-offsets,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
            }
            ch=heatmap[i,k]
            if (b.add.lines && ch!="white"){
              segments(x0 = idx[k]-0.5+offsets,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
            }
          }
        }
        if(b.add.lines && ch!="white"){
          segments(x0 = idx[k]+0.5,y0=bot,y1=tree$top,col=ch,lwd=LWD,lty=LTY)
        }
      }
    }
    
    
    if(!is.null(labels)){
      text(labels[i],x=par("usr")[1]+0.1,y=0.5*(top+bot),pos = pos,cex = cex.label,col=txtcols[i])
    }
    top=bot-gap2
  }
  tree
}

cut_tree=function(tree,nmuts){
  nh=nodeHeights(tree)
  idx=which(nh[,2]>nmuts & nh[,1]<nmuts)
  if(length(idx)>0){
    tree$edge.length[idx]=nmuts-nh[idx,1]
  }
  idx=which(nh[,1]>=nmuts)
  if(length(idx)>0){
    tree$edge.length[idx]=0
  }
  tree
}
