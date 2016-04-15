library(igraph)
edgeList=read.table("data/fb_caltech_small_edgelist.txt", sep=" ", header=FALSE)
attrList=read.csv("data/fb_caltech_small_attrlist.csv", sep=",", header=TRUE)
edgeList=edgeList+1
graph=graph.edgelist(as.matrix(edgeList), directed=FALSE)
attrnames=names(attrList)
for(j in 1:length(attrnames)){
	for(i in 1:max(edgeList)){
		graph=set_vertex_attr(graph, index=i, name=attrnames[j], value=attrList[i,j])
	}
 }
 graph=set.edge.attribute(graph, name="weight", value=1)
simil = function(m1, m2){
 # Check dimensions are the same
	if (any(dim(m1) != dim(m2)))
		stop(paste("ERROR: matrices are not the same size: ",nrow(m1), "x", ncol(m1), "vs",nrow(m2), "x", ncol(m2)))
	# Linearize the matrices
	m1 = as.vector(m1)
	m2 = as.vector(m2)
	# Cosine similarity
	similarity = (m1%*%m2)/sqrt((m1%*%m1) * (m2%*%m2))
	return(similarity)
 }
 
 get.attr.vector=function(vertex.id){
    attri=0
    for (n in 1:length(attrnames)){
        attri=c(attri, get.vertex.attribute(graph, attrnames[n], vertex.id))
    }
    attri=attri[-1]
    return (attri)
 }
 
alpha=as.numeric(commandArgs(trailingOnly = TRUE))

 times=1
 
 memb.table=data.frame(vertices=c(1:length(V(graph))), id=1:length(V(graph)))
 #set.seed(111)
 
 repeat{
 print(paste("Running : ", times, " /15 iteration"))
 #print(nrow(memb.table))
 
 # Phase 1
 nvert=length(V(graph))
 membership=1:length(V(graph))
 randomized=sample(1:nvert)
 iter=1
 
 attrList=t(sapply(1:length(V(graph)),  function(v) get.attr.vector(v)))
 repeat{
	#print (iter)
	i=randomized[iter]
	max=0
	w=i
	temp.memb=membership
	init.Q.new=modularity(graph, membership)
	
	init.Q.attr=0
	memb.i=which(membership==membership[i])
	if(length(memb.i==1)){
	    init.Q.attr=0
	    } else{
	        for(c in memb.i){
			    init.Q.attr=init.Q.attr+(as.numeric(simil(as.matrix(attrList[i,]), as.matrix(attrList[c,]))))
			}
		}
	init.Q.attr=init.Q.attr/(length(memb.i)*length(memb.i))
	
	for (j in neighbors(graph, i)){
		if(temp.memb[i]==temp.memb[j]){
		next
		}
		else{
		memb.j=which(temp.memb==temp.memb[j])
		flag=1
		temp.memb[i]=temp.memb[j]
		new.Q.new=modularity(graph, temp.memb)
		new.Q.attr=0
		
		for(k in memb.j){
		    
			new.Q.attr=new.Q.attr+(as.numeric(simil(as.matrix(attrList[i,]), as.matrix(attrList[k,]))))
			}
		new.Q.attr=new.Q.attr/(length(memb.j)*length(memb.j))
		diff=(alpha*(new.Q.new-init.Q.new))+((1-alpha)*(new.Q.attr-init.Q.attr))
		if(diff>max){
			max=diff
			w=j
		}
		}
	}
	membership[i]=membership[w]
	
	
    if(diff<=0 || iter==length(V(graph)) ){
	#if(iter==15){
		break
	}
	iter=iter+1
	
	}
	
	#Phase 2
	membership=membership+1000
	temp=1
	temp.memb.table=data.frame(vertices=0, id=0)
	for(i in 1:length(unique(membership))){
	    label=membership[temp]
	    if(label<i){
	        while(label<i){
	            temp=temp+1
	            label=membership[temp]
	            
	        }
	    }
	        
	        lst.to.change=which(membership==(label))
	        list.to.append=subset(memb.table, id==label-1000)$vertices
	        
	        for (j in lst.to.change){
	            membership[j]=i
	            list.to.append=c(list.to.append, subset(memb.table, id==j)$vertices)
	        }
	        temp.memb.table=rbind(temp.memb.table, data.frame(vertices=list.to.append, id=i))
	}
	
	memb.table=temp.memb.table
	
	graph=contract.vertices(graph, membership, vertex.attr.comb="median")
	graph=simplify(graph, edge.attr.comb="sum")
	
	
	if(times==15 || vcount(graph)==1){
	    break;
	}
	times=times+1
}

 for( comm in 1:length(V(graph))){
    text.to.append=paste(unique(subset(memb.table, id==comm)$vertices), collapse=",")
    write.table(text.to.append,"communities.txt", append=TRUE, quote=FALSE, eol="\n", col.names=FALSE, row.names=FALSE)
    
 }
  
