library(stargazer)

# Function to find, for every compound, his pathways using the 4 databases
# Input list of name, puchem, kegg (or one of the 3)
# 3 lists must be the same size even if they are empty
# Returns a list: Name, Pubchem, HMDB, KeggID, BioCyc, InchiKey, Map, KeggPathway (as a list), KeggPathID, SmpdbPathway (as a list), HmdbPathway (as a list), BiocycPathway (as a list)
FindPathway<-function(name,puchem,kegg){
	sapply(1:length(name), function(i){
		HMDB<-KeggID<-BioCyc<-InchiKey<-Map<-KeggPathway<-KeggPathID<-SmpdbPathway<-HmdbPathway<-BiocycPathway<-""

		# Check if compound is in KEGG
		s1<-SearchKEGG(puchem[i],"pubchem")
		if(!is.list(s1)){s1<-SearchKEGG(kegg[i],"kegg")}
		if(!is.list(s1)){s1<-SearchKEGG(name[i],"name")}
		
		# If it's not in Kegg	
		if(!is.list(s1)){ 
			# Check if compound is in SMPDB
			s2<-SearchSMPDB(kegg[i],"kegg")
			if(!is.list(s2)){s2<-SearchSMPDB(name[i],"name")}		
			if(is.list(s2)){
				HMDB<-s2$HMDB
				keggID<-s2$KeggID
				SmpdbPathway<-s2$SmpdbPathway
			}
		}
		else{ 
			Pubchem<-s1$Pubchem
			HMDB<-s1$HMDB
			KeggID<-s1$KeggID
			BioCyc<-s1$BioCyc
			InchiKey<-s1$InchiKey
			Map<-s1$Map
			KeggPathway<-s1$KeggPathway
			KeggPathID<-s1$KeggPathID	
			
			# Check if compound is in SMPDB
			f1<-SearchSMPDB(HMDB,"hmdb")
			if(!is.list(f1)){f1<-SearchSMPDB(KeggID,"kegg")}
			if(is.list(f1)){SmpdbPathway<-f1$SmpdbPathway}
		}
		
		# Check if compound is in HMDB
		if(nchar(as.vector(HMDB))>5){
			h1<-SearchHMDB(HMDB,"Pathways")
			if(is.vector(h1)){
				if(!h1=="error"){HmdbPathway<-h1}
			}
			
			# Check if compound is in Biocyc
			if(nchar(as.vector(BioCyc))>2){
				q1<-SearchBIOCYC(BioCyc,"Multiorganism")
				if(is.vector(q1)){
					if(!q1=="error"){MetacycPathway<-q1}
				}
			}
			# If no biocyc id, check if hmdb has it
			else{
				biocyc_id<-SearchHMDB(HMDB,"Biocyc")
				if(!biocyc_id=="error"){
					BioCyc<-biocyc_id
					q2<-SearchBIOCYC(BioCyc,"Multiorganism")
					if(is.vector(q2)){BiocycPathway<-q2}
				}
			}
		}	
	list("Name"=as.vector(name[i]),"Pubchem"=as.vector(puchem[i]),"HMDB"=as.vector(HMDB),"KeggID"=as.vector(KeggID),"BioCyc"=as.vector(BioCyc),"InchiKey"=as.vector(InchiKey),"Map"=as.vector(Map),"KeggPathway"=as.vector(KeggPathway),"KeggPathID"=as.vector(KeggPathID),"SmpdbPathway"= as.vector(SmpdbPathway),"HmdbPathway"= as.vector(HmdbPathway),"BiocycPathway"=as.vector(BiocycPathway))
		
	})	
}

# Function to format result from FindPathway to CSV (table)
FindPathwayToCsv<-function(puchem,name,kegg){
	result<-FindPathway(puchem,name,kegg)
	for(i in 1:dim(result)[2]){
		# To store list of pathway in the same string
		result[,i]$KeggPathway<-paste(result[,i]$KeggPathway,collapse="__")
		result[,i]$SmpdbPathway<-paste(result[,i]$SmpdbPathway,collapse="__")
		result[,i]$HmdbPathway<-paste(result[,i]$HmdbPathway,collapse="__")
		result[,i]$BiocycPathway<-paste(result[,i]$BiocycPathway,collapse="__")	
	}
	write.csv(t(as.matrix(result)),"result.csv",row.names=TRUE,col.names=FALSE)
}
	
# Function to find pathway wich has more metabolite using Kegg database
# Input list of name and pubchem. Has to be the same length
# Returns a matrix with: Pathways, Number (of compound in the pathway), and Example (of compound)
TopPathwayKegg<-function(name,pubchem){
	# Matrix where we store evry pathway
	Pathways<-matrix(ncol=3)
	for(i in 1:length(puchem)){
	
		# Look in Kegg
		s1<-SearchKEGG(puchem[i],"pubchem")
		if(!is.list(s1)){s1<-SearchKEGG(name[i],"name")}
		if(is.list(s1)){
			KeggPathway<-s1$KeggPathway
			# If found, look for every pathway of that compound
			for(j in 1:length(KeggPathway))
			{
				if(is.character(KeggPathway[j]) && length(KeggPathway[j])>0 ){
					# Look if pathway already exists in the matrix
					s2<-Mmatch(KeggPathway[j],Pathways[,1])
	
					if(length(s2)>0){
						position<-s2[1]
						# Alreay in list, add 1 to number
						Pathways[position,2]<-as.numeric(Pathways[position,2])+1
						}else{
							# Store pathway in matrix
						Pathways<-rbind(Pathways,c(KeggPathway[j],"1",s1$Name))
					}
				}	
			}
		}
	
	}
	Pathways<-Pathways[!is.na(Pathways[,1]),] 
	Pathways<-Pathways[order(Pathways[,2],decreasing=TRUE),]
	colnames(Pathways)<-c("Pathway","Number","Example")
	return(Pathways)
}	

# Function that output in latex the table from TopPathwayKegg
TopPathwayKeggToLatex<-function(puchem,name,top=20,index=1){
	Pathways<-TopPathKegg(puchem,name)
	Pathways<-Pathways[1:top,]
	
	output <- capture.output(stargazer(as.data.frame(Pathways),title="Pathway Analysis",font.size="scriptsize",type="latex",summary=FALSE, label=paste0("table:tab",index+1,seq=""),table.placement="h"))

	 # cat out the results 
	cat(paste(output, collapse = "\n"), "\n", file="pathways.tex",append = TRUE)
	
}	

	
	