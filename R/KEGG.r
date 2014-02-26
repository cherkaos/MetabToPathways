library(RCurl)
library(KEGGREST)



# KEGG: Function to get pathways of a compound
# organism can be: multiorganism (general map), human(hsa), ecoli(eco), bacillus(bsu), mouse (mmu), rat(rno)
# Input: kegg id of the compound


getKEGG<-function(id,organism="multiorganism"){ 
	
	compoundINFO<-keggGet(id)
	compoundPATHWAYS<-compoundINFO[[1]]$PATHWAY
	
	if(organism =="multiorganism")organism= "map"
	else if(organism =="human")organism = "hsa"
	else if(organism =="ecoli")organism = "eco"
	else if(organism =="bacillus")organism = "bsu"
	else if(organism =="mouse")organism = "mmu"
	else if(organism =="rat")organism = "rno"
	

	
	pathwaysID <-names(compoundPATHWAYS)
	pathwaysID<-sub("map", organism, pathwaysID )

	# To check if pathway exist for that organism
	pathwaysORGANISM <- sapply(1:length(pathwaysID), function(i){
		pathwaysINFO <- tryCatch(keggGet(pathwaysID[i]),error=function(e){NULL})
		if(length(pathwaysINFO)>0){
				return(compoundPATHWAYS[i])
		}
	})	
	
	pathwaysORGANISM[sapply(pathwaysORGANISM, is.null)] <- NULL
	pathwaysORGANISM<-unlist(pathwaysORGANISM)
	return(pathwaysORGANISM)
}