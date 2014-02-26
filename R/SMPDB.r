library(httr)
library(XML)
library(RCurl)

# Function for multiple matching
# Pattern is the character to match, list is where to look
# If multiple occurrence of the pattern, return all the position of the pattern in the list

Mmatch<-function(pattern,list){
	if(!(length(list)>0) || !(length(pattern)>0) ){return(NULL)}else{
		list<-fixlc(list)
		position<-lapply(1:length(list),function(i){
		if(length(grep(pattern,list[i],ignore.case = TRUE))>0) {
			if((nchar(as.vector(pattern))==nchar(as.vector(list[i]))))
			{return(i)}}
		})
		return(unlist(position))
	}	
}

# Common function to fix list
fixlc<-function(obj){
	as.character(unlist(obj))
}


getDBSMPDB<-function(url="https://gist.githubusercontent.com/cherkaos/9190435/raw/ba9859965f666f29b730cbb9ebd67f91b4557a76/SMPDB_metabolites"){
	options(warn=-1)
	DB<-tryCatch(getURL(url,ssl.verifypeer=FALSE),error=function(e){NULL})
	tmp<-strsplit(DB,"\\n")
	tmp2<-strsplit(as.character(unlist(tmp)), "\t")
	#convert to matrix
	obj<-t(do.call("cbind",sapply(tmp2,unlist)))
	#try to fix colnames
	names<-unlist(strsplit(obj[1,]," "))[1:ncol(obj)]
	tmp<-obj[-1,]
	colnames(tmp)<-names
	return(as.matrix(tmp))
}	

DBSMPDB<-getDBSMPDB()


# SMPDB: Function to get pathway of a compound
# Type of id can be : name, hmdb, kegg, chebi, drugbank
# Returns a list: Name, HMDB, KeggID, SmpdbPathway
# organism: human specific

getSMPDB<-function(id,idtype="hmdb"){ 
	position<-vector()
	if(idtype=="name"){position<-Mmatch(id,DBSMPDB[,4])}			
	if(idtype=="hmdb"){position<-Mmatch(id,DBSMPDB[,3])}	
	if(idtype=="kegg"){position<-Mmatch(id,DBSMPDB[,5])}	
	if(idtype=="chebi"){position<-Mmatch(id,DBSMPDB[,6])}	
	if(idtype=="drugbank"){position<-Mmatch(id,DBSMPDB[,7])}	
	
	
	if(!length(position)==0) # If compound found
	{
		smpdbPATHWAY<-DBSMPDB[position,2]
		names(smpdbPATHWAY) <-  DBSMPDB[position,1]
		return(smpdbPATHWAY)		
	}else(return("error"))
}
