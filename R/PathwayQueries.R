library(httr)
library(XML)
library(RCurl)

# Function to get kegg database, obtained from Ideom
IDEOM<-function(url="https://gist.github.com/dgrapov/5548790/raw/399f0958306c1018a6be846f58fd076ae83f1b78/IDEOM+small+list"){
	options(warn=-1)
	if(require(RCurl)==FALSE){
		install.packages("RCurl");library(RCurl)
	}
	 else{
		 library(RCurl)
	 }
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

# Common function to fix list
fixlc<-function(obj){
	as.character(unlist(obj))
}

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


# Lookup table for Kegg and Smpdb database, locally saved  
DBKEGG<-IDEOM()
DBSMPDB<-read.csv("http://www.smpdb.ca/public/system/current/metabolite_mapping.csv") 

# KEGG: Function to find a compound information using his id
# Type of id can be : name, pubchem, hmdb, kegg, biocyc, inchikey
# Returns a list: Name, Pubchem, HMDB, KeggID, BioCyc, InchiKey, Map, KeggPathway (as a list), KeggPathID
SearchKEGG<-function(id,type){ 
	position<-vector()
	if(type=="name"){position<-Mmatch(id,DBKEGG[,3])}	
	if(type=="pubchem"){position<-Mmatch(id,DBKEGG[,9])}		
	if(type=="hmdb"){position<-Mmatch(id,DBKEGG[,11])}	
	if(type=="kegg"){position<-Mmatch(id,DB1DBKEGG[,8])}	
	if(type=="biocyc"){position<-Mmatch(id,DBKEGG[,12])}	
	if(type=="inchikey"){position<-Mmatch(id,DBKEGG[,15])}	
	
	if(length(position)>0) # If compound found
	{
		# Store information of that compound
		Name<-fixlc(DBKEGG[position,3])
		Pubchem<-fixlc(DBKEGG[position,9])
		HMDB<-fixlc(DBKEGG[position,11])
		KeggID<-fixlc(DBKEGG[position,8]) 
		BioCyc<-fixlc(DBKEGG[position,12])
		InchiKey<-fixlc(DBKEGG[position,15])
			
		# Store info of pathway of that compound
		Map<-fixlc(DBKEGG[position,5])
		# Separate Pathways name, store it as a list
		KeggPathway<-fixlc(strsplit(DBKEGG[position,6],"__")) 
		KeggPathID<-fixlc(DBSMPDB[position,7])
		return(list("Name"=Name,"Pubchem"=Pubchem,"HMDB"=HMDB,"KeggID"=KeggID,"BioCyc"=BioCyc,"InchiKey"=InchiKey,"Map"=Map,"KeggPathway"=KeggPathway,"KeggPathID"=KeggPathID))
	}
	else(return("error"))
}

# SMPDB: Function to find a compound information using his id
# Type of id can be : name, hmdb, kegg
# Returns a list: Name, HMDB, KeggID, SmpdbPathway
SearchSMPDB<-function(id,type){ 
	position<-vector()
	if(type=="name"){position<-Mmatch(id,DBSMPDB[,4])}			
	if(type=="hmdb"){position<-Mmatch(id,DBSMPDB[,3])}	
	if(type=="kegg"){position<-Mmatch(id,DBSMPDB[,5])}	
	
	if(!length(position)==0) # If compound found
	{
	
		Name<-fixlc(DBSMPDB[position[1],4])
		HMDB<-fixlc(DBSMPDB[position[1],3])
		KeggID<-fixlc(DBSMPDB[position[1],5])
		SmpdbPathway<-fixlc(DBSMPDB[position,2])


		return(list("Name"=Name,"HMDB"=HMDB,"KeggID"=KeggID,"SmpdbPathway"=SmpdbPathway))		
	}
	else(return("error"))
}

# HMDB: Function to get info (information you want to get) of a compound
# info can be: Pathways, Biocyc, Synonyms, CellularLocation, BiofluidLocation, TissueLocation
# Must have hmdb id of the compound
SearchHMDB <- function(hmdb,info="Pathways"){

                switch(info,
                "Pathways" = "pathway/name",
                "Biocyc" = "biocyc_id",
                "Synonyms" = "synonym",
				"CellularLocation"="cellular_location",
				"BiofluidLocation"="biofluid",
				"TissueLocation"="tissue"
				)
		
	if(!hmdb=="error"|| !is.null(hmdb)){
		# Create the url
		url <- paste0("http://www.hmdb.ca/metabolites/",hmdb,".xml")
		xmlhttp<-getURL(url)
		 # To check if HMDB returns a xml
		if(grepl("metabolite",xmlhttp, ignore.case = TRUE)==TRUE) 
		{
			if(is.null(xmlhttp)){return("error")}
			else{
				doc<-(xmlParse(xmlhttp))
				# Parsing the xml file to get info wanted
				src <- tryCatch(xpathSApply(doc, paste0("//",info), xmlValue),error=function(e){NULL})
				if(length(src)>0){
					src<-unlist(src)
					return(src)
				}
				else(return("error"))		
			}
		}
		else{return("error")}
	}
	else{return("error")}	
}


# BIOCYC: Function to get pathway of a compound
# organism can be: Multiorganism (META), Human, Ecoli, Basillus
# Most have hmdb id of the compound
SearchBIOCYC <- function(biocyc,organism="Multiorganism"){
    switch(organism,
    "Multiorganism" = "META",
    "Human" = "HUMAN",
    "Ecoli" = "ECOLI",
	"Basillus"="BSUB",
	)
	
	if(!biocyc=="error"){
		url <- paste0("http://websvc.biocyc.org/apixml?fn=pathways-of-compound&id=",organism,":",biocyc)
		xmlhttp<-getURL(url)
		# To check if Biocyc resturns a xml
		if(grepl("ptools-xml",xmlhttp, ignore.case = TRUE)==TRUE)
		{
			if(is.null(xmlhttp)){return("error")}
			else{
				doc<-xmlParse(xmlhttp)
				# Paser to get common name of Pathway
				src <- tryCatch(xpathApply(doc, "//common-name", xmlValue),error=function(e){NULL})
				html<-c("&","<i>","</i>")
				for(i in 1:length(html)){src<-sub(html[i],"",src,fixed = FALSE)}
				src<-as.vector(fixlc(src))
				return(src)
			}
		}
		else{return("error")}
	}
	else{return("error")}
}


##To get the input
#vars<-data.frame(id = c(1:nrow(test.index)),p.value=test.index[,i], CID = col.data$PubChem.id, name = col.metadata$BinBase.name, kegg=col.metadata$KEGG.id )
#vars<-vars[order(test.index[,i,drop=F],decreasing=TRUE),]

#limit to variables with CIDs (do conversions database conversions seperately)
#tmp<-vars[!is.na(as.numeric(fixlc(vars$CID)))& as.numeric(vars$p.value)<=alpha,]

#input is the table of top compound, with CID and name 
#Put in a csv and send in result






                          

