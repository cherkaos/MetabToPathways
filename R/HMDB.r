library(httr)
library(XML)
library(RCurl)


# HMDB: Function to get information (ex: pathways) of a compound
# info can be: pathways, synonyms, cellularLocation, biofluidLOCATION, tissueLOCATION, class (compound<s class)
# Input: hmdb id of the compound
# organism: human specific

getHMDB <- function(id,info="pathways"){
                if(info =="pathways")info = "/pathway/name"
				else if(info =="/synonyms")info = "synonym"
				else if(info =="/cellularLOCATION")info = "tissue"
				else if(info =="class")info = "super_class"
				else if(info =="name")info = "name"
		
	if(!id=="error"|| !is.null(id)){
		# Create the url
		url <- paste0("http://www.hmdb.ca/metabolites/",id,".xml")
		xmlhttp<-getURL(url)
		 # To check if HMDB returns a xml
		if(grepl("metabolite",xmlhttp, ignore.case = TRUE)==TRUE) 
		{
			if(is.null(xmlhttp)){return("error")}
			else{
				doc<-(xmlParse(xmlhttp))
				# Parsing the xml file to get info wanted
				src <- tryCatch(xpathSApply(doc, paste0("/",info), xmlValue),error=function(e){NULL})
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
