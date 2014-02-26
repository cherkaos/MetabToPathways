library(httr)
library(XML)
library(RCurl)


# BIOCYC: Function to get pathways of a compound
# organism can be: multiorganism (metacyc), human(humancyc), ecoli(ecocyc), bacillus(bsubcyc)
# Input: biocyc id of the compound

getBIOCYC <- function(biocyc,organism="multiorganism"){

	if(organism =="multiorganism")organism= "META"
	else if(organism =="human")organism = "HUMAN"
	else if(organism =="ecoli")organism = "ECOLI"
	else if(organism =="bacillus")organism = "BSUB"
	
	
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
