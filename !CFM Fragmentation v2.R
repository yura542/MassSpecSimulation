# Fragment molecule
# USE R 3.6 OR HIGHER!!!!!!!
# USE BOTH R AND JAVA x64!!!!!!

require(rcdk)

Isomer_Molecules=read.csv("Isomer_Molecules.csv", header=TRUE, sep=',')

pdf("MS2 Fragmentation.pdf")
pb <- txtProgressBar(min = 0, max = nrow(Isomer_Molecules), style = 3)
k=0

for (index in c(1:nrow(Isomer_Molecules))) {
	Isomer_Molecules_index=Isomer_Molecules[index,]
	SMILES_str=Isomer_Molecules_index$smiles

	k=k+1
	#SMILES_str="CCCC(C(=O)C1=CC2=C(C=C1)OCO2)N3CCCC3"
	
	results=shell(paste("cfm-predict.exe ", SMILES_str, " 0.001 param_output.log param_config.txt 1",sep=""), intern=TRUE)
	if (!(TRUE %in% grepl("Unsupported",results))) {
		results_fragments=results[c((1+which(results=="")):length(results))]
		
		# determine unique fragments
		
		fragments=data.frame()
		for (i in c(1:length(results_fragments))) {
			line_i=results_fragments[i]
			line_i_split=strsplit(line_i,split=" ")[[1]]
			fragments=rbind(fragments,line_i_split)
				i_factor <- sapply(fragments, is.factor)
				fragments[i_factor] <- lapply(fragments[i_factor], as.character)
		}
		names(fragments)=c("n","m","smiles")
		fragments$i=0
		
		# assign intensities
		
		results_fragments_int=results[c(1:(which(results=="")-1))]
		results_fragments_int=results_fragments_int[-grep("ener",results_fragments_int)]
		
		for (i in c(1:length(results_fragments_int))) {
			line_i=results_fragments_int[i]
			line_i_split=strsplit(line_i,split=" ")[[1]]; line_i_split=gsub("[:(:]","",line_i_split); line_i_split=gsub("[:):]","",line_i_split)
			number_of_fragments=(length(line_i_split)-2)/2
			for (ii in c(1:number_of_fragments)) {
				fragment_n=line_i_split[2+ii]
				fragment_i=line_i_split[2+number_of_fragments+ii]
				
				fragments[fragments$n==fragment_n,]$i=fragments[fragments$n==fragment_n,]$i+as.numeric(fragment_i)
			}
			
		}
		
		#plot(fragments$m,fragments$i,ty='h')
		
		# Select most intense fragments
		
		most_intense_number=10
		fragments_most_intense=head(fragments[order(fragments$i,decreasing=TRUE),],most_intense_number)
		
		par(oma = c(1,0,0,0) + 0.1,
			mar = c(0,0.5,0,0)+1,
			mgp = c(1.7,0.5,0))
		layout(matrix(c(1:2), 1, 2, byrow = TRUE),widths=c(3,2),heights=c(1,1))
	
		plot(fragments_most_intense$m,fragments_most_intense$i,ty='h', main=SMILES_str)
		img <- view.image.2d(parse.smiles(as.character(Isomer_Molecules_index$smiles))[[1]])
		rasterImage(img, min(as.numeric(fragments_most_intense$m)),min(as.numeric(fragments_most_intense$i))+(max(as.numeric(fragments_most_intense$i))-min(as.numeric(fragments_most_intense$i)))/2, min(as.numeric(fragments_most_intense$m))+(max(as.numeric(fragments_most_intense$m))-min(as.numeric(fragments_most_intense$m)))/2,max(as.numeric(fragments_most_intense$i)))
		
		
		text(fragments_most_intense$m,fragments_most_intense$i,round(as.numeric(fragments_most_intense$m),3))
		
		plot(c(0:10),c(0:10),col=0,main=paste(as.numeric(Isomer_Molecules_index[,c(9:16)]),collapse=','))
		for (image_n in c(1:min(most_intense_number,nrow(fragments_most_intense)))) {
			img <- view.image.2d(parse.smiles(fragments_most_intense[order(as.numeric(fragments_most_intense$m)),]$smiles[image_n])[[1]])
			if (image_n<=5) {
				rasterImage(img, 0,2*(image_n-1), 5,2*image_n)
			} else {
				rasterImage(img, 5,2*(image_n-6), 10,2*(image_n-5))
			}
		}
		text(c(rep(2.5,5),rep(7.5,5)),c(rep(seq(0,8,2),2)),round(as.numeric(fragments_most_intense[order(as.numeric(fragments_most_intense$m)),]$m),3))
		fragments_most_intense$precursor_smiles=SMILES_str
		write.table(fragments_most_intense,"Fragmented spectra.txt", append=TRUE, col.names=FALSE)
	}
	setTxtProgressBar(pb, k)
}
dev.off()

#results=shell("cfm-predict.exe CCCC(C(=O)C1=CC2=C(C=C1)OCO2)N3CCCC3 0.001 param_output.log param_config.txt 1", intern=TRUE)


















