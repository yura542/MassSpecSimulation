# Fragment molecule
# USE R 3.6 OR HIGHER!!!!!!!
# USE BOTH R AND JAVA x64!!!!!!

require(rcdk)
require(RSQLite)
require(igraph)

# Connect to PubChem_2019_09_09.db

con <- dbConnect(RSQLite::SQLite(), "PubChem_2019_09_09.db")
dbListTables(con)
dbListFields(con, "molecules")

# Select molecule
Target_Molecule_Formula="C8H9NO2"
Target_Molecule_Formula="C16H21NO3"
Target_Molecule_Formula="C15H21NO"
Target_Molecule_Formula="C17H19NO3" 
Target_Molecule_Formula="C8H17NO2"
Target_Molecule_Formula="C12H12N2O3"
Target_Molecule_Formula="C24H31NO4"
Target_Molecule_Formula="C41H64O14"
Target_Molecule_Formula="C5H10N2O2S"
Target_Molecule_Formula="C12H14ClNO2"
Target_Molecule_Formula="C18H24N2O4"
Target_Molecule_Formula="C14H20ClNO2"
Target_Molecule_Formula="C50H73N15O11"
Target_Molecule_Formula="C63H98N18O13S"

Target_Molecule_Formula="C9H11NO4" # DOPA

res <- dbSendQuery(con, paste("SELECT * FROM molecules WHERE formula = '",Target_Molecule_Formula,"'",sep=""))
Isomer_Molecules=dbFetch(res)

head(Isomer_Molecules)
nrow(Isomer_Molecules)

# Determine number of functional grouops
NH =  "[#7;H1]" # 
NH2=  "[#7;H2]"
OH = "[#8;H1]"
SH = "[#16;H1]"

Carbonyle=  "[CX3]=[OX1]"
Carboxylic_Acid = "[CX3](=O)[OX2H1]"

Functional_Groups=c(NH, NH2,OH,SH, Carbonyle,Carboxylic_Acid)
Functional_Groups_Names=c("NH","NH2", "OH","SH", "Carbonyle","Carboxylic_Acid")

Isomer_Molecules_1=data.frame()

pb <- txtProgressBar(min = 0, max = nrow(Isomer_Molecules), style = 3)

for (index in c(1:nrow(Isomer_Molecules))) {
	Isomer_Molecules_index=Isomer_Molecules[index,]
	Isomer_Molecules_JAVA_molecules=sapply(Isomer_Molecules_index$smiles, parse.smiles)
	
	Functional_Groups_vector=c()
	
	if (is.null(Isomer_Molecules_JAVA_molecules[[1]])==FALSE) {
	
		for (group in Functional_Groups) {
			mappings <- matches(group, Isomer_Molecules_JAVA_molecules, TRUE)
		
			Number_Of_Groups = as.numeric(sapply(mappings, function (mappings) {a=mappings[[2]]; return (length(a))}))
			
			Functional_Groups_vector=c(Functional_Groups_vector,Number_Of_Groups)
			
		}
		Functional_Groups_vector_d=as.data.frame(t(Functional_Groups_vector))
		names(Functional_Groups_vector_d)=Functional_Groups_Names
	Functional_Groups_vector_1=cbind(Isomer_Molecules_index,Functional_Groups_vector_d)	
	Isomer_Molecules_1=rbind(Isomer_Molecules_1,Functional_Groups_vector_1)
	
	}
		
	setTxtProgressBar(pb, index)

}

nrow(Isomer_Molecules_1)


Isomer_Molecules_1$Total_L_H=Isomer_Molecules_1$NH+2*Isomer_Molecules_1$NH2+Isomer_Molecules_1$OH+Isomer_Molecules_1$SH
Isomer_Molecules_1$Total_L_O=Isomer_Molecules_1$Carbonyle+Isomer_Molecules_1$Carboxylic_Acid

Isomer_Molecules_H=Isomer_Molecules_1[Isomer_Molecules_1$Total_L_H==5,]
nrow(Isomer_Molecules_H)

Isomer_Molecules_HO=Isomer_Molecules_H[Isomer_Molecules_H$Total_L_O==2,]
nrow(Isomer_Molecules_HO)

########## Plot isomers

Number_Of_Slides=floor(nrow(Isomer_Molecules_HO)/10)+1
pdf("Isomers1.pdf")
pb <- txtProgressBar(min = 0, max = Number_Of_Slides, style = 3)


for (index in c(1:Number_Of_Slides)) {

	Isomer_Molecules_index=Isomer_Molecules_HO[c(((index-1)*10+1):(index*10)),]
	plot(c(0:10),c(0:10),col=0)
	for (image_n in c(1:nrow(Isomer_Molecules_index))) {
		img <- view.image.2d(parse.smiles(Isomer_Molecules_index$smiles[image_n])[[1]])
		if (image_n<=5) {
			rasterImage(img, 0,2*(image_n-1), 5,2*image_n)
		} else {
			rasterImage(img, 5,2*(image_n-6), 10,2*(image_n-5))
		}
	}
	text(c(rep(2.5,5),rep(7.5,5)),c(rep(seq(0,8,2),2)),c(((index-1)*10+1):(index*10)))
	setTxtProgressBar(pb, index)
}

dev.off()

########## Isomers for DOPA


#############################################

SMILES_pattern_ortho="[$(caO);H1]"       # carbon ortho to O with 1 H
SMILES_pattern_para="[$(caaaO);H1]"       # carbon para to O with 1 H

Test_list=c(SMILES_pattern_ortho,SMILES_pattern_para)
Isomer_Molecules_HO$NonLab=0

for (index in c(1:nrow(Isomer_Molecules_HO))) {
	Isomer_Molecules_index=Isomer_Molecules_HO[index,]
	Molecule_JAVA_molecules=sapply(Isomer_Molecules_index$smiles, parse.smiles)
	Atoms_list=c()
	for (pattern in Test_list) {
		mappings  <- matches(pattern, Molecule_JAVA_molecules, TRUE)
		Atoms_list_1=unlist(mappings[[1]][[2]])
		Atoms_list=c(Atoms_list,Atoms_list_1)
		
	}
	Atoms_list_unique=unique(Atoms_list)
	Isomer_Molecules_HO[index,]$NonLab=length(Atoms_list_unique)
}
Isomer_Molecules_HO[Isomer_Molecules_HO$NonLab==3,]

Molecule_smiles=Isomer_Molecules_HO[index,]$smiles  #DOPA
Molecule_JAVA_molecules=sapply(Molecule_smiles, parse.smiles)
get.atom.count(Molecule_JAVA_molecules[[1]])
get.atoms(Molecule_JAVA_molecules[[1]])
generate.2d.coordinates(Molecule_JAVA_molecules[[1]])
adjm1=get.connection.matrix(Molecule_JAVA_molecules[[1]])
g1<-graph.adjacency(adjm1);

mappings <- matches(SMILES_pattern_ortho, Molecule_JAVA_molecules, TRUE)   # PLUS 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mappings
Number_Of_Groups = as.numeric(sapply(mappings, function (mappings) {a=mappings[[2]]; return (length(a))}))
Number_Of_Groups

par(mfrow=c(1,2))
plot(c(0:10),c(0:10),col=0)
img <- view.image.2d(Molecule_JAVA_molecules[[1]])
rasterImage(img, 0,0,10,10)
plot(g1)




