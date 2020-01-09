# Fragment molecule
# USE R 3.6 OR HIGHER!!!!!!!
# USE BOTH R AND JAVA x64!!!!!!

require(rcdk)
require(RSQLite)

# Connect to PubChem_2019_09_09.db

con <- dbConnect(RSQLite::SQLite(), "PubChem_2019_09_09.db")
dbListTables(con)
dbListFields(con, "molecules")

# Select molecule
Target_Molecule_Formula="C8H9NO2"

res <- dbSendQuery(con, paste("SELECT * FROM molecules WHERE formula = '",Target_Molecule_Formula,"'",sep=""))
Isomer_Molecules=dbFetch(res)

head(Isomer_Molecules)

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

Isomer_Molecules_1$Total_L_H=Isomer_Molecules_1$NH+2*Isomer_Molecules_1$NH2+Isomer_Molecules_1$OH+Isomer_Molecules_1$SH
Isomer_Molecules_1$Total_L_O=Isomer_Molecules_1$Carbonyle+Isomer_Molecules_1$Carboxylic_Acid

Isomer_Molecules_H=Isomer_Molecules_1[Isomer_Molecules_1$Total_L_H==2,]
nrow(Isomer_Molecules_H)

Isomer_Molecules_HO=Isomer_Molecules_H[Isomer_Molecules_H$Total_L_O==1,]
nrow(Isomer_Molecules_HO)

write.csv(Isomer_Molecules_1,"Isomer_Molecules.csv")

TARGET_SMILES="CC(O)=Nc1ccc(O)cc1" # Paracetomol

# Run fragmentation for each molecule
source("!CFM Fragmentation v2.R")
dev.off()

# Predict LC retention time for each molecule
source("!Predict retention time v1.R")

# Analyze fragments
source("!Fragments Analysis v1")

#
#s=head(Isomer_Molecules,500)
#JAVA_molecules=sapply(s$smiles, parse.smiles)
#
#for (group in Functional_Groups) {
#	mappings <- matches(group, JAVA_molecules, TRUE)
#
#	Number_Of_Groups = as.numeric(sapply(mappings, function (mappings) {a=mappings[[2]]; return (length(a))}))
#	
#	s=cbind(s,Number_Of_Groups)
#}

plot(Isomer_Molecules_1$Total_L_H,Isomer_Molecules_1$Total_L_O)











