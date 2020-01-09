# Fragment molecule
# USE R 3.6 OR HIGHER!!!!!!!
# USE BOTH R AND JAVA x64!!!!!!

require(rcdk)

s=read.table("Fragmented spectra.txt")
names(s)=c('i','n','m','smiles','int','precursor')

s_precursor=s[s$precursor==TARGET_SMILES,]

Unique_fragments_m=unique(s$m)

unified_fragments=data.frame()

for(precursor in unique(s$precursor)) {
	s1=s[s$precursor==precursor,]
	unified_fragments_vec=rep(0,length(Unique_fragments_m))
	unified_fragments_vec[match(s1$m, Unique_fragments_m)]=1
	unified_fragments=rbind(unified_fragments,unified_fragments_vec)
}

Difference=apply(head(unified_fragments,5000),1,function(x) {return(as.numeric(x-unified_fragments[which (TARGET_SMILES==as.character(unique(s$precursor))),]))})
Difference_ind=apply(Difference,2,function(x) {return(10-length(x[x==-1]))})

hist(Difference_ind,col=2,31)

unique(s$precursor)[Difference_ind==10]

hist(Difference_ind,31,col=2,plot=FALSE)



###### Unresolved
Isomer_Molecules_Unresolved=read.csv("Isomer_Molecules_Unresolved.csv", header=TRUE, sep=',')

s=s[s$precursor %in% Isomer_Molecules_Unresolved$smiles,]

s=s[s$precursor %in% Isomer_Molecules_HO$smiles,]

s=rbind(s_precursor,s)