# Fragment molecule
# Launch HILIC trained.RData

require ("rcdk")
require("caret")
require("xgboost")


#' Calculate Chemical Descriptors using RCDK
#' @param x A dataframe with 3 mandatory column: "Name", "InchIKey", "SMILES".
#' @export getCD
#' @return A dataframe with calculated CD. Some compund can be missing if smiles is incorect or if molecule returns error in CDK
#' @examples
#' \donttest{
#' # RP and HILIC are previusly loaded dataset from excel with
#' # Name, InchIKey, SMILES and Retention Time
#' descs <- getCD(RP)
#' descs <- getCD(HILIC)}


getCD <- function(x){

  print(paste0("Converting SMILES..."))

  for (i in 1:nrow(x)) {
    smi <- rcdk::parse.smiles(as.character(unlist(x[i,"SMILES"]))) [[1]]
    smi1 <- rcdk::generate.2d.coordinates(smi)
    smi1 <- rcdk::get.smiles(smi,smiles.flavors(c('CxSmiles')))
    x$SMILES[i] <- smi1
    print(paste0(i," of ",nrow(x)))
  }


  # select all possible descriptors
  descNames <- rcdk::get.desc.names(type = "all")
  # select only one descriptors. This helps to remove compounds that makes errors
  descNames1 <- c('org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor')

  print(paste0("Checking for compound errors..."))

  # calculate only 1 descriptor for all the molecules
  mols_x <- rcdk::parse.smiles(as.character(unlist(x[1,"SMILES"])))
  descs1_x <- rcdk::eval.desc(mols_x, descNames1)

  for (i in 2:nrow(x)) {
    mols1 <- rcdk::parse.smiles(as.character(unlist(x[i,"SMILES"])))
    descs1_x[i,] <- rcdk::eval.desc(mols1, descNames1)
    print(paste0(i," of ",nrow(x)))
  }

  # remove molecules that have NA values with only one descriptor
  x_na <- data.frame(descs1_x,x)

  x_na_rem <- x_na[stats::complete.cases(x_na), ]

  x_na_rem <- x_na_rem [,-c(1:6)]

  # computing the whole descriptos on the good on the clean dataset
  print(paste0("Computing Chemical Descriptors 1 of ",nrow(x_na_rem)," ... Please wait"))

  mols_x1 <- rcdk::parse.smiles(as.character(unlist(x_na_rem[1,"SMILES"])))[[1]]
  rcdk::convert.implicit.to.explicit(mols_x1)
  descs_x_loop <- rcdk::eval.desc(mols_x1, descNames)


  for (i in 2:nrow(x_na_rem)) {
    mols <- rcdk::parse.smiles(as.character(unlist(x_na_rem[i,"SMILES"])))[[1]]
    rcdk::convert.implicit.to.explicit(mols)
    descs_x_loop[i,] <- rcdk::eval.desc(mols, descNames)
    print(paste0(i," of ",nrow(x_na_rem)))
  }
  datadesc <- data.frame(x_na_rem,descs_x_loop)
  return(datadesc)
}



fit.xgboost <- function(x){

      # set up train control for 10 times cross validation
      cv.ctrl <-caret::trainControl(method = "cv",number = 10)

      # These are the tune grid parameters
      xgb.grid <- base::expand.grid(nrounds=c(100,200,300,400,500,600,700),
                              max_depth = c(5),
                              eta = c(0.025,0.05),
                              gamma = c(0.01),
                              colsample_bytree = c(0.75),
                              subsample = c(0.50),
                              min_child_weight = c(0))

      print("Computing model Xgboost  ... Please wait ...")

      # Model training using the above parameters
      set.seed(101)
      model_xgb <-caret::train(RT ~.,
                        data=x,
                        method="xgbTree",
                        metric = "RMSE",
                        trControl=cv.ctrl,
                        tuneGrid=xgb.grid,
                        tuneLength = 14)




      print("End training")


      return(model_xgb)
}



cesc <- function(db_rt){
# calculate center and scale with caret function
db_rt_rt <- data.frame(db_rt$RT)
db_rt2 <- db_rt[,-1]
preProc <- caret::preProcess(db_rt2,method = c("center","scale"),rangeBounds = c(0,1))


return(preProc)

}



proc.data  <- function(x){
            # remove na columns
            db <- x[, !apply(x, 2, function(x) any(is.na(x)) )]

            set.seed(101)

            # remove columns with near zero variance values
            nzvColumns <- caret::nearZeroVar(db)

            mdrDescFH <- db[, -nzvColumns]

            db_rt <- data.frame(mdrDescFH[,-c(1:3)])

            return(db_rt)

          }



getRT.smile <- function(smile="O=C(O)C(NC(=O)C(N)C)CC(C)C",training,model=model,cesc=cesc) {
# without center and scale function activated
  if (missing(cesc)){
# convert smile into CDK object
mol <- parse.smiles(smile)
# set the descriptors that have to be computed
descNames <- get.desc.names(type = "all")
# calculate the descriptors with CDK
target <- eval.desc(mol, descNames)

# convert to matrix and use only column present in custom training set
target_desc <- as.matrix(target[names(target) %in% names(training)])
# make the RTP prediction and print in console
ncolt1 <- ncol(target_desc)
pred_target <- round(stats::predict(model, target_desc),2)

}
# The same but with center and scale function activated
  else {

  mol <- parse.smiles(smile)
  descNames <- get.desc.names(type = "all")
  target <- eval.desc(mol, descNames)

  target <- stats::predict(cesc,target)

  target_desc <- as.matrix(target[names(target) %in% names(training)])

  ncolt1 <- ncol(target_desc)
  pred_target <- round(stats::predict(model, target_desc),2)


  }

return(pred_target)

}

# !!! Uncomment if running NOT HILIC.Rdataframe
# !!! This data frame contains previously trained xgb model

#descs2 <- getCD(head(HILIC,1000))
#db_rt <- proc.data(descs2)
#
#preProc <- cesc(db_rt) #Build a model to use for center and scale a dataframe 
#db_rt <- predict(preProc,db_rt) # use the above created model for center and scale dataframe
#
#inTraining <- caret::createDataPartition(db_rt$XLogP, p = .8, list = FALSE)
#training <- db_rt[ inTraining,]
#
#xgb <- fit.xgboost(training)

Isomer_Molecules=read.csv("Isomer_Molecules.csv", header=TRUE, sep=',')
Isomers=Isomer_Molecules$smiles
Isomer_Molecules$RT=0

pb <- txtProgressBar(min = 0, max = length(Isomers), style = 3)

for (i in c(1:nrow(Isomer_Molecules))) {
	
	RT=try(getRT.smile(smile=as.character(Isomers[i]),training,model=xgb),silent=TRUE)
	if (!grepl("Error",RT)) {
		Isomer_Molecules[i,]$RT=RT
	}
	setTxtProgressBar(pb, i)
}


Isomer_Molecules$RTdif=Isomer_Molecules$RT-getRT.smile(smile=TARGET_SMILES,training,model=xgb)

Isomer_Molecules[as.character(Isomer_Molecules$smiles) == TARGET_SMILES,]

hist(Isomer_Molecules$RTdif,100)

write.csv(Isomer_Molecules[abs(Isomer_Molecules$RTdif)<5/60,],"Isomer_Molecules_Unresolved.csv")
