#Comparing spatial conservation prioritisation methods using site- versus linkage-based connectivity
#This is the underlying code for the research paper found at ____. Section 3 is for generating simulated seascapes and running Marxan spatial prioritisation on these. Prioritisation is run without using connectivity (baseline), using connectivity-based features methods (where metrics of connectivity are calculated for each planning unit and targets are set for these), or using spatial dependency methods (described in Beger et al. 2010 in Conservation Letters). Section 4 is for running Marxan spatial prioritisation on case studies of four species from two regions, where matrices describe probability of larval dispersal  for rabbitfish and mudcrab in Southeast Sulawesi, Indonesia, and coral trout and sea cucumber in the Coral Triangle
#Dominic Muenzel bsdkm@leeds.ac.uk last updated 2022-08


#=#=#=#=#=#=#=#=#=#= 1. LIBRARIES #=#=#=#=#=#=#=#=#=#=

#for running
library(foreach)
library(doParallel)
library(igraph)
library(dplyr)
library(stringr)
library(tidyr)
library(metacapa)
library(prioritizr)
library(ggforce)
library(reshape)
library(tnet)
library(foreign)
library(ConnMatTools)
library(ggConvexHull)
library(lme4)
library(lmerTest)
library(dplyr)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(grid)
library(randomcoloR)
library(ggplot2)
library(gridExtra)
library(sparcl)
library(dendextend)
library(ggrepel)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(bayesplot)
library(rstanarm)
library(rstan)
library(insight)
library(bayestestR)
library(rnaturalearth)
library(rnaturalearthdata)

#=#=#=#=#=#=#=#=#=#= 2. SET UP WORKSPACE #=#=#=#=#=#=#=#=#=#=

#Clean environment
rm(list=ls())

#Set working directory to location of this script
oDir<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(oDir)

#Set directory of Marxan executable and outputs 
sMxDir<-paste0(oDir,"/MarxanData")

#Import custom functions 
source("Comparing spatial conservation prioritisation - functions.R")

#Create sub-directories if not already present. Note that the 'Marxan.exe' and 'input.dat' need to be placed in the directory /MarxanData for everything to run.
setup_directories()

#=#=#=#=#=#=#=#=#=#= 3. PRIORITISATION OF SIMULATED NETWORKS #=#=#=#=#=#=#=#=#=#=
#First this generates 100 seascapes with a near neighbour pattern for five dispersal distances and 100 seascapes with a small-world links pattern for five dispersal distances (distances of 50,100,150,200,250,a total of 500 for each pattern). Then it runs Marxan prioritisation scenarios for each seascape, using either site- or linkage-based connectivity methods to design reserve networks for six habitat targets (5,10,15,20,25,30%)  

#=#=#=#=#=#=#=#=#=#= 3.1. SET UP SEASCAPES #=#=#=#=#=#=#=#=#=#= 

#oQ is the near neighbour pattern, oQrw0.1 is the small-world links pattern. This performs the steps for each successively.
for(dispersalpattern in c("oQ","oQrw0.1")){
  
  #rr are the 100 replications of simulated seascapes
  for(rr in 1:100){ 
    
    #i are the five dispersal distances (50,100,150,200,250) 
   for(i in 1:5){    
   
     #Creating near neighbour seascapes
     
      #Parameters describing the seascape
      parameters<-generate_parameters(N=100,mapsize=5000,disp=c(50,100,150,200,250)[i],area="same",cost="same",discretize=4,toimport=F,prioritisation="marxan3",numreps=100,topN=10,fishtype="coraltrout",FIXEDCSM=100,runversion=rr,addname=dispersalpattern) 
      
      #Generate seascape
      seascape<-generate_seascape(parameters)
      
      #habitatMat is a matrix of seascape information consisting of the following columns: unique identifier of patch (puid), x,y of patch centroids, area of habitat, cost of patch, and site-specific graph metrics of patches
      habitatMat<-seascape[["habitatMat"]] 
      
      #habitatDispersal is the probability matrix of larval dispersal between patches
      habitatDispersal<-seascape[["habitatDispersal"]]
      
      #habitatDistance is a matrix of straight line interpatch distance
      habitatDistance<-seascape[["habitatDistance"]]
    
      #Creating small-world links seascapes
      if(dispersalpattern=="oQrw0.1"){
        
        #Edit the name used to save outputs
        parameters$figname<-edit_parametersfignameREWIRE(parameters,rewirePROB=0.1)
        
        #Rewire some connections in the probability matrix to create small-world patterns
        habitatDispersal<-rewirehabitatDispersal(habitatDispersal)
      }
      
      #Save the parameters, habitat matrix, distance matrix, dispersal matrix
      save_habitatMat(parameters,habitatMat,habitatDistance,habitatDispersal,parameters$figname)
      
      #Save dispersal probability as a numeric vector 
      save_disps(habitatDispersal,parameters)

#=#=#=#=#=#=#=#=#=#= 3.2. SPATIAL PRIORITISATION OF SEASCAPES #=#=#=#=#=#=#=#=#=#= 

      #Create a dataframe of scenarios for each connectivity method and habitat target, each row is a different scenario
      scenarios<-generate_scenarios(habitatMat,parameters)
      
      #Run Marxan for each scenario (goes row by row)
      run_scenarios(scenarios,habitatMat,habitatDispersal,parameters)
      
#=#=#=#=#=#=#=#=#=#= 3.3. ASSESSMENT OF RESERVE NETWORK PERFORMANCE #=#=#=#=#=#=#=#=#=#= 
      
      #Vector of lifetime egg production values (proxy for degradation of unprotected areas)
      allLEPS<-c(0,0.25,0.5,0.75)
      
      #Create four new names for files, adding the LEP value in the name
      fignameme1<-edit_parametersfignameLEP(parameters,allLEPS[1])
      fignameme2<-edit_parametersfignameLEP(parameters,allLEPS[2])
      fignameme3<-edit_parametersfignameLEP(parameters,allLEPS[3])
      fignameme4<-edit_parametersfignameLEP(parameters,allLEPS[4])
      
      #Combine four names
      allfignames<-c(fignameme1,fignameme2,fignameme3,fignameme4)
      
      #Each of the four LEP names corresponds to a different assumption about degradation or metapopulation contribution of unprotected habitats. For each of these, assess the performance of the reserve network
      for(y in 1:length(allfignames)){
         
         #For each of the top ten solutions of each scenario, calculate equilibrium settlement, using one of the four LEP assumptions
         for(runnumber in 1:parameters$topN){
           
            #Create a dataframe where results are saved. This stores characteristics such as number of reserves, cost, ID of planning unit solutions
            results<-run_scenarios_resultsdf(scenarios,habitatMat,parameters,run=runnumber)
            
            #Calculate equilibrium settlement using the DispersalPerRecruitModel from the package ConnMatTools
            results<-run_DPR(results,habitatMat,habitatDispersal,parameters,allLEPS[y])
            
            #Create a name to save the output, adding the id of the top run (runnumber)
            fignameprint<-paste0(allfignames[y],"_r",sprintf('%0.4d', runnumber))
            
            #Convert the results into a dataframe
            resultsDF<-results_df(results,parameters,scenarios)
            
            #Save the results dataframe
            save_resultsDF(resultsDF,fignameprint)
            
            #Edit the label for methods which are spatial dependency to include the CSM value
            resultsDF$detailtype[which(resultsDF$type=="spatial")]<-paste0(resultsDF$detailtype[which(resultsDF$type=="spatial")],"csm",sprintf('%0.3f',resultsDF$csm[which(resultsDF$type=="spatial")]))
            
            #Dataframe to find the CSM value which achieves the highest equilibrium settlement at a similar cost to baseline runs. CSM is the connectivity strength modifier, equivalent to boundary length modifier (BLM) if using bound.dat in Marxan to define physically adjacent boundaries between planning units
            bestcsmfortarget<-data.frame("target"=unique(resultsDF$target),"csm"=0,"index"=0)
              
               #Find best CSM at each unique % target
               for(i in 1:length(unique(resultsDF$target))){  #i=1
                  
                  #Vector of all % targets
                  alltargets<-unique(resultsDF$target)
                  
                  #Subset results dataframe to those matching the target
                  workingresultsDF<-subset(resultsDF,resultsDF$target==alltargets[i])
                  
                  #Subset the results dataframe to those running spatial dependency method
                  workingresultsDF<-subset(workingresultsDF,workingresultsDF$type=="spatial")
                  
                  #Find the best CSM at the target that has highest equilibrium settlement at similar cost
                  bestcsmfortarget$csm[i]<-workingresultsDF$csm[which(workingresultsDF$DPRsettlers==max(workingresultsDF$DPRsettlers[which(workingresultsDF$HitTarget<1.01)]))[1]]
               }
              
            #Paste the best CSM back into the results dataframe (resultsDF) for all scenarios
            for(i in 1:nrow(bestcsmfortarget)){
              
               #Subset of the resultsDF matching the best CSM
               resultsDFsub<-resultsDF[which(resultsDF$csm==bestcsmfortarget$csm[i]),]
               
               #Subset of the resultsDF subset matching the target
               resultsDFsub<-resultsDFsub[which(resultsDFsub$target==bestcsmfortarget$target[i]),]
               
               #Index to match to correct name
               bestcsmfortarget$index[i]<-which(resultsDFsub$Name==resultsDF$Name)
            }
            
            #Save a subset of the resultsDF where only one csm value per scenario is retained. This is the best CSM value that has similar cost but highest equilibrium settlement. The other csm value scenarios are not of interest.
            resultsDFsubset<-resultsDF[c(bestcsmfortarget$index,which(resultsDF$type=="metric")),]
            resultsDFsubset$detailtype[grep("spatial",resultsDFsubset$detailtype)]<-"spatial"

            #Append "_subset" to the name of the file so the original entire resultsDF is not overwritten
            fignameprintsub<-paste0(fignameprint,"_subset")
            
            #Save the subset of the results dataframe (resultsDF)
            save_resultsDF(resultsDFsubset,fignameprintsub)
         }
      }
   }
  }
}

#=#=#=#=#=#=#=#=#=#= 4. PRIORITISATION OF CASE STUDIES #=#=#=#=#=#=#=#=#=#=
#This imports the dispersal probability matrices for the four species and then  runs Marxan prioritisation scenarios for each species, using either site- or linkage-based connectivity methods to design reserve networks for six habitat targets (5,10,15,20,25,30%) 

#=#=#=#=#=#=#=#=#=#= 4.1. IMPORT CASE STUDIES  #=#=#=#=#=#=#=#=#=#= 

#Cuke is sea cucumber, Trout is coral trout. These are both Coral Triangle species. Mudcrab and rabbitfish are Southeast Sulawesi species. This runs the steps for each species one at a time.
for(casestudy in c("Cuke","Trout","Mudcrab","rabbitfish")){  

  #Parameters describing the species
  parameters<-generate_parameters(N=100,mapsize=5000,disp=10,area="same",cost="same",discretize=10,toimport=T,prioritisation="marxan3",numreps=100,topN=10,fishtype=casestudy,FIXEDCSM=100,runversion=1,addname="")
  
  #Import the dispersal probability matrix and the x,y centroids of the planning units
  seascape<-generate_seascape(parameters)
  habitatMat<-seascape[["habitatMat"]]
  habitatDispersal<-seascape[["habitatDispersal"]]
  habitatDistance<-seascape[["habitatDistance"]]
  
#=#=#=#=#=#=#=#=#=#= 4.2. SPATIAL PRIORITISATION OF CASE STUDIES #=#=#=#=#=#=#=#=#=#= 
  
  #Create a dataframe of scenarios for each connectivity method and habitat target, each row is a different scenario
  scenarios<-generate_scenarios(habitatMat,parameters)
  
  #Run Marxan for each scenario (goes row by row)
  run_scenarios(scenarios,habitatMat,habitatDispersal,parameters)
  
#=#=#=#=#=#=#=#=#=#= 4.3. ASSESSMENT OF RESERVE NETWORK PERFORMANCE #=#=#=#=#=#=#=#=#=#= 
  
  #Vector of lifetime egg production values (proxy for degradation of unprotected areas)
  allLEPS<-c(0,0.25,0.5,0.75)
  
  #Adjust the names of the files, adding "LEP" and the value of LEP. This ensures files are not overwitten. 
  allfignames<-c(paste0("importLEP00_",casestudy,"_marxan3"),paste0("importLEP25_",casestudy,"_marxan3"),paste0("importLEP50_",casestudy,"_marxan3"),paste0("importLEP75_",casestudy,"_marxan3"))
  
  #Each of the four LEP names corresponds to a different assumption about degradation or metapopulation contribution of unprotected habitats. For each of these, assess the performance of the reserve network
  for(y in 1:length(allfignames)){
    
    #For each of the top ten solutions of each scenario, calculate equilibrium settlement, using one of the four LEP assumptions
     for(runnumber in 1:parameters$topN){
       
        #Create a dataframe where results are saved. This stores characteristics such as number of reserves, cost, ID of planning unit solutions
        results<-run_scenarios_resultsdf(scenarios,habitatMat,parameters,run=runnumber)
      
        #Calculate equilibrium settlement using the DispersalPerRecruitModel from the package ConnMatTools
        results<-run_DPR(results,habitatMat,habitatDispersal,parameters,allLEPS[y])
        
        #Create a name to save the output, adding the id of the top run (runnumber)
        fignameprint<-paste0(allfignames[y],"_r",sprintf('%0.4d', runnumber))
        
        #Convert the results into a dataframe
        resultsDF<-results_df(results,parameters,scenarios)
        
        #Save the results dataframe
        save_resultsDF(resultsDF,fignameprint)
        
        #Edit the label for methods which are spatial dependency to include the CSM value
        resultsDF$detailtype[which(resultsDF$type=="spatial")]<-paste0(resultsDF$detailtype[which(resultsDF$type=="spatial")],"csm",sprintf('%0.3f',resultsDF$csm[which(resultsDF$type=="spatial")]))
        
        #Dataframe to find the CSM value which achieves the highest equilibrium settlement at a similar cost to baseline runs. CSM is the connectivity strength modifier, equivalent to boundary length modifier (BLM) if using bound.dat in Marxan to define physically adjacent boundaries between planning units
        bestcsmfortarget<-data.frame("target"=unique(resultsDF$target),"csm"=0,"index"=0)
        
        #Find best CSM at each unique % target
        for(i in 1:length(unique(resultsDF$target))){  
          
          #Vector of all % targets
          alltargets<-unique(resultsDF$target)
          
          #Subset results dataframe to those matching the target
          workingresultsDF<-subset(resultsDF,resultsDF$target==alltargets[i])
          
          #Subset the results dataframe to those running spatial dependency method
          workingresultsDF<-subset(workingresultsDF,workingresultsDF$type=="spatial")
            
          #Find the best CSM at the target that has highest equilibrium settlement at similar cost
          bestcsmfortarget$csm[i]<-workingresultsDF$csm[which(workingresultsDF$DPRsettlers==max(workingresultsDF$DPRsettlers[which(workingresultsDF$HitTarget<1.15)]))[1]]  
         }
         
        #Paste the best CSM back into the results dataframe (resultsDF) for all scenarios
        for(i in 1:nrow(bestcsmfortarget)){
          
          #Subset of the resultsDF matching the best CSM
          resultsDFsub<-resultsDF[which(resultsDF$csm==bestcsmfortarget$csm[i]),]
          
          #Subset of the resultsDF subset matching the target
          resultsDFsub<-resultsDFsub[which(resultsDFsub$target==bestcsmfortarget$target[i]),]
          
          #Index to match to correct name
          bestcsmfortarget$index[i]<-which(resultsDFsub$Name==resultsDF$Name)
        }
        
        #Save a subset of the resultsDF where only one csm value per scenario is retained. This is the best CSM value that has similar cost but highest equilibrium settlement. The other csm value scenarios are not of interest.
        resultsDFsubset<-resultsDF[c(bestcsmfortarget$index,which(resultsDF$type=="metric")),]
        resultsDFsubset$detailtype[grep("spatial",resultsDFsubset$detailtype)]<-"spatial"
        
        #Append "_subset" to the name of the file so the original entire resultsDF is not overwritten
        fignameprintsub<-paste0(fignameprint,"_subset")
        
        #Save the subset of the results dataframe (resultsDF)
        save_resultsDF(resultsDFsubset,fignameprintsub)
     }
  }
}
#=#=#=#=#=#=#=#=#=#=#=#= END =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

