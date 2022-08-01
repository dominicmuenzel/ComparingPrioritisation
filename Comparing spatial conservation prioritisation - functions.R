#Comparing spatial conservation prioritisation methods using site- versus linkage-based connectivity
#This is the underlying code for the paper found at ____. Section 3 is for generating simulated seascapes and running Marxan spatial prioritisation on these. Prioritisation is run without using connectivity (baseline), using connectivity-based features methods (where metrics of connectivity are calculated for each planning unit and targets are set for these), or using linkage-based methods (the spatial dependency method described in Beger et al. 2010 in Conservation Letters). Section 4 is for running Marxan spatial prioritisation on case studies of four species from two regions, where matrices describe probability of larval dispersal  for rabbitfish and mudcrab in Southeast Sulawesi, Indonesia, and coral trout and sea cucumber in the Coral Triangle
#Dominic Muenzel bsdkm@leeds.ac.uk last updated 2022-04

#=#=#=#=#=#=#=#=#=#=#=#= FUNCTIONS =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

## Create directories if not already present
setup_directories<-function(){
  
  #Set up the results directory
  if(!dir.exists("results")){dir.create("results")}
  
  #Set up the habitatMat directory. Here the habitatMat matrices are stored when they are created for either simulated seascapes or case studies. The habitatMat contains rows for each planning unit and gives x,y coordinates as well as cost and area
  if(!dir.exists("results/habitatMat")){dir.create("results/habitatMat")}
  
  #Set up the resultsDF directory. This stores the resultsDF csv files which will be produced after having run Marxan, when reserve networks are asessed for their performance. 
  if(!dir.exists("results/habitatMat")){dir.create("results/resultsDF")}
  
  #Set up the stan directory. This stores the model coefficients from the Bayesian stan models. 
  if(!dir.exists("results/habitatMat")){dir.create("results/stan")}
  
  #Set up the MarxanData directory. The marxan.exe executable, as well as the input.dat file need to be placed in here
  if(!dir.exists("MarxanData")){dir.create("MarxanData")}
  
  #Set up the input directory. This is where Marxan input files will be stored (pu.dat, spec.dat, puvspr2.dat, assymbound.dat)
  if(!dir.exists("MarxanData/input")){dir.create("MarxanData/input")}
  
  #Set up the output directory. This is where Marxan output files will be stored
  if(!dir.exists("MarxanData/output")){dir.create("MarxanData/output")}
}

## Generate n colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Opposite of %in%
`%notin%` <- Negate(`%in%`)

## Calculate the Euclidean distance between patches, used for generating dispersal matrix
patch_distance<-function(habitatMat){
  
  #Get the x,y coordinates from habitatMat matrix
  xy<-habitatMat[,c("x","y")]
  
  #Create a distance matrix from the x,y coordinates 
  xy<-as.matrix(dist(xy,diag=T,upper=T))
  
  #Return matrix of distance
  return(as.matrix(xy,nrow=nrow(habitatMat),ncol=nrow(habitatMat)))
}

## Calculate distance between patches as edge to edge, used to check for overlaps
patch_distance_edge_to_edge_nodiag<-function(habitatMat,radiuses){
  
  #Get the x,y coordinates from habitatMat matrix
  xy<-habitatMat[,c("x","y")]
  
  #Create a distance matrix from the x,y coordinates 
  xy<-as.matrix(dist(xy,diag=T,upper=T))
  
  #Subtract radius from columns aka destination
  xy <- sweep(xy, 2, habitatMat[,"radius"], "-")  
  
  #Subtract radius from rows
  xy <- sweep(xy, 1, habitatMat[,"radius"], "-")      
  
  #Remove diagonals
  diag(xy)<-0
  
  #Set NAs to zero
  xy[is.na(xy)]<-0
  
  #Return matrix of distance
  return(as.matrix(xy,nrow=nrow(habitatMat),ncol=nrow(habitatMat)))
}

## Parameters describing the seascape, including number of patches (N), size of map (mapsize), dispersal decay rate (disp), whether area and cost are uniform (area, cost), number of quantiles for discretising metrics (discretize), whether seascape is simulated or an imported case study (toimport), which prioritisation is used although this is not implemented (prioritisation), number of Marxan replicates (numreps), number of top solutions to analyse (topN), name of fish for case studies (fishtype), 
generate_parameters<-function(N,mapsize,disp,area,cost,discretize,toimport,prioritisation,numreps,topN,fishtype,FIXEDCSM,runversion,addname){
  
  #File name if seascape is simulated
  if(isFALSE(toimport)){
    
    #Unique name based on parameters to avoid overwriting each time
    figname<-paste0("N",N,addname,"_","_d",discretize,"_","disp",sprintf('%0.4d',round(disp,4)),"_",prioritisation,"_rv",sprintf('%0.3d',runversion)) 
    
  #Parameters and file name if seascape is an imported case study  
  } else {
    if(isTRUE(toimport)){
      
      #Rabbitfish case study from Southeast Sulawesi
     if(fishtype %in% c("rabbitfish")){
      
       #Number of reef patches
       N<-487
       
       #Name of files
      figname<-paste0("import_",fishtype,"_",prioritisation)
     } else {
       
       #Sea cucumber and coral trout case study from Coral Triangle
         if(fishtype %in% c("Cuke","Trout")){
           
           #Number of reef patches
           N<-425
           
           #Name of file
           figname<-paste0("import_",fishtype,"_",prioritisation)
           
           #Mudcrab case study from Southeast Sulawesi
         } else {
           if(fishtype%in% c("Mudcrab")){
             
             #Number of mangrove patches
             N<-216
             
             #Name of file
             figname<-paste0("import_",fishtype,"_",prioritisation)
         }
       }
     }
    }
  }

#List of parameter values
  parameters<-list("N"=N,"mapsize"=mapsize,"disp"=disp,"area"=area,"cost"=cost,"discretize"=discretize,"toimport"=toimport,"figname"=figname,"prioritisation"=prioritisation,"numreps"=numreps,"topN"=topN,"fishtype"=fishtype,"FIXEDCSM"=FIXEDCSM,"runversion"=runversion,"addname"=addname)

  return(parameters)  
}

## create a seascape with randomly selected x and y and patch size, patches do not overlap
habitat_spawn_loop<-function(parameters){
  
  #Number of patches to simulate
  N<-parameters$N 
  
  #Size of extent
  mapsize<-parameters$mapsize  
  
  #Exponential decay rate parameter
  disp<-parameters$disp

  #Column names of habitat matrix
  habitatMatcolumns<-c("puid","x","y","area","radius","disp","cost")
  
  #Create blank matrix
  habitatMat<-matrix(NA,nrow=N,ncol=length(habitatMatcolumns))  #create a blank data frame
  colnames(habitatMat)<-habitatMatcolumns
  
  #Planning unit id is a unique identifier for each patch
  habitatMat[,"puid"]<-1:N 
  
  #Uniform area for all patches
  if(parameters$area=="same"){
    habitatMat[,"area"]<-50
  }

  #Radius of each patch
  habitatMat[,"radius"]<-round(sqrt(habitatMat[,"area"]/pi),4)

  #Generate non-overlapping patches
  for (i in 1:N){  
    
    repeat{
      
      #random sample from uniform pdf
      habitatMat[i,"x"]<-round(runif(1,0,mapsize),2)  
      
      #random sample from uniform pdf
      habitatMat[i,"y"]<-round(runif(1,0,mapsize),2)
      
      #does the ith habitat fall outside the area?
      abovex<-habitatMat[i,"x"]+habitatMat[i,"radius"]>mapsize
      underx<-habitatMat[i,"x"]-habitatMat[i,"radius"]<0
      abovey<-habitatMat[i,"y"]+habitatMat[i,"radius"]>mapsize
      undery<-habitatMat[i,"y"]-habitatMat[i,"radius"]<0
      
      #index to test any patches outside
      outside<-sum(abovex,underx,abovey,undery)  
      
      #does the ith habitat overlap with one of the preceding habitats?
      overlap<-patch_distance_edge_to_edge_nodiag(habitatMat)
      
      #index to test any overlaps
      overlap<-sum(overlap<0) 
      
      if(outside==0 & overlap==0){
        break
      }
    }
  }

  return(habitatMat)
}

## Take the dispersal probability matrix and create the Marxan assymbound input file giving the penalty between patches. This is used in the linkage-based method of using connectivity, equivalent to boundary.dat in implementations where physical shared boundary of planning units is used to create spatially clustered solutions. 
generate_assymbound<-function(habitatDispersal){
  
  #igraph graph object from dispersal matrix
  g  <- graph.adjacency(as.matrix(habitatDispersal),weighted=TRUE)
  
  #Edge list with weights of graph
  df <- get.data.frame(g)
  
  #Remove self connections
  df<- df[  -which(df$from==df$to) , ]   
  
  return(df)
  
}

## Turns the dispersal matrix into an igraph graph
habitatMat_to_graph<-function(habitatMat,habitatDispersal){
  
  #igraph graph from dispersal matrix
  g  <- graph_from_adjacency_matrix(habitatDispersal,weighted=TRUE,mode="directed")
  
  #Assign attribute of area with values taken from habitat matrix
  V(g)$area<-habitatMat[,"area"]
  return(g)
}

## Save the habitatMat matrix, distance matrix, dispersal matrix, and parameters
save_habitatMat<-function(parameters,habitatMat,habitatDistance,habitatDispersal,figname){
  setwd(oDir)
  save(parameters,habitatMat,habitatDistance,habitatDispersal,file=paste0("results/habitatMat/",figname,"_habitatMat.RData"))  
}

## Calculate graph metrics of each patch and add columns named after metrics to the habitatMat matrix
add_metrics<-function(habitatMat,habitatDispersal,parameters){

  #Flow matrix is calculated as the imported probability matrix row multiplied by area of patch 
  flowMat<-sweep(habitatDispersal, MARGIN=1, habitatMat[,"area"], `*`)
  
  #Migration matrix is the flow matrix divided by column sums
  migMat<-t(t(flowMat)/colSums(flowMat))
        
  #Dispersal probability matrix as an igraph graph  
  gprob<-habitatMat_to_graph(habitatMat,habitatDispersal)
  
  #Number of vertices
  N<-gorder(gprob)
  
  #igraph graph of flow matrix
  gflow<-habitatMat_to_graph(habitatMat,flowMat)
  
  #igraph graph of migration matrix
  gmig<-habitatMat_to_graph(habitatMat,migMat)
  
  #Get rid of NAs in the migration matrix
  migMat[is.na(migMat)]<-0

  #Calculate Google Page Rank
  page<-page_rank(gflow)$vector

  #Calculate Betweenness Centrality, igraph betweenness function can't cope with very small values, so set very small values as the same
  gminup<-gflow
  E(gminup)$weight[which(E(gminup)$weight<1e-10)]<-1e-10
  betweenness<-betweenness(gminup,normalized=F) 
  
  #Calculate Eigenvector Centrality
  eigenvector<-eigen_centrality(gmig)$vector 
  
  #Calculate Outflux  
  outflux<-sapply(1:N, function(i) sum(E(gflow)[from(i)]$weight)-as.numeric(diag(flowMat)[i]))

  #Calculate Influx  
  influx<-sapply(1:N, function(i) sum(E(gflow)[to(i)]$weight)-diag(flowMat)[i])
  
  #Calculate In Degree
  indegree <- degree(gflow,mode="in",loops=T)
  
  #Calculate Out Degree
  outdegree <- degree(gflow,mode="out",loops=T)
  
  #Calculate Degree
  degree <- degree(gflow,mode="all",loops=T)

  #Calculate Local retention
  localretention<-diag(habitatDispersal)

  #Add the newly calculated metrics to the habitatMat matrix
  habitatMat_add<-cbind(habitatMat[,1:7],page,betweenness,eigenvector,outflux,influx,indegree,outdegree,localretention)
 
  #Set any NAs to zero
  habitatMat_add[which(is.na(habitatMat_add))]<-0
  
  #Rescale all metrics going from 0-100
  for(i in 8:ncol(habitatMat_add)){
    habitatMat_add[,i]<-scales::rescale(habitatMat_add[,i],to=c(0,100))
    
  }
  return(habitatMat_add)
}

## This edits the spec.dat Marxan input file. Each time a scenario is run, the spec.dat is edited to have the correct habitat target.
adjust_spec.dat<-function(k,relativetargets,scenarios,parameters){

  #Set working directory
  setwd(oDir)
  
  #Read the spec.dat input file located in the /MarxanData/input directory
  x<-read.table(paste0("MarxanData/input/",scenarios[,"specdat"][k],".dat"),stringsAsFactors = F,header=T,sep=",")
  
  #If running the spatial dependency scenario 
  if(scenarios[,"scenario_name"][k]=="spatial_dependency"){
    
    #Take the value for habitat target "prop" from the scenarios dataframe 
    x$prop[which(x$name=="Habitat")]<-scenarios[,"target"][k]
        
    #If running the site-based connectivity metric method
      } else {
        
          #Set the target for the metric at 0.5 of the habitat target for in- and out-flow, higher values result in cost being greater than baseline runs 
          if(scenarios[,"scenario_name"][k] %in% c("outflux","influx")){
            x$prop<-relativetargets[k,]*0.5
          } else {
            
          #Set the target for the metric at 0.75 for other metrics, higher values result in cost being greater than baseline runs
            x$prop<-relativetargets[k,]*0.75
          }
          
          #Set the habitat target % labelled as "prop" in spec.dat with the value taken from the scenarios dataframe
          x$prop[which(x$name=="Habitat")]<-scenarios[,"target"][k]
        
      }
    
    #Save the spec.dat file in the /MarxanData/input/ directory with the correct habitat and metric targets 
    write.table(x,file=paste0("MarxanData/input/",scenarios[,"specdat"][k],".dat"),row.names = F, quote = F,sep=",")
}

## Function of a negative exponential dispersal kernel
dispersal_negexp <- function(alpha) {
  stopifnot(is.numeric(alpha))
  stopifnot(!is.na(alpha))
  stopifnot(length(alpha) == 1)
  stopifnot(alpha > 0)
  
  function(d) {
    exp(-alpha * abs(d))
  }
}

## Create a dataframe of prioritisation scenarios, each row is a connectivity method (site- or linkage-based) and a given habitat target (%)
generate_scenarios<-function(habitatMat,parameters){
  
  #Vector of habitat protection targets (5,10,15,20,25,30%)
  targets<-seq(from=0.05,to=.3,by = .05)   
  
  #Vector of connectivity strength modifier values (CSM) used to calibrate the spatial dependency method
  allcsms<-c(0.1,0.5,1,5,seq(from=10,to=300,by=10))
  
  #Create a dataframe of scenarios with appropriate column names.  
  scennames<-colnames(habitatMat)[8:ncol(habitatMat)]
  columnnames<-c("uniqid","scenario_name","pudat","puvspr2dat","specdat","assymbounddat","name","csm","target")
  scenariomat<-data.frame(matrix(0,ncol=length(columnnames),nrow=c(length(targets)+length(targets)*length(allcsms)+((length(scennames))*length(targets)))))
  colnames(scenariomat)<-columnnames
  
  #Give each scenario a unique id and name (baseline, spatial, or name of connectivity metric)
  scenariomat[,"uniqid"]<-1:nrow(scenariomat)
  scenariomat[,"scenario_name"]<-c(rep("baseline",length(targets)),rep("spatial",(length(targets))*length(allcsms)),rep((scennames),length(targets)))
  
  #Name the input files used by Marxan after the name of the seascape, including the pu.dat, puvspr2.dat, spec.dat, and assymbound.dat
  scenariomat[,"pudat"]<-paste0("pu_",parameters$figname)
  scenariomat[,"puvspr2dat"]<-paste0("puvspr2_",parameters$figname)
  scenariomat[,"specdat"]<-paste0("spec_",parameters$figname)

  scenariomat[,"assymbounddat"]<-c(rep(NA,(length(targets))),rep(paste0("assymbound_",parameters$figname),(length(targets)*length(allcsms))),rep(NA,(length(targets))*length(scennames)))

  #Set the connectivity strength modifier (csm) for each scenario. Only spatial dependency uses a csm, baseline and connectivity metrics do not. 
  scenariomat[,"csm"]<-c(rep(NA,(length(targets))),rep(allcsms,each=length(targets)),rep(NA,length(rep(scennames,length(targets)))))
  
  #Set the habitat protection target for each scenario
  scenariomat[,"target"]<-c(targets,rep(targets,length(allcsms)),rep(targets,each=length(scennames)))
  
  #Give each scenario a unique name, combining the % target with the method of using connectivity
  scenariomat[,"name"]<-paste0(scenariomat[,"scenario_name"],"_t",as.character(format(as.numeric(scenariomat[,"target"]),nsmall = 2)))
  scenariomat[which(!is.na(scenariomat[,"csm"])),"name"]<-paste0(scenariomat[which(!is.na(scenariomat[,"csm"])),"name"],"_csm",sprintf('%0.4f',as.numeric(scenariomat[which(!is.na(scenariomat[,"csm"])),"csm"])))
  
  return(scenariomat)
}

## Add uniform cost to habitatMat matrix
add_cost<-function(habitatMat,parameters){
  
  #If using uniform cost then all patches have same cost
  if(parameters$cost=="same"){
    habitatMat[,which(colnames(habitatMat)=="cost")]<-15
  }
    
  return(habitatMat)
}

## Creates a probability matrix of larval dispersal based on a negative exponential dispersal kernel given a decay parameter and x,y, coordinates of patches
dispersal<-function(habitatMat,parameters,localrecruitment){
  
  # a negative exponential function
  disp<-dispersal_negexp(1/parameters$disp)
  
  #distance from centroid to centroid
  habitatDistance<-patch_distance(habitatMat)
  
  #distance passed through dispersal function
  habitatDispersal<-disp(habitatDistance) # 
  
  #add local recruitment
  if(isTRUE(localrecruitment)){
    habAREA<-runif(parameters$N,min=1,max=100)
    diag(habitatDispersal) <- diag(habitatDispersal)*sqrt(habAREA)/mean(sqrt(habAREA))
  }

  #Assume that probabilities smaller than 1e-6 are ecologically insignificant  
  habitatDispersal[which(habitatDispersal<1e-6)]<-0
  
  return(habitatDispersal)

}

## Save the dispersal as numeric vector. Used to create box plots of dispersal
save_disps<-function(habitatDispersal,parameters){
  
  #Matrix to remove diagonals
  habitatDispersalnodiag<-habitatDispersal
  
  #Remove diagonals
  diag(habitatDispersalnodiag)<-0
  
  #Numeric vector of non diagonal dispersals
  nondiagonaldisps<-as.numeric(habitatDispersalnodiag)
  
  #Numeric vector of diagonal dispersals
  diagonaldisps<-as.numeric(diag(habitatDispersal))
  
  #Save both vectors
  save(nondiagonaldisps,diagonaldisps,file=paste0("results/habitatMat/",parameters$figname,"_disps.RData"))
}

## Edit the file name to add the rewire probability of 0.1 used to create small-world links pattern
edit_parametersfignameREWIRE<-function(parameters,rewirePROB){
  
  #Split the file name by "_"
  parametersfignamedit<-  str_split(parameters$figname,"_")[[1]]
  
  #Add "rw0.1" to the file name
  parametersfignamedit[1]<-paste0(parametersfignamedit[1],"rw",rewirePROB)
  
  #Collapse new file name
  parametersfignamedit<-paste(parametersfignamedit,collapse="_")
  
  return(parametersfignamedit)
}

## Edit the file name to add the value for lifetime egg production (LEP)
edit_parametersfignameLEP<-function(parameters,LEP){
  
  #Split the name of the file
  parametersfignamedit<-  str_split(parameters$figname,"_")[[1]]
  
  #Add "LEP" and LEP value to file name
  parametersfignamedit[1]<-paste0(parametersfignamedit[1],"LEP",sprintf('%0.2d',LEP*100))
  
  #Collapse file name
  parametersfignamedit<-paste(parametersfignamedit,collapse="_")
  
  #Return modified file name
  return(parametersfignamedit)
}

## Generates the dataframe of results by reading in the appropriate Marxan output files. This contains sum of reserves, cost, id of planning units, etc.
run_scenarios_resultsdf<-function(scenarios,habitatMat,parameters,run){
  
  #Set working directory
  setwd(oDir)
  
  #Give colum names
  columnames<-c("Name","SolutionCost","SumCost","SumArea","NrReserves","HitTarget","DPRrecruits","DPRsettlers","DPRsettlersfished","DPRsettlersnotfished","DPRrecruitsfished","DPRrecruitsnotfished",paste0("s",1:nrow(habitatMat)))
  
  #Create an empty dataframe and give columns names
  solutionsmat<-data.frame(matrix(NA,nrow=nrow(scenarios),ncol=length(columnames)))
  colnames(solutionsmat)<-columnames
  
  #For each scenario
  for(k in 1:nrow(scenarios)){   #k=1
    
    #Directory of Marxan outputs
    thisdir<-paste0(sMxDir,"/output/")
    
    #Read in the summary table for the scenario which is produced by Marxan
    sumtable<-read.table(paste0(thisdir,parameters$figname,"/",parameters$figname,"_",scenarios[k,"name"],"_sum.csv"),sep=",",header=T)
    
    #Re-order summary table to find best solutions, which have lowest scores
    sumtable<-sumtable[order(sumtable$Score),]
    
    #Get ID of the solution and its cost 
    topNsolution <- sumtable$Run_Number[run]
    topNsolutioncost<-sumtable$Cost[run]
    
    #Read the run file which Marxan produces for the specific solution 
    topNsolutioncsv<-read.table(paste0(thisdir,parameters$figname,"/",parameters$figname,"_",scenarios[k,"name"],"_","r",sprintf('%0.5d', topNsolution),".csv"),sep=",",header=T)

    #Order the run file in descending order of planning unit ID
    topNsolutioncsv<-topNsolutioncsv[order(topNsolutioncsv$PUID),]
    
    #Vector of planning unit designations (0 not designated, 1 designated)
    topNsolutioncsvonly<-as.numeric(topNsolutioncsv$SOLUTION)
    
    #Get name, cost, amount of habitat protected, number of reserves, habitat target set in the prioritisation, and proportion of habitat target achieved for the scenario
    solname<-paste0(parameters$figname,"_",scenarios[k,"name"])
    solcost<-topNsolutioncost
    sumcost<-sum(habitatMat[,"cost"][which(topNsolutioncsvonly==1)])
    sumarea<-sum(habitatMat[,"area"][which(topNsolutioncsvonly==1)])
    reserves<-sum(topNsolutioncsvonly)
    target<-sum(habitatMat[,"area"])*unlist(as.numeric(scenarios[k,"target"]))
    hittarget<-round(sumarea/target,2)
    
    #Put together solution dataframe 
    solutionsmat[k,]<-c(solname,solcost,sumcost,sumarea,reserves,hittarget,0,0,0,0,0,0,topNsolutioncsvonly)   
    }
  
  #Shorten cost, sum cost, and sum area to three decimal points
  for(g in c(2,3,4)){
    solutionsmat[,g]<-round(as.numeric(solutionsmat[,g]),3)
  }

return(solutionsmat)  
}

## This generates the standard input files for Marxan including pu.dat, spec.dat, puvspr2.dat, assymbound.dat. It also returns a matrix giving the relative targets for each of the metric quantiles used in the connectivity-based feature method
generate_marxan_inputfiles<-function(scenarios,habitatMat,habitatDispersal,parameters){
  
  #Set arbitrary beginning target, this doesn't matter as it is replaced
  habitattarget<-0.1
  
  #Get the connectivity metrics used in the site-based method from the scenario dataframe
  allmetrics<-unique(scenarios[,"scenario_name"])[which(unique(scenarios[,"scenario_name"])%notin%c("baseline","spatial"))]
  
  #Dataframe of planning unit id and cost
  x<-data.frame("id"=habitatMat[,"puid"],"cost"=habitatMat[,"cost"])
  
  #Number of quantiles to discretise the continuous connectivity metrics into discrete conservation features
  discretizationby<-parameters$discretize
  discretizationprobability<-seq(from=0,to=1,by=1/discretizationby)
  
  #Combination of all connectivity metrics and their quantiles
  allmetricscombos<-paste0(rep(allmetrics,each=discretizationby),"_",0:(discretizationby-1))

  #If running the connectivity metrics method, then the conservation features consist of the habitat and the connectivity features split into quantiles
  if(length(allmetrics)>1){
  features<-data.frame("id"=c(1:(length(allmetricscombos)+1)),"name"=c("Habitat",allmetricscombos))
  } else {
    
  #If running either the baseline or spatial dependency method, conservation features consist only of habitat
    features<-data.frame("id"=c(1),"name"=c("Habitat"))
  }
  
  #Dataframe containing planning unit id, species id, and amount of habitat
  rij<-data.frame("pu"=habitatMat[,"puid"],"species"=1,"amount"=habitatMat[,"area"])
  
  #If running the connectivity features method
  if(length(allmetrics)>1){
    
  #For each of the connectivity metrics
  for(i in 1:length(allmetrics)){  #i=1
    
    #List of unique species id, based on the combination of number of connectivity metrics and number of quantiles
    spec_0s<-seq(from=2,to=length(allmetricscombos)+2,by=discretizationby)[i]
    
    #Get the values of the metrics from the habitatMat matrix
    metricvalues<-unlist(as.numeric(habitatMat[,allmetrics[i]]))
    
    #If there are zeroes, it throws off the quantiles and these PUs are never selected. Instead give them a small value
    metricvalues[which(metricvalues==0)]<-runif(length(which(metricvalues==0)),min=0.00001,max=0.0001)
    
    #If the amount of the connectivity metric is the same across all planning units, each planning unit has the same value 
    if(var(metricvalues)==0){
      
      #Dataframe of planning unit ids and amount of each connectivity metric feature (termed "species")
      rijfeat<-data.frame("pu"=sample(unlist(as.numeric(habitatMat[,"puid"])),nrow(habitatMat),replace=F),"species"=rep(seq(from=spec_0s,to=(spec_0s+discretizationby-1)),each=10),"amount"=1)
      
      #Combine the connectivity metric dataframe with the habitat dataframe
      rij<-rbind(rij,rijfeat)
      
      #If the amount of the connectivity metric differs across planning units
    } else {
      
      #This is the highest quantile of the connectivity metric
      spec_0<-which(metricvalues>quantile(metricvalues,probs=discretizationprobability)[discretizationby])
      
      if(length(spec_0)>0){
        
        #If there is no feature in the highest quantile, then all features are the same
        rijfeat<-data.frame("pu"=spec_0,"species"=spec_0s,"amount"=1)
        
        #Combine the connectivity metric dataframe with the habitat dataframe
        rij<-rbind(rij,rijfeat)
      }
      
      #Add the remaining quantiles of the connectivity metric and the amount of each quantile in each planning unit
      for(k in 1:(discretizationby-2)){ #
        
        #Upper limit of the ith quantile
        upperquants<-seq(from=discretizationby,to=2)[k]
        
        #Lower limit of the ith quantile
        lowerquants<-upperquants-1
        
        #ID of planning units which contain this quantile
        pus<-which(metricvalues<=quantile(metricvalues,probs=discretizationprobability)[upperquants] & metricvalues>quantile(metricvalues,probs=discretizationprobability)[lowerquants])
        
        #If there is at least one planning unit containing this quantile, add to the conservation features dataframe
        if(length(pus)>0){
          rijfeat<-data.frame("pu"=pus,"species"=spec_0s+k,"amount"=1)
          rij<-rbind(rij,rijfeat)
        }
      }
      
      #This is the lowest quantile of the connectivity metric
      spec_12<-which(metricvalues<=quantile(metricvalues,probs=discretizationprobability)[2])
      
      #If there is at least one planning unit containing this quantile, add to the conservation features dataframe
      if(length(spec_12)>0){
        rijfeat<-data.frame("pu"=spec_12,"species"=spec_0s+k+1,"amount"=1)
        rij<-rbind(rij,rijfeat)
      }
    }
  }
  }
 
  #Create the puvspr2.dat Marxan input file, which gives the amount of each feature in each planning unit
  puvspr2.dat<-data.frame("species"=rij$species,"pu"=rij$pu,"amount"=rij$amount) 
  
  #Order the puvspr2.dat by planning unit
  puvspr2.dat<-puvspr2.dat[order(puvspr2.dat$pu),]
  
  #Write the puvpsr2.dat file in the /MarxanData/input/ directory
  write.table(puvspr2.dat,file=paste0("MarxanData/input/puvspr2_",parameters$figname,".dat"),sep=",",row.names = F,quote=F)
 
  #Create the pu.dat Marxan input file with columns for planning unit id, cost, status
  pu.dat<-data.frame("id"=habitatMat[,"puid"],"cost"=habitatMat[,"cost"],"status"=0)

  #Write the pu.dat file in the /MarxanData/input/ directory
  write.table(pu.dat,file=paste0("MarxanData/input/pu_",parameters$figname,".dat"),sep=",",row.names = F,quote=F)
  
  #Create the spec.dat Marxan input file with columns for species id, prop (target %), spf (species penalty factor), name
  spec.dat<-data.frame("id"=features$id,"prop"=c(habitattarget,rep(0,(nrow(features)-1))),"spf"=10,"name"=features$name)
  
 #Write the spec.dat file in the /MarxanData/input/ directory
  write.table(spec.dat,file=paste0("MarxanData/input/spec_",parameters$figname,".dat"),sep=",",row.names = F,quote=F)

  #If running the case studies
  if(isTRUE(parameters$toimport)){
    
    #Create the assymbound.dat file (equivalent to bound.dat if using boundary length in Marxan) which gives the penalty between patches for missing strong connections. This is only used in the linkage-based method of spatial dependency
    assymbound<-generate_assymbound(habitatDispersal)
   
    #Write the assymbound.dat into directory /MarxanData/input/
    write.table(assymbound,file=paste0("MarxanData/input/assymbound_",parameters$figname,".dat"),sep=",",row.names = F,quote=F)
    
    #If running the simulated seascapes
  } else {
    
      #The flow matrix is the dispersal probability matrix row multiplied by the habitat area in each patch (this is the same if using a uniform area)
      flowDisp<-sweep(habitatDispersal, MARGIN=1, habitatMat[,"area"], `*`)
      
      #The migration matrix is the flow matrix divided by column sums
      migDisp<-t(t(flowDisp)/colSums(flowDisp))
      migDisp[which(is.na(migDisp))]<-0
      assymbound<-generate_assymbound(migDisp)
      
      #Write the assymbound.dat into directory /MarxanData/input/
      write.table(assymbound,file=paste0("MarxanData/input/assymbound_",parameters$figname,".dat"),sep=",",row.names = F,quote=F)
  }
  
  #Create a matrix where the targets for each of the connectivity metric conservation features is stored 
  relativetargets<-matrix(0,ncol=nrow(features),nrow=nrow(scenarios))
  
  #Index of the first metric
  firstmetric<-which(scenarios[,"scenario_name"]=="spatial")[length(which(scenarios[,"scenario_name"]=="spatial"))]+1
  
  #Row index of scenario dataframe which are using this metric
  allforthismetric<-which(scenarios[,"scenario_name"]==scenarios[,"scenario_name"][firstmetric])
  
  #Calculate the target for each of the metric quantiles based on the total target of the metric
  if(length(allforthismetric)>0){
  
  #Last index of the scenario row using this metric  
  lastmetric<-firstmetric+length(allmetrics)-1
  
  for(y in firstmetric:lastmetric ){
    
    #Inxed of scenario row matching metric
    allforthismetric<-which(scenarios[,"scenario_name"]==scenarios[,"scenario_name"][y])
    
    #Index of metric
    metricindexfortarget<-which(scenarios[,"scenario_name"][y]==gsub("_.*","",features$name))
    
    #Subset of the puvspr2.dat file
    rijsubset<-subset(rij,rij$species%in%c(metricindexfortarget))
    
    #Subset of the target matrix
    extractmat<-relativetargets[allforthismetric,metricindexfortarget]
    
    #How many planning units are in each quantile
    puthatcanbeineach<-as.numeric(table(rijsubset$species))
    
    #Look for any planning units that are not in any quantile
    anymissing<-outersect(metricindexfortarget,     as.numeric(names(table(rijsubset$species))))
    
    #This calculates how many planning units are in each quantile for a given total target of the metric
    if(length(anymissing)>0){
      for(m in 1:length(anymissing)){  
        position<-which(anymissing[m] == metricindexfortarget)
        if(position==1){
          puthatcanbeineach<-c(0,puthatcanbeineach)
        } else{
          puthatcanbeineach<-c(puthatcanbeineach[1:(position-1)],0,puthatcanbeineach[position:length(puthatcanbeineach)])
        }
      }
    }
    
      #Set the number of planning units required in each quantile to meet the target
      extractmat[,1]<-as.numeric(scenarios[,"target"][allforthismetric])*parameters$N
      for(h in 1:nrow(extractmat)){ 
        for(g in 1:ncol(extractmat)){ 
          if(extractmat[h,g]>puthatcanbeineach[g]){
            extractmat[h,g+1]<-extractmat[h,g]-puthatcanbeineach[g]
            extractmat[h,g]<-puthatcanbeineach[g]
          }
        }
      }
      
      #Divide the number of planning units by total possible number in the quantile so it becomes a percentage target
      for(k in 1:ncol(extractmat)){
        if(puthatcanbeineach[k]!=0){
          extractmat[,k]<-extractmat[,k]/puthatcanbeineach[k]
        }
      }
    }
    
    #Fill out relevant section of the matrix storing all of the targets for all of the quantiles 
    relativetargets[allforthismetric,metricindexfortarget]<-extractmat
    
  }
  
  #Set habitat area target
  relativetargets[,1]<-habitattarget
  
  return(relativetargets)
}

## This takes the original near neighbour probability matrix of dispersal and rewires some non-diagonal connections to create a small-world links pattern
rewirehabitatDispersal<-function(habitatDispersal){
  
  #New matrix without diagonals
  habitatDispersalnodiag<-habitatDispersal
  
  #Set diagonals to zero
  diag(habitatDispersalnodiag)<-0
  
  #Create graph object
  g1<-graph_from_adjacency_matrix(habitatDispersalnodiag,weighted=TRUE,mode="directed")
  
  #Rewire graph
  g2<-rewire(g1,each_edge(rewirePROB=0.1,loops=F,multiple=F))
  
  #Graph as matrix
  habitatDispersalrewire<-as.matrix(as_adjacency_matrix(g2, attr="weight"))
  
  #Restore original diagonal values
  diag(habitatDispersalrewire)<-diag(habitatDispersal)
  
  #Export the small-world links probability matrix
  habitatDispersal<-habitatDispersalrewire
  return(habitatDispersal)
}

## Runs Marxan for each scenario of connectivity method and habitat protection target
run_scenarios<-function(scenarios,habitatMat,habitatDispersal,parameters){
  
  #Set working directory
  setwd(oDir)
  
  #Create a matrix of targets for each quantile of each connectivity metric (in this case these are the conservation features added to the habitat target). Also generate the pu.dat, spec.dat, puvspr2.dat, and assymbound.dat file in the directory /MarxanData/input/
  relativetargets<-generate_marxan_inputfiles(scenarios,habitatMat,habitatDispersal,parameters)

  #Set working directory
  setwd(oDir)
        
        #Loop running Marxan for each scenario (row of the scenario dataframe)
        for(k in 1:nrow(scenarios)){
          
          #Set number of Marxan replicates 
          numreps<-parameters$numreps
          
          #If running either baseline or connectivity metrics methods, then turn off use of boundary file
          if(scenarios[,"csm"][k] %in% c(0,NA)){
            assymcon<-0
            FIXEDCSM<-0
          } else {
            
            #If running linkage-based spatial dependency methods, then turn on boundary file and take connectivity strength modifier (CSM) value from the scenario dataframe
            assymcon<-1 
            FIXEDCSM<-scenarios[,"csm"][k]
          }
          
          #Adjust the spec.dat file in the directory /MarxanData/input/ to have the correct targets for features 
          adjust_spec.dat(k,relativetargets,scenarios,parameters)  #this will change spec.dat
          
          #Set working directory to /MarxanData subdirectory
          setwd(sMxDir)
          
          #Load input.dat
          input.dat <- readLines("input.dat")
          
          #Set scenario name
          iSCENNAMEParam <- which(regexpr("SCENNAME",input.dat)==1)
          input.dat[iSCENNAMEParam] <- paste0("SCENNAME ",parameters$figname,"_",scenarios[,"name"][k])
          
          #Set spec.dat file
          iSPECNAMEParam <- which(regexpr("SPECNAME",input.dat)==1)
          input.dat[iSPECNAMEParam]<-paste0("SPECNAME ",scenarios[,"specdat"][k],".dat")
          
          #Set pu.dat file
          iPUNAMEParam <- which(regexpr("PUNAME",input.dat)==1)
          input.dat[iPUNAMEParam] <- paste0("PUNAME ",scenarios[,"pudat"][k],".dat")
          
          #Set puvspr2.dat file
          iPUVSPRNAMEParam <- which(regexpr("PUVSPRNAME",input.dat)==1)
          input.dat[iPUVSPRNAMEParam] <- paste0("PUVSPRNAME ",scenarios[,"puvspr2dat"][k],".dat")
          
          #Set assymbound.dat file
          iBOUNDNAMEParam <- which(regexpr("BOUNDNAME",input.dat)==1)
          input.dat[iBOUNDNAMEParam] <- paste0("BOUNDNAME ",scenarios[,"assymbounddat"][k],".dat")
          
          #Set connectivity strength modifier (CSM) value
          iBLMParam <- which(regexpr("BLM",input.dat)==1)
          input.dat[iBLMParam] <- paste0("BLM ",FIXEDCSM)
          
          #Set number of Marxan replicates
          iNUMREPSParam <- which(regexpr("NUMREPS",input.dat)==1)
          input.dat[iNUMREPSParam] <- paste0("NUMREPS ",numreps)
          
          #Turn on or off use of boundary file
          iASYMMETRICCONNECTIVITYParam <- which(regexpr("ASYMMETRICCONNECTIVITY",input.dat)==1)
          input.dat[iASYMMETRICCONNECTIVITYParam] <- paste0("ASYMMETRICCONNECTIVITY ",assymcon)
          
          #Save input.dat
          writeLines(input.dat,con="input.dat")
          
          #run Marxan
          system("marxan.exe -s input.dat")
        }

        #Save Marxan outputs in directories named after scenarios. Any existing directories are overwritten. 
        if (file.exists(paste0("output/",parameters$figname))){unlink(paste0("output/",parameters$figname),recursive=T)}   
        dir.create(paste0("output/",parameters$figname))
        
        #Save each scenario
        for(u in 1:nrow(scenarios)){
          
          #Grab the files matching the scenario name
          myfiles<-list.files(pattern=paste0(parameters$figname,"_",scenarios[u,"name"]),path="output/")
          
          #Move the files into the newly created directory named after the scenario
          file.copy(paste0("output/",myfiles),to=paste0("output/",parameters$figname,"/",myfiles))
          options(warn=-1)
          file.remove(paste0("output/",myfiles))
          options(warn=0)
        }
        
        #Set working directory to root directory
        setwd(oDir)
}

## Creates a dataframe of Marxan results, including id of solutions (0 not designated, 1 designated), cost of solutions, sum area protected, number of patches selected, equilibrium settlement 
results_df<-function(results,parameters,scenarios){
  
  #Take results dataframe and get column names
  resultsdf<-results
  colnames(resultsdf)<-colnames(results)
  
  #Set relevant columns to numeric
  class(resultsdf$SolutionCost)<-"numeric"
  class(resultsdf$SumCost)<-"numeric"
  class(resultsdf$SumArea)<-"numeric"
  class(resultsdf$NrReserves)<-"numeric"
  class(resultsdf$HitTarget)<-"numeric"
  
  
  #Type is whether the method of using connectivity is spatial (linkage-based) or connectivity metric (site-based)
  resultsdf$type<-0
  resultsdf$type[grep("spatial",resultsdf$Name)]<-"spatial"
  resultsdf$type[which(resultsdf$type==0)]<-"metric"
  
  #New empty column to save the name of the connectivity metric if 'type' is 'metric'
  resultsdf$detailtype<-0
  
  #Set the name of the connectivity metric for case studies
  if(isTRUE(parameters$toimport)){
  allnames<-strsplit(as.character(resultsdf$Name), "_")
  for(u in 1:length(allnames)){  
    workingname<-allnames[[u]]
    resultsdf$detailtype[u]<-workingname[grep(parameters$prioritisation,workingname)+1]
  }
  
  #Set the name of the connectivity metric for simulated seascapes
  } else {
      allnames<-strsplit(as.character(resultsdf$Name), "_")
      for(u in 1:length(allnames)){
        workingname<-allnames[[u]]
        resultsdf$detailtype[u]<-workingname[grep(parameters$prioritisation,workingname)+2]
      }  
  }
  
  #Set the % habitat target for the scenario 
  resultsdf$target<-scenarios[,"target"]

  #Set the connectivity strenth modifier (CSM) value of the scenario
  resultsdf$csm<-scenarios[,"csm"]
    
  return(resultsdf)
}

## Save the results dataframe 
save_resultsDF<-function(resultsDF,figname){
  setwd(oDir)
  write.table(resultsDF,paste0("results/resultsDF/",figname,".csv"),sep=",",row.names = F)  
}

## Opposite of intersect
outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

## Maximum of columns
colMax <- function(data) sapply(data, max, na.rm = TRUE)

## Generates the simulated seascape or imports habitat patches and dispersal matrices of case studies
generate_seascape<-function(parameters){

  #Generating simulated sescapes
  if(isFALSE(parameters$toimport)){
    
      #1. Generate seascape
      habitatMat<-habitat_spawn_loop(parameters) 
      
      #2. Add cost to seascape matrix    
      habitatMat<-add_cost(habitatMat,parameters)           
      
      #3. Calculate straight line distance between patches
      habitatDistance<-patch_distance(habitatMat)
      
      #4. Calculate dispersal as probability matrix
      habitatDispersal<-dispersal(habitatMat,parameters,localrecruitment=T) 
      
   #Importing case study dispersal matrices and x,y coordinates of habitat patches
  } else {
    
      #1. Import x,y of habitat patches
      habitatMat<-import_hmat(parameters)
      
      #2. Import probability matrix of larval dispersal
      habitatDispersal<-import_pmat(parameters)
      
      #3. Calculate straight line distance between patches
      habitatDistance<-patch_distance(habitatMat)
     
  }
  
  #5. Add network metrics to seascape matrix 
  habitatMat<-add_metrics(habitatMat,habitatDispersal,parameters)
  
  #6. Save items
  save_habitatMat(parameters,habitatMat,habitatDistance,habitatDispersal,parameters$figname)
  
 #List of habitat matrix, distance matrix, dispersal matrix
  l1<-list("habitatMat"=habitatMat,"habitatDistance"=habitatDistance,"habitatDispersal"=habitatDispersal)
  
  return(l1)
}




## Run the dispersal-per-recruit model from the package ConnMatTools to calculate equilibrium settlement of larvae as a measure of reserve performance.
run_DPR<-function(results,habitatMat,habitatDispersal,parameters,LEPvalue){
  
  #Calculate the equilibrium settlement for each row of the results dataframe
  for(i in 1:nrow(results)){
  
    #Vector of planning unit solutions where 0 is not designated, 1 is designated
    solution_1<-as.numeric(results[i,grep("s1",colnames(results))[1]:grep(paste0("s",parameters$N),colnames(results))])
    
    #ID of planning units which are designated
    solution_1INDEX<-which(solution_1==1)
    
    #ID of planning units which are not designated
    solution_0INDEX<-which(solution_1==0)
      
    #Get appropriate collapse slope
    slope <- settlerRecruitSlopeCorrection(habitatDispersal,use.arpack = F)
  
    #Number of patches
    n <- dim(habitatDispersal)[2]
    
    #Vector of lifetime egg production values
    LEP <- rep(LEPvalue,n)
    
    #Set LEP of reserves to 1  
    LEP[solution_1INDEX] <- 1
    
    #Run DPR model from ConnMatTools
    ret <- DispersalPerRecruitModel(LEP,habitatDispersal,recruits0=rep(Rmax,n),timesteps=1:200,slope=slope,Rmax=1,settler.recruit.func=BevertonHolt)
  
    #Dataframe of time steps, settlers, recruits, eggs
    popDF<-data.frame("year"=1:ncol(ret$settlers),"settlers"=colSums(ret$settlers),"recruits"=colSums(ret$recruits),"eggs"=colSums(ret$eggs))
    
    #Save the values of settlers and recruits across all planning units, in designated planning units, and in non-designated planning units
    results$DPRsettlers[i]<-round(as.numeric(popDF$settlers[timetorun]),4)
    results$DPRrecruits[i]<-round(as.numeric(popDF$recruits[timetorun]),4)
    results$DPRsettlersfished[i]<-round(as.numeric(sum(ret$settlers[solution_0INDEX,timetorun])),4)
    results$DPRsettlersnotfished[i]<-round(as.numeric(sum(ret$settlers[solution_1INDEX,timetorun])),4)
    results$DPRrecruitsfished[i]<-round(as.numeric(sum(ret$recruits[solution_0INDEX,timetorun])),4)
    results$DPRrecruitsnotfished[i]<-round(as.numeric(sum(ret$recruits[solution_1INDEX,timetorun])),4)
      
  }
  
  #Return dataframe of results
  return(results) 
}

## Import the habitat patch x and y coordinates for the case studies
import_hmat<-function(parameters){  
  
  #Set working directory
  setwd(oDir)
  
  #Read the imported planning unit table containing x,y, centroids and planning unit id
  putable<-read.table("imported habitat matrix.csv",stringsAsFactors = F,header=T,sep=",")
  
  #Subset to the species
  putable<-putable[grep(parameters$fishtype, putable$species),]

  #Build the hMat matrix containing planning unit id, x,y, centroids, uniform area, radius (not used), disp (dispersal not used), and cost
  hMat<-matrix(0,ncol=7,nrow=nrow(putable))
  colnames(hMat)<-c("puid","x","y","area","radius","disp","cost")
  hMat[,"puid"]<-putable$puid
  hMat[,"x"]<-putable$x_coord_epsg4326
  hMat[,"y"]<-putable$y_coord_epsg4326
  hMat[,"area"]<-10
  hMat[,"radius"]<-sqrt(hMat[,"area"]/pi)
  hMat[,"disp"]<-NA
  hMat[,"cost"]<-10
        
  return(hMat)
}

## Import the dispersal probability matrix for the case study species
import_pmat<-function(parameters){   
  
  #Read in the dispersal probability connectivity matrix
  cxmat<-read.table(paste0("imported probability matrix.csv"),stringsAsFactors = F,header=F,sep=",")
  
  #Subset to the appropriate species
  cxmat<-cxmat[grep(parameters$fishtype,cxmat$V1),]
  
  #Remove NA columns
  cxmat <- cxmat[,colSums(is.na(cxmat))<nrow(cxmat)]
  cxmat<-cxmat[,-1]
  
  #Object as matrix
  cxmat2<-as.matrix(cxmat)
  
  #Name the rows and columns after the sequential index of the habitat patch
  colnames(cxmat2)<-1:ncol(cxmat2)
  rownames(cxmat2)<-1:nrow(cxmat2)
  
  #Check that the dimensions of the square matrix match the number of habitat patches defined in 'parameters' 
  dim(cxmat2)==parameters$N
  
  return(cxmat2)
}

## Create figure 8 showing the histogram of local retention and box plot of equilibrium settlement in solutions
plot_figure_8<-function(){
  
#=#=#=#= Calculate the equilibrium settlement in each planning unit to create the box plots in figure 8, planning units are differentiated between those which are selected in a solution when using the connectivity-based features method with a 20% target for local retention and those that are not selected. Also grab the local retention value of each species for each planning unit and store in dataframe to create the histograms.
  
  #Blank dataframe in which to store the equilibrium settlement of each planning unit
  puidpopdfapp<-data.frame("puid"=c(),"pop"=c(),"localsolution20"=c(),"fish"=c())
  
  #Blank dataframe in which to store all local retention values
  allfishLRdfapp<-data.frame("species"=c(),"LR"=c())
  
  #List of all case studies
  allcasestudies<-c("Cuke","Trout","Mudcrab","rabbitfish")
  
  #For each species, grab the equilibrium settlement in each planning unit
  for(qw in 1:4){  
    
    #Load the results for the species with LEP = 0
    theseresults<-read.table(paste0("resultsDF/importLEP00_",allcasestudies[qw],"_marxan3_r0001_subset.csv"))
    
    #Column id which is "s1"
    colidsolFIRST<-which(colnames(theseresults)=="s1")
    
    #Column id which is the solution of the number of patches
    colidsolLAST<-which(colnames(theseresults)==paste0("s",parameters$N))
    
    #Grab the solutions for method using connectivity metrics with a 20% target for local retention
    solutionloc<-as.numeric(unname(theseresults[grep("localretention_t0.20",theseresults$Name),colidsolFIRST:colidsolLAST]))
    
    #Set up parameter values for case study
    parameters<-generate_parameters(N=100,mapsize=5000,disp=10,area="same",cost="same",discretize=10,toimport=T,prioritisation="marxan3",numreps=100,topN=10,fishtype=allcasestudies[qw],FIXEDCSM=100,runversion=1,addname="")
    
    #Generate habitatMat, dispersal, and distance matrices
    seascape<-generate_seascape(parameters)
    habitatMat<-seascape[["habitatMat"]]
    habitatDispersal<-seascape[["habitatDispersal"]]
    habitatDistance<-seascape[["habitatDistance"]]
    
    #Dataframe of local retention values
    allfishLRdf<-data.frame("species"=allcasestudies[qw],"LR"=diag(habitatDispersal))
    
    #Attach local retention values to dataframe of all species
    allfishLRdfapp<-rbind(allfishLRdfapp,allfishLRdf)
    
    #Re-run the DispersalPerRecruitModel from the package ConnMatTools
    #Get appropriate collapse slope
    slope <- settlerRecruitSlopeCorrection(habitatDispersal,use.arpack = F)
    
    #Vector of lifetime egg production values, here all planning units have the same value as we are assessing equilibrium settlement before implementation of reserves
    LEP <- 1
    
    #Run DPR model from ConnMatTools
    ret <- DispersalPerRecruitModel(LEP,habitatDispersal,recruits0=rep(Rmax,parameters$N),timesteps=1:200,slope=slope,Rmax=1,settler.recruit.func=BevertonHolt)
    
    #Dataframe to store equilibrium settlement in each planning unit
    puidpopdf<-data.frame("puid"=1:nrow(habitatDispersal),"pop"=ret$settlers[,200],"localsolution20"=solutionloc)
    
    #Subset of planning units which are selected as solutions
    LRsolspopCuke<-puidpopdf$pop[which(puidpopdf$localsolution20==1)]
    
    #Subset of planning units which are not selected as solutions
    nonLRsolspopCuke<-puidpopdf$pop[which(puidpopdf$localsolution20==0)]
    
    #T-test to see whether the two are different
    print(t.test(LRsolspopCuke,nonLRsolspopCuke,var.equal=F))
    
    #Create a new column to store the name of the species
    puidpopdf$fish<-allcasestudies[qw]
    
    #Append this species dataframe to the dataframe containing all species
    puidpopdfapp<-rbind(puidpopdfapp,puidpopdf)
  }
  
  #Turn the species column into a factor with levels
  puidpopdfapp$fish<-factor(puidpopdfapp$fish,levels=c("Trout","Cuke","rabbitfish","Mudcrab"))
  
  #Dataframe of labels to put into figure 8 box plots
  boxtextdata<-data.frame(labels=c("b","d","f","h"),"localsolution20"=0.5,pop=c(3,2,.2,.2),"fish"=factor(c("Trout","Cuke","rabbitfish","Mudcrab")))
  
  #Box plots on the right-hand side of figure 8
  pBOXPLOTS<-ggplot(puidpopdfapp,aes(x=as.factor(localsolution20),group=localsolution20,y=pop))+geom_boxplot()+theme_classic()+scale_x_discrete(labels=c("Non-reserves","Reserves"),name="")+ylab("Population size")+facet_wrap(~fish,ncol=1,scales="free_y")+theme(strip.text = element_blank(),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12))+geom_text(data=boxtextdata,aes(x=localsolution20,y=pop,label=labels),size=5)
  
  #Dataframe of labels to put into figure 8 histograms
  histtextdata<-data.frame(labels=c("c","a","g","e"),"species"=c("bCuke","aTrout","dMudcrab","crabbitfish"),"Count"=30,LR=0.05)
  
  #Re-label species to the correct order
  allfishLRdfapp$species<-recode(allfishLRdfapp$species,"Cuke"="bCuke",
                              "Trout"="aTrout",
                              "Mudcrab"="dMudcrab",
                              "rabbitfish"="crabbitfish")
  
  #Get the mean for creating a vertical line
  mu <- plyr::ddply(allfishLRdfapp, "species", summarise, grp.mean=mean(LR))
  
  #Histogram of local retention on the left-hand side of figure 8
  pHISTOGRAMS<-ggplot(allfishLRdfapp, aes(x=LR)) +geom_density()+  geom_vline(data=mu, aes(xintercept=grp.mean),linetype="dashed")+facet_wrap(~species,ncol=1)+geom_histogram(aes(y=..density..), alpha=0.5, position="identity",binwidth=0.025)+theme_classic()+xlab("Local retention")+ylab("Count")+theme(strip.text = element_blank(),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12))+geom_text(data=histtextdata,aes(x=LR,y=Count,label=labels),size=5)+scale_y_continuous(breaks=c(0,15,30),labels=c(0,"",30))+scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))
  
}

## Compare the spatial distribution of the solutions, creates maps in appendix
plot_appendix_maps<-function(){
  
  #% habitat target for which to create the map
  thistarget<-0.2

  #List of all case studies
  allcasestudies<-c("Cuke","Trout","Mudcrab","rabbitfish")
  
  #For each species, make maps
  for(qw in 1:4){  
  
    #Set up parameters for imported case study
    parameters<-generate_parameters(N=100,mapsize=5000,disp=10,area="same",cost="same",discretize=10,toimport=T,prioritisation="marxan3",numreps=100,topN=10,fishtype=allcasestudies[qw],FIXEDCSM=100,runversion=1,addname="")
  
    #Generate habitatMat, dispersal, and distance matrices
    seascape<-generate_seascape(parameters)
    habitatMat<-seascape[["habitatMat"]]
    habitatDispersal<-seascape[["habitatDispersal"]]
    habitatDistance<-seascape[["habitatDistance"]]
    
    #Read in the resultsDF dataframe of results so we can colour in patches by reserve designation (0 or 1)
    resultsDFread<-read.table(paste0("resultsDF/importLEP00_",allcasestudies[qw],"_marxan3_r0001_subset.csv"))
    
    #xmin, xmax, ymin, ymax to plot in maps for Southeast Sulawesi species
    if(allcasestudies[qw]%in%c("rabbitfish","Mudcrab")){
      mapextents<-c(119.8177,124.615,-7.80792,-2.486466) 
    }
    
    #xmin, xmax, ymin, ymax to plot in maps for Coral Triangle species
    if(allcasestudies[qw]%in%c("Trout","Cuke")){
      mapextents<-c(95.00933,168.8435,-12.28981,21.12061)
    }
    
    #Column id which is "s1"
    colidsolFIRST<-which(colnames(resultsDFread)=="s1")
    
    #Column id which is the solution of the number of patches
    colidsolLAST<-which(colnames(resultsDFread)==paste0("s",parameters$N))
    
    #Combine first and last column index
    colindex<-c(colidsolFIRST:colidsolLAST)
    
    #Set what % target to map
    rowindex<-which(resultsDFread$target==thistarget)
    
    #Create a matrix of only the solutions (0s and 1s for all planning units for all scenarios)
    summedsolnsmat<-resultsDFread[rowindex,colindex]
    
    #Empty vector
    colmeinbase<-c()
    
    #For each scenario get a numeric vector to color in by
    for(i in 1:nrow(summedsolnsmat)){
    
      #If designated, 1, if not, 0
      colmein<-as.numeric(summedsolnsmat[i,]==1)
      
      #Attach to empty vector 
      colmeinbase<-c(colmeinbase,colmein)
    }
    
  #Get a shapefile for countries
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  #Repeat the habitatMat for each of the method types (of which there are 10), used to facet_wrap the map
  habitatMatrepeat<-do.call("rbind", replicate(length(unique(resultsDFread$detailtype)), data.frame(habitatMat), simplify = FALSE))
  
  #Give correct name of method
  habitatMatrepeat$detailtype<-rep(resultsDFread$detailtype[rowindex],each=nrow(habitatMat))
  
  #Create a set of maps showing which planning units are designated and which are not  
  pMAPAPPENDIX<-ggplot()+geom_point(data=habitatMatrepeat,aes(x=x,y=y,fill=as.factor(colmeinbase),color=as.factor(colmeinbase)),pch=21,size=1)+geom_sf(data=world,fill="gray90")+xlim(mapextents[1],mapextents[2])+ylim(mapextents[3],mapextents[4])+theme_classic()+scale_color_manual(values=c("grey40","grey40"),name="",labels=c("Non-reserve","Reserve"))+scale_fill_manual(values=c("white","black"),name="",labels=c("Non-reserve","Reserve"))+ylab("Latitude")+xlab("Longitude")+facet_wrap(~detailtype,nrow=2)+theme(strip.background = element_blank(),legend.position = "none")+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf)+annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf)+annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  }
}

## This changes the names of the files so they are not overwritten. It adds a label for LEP (lifetime egg production) and the value of the LEP to the file.
givemenewfignames<-function(newtype,LEPval,casetype){
  
    #New file names if running sea cucumber
    if(newtype=="Cuke"){
      if(LEPval==0){giveout<-"importLEP00_Cuke_marxan3"  } else {
        if(LEPval==25){giveout<-"importLEP25_Cuke_marxan3"} else {
          if(LEPval==50){giveout<-"importLEP50_Cuke_marxan3"} else {
            if(LEPval==75){giveout<-"importLEP75_Cuke_marxan3"}
          }
        }
      }
    }
    
    #New file names if running coral trout
    if(newtype=="Trout"){
      if(LEPval==0){giveout<-"importLEP00_Trout_marxan3"  } else {
        if(LEPval==25){giveout<-"importLEP25_Trout_marxan3"} else {
          if(LEPval==50){giveout<-"importLEP50_Trout_marxan3"} else {
            if(LEPval==75){giveout<-"importLEP75_Trout_marxan3"}
          }
        }
      }
    }
    
    #New file names if running rabbitfish
    if(newtype=="rabbitfish"){
      if(LEPval==0){giveout<-"importLEP00_rabbitfish_marxan3"  } else {
        if(LEPval==25){giveout<-"importLEP25_rabbitfish_marxan3"} else {
          if(LEPval==50){giveout<-"importLEP50_rabbitfish_marxan3"} else {
            if(LEPval==75){giveout<-"importLEP75_rabbitfish_marxan3"} 
          }
        }
      }
    }
    
    #New file names if running mudcrab
    if(newtype=="Mudcrab"){
      if(LEPval==0){giveout<-"importLEP00_Mudcrab_marxan3"  } else {
        if(LEPval==25){giveout<-"importLEP25_Mudcrab_marxan3"} else {
          if(LEPval==50){giveout<-"importLEP50_Mudcrab_marxan3"} else {
            if(LEPval==75){giveout<-"importLEP75_Mudcrab_marxan3"}
          }
        }
      }
    }
  
  return(giveout)
  
}

## STAN Bayesian model to estimate the effects of connectivity method on equilibrium settlement at each % target
STANmodel<-function(){
  
  #Run Bayesian model for simulated seascapes with near neighbour pattern (simulatedNN) or small-world links (simulatedSWL) or for the case studies (case)
  for(casetype in c("simulatedNN","simulatedSWL","case")){
    
    #If running case studies
    if(casetype=="case"){
      
      #List of all species in case studies
      allcasestudies<-c("Cuke","Trout","Mudcrab","rabbitfish")
    }
    
    #If running simulated seascapes
    if(casetype%in%c("simulatedNN","simulatedSWL")){
      
      #All of the dispersal distances in the simulated seascapes
      allcasestudies<-c(50,100,150,200,250)
    }
    
    #Empty dataframe to store the results
    coefDFapp<-data.frame("type"=c(),"coef"=c(),"stderror"=c(),"target"=c(),"stanmedianfished"=c(),"stanmeanfished"=c(),"stanstderrorupperfished"=c(),"stanstderrorlowerfished"=c(),"stanmedianNOTfished"=c(),"stanmeanNOTfished"=c(),"stanstderrorupperNOTfished"=c(),"stanstderrorlowerNOTfished"=c(),"fish"=c(),"LEP"=c(),"stanmediansum"=c(),"stanmeansum"=c(),"stanstderroruppersum"=c(),"stanstderrorlowersum"=c(),"stanmedianfished.r"=c(),"stanmeanfished.r"=c(),"stanstderrorupperfished.r"=c(),"stanstderrorlowerfished.r"=c(),"stanmedianNOTfished.r"=c(),"stanmeanNOTfished.r"=c(),"stanstderrorupperNOTfished.r"=c(),"stanstderrorlowerNOTfished.r"=c(),"stanmediansum.r"=c(),"stanmeansum.r"=c(),"stanstderroruppersum.r"=c(),"stanstderrorlowersum.r"=c())
    
    #Loop over all of the species or dispersal distances
    for(ra2 in 1:length(allcasestudies)){
      
      #Loop over all of the lifetime egg production values (LEP)
      for(ra3 in c(0,25,50,75)){
    
        #If running the case studies
        if(casetype=="case"){
          
          #Grab the modified file name with LEP value and species name
          figname<-givemenewfignames(allcasestudies[ra2],LEPval=ra3,casetype)
            
            #Read in the results dataframe of the top solution (*r0001.csv)
            for(runnumber in 1){
              readresultsDF<-read.table(paste0("results/resultsDF/",figname,formatC(runnumber, width = 4, format = "d", flag = "0"),"_subset.csv"),sep=",",header=T)
              
              #Define a column with the Marxan run
              readresultsDF$run<-runnumber
            }
            
           #Read in the results dataframe of the second best to tenth best solutions (*r0002.csv - *r0010.csv)
            for(runnumber in 2:10){
              readresultsDF1<-read.table(paste0("results/resultsDF/",figname,formatC(runnumber, width = 4, format = "d", flag = "0"),"_subset.csv"),sep=",",header=T)
              
              #Define a column with the Marxan run number
              readresultsDF1$run<-runnumber
              
              #Append to the other results dataframes
              readresultsDF<-rbind(readresultsDF,readresultsDF1)
            }
        }
        
        #If running the simulated seascapes
        if(casetype%in%c("simulatedNN","simulatedSWL")){
          
          #Format the csm value
          lpv<-formatC(ra3, width = 2, format = "d", flag = "0")
          
          #Files to read in if running the near neighbour pattern
          if(casetype=="simulatedNN"){
            allfilesmain<-list.files(path="results/resultsDF/",pattern=paste0("N100oQLEP",lpv))
          }
          
          #Files to read in if running the small world links pattern
          if(casetype=="simulatedSWL"){
            allfilesmain<-list.files(path="results/resultsDF/",pattern=paste0("N100oQrw0.1LEP",lpv))
          }
          
          #Subset to the files which have "_subset" in them as we only want the best csm
          allfilesmain<-allfilesmain[grep("_subset.csv",allfilesmain)]
          
          #Subset to the files matching the dispersal
          allfilesmain<-allfilesmain[grep(paste0("disp",formatC(allcasestudies[ra2], width = 4, format = "d", flag = "0")),allfilesmain)]
          
          #Sort the files alphabetically
          allfilesmain<-sort(allfilesmain)
          
          #Read in the results dataframe of the first results file
          for(readrow in 1){   
            readresultsDF<-read.table(paste0("results/resultsDF/",allfilesmain[readrow]),sep=",",header=T)
            
            #Grab the id of the seascape replicate (1-100), this is different from the Marxan solution (1-10)
            rv<-sub(".*_rv","",allfilesmain[readrow])
            rv<-sub("_r.*","",rv)
            
            #Create a new column and give it the value of the seascape replicate
            readresultsDF$rv<-rv
          }
          
          #Read in the results dataframe of the remaining results files
          for(readrow in 2:length(allfilesmain)){
            readresultsDF1<-read.table(paste0("results/resultsDF/",allfilesmain[readrow]),sep=",",header=T)
            
            #Grab the id of the seascape replicate (1-100), this is different from the Marxan solution (1-10)
            rv<-sub(".*_rv","",allfilesmain[readrow])
            rv<-sub("_r.*","",rv)
            
            #Create a new column and give it the value of the seascape replicate
            readresultsDF1$rv<-rv
            
            #Append to the other results dataframes
            readresultsDF<-rbind(readresultsDF,readresultsDF1)
          }
        }
        
     #New column in the results dataframe (readresultsDF) that is the sum of settler and recruit settlement in designated and not designated planning units
     readresultsDF$DPRsettlerssum<-readresultsDF$DPRsettlersfished+readresultsDF$DPRsettlersnotfished
     readresultsDF$DPRrecruitssum<-readresultsDF$DPRrecruitsfished+readresultsDF$DPRrecruitsnotfished
     
    #Vector of all habitat target percentages 
    alltargets<-unique(readresultsDF$target)
    
    #Empty dataframe in which to store model coefficients
    coefDF<-data.frame("type"=rep(c("abaseline","betweenness","eigenvector","indegree","influx","localretention","outdegree","outflux","page","spatial"),length(alltargets)),"coef"=0,"stderror"=0,"target"=rep(alltargets,each=10),"stanmedianfished"=0,"stanmeanfished"=0,"stanstderrorupperfished"=0,"stanstderrorlowerfished"=0,"stanmedianNOTfished"=0,"stanmeanNOTfished"=0,"stanstderrorupperNOTfished"=0,"stanstderrorlowerNOTfished"=0,"stanmediansum"=0,"stanmeansum"=0,"stanstderroruppersum"=0,"stanstderrorlowersum"=0,"stanmedianfished.r"=0,"stanmeanfished.r"=0,"stanstderrorupperfished.r"=0,"stanstderrorlowerfished.r"=0,"stanmedianNOTfished.r"=0,"stanmeanNOTfished.r"=0,"stanstderrorupperNOTfished.r"=0,"stanstderrorlowerNOTfished.r"=0,"stanmediansum.r"=0,"stanmeansum.r"=0,"stanstderroruppersum.r"=0,"stanstderrorlowersum.r"=0)
    
      #Loop for each of the habitat targets, run the Bayesian model to obtain model coefficients
      for(targ in 1:length(alltargets)){
        
        #Index of first row to read
        dfindex<-targ+9*(targ-1)
        
        #Index of last row to read
        dfindexend<-dfindex+9       
       
        #MODELS FOR SETTLERS
        #STAN model with settlers in non-designated planning units as a function of connectivity method used (detailtype)
        stfished<-stan_glm( DPRsettlersfished ~ detailtype, data = subset(readresultsDF,readresultsDF$target==alltargets[targ]))
        
        #STAN model with settlers in designated planning units as a function of connectivity method used (detailtype)
        stNOTfished<-stan_glm( DPRsettlersnotfished ~ detailtype, data = subset(readresultsDF,readresultsDF$target==alltargets[targ]))
        
        #STAN model with settlers in both designated and non-designated planning units as a function of connectivity method used (detailtype)
        stsum<-stan_glm( DPRsettlerssum ~ detailtype, data = subset(readresultsDF,readresultsDF$target==alltargets[targ]))
      
        #Get parameters of the model for non-designated planning units
        posteriorsfished <- insight::get_parameters(stfished)
        
        #Empty vector to store upper interval
        aupperfished<-c()
        
        #Empty vector to store lower interval
        alowerfished<-c()
        
        #Calculate 89% credible interval using highest density interval
        for(i in 1:ncol(posteriorsfished)){   
          a1fished<-hdi(posteriorsfished[,i],ci=0.89)
          aupperfished[i]<-a1fished$CI_high
          alowerfished[i]<-a1fished$CI_low
        }
        
        #Get parameters of the model for designated planning units
        posteriorsNOTfished <- insight::get_parameters(stNOTfished)
        
        #Empty vector to store upper interval
        aupperNOTfished<-c()
        
        #Empty vector to store lower interval
        alowerNOTfished<-c()
        
        #Calculate 89% credible interval using highest density interval
        for(i in 1:ncol(posteriorsNOTfished)){   
          a1NOTfished<-hdi(posteriorsNOTfished[,i],ci=0.89)
          aupperNOTfished[i]<-a1NOTfished$CI_high
          alowerNOTfished[i]<-a1NOTfished$CI_low
        }
        
        #Get parameters of the model for both designated and non-designated planning units
        posteriorssum <- insight::get_parameters(stsum)
        
        #Empty vector to store upper interval
        auppersum<-c()
        
        #Empty vector to store lower interval
        alowersum<-c()
        
        #Calculate 89% credible interval using highest density interval
        for(i in 1:ncol(posteriorssum)){   
          a1sum<-hdi(posteriorssum[,i],ci=0.89)
          auppersum[i]<-a1sum$CI_high
          alowersum[i]<-a1sum$CI_low
        }
        
        #Save the model parameters in the dataframe for non-designated planning units
        coefDF$stanmedianfished[dfindex:dfindexend]<-  c(apply(posteriorsfished,2,median))
        coefDF$stanmeanfished[dfindex:dfindexend]<-  c(apply(posteriorsfished,2,mean))
        coefDF$stanstderrorupperfished[dfindex:dfindexend]<-c(aupperfished-apply(posteriorsfished,2,median))
        coefDF$stanstderrorlowerfished[dfindex:dfindexend]<-c(apply(posteriorsfished,2,median)-alowerfished)
        
        #Save the model parameters in the dataframe for designated planning units
        coefDF$stanmedianNOTfished[dfindex:dfindexend]<-  c(apply(posteriorsNOTfished,2,median))
        coefDF$stanmeanNOTfished[dfindex:dfindexend]<-  c(apply(posteriorsNOTfished,2,mean))
        coefDF$stanstderrorupperNOTfished[dfindex:dfindexend]<-c(aupperNOTfished-apply(posteriorsNOTfished,2,median))
        coefDF$stanstderrorlowerNOTfished[dfindex:dfindexend]<-c(apply(posteriorsNOTfished,2,median)-alowerNOTfished)
        
        #Save the model parameters in the dataframe for both designated and non-designated planning units
        coefDF$stanmediansum[dfindex:dfindexend]<-  c(apply(posteriorssum,2,median))
        coefDF$stanmeansum[dfindex:dfindexend]<-  c(apply(posteriorssum,2,mean))
        coefDF$stanstderroruppersum[dfindex:dfindexend]<-c(auppersum-apply(posteriorssum,2,median))
        coefDF$stanstderrorlowersum[dfindex:dfindexend]<-c(apply(posteriorssum,2,median)-alowersum)
        
        #MODELS FOR RECRUITS
        #STAN model with recruits in non-designated planning units as a function of connectivity method used (detailtype)
        stfished.r<-stan_glm(DPRrecruitsfished ~ detailtype, data = subset(readresultsDF,readresultsDF$target==alltargets[targ]))
        
        #STAN model with recruits in designated planning units as a function of connectivity method used (detailtype)
        stNOTfished.r<-stan_glm(DPRrecruitsnotfished ~ detailtype, data = subset(readresultsDF,readresultsDF$target==alltargets[targ]))
        
        #STAN model with recruits in both designated and non-designated planning units as a function of connectivity method used (detailtype)
        stsum.r<-stan_glm(DPRrecruitssum ~ detailtype, data = subset(readresultsDF,readresultsDF$target==alltargets[targ]))
        
        #Get parameters of the model for non-designated planning units
        posteriorsfished.r <- insight::get_parameters(stfished.r)
        
        #Empty vector to store upper interval
        aupperfished.r<-c()
        
        #Empty vector to store lower interval
        alowerfished.r<-c()
        
        #Calculate 89% credible interval using highest density interval
        for(i in 1:ncol(posteriorsfished.r)){   
          a1fished.r<-hdi(posteriorsfished.r[,i],ci=0.89)
          aupperfished.r[i]<-a1fished.r$CI_high
          alowerfished.r[i]<-a1fished.r$CI_low
        }
        
        #Get parameters of the model for designated planning units
        posteriorsNOTfished.r <- insight::get_parameters(stNOTfished.r)
        
        #Empty vector to store upper interval
        aupperNOTfished.r<-c()
        
        #Empty vector to store lower interval
        alowerNOTfished.r<-c()
        
        #Calculate 89% credible interval using highest density interval
        for(i in 1:ncol(posteriorsNOTfished.r)){   
          a1NOTfished.r<-hdi(posteriorsNOTfished.r[,i],ci=0.89)
          aupperNOTfished.r[i]<-a1NOTfished.r$CI_high
          alowerNOTfished.r[i]<-a1NOTfished.r$CI_low
        }
        
        #Get parameters of the model for both designated and non-designated planning units
        posteriorssum.r <- insight::get_parameters(stsum.r)
        
        #Empty vector to store upper interval
        auppersum.r<-c()
        
        #Empty vector to store lower interval
        alowersum.r<-c()
        
        #Calculate 89% credible interval using highest density interval
        for(i in 1:ncol(posteriorssum.r)){   
          a1sum.r<-hdi(posteriorssum.r[,i],ci=0.89)
          auppersum.r[i]<-a1sum.r$CI_high
          alowersum.r[i]<-a1sum.r$CI_low
        }
        
        #Save the model parameters in the dataframe for non-designated planning units
        coefDF$stanmedianfished.r[dfindex:dfindexend]<-  c(apply(posteriorsfished.r,2,median))
        coefDF$stanmeanfished.r[dfindex:dfindexend]<-  c(apply(posteriorsfished.r,2,mean))
        coefDF$stanstderrorupperfished.r[dfindex:dfindexend]<-c(aupperfished.r-apply(posteriorsfished.r,2,median))
        coefDF$stanstderrorlowerfished.r[dfindex:dfindexend]<-c(apply(posteriorsfished.r,2,median)-alowerfished.r)
        
        #Save the model parameters in the dataframe for designated planning units
        coefDF$stanmedianNOTfished.r[dfindex:dfindexend]<-  c(apply(posteriorsNOTfished.r,2,median))
        coefDF$stanmeanNOTfished.r[dfindex:dfindexend]<-  c(apply(posteriorsNOTfished.r,2,mean))
        coefDF$stanstderrorupperNOTfished.r[dfindex:dfindexend]<-c(aupperNOTfished.r-apply(posteriorsNOTfished.r,2,median))
        coefDF$stanstderrorlowerNOTfished.r[dfindex:dfindexend]<-c(apply(posteriorsNOTfished.r,2,median)-alowerNOTfished.r)
        
        #Save the model parameters in the dataframe for both designated and non-designated planning units
        coefDF$stanmediansum.r[dfindex:dfindexend]<-  c(apply(posteriorssum.r,2,median))
        coefDF$stanmeansum.r[dfindex:dfindexend]<-  c(apply(posteriorssum.r,2,mean))
        coefDF$stanstderroruppersum.r[dfindex:dfindexend]<-c(auppersum.r-apply(posteriorssum.r,2,median))
        coefDF$stanstderrorlowersum.r[dfindex:dfindexend]<-c(apply(posteriorssum.r,2,median)-alowersum.r)
      }
    
      #Create a column and assign species name if case study or dispersal distance if simulated seascape
      coefDF$fish<-allcasestudies[ra2]
      
      #Create a column and assign LEP value
      coefDF$LEP<-ra3
      
      #Append results dataframe
      coefDFapp<-rbind(coefDFapp,coefDF)
    }
    }
    
    #Save the results dataframe
    write.table(coefDFapp,file=paste0("results/stan/STAN results ",casetype,".csv"),sep=",",row.names = F)
  } 
}

## Create figures 4, 5, and 7. Bar charts showing the equilibrium settlement in reserve networks designed using different connectivity methods at 20% target. Settlement is differentiated between inside reserves and outside. The cumulative length of each bar is the sum of both inside and outside. 
plot_figures_4_5_7<-function(){
 
  #Create bar charts for simulated seascapes with near neighbour pattern (simulatedNN) or small-world links (simulatedSWL) or for the case studies (case)
  for(casetype in c("simulatedNN","simulatedSWL","case")){
    
    #Read in the results file containing the model coefficients 
    coefDFapp<-read.table(paste0("results/stan/STAN results ",casetype,".csv"),header=T,sep=",")
    
    #Rescale so that the intercept of 1 is the baseline, a value of 2 is twice of 1
    coefDFapp<-coefDFapprescale(coefDFapp)
     
    #Subset of results to only 20% target
    coefDFappsubsub<-subset(coefDFappsub,coefDFappsub$target==0.2)
    
    #Re-label the connectivity methods
    coefDFappsubsub$type<-dplyr::recode(coefDFappsubsub$type,
                                           "abaseline"="Baseline",
                                           "spatial" = "Spatial dependency",
                                           "page" = "Google PageRank",
                                           "betweenness" = "Betweenness centrality",
                                           "eigenvector" = "Eigenvector centrality",
                                           "outflux" = "Out-flow",
                                           "outdegree" = "Out-degree",
                                           "localretention" = "Local retention",
                                           "indegree"="In-degree",
                                           "influx"="In-flow")
    
    #Change the connectivity methods column to a factor with defined levels
    coefDFappsubsub$type<-factor(coefDFappsubsub$type,levels=c("Baseline","Betweenness centrality","Eigenvector centrality","Google PageRank","Local retention","In-degree","In-flow","Out-degree","Out-flow","Spatial dependency"))
    #Subset of the results to required columns
    amm2<-coefDFappsubsub[,c("stanmedianfishedSCL2","stanmedianNOTfishedSCL2","type","fish","LEP")]
    
    #Results to long format
    ammDFwide<-amm2 %>% gather(stanmedianfishedSCL2,stanmedianNOTfishedSCL2,-c(type,fish,LEP))
    
    #Rename 4th and 5th column of long dataframe
    colnames(ammDFwide)[4:5]<-c("where","effect")
    
    #Add in the upper and lower limits of model coefficients
    ammDFwide$upper<-c(coefDFappsubsub$stanstderrorupperfishedSCL2,coefDFappsubsub$stanstderrorupperNOTfishedSCL2)
    ammDFwide$upper<-ammDFwide$effect+ammDFwide$upper
    ammDFwide$lower<-c(coefDFappsubsub$stanstderrorlowerfishedSCL2,coefDFappsubsub$stanstderrorlowerNOTfishedSCL2)
    ammDFwide$lower<-ammDFwide$effect-ammDFwide$lower
    
    #Combine the lower limits in designated and non-designated so they align in the bar chart
    ammDFwide$lower[ammDFwide$where == "stanmedianfishedSCL2"] <- with(ammDFwide,lower[where == "stanmedianfishedSCL2"] +lower[where == "stanmedianNOTfishedSCL2"])
    
    #Combine the upper limits in designated and non-designated so they align in the bar chart
    ammDFwide$upper[ammDFwide$where == "stanmedianfishedSCL2"] <- with(ammDFwide,upper[where == "stanmedianfishedSCL2"] +upper[where == "stanmedianNOTfishedSCL2"])
    
    #Assign labels of methods from the results dataframe
    ammDFwide$labels<-ammDFwide$type
    
    #Re-label whether equilibrium settlement is in designated or non-designated planning units
    ammDFwide$where<-recode(ammDFwide$where,"stanmedianNOTfishedSCL2"="Reserve",
                            "stanmedianfishedSCL2"="Non-reserve")
    
    #Column as character
    ammDFwide$type<-as.character(ammDFwide$type)
    
    #Want to remove some axis labels later on when using facet_wrap
    ammDFwide$type[which(ammDFwide$LEP%in%c(25,50,75))]<-paste0(  ammDFwide$type[which(ammDFwide$LEP%in%c(25,50,75))],"no_display")
    
    #For simulated seascapes
    if(casetype%in%c("simulatedNN","simulatedSWL")){
    
      #For near neighbour pattern
      if(casetype=="simulatedNN"){ 
         
         #Grab the dispersal distances 50, 100, 150, 150, 200, 250
         ammDFwide<-ammDFwide[which(ammDFwide$fish %in%c("50","100","150","200","250")),]
         
         #Set dispersal distance as factor with levels
         ammDFwide$fish<-factor(ammDFwide$fish,levels=rev(c("50","100","150","200","250")))
      }
      
      #For small-world links pattern
      if(casetype=="simulatedSWL"){
        
        #Grab the dispersal distances 50, 100, 150, 150, 200, 250
        ammDFwide<-ammDFwide[which(ammDFwide$fish %in%c("50rw0.1","100rw0.1","150rw0.1","200rw0.1","250rw0.1")),]
        
        #Set dispersal distance as factor with levels
        ammDFwide$fish<-factor(ammDFwide$fish,levels=rev(c("50rw0.1","100rw0.1","150rw0.1","200rw0.1","250rw0.1")))
        }
    
      #Remove out-flow and out-degree as these are identical to in-flow and in-degree  
      ammDFwide<-ammDFwide[-grep("Out",ammDFwide$type),]
      
      #Re-label connectivity methods, no_display will be removed in facet_wrap
      ammDFwide$type[which(ammDFwide$type=="In-flow")]<-"Flow"
      ammDFwide$type[which(ammDFwide$type=="In-degree")]<-"Degree"
      ammDFwide$type[which(ammDFwide$type=="In-flowno_display")]<-"Flowno_display"
      ammDFwide$type[which(ammDFwide$type=="In-degreeno_display")]<-"Degreeno_display"
      
      #Set the connectivity method as a factor with levels
      ammDFwide$type<-factor(ammDFwide$type,levels=c("Baseline","Betweenness centrality","Eigenvector centrality","Google PageRank","Local retention","Degree","Flow","Spatial dependency","Baselineno_display","Betweenness centralityno_display","Eigenvector centralityno_display","Google PageRankno_display","Local retentionno_display","Degreeno_display","Flowno_display","Spatial dependencyno_display"))
    }
    
    #If a case study
    if(casetype=="case"){
      
      #Set the species as a factor with levels
      ammDFwide$fish<-factor(ammDFwide$fish,levels=c("Trout","Cuke","rabbitfish","Mudcrab"))
      
      #Set the connectivity method as a factor with levels
      ammDFwide$type<-factor(ammDFwide$type,levels=c("Baseline","Betweenness centrality","Eigenvector centrality","Google PageRank","Local retention","In-degree","Out-degree","In-flow","Out-flow","Spatial dependency","Baselineno_display","Betweenness centralityno_display","Eigenvector centralityno_display","Google PageRankno_display","Local retentionno_display","In-degreeno_display","Out-degreeno_display","In-flowno_display","Out-flowno_display","Spatial dependencyno_display"))
    }
    
  #Bar chart of equilibrium settlement in reserve networks created using each of the connectivity methods, split into inside and outside reserves
  pBAR<-ggplot(ammDFwide, aes(fill=where, x=effect, y=type)) + geom_bar( stat="identity")+theme_classic()+facet_wrap(~fish+LEP,scales="free",ncol=4)+geom_errorbar(aes(xmin=lower,xmax=upper),width=0.2)+scale_y_discrete(label = delete_no_display,name="")+scale_fill_manual(name="",values=c("gray78","gray38"))+theme(strip.text = element_blank(),legend.position="none",strip.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_text(size=10),axis.text.x=element_text(size=12),axis.title.x=element_text(size=12))+xlab("Relative equilibrium settlement")
  
    #Function for manual scale 
    scale_override <- function(which, scale) {
      if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
        stop("which must be an integer of length 1")
      }
      
      if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
        stop("scale must be an x or y position scale")
      }
      
      structure(list(which = which, scale = scale), class = "scale_override")
    }
    
    #Function for manual facetwrap
    CustomFacetWrap <- ggproto(
      "CustomFacetWrap", FacetWrap,
      init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
        # make the initial x, y scales list
        scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
        
        if(is.null(params$scale_overrides)) return(scales)
        
        max_scale_x <- length(scales$x)
        max_scale_y <- length(scales$y)
        
        # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
        for(scale_override in params$scale_overrides) {
          which <- scale_override$which
          scale <- scale_override$scale
          
          if("x" %in% scale$aesthetics) {
            if(!is.null(scales$x)) {
              if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
              scales$x[[which]] <- scale$clone()
            }
          } else if("y" %in% scale$aesthetics) {
            if(!is.null(scales$y)) {
              if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
              scales$y[[which]] <- scale$clone()
            }
          } else {
            stop("Invalid scale")
          }
        }
        
        # return scales
        scales
      }
    )
    
    #Function for custom facet wrap
    facet_wrap_custom <- function(..., scale_overrides = NULL) {
      # take advantage of the sanitizing that happens in facet_wrap
      facet_super <- facet_wrap(...)
      
      # sanitize scale overrides
      if(inherits(scale_overrides, "scale_override")) {
        scale_overrides <- list(scale_overrides)
      } else if(!is.list(scale_overrides) || 
                !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
        stop("scale_overrides must be a scale_override object or a list of scale_override objects")
      }
      
      facet_super$params$scale_overrides <- scale_overrides
      
      ggproto(NULL, CustomFacetWrap,
              shrink = facet_super$shrink,
              params = facet_super$params
      )
    }
    
    #Function for custom breaks
    count <- 0
    breaks_fun <- function(x) {
      count <<- count + 1L
      switch(
        count,
        c(0.025, 0.125),
        c(0.05, 0.054),
        c(0.05, 0.150),
        c(0.138, 0.144)
      )
    }
    
    #Manual breaks for the near neighbour bar chart
    if(casetype=="simulatedNN"){
      (pBARmanual<-pBAR +
         facet_wrap_custom(~fish+LEP, scales = "free", ncol = 4, scale_overrides = list(
           scale_override(17, scale_x_continuous(breaks = c(0,5,10, 15),labels=c(0,"",10,""))),
           scale_override(18, scale_x_continuous(breaks = c(0,2,4,6),labels=c(0,"",4,""))),
           scale_override(19, scale_x_continuous(breaks = c(0,.5,1,1.5,2),labels=c(0,"",1,"",2))),
           scale_override(20, scale_x_continuous(breaks = c(0,.5,1),labels=c(0,"",1))),
           scale_override(13, scale_x_continuous(breaks = c(0,5,10,15),labels=c(0,"",10,""))),
           scale_override(14, scale_x_continuous(breaks = c(0,2,4,6),labels=c(0,"",4,""))),
           scale_override(15, scale_x_continuous(breaks = c(0, .5,1,1.5,2),labels=c(0,"",1,"",2))),
           scale_override(16, scale_x_continuous(breaks = c(0, .5,1),labels=c(0,"",1))),
           scale_override(9, scale_x_continuous(breaks = c(0,5,10,15),labels=c(0,"",10,""))),
           scale_override(10, scale_x_continuous(breaks = c(0,2,4,6),labels=c(0,"",4,""))),
           scale_override(11, scale_x_continuous(breaks = c(0, .5,1,1.5,2),labels=c(0,"",1,"",2))),
           scale_override(12, scale_x_continuous(breaks = c(0,.5,1),labels=c(0,"",1))),
           scale_override(5, scale_x_continuous(breaks = c(0,10,20),labels=c(0,"",20))),
           scale_override(6, scale_x_continuous(breaks = c(0,2,4),labels=c(0,"",4))),
           scale_override(7, scale_x_continuous(breaks = c(0, .5,1,1.5),labels=c(0,"",1,""))),
           scale_override(8, scale_x_continuous(breaks = c(0, .25,.5,.75,1,1.25),labels=c(0,"",.5,"",1,""))),
           scale_override(1, scale_x_continuous(breaks = c(0, 10,20,30,40),labels=c(0,"",20,"",40))),
           scale_override(2, scale_x_continuous(breaks = c(0, 2,4),labels=c(0,"",4))),
           scale_override(3, scale_x_continuous(breaks = c(0, .5,1,1.5),labels=c(0,"",1,""))),
           scale_override(4, scale_x_continuous(breaks = c(0, .5,1),labels=c(0,"",1)))
         )))
    }
    
    #Manual breaks for the small-world links bar chart
    if(casetype=="simulatedSWL"){
      (pBARmanual<-pBAR +
         facet_wrap_custom(~fish+LEP, scales = "free", ncol = 4, scale_overrides = list(
           scale_override(17, scale_x_continuous(breaks = c(0,5,10),labels=c(0,"",10))),
           scale_override(18, scale_x_continuous(breaks = c(0,2,4,6),labels=c(0,"",4,""))),
           scale_override(19, scale_x_continuous(breaks = c(0,.5,1,1.5,2,2.5),labels=c(0,"",1,"",2,""))),
           scale_override(20, scale_x_continuous(breaks = c(0,.5,1),labels=c(0,"",1))),
           scale_override(13, scale_x_continuous(breaks = c(0,4,8,12),labels=c(0,"",8,""))),
           scale_override(14, scale_x_continuous(breaks = c(0,2,4,6),labels=c(0,"",4,""))),
           scale_override(15, scale_x_continuous(breaks = c(0, .5,1,1.5,2),labels=c(0,"",1,"",2))),
           scale_override(16, scale_x_continuous(breaks = c(0, .5,1),labels=c(0,"",1))),
           scale_override(9, scale_x_continuous(breaks = c(0,5,10,15),labels=c(0,"",10,""))),
           scale_override(10, scale_x_continuous(breaks = c(0,2,4),labels=c(0,"",4))),
           scale_override(11, scale_x_continuous(breaks = c(0, .5,1,1.5),labels=c(0,"",1,""))),
           scale_override(12, scale_x_continuous(breaks = c(0,.25,.5,.75,1,1.25),labels=c(0,"",.5,"",1,""))),
           scale_override(5, scale_x_continuous(breaks = c(0,5,10,15),labels=c(0,"",10,""))),
           scale_override(6, scale_x_continuous(breaks = c(0,1,2,3,4),labels=c(0,"",2,"",4))),
           scale_override(7, scale_x_continuous(breaks = c(0, .5,1,1.5),labels=c(0,"",1,""))),
           scale_override(8, scale_x_continuous(breaks = c(0, .25,.5,.75,1),labels=c(0,"",.5,"",1))),
           scale_override(1, scale_x_continuous(breaks = c(0, 10,20,30),labels=c(0,"",20,""))),
           scale_override(2, scale_x_continuous(breaks = c(0, 1,2,3,4),labels=c(0,"",2,"",4))),
           scale_override(3, scale_x_continuous(breaks = c(0, .5,1,1.5),labels=c(0,"",1,""))),
           scale_override(4, scale_x_continuous(breaks = c(0, .5,1),labels=c(0,"",1)))
         )))
    }
    
    #Manual breaks for the case studies bar chart
    if(casetype=="case"){
      (pBARmanual<-pBAR +
         facet_wrap_custom(~fish+LEP, scales = "free", ncol = 4, scale_overrides = list(
           scale_override(5, scale_x_continuous(breaks = c(0,1,2,3,4),labels=c(0,"",2,"",4))),
           scale_override(6, scale_x_continuous(breaks = c(0,1,2),labels=c(0,"",2))),
           scale_override(7, scale_x_continuous(breaks = c(0,.5,1,1.5),labels=c(0,"",1,""))),
           scale_override(8, scale_x_continuous(breaks = c(0,.25,.5,.75,1,1.25),labels=c(0,"",.5,"",1,""))),
           scale_override(1, scale_x_continuous(breaks = c(0,5,10),labels=c(0,"",10))),
           scale_override(2, scale_x_continuous(breaks = c(0,1,2,3,4),labels=c(0,"",2,"",4))),
           scale_override(3, scale_x_continuous(breaks = c(0, .5,1,1.5),labels=c(0,"",1,""))),
           scale_override(4, scale_x_continuous(breaks = c(0, .5,1),labels=c(0,"",1))),
           scale_override(13, scale_x_continuous(breaks = c(0,1,2,3,4),labels=c(0,"",2,"",4))),
           scale_override(14, scale_x_continuous(breaks =c(0,1,2,3,4),labels=c(0,"",2,"",4))),
           scale_override(15, scale_x_continuous(breaks = c(0, 1,2),labels=c(0,"",2))),
           scale_override(16, scale_x_continuous(breaks = c(0,.5,1,1.5),labels=c(0,"",1,""))),
           scale_override(9, scale_x_continuous(breaks = c(0,1,2,3,4),labels=c(0,"",2,"",4))),
           scale_override(10, scale_x_continuous(breaks = c(0,1,2,3),labels=c(0,"",2,""))),
           scale_override(11, scale_x_continuous(breaks = c(0, .5,1,1.5),labels=c(0,"",1,""))),
           scale_override(12, scale_x_continuous(breaks = c(0,.5,1),labels=c(0,"",1)))
         )))
    }
  }
}

## Function to suppress labels in plot
delete_no_display <- function(v) {
  if_else(str_detect(v, 'no_display'), '', v)
}

## Rescales the equilibrium settlement so that a value of 1 is the baseline. A value of 2 is twice as high as 1
coefDFapprescale<-function(coefDFapp){
  
  #Duplicate dataframe storing the Bayesian model coefficients
  coefDFapp2<-coefDFapp
  
  #Unique row id to dataframe
  coefDFapp2$rowid<-1:nrow(coefDFapp2)
  
  #Dataframe to store rescaled coefficients
  DFappend<-data.frame("stanmedianfishedSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorupperfishedSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorlowerfishedSCL"=rep(NA,nrow(coefDFapp2)),"stanmedianNOTfishedSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorupperNOTfishedSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorlowerNOTfishedSCL"=rep(NA,nrow(coefDFapp2)), "stanmedianfishedSCL2"=rep(NA,nrow(coefDFapp2)),"stanstderrorupperfishedSCL2"=rep(NA,nrow(coefDFapp2)),"stanstderrorlowerfishedSCL2"=rep(NA,nrow(coefDFapp2)),"stanmedianNOTfishedSCL2"=rep(NA,nrow(coefDFapp2)),"stanstderrorupperNOTfishedSCL2"=rep(NA,nrow(coefDFapp2)),"stanstderrorlowerNOTfishedSCL2"=rep(NA,nrow(coefDFapp2)),"stanmediansumSCL"=rep(NA,nrow(coefDFapp2)),"stanstderroruppersumSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorlowersumSCL"=rep(NA,nrow(coefDFapp2)),"stanmedianfished.rSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorupperfished.rSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorlowerfished.rSCL"=rep(NA,nrow(coefDFapp2)),"stanmedianNOTfished.rSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorupperNOTfished.rSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorlowerNOTfished.rSCL"=rep(NA,nrow(coefDFapp2)),"stanmediansum.rSCL"=rep(NA,nrow(coefDFapp2)),"stanstderroruppersum.rSCL"=rep(NA,nrow(coefDFapp2)),"stanstderrorlowersum.rSCL"=rep(NA,nrow(coefDFapp2)))
  coefDFapp2<-cbind(coefDFapp2,DFappend)
  
  #Vector of habitat targets
  alltargets<-unique(coefDFapp2$target)
  
  #Vector of lifetime egg production values
  allLEPS<-unique(coefDFapp2$LEP)
  
  #Vector of unique species or dispersal distances
  allfish<-unique(coefDFapp2$fish)
  
  #Loop over all targets
  for(i in 1:length(alltargets)){   
    
    #Loop over all lifetime egg production values
    for(k in 1:length(allLEPS)){  
      
      #Loop over all species or dispersal distances
      for(j in 1:length(allfish)){   
        
        #Subset to the model coefficients for a given target, LEP, and species or dispersal distance
        coefDFapp2sub<-coefDFapp2[which(coefDFapp2$target==alltargets[i] & coefDFapp2$LEP==allLEPS[k] & coefDFapp2$fish==allfish[j]),]
        
        #SETTLERS
        
        #Add model coefficient of connectivity method to the intercept, then divide by intercept. Settlement in non-designated planning units.
        coefDFapp2sub$stanmedianfishedSCL<-(coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])/coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")]
        
        #Add model coefficient of connectivity method to the intercept, then divide by intercept. Settlement in designated planning units.
        coefDFapp2sub$stanmedianNOTfishedSCL<-(coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])/coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")]
        
        #Add model coefficient of connectivity method to the intercept, then divide by intercept. Settlement in both designated and non-designated planning units.
        coefDFapp2sub$stanmediansumSCL<-(coefDFapp2sub$stanmediansum+coefDFapp2sub$stanmediansum[which(coefDFapp2sub$type=="abaseline")])/coefDFapp2sub$stanmediansum[which(coefDFapp2sub$type=="abaseline")]
        
        #Scale the upper 89% confidence intervals for non-designated planning units
        coefDFapp2sub$stanstderrorupperfishedSCL<-coefDFapp2sub$stanmedianfishedSCL*((coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderrorupperfished)/(coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmedianfishedSCL
        
        #Scale the upper 89% confidence intervals for designated planning units
        coefDFapp2sub$stanstderrorupperNOTfishedSCL<-coefDFapp2sub$stanmedianNOTfishedSCL*((coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderrorupperNOTfished)/(coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmedianNOTfishedSCL
        
        #Scale the upper 89% confidence intervals for both designated and non-designated planning units
        coefDFapp2sub$stanstderroruppersumSCL<-coefDFapp2sub$stanmediansumSCL*((coefDFapp2sub$stanmediansum+coefDFapp2sub$stanmediansum[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderroruppersum)/(coefDFapp2sub$stanmediansum+coefDFapp2sub$stanmediansum[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmediansumSCL
        
        #Scale the lower 89% confidence intervals for non-designated planning units
        coefDFapp2sub$stanstderrorlowerfishedSCL<-coefDFapp2sub$stanmedianfishedSCL-coefDFapp2sub$stanmedianfishedSCL*((coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanstderrorlowerfished)/(coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])
        
        #Scale the lower 89% confidence intervals for designated planning units
        coefDFapp2sub$stanstderrorlowerNOTfishedSCL<-coefDFapp2sub$stanmedianNOTfishedSCL-coefDFapp2sub$stanmedianNOTfishedSCL*((coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanstderrorlowerNOTfished)/(coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])
        
        #Scale the lower 89% confidence intervals for both designated and non-designated planning units
        coefDFapp2sub$stanstderrorlowersumSCL<-coefDFapp2sub$stanmediansumSCL-coefDFapp2sub$stanmediansumSCL*((coefDFapp2sub$stanmediansum+coefDFapp2sub$stanmediansum[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanstderrorlowersum)/(coefDFapp2sub$stanmediansum+coefDFapp2sub$stanmediansum[which(coefDFapp2sub$type=="abaseline")])
        
        #RECRUITS
        
        #Add model coefficient of connectivity method to the intercept, then divide by intercept. Recruitment in non-designated planning units.
        coefDFapp2sub$stanmedianfished.rSCL<-(coefDFapp2sub$stanmedianfished.r+coefDFapp2sub$stanmedianfished.r[which(coefDFapp2sub$type=="abaseline")])/coefDFapp2sub$stanmedianfished.r[which(coefDFapp2sub$type=="abaseline")]
        
        #Add model coefficient of connectivity method to the intercept, then divide by intercept. Recruitment in designated planning units.
        coefDFapp2sub$stanmedianNOTfished.rSCL<-(coefDFapp2sub$stanmedianNOTfished.r+coefDFapp2sub$stanmedianNOTfished.r[which(coefDFapp2sub$type=="abaseline")])/coefDFapp2sub$stanmedianNOTfished.r[which(coefDFapp2sub$type=="abaseline")]

        #Add model coefficient of connectivity method to the intercept, then divide by intercept. Recruitment in both designated and non-designated planning units.
        coefDFapp2sub$stanmediansum.rSCL<-(coefDFapp2sub$stanmediansum.r+coefDFapp2sub$stanmediansum.r[which(coefDFapp2sub$type=="abaseline")])/coefDFapp2sub$stanmediansum.r[which(coefDFapp2sub$type=="abaseline")]
        
        #Scale the upper 89% confidence intervals for non-designated planning units
        coefDFapp2sub$stanstderrorupperfished.rSCL<-coefDFapp2sub$stanmedianfished.rSCL*((coefDFapp2sub$stanmedianfished.r+coefDFapp2sub$stanmedianfished.r[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderrorupperfished.r)/(coefDFapp2sub$stanmedianfished.r+coefDFapp2sub$stanmedianfished.r[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmedianfished.rSCL
        
        #Scale the upper 89% confidence intervals for designated planning units
        coefDFapp2sub$stanstderrorupperNOTfished.rSCL<-coefDFapp2sub$stanmedianNOTfished.rSCL*((coefDFapp2sub$stanmedianNOTfished.r+coefDFapp2sub$stanmedianNOTfished.r[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderrorupperNOTfished.r)/(coefDFapp2sub$stanmedianNOTfished.r+coefDFapp2sub$stanmedianNOTfished.r[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmedianNOTfished.rSCL
        
        #Scale the upper 89% confidence intervals for both designated and non-designated planning units
        coefDFapp2sub$stanstderroruppersum.rSCL<-coefDFapp2sub$stanmediansum.rSCL*((coefDFapp2sub$stanmediansum.r+coefDFapp2sub$stanmediansum.r[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderroruppersum.r)/(coefDFapp2sub$stanmediansum.r+coefDFapp2sub$stanmediansum.r[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmediansum.rSCL
        
        #Scale the lower 89% confidence intervals for non-designated planning units
        coefDFapp2sub$stanstderrorlowerfished.rSCL<-coefDFapp2sub$stanmedianfished.rSCL-coefDFapp2sub$stanmedianfished.rSCL*((coefDFapp2sub$stanmedianfished.r+coefDFapp2sub$stanmedianfished.r[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanstderrorlowerfished.r)/(coefDFapp2sub$stanmedianfished.r+coefDFapp2sub$stanmedianfished.r[which(coefDFapp2sub$type=="abaseline")])
        
        #Scale the lower 89% confidence intervals for designated planning units
        coefDFapp2sub$stanstderrorlowerNOTfished.rSCL<-coefDFapp2sub$stanmedianNOTfished.rSCL-coefDFapp2sub$stanmedianNOTfished.rSCL*((coefDFapp2sub$stanmedianNOTfished.r+coefDFapp2sub$stanmedianNOTfished.r[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanstderrorlowerNOTfished.r)/(coefDFapp2sub$stanmedianNOTfished.r+coefDFapp2sub$stanmedianNOTfished.r[which(coefDFapp2sub$type=="abaseline")])
        
        #Scale the lower 89% confidence intervals for both designated and non-designated planning units
        coefDFapp2sub$stanstderrorlowersum.rSCL<-coefDFapp2sub$stanmediansum.rSCL-coefDFapp2sub$stanmediansum.rSCL*((coefDFapp2sub$stanmediansum.r+coefDFapp2sub$stanmediansum.r[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanstderrorlowersum.r)/(coefDFapp2sub$stanmediansum.r+coefDFapp2sub$stanmediansum.r[which(coefDFapp2sub$type=="abaseline")])
        
        
        #Divide model coefficient in designated and non-designated by the coefficient in both. For baseline 
        for(trf in 1){
          coefDFapp2sub$stanmedianfishedSCL2[trf]<-coefDFapp2sub$stanmedianfished[trf]/coefDFapp2sub$stanmediansum[trf]
          coefDFapp2sub$stanmedianNOTfishedSCL2[trf]<-coefDFapp2sub$stanmedianNOTfished[trf]/coefDFapp2sub$stanmediansum[trf]

        }
        
        #Divide model coefficient in designated and non-designated by the coefficient in both. For connectivity runs
        for(trf in 2:nrow(coefDFapp2sub)){
          coefDFapp2sub$stanmedianfishedSCL2[trf]<-(coefDFapp2sub$stanmedianfished[trf]+coefDFapp2sub$stanmedianfished[1])/(coefDFapp2sub$stanmediansum[trf]+coefDFapp2sub$stanmediansum[1])*coefDFapp2sub$stanmediansumSCL[trf]
          coefDFapp2sub$stanmedianNOTfishedSCL2[trf]<-(coefDFapp2sub$stanmedianNOTfished[trf]+coefDFapp2sub$stanmedianNOTfished[1])/(coefDFapp2sub$stanmediansum[trf]+coefDFapp2sub$stanmediansum[1])*coefDFapp2sub$stanmediansumSCL[trf]
        }
        
        #Adjust the upper confidence interval in non-designated planning units to the new baseline
        coefDFapp2sub$stanstderrorupperfishedSCL2<-coefDFapp2sub$stanmedianfishedSCL2*((coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderrorupperfished)/(coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmedianfishedSCL2

        #Adjust the upper confidence interval in designated planning units to the new baseline
        coefDFapp2sub$stanstderrorupperNOTfishedSCL2<-coefDFapp2sub$stanmedianNOTfishedSCL2*((coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderrorupperNOTfished)/(coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmedianNOTfishedSCL2
        
        #Adjust the lower confidence interval in non-designated planning units to the new baseline
        coefDFapp2sub$stanstderrorlowerfishedSCL2<-coefDFapp2sub$stanmedianfishedSCL2*((coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderrorlowerfished)/(coefDFapp2sub$stanmedianfished+coefDFapp2sub$stanmedianfished[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmedianfishedSCL2

        #Adjust the lower confidence interval in designated planning units to the new baseline
        coefDFapp2sub$stanstderrorlowerNOTfishedSCL2<-coefDFapp2sub$stanmedianNOTfishedSCL2*((coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])+coefDFapp2sub$stanstderrorlowerNOTfished)/(coefDFapp2sub$stanmedianNOTfished+coefDFapp2sub$stanmedianNOTfished[which(coefDFapp2sub$type=="abaseline")])-coefDFapp2sub$stanmedianNOTfishedSCL2
     
        #Place the adjusted values in the correct location of the master results file
        coefDFapp2[coefDFapp2sub$rowid,]<-coefDFapp2sub
      }
    }
  }
  
  #Return results dataframe with adjusted values
  return(coefDFapp2)
}

## Plot figures 3 and 6 showing the connectivity methods which achieved the highest equilibrium settlement.
plot_topmethodssummary<-function(){
  
  #Create plots for simulated seascapes with near neighbour pattern (simulatedNN) or small-world links (simulatedSWL) or for the case studies (case)
  for(casetype in c("simulated","case")){
    
    #If case study
    if(casetype=="case"){
      #Read in the results file containing the model coefficients 
      coefDFapp<-read.table(paste0("results/stan/STAN results ",casetype,".csv"),header=T,sep=",")
      
      #Rescale so that the intercept of 1 is the baseline, a value of 2 is twice of 1
      coefDFapp<-coefDFapprescale(coefDFapp)
    }
    
    #If simulated seascapes
    if(casetype=="simulated"){
      
      #Read in the near neighbour results  
      coefDFapp<-read.table(paste0("results/stan/STAN results simulatedNN.csv"),header=T,sep=",")
      
      #Rescale so that the intercept of 1 is the baseline, a value of 2 is twice of 1
      coefDFapp<-coefDFapprescale(coefDFapp)
      
      #Read in the small-world links results  
      coefDFappSWL<-read.table(paste0("results/stan/STAN results simulatedSWL.csv"),header=T,sep=",")
      
      #Add "rw0.1" to the dispersal distance columns so we can differentiate between near neighbour and small-world links 
      coefDFappSWL$fish<-paste0(coefDFappSWL$fish,"rw0.1")
      
      #Rescale so that the intercept of 1 is the baseline, a value of 2 is twice of 1
      coefDFappSWL<-coefDFapprescale(coefDFappSWL)
      
      #Combine the results of near neighbour and small-world links pattern
      coefDFapp<-rbind(coefDFapp,coefDFappoQrw0.1)
    }
    
    #Remove baselines
    coefDFappsub<-coefDFappsub[-which(coefDFappsub$type=="abaseline"),]
    
    #Re-label connectivity methods 
    coefDFappsub$type<-dplyr::recode(coefDFappsub$type,
                                     "spatial" = "Spatial dependency",
                                     "page" = "Google PageRank",
                                     "betweenness" = "Betweenness centrality",
                                     "eigenvector" = "Eigenvector centrality",
                                     "outflux" = "Out-flow",
                                     "outdegree" = "Out-degree",
                                     "localretention" = "Local retention",
                                     "indegree"="In-degree",
                                     "influx"="In-flow")
    
    #Set connectivity methods as factor with levels
    coefDFappsub$type<-factor(coefDFappsub$type,levels=c("Betweenness centrality","Eigenvector centrality","Google PageRank","Local retention","In-degree","In-flow","Out-degree","Out-flow","Spatial dependency"))
    
    #If case study
    if(casetype=="case"){ 
      
      #Species as factor with levels
      coefDFappsub$fish<-factor(coefDFappsub$fish,levels=c("Cuke","Trout","Mudcrab","rabbitfish"))
    }
    
    #If simulated seascapes
    if(casetype=="simulated"){
      
      #Remove in-flow
      coefDFappsub<-coefDFappsub[-which(coefDFappsub$type=="In-flow"),]
      
      #Remove in-degree
      coefDFappsub<-coefDFappsub[-which(coefDFappsub$type=="In-degree"),]

      #Dispersal distances as factor with levels
      coefDFappsub$fish<-factor(coefDFappsub$fish,levels=c("50","100","150","200","250","50rw0.1","100rw0.1","150rw0.1","200rw0.1","250rw0.1"))
    }
    
    #Get the model coefficients for settlement in both designated and non-designated planning units 
    coefDFappsub$stanmedian<-coefDFappsub$stanmediansumSCL
    coefDFappsub$stanstderrorlower<-coefDFappsub$stanstderrorlowersumSCL
    coefDFappsub$stanstderrorupper<-coefDFappsub$stanstderroruppersumSCL
    
    #Empty dataframe to store the names of the best performing methods
    targDF<-data.frame("target"=c(),"LEP"=c(),"fish"=c(),"maxtype"=c())
    
    #Vector of unique species or dispersal distances
    allfish<-unique(coefDFappsub$fish)
    
    #Vector of unique lifetime egg production values
    allLEP<-unique(coefDFappsub$LEP)
    
    #Vector of habitat protection targets
    alltargets<-unique(coefDFappsub$target)
    
    #Loop over each species or dispersal distance
    for(g in 1:length(allfish)){
      
      #Loop over each lifetime egg production value
      for(h in 1:length(allLEP)){
        
        #Loop over each habitat target
        for(f in 1:length(alltargets)){
          
          #Subset to the results of a given species or dispersal distance, lifetime egg produciton, and habitat target
          coefDFappsubMMT<-coefDFappsub[which(coefDFappsub$fish==allfish[g] & coefDFappsub$LEP==allLEP[h] & coefDFappsub$target==alltargets[f]),]
          
          #Get row index of the method having the highest coefficient
          rowid<-which(coefDFappsubMMT$stanmediansumSCL==max(coefDFappsubMMT$stanmediansumSCL))
   
          #Find all methods whose confidence interval overlaps with the confidence interval of the maximum
          alltops<-as.character(coefDFappsubMMT$type[which((coefDFappsubMMT$stanmediansumSCL+coefDFappsubMMT$stanstderroruppersumSCL)>(coefDFappsubMMT$stanmediansumSCL[rowid]-coefDFappsubMMT$stanstderrorlowersumSCL[rowid]))])
          
          #Fill in the dataframe giving the name of the top methods
          targDFapp<-data.frame("target"=rep(c(alltargets[f]),length(alltops)),"LEP"=rep(c(allLEP[h]),length(alltops)),"fish"=rep(c(allfish[g]),length(alltops)),"maxtype"=alltops)
  
          #Append to the cumulative dataframe giving which method or methods were best for a given lifetime egg production value and target
          targDF<-rbind(targDF,targDFapp)
        }
      }
    }
    
    #Create the plot if using case studies
    if(casetype=="case"){
      
      #New column giving location of case study as Coral Triangle
      targDF$casestudy<-"CTI"
      
      #Label location for species in Southeast Sulawesi
      targDF$casestudy[which(targDF$fish%in%c("Mudcrab","rabbitfish"))]<-"SeSU"
      
      #Dataframe to the best performing methods
      finDF<-data.frame("target"=c(),"LEP"=c(),"maxtype"=c(),"howmany"=c(),"casestudy"=c(),"tothisfish"=c())
      
      #Vector of case studies
      allcases<-c("CTI","SeSU")
    
      #Loop over each lifetime egg production value
      for(h in 1:length(allLEP)){
        
        #Loop over each habitat target
        for(f in 1:length(alltargets)){    
          
            #Loop over each location of the case study
            for(x in 1:length(allcases)){
              
              #Subset of the results for a given lifetime egg production, target, and case study
              targDFRRRR<-targDF[which(targDF$LEP==allLEP[h] & targDF$target==alltargets[f] & targDF$casestudy==allcases[x]),]
              
              #Count the occurrence of each method as top
              d1<-data.frame(table(targDFRRRR$maxtype))
              
              #Order the count table
              d1<-d1[order(d1$Var1),]
              
              #Fill in the results dataframe with how many times each method was top
              finDFapp<-data.frame("target"=c(alltargets[f]),"LEP"=c(allLEP[h]),"maxtype"=d1$Var1,"howmany"=d1$Freq,"casestudy"=allcases[x],"tothisfish"="none")
              
              #In which of the species was the method top
              for(hwmany in which(finDFapp$howmany==1)){
                finDFapp$tothisfish[hwmany]<-as.character(targDFRRRR$fish[which(targDFRRRR$maxtype==finDFapp$maxtype[hwmany])])
              }
      
              #Append results to master dataframe
              finDF<-rbind(finDF,finDFapp)
            }
        }
      }
    
      #Shorter labels for connectivity methods
      finDF$maxtype<-dplyr::recode(finDF$maxtype,"Spatial dependency"="SpaDep",
                            "Google PageRank"="GooPag",
                            "Local retention"="LocRet",
                            "Eigenvector centrality"="EigCen",
                            "Out-flow"="OutFlo",
                            "In-degree"="In-Deg",
                            "In-flow"="In-Flo",
                            "Betweenness centrality"="BetCen")
    
      #Rescale how many times each method occurred
      finDF$howmany2<-scales::rescale(finDF$howmany,to=c(1,1.02))
    
      #Column giving unique row di
      finDF$uniqrow<-1:nrow(finDF)
      
      #Location to place LEP label
      finDF$LEPplacement<--999
     
      #How far apart method labels are to be
      stepsize<-5
      
      #Species as codes
      finDF$tothisfishcoded<-recode(finDF$tothisfish,"Trout"=1,
                                    "none"=2,
                                    "Cuke"=3,
                                    "rabbitfish"=1,
                                    "Mudcrab"=3)
      
      #Loop for each case study
      for(x in 1:2){
        
        #Loop for each lifetime egg production value
        for(h in 1:length(allLEP)){
          
          #Loop for each habitat protection target
          for(f in 1:length(alltargets)){
            finDFsub<-finDF[which(finDF$LEP==allLEP[h] & finDF$target==alltargets[f] & finDF$casestudy==allcases[x]),]
    
            #If there is only one top method then the placement is the same as the lifetime egg production value
            if(nrow(finDFsub)==1){
              finDF$LEPplacement[which(finDF$uniqrow==finDFsub$uniqrow)]<-finDFsub$LEP
            } else {
    
               #Re-order by code of species
               finDFsub<-finDFsub[order(finDFsub$tothisfishcoded),]
    
               #If there is an even number of top methods
               if(nrow(finDFsub)%%2==0){
                 
                 #Get number of rows in results
                 nrows<-nrow(finDFsub)
                 
                 #Placement above the mean
                 upper<-finDFsub$LEP[1]+stepsize*nrows/2-stepsize/2
                 
                 #Placement below the mean
                 lower<-finDFsub$LEP[1]-stepsize*nrows/2+stepsize/2
    
                 #Create a vector of placement
                 finDFsub$LEPplacement<-seq(from=lower,to=upper,by=stepsize)
                 
                 #Assign appropriate placement to each method
                 for(i in 1:nrow(finDFsub)){ 
                   finDF$LEPplacement[which(finDF$uniqrow==finDFsub$uniqrow[i])]<-finDFsub$LEPplacement[i]
                 }
    
                 #If there is an odd number of top methods
               } else {
    
                 #Get number of rows in results
                 nrows<-nrow(finDFsub)
                 
                 #Placement above the mean
                 upper<-finDFsub$LEP[1]+stepsize*nrows/2-stepsize/2
                 
                 #Placement below the mean
                 lower<-finDFsub$LEP[1]-stepsize*nrows/2+stepsize/2
    
                 #Create a vector of placement
                 finDFsub$LEPplacement<-seq(from=lower,to=upper,by=stepsize)
                 
                 #Assign appropriate placement to each method
                 for(i in 1:nrow(finDFsub)){ #  i=1
                   finDF$LEPplacement[which(finDF$uniqrow==finDFsub$uniqrow[i])]<-finDFsub$LEPplacement[i]
                 }
               }
             }
           }
         }
       }
    
     
      #Subset to relevant columns
      finDF<-finDF[,c("target","LEP","tothisfish","casestudy","maxtype")]
      
      #Long dataframe of top methods
      finDFlong<-data.frame("target"=c(),"LEP"=c(),"forwho1"=c(),"casestudy1"=c(),"maxtype"=c(),"targetplace"=c(),"placeLEP"=c())
      
      #Re-label the name of the top method
      finDF$maxtype<-recode(finDF$maxtype,
                                   "EigCen"="aEigCen",
                                   "GooPag"= "bGooPag",
                                   "LocRet"="cLocRet",
                                   "In-Deg"="dDeg",
                                   "In-Flo"="eFlow",
                                   "SpaDep"="fSpaDep")
      
      #Determines the spacing between methods on x-axis
      stepsize<-0.007
      
      #Determines spacing on y-axis
      steplep<-6
      
      #Order the dataframe by name of method
      finDF<-finDF[order(finDF$maxtype),]
      
      #Loop over each row of the dataframe
      for(i in 1:nrow(finDF)){
        
        #Subset of the results for a given target, lifetime egg production, and case study
        finDFlongsub<-finDF[which(finDF$target==finDF$target[i]&finDF$LEP==finDF$LEP[i]&finDF$casestudy==finDF$casestudy[i]),]
        
        #Vector of unique top methods
        allthese<-unique(finDFlongsub$maxtype)
        
        #Vector of species having top methods
        forwhos<-finDF$tothisfish[i]
        
        #Assign the name of the species based on the name of the case study
        if(forwhos=="none" & finDFlongsub$casestudy=="CTI"){
          forwhos<-c("Cuke","Trout")
        } else {
          if(forwhos=="none" & finDFlongsub$casestudy=="SeSU"){
            forwhos<-c("Mudcrab","rabbitfish")
          }
        }
        
        #If there is only one top method, the placement on x-axis is the target amount
        if(length(allthese)==1){
          targetplace<-finDF$target[i]
        } else {
          
          #If there is an even number of top methods, then place on x-axis using spacing
          if(length(allthese)%%2==0){
            
            #Number of rows of methods
            nrows<-length(allthese)
            
            #Upper limit of placement on x-axis
            upper<-finDF$target[i]+stepsize*nrows/2-stepsize/2

            #Lower limit of placement on x-axis
            lower<-finDF$target[i]-stepsize*nrows/2+stepsize/2
            
            #Sequence of placement on x-axis
            targetplace<-seq(from=lower,to=upper,by=stepsize)
            
            #If there is an odd number of top methods, then place on x-axis using spacing
          } else {
            
            #Number of rows of methods
            nrows<-length(allthese)
            
            #Upper limit of placement on x-axis
            upper<-finDF$target[i]+stepsize*nrows/2-stepsize/2
            
            #Lower limit of placement on x-axis
            lower<-finDF$target[i]-stepsize*nrows/2+stepsize/2
            
            #Sequence of placement on x-axis
            targetplace<-seq(from=lower,to=upper,by=stepsize)
          }   
        }
        
        #Assign placement on x-axis by method type
        targetplace<-targetplace[which(finDF$maxtype[i]==allthese)]
        
        #Vector of placement
        theseplaces<-seq(from=finDFlongsub$target[1]-stepsize*3,to=finDFlongsub$target[1]+stepsize*3,length.out=6)
        
        #This ensures that the methods are always in vertical alignment with themselves
        if(finDF$maxtype[i]=="aBetCen"){
          targetplace<-theseplaces[1]
        } else {
          if(finDF$maxtype[i]=="aEigCen"){
            targetplace<-theseplaces[1]
          } else {
            if(finDF$maxtype[i]=="bGooPag"){
              targetplace<-theseplaces[2]
            } else {
              if(finDF$maxtype[i]=="cLocRet"){
                targetplace<-theseplaces[3]
              } else {
                if(finDF$maxtype[i]=="dDeg"){
                  targetplace<-theseplaces[4]
                } else {
                  if(finDF$maxtype[i]=="eFlow"){
                    targetplace<-theseplaces[5]
                  } else {
                    if(finDF$maxtype[i]=="fSpaDep"){
                      targetplace<-theseplaces[6]
                    }
                  }
                }
              }
            }
          }
        }
        
        #Dataframe storing the placement of the methods
        finDFlongapp<-data.frame("target"=rep(finDF$target[i],length(forwhos)),"LEP"=rep(finDF$LEP[i],length(forwhos)),"forwho"=forwhos,"casestudy1"=rep(finDF$casestudy[i],length(forwhos)),"maxtype"=rep(finDF$maxtype[i],length(forwhos)),"targetplace"=targetplace,"placeLEP"=NA)
        
        #Loop over each row of the results
        for(qrt in 1:nrow(finDFlongapp)){
          
          #Upper limit of lifetime egg production placement
          lepup<-finDFlongapp$LEP[1]+steplep

          #Lower limit of lifetime egg production placement
          lepdown<-finDFlongapp$LEP[1]-steplep
          
          #Sequence of placement for lifetime egg production 
          thisseq<-seq(from=lepdown,to=lepup,length.out=2)
          
          #This ensures that the species are always in alignment with themselves
          if(finDFlongapp$forwho[qrt]=="Cuke"){finDFlongapp$placeLEP[qrt]<-thisseq[1]}
          if(finDFlongapp$forwho[qrt]=="Trout"){finDFlongapp$placeLEP[qrt]<-thisseq[2]}
          if(finDFlongapp$forwho[qrt]=="Mudcrab"){finDFlongapp$placeLEP[qrt]<-thisseq[1]}
          if(finDFlongapp$forwho[qrt]=="rabbitfish"){finDFlongapp$placeLEP[qrt]<-thisseq[2]}
        }
        
        #Attach the results to the master file
        finDFlong<-rbind(finDFlong,finDFlongapp)
      }
      
      #Re-label the connectivity methods
      finDFlong$maxtype<-factor(finDFlong$maxtype,levels=c("aEigCen","bGooPag","cLocRet","dDeg","eFlow","fSpaDep"))
      
      #Plot showing which connectivity methods performed best at a given habitat target (x-axis) and value for lifetime egg production (column) 
      pTOPMETHOD<-ggplot()+geom_point(data=finDFlong,aes(x=targetplace,y=placeLEP,shape=maxtype,color=maxtype),size=2)+geom_hline(yintercept=seq(from=0,to=75,by=25),lty=2,color="gray88")+geom_vline(xintercept=seq(from=0.075,to=.35,by=.05),lty=1,color="gray58")+geom_hline(yintercept=seq(from=12.5,to=87.5,by=25),lty=1,color="gray58")+ scale_x_continuous(expand = c(0, 0),limits=c(0.025,0.325),breaks=c(0.05,0.1,0.15,0.2,0.25,0.3),labels=seq(from=5,to=30,by=5)) + scale_y_continuous(expand = c(0, 0),limits=c(-12.5,87.5),breaks=c(0,25,50,75))+theme_classic()+xlab("Habitat protection (%)")+ylab("LEP of non-reserves (%)")+scale_color_manual(values=gg_color_hue(6),labels=c("Eigenvector centrality","Google PageRank","Local retention","In-degree","In-Flow","Spatial dependency"),name="Method")+theme(legend.position="bottom", strip.text = element_blank(),strip.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text=element_text(size=10))+annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf)+annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf)+facet_wrap(~casestudy1)+annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+scale_shape_manual(values=c(0,15,1,17,2,19),labels=c("Eigenvector centrality","Google PageRank","Local retention","In-degree","In-Flow","Spatial dependency"),name="Method")+geom_line(data=finDFlong,aes(x=targetplace,y=placeLEP,color=maxtype,group=interaction(maxtype, target,LEP)),size=1)
    }
   
    #Create the plot if using simulated seascapes 
    if(casetype=="simulated"){
      
      #Label as nearest neighbour (oQ)
      targDFs$casestudy<-"oQ"
      
      #Label the small-world links runs (oQrw0.1)
      targDFs$casestudy[grep("rw0.1",targDFs$fish)]<-"oQrw0.1"
      
      #Vector of both dispersal patterns
      allcases<-c("oQ","oQrw0.1")
      
      #Empty dataframe to store the top methods
      finDFs<-data.frame("target"=c(),"LEP"=c(),"maxtype"=c(),"forwho"=c(),casestudy=c())
      
      #Loop over each dispersal pattern
      for(x in 1:2){   
        
        #Loop over each lifetime egg production value
        for(h in 1:length(allLEP)){
          
          #Loop over each target
          for(f in 1:length(alltargets)){
            
            #Subset of the results for a given lifetime egg production value, target, and dispersal pattern
            targDFssub<-targDFs[which(targDFs$target==alltargets[f] & targDFs$LEP==allLEP[h] & targDFs$casestudy==allcases[x]),]
            
            #Vector of unique top methods
            alluniquetops<-unique(targDFssub$maxtype)
            
            #Loop over each top method
            for(co in 1:length(alluniquetops)){
              
              #If all dispersal distances have the same top method, then assign "all"
              if(length(which(targDFssub$maxtype==alluniquetops[co]))==length(unique(targDFs$fish))){
                finDFsapp<-data.frame("target"=c(alltargets[f]),"LEP"=c(allLEP[h]),"maxtype"=c(alluniquetops[co]),"forwho"=c("all"))
              } else {
                
                #If different dispersal distances have different top methods, then assign the name of each
                forthese<-targDFssub$fish[which(targDFssub$maxtype==alluniquetops[co])]   
                finDFsapp<-data.frame("target"=rep(alltargets[f],length(forthese)),"LEP"=rep(allLEP[h],length(forthese)),"maxtype"=rep(alluniquetops[co],length(forthese)),"forwho"=paste(as.character(forthese),sep=" ",collapse=","),"casestudy"=allcases[x])  
              }
              
              #Append results to master file
              finDFs<-rbind(finDFs,finDFsapp)
            }
          }
        }
      }
      
      #Vector of unique results
      finDFsunique<-unique(finDFs)
      
      #Get rid of "rw0.1" in names
      finDFsunique$forwho<-    gsub("rw0.1","",    finDFsunique$forwho)
      
      #Re-label the connectivity method names
      finDFsunique$maxtype<-dplyr::recode(finDFsunique$maxtype,"Spatial dependency"="SpaDep",
                                          "Google PageRank"="GooPag",
                                          "In-flow"="InFlo",
                                          "Local retention"="LocRet",
                                          "Eigenvector centrality"="EigCen",
                                          "Out-flow"="Flow",
                                          "In-degree"="InDeg",
                                          "In-flow"="InFlo",
                                          "Out-degree"="Deg",
                                          "Betweenness centrality"="BetCen")
      
      #New column for different sorting of labels
      finDFsunique$forwho2<-finDFsunique$forwho
      
      #New column for different sorting of labels
      finDFsunique$lowest<-finDFsunique$forwho
      
      #Loop for each row fo result
      for(i in 1:nrow(finDFsunique)){  
        
        #The dispersal distance
        sub2<-finDFsunique$forwho[i]
        
        #Split the dispersal distance if there is more than one
        splme<-str_split(sub2,",")[[1]]
        
        #If more than one dispersal distance
        if(length(splme)>1){
          
          #If there is "rw" in names
          if(length(grep("rw",splme))!=0){
            splme2<-splme
            
            #Remove label "rw"
            splme2<-gsub("rw0.1","",splme2)
            
            #Dispersal distances as numeric vector
            splme<-as.numeric(splme2)
            
          } else {
            #If no "rw" in names, dispersal distances as numeric vector
            splme<-as.numeric(splme)
          }
          
          #If top method is same for all distances
          if((splme[length(splme)]-splme[1])==((length(splme)-1)*50)){
            
            #Then label as "50-250"
            finDFsunique$forwho2[i]<-paste0(splme[1],"-",splme[length(splme)])
          }
          
          #Assign the name of method
          finDFsunique$lowest[i]<-splme[1]
        }
      }
      
      #Modify the label to show in parentheses for which dispersal distance the method is the best
      finDFsunique$label<-paste0(finDFsunique$maxtype," (",finDFsunique$forwho2,")")
      
      #New column for placement on x-axis
      finDFsunique$LEPplacement<-NA
      
      #ID of unique rows
      finDFsunique$uniqrow<-1:nrow(finDFsunique)
      
      #Dataframe to store vectors
      dfapp1<-data.frame("target"=c(),"LEP"=c(),"maxtype"=c(),"forwho"=c(),"colormetric"=c(),"combo"=c(),"forwho2"=c(),"lowest"=c(),"label"=c(),"LEPplacement"=c(),"uniqrow"=c(),"casestudy"=c(),"LEPplacement2"=c())
      
      #Loop over each dispersal pattern
      for(x in 1:2){
        
        #Loop over each lifetime egg production value
        for(h in 1:length(allLEP)){
          
          #Loop over each habitat target
          for(f in 1:length(alltargets)){
            
            #Subset to the target, lifetime egg production, and dispersal pattern
            finDFsuniquesub<-finDFsunique[which(finDFsunique$target==alltargets[f] & finDFsunique$LEP==allLEP[h] & finDFsunique$casestudy==allcases[x]),]
            
            #Column to store mean dispersal distance
            finDFsuniquesub$average<-NA
            
            #If there is only one top method
            if(nrow(finDFsuniquesub)==1){
              
              #Placement on x-axis is LEP value
              finDFsuniquesub$LEPplacement<-finDFsuniquesub$LEP
              finDFsuniquesub$LEPplacement2<-finDFsuniquesub$LEP
              finDFsuniquesub$average<-mean(c(50,100,150,200,250))
            } else {
              
              #Sets distance between lines
              stepsize=3.5
              
              #Loop over each top method and get the average dispersal distance for y-axis position
              for(kvf in 1:nrow(finDFsuniquesub)){
                finDFsuniquesub$average[kvf]<-mean(as.numeric(str_split(finDFsuniquesub$forwho[kvf],",")[[1]]))
              }
              
              #Order the methods by mean of dispersal distance
              finDFsuniquesub<-finDFsuniquesub[order(finDFsuniquesub$average),]
              
              #If there is an even number of top methods
              if(nrow(finDFsuniquesub)%%2==0){
                
                #Number of rows in results
                nrows<-nrow(finDFsuniquesub)
                
                #Upper limit for x-axis placement
                upper<-finDFsuniquesub$LEP[1]+stepsize*nrows/2-stepsize/2

                #Lower limit for x-axis placement
                lower<-finDFsuniquesub$LEP[1]-stepsize*nrows/2+stepsize/2
              
                #Sequence of x-axis placement  
                finDFsuniquesub$LEPplacement2<-seq(from=lower,to=upper,by=stepsize)
                
                #Loop over each result
                for(i in 1:nrow(finDFsuniquesub)){
                  
                  #Assign the placement on x-axis
                  finDFsunique$LEPplacement[which(finDFsunique$uniqrow==finDFsuniquesub$uniqrow[i])]<-finDFsuniquesub$LEPplacement[i]
                }
               
                #If there is an odd number of top methods
              } else {
                #Number of rows in results
                nrows<-nrow(finDFsuniquesub)
                
                #Upper limit for x-axis placement
                upper<-finDFsuniquesub$LEP[1]+stepsize*nrows/2-stepsize/2
                
                #Lower limit for x-axis placement
                lower<-finDFsuniquesub$LEP[1]-stepsize*nrows/2+stepsize/2
                
                #Sequence of x-axis placement  
                finDFsuniquesub$LEPplacement2<-seq(from=lower,to=upper,by=stepsize)
              }   
            }
            
            #Append the results to the master file
            dfapp1<-rbind(dfapp1,finDFsuniquesub)
          }
        }
      }
      
      #Column to count for how many dispersal distances a method was best
      dfapp1$howmany<-NA
      
      #Loop over each row of result
      for(i in 1:nrow(dfapp1)){ 
        
        #Split the dispersal distances
        splitme<-str_split(dfapp1$forwho[i],",")[[1]]
        
        #If method was best for all
        if(all(splitme=="all")){
          
          #Method was best for all five dispersal distances
          dfapp1$howmany[i]<-5
        } else {
          
          #If not best for all, count how many it was best for
          dfapp1$howmany[i]<-length(splitme)
        }
      }
    }
 
    #Copy of results dataframe
    dfapp1redone<-dfapp1
    
    #Subset to relevant columns  
    dfapp1redone<-dfapp1redone[,c("target","LEP","forwho","casestudy","maxtype")]
    
    #Empty dataframe to store placement of labels  
    dfapp1redonelong<-data.frame("target"=c(),"LEP"=c(),"forwho1"=c(),"casestudy1"=c(),"maxtype"=c(),"targetplace"=c(),"placeLEP"=c())
    
    #Re-label the best method name  
    dfapp1redone$maxtype<-recode(dfapp1redone$maxtype,
                                   "EigCen"="aEigCen",
                                   "GooPag"= "bGooPag",
                                   "Flow"="cFlow",
                                   "SpaDep"="dSpaDep")
      
    #This determines the spacing in the x-dimension
    stepsize<-0.007
    
    #This determines the spacing in the y-dimension
    steplep<-10
    
    #Re-order the results file
    dfapp1redone<-dfapp1redone[order(dfapp1redone$maxtype),]
    
    #Loop over each row of results
    for(i in 1:nrow(dfapp1redone)){

      #Subset to a given target, lifetime egg production, and case study
      dfapp1redonelongsub<-dfapp1redone[which(dfapp1redone$target==dfapp1redone$target[i]&dfapp1redone$LEP==dfapp1redone$LEP[i]&dfapp1redone$casestudy==dfapp1redone$casestudy[i]),]
      
      #Vector of unique top methods
      allthese<-unique(dfapp1redonelongsub$maxtype)
      
      #Which dispersal distance the method was best for
      forwhos<-str_split(dfapp1redone$forwho[i],",")[[1]]
      
      #If only one top method
      if(length(allthese)==1){
        
        #Then placement on x-axis is the target amount
        targetplace<-dfapp1redone$target[i]
      } else {
        
        #If there is an even number of top methods
        if(length(allthese)%%2==0){
          
          #Number of rows of methods
          nrows<-length(allthese)
          
          #Upper limit for placement on x-axis
          upper<-dfapp1redone$target[i]+stepsize*nrows/2-stepsize/2

          #Lower limit for placement on x-axis
          lower<-dfapp1redone$target[i]-stepsize*nrows/2+stepsize/2
          
          #Vector of placement for x-axis
          targetplace<-seq(from=lower,to=upper,by=stepsize)
        
          #If there is an odd number of top methods
        } else {
          
          #Number of rows of methods
          nrows<-length(allthese)
          
          #Upper limit for placement on x-axis
          upper<-dfapp1redone$target[i]+stepsize*nrows/2-stepsize/2
          
          #Lower limit for placement on x-axis
          lower<-dfapp1redone$target[i]-stepsize*nrows/2+stepsize/2
          
          #Vector of placement for x-axis
          targetplace<-seq(from=lower,to=upper,by=stepsize)
        }   
      }
    
      #Assign where the top method label should be placed
      targetplace<-targetplace[which(dfapp1redone$maxtype[i]==allthese)]
      
      #If instead we want to align the methods vertically so they are always matching, use this vector for placement
      theseplaces<-seq(from=dfapp1redonelongsub$target[1]-stepsize*3,to=dfapp1redonelongsub$target[1]+stepsize*3,length.out=length(unique(dfapp1redone$maxtype)))
      
      #This sets the placement of the label based on what the top method is
      if(dfapp1redone$maxtype[i]=="aBetCen"){
        targetplace<-theseplaces[1]
      } else {
        if(dfapp1redone$maxtype[i]=="aEigCen"){
          targetplace<-theseplaces[1]
        } else {
          if(dfapp1redone$maxtype[i]=="bGooPag"){
            targetplace<-theseplaces[2]
          } else {
            if(dfapp1redone$maxtype[i]=="dLocRet"){
              targetplace<-theseplaces[4]
            } else {
              if(dfapp1redone$maxtype[i]=="eDeg"){
                targetplace<-theseplaces[5]
              } else {
                if(dfapp1redone$maxtype[i]=="cFlow"){
                  targetplace<-theseplaces[3]
                } else {
                  if(dfapp1redone$maxtype[i]=="dSpaDep"){
                    targetplace<-theseplaces[4]
                  }
                }
              }
            }
          }
        }
      }
      
      #Long dataframe storing the placement of labels
      dfapp1redonelongapp<-data.frame("target"=rep(dfapp1redone$target[i],length(forwhos)),"LEP"=rep(dfapp1redone$LEP[i],length(forwhos)),"forwho"=forwhos,"casestudy1"=rep(dfapp1redone$casestudy[i],length(forwhos)),"maxtype"=rep(dfapp1redone$maxtype[i],length(forwhos)),"targetplace"=targetplace,"placeLEP"=NA)
      
      #Loop over each row of dataframe
      for(qrt in 1:nrow(dfapp1redonelongapp)){
        
        #Upper limit of y-axis placement
        lepup<-dfapp1redonelongapp$LEP[1]+steplep
      
        #Lower limit of y-axis placement
        lepdown<-dfapp1redonelongapp$LEP[1]-steplep
        
        #Vector of y-axis placement
        thisseq<-seq(from=lepdown,to=lepup,length.out=5)
        
        #Assign the y-axis position based on what the dispersal distance is
        if(dfapp1redonelongapp$forwho[qrt]==50){dfapp1redonelongapp$placeLEP[qrt]<-thisseq[1]}
        if(dfapp1redonelongapp$forwho[qrt]==100){dfapp1redonelongapp$placeLEP[qrt]<-thisseq[2]}
        if(dfapp1redonelongapp$forwho[qrt]==150){dfapp1redonelongapp$placeLEP[qrt]<-thisseq[3]}
        if(dfapp1redonelongapp$forwho[qrt]==200){dfapp1redonelongapp$placeLEP[qrt]<-thisseq[4]}
        if(dfapp1redonelongapp$forwho[qrt]==250){dfapp1redonelongapp$placeLEP[qrt]<-thisseq[5]}
      }
      
      #Append the results to the master file
      dfapp1redonelong<-rbind(dfapp1redonelong,dfapp1redonelongapp)
    }
    
    #Vector of colours for points
    colormethis<-gg_color_hue(6)
    
    #Don't want the 3rd and 4th colours, so this matches the figure for case studies
    colormethis<-colormethis[-c(3,4)]
    
    #Labels for methods
    labelmelikethis<-c("Eigenvector centrality","Google PageRank","Flow","Spatial dependency")
    
    #Methods as factor with levels
    dfapp1redonelong$maxtype<-factor(dfapp1redonelong$maxtype,levels=c("aEigCen", "bGooPag", "cFlow",   "dSpaDep"))
    
    #Plot showing which connectivity methods performed best at a given habitat target (x-axis) and value for lifetime egg production (column) 
    pTOPMETHOD<-ggplot()+geom_point(data=dfapp1redonelong,aes(x=targetplace,y=placeLEP,shape=maxtype,color=maxtype),size=2)+geom_hline(yintercept=seq(from=-12.5,to=87.5,by=5),lty=2,color="gray88")+geom_vline(xintercept=seq(from=0.075,to=.35,by=.05),lty=1,color="gray58")+geom_hline(yintercept=seq(from=12.5,to=87.5,by=25),lty=1,color="gray58")+ scale_x_continuous(expand = c(0, 0),limits=c(0.025,0.325),breaks=c(0.05,0.1,0.15,0.2,0.25,0.3),labels=seq(from=5,to=30,by=5)) + scale_y_continuous(expand = c(0, 0),limits=c(-12.5,87.5),breaks=c(0,25,50,75))+theme_classic()+xlab("Habitat protection (%)")+ylab("LEP of non-reserves (%)")+scale_color_manual(values=colormethis,name="Method",labels=labelmelikethis)+theme(legend.position="bottom", strip.text = element_blank(),strip.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text=element_text(size=10))+annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf)+annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf)+facet_wrap(~casestudy1)+annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+scale_shape_manual(values=c(0,15,2,19),labels=labelmelikethis,name="Method")+geom_line(data=dfapp1redonelong,aes(x=targetplace,y=placeLEP,color=maxtype,group=interaction(maxtype, target,LEP,groupagain)),size=1)+guides(shape=guide_legend(ncol=2),color=guide_legend(ncol=2))
  } 
}

#=#=#=#=#=#=#=#=#=#=#=#= END =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
