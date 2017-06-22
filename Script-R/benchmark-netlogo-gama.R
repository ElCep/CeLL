library(gamar)
defpath("~/app/Gama1.7.2RC/")



experiment1 <- getmodelparameter("~/gama_workspace/CeLL/models/Lobesia.gaml","Rexperiment")

getoutputnames(experiment1)

experiment1 <- setparametervalue(experiment1,"nblobesias",100)
experiment1 <- setparametervalue(experiment1,"visibility_other",1)
experiment1 <- setparametervalue(experiment1,"rotation",90)
experiment1 <- setparametervalue(experiment1,"ageOfDie",16)
experiment1 <- setparametervalue(experiment1,"nb_egg",10)
experiment1 <- setparametervalue(experiment1,"ageOfDie",16)

experiment1 <- setparametervalue(experiment1,"rational_traps",0)

# //Parameter agents behaviors
# parameter "Number of L. botrana" var:nblobesias  min: 10 max: 1000;
# parameter "Dis. Vision" var:visibility_other <- 1.0 min: 1.0 max: 5.0;
# parameter "rotation" var:rotation <- 170 min: 20 max: 180;
# parameter "Adult age of die" var:ageOfDie <- 16 min: 5 max: 20;
# parameter "Num. of egg laid" var:nb_egg <- 10 min: 5 max: 40;
# //Parameter trap 
# parameter "Senarii" var:rational_traps min: 0 max: 4;
# parameter "Confusion" var:confusion <- true;
# parameter "Init. Disp" var:radius_hotSpots <- 8.0 min: 5.0 max: 20.0;


experiment1 <- setfinalstep(experiment1,180)



experimentplan <- addtoexperimentplan(experiment1)

time_init <- Sys.time()
output <- runexpplan(experimentplan,1)
time_exp <- Sys.time() - time_init

with(output[[1]],{
  plot(nb_lobesias,step,type="l",lwd=2,col="red",
       ylab="number of lobesia",xlab="Step",
       ylim=c(0,70000),xlim=c(70,175), main = paste("ex. time", round(time_exp,digits = 2),"sec",sep = " "))
})
# 
# with(output[[2]],{
#   plot(step,I,type="l",lwd=2,col="red",
#        ylab="number of infected",ylim=c(0,1000))
#   lines(step,S,lwd=2,col="blue")
#   lines(step,R,lwd=2,col="green")
# })

#########################################################################################################################
## Netlogo proc
#########################################################################################################################
Sys.setenv(JAVA_HOME = '/usr/lib/jvm/jdk1.8.0_65')
library(RNetLogo)

#localisaer l'installation de netlogo
nl.path <- "/home/delaye/app/NetLogo_6.0.1/app/"
# nl.path <- "/opt/netlogo-5.3.1-64/app/"
NLStart(nl.path, gui = FALSE, nl.jarname = "netlogo-6.0.1.jar") #lance netlogo avec (TRUE) ou non (FALSE) une gui
##Definition du chemin du modèle
model.path <- "/home/delaye/github/these_ed/models/CeLL/nl602_lobesia0.6.2.nlogo"
## chargement du modele dans netlogo
NLLoadModel(model.path)


NLQuit()
