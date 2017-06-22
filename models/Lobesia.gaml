/**
* Name: CeLL 
* Author: Patrick Taillandie et Etienne Delay a partir du modele developper et implenter 
* 		par Etienne Delay sous netlogo  l'aide Amélia Caffara.
* Description: Describe here the model and its experiments
* Tags: Tag1, Tag2, TagN
*/

model LobesiaModel 

global {
	list<lobesia> lobesias -> {lobesia where not dead(each)};
	//file temperature_raster <- file("../includes/gis/temperature_zone.tif");
	file parcelles_viti <- file("../includes/gis/parcelles_viti_zone2.shp");
	file parcelles_confu <- file("../includes/gis/parcelles2.shp");
	file tempFile <- file("../includes/temperature_banyuls/rimbau2012_delta.csv");
	file timeFile <- csv_file("../includes/temperature_banyuls/date2012.csv", " ",string);
	list<float> TempData;
	list<string> timeData;
	//definition of the date of begining of the simulation
	date starting_date <- date([2012,3]);
	//be careful, when real dates are used, modelers should not use the #month and #year values that are not consistent with them
	float step <- 1 #day; 
	string toDay ;//<- (date('2012-03-09') to date('2012-08-30')) every (#day); 
	geometry shape <- envelope(parcelles_viti);
	int beginPlantDev <- 100;
	bool climatChange <- false;
	int rotation <- 90;
	float T <- 17.0;
	float pheromone_left <- 100.0;
	int passoEnergy <- 4;
	float reprod_distance <-  0.5;
	float trapdistance <- 2.0;
	float visibility_other <- 3.0;
	float distMove <- 16.0;
	int nb_egg <- 10;
	int ageOfDie <- 16;
	float downPherom <- 2.0;
	
	
	float coeff_dist <- shape.width/250;
	float bufferTrapping <- 10.0 * coeff_dist;
	//float aggregation_factor <- 20.0 *coeff_dist;
	float setupDistTrap <- 10.0 * coeff_dist;
	float radius_trap <- 10.0 * coeff_dist;
	float distEffectiveTrap <- 15.0 * coeff_dist;
	int nblobesias <- 200;
	float radius_hotSpots <- 9.0;
	int nbcomptage;
	int fisrtComptage;
	int fisrtticks;
	bool lobesias_HotSpot <- true;
	float coeff_diff <- 0.8;
	float coeff_diff_trap <- 0.6;
	list<cell> active_cells;
	bool confusion <- true;
	int rational_traps <- 0 min: 0 max: 4;
	int nb_traps <- 0;
	int nbClandestin <- 5;
	
	map<string,list> generation_dist;
	init {
		ask cell{
			temperature_here <- 14.0 - rnd(1.0)  - 2.504503;
			//temperature_here <-max([0, first(temp_cell overlapping location).grid_value - 2.504503]) ;
		} 
		active_cells <- cell where (each.temperature_here >= 0.0);
		create parcelle from: parcelles_viti {
			cells <-  (active_cells inside self) where (each.location overlaps self);
			ask cells{
				viticole <- true;
				my_parcelle <- myself;
			}
		}
		create parcelle_conf from: parcelles_confu {
			cells <-  (active_cells inside self) where (each.location overlaps self);
			ask cells{
				confusion <- true;
				my_SubParcelle <- myself;
			}
		}
		TempData <- tempFile.contents collect (float(each));
		timeData <- timeFile.contents collect (string(each));
		//step <- (date('2012-03-09') to date('2012-08-30')) every (#day);
		nbcomptage <- 0;
  		fisrtComptage <- 0;
  		fisrtticks <- 0;
	
  		if lobesias_HotSpot {
  			list<cell> vitipatches <- cell where (each.viticole);
    		create seeds with:[location::one_of(vitipatches).location];
  		} else {
  			create lobesia number: nblobesias {
				do difference_Lobesias;	
			}
  		}
  		
		do cure_pheromone;
  		do integreted_lutte;
  		
  		ask active_cells where empty(trap overlapping (each.location buffer distEffectiveTrap)) {
       		outOftrapsInfluence <- true;
       	}
       	do update_stat;
	}
	
	user_command "creation d'un piege" {
		create trap with:[location::#user_location];
		ask active_cells {outOftrapsInfluence <- false;}
		ask active_cells where empty(trap overlapping (each.location buffer distEffectiveTrap)) {
       		outOftrapsInfluence <- true;
       	}
	}
	//=====================================================================================
	// we left the init part
	//=====================================================================================
	reflex dynamique {
		do temperatureEvolution;
  		
		ask active_cells  parallel: true{do maj_temp;}
	  	ask lobesia parallel: true{
			do dynamique1;
		}
		ask lobesia {
			do reproduction;
		}
		ask lobesia parallel: true{
			do dynamique2;
		}
		do cure_pheromone ;
		do update_stat;
		if (empty(lobesia)) {do pause;}
	}
	
	action update_stat {
		generation_dist["legend"] <- [];
		generation_dist["values"] <- [];
		loop g over: remove_duplicates(lobesias collect each.generation)
		{
			generation_dist["legend"] << string(g);
			generation_dist["values"] << lobesias count (each.generation = g);
		}  	
	}
	
	action integreted_lutte {
		 if confusion {
		 	//if not empty(trap where not each.observer) {
		 		do install_trapping;
		 		ask active_cells where not empty(trap overlapping (each.location buffer distEffectiveTrap)) {
       				outOftrapsInfluence <- false;
  	 			} 
		 	//}     
		 }
       }
	
	action install_trapping {
		  list<cell> vitipatches <- cell where each.viticole;
		  switch rational_traps {
		  	match 0 {
		  		create trap number: nb_traps with: [location::one_of (vitipatches).location];
		  	}
		  	match 1 {
		  		ask one_of(seeds) {
		  			list<cell> tappingZone <- cell at_distance (radius_hotSpots + bufferTrapping); 
		  			ask nb_traps among tappingZone {
		  				if empty(trap overlapping self) {
			  				create trap number: nb_traps with: [location::location]; 	
		  				}
		  			}
		  		}
		  	}
		  	match 2 {
		  		int surface_lute <- length(vitipatches);
		  		float surface_trap <- #pi * (radius_trap ^ 2);
		  		list<cell> myparcelle <-  vitipatches where empty(trap overlapping (each.location buffer distEffectiveTrap)) ;
		  		loop while: surface_lute > surface_trap and not empty(myparcelle) {
		  			cell parc <- one_of(myparcelle);
		  			create trap with: [location:: parc.location]{
		  				myparcelle <- myparcelle - (cell at_distance distEffectiveTrap);
		  			}
		  			surface_lute <- length(vitipatches);
		  		}
		  	}
		  	match 3 {
		  		 list<parcelle> list_parcelles <- list<parcelle>(remove_duplicates(vitipatches collect each.my_parcelle));
		  		 loop p over: list_parcelles {
		  		 	int surface_lute <- length(p.cells);
		  			float surface_trap <- #pi * (radius_trap ^ 2);
		  			list<cell> my_cells <- copy(p.cells);
		  			loop while: surface_lute > surface_trap and not empty(my_cells) {
			  			cell parc <- one_of(my_cells);
			  			create trap with: [location:: parc.location]{
			  				my_cells <- my_cells - (cell at_distance distEffectiveTrap);
			  			}
			  			surface_lute <- length(vitipatches);
			  		}
		  		}	
		  	}
		  	match 4 {
		  		int surface_lute <- length(vitipatches);
		  		float surface_trap <- #pi * (radius_trap ^ 2);
		  		list<cell> myparcelle <-  vitipatches where empty(trap overlapping (each.location buffer distEffectiveTrap)) ;
		  		loop while: surface_lute > surface_trap and not empty(myparcelle) {
		  			cell parc <- one_of(myparcelle);
		  			create trap with: [location:: parc.location]{
		  				myparcelle <- myparcelle - (cell at_distance distEffectiveTrap);
		  			}
		  			surface_lute <- length(vitipatches);
		  		}
		  		 list<parcelle> list_parcelles <- list<parcelle>(remove_duplicates(vitipatches collect each.my_parcelle));
		  		 ask nbClandestin among list_parcelles {
		  		 	ask cells {
		  		 		clandestin <- true;
		  		 	}
		  		 	ask trap overlapping self { do die;}
		  		 }
		  	}	
		 }
	}
  		
 
	action cure_pheromone {
		ask trap where (each.observer) {
			list<lobesia> lobesias <- (lobesia at_distance trapdistance) where (each.phase = "adult" and sex = 1);
			if (not empty(lobesias)) {
				lobesias_traps <- length(lobesias);
				ask lobesias {do die;}
			}
		}
     
	  	ask trap {
	    	ask cell(location) {
	    		Ppheromone_trap <- pheromone_left;
	    	}
    	}
    	ask active_cells parallel: true{
    		Ppheromone_tmp <- 0.0;
    		Ppheromone_trap_tmp <- 0.0;
    	}
    	
    //	diffuse var:Ppheromone on:cell proportion: coeff_diff  ;
    //	diffuse var:Ppheromone_trap on:cell proportion: coeff_diff_trap radius:1 ;
    	
    	ask active_cells parallel: true{
    		if (Ppheromone != 0) {
    			float diff <- coeff_diff * Ppheromone;
    			Ppheromone_tmp <- Ppheromone_tmp + Ppheromone - diff; 
    			diff <- diff / length(neighbors);
    			ask neighbors {Ppheromone_tmp <- Ppheromone_tmp + diff; }
    		}
    		if (Ppheromone_trap != 0) {
    			float diff <- coeff_diff_trap * Ppheromone_trap;
    			Ppheromone_trap_tmp <- Ppheromone_trap_tmp + Ppheromone_trap - diff; 
    			diff <- diff / length(neighbors);
    			ask neighbors {Ppheromone_trap_tmp <- Ppheromone_trap_tmp + diff; }
    		}
    	}
    	ask active_cells parallel: true {
    		Ppheromone <- (Ppheromone_tmp / downPherom) with_precision 2;
    		//Ppheromone <- Ppheromone / downPherom;
    		Ppheromone_trap <- Ppheromone_trap_tmp with_precision 2;
    		
    		if Ppheromone != 0 or Ppheromone_trap != 0 {
    			mixed_pheromone <- max([Ppheromone, Ppheromone_trap]);
    		} else {
    			mixed_pheromone <- 0.0;
    		}
    	}
    	
    	ask active_cells where (each.outOftrapsInfluence) {
    		Ppheromone_trap <- 0.0;
    	}
	}
	
	action temperatureEvolution {
		if (empty(TempData)) {do pause;}
		float deltaDay <- first(TempData);
		TempData >> deltaDay;
  		ask cell {
  			temperature_here <- temperature_here + deltaDay + (climatChange ? 0.5 : 0.0);
  		}
    	toDay <- first(timeData);
    	timeData >> toDay;
    	list<string> v <- toDay split_with "-";
    	toDay <- v[2]+ "-"+v[1];
	}
}


//grid temp_cell file: temperature_raster frequency: 0;

grid cell width: 250 height: 159 neighbors: 8 frequency: 0{
	float termicalAcum <- 0.0;
	float mixed_pheromone;
	float temperature_here;
	bool flower <- false;
	float  Ppheromone_trap ;
	float unitaTermicaForch;
	float flower_time;
	float  Ppheromone <- 0.0;
	float Ppheromone_tmp;
	float  Ppheromone_trap_tmp ;
	bool viticole <- false;
	bool clandestin <- false;
	rgb color <- #black;
	parcelle my_parcelle;
	parcelle_conf my_SubParcelle;
	bool outOftrapsInfluence <- false;
	bool devplant_dyn <- false;
	
	aspect default {
		if flower {
			draw pyramid(3) color: #red;
		}
		else if (mixed_pheromone > 0.2) {
			draw shape color: color depth: 0.2 border: #black;
		}
	}
	action maj_temp {
		termicalAcum <- termicalAcum + temperature_here;
   		if devplant_dyn {do devPlant;}
   		do color_patches;
   	}
   	
   	action devPlant {
   		if cycle > beginPlantDev  {
   			if unitaTermicaForch >= flower_time {
   		  		flower <- true;
   		  	}
   		  }
   		  if unitaTermicaForch > flower_time + 5 {
	    	flower <- false;
	     }
	  }
	   	
   	action color_patches {
   		int val <- round(mixed_pheromone * 10 );
   		color <- rgb(0,0,val);
   	}
}

species parcelle frequency: 0{
	list<cell> cells;
	aspect default {
		draw shape color: #gray border: #black depth: 0.1;
	}
	aspect contours {
		draw shape.contour color: #red depth:10;
	}
}

species parcelle_conf frequency: 0{
	list<cell> cells;
	aspect default {
		draw shape color: #red empty: true width: 2 border: #red;
	}
}

species seeds {
	init {
		//create trap with:[location::location];
		list<cell> cells_neighbors <- cell at_distance  radius_hotSpots;
		create lobesia number: nblobesias with:[location :: any_location_in(one_of(cells_neighbors))] {
			do difference_Lobesias;	
		}
	
      	
	}
	aspect default {
		draw triangle(2) color: #violet border: #black;
	}
}

species trap frequency: 0{
	bool observer <- false;
	int lobesias_traps;
	aspect default {
		draw pyramid(3) color: #red border: #black;
	}
}

species lobesia skills: [moving] frequency: 0{
	rgb color <- #red;
    string phase <- "egg";//gg larva crisalide adageOfDieulte
    int phaseNumber <- 1;
    int energy <- 50; // how much I'm hangry
    int age <- 0 ; // max 20 jours
    int sex <- rnd(1); //  ;;male = 0 femmelle = 1
    int gestationTime <- 0;
    int heatrequired <- 0;
    float a <- 0.397370;
    float b <- 0.183374;
    float c <- 0.187975;
    int generation <- 0;
    float unitaDev;
    float t;
    bool reproductionFlag;
    float nb_mating;
    
    aspect default {
    	if (phaseNumber >= 4) {
    		draw file("../images/bug-1.gif") size: 5 at: location + {0,0,1} ;
			//draw obj_file("../images/Insect.obj", -90::{1,0,0}) color: color size: 1 at: location + {0,0,1} ;
    	}
    	else {
    		draw sphere(0.5) color: color;
    	} 
    }
    action difference_Lobesias {
    	float normalDist <- gauss(50, 47.5 / 2);
    	if normalDist <= 5 {
 	    	heatrequired <- 351;
    	} else if normalDist > 5 and  normalDist <= 25 {
    		heatrequired <- 374;
    	} else if normalDist > 25 and  normalDist <= 50 {
    		heatrequired <- 401;
    	} else if normalDist > 50 and  normalDist <= 75 {
    		heatrequired <- 442;
    	} else if normalDist > 75 and  normalDist <= 95 {
    		heatrequired <- 518;
    	} else if normalDist > 95 {
    		heatrequired <- 600;
    	} 
    }
    action goTo {
    	if phase = "adulte" {
    		if sex = 1 {
    			loop times: distMove{
            		do wander amplitude: rotation speed: 1.0;
            		ask cell(location) {
       					 Ppheromone <- pheromone_left;
     				}
            	}
             }else {
    			loop times: distMove{
    				cell cible <- (cell(location) neighbors_at visibility_other) with_max_of (each.mixed_pheromone);
    				do move speed: 1.0 heading: self towards cible;
            	}
    		}
    	}
    } 
    action phaseEvolution {
    	if (cell(location).termicalAcum >= heatrequired) and (generation = 0) {
    		if phase = "egg" {
    			color <- #orange;
    			phase <- "larva";
    			phaseNumber <- 2;
    			a <- 0.335958;
    			b <- 0.195681;
    			c <- 0.197009;
    			t <- cell(location).temperature_here;
    			unitaDev <- 0.0;
    		}
    	}
    	if  phase = "egg" and generation > 0{
    		phase <- "larva";
    		color <- #orange;
    		phaseNumber <- 2;
    		a <- 0.335958;
    		b <- 0.195681;
    		c <- 0.197009;
    		t <- cell(location).temperature_here;
    		unitaDev <- 0.0;
    	}
    	if  phase = "larva" and unitaDev >= 1{
    		phase <- "crisalide";
    		color <- #brown;
    		phaseNumber <- 3;
    		a <- 0.439051;
    		b <- 0.311930;
    		c <- 0.313915;
    		t <- cell(location).temperature_here;
    		unitaDev <- 0.0;
    	}
    	if  phase = "crisalide" and unitaDev >= 1{
    		phase <- "adulte";
    		color <- #yellow;
    		phaseNumber <- 4;
    	}
    		
    }
   
    float foncDev(float tt, float aa, float bb, float cc) {
    	return aa * (exp(bb * (tt - 10)) -exp(bb * (35 - 10) - cc * (35 - tt)));
    }
    
    bool mating (float x){
    	float mating <- rnd(100.0);
  		float proba_mating <- 100 / (x + 1);
  		return mating <= proba_mating;
    }
    
    action reproduction {
    	if reproductionFlag and sex = 1 and cell(location).viticole{
    		list<lobesia> voisins_males <- (lobesia at_distance reprod_distance) where ((each.sex = 0) and (phase = "adulte") and mating(nb_mating));
    	 	if (not empty(voisins_males)) {
    	 		int myheatrequired <- heatrequired;
    	 		int mygeneration <- generation;
    	 		using topology(cell) {
    	 			list<cell> where_lay <- 5 among (cell at_distance 2);
    	 			ask where_lay with_min_of length(lobesia inside each) {
    	 				create lobesia number: nb_egg {
    	 					age <- 0;
    	 					sex <- rnd(1);//   ;;male = 0 femmelle = 1
     					   	heatrequired <- myheatrequired;
     					   	phase <- "egg";
     					   	phaseNumber <- 1;
     					   	generation <- mygeneration + 1;
           					a <- 0.335958;
							b <- 0.195681;
							c <- 0.197009;
							t <- myself.temperature_here;
        					unitaDev <- 0.0;
        				 	location <- any_location_in(myself);
    	 				}
    	 			}
    	 		}
    	 	
	    	reproductionFlag <- false;
	    	gestationTime  <- 0;
	    	nb_mating <- nb_mating + 1;
	    	}
	    }
    	else {
    		gestationTime <- gestationTime + 1;
   	 		if gestationTime > 2 {
    			reproductionFlag <- true;
    
   			}
   		}
   	}
   	
    action dynamique1 {
    	do goTo;
    	do phaseEvolution;
    	unitaDev <- unitaDev + foncDev(t, a, b, c);
    	if phase = "larva" {
  			 energy <- passoEnergy;
  		}
  	}
    
    action dynamique2 {
    	if phase = "adulte" {
    	 	age <- age + 1;
    	 	if (age > ageOfDie) {
    	 		do die;
    	 	}
    	 }
    	 do difference_Lobesias;
    }
	
}



experiment main type: gui {
	//Parameter agents behaviors
	parameter "Number of L. botrana" var:nblobesias  min: 10 max: 1000;
	parameter "Dis. Vision" var:visibility_other <- 1.0 min: 1.0 max: 5.0;
	parameter "rotation" var:rotation <- 170 min: 20 max: 180;
	parameter "Adult age of die" var:ageOfDie <- 16 min: 5 max: 20;
	parameter "Num. of egg laid" var:nb_egg <- 10 min: 5 max: 40;
	//Parameter trap 
	parameter "Senarii" var:rational_traps min: 0 max: 4;
	parameter "Confusion" var:confusion <- true;
	parameter "Init. Disp" var:radius_hotSpots <- 8.0 min: 5.0 max: 20.0;
	
	
	output {
		display carte type: opengl {
			image "../images/background.png" refresh: false;
			//grid cell lines: #black transparency: 0.0;
			species cell transparency: 0.5;
			species parcelle transparency: 0.5 refresh: false;
			species parcelle aspect: contours refresh: false;
			species trap ;
			species lobesia;
		}
		
		display charts {
			chart "pop lobesias" type: series size: {0.5,0.5} x_serie_labels: toDay{
				data "nb of lobesias" color: #gray value: length(lobesia);
				data "nb of eggs" color: #red value: lobesias count (each.phaseNumber = 1);
				data "nb of larva" color: #orange value: lobesias count (each.phaseNumber = 2);
				data "nb of crisalide" color: #brown value: lobesias count (each.phaseNumber = 3);
				data "nb of adults" color: #yellow value: lobesias count (each.phaseNumber = 4);
			}
			chart "phase evolution" type: histogram size: {0.5,0.5} position: {0.5,0.0}{
				data "nb of eggs" color: #red value: lobesias count (each.phaseNumber = 1);
				data "nb of larva" color: #orange value: lobesias count (each.phaseNumber = 2);
				data "nb of crisalide" color: #brown value: lobesias count (each.phaseNumber = 3);
				data "nb of adults" color: #yellow value: lobesias count (each.phaseNumber = 4);
			}
			chart "generation histogram" type: histogram size: {0.5,0.5} position: {0.0,0.5}{
				datalist list(generation_dist at "legend") value: list(generation_dist at "values");
			}
			chart "male/femelle" type: series size: {0.5,0.5} position: {0.5,0.5} x_serie_labels: toDay{
				data "nb of male" color: #blue value: lobesias count (each.sex = 0);
				data "nb of femelle" color: #violet value: lobesias count (each.sex = 1);
			}
			
		}
		display obs {
			chart "Temperature" type: series size: {0.5,0.5} x_serie_labels: toDay{
				data "nb of lobesias" color: #gray value: mean(cell collect (each.temperature_here));
			}
			chart "phase evolution" type: histogram size: {0.5,0.5} position: {0.5,0.0}{
				data "nb of eggs" color: #red value: lobesias count (each.phaseNumber = 1);
				data "nb of larva" color: #orange value: lobesias count (each.phaseNumber = 2);
				data "nb of crisalide" color: #brown value: lobesias count (each.phaseNumber = 3);
				data "nb of adults" color: #yellow value: lobesias count (each.phaseNumber = 4);
			}
		
		}
	}
}

experiment Rexperiment type: gui  {
	//Parameter agents behaviors
	parameter "Number of L. botrana" var:nblobesias  min: 10 max: 1000;
	parameter "Dis. Vision" var:visibility_other <- 1.0 min: 1.0 max: 5.0;
	parameter "rotation" var:rotation <- 170 min: 20 max: 180;
	parameter "Adult age of die" var:ageOfDie <- 16 min: 5 max: 20;
	parameter "Num. of egg laid" var:nb_egg <- 10 min: 5 max: 40;
	//Parameter trap 
	parameter "Senarii" var:rational_traps min: 0 max: 4;
	parameter "Confusion" var:confusion <- true;
	parameter "Init. Disp" var:radius_hotSpots <- 8.0 min: 5.0 max: 20.0;
	
	output {
		monitor "nb_lobesias" value: length(lobesia);
	}
	
	
	
	
}
