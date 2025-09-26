/**
* Name: camisol
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model cammisol

import "environment/grid.gaml"
import "nematode/nematode.gaml"

global {
	int nematodes_count <- 20;
		
	// 5E8 -> 5E9 bacterie / gramme de sol
	/*
	 * A modifier................................................ 
	 */
	 // TODO: poids total de bactérie par gramme de sol * surface du modèle
	/**
	 * Total bacteria weight in the model.
	 */
	float total_initial_bacteria_weight <-  0.05*1.5#gram/(#cm*#cm)*world.shape.area;
	// TODO: what rates?
	float init_O_rate <- 0.1;
	float init_F_rate <- 0.2;
	float init_M_rate <- 0.7;
	
	float rain_diffusion_rate <- 0.1;
	float rain_period <- 7#days;
	
	int simulation_cycle_end <- 302;
	
	bool simulationTerminee <- false;
	bool distributed_simulation <- false;
	
	init {
		
		seed <- 10.0;
		
		write("distributed_simulation ?? " + distributed_simulation);
		
		do init_grid;
		do init_enzymatic_optimisation;
		do init_enzymes;
		// Counts the number of PORES after the initialization of the grid
		int pores_count <- length(PoreParticle);
		ask PoreParticle {
			// The carrying capacity of each pore is equal to 10 times the initial bacteria population
			carrying_capacity <- 10 * total_initial_bacteria_weight / pores_count;
			
			create O_Strategist with: [C::init_O_rate * total_initial_bacteria_weight / length(PoreParticle)] {
				add self to:myself.populations;
			}
			create F_Strategist with: [C::init_F_rate * total_initial_bacteria_weight / length(PoreParticle)]{
				add self to:myself.populations;
			}
			create M_Strategist with: [C::init_M_rate * total_initial_bacteria_weight / length(PoreParticle)]{
				add self to:myself.populations;
			}
		}
		
		create Nematode number: nematodes_count
		{
			current_pore <- one_of(PoreParticle); 
			location <- any_location_in(current_pore);
		}
		
		ask PoreParticle {
			ask populations {
				do update;
				//write "Initial enzymes optimization for " + self;
				//do optimize_enzymes(myself.dam, myself.accessible_organics);
			}
		}
		
		/*int s1 <- length(MineralParticle);
		int s2 <- length(PoreParticle);
		int s3 <- (length(OrganicParticle) - length(PoreParticle));
		
		write("lenght MineralParticle " + s1);
		write("lenght PoreParticle " + s2);
		write("lenght OrganicParticle " + s3);
		
		write("lenght TOTAL " + (s1+s2+s3));
		write("Nematode " + length(Nematode));
		write("grid_size * grid_size = " + grid_size * grid_size);*/
	}

	reflex { // main loop of the model
		ask shuffle(Nematode) {
			if(scheduled) // execute only the Nematode if they are scheduled
			{			
				//write("self " + self);
				do life;
			}else
			{
				//write("" + self + " is not scheduled");
			}
		}
		ask shuffle(PoreParticle) {
			if(scheduled)
			{
				//write("self " + self);
				ask populations {
					if flip(local_step / enzymes_optimization_period) {
						do update;
						do optimize_enzymes(myself.dam, myself.accessible_organics);
					}
				}
				do decompose;
				do microbe_life;
			}else
			{
				//write("" + self + " is not scheduled");
			}
		}
	}
	
	reflex
	{
		write("running cycle " + cycle);
	} 
	
	reflex when: cycle = simulation_cycle_end
	{
		
 		//write("distribution step : --------------------------------------" + cycle);
 		//write("total_duration " + float(total_duration)/1000 + "s");
 		//write("duration " + float(duration)/1000 + "s");
 		simulationTerminee <- true;
 		write("simulation_cycle_end REACHED ");
 		write("total_duration " + float(total_duration)/1000 + "s");
		//do die;
	}
}

experiment base_cammisol_output {
	parameter "Grid size" category: "Environment" var:grid_size;
	parameter "Nematodes count" category: "Environment" var:nematodes_count;
	parameter "Organic particle rate" category: "Environment" var:organic_rate;
	parameter "Mineral particle rate" category: "Environment" var:mineral_rate;
	parameter "Init C (FOM)" category: "Environment" var:C_concentration_in_pom;
	parameter "Init N (FOM)" category: "Environment" var:N_concentration_in_pom;
	parameter "Init P (FOM)" category: "Environment" var:P_concentration_in_pom;
	parameter "Init labile rate (FOM)" category: "Environment" var:labile_rate_pom;
	
	parameter "Nematodes predation rate" category: "Nematode" var:nematode_predation_rate;
	parameter "Nematodes CUE" category: "Nematode" var:nematode_CUE;
	parameter "Nematodes C/N" category: "Nematode" var:nematode_C_N;
	parameter "Nematodes C/P" category: "Nematode" var:nematode_C_P;
	
	parameter "Dividing time (O)" category: "O Strategists" var:dividing_time_O;
	parameter "CUE (O)" category: "O Strategists" var:carbon_use_efficiency_O;
	parameter "Minimum active rate (O)" category: "O Strategists" var:minimum_active_rate_O;
	
	parameter "Dividing time (F)" category: "F Strategists" var:dividing_time_F;
	parameter "CUE (F)" category: "F Strategists" var:carbon_use_efficiency_F;
	parameter "Minimum active rate (F)" category: "F Strategists" var:minimum_active_rate_F;
	
	parameter "Dividing time (M)" category: "M Strategists" var:dividing_time_M;
	parameter "CUE (M)" category: "M Strategists" var:carbon_use_efficiency_M;
	parameter "Minimum active rate (M)" category: "M Strategists" var:minimum_active_rate_M;
	
	reflex update_particle_color {
		float max_population <- 0.0;
		ask PoreParticle {
			float population <- sum(populations collect (each.C + each.cytosol_C));
			if(population > max_population) {
				max_population <- population;
			}
		}
		if(max_population > 0.0) {
			ask Particle {
				if(type = PORE) {
					color <- rgb(0, 0, 255 * sum(PoreParticle(particle).populations collect (each.C + each.cytosol_C))/max_population);
				}
			}	
		}
	}
	
	reflex save_charts when: !distributed_simulation
	{
		if(cycle mod 10 = 0)
		{
			ask simulation 
			{	
				save (snapshot("Awoken population")) to: "../output.log/snapshot/central/awoken_" + cycle + ".png" rewrite: true;
				save (snapshot("dam")) to: "../output.log/snapshot/central/dam_" + cycle + ".png" rewrite: true;
				save (snapshot("organics")) to: "../output.log/snapshot/central/organics_" + cycle + ".png" rewrite: true;
				save (snapshot("populations")) to: "../output.log/snapshot/central/pop_" + cycle + ".png" rewrite: true;
			}	
		}
			
	}
	
	reflex save_result when: !distributed_simulation
	{
		if(cycle = 0)
		{
			save "N_dom; P_dom; C_dom; N_dim; P_dim" to: "../output.log/results_central/dam.csv" format: 'csv' rewrite: true;
	 		save "Cl; Nl; Pl; Cr; Nr; Pr" to: "../output.log/results_central/organics.csv" format: 'csv' rewrite: true;
	 		save "O_c; F_c; M_c" to: "../output.log/results_central/bacteria.csv" format: 'csv' rewrite: true;
	 		save "nematode_CO2_emissions" to: "../output.log/results_central/CO2.csv" format: 'csv' rewrite: true;
		}
		if(cycle mod 10 = 0)
		{
			ask simulation 
			{
				/*save "" + sum(Dam collect each.dom[0])/#gram +
	 			";" + sum(Dam collect each.dom[1])/#gram +
	 			";" + sum(Dam collect each.dom[2])/#gram +
	 			";" + sum(Dam collect each.dim[0])/#gram + 
	 			";" + sum(Dam collect each.dim[1])/#gram 
	 			to: "../output.log/results_central/dam.csv" format: 'csv' rewrite: false;
	 			
				save "" + sum(OrganicParticle collect each.C_labile)/#gram + 
				";" +sum(OrganicParticle collect each.N_labile)/#gram + 
				";" + sum(OrganicParticle collect each.P_labile)/#gram + 
				";" + sum(OrganicParticle collect each.C_recalcitrant)/#gram +
				";" +sum(OrganicParticle collect each.N_recalcitrant)/#gram +
				";" +sum(OrganicParticle collect each.P_recalcitrant)/#gram 
				to: "../output.log/results_central/organics.csv" format: 'csv' rewrite: false;
				
				
				let O <- PoreParticle collect each.populations[0].C;
				let F <- PoreParticle collect each.populations[1].C;
				let M <- PoreParticle collect each.populations[2].C;
				
				write("O " + O);
				write("F " + F);
				write("M " + M);
				
				save "" + sum(O)/#gram + 
				";" + sum(F)/#gram + 
				";" + sum(M)/#gram 
				to: "../output.log/results_central/bacteria.csv" format: 'csv' rewrite: false;	*/
	 			save nematode_CO2_emissions  to: "../output.log/results_central/CO2.csv" format: 'csv' rewrite: false;
			}	
		}
}
	
	output {
		display grid {
			grid Particle;
			species Nematode aspect: red_dot;
		}
		
		display "Awoken population" type: java2D {
			chart "Awoken population" type: series {
				if (length(PoreParticle) > 0) {
					data "O awake (%)" value: sum(O_Strategist collect (each.active_rate))/length(O_Strategist) * 100 style:spline color: #red marker:false thickness:3;
					data "A awake (%)" value: sum(F_Strategist collect (each.active_rate))/length(F_Strategist) * 100 style:spline color: #green marker:false thickness:3;
					data "M awake (%)" value: sum(M_Strategist collect (each.active_rate))/length(M_Strategist) * 100 style:spline color: #blue marker:false thickness:3;
				}
				if(nematodes_count > 0) {
					data "Nematode awake (%)" value: (sum(Nematode collect (each.awake as int)) / length(Nematode)) * 100 style:spline color: #yellow marker:false thickness:3;				
				}
			}
		}
			
		display "dam" type: java2D {
			chart "dam" type:series {
				data "N dom (g)" value: sum(Dam collect each.dom[0])/#gram style:spline marker:false thickness:3;
				data "P dom (g)" value: sum(Dam collect each.dom[1])/#gram style:spline marker:false thickness:3;
				data "C dom (g)" value: sum(Dam collect each.dom[2])/#gram style:spline marker:false thickness:3;
				data "N dim (g)" value: sum(Dam collect each.dim[0])/#gram style:spline marker:false thickness:3;
				data "P dim (g)" value: sum(Dam collect each.dim[1])/#gram style:spline marker:false thickness:3;
			}
		}
		
		display "organics" type:java2D {
			chart "Organics composition" type:series {
				data "C labile (g)" value: sum(OrganicParticle collect each.C_labile)/#gram style:spline marker:false thickness:3;
				data "N labile (g)" value: sum(OrganicParticle collect each.N_labile)/#gram style:spline marker:false thickness:3;
				data "P labile (g)" value: sum(OrganicParticle collect each.P_labile)/#gram style:spline marker:false thickness:3;
				data "C recalcitrant (g)" value: sum(OrganicParticle collect each.C_recalcitrant)/#gram style:spline marker:false thickness:3;
				data "N recalcitrant (g)" value: sum(OrganicParticle collect each.N_recalcitrant)/#gram style:spline marker:false thickness:3;
				data "P recalcitrant (g)" value: sum(OrganicParticle collect each.P_recalcitrant)/#gram style:spline marker:false thickness:3;
			}
		}
		
		display "populations" type:java2D {
			chart "Bacteria populations" type:series {
				data "O (g)" value: sum(O_Strategist collect each.C)/#gram style:spline color: #red marker:false thickness:3;
				data "F (g)" value: sum(F_Strategist collect each.C)/#gram style:spline color: #green marker:false thickness:3;
				data "M (g)" value: sum(M_Strategist collect each.C)/#gram style:spline color: #blue marker:false thickness:3;
			}
		}
	}
}

experiment test_microbes parent:base_cammisol_output {
	float C_soil <- 1#gram/(#cm#cm);
	float N_soil <- 0.1#gram/(#cm#cm);
	float P_soil <- 0.05#gram/(#cm#cm);
	
	parameter "Grid size" var:grid_size;
	parameter "Nematodes" var:nematodes_count;
	parameter "C soil" var:C_soil;
	parameter "N soil" var:N_soil;
	parameter "P soil" var:P_soil;
	parameter "Initial O population rate" var:init_O_rate;
	parameter "Initial F population rate" var:init_F_rate;
	parameter "Initial M population rate" var:init_M_rate;
	
	parameter "Enzymes optimization period" var:enzymes_optimization_period init:local_step;
	
	init {
		ask simulation {
			do init_POM(myself.C_soil, myself.N_soil, myself.P_soil);		
		}
	}
}

experiment display_grid {
	bool show_nematodes;
	parameter "Show nematodes" var: show_nematodes <- false;
	parameter "Grid size" category: "Environment" var:grid_size;
	
	output {
		display grid type:2d axes:false {
			grid Particle;
			graphics legend {
				draw square(cell_size) at: {soil_size + cell_size, 1.5*cell_size, 0} color: #yellow;
				draw square(cell_size) at: {soil_size + cell_size, 2.8*cell_size, 0} color: #green;
				draw square(cell_size) at: {soil_size + cell_size, 4.1*cell_size, 0} color: #black;
				draw "Mineral particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 1.7*cell_size, 0} color: #black;
				draw "Organic particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 3*cell_size, 0} color: #black;
				draw "Pore particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 4.3*cell_size, 0} color: #black;
				if show_nematodes {
					draw circle(cell_size/4) at: {soil_size + cell_size, 5.4*cell_size, 0} color: #red;
					draw "Nematode" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 5.5*cell_size, 0} color: #black;
				}
			}
			species Nematode aspect: red_dot visible: show_nematodes;
		}
	}
}

experiment cammisol_no_output {

	reflex state when: local_cycle mod 100 = 0 {
		ask simulation {
			write "Time: " + time/#day;
			write "Dom: N=" + sum(Dam collect each.dom[0]) + ", P=" + sum(Dam collect each.dom[1]) + ", C=" + sum(Dam collect each.dom[2]);
			write "Dim: N=" + sum(Dam collect each.dim[0]) + ", P=" + sum(Dam collect each.dim[1]);
			write "Active nematodes: " + length(Nematode where each.awake);
			write "";
		}
	}
	
	reflex when: (local_time>6#month) {
		ask simulation {
			do pause;
		}
	}
}

experiment cammisol parent:base_cammisol_output {
	// Test to reactivate bacterias
	reflex add_N_P when: local_cycle mod 300 = 0 {
		ask simulation {
			// Redistributes the initial quantity of C/N/P in all pores.
			// Notice that this has no scientific meaning, and is used only for test purpose.
			// This might however represent a kind of fertilizer added to the soil.
			int pores_count <- length(PoreParticle);
			float carbon_in_all_pore <- soil_weight * C_concentration_in_dom; 
			float azote_in_all_pore <- soil_weight * N_concentration_in_dom;
			float phosphore_in_all_pore <- soil_weight * P_concentration_in_dom;
			
			ask PoreParticle {
				ask dam {
					float carbon_in_pore <- carbon_in_all_pore / pores_count;
					float azote_in_pore <- azote_in_all_pore / pores_count;
					float phosphore_in_pore <- phosphore_in_all_pore / pores_count;
					
					dom[0] <- dom[0] + azote_in_pore;
					dom[1] <- dom[1] + phosphore_in_pore;
					dom[2] <- dom[2] + carbon_in_pore;
				}
			}
		}
	}
}

