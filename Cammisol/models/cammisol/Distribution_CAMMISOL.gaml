/**
* Name: DistributionCAMMISOL
* Distribution model of CAMMISOL using GRID partitionning. 
* Author: Lucas Grosjean
* Tags: Distribution, High Performance Computing, CAMMISOL, GRID
*/

model DistributionCAMMISOL

import "cammisol.gaml" as Thematic

global
{
	int MPI_RANK <- 0;								// MPI RANK of the current  model instance
	int MPI_SIZE;									// number of MPI rank on the network
	
	int cluster_number;						// number of cluster wanted
	int grid_size <- 20;							// size of the grid
	int nematodes_count <- 20;						// number of nematodes
	
	list<list<int>> clusters;						// clusters of cells
	
	list<int> index_of_pores;						// list of the index of pores
	
	map<int,list<int>> pores_by_clusters;			// key : ID of cluster; value : list of the pores for that cluster 
	map<int,int> cell_to_cluster; 					// key : index of the cell; value : cluster ID;
	
	list<int> pores_neigh_border; 					// index of pores simulated by other proc adjacent to one of my pore
	list<int> organics_neigh_border; 				// index of organics simulated by other proc adjacent to one of my pore
	
	list<int> my_adjacent_pores;					// list of the adjacent pores simulated by this instance
	list<int> my_border_organics;					// list of the border cells simulated by this instance
	
	list<int> my_cells;								// list of the cells simulated by this instance
	list<int> my_pores;								// list of the pores simulated by this instance
	list<Nematode> my_nematodes <- list<Nematode>([]);							// list of the nematodes simulated by this instance
	
	map<int,list<int>> my_border_organics_neighbors; // key : MPI RANK; value : list of organics to send
	
	map<Nematode,int> nematode_to_cell;				// key : Nematode agent; value : cell of the key nematode
	
	list<int> grid_score;
		
	init
 	{	
 		write("seed " + seed);
 		seed <- 33.0;
 		 
 		create Communication_Agent_MPI; 					// init of the communication agent
 		MPI_RANK <- Communication_Agent_MPI[0].MPI_RANK;	// get the MPI Rank of this instance
 		MPI_SIZE <- Communication_Agent_MPI[0].MPI_SIZE;	// get the size of the MPI Network
 		
 		create Partitionning_Agent;
 		create Synchronization_Agent;
 		
 		cluster_number <- MPI_SIZE; // number of cluster = MPI_SIZE
 		
 		write("MPI_RANK " + MPI_RANK);
 		write("MPI_RANK " + MPI_SIZE);
 		
 		create Thematic.cammisol with: [grid_size:grid_size, nematodes_count:nematodes_count, seed:seed, distributed_simulation: true];
 		
 		ask Thematic.cammisol[0]
 		{
 			write("..................INIT CAMMISOL........................");
 			myself.grid_size <- grid_size;
 			grid_score <- Particle collect each.score;
 			loop particle over: list(Particle)
 			{
				myself.particle_color[particle.index].type <- particle.type; // cammisol coloring
 				if(particle.type = "pore")
 				{
 					index_of_pores << particle.index; 
 				}
 			}
 		}
 		
 		ask Partitionning_Agent
 		{
			//clusters <- grid_grid_partitioning(grid_size, grid_size, cluster_number,4);
			//clusters <- grid_KMEAN_partitionning(grid_size, grid_size, cluster_number, 4); // todo only cluster by 0 + scatter
			clusters <- grid_score_partitioning(grid_size, grid_size, grid_score, cluster_number);
			
	 		//clusters <- grid_grid_partitioning(grid_size, grid_size, cluster_number, 4);
	 		my_cells <- clusters[MPI_RANK];
	 		
			do coloring(clusters);					
			write("grid_KMEAN_paritionning");
 			my_pores <- pores_by_clusters[MPI_RANK];
 			
 			do color_nematode_cell;
			do border_cell;
 			do unschedule;
	 		do unschedule_color;
	 		
			write("pores_by_clusters " + pores_by_clusters);
			write("cell_to_cluster " + cell_to_cluster);
			write("pores_neigh_border " + pores_neigh_border);
			write("organics_neigh_border " + organics_neigh_border);
			
			write("my_border_organics " + my_border_organics );
			write("my_adjacent_pores " + my_adjacent_pores);
	 		write("my_cells" + my_cells);
	 		write("my_pores" + my_pores);
	 		write("my_nematodes "  + my_nematodes);
	 		
	 		
			write("pourcent of border_cells " + (length(organics_neigh_border)/(grid_size * grid_size))*100 );
	 		write("my_border_cells_neighbors " + my_border_organics_neighbors);
		}
 	}
 	
 	action end_simulations
	{
 		if(check_end_simulation())
 		{
 			//write("clusters : " + sample(clusters));
 			
 			write("DISTRIBUTION MODEL IS OVER");
 			write("total_duration " + float(total_duration)/1000 + "s");
 			
 			ask Communication_Agent_MPI
 			{
 				do MPI_BARRIER();
 			}
 			ask Thematic.cammisol[0].simulation
 			{
 				do die;
 			}
	 		do die;
	 		
 		}
	}
	
	reflex distributed_main
 	{ 
 		do run_thematic_model(); 		// run CAMMISOL
 		do end_simulations();			// check end of the model
 		
 		
 		ask Synchronization_Agent
 		{
	 		ask Thematic.cammisol[0]
	 		{		
		 		write("lenght organics before " + length(OrganicParticle));
	 		} 
	 		if(cycle mod 10 = 0)
	 		{ 			
	 			do send_my_border_cell;
	 		}
	 		ask Thematic.cammisol[0]
	 		{		
		 		write("lenght organics after " + length(OrganicParticle)); // todo fix this
	 		} 
	 		
	 		do migrate_nematode;
 		}
 	}
 	
 	action run_thematic_model // run a cycle of the CAMMISOL model
 	{
 		write("" + MPI_RANK + " distribution step : --------------------------------------" + cycle);
 		write("total_duration " + float(total_duration)/1000 + "s");
 		write("duration " + float(duration)/1000 + "s");
 		if(check_end_simulation())
 		{
 			return; // CAMMISOL is over we don't need to run step again
 		}
 		write("RUNNING THEMATIC");
 		ask Thematic.cammisol[0].simulation
 		{
 			do _step_; // run a step of sub model
 		}
 	}
 	
	bool check_end_simulation
 	{
 		bool end_simu <- false;
 		ask Thematic.cammisol[0].simulation
 		{
	 		if(simulationTerminee) // check if sub model is done simulating
	 		{	
	 			end_simu <- true;
	 			write(" ");
	 			write("CAMMISOL SIMULATION IS OVER");
	 		}	
 		}
 		return end_simu;
 	}
}

// fake grid for coloring and cammisol data manipulation
grid particle_color width: grid_size height: grid_size neighbors: 4
{
	rgb color_border <- #black;
	string type;
	bool nematode <- false;
	bool scheduled <- true;
	
	bool neighbor_pore_cell <- false;
	bool my_border_pore_cell <- false;
	
	bool my_border_cell <- false;
	
	/**
	 * Each pore cell check the neighborhood and verify that they are on the same cluster, if not color them and fill containes 
	  */
	aspect distributed_global
	{
		draw self color: color border: color_border;
		
		if(nematode)
		{
			draw circle(0.5) at: self.location color: #red;
		}
		
		draw ""+index color: #white font: font('Default', 8, #bold);
	}
	aspect distributed
	{
		if(neighbor_pore_cell)
		{		
			draw self color: #black border: color_border;
		}else if(my_border_pore_cell)
		{
			draw self color: #brown border: color_border;
		}else
		{
			draw self color: color border: color_border;
		}
		
		if(nematode)
		{
			draw circle(0.5) at: self.location color: #red;
		}
		draw ""+index color: #white font: font('Default', 7, #bold) at: {self.location.x-1,self.location.y};
	}
	
	aspect cammisole{
		switch type {
			match "pore" {draw self color: #black border: color_border;}
			match "mineral" {draw self color: #yellow border: color_border;}
			match "organic" {draw self color: #green border: color_border;}
		}
		draw ""+index color: #white font: font('Default', 7, #bold) at: {self.location.x-1,self.location.y};
	}
	
	aspect UNSCHEDULED{ draw self color: scheduled ? #white : #black;}
}

species Partitionning_Agent
{
	action border_cell
	{
		loop cell over: list(particle_color)
		{
			if(my_cells contains cell.index) // my cells
			{
				int my_cluster <- cell_to_cluster[cell.index];
				if(cell.type = "pore") // if neighbor on different cluster = cell that will need to be updated from other proc
				{
					loop neigh over: cell.neighbors
					{
						int neigh_cluster <- cell_to_cluster[neigh.index];
						if(my_cluster != neigh_cluster) 		// we are on different cluster
						{
							
							switch neigh.type { // type of neighbor cell
								match "pore" {
									add neigh.index to: pores_neigh_border;		// neighbor = adjacent pore (nematode move !)
									add cell.index to: my_adjacent_pores;		// neighbor = adjacent pore (nematode move !)
									neigh.neighbor_pore_cell <- true;		// neighbor is a pore
									cell.my_border_pore_cell <- true;		// my cell is a pore
								}
								match "organic" {
									add neigh.index to: organics_neigh_border;		// will need to receive update from this cell
								}
							}
						}
					}
				}else if(cell.type = "organic") // if neighbor on different cluster = cell that we will update
				{
					loop neigh over: cell.neighbors
					{
						int neigh_cluster <- cell_to_cluster[neigh.index];
						if(my_cluster != neigh_cluster) 	// we are on different cluster
						{
							switch neigh.type {
								match "pore" {
									write("self : " + self);
									add cell.index to: my_border_organics;
									
									// important use case : the current cell is an organics and is next to a pore on another cluster, we need to update that cell
									do add_my_border_cells_neighbors(cell, neigh_cluster);
								}
								match "organic" {
									write("self : " + cell);
								}
							}
						}
					}
				}
			}
		}
	}
	
	action add_my_border_cells_neighbors(particle_color cell, int neigh_cluster)
	{
		if( my_border_organics_neighbors[neigh_cluster] = nil)
		{
			my_border_organics_neighbors[neigh_cluster] <- list<int>(cell.index);
		}else
		{
			if(!(my_border_organics_neighbors[neigh_cluster] contains cell.index))
			{
				my_border_organics_neighbors[neigh_cluster] << cell.index;
			}
		}
	}
	
	action coloring(list<list<int>> clusters_to_color)
	{
		list colors <- [#green, #darkcyan, #red, #purple, #antiquewhite, #aqua, #aquamarine, #azure, #beige, #bisque, #black, #blanchedalmond, #blue, #blueviolet, #brown, #burlywood, #cadetblue, #chartreuse, #chocolate, #coral, #cornflowerblue, #cornsilk, #crimson, #cyan, #darkblue, #darkcyan, #darkgoldenrod, #darkgray, #darkgreen, #darkkhaki, #darkmagenta, #darkolivegreen, #darkorange, #darkorchid, #darkred, #darksalmon, #darkseagreen, #darkslateblue, #darkslategray, #darkturquoise, #darkviolet, #deeppink, #deepskyblue, #dimgray, #dodgerblue, #firebrick, #floralwhite, #forestgreen, #fuchsia, #gainsboro, #ghostwhite, #gold, #goldenrod, #gray, #green, #greenyellow, #honeydew, #hotpink, #indianred, #indigo, #ivory, #khaki, #lavender, #lavenderblush, #lawngreen, #lemonchiffon, #lightblue, #lightcoral, #lightcyan, #lightgoldenrodyellow, #lightgray, #lightgreen, #lightpink, #lightsalmon, #lightseagreen, #lightskyblue, #lightslategray, #lightsteelblue, #lightyellow, #lime, #limegreen, #linen, #magenta, #maroon, #mediumaquamarine, #mediumblue, #mediumorchid, #mediumpurple, #mediumseagreen, #mediumslateblue, #mediumspringgreen, #mediumturquoise, #mediumvioletred, #midnightblue, #mintcream, #mistyrose, #moccasin, #navajowhite, #navy, #oldlace, #olive, #olivedrab, #orange, #orangered, #orchid, #palegoldenrod, #palegreen, #paleturquoise, #palevioletred, #papayawhip, #peachpuff, #peru, #pink, #plum, #powderblue, #purple, #red, #rosybrown, #royalblue, #saddlebrown, #salmon, #sandybrown, #seagreen, #seashell, #sienna, #silver, #skyblue, #slateblue, #slategray, #snow, #springgreen, #steelblue, #tan, #teal, #thistle, #tomato, #turquoise, #violet, #wheat, #white, #whitesmoke, #yellow, #yellowgreen];
		int index_color <- 0;
		
		loop cluster over: clusters_to_color
		{
			loop cell over: cluster
			{
				particle_color[cell].color <- colors[index_color];
				
				if(index_of_pores contains cell)
				{
					//particle_color[cell].color <- #white; todo
					
					particle_color[cell].color <- colors[index_color];
					//particle_color[cell].color_border <- colors[index_color];
					
					if(pores_by_clusters[index_color] = nil)
					{
						pores_by_clusters[index_color]  <- list(cell);
					}else
					{
						pores_by_clusters[index_color]  << cell;
					}
				}
				cell_to_cluster[cell] <- index_color;
			}
			index_color <- index_color + 1;	
		}		
	}
	
	/**
 	 * Unschedule nematodes and cell not in my clusters
 	 */
 	action unschedule
 	{
 		loop nematode over: nematode_to_cell.pairs 
 		{
 			write("nematode might be unschedleudl ? " + nematode);
 			if(!(my_pores contains nematode.value))
 			{
 				write("UNSCHELUDED " + nematode);
 				nematode.key.scheduled <- false;
 			}
 		}
 		
 		int index_pore_particle <- 0;
 		list<int> index_of_pore_to_unschedule;
 		write("my_pores index_of_pore_to_unschedule " + my_pores);
 		loop pore over: index_of_pores
 		{
 			if(!(my_pores contains pore))
 			{
 				index_of_pore_to_unschedule << index_pore_particle;
 			}
 			index_pore_particle <- index_pore_particle + 1;
 		}
 		
		ask Thematic.cammisol[0]
 		{
 			loop index_unschedule over: index_of_pore_to_unschedule
 			{
 				PoreParticle[index_unschedule].scheduled <- false;
 				//write("UNSCHELUDED " + PoreParticle[index_unschedule]);
 			}
 		}
 	}
 
 	/**
  	* color Nematode position and store them in  nematode_to_cell
  	*/
	action color_nematode_cell
	{
		write("my_poresmy_poresmy_poresmy_pores " + my_pores);
 		ask Thematic.cammisol[0]
 		{	
 			let tmp_my_pores <- my_pores;
			ask Nematode
			{
				nematode_to_cell[self] <- index_of_pores[current_pore.index]; // can't directly access particle_color in this context
				write("" + self + "current_pore.index " + current_pore.index);
				write(" " + self + " currurur " + index_of_pores[current_pore.index]);
				if(tmp_my_pores contains index_of_pores[current_pore.index])
				{
					my_nematodes << self;
				}
			}
			write("tmp_my_pores " + tmp_my_pores);
		}
		loop cell over: nematode_to_cell
		{
			particle_color[cell].nematode <- true;	
		}
		
		write("init my_nematodes " + my_nematodes);
	}
	
	/**
	 * Cell not in my cluster are unscheduled
	  */
	action unschedule_color
	{
		loop index_cell from: 0 to: length(particle_color) { 
			if(!(my_cells contains index_cell))
			{
				particle_color[index_cell].scheduled <- false;
			}else
			{
				particle_color[index_cell].scheduled <- true;
			}
		}
	}
}

species Communication_Agent_MPI skills:[MPI_SKILL]{} // communication agent with MPI Librairy

species Synchronization_Agent
{	
	map<int,list<OrganicParticle>> organics_to_send; 	// key : MPI RANK; value : list of OrganicParticle agent to send
	map<int,list<OrganicParticle>> updated_organics;	// key : MPI RANK; value : list of OrganicParticle agent with data
	
	action send_my_border_cell
	{
		// send my_border_cells
		loop cell over: my_border_organics_neighbors.pairs 
		{
			write("MPI RANK " + cell.key + " need to be receive organic " + cell.value);
			ask Thematic.cammisol[0]
 			{
				loop organic over: cell.value
				{
					write("Particle " + Particle[organic] + " :: " + Particle[organic].particle);
					
					if(myself.organics_to_send[cell.key] = nil)
					{
						myself.organics_to_send[cell.key]  <- list<OrganicParticle>((Particle[organic].particle as OrganicParticle));
					}else
					{
						myself.organics_to_send[cell.key]  << Particle[organic].particle as OrganicParticle;
					}
				}	
			}
		}
		
		write("organics_to_send " + organics_to_send);
		
		ask Communication_Agent_MPI
		{
	    	myself.updated_organics <- MPI_ALLTOALL(myself.organics_to_send); // MPI all to all	
	    	write("data_recv " + myself.updated_organics);
		}
		
		do update_border_cell; // update the cell with the new values
		do kill_updated_cell;
		
		organics_to_send <- nil; // emptying container
		updated_organics <- nil; // emptying container
		
	}
	
	action kill_updated_cell
	{
		loop killer over: updated_organics
		{
			write("killer " + killer);
			loop agent_to_kill over: killer
			{
				ask agent_to_kill
				{
					write("self " + self);
					do die;
				}
			}
		}
		// we don't need them anymore -> go kill yourselff
	}
	action update_border_cell
	{
		// update the border_cells
		write("update_border_cell " + updated_organics);
		
		loop updated_organic_pair over: updated_organics.pairs
		{
			write("From " + updated_organic_pair.key + " " + updated_organic_pair.value);
			loop updated_organic over: updated_organic_pair.value
			{
				write("updated_organic neigh " + updated_organic.organic_neighbors);
				write("updated_organic index " + updated_organic.cell_index);
				
				ask Thematic.cammisol[0]
 				{	
 					OrganicParticle OrganicToUpdate <- Particle[updated_organic.cell_index].particle as OrganicParticle;
 					write("OrganicToUpdate.N_labile " + OrganicToUpdate.N_labile);
 					write("OrganicToUpdate.C_labile " + OrganicToUpdate.C_labile);
 					write("OrganicToUpdate.P_labile " + OrganicToUpdate.P_labile);
 					write("OrganicToUpdate.C_recalcitrant " + OrganicToUpdate.C_recalcitrant);
 					write("OrganicToUpdate.N_recalcitrant " + OrganicToUpdate.N_recalcitrant);
 					write("OrganicToUpdate.P_recalcitrant " + OrganicToUpdate.P_recalcitrant);
 					
 					OrganicToUpdate.N_labile <- updated_organic.N_labile;
 					OrganicToUpdate.C_labile <- updated_organic.C_labile;
 					OrganicToUpdate.P_labile <- updated_organic.P_labile;
 					OrganicToUpdate.C_recalcitrant <- updated_organic.C_recalcitrant;
 					OrganicToUpdate.N_recalcitrant <- updated_organic.N_recalcitrant;
 					OrganicToUpdate.P_recalcitrant <- updated_organic.P_recalcitrant;
 					
 					write("OrganicToUpdate.N_labile " + OrganicToUpdate.N_labile);
 					write("OrganicToUpdate.C_labile " + OrganicToUpdate.C_labile);
 					write("OrganicToUpdate.P_labile " + OrganicToUpdate.P_labile);
 					write("OrganicToUpdate.C_recalcitrant " + OrganicToUpdate.C_recalcitrant);
 					write("OrganicToUpdate.N_recalcitrant " + OrganicToUpdate.N_recalcitrant);
 					write("OrganicToUpdate.P_recalcitrant " + OrganicToUpdate.P_recalcitrant);
 					
 					Particle[updated_organic.cell_index].particle <- OrganicToUpdate;
 				}
			}
		}
	}
	
	action migrate_nematode
	{
		write("migrate_nematode start " + my_nematodes);
		write("my_pores " + my_pores);
		map<int, list<Nematode>> nematodes_to_migrate;
		loop nematode over: my_nematodes
		{
			write("" + nematode + " index_of_pores[nematode.current_pore.index] " + index_of_pores[nematode.current_pore.index]);
			
			if(!(my_pores contains index_of_pores[nematode.current_pore.index])) // nematode moved to a pore that I don't simulate
			{
				write(" ?? ??? ?? ? ? ? ? ?  " + nematode);
				// migrate the nematode to the right instace
				int cluster <- cell_to_cluster[nematode.current_pore.index];
				
				if(nematodes_to_migrate[cluster] = nil)
				{
					nematodes_to_migrate[cluster] <- list<Nematode>(nematode);
				}else
				{
					nematodes_to_migrate[cluster] << nematode;
				}
			}
		}
		
		write("nematodes_to_migrate " + nematodes_to_migrate);
	}
}

experiment distribution type: MPI_EXP
{

	float N_dom <- 0.0;
	float P_dom <- 0.0;
	float C_dom <- 0.0;
	
	float N_dim <- 0.0;
	float P_dim <- 0.0;
	
	reflex snapshot // take a snapshot of the current distribution model instance
	{
		write("SNAPPING___________________________________ " + cycle);
		int mpi_id <- MPI_RANK;
		
		if(cycle = 0)
		{
			ask simulation
			{	
				if(mpi_id = 0)
				{
					save (snapshot("cammisole")) to: "../output.log/snapshot/cammisol.png" rewrite: true;	
					save (snapshot("global distributed")) to: "../output.log/snapshot/distributed_cammisol_global.png" rewrite: true;
				}
				save (snapshot("distributed")) to: "../output.log/snapshot/distributed_cammisol_" + mpi_id + ".png" rewrite: true;	
					
				save (snapshot("UNSCHEDULED")) to: "../output.log/snapshot/UNSCHEDULED_" + mpi_id + ".png" rewrite: true;
				
				
	 			//save "N_dom; P_dom; C_dom; N_dim; P_dim" to: "../output.log/results/"+mpi_id+"/dam_" + mpi_id + ".csv" format: 'csv' rewrite: true;
	 			//save "Cl; Nl; Pl; Cr; Nr; Pr" to: "../output.log/results/"+mpi_id+"/organics_" + mpi_id + ".csv" format: 'csv' rewrite: true;
	 			//save "O_c; F_c; M_c" to: "../output.log/results/"+mpi_id+"/bacteria_" + mpi_id + ".csv" format: 'csv' rewrite: true;
	 			
	 			save "nematode_CO2_emissions"  to: "../output.log/results/CO2/CO2_" + mpi_id + ".csv" format: 'csv' rewrite: true;
			}
		}
		if(cycle mod 10 = 0)
		{
			ask Thematic.cammisol[0]
			{	
				/*list<PoreParticle> scheduled <- PoreParticle where each.scheduled;
				
	 			save "" + sum(scheduled collect each.dam.dom[0])/#gram +
	 			";" + sum(scheduled collect each.dam.dom[1])/#gram +
	 			";" + sum(scheduled collect each.dam.dom[2])/#gram +
	 			";" + sum(scheduled collect each.dam.dim[0])/#gram + 
	 			";" + sum(scheduled collect each.dam.dim[1])/#gram 
	 			to: "../output.log/results/"+mpi_id+"/dam_" + mpi_id + ".csv" format: 'csv' rewrite: false;
	 			
				save "" + sum(OrganicParticle collect each.C_labile)/#gram + 
				";" +sum(OrganicParticle collect each.N_labile)/#gram + 
				";" + sum(OrganicParticle collect each.P_labile)/#gram + 
				";" + sum(OrganicParticle collect each.C_recalcitrant)/#gram +
				";" +sum(OrganicParticle collect each.N_recalcitrant)/#gram +
				";" +sum(OrganicParticle collect each.P_recalcitrant)/#gram 
				to: "../output.log/results/"+mpi_id+"/organics_" + mpi_id + ".csv" format: 'csv' rewrite: false;
				
				let t <- list<MicrobePopulation>((scheduled collect each.populations));
				
				let O <- scheduled collect each.populations[0].C;
				let F <- scheduled collect each.populations[1].C;
				let M <- scheduled collect each.populations[2].C;
				
				save "" + sum(O)/#gram + 
				";" + sum(F)/#gram + 
				";" + sum(M)/#gram 
				 to: "../output.log/results/"+mpi_id+"/bacteria_" + mpi_id + ".csv" format: 'csv' rewrite: false;*/
				 
				 
	 			save nematode_CO2_emissions  to: "../output.log/results/CO2/CO2_" + mpi_id + ".csv" format: 'csv' rewrite: false;
			}
			
			/*ask simulation
			{
				save (snapshot("dam_distributed")) to: "../output.log/snapshot/dam_distributed/dam_" + cycle + "__" + mpi_id + ".png" rewrite: true;
			}*/	
		}
				
			
		/*}else
		{
		}*/
	}
	output
	{
		display 'global distributed'
		{
			species particle_color aspect: distributed_global;
		}
		display 'distributed'
		{
			species particle_color aspect: distributed;
		}	
		display 'cammisole'
		{
			species particle_color aspect: cammisole;
		}	
		display 'UNSCHEDULED'
		{
			species particle_color aspect: UNSCHEDULED;
		}
		
		/*display "dam_distributed" type: java2D {
			chart "dam" type:series {
				data "N dom (g)" value: N_dom/#gram style:spline marker:false thickness:3;
				data "P dom (g)" value: P_dom/#gram style:spline marker:false thickness:3;
				data "C dom (g)" value: C_dom/#gram style:spline marker:false thickness:3;
				data "N dim (g)" value: N_dim/#gram style:spline marker:false thickness:3;
				data "P dim (g)" value: P_dim/#gram style:spline marker:false thickness:3;
			}
		}*/
	}
	
} 
