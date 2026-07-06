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

	// Cycles between two border-organics synchronizations. Must match cammisol.gaml's
	// results_save_period and compare_co2.py's --step default, since they all assume
	// the same save/sync cadence.
	int sync_period <- 10;

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
 		seed <- 38.0; // seed for this model's own RNG (partitioning); the cammisol submodel fixes its own seed independently, see cammisol.gaml

 		create Communication_Agent_MPI; 					// init of the communication agent
 		MPI_RANK <- Communication_Agent_MPI[0].MPI_RANK;	// get the MPI Rank of this instance
 		MPI_SIZE <- Communication_Agent_MPI[0].MPI_SIZE;	// get the size of the MPI Network

 		create Partitionning_Agent;
 		create Synchronization_Agent;

 		cluster_number <- MPI_SIZE; // number of cluster = MPI_SIZE

 		write("rank " + MPI_RANK + "/" + MPI_SIZE + " starting");

 		create Thematic.cammisol with: [grid_size:grid_size, nematodes_count:nematodes_count, distributed_simulation: true];

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
			clusters <- grid_score_partitioning(grid_size, grid_size, grid_score, cluster_number);
	 		my_cells <- clusters[MPI_RANK];
			do coloring(clusters);
 			my_pores <- pores_by_clusters[MPI_RANK];

 			do color_nematode_cell;
			do border_cell;
 			do unschedule;
	 		do unschedule_color;

			write("rank " + MPI_RANK + ": " + length(my_cells) + " cells, " + length(my_pores) + " pores, " + length(my_nematodes) + " nematodes");
		}
 	}

 	action end_simulations
	{
 		if(check_end_simulation())
 		{
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
	 		if(cycle mod sync_period = 0)
	 		{
	 			do send_my_border_cell;
	 		}
	 		do migrate_nematode;
 		}
 	}

 	action run_thematic_model // run a cycle of the CAMMISOL model
 	{
 		write("" + MPI_RANK + " distribution step : --------------------------------------" + cycle);
 		//write("total_duration " + float(total_duration)/1000 + "s");
 		//write("duration " + float(duration)/1000 + "s");
 		if(check_end_simulation())
 		{
 			return; // CAMMISOL is over we don't need to run step again
 		}
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
	 			write("" + MPI_RANK + "CAMMISOL SIMULATION IS OVER");
				write("" + MPI_RANK + "total_duration " + float(total_duration)/1000 + "s");
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
									add cell.index to: my_border_organics;

									// important use case : the current cell is an organics and is next to a pore on another cluster, we need to update that cell
									do add_my_border_cells_neighbors(cell, neigh_cluster);
								}
								match "organic" {}
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
 			if(!(my_pores contains nematode.value))
 			{
 				nematode.key.scheduled <- false;
 			}
 		}

 		int index_pore_particle <- 0;
 		list<int> index_of_pore_to_unschedule;
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
 			}
 		}
 	}

 	/**
  	* color Nematode position and store them in  nematode_to_cell
  	*/
	action color_nematode_cell
	{
 		ask Thematic.cammisol[0]
 		{
 			let tmp_my_pores <- my_pores;
			ask Nematode
			{
				nematode_to_cell[self] <- index_of_pores[current_pore.index]; // can't directly access particle_color in this context
				if(tmp_my_pores contains index_of_pores[current_pore.index])
				{
					my_nematodes << self;
				}
			}
		}
		loop cell over: nematode_to_cell
		{
			particle_color[cell].nematode <- true;
		}
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

// Transient carrier used to hand a Nematode's state over to the rank that now
// owns the cluster it wandered into (see Synchronization_Agent.migrate_nematode).
// It is only a serializable payload for MPI_ALLTOALL: it never acts on its own
// and is killed as soon as its fields have been copied onto the receiving
// rank's local Nematode with the same nematode_id.
species NematodeMigration
{
	int nematode_id;
	float stomack_C;
	float stomack_N;
	float stomack_P;
	bool awake;
	int dest_cell_index; // grid cell index (Particle.index) of the pore the nematode moved into
}

species Synchronization_Agent
{
	map<int,list<OrganicParticle>> organics_to_send; 	// key : MPI RANK; value : list of OrganicParticle agent to send
	map<int,list<OrganicParticle>> updated_organics;	// key : MPI RANK; value : list of OrganicParticle agent with data
	map<int,list<NematodeMigration>> migrations_to_send;	// key : MPI RANK; value : list of NematodeMigration payloads to send
	map<int,list<NematodeMigration>> migrations_received;	// key : MPI RANK; value : list of NematodeMigration payloads received

	action send_my_border_cell
	{
		// send my_border_cells
		loop cell over: my_border_organics_neighbors.pairs
		{
			ask Thematic.cammisol[0]
 			{
				loop organic over: cell.value
				{
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

		ask Communication_Agent_MPI
		{
	    	myself.updated_organics <- MPI_ALLTOALL(myself.organics_to_send); // MPI all to all
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
			loop agent_to_kill over: killer
			{
				ask agent_to_kill
				{
					do die;
				}
			}
		}
		// we don't need them anymore -> go kill yourself
	}
	action update_border_cell
	{
		// update the border_cells
		loop updated_organic_pair over: updated_organics.pairs
		{
			loop updated_organic over: updated_organic_pair.value
			{
				ask Thematic.cammisol[0]
 				{
 					OrganicParticle OrganicToUpdate <- Particle[updated_organic.cell_index].particle as OrganicParticle;

 					OrganicToUpdate.N_labile <- updated_organic.N_labile;
 					OrganicToUpdate.C_labile <- updated_organic.C_labile;
 					OrganicToUpdate.P_labile <- updated_organic.P_labile;
 					OrganicToUpdate.C_recalcitrant <- updated_organic.C_recalcitrant;
 					OrganicToUpdate.N_recalcitrant <- updated_organic.N_recalcitrant;
 					OrganicToUpdate.P_recalcitrant <- updated_organic.P_recalcitrant;

 					Particle[updated_organic.cell_index].particle <- OrganicToUpdate;
 				}
			}
		}
	}

	action migrate_nematode
	{
		list<Nematode> nematodes_that_left;

		loop nematode over: my_nematodes
		{
			int dest_cell_index <- index_of_pores[nematode.current_pore.index];

			if(!(my_pores contains dest_cell_index)) // nematode moved to a pore that I don't simulate
			{
				int dest_cluster <- cell_to_cluster[dest_cell_index]; // cluster id == owning MPI rank
				write("rank " + MPI_RANK + ": nematode " + nematode + " leaving my cluster for cluster " + dest_cluster);

				create NematodeMigration with: [
					nematode_id: nematode.nematode_id,
					stomack_C: nematode.stomack_C,
					stomack_N: nematode.stomack_N,
					stomack_P: nematode.stomack_P,
					awake: nematode.awake,
					dest_cell_index: dest_cell_index
				] returns: created_migration;

				if(migrations_to_send[dest_cluster] = nil)
				{
					migrations_to_send[dest_cluster] <- list<NematodeMigration>(created_migration[0]);
				}else
				{
					migrations_to_send[dest_cluster] << created_migration[0];
				}
				nematodes_that_left << nematode;
			}
		}

		// I no longer simulate the nematodes that left my cluster: freeze them here
		// (their state is now carried by the migration payload) until/unless they
		// wander back into my cluster later, at which point they'll be re-adopted
		// below just like any other incoming nematode.
		loop nematode over: nematodes_that_left
		{
			remove nematode from: my_nematodes;
			nematode.scheduled <- false;
		}

		ask Communication_Agent_MPI
		{
			myself.migrations_received <- MPI_ALLTOALL(myself.migrations_to_send);
		}

		// Adopt every nematode that just moved into my cluster: find my own local
		// copy of that nematode_id (every rank already has one, created at init
		// with the same seed) and bring it up to date, then start scheduling it.
		loop incoming_pair over: migrations_received.pairs
		{
			loop migration over: incoming_pair.value
			{
				ask Thematic.cammisol[0]
				{
					Nematode target <- Nematode first_with (each.nematode_id = migration.nematode_id);
					if(target != nil)
					{
						target.current_pore <- Particle[migration.dest_cell_index].particle as PoreParticle;
						target.location <- any_location_in(target.current_pore);
						target.stomack_C <- migration.stomack_C;
						target.stomack_N <- migration.stomack_N;
						target.stomack_P <- migration.stomack_P;
						target.awake <- migration.awake;
						target.scheduled <- true;
						nematode_to_cell[target] <- migration.dest_cell_index;
						my_nematodes << target;
					}
				}
				// the payload was only a carrier for the state above, discard it
				ask migration
				{
					do die;
				}
			}
		}

		migrations_to_send <- nil;
		migrations_received <- nil;
	}
}

experiment distribution type: MPI_EXP
{
	reflex snapshot // take a snapshot of the current distribution model instance
	{
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

	 			save "nematode_CO2_emissions"  to: "../output.log/results/CO2/CO2_" + mpi_id + ".csv" format: 'csv' rewrite: true;
			}
		}
		if(cycle mod sync_period = 0)
		{
			ask Thematic.cammisol[0]
			{
	 			save nematode_CO2_emissions  to: "../output.log/results/CO2/CO2_" + mpi_id + ".csv" format: 'csv' rewrite: false;
			}
		}
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
	}

}
