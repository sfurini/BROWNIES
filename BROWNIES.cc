/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Copyright (C) 2014 Claudio Berti
-----------------------------------------------------------------------------
Claudio Berti, Ph.D.,
Department of Molecular Biophysics and Physiology,
Rush University Medical Center
and
University of Bologna, Bologna, Italy

Tel: 312.942.6756
Fax: 312.942.8711
E-Mail: brt.cld@gmail.com
-----------------------------------------------------------------------------		

This file is part of BROWNIES.

BROWNIES is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Private, research, and institutional use is free under the condition of
acknowledge the author and contributors citing the following papers:

- Berti, C.; Furini, S.; Gillespie, PACO: PArticle COunting Method To Enforce
Concentrations in Dynamic Simulations. J. Chem. Theory Comput. (in press)

- Berti, C.; Furini, S.; Gillespie, D.; Boda, D.; Eisenberg, R. S.; 
Sangiorgi, E; Fiegna, A 3-D Brownian dynamics simulator for the study of ion 
permeation through membrane pores J. Chem. Theory Comput. 2014, 10, 2911-2926

- Berti, C.; Gillespie, D.; Bardhan, J. P.; Eisenberg, R. S.; Fiegna, C.,
Comparison of three-dimensional Poisson solution methods for particle-based 
simulation and inhomogeneous dielectrics. Phys. Rev. E 2012, 86, 011912

- Berti, C.; Gillespie, D.; Eisenberg, R.; Fiegna, C. Particle-based 
simulation of charge transport in discrete-charge nano-scale systems: the 
electrostatic problem. Nanoscale Research Letters 2012, 7, 135



You may distribute modified versions of this code UNDER THE CONDITION THAT
THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER 
COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE 
FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE 
MODIFICATIONS.  Distribution of this code as part of a commercial system
is permissible ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR. (If you are 
not directly supplying this code to a customer, and you are instead 
telling them how they can obtain it for free, then you are not required
to make any arrangement with me.) 

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <sstream> 
#include <cstring>
#include <algorithm>
#include <dirent.h>
#include <map>
#include <cstdio>
#include <cctype>
#include <cmath>
#include <ctime>
#include <csignal>
#include <limits>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>
#include <boost/random/linear_congruential.hpp> 
#include <boost/random/lagged_fibonacci.hpp>	
#include <boost/random/uniform_int.hpp>		
#include <boost/random/uniform_real.hpp>	
#include <boost/random/normal_distribution.hpp>	
#include <boost/random/variate_generator.hpp>	

#include "constants.h"
#include "ions_properties.h"
#include "utils.h"
#include "classes.h"
#include "physics_functions.h"
#include "file_functions.h"
#include "input_output.h"
#include "sim_domain.h"
#include "ions_functions.h"
#include "PACO.h"
#include "icc.h"
#include "retrace.h"
#include "sim_structures.h"

using namespace std;

int main(int argc,char *argv[]){
	
	srand(time(NULL));
	
//###############################################
//		PARSE AND CHECK ARGUMENTS
//###############################################	
	vector<string> arguments;
	for(int i=0; i<argc; i++){
		arguments.push_back(string(argv[i]));
	}
	if(arguments.size()!=3){
		cerr << "USAGE:"<<endl;
		cerr << "1) ./BROWNIES -h\t\t\t\t\tfor help"<<endl;
		cerr << "2) ./BROWNIES -i path to .conf file\t\t\tto run a simulation"<<endl;
		exit(1);
	}
	if(arguments.at(1).compare("-i")!=0){
		cerr << "USAGE:"<<endl;
		cerr << "1) ./BROWNIES -h\t\t\t\t\tfor help"<<endl;
		cerr << "2) ./BROWNIES -i path to .conf file\t\t\tto run a simulation"<<endl;
		exit(1);
	}
	retrieve_parameters(arguments.at(2));
	check_parameters();
	
//###############################################
//		CREATE SIMULATION DOMAIN
//###############################################	
	create_simulation_domain();
	create_EPFFs();
	if(PRM.SIM_TYPE.compare("BULK")){
		create_surfaces();
		MATRIX_A.clear();
		adjust_surfaces();
		compute_SR_dielectric_boundary();
		PRM.NUM_OF_SUB_DIV=subSurfaces.at(0).size();
	}
	
//###############################################
//		SURFACES HANDLING
//###############################################	
	NUM_OF_SURFACES=surfaces.size();
	if(NUM_OF_SURFACES>1200){
		cout << "THE MAXIMUM NUMBER OF SURFACES IS 1200. Change the .conf file";
		exit(1000);
	}
	for(int i=0; i<NUM_OF_SURFACES; i++){
		SURFACES[i].center_x=surfaces.at(i).center.x;
		SURFACES[i].center_y=surfaces.at(i).center.y;
		SURFACES[i].center_z=surfaces.at(i).center.z;
		SURFACES[i].normal[0]=surfaces.at(i).normal[0];
		SURFACES[i].normal[1]=surfaces.at(i).normal[1];
		SURFACES[i].normal[2]=surfaces.at(i).normal[2];
		SURFACES[i].area=surfaces.at(i).area;
	}
	posix_memalign((void **)(&vector_c), 128, sizeof(double)*surfaces.size());
	posix_memalign((void **)(&vector_h), 128, sizeof(double)*surfaces.size());
	posix_memalign((void **)(&vector_c_MD_map), 128, sizeof(double)*surfaces.size());
	posix_memalign((void **)(&vector_h_MD_map), 128, sizeof(double)*surfaces.size());
	
//###############################################
//		COMPUTE C_BASE
//###############################################	
	compute_c_base();		
	print_vmd_channel_files();
	if(PRM.SIM_TYPE.compare("PORE")==0 && PRM.DIFF_COEFF_IN_CHANNEL!=1){
		compute_diffusion_coefficients();
	}

//###############################################
//		OUTPUT PARAMETERS
//###############################################	
	cout << PRM;
	
//###############################################
//		INITIALIZE CONTROL CELLS
//###############################################	
	initialize_left_control_cell();
	initialize_right_control_cell();
	Control_cell left_cell0(left_cell.back());
	Control_cell right_cell0(right_cell.back());
	
//###############################################
//		INITIALIZE ION BOXES
//###############################################	
	create_ion_boxes();
	create_SR_ion_boxes();
	
//###############################################
//		PREPARE STATISTICS STRUCTURES 
//###############################################
	STAT.reset_statistics();
	
//###############################################
//		RE-ADJUST EXTERNAL ELECTRIC FIELD 
//###############################################	
	compute_potential_on_trajectories();
	
//###############################################
//		MAKE FIRST STEP
//###############################################
	int MAX_NUMBER_OF_IONS_PER_STEP=estimate_max_number_of_ions();
	cout << "MAX_NUMBER_OF_IONS_PER_STEP:\t" << MAX_NUMBER_OF_IONS_PER_STEP << endl;
	for(int i=0; i<PRM.HISTORY_SIZE; i++){
		STEPS[i]=0;
		NUM_OF_IONS_IN_STEP[i]=0;
	}
	LEFT_CELLS[0]=left_cell.back();
	RIGHT_CELLS[0]=right_cell.back();
	do_initial_step();
	//print initial ion configuration
	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){
		cout << ion_index<< " " << IONS[INDEX_LAST_STEP][ion_index].name<< " "<< IONS[INDEX_LAST_STEP][ion_index].kind<< " "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].x <<" "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].y<<" "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].z<<endl;
	}
	
//###############################################
//		SIMULATION
//###############################################	
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	cout << "SIMULATION STARTS!\t-\t" << asctime (timeinfo) <<"\t"<<flush;
	retrace_window=1;
	last_step_statistics=0;
	last_step_output=-1;
	step_cap=-10000000;
	local_step_cap1=-10000000;
	local_step_cap2=-10000000;
	local_step_cap3=-10000000;
	local_step_cap4=-10000000;
	local_step_cap5=-10000000;
	attempts=0;		
	double MAX_STATS_OUT_EVERY=5.1e-9;		// data every 51 ns at most
	PRM.STATS_OUT_FREQ=1;
	int new_freq_out=1;
	int den_freq_out=1;
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		sample_ions[is]=Ion(is);
	}
	if(PRM.SIM_TYPE.compare("PORE")==0){		
		createFile(PRM.PREFIX+".filter.pdb");
	}
	
//=================================================================//		
//#################################################################//		
//####		SIMULATION STARTS HERE! 		       ####//
//#################################################################//
//=================================================================//
	cout << flush << endl;
	INDEX_LAST_STEP=(-1+50)%50;
	INDEX_STAT_STEP=(2)%50;
	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[0]; ion_index++){
		IONS[INDEX_LAST_STEP][ion_index]=IONS[0][ion_index];
		IONS[INDEX_STAT_STEP][ion_index]=IONS[0][ion_index];
	}
	NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=NUM_OF_IONS_IN_STEP[0];
	NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]=NUM_OF_IONS_IN_STEP[0];
	LEFT_CELLS[INDEX_LAST_STEP]=LEFT_CELLS[0];
	LEFT_CELLS[INDEX_STAT_STEP]=LEFT_CELLS[0];
	RIGHT_CELLS[INDEX_LAST_STEP]=RIGHT_CELLS[0];
	RIGHT_CELLS[INDEX_STAT_STEP]=RIGHT_CELLS[0];
	cout << "\nINDEX_LAST_STEP:\t" << INDEX_LAST_STEP << endl;
	cout << "INDEX_STAT_STEP:\t" << INDEX_STAT_STEP << endl;
	STEPS[INDEX_LAST_STEP]=PRM.FIRST_STEP-1;
	STEPS[INDEX_STAT_STEP]=PRM.FIRST_STEP-1;
	indexesGap=0;
	  
	for(step=PRM.FIRST_STEP; step<=PRM.SIM_STEPS; step++){

		if ((PRM.trajectory != 0) && (step > PRM.PREP_STEPS)) {
			print_ions_trajectories();
		}
  
//======================================================================================
//		01		checks and output handling	
//======================================================================================			
		if(step==1){
			attempts=0.00;
			STAT.ERROR_nan=0.00;
			STAT.ERROR_long_jump=0.00;
			STAT.ERROR_into_membrane=0.00;
			STAT.ERROR_ion_boxes=0.00;
		}
		attempts=attempts+double(1.00);	
//======================================================================================
//		02		duplicate last configuration	
//======================================================================================		
		duplicate_last_configuration();
//======================================================================================
//		03		langevin motion	
//======================================================================================		
		if(!PRM.CONSTANT_DIFF_COEFF){
			for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
				if(IONS[INDEX_LAST_STEP][ion_index].name.substr(0, 1).compare("J")){
					int index_z=IONS[INDEX_LAST_STEP][ion_index].z-PRM.MIN_Z;
					IONS[INDEX_LAST_STEP][ion_index].diffusion_coeff=diffusion_coefficients.at(IONS[INDEX_LAST_STEP][ion_index].kind).at(index_z);
					IONS[INDEX_LAST_STEP][ion_index].set_brownian_parameters();
				}
			}
		}		
		langevin_compute_new_position_and_velocity();
//======================================================================================
//		04		boundary conditions for ion concentrations	
//======================================================================================		
		PACO_update_left_CC();
		PACO_update_right_CC();
//======================================================================================
//		05		manage simulation history
//======================================================================================		
		bool CORRECT=true;
		CORRECT=motion_check_correctness();
		if(CORRECT){
			motion_prepare_next_step();
			if(INDEX_LAST_STEP==(INDEX_STAT_STEP+50-2)%50){
				new_freq_out=STEPS[INDEX_STAT_STEP]/den_freq_out;
				if(new_freq_out==10){
					den_freq_out=den_freq_out*10;
					if(PRM.DELTA_T*den_freq_out<=MAX_STATS_OUT_EVERY){
						PRM.STATS_OUT_FREQ=den_freq_out;
					}
				}
				if(STEPS[INDEX_STAT_STEP]>last_step_statistics){
					STAT.update_statistics();
					last_step_statistics=STEPS[INDEX_STAT_STEP];	
				}
				if(STEPS[INDEX_STAT_STEP]>=0 && STEPS[INDEX_STAT_STEP]>last_step_output && fmod(STEPS[INDEX_STAT_STEP],PRM.STATS_OUT_FREQ)==0){
					STAT.print_statistics();
					time_t rawtime;
					struct tm * timeinfo;
					time ( &rawtime );
					timeinfo = localtime ( &rawtime );
					cout <<"STEP " <<  STEPS[INDEX_STAT_STEP] << " DONE!\t-\t" << asctime (timeinfo) <<flush;
					double sim_time=STEPS[INDEX_STAT_STEP]*PRM.DELTA_T;			
					cout << "Simulated " << getNumberWithScale(sim_time) << "s"<<endl;
					cout << "Efficiency: " << double(step)/(attempts) <<endl <<endl;
					last_step_output=STEPS[INDEX_STAT_STEP];
				}
				INDEX_STAT_STEP=(INDEX_STAT_STEP+1)%50;
			}
		}
		else{
			if(step==step_cap){
				retrace_window++;
			}
			else if(step>step_cap){
				retrace_window=1;
				step_cap=step;
			}
			else{
				if(step==local_step_cap1){
					retrace_window++;
				}
				else if(step>local_step_cap1){
					retrace_window=1;
					local_step_cap1=step;
				}
				else{
					if(step==local_step_cap2){
						retrace_window++;
					}
					else if(step>local_step_cap2){
						retrace_window=1;
						local_step_cap2=step;
					}
					else{
						if(step==local_step_cap3){
							retrace_window++;
						}
						else if(step>local_step_cap3){
							retrace_window=1;
							local_step_cap3=step;
						}
						else{
							if(step==local_step_cap4){
								retrace_window++;
							}
							else if(step>local_step_cap4){
								retrace_window=1;
								local_step_cap4=step;
							}
							else{
								if(step==local_step_cap5){
									retrace_window++;
								}
								else if(step>local_step_cap5){
									retrace_window=1;
									local_step_cap5=step;
								}
							}
						}
					}
				}
			}			
			bool restart_from_the_beginning=false;
			for(int iii=0; iii<retrace_window; iii++){
				step=step-1;
				INDEX_LAST_STEP=(INDEX_LAST_STEP-1+50)%50;
				if(step<PRM.FIRST_STEP-1 || INDEX_LAST_STEP==INDEX_STAT_STEP){
					restart_from_the_beginning=true;
				}
			}
			if(restart_from_the_beginning){
				cout << "RESTART FROM THE BEGINNING!!" << endl;
				INDEX_LAST_STEP=(-1+50)%50;
				INDEX_STAT_STEP=(2)%50;
				STEPS[INDEX_LAST_STEP]=PRM.FIRST_STEP-1;
				STEPS[INDEX_STAT_STEP]=PRM.FIRST_STEP-1;
				LEFT_CELLS[INDEX_LAST_STEP]=left_cell0;
				RIGHT_CELLS[INDEX_LAST_STEP]=right_cell0;
				do_initial_step();
				STAT.reset_after_restarting_from_the_beginning();
				retrace_window=1;
				last_step_statistics=0;
				last_step_output=-1;
				step_cap=-10000000;
				local_step_cap1=-10000000;
				local_step_cap2=-10000000;
				local_step_cap3=-10000000;
				local_step_cap4=-10000000;
				local_step_cap5=-10000000;
				attempts=0;	
				step=PRM.FIRST_STEP-1;	
				retrace_window=1;
				new_freq_out=1;
				den_freq_out=1;
				cout << "let's see.." << endl;
				string yyy; getline(cin, yyy);
			}
		}
	}
//=================================================================//		
//#################################################################//		
//####		     SIMULATION STOPS HERE!	               ####//
//#################################################################//
//=================================================================//		
	cout << "simulation ended"<<endl;
	if(!surfaces.empty()){
		free(STAT.vector_h_RAMO);
	}	
	return 0;
}
