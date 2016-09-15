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

#include "retrace.h"
#include "constants.h"
#include "utils.h"
#include "classes.h"
#include "file_functions.h"
#include "input_output.h"
#include "sim_domain.h"
#include "ions_functions.h"
#include "sim_structures.h"


void langevin_step_0(){

	

	cout << "langevin_step_0..............."<<flush;
	
	compute_force_on_ions();
	
	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){
		initialize_ion_velocity(IONS[INDEX_LAST_STEP][ion_index]);
		
		IONS[INDEX_LAST_STEP][ion_index].X_n[0]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_X;
		IONS[INDEX_LAST_STEP][ion_index].X_n[1]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_X;
		IONS[INDEX_LAST_STEP][ion_index].X_n[2]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_X;
		
		IONS[INDEX_LAST_STEP][ion_index].x_next=IONS[INDEX_LAST_STEP][ion_index].x+IONS[INDEX_LAST_STEP][ion_index].velocity[0]*IONS[INDEX_LAST_STEP][ion_index].K_step0_1+IONS[INDEX_LAST_STEP][ion_index].force[0]*IONS[INDEX_LAST_STEP][ion_index].K_step0_2+IONS[INDEX_LAST_STEP][ion_index].X_n[0];
		IONS[INDEX_LAST_STEP][ion_index].y_next=IONS[INDEX_LAST_STEP][ion_index].y+IONS[INDEX_LAST_STEP][ion_index].velocity[1]*IONS[INDEX_LAST_STEP][ion_index].K_step0_1+IONS[INDEX_LAST_STEP][ion_index].force[1]*IONS[INDEX_LAST_STEP][ion_index].K_step0_2+IONS[INDEX_LAST_STEP][ion_index].X_n[1];
		IONS[INDEX_LAST_STEP][ion_index].z_next=IONS[INDEX_LAST_STEP][ion_index].z+IONS[INDEX_LAST_STEP][ion_index].velocity[2]*IONS[INDEX_LAST_STEP][ion_index].K_step0_1+IONS[INDEX_LAST_STEP][ion_index].force[2]*IONS[INDEX_LAST_STEP][ion_index].K_step0_2+IONS[INDEX_LAST_STEP][ion_index].X_n[2];
		
		IONS[INDEX_LAST_STEP][ion_index].x_prev=IONS[INDEX_LAST_STEP][ion_index].x;
		IONS[INDEX_LAST_STEP][ion_index].y_prev=IONS[INDEX_LAST_STEP][ion_index].y;
		IONS[INDEX_LAST_STEP][ion_index].z_prev=IONS[INDEX_LAST_STEP][ion_index].z;
		
		
	}
	
	
	cout << "OK!\n"<<flush;

	
	return;
}

bool motion_check_step_0(){
	
	
	
	cout << "motion_check_step_0..............."<<flush;
	
// check flight length	
	double MAX_VEL=1.5e4;	// 2e3
	double MAX_FLIGHT=MAX_VEL*PRM.DELTA_T;
	
	bool CORRECT=true;
	
	vector <int> long_flights;
	long_flights.clear();
	
	int jumps=0;
	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){		
		double flight=get_flight_length(IONS[INDEX_LAST_STEP][ion_index]);
		double velocity=get_ion_velocity(IONS[INDEX_LAST_STEP][ion_index]);
		if(flight>MAX_FLIGHT || velocity>MAX_VEL){
			
			PRM.long_flight++;
			long_flights.push_back(1);
		}
		else{
			long_flights.push_back(0);
		}
	}
	cout << "a" <<endl;
	for(int ion_index_1=0; ion_index_1<long_flights.size(); ion_index_1++){
		if(long_flights.at(ion_index_1)==1){
			
			cout << "long flight..."<<flush;
			reinitialize_ion_first_step(ion_index_1);
			CORRECT=false;
			cout << "OK!\n" << flush;
		}
	}
	cout << "b" <<endl;
	if(PRM.SIM_TYPE.compare("BULK")==0){
		
	}
	else{
	
//membrane check...
		for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].z_next>=PRM.Z_MOUTH_LEFT && 1e12*IONS[INDEX_LAST_STEP][ion_index].z_next<=PRM.Z_MOUTH_RIGHT){
				int ind_on_lim=1e12*IONS[INDEX_LAST_STEP][ion_index].z_next-PRM.Z_MOUTH_LEFT;
				double dist_from_axis=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index].x_next, IONS[INDEX_LAST_STEP][ion_index].y_next, 0.00, 0.00, 0.00, 0.00);	
				
				if(dist_from_axis>=limits.at(ind_on_lim)-IONS[INDEX_LAST_STEP][ion_index].pm_radius){
					cout << "ion into the membrane..."<<flush;
					reinitialize_ion_first_step(ion_index);
					CORRECT=false;
					cout << "OK!\n" << flush;
				}
			}
		}
		cout << "c" <<endl;
//ion boxes error check..		
		for(int ib=0; ib<ion_boxes.size(); ib++){
			for(int it=0; it<ion_boxes.at(ib).ion_indexes.size(); it++){
				if(1e12*IONS[INDEX_LAST_STEP][ion_boxes.at(ib).ion_indexes.at(it)].z_next<ion_boxes.at(ib).MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_boxes.at(ib).ion_indexes.at(it)].z_next>ion_boxes.at(ib).MAX_Z){
					
					cout << 1e12*IONS[INDEX_LAST_STEP][ion_boxes.at(ib).ion_indexes.at(it)].z_next << "\t" << ion_boxes.at(ib).MIN_Z << "\t" << ion_boxes.at(ib).MAX_Z <<endl;
					cout << "ion out of ion_box "<<ion_boxes.at(ib).index<<"..."<<flush;
					
					cout << "indice ione e': " << ion_boxes.at(ib).ion_indexes.at(it) << endl;
				
					reinitialize_ion_first_step(ion_boxes.at(ib).ion_indexes.at(it));
					CORRECT=false;
					cout << "OK!\n" << flush;
					string jjj;getline(cin, jjj);
				}
			}
		}
	}
cout << "d" <<endl;
	return CORRECT;	

}

void reinitialize_ion_first_step(int ion_index_1){
	
	
	
	bool found=false;
	for(int ib=0; ib<ion_boxes.size(); ib++){
		for(int it=0; it<ion_boxes.at(ib).ion_indexes.size(); it++){
			
			if(ion_boxes.at(ib).ion_indexes.at(it)==ion_index_1){
				found=true;
				
				double min_z=ion_boxes.at(ib).MIN_Z+double(1.00)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius;
				double width_z=(ion_boxes.at(ib).MAX_Z-ion_boxes.at(ib).MIN_Z)-double(2.00)*double(1.00)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius;
				
				bool goon=true;
				do{
					goon=true;

					IONS[INDEX_LAST_STEP][ion_index_1].z=min_z+width_z*((double)rand()/((double)(RAND_MAX)+(double)(1)));
				
					double R=PRM.SIM_DOMAIN_WIDTH_X/double(2.00)-double(1.00)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius;
					if(IONS[INDEX_LAST_STEP][ion_index_1].z>PRM.Z_MOUTH_LEFT && IONS[INDEX_LAST_STEP][ion_index_1].z<PRM.Z_MOUTH_RIGHT){
						int ind_on_lim=IONS[INDEX_LAST_STEP][ion_index_1].z-PRM.Z_MOUTH_LEFT;
						R=limits.at(ind_on_lim)-double(1.00)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius;
					}
					double r=(double)rand()/((float)RAND_MAX/R);
					double theta=(double)rand()/((double)RAND_MAX/(double(2.00)*M_PI));
					IONS[INDEX_LAST_STEP][ion_index_1].x=r*cos(theta);
					IONS[INDEX_LAST_STEP][ion_index_1].y=r*sin(theta);

					for(int j=0; j<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; j++){
						if(ion_index_1!=j){
							goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, double(1.00)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, double(1.00)*IONS[INDEX_LAST_STEP][j].pm_radius);
						}
					}
				}
				while(!goon);
				
				IONS[INDEX_LAST_STEP][ion_index_1].x=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].x;
				IONS[INDEX_LAST_STEP][ion_index_1].y=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].y;
				IONS[INDEX_LAST_STEP][ion_index_1].z=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].z;
			}
		}
	}
	
	if(!found){
		
		bool goon=true;
		do{
			goon=true;

			IONS[INDEX_LAST_STEP][ion_index_1].x=PRM.MIN_X+PRM.SIM_DOMAIN_WIDTH_X*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][ion_index_1].y=PRM.MIN_Y+PRM.SIM_DOMAIN_WIDTH_Y*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][ion_index_1].z=PRM.LEFT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			if(IONS[INDEX_LAST_STEP][ion_index_1].side==1){
				IONS[INDEX_LAST_STEP][ion_index_1].z=PRM.RIGHT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			}

			for(int j=0; j<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; j++){
				if(ion_index_1!=j){
					goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, double(1.00)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, double(1.00)*IONS[INDEX_LAST_STEP][j].pm_radius);
				}
			}
		}
		while(!goon);

		IONS[INDEX_LAST_STEP][ion_index_1].x=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].x;
		IONS[INDEX_LAST_STEP][ion_index_1].y=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].y;
		IONS[INDEX_LAST_STEP][ion_index_1].z=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].z;
	}
	
	return;
}

bool motion_check_correctness_step_0(){
	
	
	
// check flight length	
	double MAX_VEL=1.5e4;	// 2e3
	double MAX_FLIGHT=MAX_VEL*PRM.DELTA_T;
	
	
	bool CORRECT=true;
	
	vector <int> long_flights;
	long_flights.clear();
	
	int jumps=0;
	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){		
		double flight=get_flight_length(IONS[INDEX_LAST_STEP][ion_index]);
		double velocity=get_ion_velocity(IONS[INDEX_LAST_STEP][ion_index]);
		if(flight>MAX_FLIGHT || velocity>MAX_VEL){
			
			PRM.long_flight++;
			long_flights.push_back(1);
		}
		else{
			long_flights.push_back(0);
		}
	}
	
	for(int ion_index_1=0; ion_index_1<long_flights.size(); ion_index_1++){
		if(long_flights.at(ion_index_1)==1){
			
			cout << "long flight..."<<flush;
			brownian_reinitialize_ion(ion_index_1);
			CORRECT=false;
			cout << "OK!\n"<<flush;
		}
	}
	
	if(PRM.SIM_TYPE.compare("BULK")==0){
		
	}
	else{
	
//membrane check...
		for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].z_next>=PRM.Z_MOUTH_LEFT && 1e12*IONS[INDEX_LAST_STEP][ion_index].z_next<=PRM.Z_MOUTH_RIGHT){
				int ind_on_lim=1e12*IONS[INDEX_LAST_STEP][ion_index].z_next-PRM.Z_MOUTH_LEFT;
				double dist_from_axis=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index].x_next, IONS[INDEX_LAST_STEP][ion_index].y_next, 0.00, 0.00, 0.00, 0.00);	
				
				if(dist_from_axis>=limits.at(ind_on_lim)-IONS[INDEX_LAST_STEP][ion_index].pm_radius){
					cout << "ion into the membrane..."<<flush;
					brownian_reinitialize_ion(ion_index);
					CORRECT=false;
					cout << "OK!\n" << flush;
				}
			}
		}
		
//ion boxes error check..		
		for(int ib=0; ib<ion_boxes.size(); ib++){
			for(int it=0; it<ion_boxes.at(ib).ion_indexes.size(); it++){
				if(1e12*IONS[INDEX_LAST_STEP][ion_boxes.at(ib).ion_indexes.at(it)].z_next<ion_boxes.at(ib).MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_boxes.at(ib).ion_indexes.at(it)].z_next>ion_boxes.at(ib).MAX_Z){
					cout << "ion out of ion_box "<<ion_boxes.at(ib).index<<"..."<<flush;
					brownian_reinitialize_ion(ion_boxes.at(ib).ion_indexes.at(it));
					CORRECT=false;
					cout << "OK!\n" << flush;
				}
			}
		}
	}

	return CORRECT;	
}

bool motion_check_correctness(){
	
	bool CORRECT=true;
	
	//nan check....	
	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
		if(!isnormal(IONS[INDEX_LAST_STEP][ion_index].x_next) || !isnormal(IONS[INDEX_LAST_STEP][ion_index].y_next) || !isnormal(IONS[INDEX_LAST_STEP][ion_index].z_next)){
			CORRECT=false;
			STAT.ERROR_nan+=1.00;
			
			return CORRECT;
		}
	}
	
	//jump check....
	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){
		double flight=get_flight_length(IONS[INDEX_LAST_STEP][ion_index]);
		double velocity=get_ion_velocity(IONS[INDEX_LAST_STEP][ion_index]);
		if(flight>PRM.MAX_FLIGHT || velocity>PRM.MAX_VEL){
			CORRECT=false;
			STAT.ERROR_long_jump+=1.00;
			
			return CORRECT;
		}
		else{
			
		}
	}	
	if(PRM.SIM_TYPE.compare("BULK")==0){
	
		//out of simulation domain check....
		for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].x_next<double(1.5)*PRM.MIN_X || 1e12*IONS[INDEX_LAST_STEP][ion_index].x_next>double(1.5)*PRM.MAX_X || 1e12*IONS[INDEX_LAST_STEP][ion_index].y_next<double(1.5)*PRM.MIN_Y || 1e12*IONS[INDEX_LAST_STEP][ion_index].y_next>double(1.5)*PRM.MAX_Y || 1e12*IONS[INDEX_LAST_STEP][ion_index].z_next<double(1.5)*PRM.MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_index].z_next>double(1.5)*PRM.MAX_Z){
				CORRECT=false;
				STAT.ERROR_into_membrane+=1.00;
			
				return CORRECT;
			}
		}
	}
	else{
		//membrane check...
		for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].z_next>=PRM.Z_MOUTH_LEFT && 1e12*IONS[INDEX_LAST_STEP][ion_index].z_next<=PRM.Z_MOUTH_RIGHT){
				int ind_on_lim=1e12*IONS[INDEX_LAST_STEP][ion_index].z_next-PRM.Z_MOUTH_LEFT;
				double dist_from_axis=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index].x_next, IONS[INDEX_LAST_STEP][ion_index].y_next, 0.00, 0.00, 0.00, 0.00);	
				
				if(dist_from_axis>=limits.at(ind_on_lim)-IONS[INDEX_LAST_STEP][ion_index].pm_radius){
					CORRECT=false;
					STAT.ERROR_into_membrane+=1.00;
					
					return CORRECT;
				}
			}
		}
			
		//ion boxes error check..		
		for(int ib=0; ib<ion_boxes.size(); ib++){
			for(int it=0; it<ion_boxes.at(ib).ion_indexes.size(); it++){
				if(1e12*IONS[INDEX_LAST_STEP][ion_boxes.at(ib).ion_indexes.at(it)].z_next<ion_boxes.at(ib).MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_boxes.at(ib).ion_indexes.at(it)].z_next>ion_boxes.at(ib).MAX_Z){
					CORRECT=false;
					STAT.ERROR_ion_boxes+=1.00;
				
					return CORRECT;
				}

			}
		}
	}
	
	

	return CORRECT;	
}

void brownian_reinitialize_ion(int ion_index_1){
	
	
	
	bool found=false;
	for(int ib=0; ib<ion_boxes.size(); ib++){
		for(int it=0; it<ion_boxes.at(ib).ion_indexes.size(); it++){
			
			if(ion_boxes.at(ib).ion_indexes.at(it)==ion_index_1){
				found=true;
				
				double min_z=ion_boxes.at(ib).MIN_Z+double(1.25)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius;
				double width_z=(ion_boxes.at(ib).MAX_Z-ion_boxes.at(ib).MIN_Z)-double(2.00)*double(1.25)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius;
				
				bool goon=true;
				do{
					goon=true;

					IONS[INDEX_LAST_STEP][ion_index_1].z=min_z+width_z*((double)rand()/((double)(RAND_MAX)+(double)(1)));
				
					double R=PRM.SIM_DOMAIN_WIDTH_X/double(2.00)-double(1.25)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius;
					if(IONS[INDEX_LAST_STEP][ion_index_1].z>PRM.Z_MOUTH_LEFT && IONS[INDEX_LAST_STEP][ion_index_1].z<PRM.Z_MOUTH_RIGHT){
						int ind_on_lim=IONS[INDEX_LAST_STEP][ion_index_1].z-PRM.Z_MOUTH_LEFT;
						R=limits.at(ind_on_lim)-double(1.25)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius;
					}
					double r=(double)rand()/((float)RAND_MAX/R);
					double theta=(double)rand()/((double)RAND_MAX/(double(2.00)*M_PI));
					IONS[INDEX_LAST_STEP][ion_index_1].x=r*cos(theta);
					IONS[INDEX_LAST_STEP][ion_index_1].y=r*sin(theta);

					for(int j=0; j<ions.back().size(); j++){
						if(ion_index_1!=j){
							goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, double(1.1)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, double(1.1)*IONS[INDEX_LAST_STEP][j].pm_radius);
						}
					}
				}
				while(!goon);
				
				IONS[INDEX_LAST_STEP][ion_index_1].x=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].x;
				IONS[INDEX_LAST_STEP][ion_index_1].y=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].y;
				IONS[INDEX_LAST_STEP][ion_index_1].z=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].z;
			}
		}
	}
	if(!found){
		bool goon=true;
		do{
			goon=true;

			IONS[INDEX_LAST_STEP][ion_index_1].x=PRM.MIN_X+PRM.SIM_DOMAIN_WIDTH_X*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][ion_index_1].y=PRM.MIN_Y+PRM.SIM_DOMAIN_WIDTH_Y*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][ion_index_1].z=PRM.LEFT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			if(IONS[INDEX_LAST_STEP][ion_index_1].side==1){
				IONS[INDEX_LAST_STEP][ion_index_1].z=PRM.RIGHT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			}

			for(int j=0; j<ions.back().size(); j++){
				if(ion_index_1!=j){
					goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, double(1.25)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, double(1.25)*IONS[INDEX_LAST_STEP][j].pm_radius);
				}
			}
		}
		while(!goon);

		IONS[INDEX_LAST_STEP][ion_index_1].x=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].x;
		IONS[INDEX_LAST_STEP][ion_index_1].y=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].y;
		IONS[INDEX_LAST_STEP][ion_index_1].z=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].z;
		
	}
	
	return;
}

void duplicate_last_configuration(){

	int NEXT_STEP_INDEX=(INDEX_LAST_STEP+1)%50;
	LEFT_CELLS[NEXT_STEP_INDEX]=LEFT_CELLS[INDEX_LAST_STEP];
	RIGHT_CELLS[NEXT_STEP_INDEX]=RIGHT_CELLS[INDEX_LAST_STEP];
	
	for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
		IONS[NEXT_STEP_INDEX][i]=IONS[INDEX_LAST_STEP][i];
	}	

	NUM_OF_IONS_IN_STEP[NEXT_STEP_INDEX]=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP];
	INDEX_LAST_STEP=NEXT_STEP_INDEX;
	STEPS[INDEX_LAST_STEP]=step;

	return;
}

void duplicate_first_configuration(){
	
	
	
	for(int i=ions.size()-1; i>0; i--){	
		ions.erase(ions.begin()+i);
	}
		
	
	
	vector <Ion> aux_ions;
	
	for(int i=1; i<ions.back().size(); i++){
		Ion backup_ion(ions.back().at(i).kind);
		copy_ion(ions.back().at(i), backup_ion);
		aux_ions.push_back(backup_ion);
	}
	
	ions.push_back(aux_ions);
	
	return;
}

void do_initial_step(){
	
	bool CORRECT=true;
	do{
		NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=0;
		ions.clear();
		CORRECT=initialize_ion_boxes();
		
		//~ cout << " 1 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		//~ for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){
			//~ cout << ion_index<< " " << IONS[INDEX_LAST_STEP][ion_index].name<< " "<< IONS[INDEX_LAST_STEP][ion_index].kind<< " "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].x <<" "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].y<<" "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].z<<endl;
		//~ }
		//~ cout << " 1 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	
	string hhh="\nOK FATTO "; 	
	string hhh2="\nOK FATTO "; 	
		
		if(CORRECT){
			initialize_ions_left_cell();
			initialize_ions_right_cell();
			initialize_ions_left_bath();
			initialize_ions_right_bath();
			//~ cout << " 2 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
			//~ for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){
				//~ cout << ion_index<< " " << IONS[INDEX_LAST_STEP][ion_index].name<< " "<< IONS[INDEX_LAST_STEP][ion_index].kind<< " "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].x <<" "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].y<<" "<< 1e10*IONS[INDEX_LAST_STEP][ion_index].z<<endl;
			//~ }
			//~ cout << " 2 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
			//~ getline(cin, hhh);
			
			for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
				if(IONS[INDEX_LAST_STEP][ion_index].name.substr(0, 1).compare("J")){
					IONS[INDEX_LAST_STEP][ion_index].set_brownian_parameters();
				}
			}
			cout << "andata..." <<endl;
			cout << "number of ions: " << NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP] <<endl;
			langevin_step_0();
			cout << "1..." <<endl;
			CORRECT=CORRECT*motion_check_step_0();
			cout << "2..." <<endl;
		}
	}
	while(!CORRECT);
	cout << "3..." <<endl;

	motion_check_correctness_step_0();
	cout << "4..." <<endl;
	motion_prepare_next_step();
	cout << "5..." <<endl;
	cout << "OK!\n"<<flush;
	string ions_pdb_file=PRM.PREFIX+".pdb";
	createFile(ions_pdb_file);	

	return;
}
