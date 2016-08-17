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




#include "constants.h"
#include "utils.h"
#include "classes.h"
#include "file_functions.h"
#include "input_output.h"
#include "sim_domain.h"
#include "ions_functions.h"
#include "PACO.h"
#include "sim_structures.h"




void PACO_update_left_CC(){
	
	//eliminate ions from the other control cell
	for(int i=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]-1; i>=0; i--){
		if(IONS[INDEX_LAST_STEP][i].side!=LEFT_CELLS[INDEX_LAST_STEP].SIDE){
			if((IONS[INDEX_LAST_STEP][i].z<1e-12*PRM.LEFT_CELL_MIN_Z && IONS[INDEX_LAST_STEP][i].z_next>=1e-12*PRM.LEFT_CELL_MIN_Z) || (IONS[INDEX_LAST_STEP][i].z>=1e-12*PRM.LEFT_CELL_MAX_Z && IONS[INDEX_LAST_STEP][i].z_next<1e-12*PRM.LEFT_CELL_MAX_Z)){
				if(IONS[INDEX_LAST_STEP][i].name.substr(0, 1).compare("J")){
					NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=remove_ion(INDEX_LAST_STEP, i);
				}
			}
		}
	}	
	
	for(int is=0; is<31; is++){
		if(PRM.ions_to_simulate.at(is)){

			//update SumEta
			if(LEFT_CELLS[INDEX_LAST_STEP].SumEta.at(is)>1e9){
				LEFT_CELLS[INDEX_LAST_STEP].SumEta.at(is)=LEFT_CELLS[INDEX_LAST_STEP].SumEta.at(is)-1e9;
				LEFT_CELLS[INDEX_LAST_STEP].SumNi.at(is)=LEFT_CELLS[INDEX_LAST_STEP].SumNi.at(is)-1e9;
			}
			LEFT_CELLS[INDEX_LAST_STEP].SumEta.at(is)+=LEFT_CELLS[INDEX_LAST_STEP].eta.at(is);
			
			
			//update SumNi
			for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
				if(IONS[INDEX_LAST_STEP][i].kind==is && IONS[INDEX_LAST_STEP][i].z_next>=1e-12*PRM.LEFT_CELL_MIN_Z && IONS[INDEX_LAST_STEP][i].z_next<1e-12*PRM.LEFT_CELL_MAX_Z){
					LEFT_CELLS[INDEX_LAST_STEP].SumNi.at(IONS[INDEX_LAST_STEP][i].kind)=LEFT_CELLS[INDEX_LAST_STEP].SumNi.at(IONS[INDEX_LAST_STEP][i].kind)+1;
				}
			}
			
			//evaluate Xi
			LEFT_CELLS[INDEX_LAST_STEP].Xi.at(is)=LEFT_CELLS[INDEX_LAST_STEP].SumEta.at(is)-LEFT_CELLS[INDEX_LAST_STEP].SumNi.at(is);
			
			//inserting or removing ions
			if(LEFT_CELLS[INDEX_LAST_STEP].Xi.at(is)<0){
				int nitr=-ceil(LEFT_CELLS[INDEX_LAST_STEP].Xi.at(is));
				
				for(int ind_nitr=0; ind_nitr<nitr; ind_nitr++){
				
					bool found=false;
					int ion_index=-1;
					for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
						if(IONS[INDEX_LAST_STEP][i].side!=LEFT_CELLS[INDEX_LAST_STEP].SIDE && IONS[INDEX_LAST_STEP][i].z_next>=1e-12*PRM.LEFT_CELL_MIN_Z && IONS[INDEX_LAST_STEP][i].z_next<1e-12*PRM.LEFT_CELL_MAX_Z && IONS[INDEX_LAST_STEP][i].kind==is && !found){
							found=true;
							ion_index=i;
						}
					}
					
					if(ion_index!=-1){
						NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=remove_ion(INDEX_LAST_STEP, ion_index);
					}
				}
			}
			else{
				int nitr=floor(LEFT_CELLS[INDEX_LAST_STEP].Xi.at(is));
				
				for(int ind_nitr=0; ind_nitr<nitr; ind_nitr++){
					PACO_step_0_left_CC(is);
				}
				
			}
			
		}
	}

	
	return;
	
}

void PACO_step_0_left_CC(int is){
	


	NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=insert_ion(is, INDEX_LAST_STEP);
	
	int ion_index_1=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]-1;
	IONS[INDEX_LAST_STEP][ion_index_1].side=LEFT_CELLS[INDEX_LAST_STEP].SIDE;
	

	int attempts_counter=0;
	bool goon=true;
	do{
	
		goon=true;
		
		IONS[INDEX_LAST_STEP][ion_index_1].x=PRM.MIN_X+PRM.SIM_DOMAIN_WIDTH_X*((float)rand()/((float)(RAND_MAX)+(float)(1)));
		IONS[INDEX_LAST_STEP][ion_index_1].y=PRM.MIN_Y+PRM.SIM_DOMAIN_WIDTH_Y*((float)rand()/((float)(RAND_MAX)+(float)(1)));
		IONS[INDEX_LAST_STEP][ion_index_1].z=PRM.LEFT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH*((double)rand()/((double)(RAND_MAX)+(double)(1)));

	
		for(int j=0; j<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; j++){
			goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, IONS[INDEX_LAST_STEP][ion_index_1].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, IONS[INDEX_LAST_STEP][j].pm_radius);
		}
	
		
		for(int ion_index_2=0; ion_index_2<ion_index_1; ion_index_2++){
			goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, double(1.25)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius, 1e12*IONS[INDEX_LAST_STEP][ion_index_2].x, 1e12*IONS[INDEX_LAST_STEP][ion_index_2].y, 1e12*IONS[INDEX_LAST_STEP][ion_index_2].z, double(1.25)*IONS[INDEX_LAST_STEP][ion_index_2].pm_radius);
		}
		
	
		attempts_counter++;
		if(attempts_counter==100){
			cout << "RIGHT CC attempts_counter:\t=100" << endl;
			goon=true;
		}
	}
	while(!goon);
	
	if(attempts_counter<100){
	
		IONS[INDEX_LAST_STEP][ion_index_1].x=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].x;
		IONS[INDEX_LAST_STEP][ion_index_1].y=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].y;
		IONS[INDEX_LAST_STEP][ion_index_1].z=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].z;
		initialize_ion_velocity(IONS[INDEX_LAST_STEP][ion_index_1]);

		IONS[INDEX_LAST_STEP][ion_index_1].force[0]=0.00;
		IONS[INDEX_LAST_STEP][ion_index_1].force[1]=0.00;
		IONS[INDEX_LAST_STEP][ion_index_1].force[2]=0.00;

		
		// mobile ions
		for(int ion_index_2=0; ion_index_2<ion_index_1; ion_index_2++){
			double ion2_x=IONS[INDEX_LAST_STEP][ion_index_2].x;
			double ion2_y=IONS[INDEX_LAST_STEP][ion_index_2].y;
			double ion2_z=IONS[INDEX_LAST_STEP][ion_index_2].z;			
			apply_periodic_boundary(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, ion2_x, ion2_y, ion2_z);
			double distance=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, ion2_x, ion2_y, ion2_z);			
			
			if(distance<20000){
				double force_C=0.00;
				double force_SR=0.00;
				double dx=0.00;
				double dy=0.00;
				double dz=0.00;
				
				dx=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].x-ion2_x)/distance;
				dy=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].y-ion2_y)/distance;
				dz=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].z-ion2_z)/distance;
				
				force_C=IONS[INDEX_LAST_STEP][ion_index_1].valence*IONS[INDEX_LAST_STEP][ion_index_2].DW_valence*FORCE_C[int(distance)];
				
				if(distance<4000){
					force_SR=FORCE_SR[IONS[INDEX_LAST_STEP][ion_index_1].kind][IONS[INDEX_LAST_STEP][ion_index_2].kind].at(int(distance));
				}
				
				IONS[INDEX_LAST_STEP][ion_index_1].force[0]+=(force_C+force_SR)*dx;
				IONS[INDEX_LAST_STEP][ion_index_1].force[1]+=(force_C+force_SR)*dy;
				IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=(force_C+force_SR)*dz;
			}	
			else{
				// cut-off
			}			
		}
		
		

		//step 0	
		IONS[INDEX_LAST_STEP][ion_index_1].X_n[0]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index_1].K_X;
		IONS[INDEX_LAST_STEP][ion_index_1].X_n[1]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index_1].K_X;
		IONS[INDEX_LAST_STEP][ion_index_1].X_n[2]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index_1].K_X;

		IONS[INDEX_LAST_STEP][ion_index_1].x_next=IONS[INDEX_LAST_STEP][ion_index_1].x+IONS[INDEX_LAST_STEP][ion_index_1].velocity[0]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_1+IONS[INDEX_LAST_STEP][ion_index_1].force[0]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_2+IONS[INDEX_LAST_STEP][ion_index_1].X_n[0];
		IONS[INDEX_LAST_STEP][ion_index_1].y_next=IONS[INDEX_LAST_STEP][ion_index_1].y+IONS[INDEX_LAST_STEP][ion_index_1].velocity[1]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_1+IONS[INDEX_LAST_STEP][ion_index_1].force[1]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_2+IONS[INDEX_LAST_STEP][ion_index_1].X_n[1];
		IONS[INDEX_LAST_STEP][ion_index_1].z_next=IONS[INDEX_LAST_STEP][ion_index_1].z+IONS[INDEX_LAST_STEP][ion_index_1].velocity[2]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_1+IONS[INDEX_LAST_STEP][ion_index_1].force[2]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_2+IONS[INDEX_LAST_STEP][ion_index_1].X_n[2];
		
			
		//prepare next step
		IONS[INDEX_LAST_STEP][ion_index_1].x_prev=IONS[INDEX_LAST_STEP][ion_index_1].x;
		IONS[INDEX_LAST_STEP][ion_index_1].x=IONS[INDEX_LAST_STEP][ion_index_1].x_next;
		IONS[INDEX_LAST_STEP][ion_index_1].force_old[0]=IONS[INDEX_LAST_STEP][ion_index_1].force[0];
		IONS[INDEX_LAST_STEP][ion_index_1].X_n_old[0]=IONS[INDEX_LAST_STEP][ion_index_1].X_n[0];
		
		IONS[INDEX_LAST_STEP][ion_index_1].y_prev=IONS[INDEX_LAST_STEP][ion_index_1].y;
		IONS[INDEX_LAST_STEP][ion_index_1].y=IONS[INDEX_LAST_STEP][ion_index_1].y_next;
		IONS[INDEX_LAST_STEP][ion_index_1].force_old[1]=IONS[INDEX_LAST_STEP][ion_index_1].force[1];
		IONS[INDEX_LAST_STEP][ion_index_1].X_n_old[1]=IONS[INDEX_LAST_STEP][ion_index_1].X_n[1];
		
		IONS[INDEX_LAST_STEP][ion_index_1].z_prev=IONS[INDEX_LAST_STEP][ion_index_1].z;
		IONS[INDEX_LAST_STEP][ion_index_1].z=IONS[INDEX_LAST_STEP][ion_index_1].z_next;
		IONS[INDEX_LAST_STEP][ion_index_1].force_old[2]=IONS[INDEX_LAST_STEP][ion_index_1].force[2];
		IONS[INDEX_LAST_STEP][ion_index_1].X_n_old[2]=IONS[INDEX_LAST_STEP][ion_index_1].X_n[2];
		
	}
	else{
		cout << "RIGHT CC attempts_counter:\tREMOVING ION" << endl << endl;
		NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=remove_ion(INDEX_LAST_STEP, ion_index_1);
	}
	
	return;
}

void PACO_update_right_CC(){

	//eliminate ions from the other control cell
	for(int i=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]-1; i>=0; i--){
		if(IONS[INDEX_LAST_STEP][i].side!=RIGHT_CELLS[INDEX_LAST_STEP].SIDE){
			if((IONS[INDEX_LAST_STEP][i].z<1e-12*PRM.RIGHT_CELL_MIN_Z && IONS[INDEX_LAST_STEP][i].z_next>=1e-12*PRM.RIGHT_CELL_MIN_Z) || (IONS[INDEX_LAST_STEP][i].z>=1e-12*PRM.RIGHT_CELL_MAX_Z && IONS[INDEX_LAST_STEP][i].z_next<1e-12*PRM.RIGHT_CELL_MAX_Z)){
				if(IONS[INDEX_LAST_STEP][i].name.substr(0, 1).compare("J")){
					NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=remove_ion(INDEX_LAST_STEP, i);
				}
			}
		}
	}
	
	for(int is=0; is<31; is++){
		if(PRM.ions_to_simulate.at(is)){

			//update SumEta
			if(RIGHT_CELLS[INDEX_LAST_STEP].SumEta.at(is)>1e9){
				RIGHT_CELLS[INDEX_LAST_STEP].SumEta.at(is)=RIGHT_CELLS[INDEX_LAST_STEP].SumEta.at(is)-1e9;
				RIGHT_CELLS[INDEX_LAST_STEP].SumNi.at(is)=RIGHT_CELLS[INDEX_LAST_STEP].SumNi.at(is)-1e9;
			}
			RIGHT_CELLS[INDEX_LAST_STEP].SumEta.at(is)+=RIGHT_CELLS[INDEX_LAST_STEP].eta.at(is);
			
			
			//update SumNi
			for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
				if(IONS[INDEX_LAST_STEP][i].kind==is && IONS[INDEX_LAST_STEP][i].z_next>=1e-12*PRM.RIGHT_CELL_MIN_Z && IONS[INDEX_LAST_STEP][i].z_next<1e-12*PRM.RIGHT_CELL_MAX_Z){
					RIGHT_CELLS[INDEX_LAST_STEP].SumNi.at(IONS[INDEX_LAST_STEP][i].kind)=RIGHT_CELLS[INDEX_LAST_STEP].SumNi.at(IONS[INDEX_LAST_STEP][i].kind)+1;
				}
			}
			
			//evaluate Xi
			RIGHT_CELLS[INDEX_LAST_STEP].Xi.at(is)=RIGHT_CELLS[INDEX_LAST_STEP].SumEta.at(is)-RIGHT_CELLS[INDEX_LAST_STEP].SumNi.at(is);
			
			//inserting or removing ions
			if(RIGHT_CELLS[INDEX_LAST_STEP].Xi.at(is)<0){
				int nitr=-ceil(RIGHT_CELLS[INDEX_LAST_STEP].Xi.at(is));
				
				for(int ind_nitr=0; ind_nitr<nitr; ind_nitr++){
				
					
					//removing one ion
					bool found=false;
					int ion_index=-1;
					for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
						if(IONS[INDEX_LAST_STEP][i].side!=RIGHT_CELLS[INDEX_LAST_STEP].SIDE && IONS[INDEX_LAST_STEP][i].z_next>=1e-12*PRM.RIGHT_CELL_MIN_Z && IONS[INDEX_LAST_STEP][i].z_next<1e-12*PRM.RIGHT_CELL_MAX_Z && IONS[INDEX_LAST_STEP][i].kind==is && !found){
							found=true;
							ion_index=i;
						}
					}
					
					if(ion_index!=-1){
						NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=remove_ion(INDEX_LAST_STEP, ion_index);
					}
				}
			}
			else{
				int nitr=floor(RIGHT_CELLS[INDEX_LAST_STEP].Xi.at(is));
				
				for(int ind_nitr=0; ind_nitr<nitr; ind_nitr++){
					PACO_step_0_right_CC(is);
				}
			}
			
		}
	}
	
	
	return;
	
}

void PACO_step_0_right_CC(int is){
	


	NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=insert_ion(is, INDEX_LAST_STEP);
	
	int ion_index_1=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]-1;
	IONS[INDEX_LAST_STEP][ion_index_1].side=RIGHT_CELLS[INDEX_LAST_STEP].SIDE;
	
	int attempts_counter=0;
	bool goon=true;
	do{
	
		goon=true;
		
		IONS[INDEX_LAST_STEP][ion_index_1].x=PRM.MIN_X+PRM.SIM_DOMAIN_WIDTH_X*((float)rand()/((float)(RAND_MAX)+(float)(1)));
		IONS[INDEX_LAST_STEP][ion_index_1].y=PRM.MIN_Y+PRM.SIM_DOMAIN_WIDTH_Y*((float)rand()/((float)(RAND_MAX)+(float)(1)));
		IONS[INDEX_LAST_STEP][ion_index_1].z=PRM.RIGHT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH*((double)rand()/((double)(RAND_MAX)+(double)(1)));

	
		for(int j=0; j<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; j++){
			goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, IONS[INDEX_LAST_STEP][ion_index_1].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, IONS[INDEX_LAST_STEP][j].pm_radius);
		}
	
		
		for(int ion_index_2=0; ion_index_2<ion_index_1; ion_index_2++){
			goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, double(1.25)*IONS[INDEX_LAST_STEP][ion_index_1].pm_radius, 1e12*IONS[INDEX_LAST_STEP][ion_index_2].x, 1e12*IONS[INDEX_LAST_STEP][ion_index_2].y, 1e12*IONS[INDEX_LAST_STEP][ion_index_2].z, double(1.25)*IONS[INDEX_LAST_STEP][ion_index_2].pm_radius);
		}
		
	
		attempts_counter++;
		if(attempts_counter==100){
			cout << "RIGHT CC attempts_counter:\t=100" << endl;
			goon=true;
		}
	}
	while(!goon);
	
	if(attempts_counter<100){
	
		IONS[INDEX_LAST_STEP][ion_index_1].x=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].x;
		IONS[INDEX_LAST_STEP][ion_index_1].y=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].y;
		IONS[INDEX_LAST_STEP][ion_index_1].z=1e-12*IONS[INDEX_LAST_STEP][ion_index_1].z;
		initialize_ion_velocity(IONS[INDEX_LAST_STEP][ion_index_1]);

		IONS[INDEX_LAST_STEP][ion_index_1].force[0]=0.00;
		IONS[INDEX_LAST_STEP][ion_index_1].force[1]=0.00;
		IONS[INDEX_LAST_STEP][ion_index_1].force[2]=0.00;

		
		// mobile ions
		for(int ion_index_2=0; ion_index_2<ion_index_1; ion_index_2++){
			double ion2_x=IONS[INDEX_LAST_STEP][ion_index_2].x;
			double ion2_y=IONS[INDEX_LAST_STEP][ion_index_2].y;
			double ion2_z=IONS[INDEX_LAST_STEP][ion_index_2].z;			
			apply_periodic_boundary(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, ion2_x, ion2_y, ion2_z);
			double distance=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, ion2_x, ion2_y, ion2_z);			
			
			if(distance<20000){
				double force_C=0.00;
				double force_SR=0.00;
				double dx=0.00;
				double dy=0.00;
				double dz=0.00;
				
				dx=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].x-ion2_x)/distance;
				dy=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].y-ion2_y)/distance;
				dz=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].z-ion2_z)/distance;
				
				force_C=IONS[INDEX_LAST_STEP][ion_index_1].valence*IONS[INDEX_LAST_STEP][ion_index_2].DW_valence*FORCE_C[int(distance)];
				
				if(distance<4000){
					force_SR=FORCE_SR[IONS[INDEX_LAST_STEP][ion_index_1].kind][IONS[INDEX_LAST_STEP][ion_index_2].kind].at(int(distance));
				}
				
				IONS[INDEX_LAST_STEP][ion_index_1].force[0]+=(force_C+force_SR)*dx;
				IONS[INDEX_LAST_STEP][ion_index_1].force[1]+=(force_C+force_SR)*dy;
				IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=(force_C+force_SR)*dz;
			}	
			else{
				// cut-off
			}			
		}
		
		

		//step 0	
		IONS[INDEX_LAST_STEP][ion_index_1].X_n[0]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index_1].K_X;
		IONS[INDEX_LAST_STEP][ion_index_1].X_n[1]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index_1].K_X;
		IONS[INDEX_LAST_STEP][ion_index_1].X_n[2]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index_1].K_X;

		IONS[INDEX_LAST_STEP][ion_index_1].x_next=IONS[INDEX_LAST_STEP][ion_index_1].x+IONS[INDEX_LAST_STEP][ion_index_1].velocity[0]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_1+IONS[INDEX_LAST_STEP][ion_index_1].force[0]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_2+IONS[INDEX_LAST_STEP][ion_index_1].X_n[0];
		IONS[INDEX_LAST_STEP][ion_index_1].y_next=IONS[INDEX_LAST_STEP][ion_index_1].y+IONS[INDEX_LAST_STEP][ion_index_1].velocity[1]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_1+IONS[INDEX_LAST_STEP][ion_index_1].force[1]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_2+IONS[INDEX_LAST_STEP][ion_index_1].X_n[1];
		IONS[INDEX_LAST_STEP][ion_index_1].z_next=IONS[INDEX_LAST_STEP][ion_index_1].z+IONS[INDEX_LAST_STEP][ion_index_1].velocity[2]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_1+IONS[INDEX_LAST_STEP][ion_index_1].force[2]*IONS[INDEX_LAST_STEP][ion_index_1].K_step0_2+IONS[INDEX_LAST_STEP][ion_index_1].X_n[2];
		
			
		//prepare next step
		IONS[INDEX_LAST_STEP][ion_index_1].x_prev=IONS[INDEX_LAST_STEP][ion_index_1].x;
		IONS[INDEX_LAST_STEP][ion_index_1].x=IONS[INDEX_LAST_STEP][ion_index_1].x_next;
		IONS[INDEX_LAST_STEP][ion_index_1].force_old[0]=IONS[INDEX_LAST_STEP][ion_index_1].force[0];
		IONS[INDEX_LAST_STEP][ion_index_1].X_n_old[0]=IONS[INDEX_LAST_STEP][ion_index_1].X_n[0];
		
		IONS[INDEX_LAST_STEP][ion_index_1].y_prev=IONS[INDEX_LAST_STEP][ion_index_1].y;
		IONS[INDEX_LAST_STEP][ion_index_1].y=IONS[INDEX_LAST_STEP][ion_index_1].y_next;
		IONS[INDEX_LAST_STEP][ion_index_1].force_old[1]=IONS[INDEX_LAST_STEP][ion_index_1].force[1];
		IONS[INDEX_LAST_STEP][ion_index_1].X_n_old[1]=IONS[INDEX_LAST_STEP][ion_index_1].X_n[1];
		
		IONS[INDEX_LAST_STEP][ion_index_1].z_prev=IONS[INDEX_LAST_STEP][ion_index_1].z;
		IONS[INDEX_LAST_STEP][ion_index_1].z=IONS[INDEX_LAST_STEP][ion_index_1].z_next;
		IONS[INDEX_LAST_STEP][ion_index_1].force_old[2]=IONS[INDEX_LAST_STEP][ion_index_1].force[2];
		IONS[INDEX_LAST_STEP][ion_index_1].X_n_old[2]=IONS[INDEX_LAST_STEP][ion_index_1].X_n[2];
		
	}
	else{
		cout << "RIGHT CC attempts_counter:\tREMOVING ION" << endl << endl;
		NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=remove_ion(INDEX_LAST_STEP, ion_index_1);
	}
	
	return;
}
