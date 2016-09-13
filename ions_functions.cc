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

#include <cstring>
#include "constants.h"
#include "utils.h"
#include "classes.h"
#include "file_functions.h"
#include "input_output.h"
#include "sim_domain.h"
#include "ions_functions.h"
#include "retrace.h"
#include "sim_structures.h"



int estimate_max_number_of_ions(){
	
	int MAX_NUMBER_OF_IONS=0;
	
// 1 - ion boxes
	vector <Ion> aux_ions;
	
	for(int ib=0; ib<ion_boxes.size(); ib++){
		for(int iv=0; iv<ion_boxes.at(ib).is.size(); iv++){
			for(int ino=0; ino<ion_boxes.at(ib).in.at(iv); ino++){
				for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
					Ion ion(is);
					if(ion_boxes.at(ib).is.at(iv).compare(ion.name)==0){
						MAX_NUMBER_OF_IONS++;
					}
				}
			}
		}
	}
	
//2 - left side
	double bulk_volume=1e-33*PRM.SIM_DOMAIN_WIDTH_X*PRM.SIM_DOMAIN_WIDTH_Y*(PRM.CONTROL_CELL_WIDTH+PRM.BATH_WIDTH);
	vector <int> Nis;
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		Nis.push_back(0);
	}
	
	Nis.at(4)=round(PRM.CONC_LEFT_KCL*AVOGADRO*bulk_volume); 
	Nis.at(17)=round((PRM.CONC_LEFT_KCL+double(2.00)*PRM.CONC_LEFT_CACL2+double(2.00)*PRM.CONC_LEFT_MGCL2+PRM.CONC_LEFT_NACL)*AVOGADRO*bulk_volume); 
	Nis.at(9)=round(PRM.CONC_LEFT_CACL2*AVOGADRO*bulk_volume);
	Nis.at(8)=round(PRM.CONC_LEFT_MGCL2*AVOGADRO*bulk_volume);	
	Nis.at(3)=round(PRM.CONC_LEFT_NACL*AVOGADRO*bulk_volume); 
	
	
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		for(int ii=0; ii<Nis.at(is); ii++){			
			MAX_NUMBER_OF_IONS++;
		}
	}	
	
//2 - right side
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		Nis.at(is)=0;
	}
	
	Nis.at(4)=round(PRM.CONC_RIGHT_KCL*AVOGADRO*bulk_volume); 
	Nis.at(17)=round((PRM.CONC_RIGHT_KCL+double(2.00)*PRM.CONC_RIGHT_CACL2+double(2.00)*PRM.CONC_RIGHT_MGCL2+PRM.CONC_RIGHT_NACL)*AVOGADRO*bulk_volume); 
	Nis.at(9)=round(PRM.CONC_RIGHT_CACL2*AVOGADRO*bulk_volume);
	Nis.at(8)=round(PRM.CONC_RIGHT_MGCL2*AVOGADRO*bulk_volume);	
	Nis.at(3)=round(PRM.CONC_RIGHT_NACL*AVOGADRO*bulk_volume); 
	
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		for(int ii=0; ii<Nis.at(is); ii++){			
			MAX_NUMBER_OF_IONS++;
		}
	}		
	
	MAX_NUMBER_OF_IONS=ceil(double(1.40)*double(MAX_NUMBER_OF_IONS));
	
	
	return MAX_NUMBER_OF_IONS;
	
}

int insert_ion(int type, int stepIndex){
	
	IONS[stepIndex][NUM_OF_IONS_IN_STEP[stepIndex]]=Ion(type);	
	
	return (NUM_OF_IONS_IN_STEP[stepIndex]+1);
}



int remove_ion(int stepIndex, int ionIndex){
	
	for(int i=ionIndex; i<NUM_OF_IONS_IN_STEP[stepIndex]-1; i++){
		//~ cout << i << "\t" << i+1 << "\t" << ionIndex << "\t"
		IONS[stepIndex][i]=IONS[stepIndex][i+1];
	}

	return NUM_OF_IONS_IN_STEP[stepIndex]-1;
}



void initialize_left_control_cell(){
	
	Control_cell cell;
	left_cell.push_back(cell);
	
	left_cell.back().reset_control_cell();

	left_cell.back().SIDE=0;
	
//===================================================	
// to be updated with different ionic solutions		
																	/*

	0		O		OXYGEN
	1		H		HYDROGEN
	2		LI		LITHIUM
	3		NA		SODIUM
	4		K		POTASSIUM
	5		RB		RUBIDIUM
	6		CS		CAESIUM
	7		TL		THALLIUM
	8		MG		MAGNESIUM
	9		CA		CALCIUM
	10		SR		STRONTIUM
	11		BA		BARIUM
	12		MN		MANGANESE
	13		CO		COBALT
	14		NI		NICKEL
	15		ZN		ZINC
	16		F		FLUORINE
	17		CL		CHLORINE
	18		BR		BROMINE
	19		I		IODINE
																	*/
//===================================================

	left_cell.back().eta.at(3)=(PRM.CONC_LEFT_NACL*AVOGADRO*PRM.CONTROL_CELL_VOLUME);
	left_cell.back().eta.at(4)=(PRM.CONC_LEFT_KCL*AVOGADRO*PRM.CONTROL_CELL_VOLUME);
	left_cell.back().eta.at(9)=(PRM.CONC_LEFT_CACL2*AVOGADRO*PRM.CONTROL_CELL_VOLUME);
	left_cell.back().eta.at(8)=(PRM.CONC_LEFT_MGCL2*AVOGADRO*PRM.CONTROL_CELL_VOLUME);
	
	left_cell.back().eta.at(17)=(PRM.CONC_LEFT_NACL+PRM.CONC_LEFT_KCL+double(2.00)*PRM.CONC_LEFT_CACL2+double(2.00)*PRM.CONC_LEFT_MGCL2)*AVOGADRO*PRM.CONTROL_CELL_VOLUME;
	
	
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		if(left_cell.back().eta.at(is)>1e-10){
			left_cell.back().in_control_cell.at(is)=true;
		}	
	}

	cout << endl << "===========================================" <<endl;
	cout << left_cell.back().SIDE << " CELL"<< endl;
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		if(left_cell.back().in_control_cell.at(is)){
			Ion ion_curr(is);
			cout << "\n" << ion_curr.name << ": " << endl;
			cout << "eta: " << left_cell.back().eta.at(is) << endl;
			
			cout << "REF_CONC: " << left_cell.back().eta.at(is)*(100.00)/PRM.CONTROL_CELL_WIDTH << endl<<endl;
		}
	}
	cout << "===========================================" << endl << endl;
	
	return;
	
}

void initialize_right_control_cell(){
	
	Control_cell cell;
	right_cell.push_back(cell);
	
	right_cell.back().reset_control_cell();

	right_cell.back().SIDE=1;
	
//===================================================	
// to be updated with different ionic solutions		
																	/*

	0		O		OXYGEN
	1		H		HYDROGEN
	2		LI		LITHIUM
	3		NA		SODIUM
	4		K		POTASSIUM
	5		RB		RUBIDIUM
	6		CS		CAESIUM
	7		TL		THALLIUM
	8		MG		MAGNESIUM
	9		CA		CALCIUM
	10		SR		STRONTIUM
	11		BA		BARIUM
	12		MN		MANGANESE
	13		CO		COBALT
	14		NI		NICKEL
	15		ZN		ZINC
	16		F		FLUORINE
	17		CL		CHLORINE
	18		BR		BROMINE
	19		I		IODINE
																	*/
//===================================================

	right_cell.back().eta.at(3)=(PRM.CONC_RIGHT_NACL*AVOGADRO*PRM.CONTROL_CELL_VOLUME);
	right_cell.back().eta.at(4)=(PRM.CONC_RIGHT_KCL*AVOGADRO*PRM.CONTROL_CELL_VOLUME);
	right_cell.back().eta.at(9)=(PRM.CONC_RIGHT_CACL2*AVOGADRO*PRM.CONTROL_CELL_VOLUME);
	right_cell.back().eta.at(8)=(PRM.CONC_RIGHT_MGCL2*AVOGADRO*PRM.CONTROL_CELL_VOLUME);
	
	right_cell.back().eta.at(17)=(PRM.CONC_RIGHT_NACL+PRM.CONC_RIGHT_KCL+double(2.00)*PRM.CONC_RIGHT_CACL2+double(2.00)*PRM.CONC_RIGHT_MGCL2)*AVOGADRO*PRM.CONTROL_CELL_VOLUME;
	
	
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		if(right_cell.back().eta.at(is)>1e-10){
			right_cell.back().in_control_cell.at(is)=true;
		}	
	}

	cout << endl << "===========================================" <<endl;
	cout << right_cell.back().SIDE << " CELL"<< endl;
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		if(right_cell.back().in_control_cell.at(is)){
			Ion ion_curr(is);
			cout << "\n" << ion_curr.name << ": " << endl;
			cout << "eta: " << right_cell.back().eta.at(is) << endl;
			
			cout << "REF_CONC: " << right_cell.back().eta.at(is)*(100.00)/PRM.CONTROL_CELL_WIDTH << endl<<endl;
		}
	}
	cout << "===========================================" << endl << endl;
	
	return;
}
	
void initialize_ions_left_cell(){
	cout << "initialize_ions_left_cell - begin" <<endl;	

	int start_index=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP];

	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		int N0=round(left_cell.back().eta.at(is));
		
		for(int ii=0; ii<N0; ii++){			
			NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=insert_ion(is, INDEX_LAST_STEP);
		}
	}

// add other ionic solutions	
	for(int i=start_index; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
		
		bool goon=true;
		do{
			goon=true;

			IONS[INDEX_LAST_STEP][i].x=PRM.MIN_X+PRM.SIM_DOMAIN_WIDTH_X*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][i].y=PRM.MIN_Y+PRM.SIM_DOMAIN_WIDTH_Y*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][i].z=PRM.LEFT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			
			for(int j=0; j<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; j++){
				if(i!=j){
					goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][i].x, IONS[INDEX_LAST_STEP][i].y, IONS[INDEX_LAST_STEP][i].z, double(2.00)*IONS[INDEX_LAST_STEP][i].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, double(2.00)*IONS[INDEX_LAST_STEP][j].pm_radius);
				}
			}
			if(isnormal(IONS[INDEX_LAST_STEP][i].x) && isnormal(IONS[INDEX_LAST_STEP][i].y) && isnormal(IONS[INDEX_LAST_STEP][i].z)){
				goon=true;
			}
			else{
				goon=false;
			}
		}
		while(!goon);

		IONS[INDEX_LAST_STEP][i].x=1e-12*IONS[INDEX_LAST_STEP][i].x;
		IONS[INDEX_LAST_STEP][i].y=1e-12*IONS[INDEX_LAST_STEP][i].y;
		IONS[INDEX_LAST_STEP][i].z=1e-12*IONS[INDEX_LAST_STEP][i].z;
		
		//set the initial velocity
		initialize_ion_velocity(IONS[INDEX_LAST_STEP][i]);
	}

	cout << "initialize_ions_left_cell - end" <<endl;
	return;
}
	
void initialize_ions_right_cell(){
	cout << "initialize_ions_right_cell - begin" <<endl;

	int start_index=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP];
	
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		int N0=round(right_cell.back().eta.at(is));
		
		for(int ii=0; ii<N0; ii++){			
			NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=insert_ion(is, INDEX_LAST_STEP);
		}
	}

// add other ionic solutions	
	for(int i=start_index; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
		
		bool goon=true;
		do{
			goon=true;

			IONS[INDEX_LAST_STEP][i].x=PRM.MIN_X+PRM.SIM_DOMAIN_WIDTH_X*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][i].y=PRM.MIN_Y+PRM.SIM_DOMAIN_WIDTH_Y*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][i].z=PRM.RIGHT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			
			for(int j=0; j<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; j++){
				if(i!=j){
					goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][i].x, IONS[INDEX_LAST_STEP][i].y, IONS[INDEX_LAST_STEP][i].z, double(2.00)*IONS[INDEX_LAST_STEP][i].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, double(2.00)*IONS[INDEX_LAST_STEP][j].pm_radius);
				}
			}
			if(isnormal(IONS[INDEX_LAST_STEP][i].x) && isnormal(IONS[INDEX_LAST_STEP][i].y) && isnormal(IONS[INDEX_LAST_STEP][i].z)){
				goon=true;
			}
			else{
				goon=false;
			}
		}
		while(!goon);

		IONS[INDEX_LAST_STEP][i].x=1e-12*IONS[INDEX_LAST_STEP][i].x;
		IONS[INDEX_LAST_STEP][i].y=1e-12*IONS[INDEX_LAST_STEP][i].y;
		IONS[INDEX_LAST_STEP][i].z=1e-12*IONS[INDEX_LAST_STEP][i].z;
		
		//set the initial velocity
		initialize_ion_velocity(IONS[INDEX_LAST_STEP][i]);
	}
	
	cout << "initialize_ions_right_cell - end" <<endl;
	return;
}

void initialize_ions_left_bath(){
	cout << "initialize_ions_left_bath - begin" <<endl;

	Control_cell buffer;
	buffer.reset_control_cell();
	
	buffer.SIDE=-1;
	
	int start_index=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP];
	
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		int N0=0;
		if(left_cell.back().eta.at(is)>0){
			N0=round(left_cell.back().eta.at(is)*PRM.BATH_WIDTH/PRM.CONTROL_CELL_WIDTH);
			if(N0<3){
				N0=3;
			}
			
			for(int ii=0; ii<N0; ii++){			
				NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=insert_ion(is, INDEX_LAST_STEP);
			}
		}
	}

// add other ionic solutions	
	for(int i=start_index; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
		
		bool goon=true;
		do{
			goon=true;

			IONS[INDEX_LAST_STEP][i].x=PRM.MIN_X+PRM.SIM_DOMAIN_WIDTH_X*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][i].y=PRM.MIN_Y+PRM.SIM_DOMAIN_WIDTH_Y*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][i].z=double(300)+PRM.LEFT_CELL_MAX_Z+(PRM.BATH_WIDTH-double(600))*((double)rand()/((double)(RAND_MAX)+(double)(1)));

			
			for(int j=0; j<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; j++){
				if(i!=j){
					goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][i].x, IONS[INDEX_LAST_STEP][i].y, IONS[INDEX_LAST_STEP][i].z, double(2.00)*IONS[INDEX_LAST_STEP][i].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, double(2.00)*IONS[INDEX_LAST_STEP][j].pm_radius);
				}
			}
			if(isnormal(IONS[INDEX_LAST_STEP][i].x) && isnormal(IONS[INDEX_LAST_STEP][i].y) && isnormal(IONS[INDEX_LAST_STEP][i].z)){
				goon=true;
			}
			else{
				goon=false;
			}
		}
		while(!goon);

		IONS[INDEX_LAST_STEP][i].x=1e-12*IONS[INDEX_LAST_STEP][i].x;
		IONS[INDEX_LAST_STEP][i].y=1e-12*IONS[INDEX_LAST_STEP][i].y;
		IONS[INDEX_LAST_STEP][i].z=1e-12*IONS[INDEX_LAST_STEP][i].z;
		
		//set the initial velocity
		initialize_ion_velocity(IONS[INDEX_LAST_STEP][i]);
	}

	cout << "initialize_ions_left_buffer - end" <<endl;
	return;
}

void initialize_ions_right_bath(){
	cout << "initialize_ions_right_bath - begin" <<endl;
//=============================================
// create ions in the left cell					
															/*

	0		O		OXYGEN
	1		H		HYDROGEN
	2		LI		LITHIUM
	3		NA		SODIUM
	4		K		POTASSIUM
	5		RB		RUBIDIUM
	6		CS		CAESIUM
	7		TL		THALLIUM
	8		MG		MAGNESIUM
	9		CA		CALCIUM
	10		SR		STRONTIUM
	11		BA		BARIUM
	12		MN		MANGANESE
	13		CO		COBALT
	14		NI		NICKEL
	15		ZN		ZINC
	16		F		FLUORINE
	17		CL		CHLORINE
	18		BR		BROMINE
	19		I		IODINE
														*/
//===========================================

	Control_cell buffer;
	buffer.reset_control_cell();
	
	buffer.SIDE=-1;
	
	int start_index=NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP];

	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		int N0=0;
		if(right_cell.back().eta.at(is)>0){
			N0=round(right_cell.back().eta.at(is)*PRM.BATH_WIDTH/PRM.CONTROL_CELL_WIDTH);
			if(N0<3){
				N0=3;
			}
			
			for(int ii=0; ii<N0; ii++){			
				NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]=insert_ion(is, INDEX_LAST_STEP);
			}
		}
	}

// add other ionic solutions	
	for(int i=start_index; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
		
		bool goon=true;
		do{
			goon=true;

			IONS[INDEX_LAST_STEP][i].x=PRM.MIN_X+PRM.SIM_DOMAIN_WIDTH_X*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][i].y=PRM.MIN_Y+PRM.SIM_DOMAIN_WIDTH_Y*((double)rand()/((double)(RAND_MAX)+(double)(1)));
			IONS[INDEX_LAST_STEP][i].z=PRM.RIGHT_CELL_MIN_Z-(PRM.BATH_WIDTH-double(300))+(PRM.BATH_WIDTH-double(600))*((double)rand()/((double)(RAND_MAX)+(double)(1)));			
			
			for(int j=0; j<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; j++){
				if(i!=j){
					goon=goon*verifyDistance(IONS[INDEX_LAST_STEP][i].x, IONS[INDEX_LAST_STEP][i].y, IONS[INDEX_LAST_STEP][i].z, double(2.00)*IONS[INDEX_LAST_STEP][i].pm_radius, 1e12*IONS[INDEX_LAST_STEP][j].x, 1e12*IONS[INDEX_LAST_STEP][j].y, 1e12*IONS[INDEX_LAST_STEP][j].z, double(2.00)*IONS[INDEX_LAST_STEP][j].pm_radius);
				}
			}
			if(isnormal(IONS[INDEX_LAST_STEP][i].x) && isnormal(IONS[INDEX_LAST_STEP][i].y) && isnormal(IONS[INDEX_LAST_STEP][i].z)){
				goon=true;
			}
			else{
				goon=false;
			}
		}
		while(!goon);

		IONS[INDEX_LAST_STEP][i].x=1e-12*IONS[INDEX_LAST_STEP][i].x;
		IONS[INDEX_LAST_STEP][i].y=1e-12*IONS[INDEX_LAST_STEP][i].y;
		IONS[INDEX_LAST_STEP][i].z=1e-12*IONS[INDEX_LAST_STEP][i].z;
		
		//set the initial velocity
		initialize_ion_velocity(IONS[INDEX_LAST_STEP][i]);
	}
	
	cout << "initialize_ions_right_buffer - end" <<endl;
	
	return;
}




void initialize_ion_velocity(Ion& ion){
	
	double temp=sqrt((BOLTZMANN_K*PRM.TEMPERATURE)/ion.mass);
	ion.velocity[0]=temp*gauss_cut();
	ion.velocity[1]=temp*gauss_cut();
	ion.velocity[2]=temp*gauss_cut();
	
	return;
}

void compute_diffusion_coefficients(){
	cout << "compute diffusion coefficient" <<endl;
	double R_MIN=10000;
	double R_MAX=-10000;
	for(int i=0; i<limits.size(); i++){
		if(limits.at(i)<=R_MIN){
			R_MIN=limits.at(i);
		}
		if(limits.at(i)>=R_MAX){
			R_MAX=limits.at(i);
		}
	}
	vector <double> aux_vec;
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		Ion ion_i(i);
		aux_vec.clear();
		if(PRM.ions_to_simulate.at(i) && ion_i.name.substr(0,1).compare("J")!=0){
			for(int z=PRM.MIN_Z; z<=PRM.MAX_Z; z++){
				double diff_coeff=0;
				if(z<PRM.Z_MOUTH_LEFT || z>PRM.Z_MOUTH_RIGHT){
					diff_coeff=ion_i.diffusion_coeff;
				}
				else{
					int ind_on_lim=z-PRM.Z_MOUTH_LEFT;
					diff_coeff=ion_i.diffusion_coeff*(PRM.DIFF_COEFF_IN_CHANNEL+(double(1.00)-PRM.DIFF_COEFF_IN_CHANNEL)*(limits.at(ind_on_lim)-R_MIN)/(R_MAX-R_MIN));
				}
				aux_vec.push_back(diff_coeff);
			}
		}
		diffusion_coefficients.push_back(aux_vec);
		aux_vec.clear();
	}
	return;
}

void langevin_compute_new_position_and_velocity(){
//========	langevin equation of motion is integrated with Van Gunsteren and Berendsen method

	compute_force_on_ions();

	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
		
		double Fx_prime=(IONS[INDEX_LAST_STEP][ion_index].force[0]-IONS[INDEX_LAST_STEP][ion_index].force_old[0])/(PRM.DELTA_T);
		double Fy_prime=(IONS[INDEX_LAST_STEP][ion_index].force[1]-IONS[INDEX_LAST_STEP][ion_index].force_old[1])/(PRM.DELTA_T);
		double Fz_prime=(IONS[INDEX_LAST_STEP][ion_index].force[2]-IONS[INDEX_LAST_STEP][ion_index].force_old[2])/(PRM.DELTA_T);
		
		double Yx=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_Y;
		double Yy=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_Y;
		double Yz=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_Y;
		
		IONS[INDEX_LAST_STEP][ion_index].X_n_minus[0]=IONS[INDEX_LAST_STEP][ion_index].X_n_old[0]*IONS[INDEX_LAST_STEP][ion_index].G/IONS[INDEX_LAST_STEP][ion_index].C+Yx;
		IONS[INDEX_LAST_STEP][ion_index].X_n_minus[1]=IONS[INDEX_LAST_STEP][ion_index].X_n_old[1]*IONS[INDEX_LAST_STEP][ion_index].G/IONS[INDEX_LAST_STEP][ion_index].C+Yy;
		IONS[INDEX_LAST_STEP][ion_index].X_n_minus[2]=IONS[INDEX_LAST_STEP][ion_index].X_n_old[2]*IONS[INDEX_LAST_STEP][ion_index].G/IONS[INDEX_LAST_STEP][ion_index].C+Yz;
		
		IONS[INDEX_LAST_STEP][ion_index].X_n[0]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_X;
		IONS[INDEX_LAST_STEP][ion_index].X_n[1]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_X;
		IONS[INDEX_LAST_STEP][ion_index].X_n[2]=gauss_cut()*IONS[INDEX_LAST_STEP][ion_index].K_X;
				
		IONS[INDEX_LAST_STEP][ion_index].x_next=IONS[INDEX_LAST_STEP][ion_index].x*IONS[INDEX_LAST_STEP][ion_index].K_pos_1+IONS[INDEX_LAST_STEP][ion_index].x_prev*IONS[INDEX_LAST_STEP][ion_index].K_pos_2+IONS[INDEX_LAST_STEP][ion_index].force[0]*IONS[INDEX_LAST_STEP][ion_index].K_pos_3+Fx_prime*IONS[INDEX_LAST_STEP][ion_index].K_pos_4+IONS[INDEX_LAST_STEP][ion_index].X_n[0]+IONS[INDEX_LAST_STEP][ion_index].X_n_minus[0]*IONS[INDEX_LAST_STEP][ion_index].K_pos_5;
		IONS[INDEX_LAST_STEP][ion_index].y_next=IONS[INDEX_LAST_STEP][ion_index].y*IONS[INDEX_LAST_STEP][ion_index].K_pos_1+IONS[INDEX_LAST_STEP][ion_index].y_prev*IONS[INDEX_LAST_STEP][ion_index].K_pos_2+IONS[INDEX_LAST_STEP][ion_index].force[1]*IONS[INDEX_LAST_STEP][ion_index].K_pos_3+Fy_prime*IONS[INDEX_LAST_STEP][ion_index].K_pos_4+IONS[INDEX_LAST_STEP][ion_index].X_n[1]+IONS[INDEX_LAST_STEP][ion_index].X_n_minus[1]*IONS[INDEX_LAST_STEP][ion_index].K_pos_5;
		IONS[INDEX_LAST_STEP][ion_index].z_next=IONS[INDEX_LAST_STEP][ion_index].z*IONS[INDEX_LAST_STEP][ion_index].K_pos_1+IONS[INDEX_LAST_STEP][ion_index].z_prev*IONS[INDEX_LAST_STEP][ion_index].K_pos_2+IONS[INDEX_LAST_STEP][ion_index].force[2]*IONS[INDEX_LAST_STEP][ion_index].K_pos_3+Fz_prime*IONS[INDEX_LAST_STEP][ion_index].K_pos_4+IONS[INDEX_LAST_STEP][ion_index].X_n[2]+IONS[INDEX_LAST_STEP][ion_index].X_n_minus[2]*IONS[INDEX_LAST_STEP][ion_index].K_pos_5;

		IONS[INDEX_LAST_STEP][ion_index].velocity[0]=IONS[INDEX_LAST_STEP][ion_index].K_vel_3*(IONS[INDEX_LAST_STEP][ion_index].x_next-IONS[INDEX_LAST_STEP][ion_index].x_prev+IONS[INDEX_LAST_STEP][ion_index].force[0]*IONS[INDEX_LAST_STEP][ion_index].K_vel_1+Fx_prime*IONS[INDEX_LAST_STEP][ion_index].K_vel_2+IONS[INDEX_LAST_STEP][ion_index].X_n_minus[0]-IONS[INDEX_LAST_STEP][ion_index].X_n[0]);
		IONS[INDEX_LAST_STEP][ion_index].velocity[1]=IONS[INDEX_LAST_STEP][ion_index].K_vel_3*(IONS[INDEX_LAST_STEP][ion_index].y_next-IONS[INDEX_LAST_STEP][ion_index].y_prev+IONS[INDEX_LAST_STEP][ion_index].force[1]*IONS[INDEX_LAST_STEP][ion_index].K_vel_1+Fy_prime*IONS[INDEX_LAST_STEP][ion_index].K_vel_2+IONS[INDEX_LAST_STEP][ion_index].X_n_minus[1]-IONS[INDEX_LAST_STEP][ion_index].X_n[1]);
		IONS[INDEX_LAST_STEP][ion_index].velocity[2]=IONS[INDEX_LAST_STEP][ion_index].K_vel_3*(IONS[INDEX_LAST_STEP][ion_index].z_next-IONS[INDEX_LAST_STEP][ion_index].z_prev+IONS[INDEX_LAST_STEP][ion_index].force[2]*IONS[INDEX_LAST_STEP][ion_index].K_vel_1+Fz_prime*IONS[INDEX_LAST_STEP][ion_index].K_vel_2+IONS[INDEX_LAST_STEP][ion_index].X_n_minus[2]-IONS[INDEX_LAST_STEP][ion_index].X_n[2]);
		
	}
	
	return;
}

void motion_prepare_next_step(){

	for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){	
		if(PRM.ION_RECYCLING.substr(0, 2).compare("00")==0){
//reflective boundary, NO ion recycling				
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].x_next<PRM.MIN_X || 1e12*IONS[INDEX_LAST_STEP][ion_index].x_next>PRM.MAX_X){
				IONS[INDEX_LAST_STEP][ion_index].x_next=IONS[INDEX_LAST_STEP][ion_index].x;
				IONS[INDEX_LAST_STEP][ion_index].x=IONS[INDEX_LAST_STEP][ion_index].x_prev;
				IONS[INDEX_LAST_STEP][ion_index].x_prev=IONS[INDEX_LAST_STEP][ion_index].x_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[0]=-IONS[INDEX_LAST_STEP][ion_index].force[0];	
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[0]=-IONS[INDEX_LAST_STEP][ion_index].X_n[0];				
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].x_prev=IONS[INDEX_LAST_STEP][ion_index].x;
				IONS[INDEX_LAST_STEP][ion_index].x=IONS[INDEX_LAST_STEP][ion_index].x_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[0]=IONS[INDEX_LAST_STEP][ion_index].force[0];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[0]=IONS[INDEX_LAST_STEP][ion_index].X_n[0];
			}
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].y_next<PRM.MIN_Y || 1e12*IONS[INDEX_LAST_STEP][ion_index].y_next>PRM.MAX_Y){
				IONS[INDEX_LAST_STEP][ion_index].y_next=IONS[INDEX_LAST_STEP][ion_index].y;
				IONS[INDEX_LAST_STEP][ion_index].y=IONS[INDEX_LAST_STEP][ion_index].y_prev;
				IONS[INDEX_LAST_STEP][ion_index].y_prev=IONS[INDEX_LAST_STEP][ion_index].y_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[1]=-IONS[INDEX_LAST_STEP][ion_index].force[1];	
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[1]=-IONS[INDEX_LAST_STEP][ion_index].X_n[1];				
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].y_prev=IONS[INDEX_LAST_STEP][ion_index].y;
				IONS[INDEX_LAST_STEP][ion_index].y=IONS[INDEX_LAST_STEP][ion_index].y_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[1]=IONS[INDEX_LAST_STEP][ion_index].force[1];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[1]=IONS[INDEX_LAST_STEP][ion_index].X_n[1];
			}
		}
		else if(PRM.ION_RECYCLING.substr(0, 2).compare("11")==0){ 
//open boundary, ion recycling			
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].x_next<PRM.MIN_X){				
				IONS[INDEX_LAST_STEP][ion_index].x_prev=1e-12*PRM.SIM_DOMAIN_WIDTH_X+IONS[INDEX_LAST_STEP][ion_index].x;
				IONS[INDEX_LAST_STEP][ion_index].x=1e-12*PRM.SIM_DOMAIN_WIDTH_X+IONS[INDEX_LAST_STEP][ion_index].x_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[0]=IONS[INDEX_LAST_STEP][ion_index].force[0];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[0]=IONS[INDEX_LAST_STEP][ion_index].X_n[0];
			}
			else if(1e12*IONS[INDEX_LAST_STEP][ion_index].x_next>PRM.MAX_X){
				IONS[INDEX_LAST_STEP][ion_index].x_prev=IONS[INDEX_LAST_STEP][ion_index].x-1e-12*PRM.SIM_DOMAIN_WIDTH_X;
				IONS[INDEX_LAST_STEP][ion_index].x=IONS[INDEX_LAST_STEP][ion_index].x_next-1e-12*PRM.SIM_DOMAIN_WIDTH_X;
				IONS[INDEX_LAST_STEP][ion_index].force_old[0]=IONS[INDEX_LAST_STEP][ion_index].force[0];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[0]=IONS[INDEX_LAST_STEP][ion_index].X_n[0];
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].x_prev=IONS[INDEX_LAST_STEP][ion_index].x;
				IONS[INDEX_LAST_STEP][ion_index].x=IONS[INDEX_LAST_STEP][ion_index].x_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[0]=IONS[INDEX_LAST_STEP][ion_index].force[0];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[0]=IONS[INDEX_LAST_STEP][ion_index].X_n[0];
			}
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].y_next<PRM.MIN_Y){
				IONS[INDEX_LAST_STEP][ion_index].y_prev=1e-12*PRM.SIM_DOMAIN_WIDTH_Y+IONS[INDEX_LAST_STEP][ion_index].y;
				IONS[INDEX_LAST_STEP][ion_index].y=1e-12*PRM.SIM_DOMAIN_WIDTH_Y+IONS[INDEX_LAST_STEP][ion_index].y_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[1]=IONS[INDEX_LAST_STEP][ion_index].force[1];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[1]=IONS[INDEX_LAST_STEP][ion_index].X_n[1];
			}
			else if(1e12*IONS[INDEX_LAST_STEP][ion_index].y_next>PRM.MAX_Y){
				IONS[INDEX_LAST_STEP][ion_index].y_prev=IONS[INDEX_LAST_STEP][ion_index].y-1e-12*PRM.SIM_DOMAIN_WIDTH_Y;
				IONS[INDEX_LAST_STEP][ion_index].y=IONS[INDEX_LAST_STEP][ion_index].y_next-1e-12*PRM.SIM_DOMAIN_WIDTH_Y;
				IONS[INDEX_LAST_STEP][ion_index].force_old[1]=IONS[INDEX_LAST_STEP][ion_index].force[1];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[1]=IONS[INDEX_LAST_STEP][ion_index].X_n[1];
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].y_prev=IONS[INDEX_LAST_STEP][ion_index].y;
				IONS[INDEX_LAST_STEP][ion_index].y=IONS[INDEX_LAST_STEP][ion_index].y_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[1]=IONS[INDEX_LAST_STEP][ion_index].force[1];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[1]=IONS[INDEX_LAST_STEP][ion_index].X_n[1];
			}
		}
		else if(PRM.ION_RECYCLING.substr(0, 2).compare("22")==0){ 
//open boundary, ion deletion			
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].x_next<PRM.MIN_X || 1e12*IONS[INDEX_LAST_STEP][ion_index].x_next>PRM.MAX_X){
				ions.back().erase(ions.back().begin()+ion_index);
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].x_prev=IONS[INDEX_LAST_STEP][ion_index].x;
				IONS[INDEX_LAST_STEP][ion_index].x=IONS[INDEX_LAST_STEP][ion_index].x_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[0]=IONS[INDEX_LAST_STEP][ion_index].force[0];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[0]=IONS[INDEX_LAST_STEP][ion_index].X_n[0];
			}
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].y_next<PRM.MIN_Y || 1e12*IONS[INDEX_LAST_STEP][ion_index].y_next>PRM.MAX_Y){
				ions.back().erase(ions.back().begin()+ion_index);
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].y_prev=IONS[INDEX_LAST_STEP][ion_index].y;
				IONS[INDEX_LAST_STEP][ion_index].y=IONS[INDEX_LAST_STEP][ion_index].y_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[1]=IONS[INDEX_LAST_STEP][ion_index].force[1];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[1]=IONS[INDEX_LAST_STEP][ion_index].X_n[1];
			}
		}
		
		if(PRM.ION_RECYCLING.substr(2, 1).compare("0")==0){
//reflective boundary, NO ion recycling				
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].z_next<PRM.MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_index].z_next>PRM.MAX_Z){
				IONS[INDEX_LAST_STEP][ion_index].z_next=IONS[INDEX_LAST_STEP][ion_index].z;
				IONS[INDEX_LAST_STEP][ion_index].z=IONS[INDEX_LAST_STEP][ion_index].z_prev;
				IONS[INDEX_LAST_STEP][ion_index].z_prev=IONS[INDEX_LAST_STEP][ion_index].z_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[2]=-IONS[INDEX_LAST_STEP][ion_index].force[2];	
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[2]=-IONS[INDEX_LAST_STEP][ion_index].X_n[2];				
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].z_prev=IONS[INDEX_LAST_STEP][ion_index].z;
				IONS[INDEX_LAST_STEP][ion_index].z=IONS[INDEX_LAST_STEP][ion_index].z_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[2]=IONS[INDEX_LAST_STEP][ion_index].force[2];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[2]=IONS[INDEX_LAST_STEP][ion_index].X_n[2];
			}		
		}
		else if(PRM.ION_RECYCLING.substr(2, 1).compare("1")==0){ 
//open boundary, ion recycling			
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].z_next<PRM.MIN_Z){
				IONS[INDEX_LAST_STEP][ion_index].z_prev=1e-12*PRM.SIM_DOMAIN_WIDTH_Z+IONS[INDEX_LAST_STEP][ion_index].z;
				IONS[INDEX_LAST_STEP][ion_index].z=1e-12*PRM.SIM_DOMAIN_WIDTH_Z+IONS[INDEX_LAST_STEP][ion_index].z_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[2]=IONS[INDEX_LAST_STEP][ion_index].force[2];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[2]=IONS[INDEX_LAST_STEP][ion_index].X_n[2];
				//~ IONS[INDEX_LAST_STEP][ion_index].side="RIGHT";
			}
			else if(1e12*IONS[INDEX_LAST_STEP][ion_index].z_next>PRM.MAX_Z){
				IONS[INDEX_LAST_STEP][ion_index].z_prev=IONS[INDEX_LAST_STEP][ion_index].z-1e-12*PRM.SIM_DOMAIN_WIDTH_Z;
				IONS[INDEX_LAST_STEP][ion_index].z=IONS[INDEX_LAST_STEP][ion_index].z_next-1e-12*PRM.SIM_DOMAIN_WIDTH_Z;				
				IONS[INDEX_LAST_STEP][ion_index].force_old[2]=IONS[INDEX_LAST_STEP][ion_index].force[2];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[2]=IONS[INDEX_LAST_STEP][ion_index].X_n[2];
				//~ IONS[INDEX_LAST_STEP][ion_index].side="LEFT";
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].z_prev=IONS[INDEX_LAST_STEP][ion_index].z;
				IONS[INDEX_LAST_STEP][ion_index].z=IONS[INDEX_LAST_STEP][ion_index].z_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[2]=IONS[INDEX_LAST_STEP][ion_index].force[2];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[2]=IONS[INDEX_LAST_STEP][ion_index].X_n[2];
			}
		}
		else if(PRM.ION_RECYCLING.substr(2, 1).compare("2")==0){ 
//open boundary, ion deletion	
			if(1e12*IONS[INDEX_LAST_STEP][ion_index].z_next<PRM.MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_index].z_next>PRM.MAX_Z){
				ions.back().erase(ions.back().begin()+ion_index);
			}
			else{
				IONS[INDEX_LAST_STEP][ion_index].z_prev=IONS[INDEX_LAST_STEP][ion_index].z;
				IONS[INDEX_LAST_STEP][ion_index].z=IONS[INDEX_LAST_STEP][ion_index].z_next;
				IONS[INDEX_LAST_STEP][ion_index].force_old[2]=IONS[INDEX_LAST_STEP][ion_index].force[2];
				IONS[INDEX_LAST_STEP][ion_index].X_n_old[2]=IONS[INDEX_LAST_STEP][ion_index].X_n[2];
			}
		}
	}
	
	return;
}

double get_ion_velocity(Ion i){
	
	double vx_2=i.velocity[0]*i.velocity[0];
	double vy_2=i.velocity[1]*i.velocity[1];
	double vz_2=i.velocity[2]*i.velocity[2];
	double velocity=sqrt(vx_2+vy_2+vz_2);
	
	return velocity;
}

double get_flight_length(Ion i){
	
	double ion2_x=i.x_next;
	double ion2_y=i.y_next;
	double ion2_z=i.z_next;			
	apply_periodic_boundary(i.x, i.y, i.z, ion2_x, ion2_y, ion2_z);
	
	
	double dx=ion2_x-i.x;
	double dy=ion2_y-i.y;
	double dz=ion2_z-i.z;
	double flight=sqrt(dx*dx+dy*dy+dz*dz);
	
	return flight;
}

void get_parameters_SR_PMF(Ion i1, Ion i2, vector<double>& SR_PMF_PRM){
	if         ((i1.name.compare("K") == 0) && (i2.name.compare("K")  == 0)) {
		SR_PMF_PRM.push_back(4.007e-18);
		SR_PMF_PRM.push_back(-4.580e+02);
		SR_PMF_PRM.push_back(1.304e+02);
		SR_PMF_PRM.push_back(7.906e-03);
		SR_PMF_PRM.push_back(-7.062e-22);
	} else if (((i1.name.compare("K") == 0) && (i2.name.compare("NA") == 0)) || ((i1.name.compare("NA") == 0) && (i2.name.compare("K") == 0))) {
		SR_PMF_PRM.push_back(1.834e-18);
		SR_PMF_PRM.push_back(-6.248e+02);
		SR_PMF_PRM.push_back(1.665e+02);
		SR_PMF_PRM.push_back(6.982e-03);
		SR_PMF_PRM.push_back(-3.716e-23);
	//} else if (((i1.name.compare("K") == 0) && (i2.name.compare("CA") == 0)) || ((i1.name.compare("CA") == 0) && (i2.name.compare("K") == 0))) {
	} else if (((i1.name.compare("K") == 0) && (i2.name.compare("CL") == 0)) || ((i1.name.compare("CL") == 0) && (i2.name.compare("K") == 0))) {
		SR_PMF_PRM.push_back(3.541e-20);
		SR_PMF_PRM.push_back(-1.872e+02);
		SR_PMF_PRM.push_back(2.283e+02);
		SR_PMF_PRM.push_back(9.885e-03);
		SR_PMF_PRM.push_back(-3.960e-20);
	//} else if (((i1.name.compare("K") == 0) && (i2.name.compare("O")  == 0)) || ((i1.name.compare("O") == 0)  && (i2.name.compare("K") == 0))) {
	//} else if (((i1.name.compare("K") == 0) && (i2.name[0] == "X")) || ((i1.name[0] == "X") && (i2.name.compare("K") == 0))) {
	//} else if (((i1.name.compare("K") == 0) && (i2.name[0] == "Q")) || ((i1.name[0] == "Q") && (i2.name.compare("K") == 0))) {
	} else if  ((i1.name.compare("NA") == 0) && (i2.name.compare("NA") == 0)) {
		SR_PMF_PRM.push_back(2.264e-11);
		SR_PMF_PRM.push_back(-5.244e+01);
		SR_PMF_PRM.push_back(1.503e+01);
		SR_PMF_PRM.push_back(2.695e-03);
		SR_PMF_PRM.push_back(-1.454e-13);
	//} else if (((i1.name.compare("NA") == 0) && (i2.name.compare("CA") == 0)) || ((i1.name.compare("CA") == 0) && (i2.name.compare("NA") == 0))) {
	} else if (((i1.name.compare("NA") == 0) && (i2.name.compare("CL") == 0)) || ((i1.name.compare("CL") == 0) && (i2.name.compare("NA") == 0))) {
		SR_PMF_PRM.push_back(1.909e-18);
		SR_PMF_PRM.push_back(-5.860e+02);
		SR_PMF_PRM.push_back(1.829e+02);
		SR_PMF_PRM.push_back(6.299e-03);
		SR_PMF_PRM.push_back(-6.939e-24);
	//} else if (((i1.name.compare("NA") == 0) && (i2.name.compare("O")  == 0)) || ((i1.name.compare("O") == 0)  && (i2.name.compare("NA") == 0))) {
	//} else if (((i1.name.compare("NA") == 0) && (i2.name[0] == "X")) || ((i1.name[0] == "X") && (i2.name.compare("NA") == 0))) {
	//} else if (((i1.name.compare("NA") == 0) && (i2.name[0] == "Q")) || ((i1.name[0] == "Q") && (i2.name.compare("NA") == 0))) {
	//} else if   (i1.name.compare("CA") == 0) && (i2.name.compare("CA") == 0) {
	//} else if (((i1.name.compare("CA") == 0) && (i2.name.compare("CL") == 0)) || ((i1.name.compare("CL") == 0) && (i2.name.compare("CA") == 0))) {
	//} else if (((i1.name.compare("CA") == 0) && (i2.name.compare("O")  == 0)) || ((i1.name.compare("O") == 0)  && (i2.name.compare("CA") == 0))) {
	//} else if (((i1.name.compare("CA") == 0) && (i2.name[0] == "X")) || ((i1.name[0] == "X") && (i2.name.compare("CA") == 0))) {
	//} else if (((i1.name.compare("CA") == 0) && (i2.name[0] == "Q")) || ((i1.name[0] == "Q") && (i2.name.compare("CA") == 0))) {
	} else if   ((i1.name.compare("CL") == 0) && (i2.name.compare("CL") == 0)) {
		SR_PMF_PRM.push_back(2.472e-16);
		SR_PMF_PRM.push_back(-5.274e+02);
		SR_PMF_PRM.push_back(9.312e+01);
		SR_PMF_PRM.push_back(6.810e-03);
		SR_PMF_PRM.push_back(-3.591e-22);
	//} else if (((i1.name.compare("CL") == 0) && (i2.name.compare("O") == 0)) || ((i1.name.compare("O") == 0) && (i2.name.compare("CL") == 0))) {
	//} else if (((i1.name.compare("CL") == 0) && (i2.name[0] == "X")) || ((i1.name[0] == "X") && (i2.name.compare("CL") == 0))) {
	//} else if (((i1.name.compare("CL") == 0) && (i2.name[0] == "Q")) || ((i1.name[0] == "Q") && (i2.name.compare("CL") == 0))) {
	}
	return;
}
// -END
