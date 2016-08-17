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
#include "sim_structures.h"
#include "ions_properties.h"


void create_simulation_domain(){
	
	
	PRM.MAX_VEL=1.5e4;	// 2e3
	PRM.MAX_FLIGHT=PRM.MAX_VEL*PRM.DELTA_T;
		
	
	
	if(PRM.SEED!=0){
		srand(PRM.SEED);
	}
	else{
		srand (time(NULL));
	}

	PRM.MIN_X=-double(0.5)*PRM.SIM_DOMAIN_WIDTH_X;
	PRM.MAX_X=double(0.5)*PRM.SIM_DOMAIN_WIDTH_X;	
	PRM.MIN_Y=-double(0.5)*PRM.SIM_DOMAIN_WIDTH_Y;
	PRM.MAX_Y=double(0.5)*PRM.SIM_DOMAIN_WIDTH_Y;	
	
	if(PRM.SIM_TYPE.compare("BULK")==0){
		PRM.SIM_DOMAIN_WIDTH_Z=double(2.00)*PRM.CONTROL_CELL_WIDTH+double(2.00)*PRM.OUTER_REGION_WIDTH+double(2.00)*PRM.BATH_WIDTH;
	}
	else{
		PRM.SIM_DOMAIN_WIDTH_Z=double(2.00)*PRM.CONTROL_CELL_WIDTH+double(2.00)*PRM.OUTER_REGION_WIDTH+double(2.00)*PRM.BATH_WIDTH+PRM.MEMBRANE_WIDTH;
	}
	
	PRM.MIN_Z=-double(0.5)*PRM.SIM_DOMAIN_WIDTH_Z;
	PRM.MAX_Z=double(0.5)*PRM.SIM_DOMAIN_WIDTH_Z;
	
	//volumes are in liters
	PRM.SIM_DOMAIN_VOLUME=PRM.SIM_DOMAIN_WIDTH_X*PRM.SIM_DOMAIN_WIDTH_Y*PRM.SIM_DOMAIN_WIDTH_Z;
	PRM.SIM_DOMAIN_VOLUME=1.e-33*PRM.SIM_DOMAIN_VOLUME;
	
	PRM.LEFT_CELL_MIN_Z=PRM.MIN_Z+PRM.OUTER_REGION_WIDTH;
	PRM.LEFT_CELL_MAX_Z=PRM.LEFT_CELL_MIN_Z+PRM.CONTROL_CELL_WIDTH;
	PRM.RIGHT_CELL_MAX_Z=PRM.MAX_Z-PRM.OUTER_REGION_WIDTH;
	PRM.RIGHT_CELL_MIN_Z=PRM.RIGHT_CELL_MAX_Z-PRM.CONTROL_CELL_WIDTH;	
	
	PRM.CONTROL_CELL_VOLUME=PRM.SIM_DOMAIN_WIDTH_X*PRM.SIM_DOMAIN_WIDTH_Y*PRM.CONTROL_CELL_WIDTH;
	PRM.CONTROL_CELL_VOLUME=1e-33*PRM.CONTROL_CELL_VOLUME;
	
	//potential is in V and electric field is in V/m
	PRM.APPLIED_FIELD=PRM.APPLIED_POTENTIAL/PRM.SIM_DOMAIN_WIDTH_Z;
	PRM.APPLIED_FIELD=1e12*PRM.APPLIED_FIELD;
	
	cout << "PRM.APPLIED_POTENTIAL: " << PRM.APPLIED_POTENTIAL <<endl;
	cout << "PRM.APPLIED_FIELD: " << PRM.APPLIED_FIELD <<endl;
	
	/*================================================================================
	the applied potential and the applied field are rescaled in order to have the 
	desired potential drop in the transport region. 
	================================================================================*/
	double potential_drop_region_width=PRM.SIM_DOMAIN_WIDTH_Z-double(2.00)*PRM.OUTER_REGION_WIDTH+PRM.CONTROL_CELL_WIDTH;
	PRM.APPLIED_POTENTIAL=PRM.APPLIED_POTENTIAL*PRM.SIM_DOMAIN_WIDTH_Z/potential_drop_region_width;
	
	PRM.APPLIED_FIELD=1e12*PRM.APPLIED_POTENTIAL/PRM.SIM_DOMAIN_WIDTH_Z;
	
	cout << "Applied potential and electric field automatically rescaled."<<endl;
	cout << "PRM.APPLIED_POTENTIAL: " << PRM.APPLIED_POTENTIAL <<endl;
	cout << "PRM.APPLIED_FIELD: " << PRM.APPLIED_FIELD <<endl;


	
	
	PRM.delta_epsilon=PRM.EPS_MEM-PRM.EPS_W;
	PRM.mean_epsilon=(PRM.EPS_MEM+PRM.EPS_W)/double(2.00);		
	double den_for_K=double(4.00)*M_PI*PRM.mean_epsilon;
	PRM.dielectrics_weight=-PRM.delta_epsilon/den_for_K;
	
	
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
	
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		PRM.ions_to_simulate.push_back(false);
	}
	
	
	
	if(isnormal(PRM.CONC_LEFT_NACL) || isnormal(PRM.CONC_RIGHT_NACL)){
		PRM.ions_to_simulate[3]=true;
	}
	else{
		if(!isnormal(PRM.CONC_LEFT_NACL)){
			PRM.CONC_LEFT_NACL=0;
		}
		else if (!isnormal(PRM.CONC_RIGHT_NACL)){
			PRM.CONC_RIGHT_NACL=0;
		}
	}
	if(isnormal(PRM.CONC_LEFT_KCL) || isnormal(PRM.CONC_RIGHT_KCL)){
		PRM.ions_to_simulate[4]=true;
	}
	else{
		if(!isnormal(PRM.CONC_LEFT_KCL)){
			PRM.CONC_LEFT_KCL=0;
		}
		else if (!isnormal(PRM.CONC_RIGHT_KCL)){
			PRM.CONC_RIGHT_KCL=0;
		}
	}
	if(isnormal(PRM.CONC_LEFT_CACL2) || isnormal(PRM.CONC_RIGHT_CACL2)){
		PRM.ions_to_simulate[9]=true;
	}
	else{
		if(!isnormal(PRM.CONC_LEFT_CACL2)){
			PRM.CONC_LEFT_CACL2=0;
		}
		else if(!isnormal(PRM.CONC_RIGHT_CACL2)){
			PRM.CONC_RIGHT_CACL2=0;
		}
	}
	if(isnormal(PRM.CONC_LEFT_MGCL2) || isnormal(PRM.CONC_RIGHT_MGCL2)){
		PRM.ions_to_simulate[8]=true;
	}
	else{
		if(!isnormal(PRM.CONC_LEFT_MGCL2)){
			PRM.CONC_LEFT_MGCL2=0;
		}
		else if(!isnormal(PRM.CONC_RIGHT_MGCL2)){
			PRM.CONC_RIGHT_MGCL2=0;
		}
	}
	
	
	if(isnormal(PRM.CONC_LEFT_NACL) || isnormal(PRM.CONC_RIGHT_NACL) || isnormal(PRM.CONC_LEFT_KCL) || isnormal(PRM.CONC_RIGHT_KCL) || isnormal(PRM.CONC_LEFT_CACL2) || isnormal(PRM.CONC_RIGHT_CACL2) || isnormal(PRM.CONC_LEFT_MGCL2) || isnormal(PRM.CONC_RIGHT_MGCL2)){
		PRM.ions_to_simulate[17]=true;
	}
	
	if(PRM.BOX1_MAX_Z-PRM.BOX1_MIN_Z>0){
		if(PRM.BOX1_IS1.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX1_IS1.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX1_IS2.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX1_IS2.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX1_IS3.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX1_IS3.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX1_IS4.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX1_IS4.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX1_IS5.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX1_IS5.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX1_IS6.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX1_IS6.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
	}
	
	if(PRM.BOX2_MAX_Z-PRM.BOX2_MIN_Z>0){
		if(PRM.BOX2_IS1.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX2_IS1.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX2_IS2.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX2_IS2.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX2_IS3.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX2_IS3.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX2_IS4.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX2_IS4.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX2_IS5.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX2_IS5.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX2_IS6.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX2_IS6.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
	}
	
	if(PRM.BOX3_MAX_Z-PRM.BOX3_MIN_Z>0){
		if(PRM.BOX3_IS1.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX3_IS1.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX3_IS2.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX3_IS2.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX3_IS3.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX3_IS3.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX3_IS4.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX3_IS4.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX3_IS5.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX3_IS5.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX3_IS6.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX3_IS6.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
	}
	
	if(PRM.BOX4_MAX_Z-PRM.BOX4_MIN_Z>0){
		if(PRM.BOX4_IS1.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX4_IS1.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX4_IS2.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX4_IS2.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX4_IS3.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX4_IS3.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX4_IS4.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX4_IS4.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX4_IS5.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX4_IS5.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX4_IS6.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX4_IS6.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
	}
	
	if(PRM.BOX5_MAX_Z-PRM.BOX5_MIN_Z>0){
		if(PRM.BOX5_IS1.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX5_IS1.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX5_IS2.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX5_IS2.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX5_IS3.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX5_IS3.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX5_IS4.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX5_IS4.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX5_IS5.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX5_IS5.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX5_IS6.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX5_IS6.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
	}
	
	if(PRM.BOX6_MAX_Z-PRM.BOX6_MIN_Z>0){
		if(PRM.BOX6_IS1.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX6_IS1.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX6_IS2.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX6_IS2.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX6_IS3.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX6_IS3.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX6_IS4.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX6_IS4.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX6_IS5.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX6_IS5.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
		if(PRM.BOX6_IS6.compare("NULL")){
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				Ion iii(is);
				if(PRM.BOX6_IS6.compare(iii.name)==0){
					PRM.ions_to_simulate.at(is)=true;
				}
			}
		}
	}
	
	PRM.FIRST_STEP=-(PRM.PREP_STEPS+PRM.HISTORY_SIZE);
	
	
	create_membrane_charges();
	
	
	return;
}

bool initialize_ion_boxes(){

	vector <Ion> aux_ions;
	
	bool CORRECT=true;
	int run=0;
	
	for(int ib=0; ib<ion_boxes.size(); ib++){
		ion_boxes.at(ib).ion_indexes.clear();
		for(int iv=0; iv<ion_boxes.at(ib).is.size(); iv++){
			for(int ino=0; ino<ion_boxes.at(ib).in.at(iv); ino++){
				
				for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
					Ion ion(is);
					
					if(ion_boxes.at(ib).is.at(iv).compare(ion.name)==0){
						
						double min_z=ion_boxes.at(ib).MIN_Z+double(1.10)*ion.pm_radius;
						double width_z=(ion_boxes.at(ib).MAX_Z-ion_boxes.at(ib).MIN_Z)-double(2.00)*double(1.10)*ion.pm_radius;
						
						if(aux_ions.empty()){
							ion.z=min_z+width_z*((double)rand()/((double)(RAND_MAX)+(double)(1)));
							
							double R=PRM.SIM_DOMAIN_WIDTH_X/double(2.00)-double(1.10)*ion.pm_radius;
							if(ion.z>PRM.Z_MOUTH_LEFT && ion.z<PRM.Z_MOUTH_RIGHT){
								int ind_on_lim=ion.z-PRM.Z_MOUTH_LEFT;
								R=limits.at(ind_on_lim)-double(1.10)*ion.pm_radius;
							}
							double r=(double)rand()/((float)RAND_MAX/R);
							double theta=(double)rand()/((double)RAND_MAX/(double(2.00)*M_PI));
							ion.x=r*cos(theta);
							ion.y=r*sin(theta);
						}
						else{
							run=0;
							
							bool goon=true;
							do{
								if(run<10000){
									goon=true;

									ion.z=min_z+width_z*((double)rand()/((double)(RAND_MAX)+(double)(1)));
								
									double R=PRM.SIM_DOMAIN_WIDTH_X/double(2.00)-double(1.10)*ion.pm_radius;
									if(ion.z>PRM.Z_MOUTH_LEFT && ion.z<PRM.Z_MOUTH_RIGHT){
										int ind_on_lim=ion.z-PRM.Z_MOUTH_LEFT;
										R=limits.at(ind_on_lim)-double(1.10)*ion.pm_radius;
									}
									double r=(double)rand()/((float)RAND_MAX/R);
									double theta=(double)rand()/((double)RAND_MAX/(double(2.00)*M_PI));
									ion.x=r*cos(theta);
									ion.y=r*sin(theta);

									for(int j=0; j<aux_ions.size(); j++){
										goon=goon*verifyDistance(ion.x, ion.y, ion.z, double(1.1)*ion.pm_radius, 1e12*aux_ions.at(j).x, 1e12*aux_ions.at(j).y, 1e12*aux_ions.at(j).z, double(1.1)*aux_ions.at(j).pm_radius);
									}
								}
								else{
									CORRECT=false;
									return CORRECT;
								}
								run++;
							}
							while(!goon);
						}
						
						ion.x=1e-12*ion.x;
						ion.y=1e-12*ion.y;
						ion.z=1e-12*ion.z;
						ion.side=-1;
						
						aux_ions.push_back(ion);		
						ion_boxes.at(ib).ion_indexes.push_back(aux_ions.size()-1);
					}
				}
			}
		}
	}
	
	ions.push_back(aux_ions);
	
	
	for(int i=0; i<aux_ions.size(); i++){
		IONS[INDEX_LAST_STEP][i]=aux_ions.at(i);
		NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]++;
	}
	
	return CORRECT;
}



void create_ion_boxes(){
	
	Ion_box box;
	ion_boxes.clear();
	
	if(PRM.BOX1_MAX_Z-PRM.BOX1_MIN_Z>0){
		box.reset_ion_box();
		box.index=1;
		
		//~ box.MIN_X=PRM.BOX1_MIN_X;
		//~ box.MAX_X=PRM.BOX1_MAX_X;
		//~ box.MIN_Y=PRM.BOX1_MIN_Y;
		//~ box.MAX_Y=PRM.BOX1_MAX_Y;
		box.MIN_Z=PRM.BOX1_MIN_Z;
		box.MAX_Z=PRM.BOX1_MAX_Z;
		
		box.TOT_CHARGE=0.0;
		
		if(PRM.BOX1_IS1.compare("NULL")){
			box.is.push_back(PRM.BOX1_IS1);
			box.in.push_back(PRM.BOX1_N1);
			
			box.TOT_CHARGE+=PRM.BOX1_N1*BOX1_ION1_VALENCE;
		}
		if(PRM.BOX1_IS2.compare("NULL")){
			box.is.push_back(PRM.BOX1_IS2);
			box.in.push_back(PRM.BOX1_N2);
			
			box.TOT_CHARGE+=PRM.BOX1_N2*BOX1_ION2_VALENCE;
		}
		if(PRM.BOX1_IS3.compare("NULL")){
			box.is.push_back(PRM.BOX1_IS3);
			box.in.push_back(PRM.BOX1_N3);
			
			box.TOT_CHARGE+=PRM.BOX1_N3*BOX1_ION3_VALENCE;
		}
		if(PRM.BOX1_IS4.compare("NULL")){
			box.is.push_back(PRM.BOX1_IS4);
			box.in.push_back(PRM.BOX1_N4);
			
			box.TOT_CHARGE+=PRM.BOX1_N4*BOX1_ION4_VALENCE;
		}
		if(PRM.BOX1_IS5.compare("NULL")){
			box.is.push_back(PRM.BOX1_IS5);
			box.in.push_back(PRM.BOX1_N5);
			
			box.TOT_CHARGE+=PRM.BOX1_N5*BOX1_ION5_VALENCE;
		}
		if(PRM.BOX1_IS6.compare("NULL")){
			box.is.push_back(PRM.BOX1_IS6);
			box.in.push_back(PRM.BOX1_N6);
			
			box.TOT_CHARGE+=PRM.BOX1_N6*BOX1_ION6_VALENCE;
		}
		ion_boxes.push_back(box);
	}
	if(PRM.BOX2_MAX_Z-PRM.BOX2_MIN_Z>0){
		box.reset_ion_box();
		box.index=2;
		
		//~ box.MIN_X=PRM.BOX2_MIN_X;
		//~ box.MAX_X=PRM.BOX2_MAX_X;
		//~ box.MIN_Y=PRM.BOX2_MIN_Y;
		//~ box.MAX_Y=PRM.BOX2_MAX_Y;
		box.MIN_Z=PRM.BOX2_MIN_Z;
		box.MAX_Z=PRM.BOX2_MAX_Z;
		
		box.TOT_CHARGE=0.0;
		
		if(PRM.BOX2_IS1.compare("NULL")){
			box.is.push_back(PRM.BOX2_IS1);
			box.in.push_back(PRM.BOX2_N1);
			
			box.TOT_CHARGE+=PRM.BOX2_N1*BOX2_ION1_VALENCE;
		}
		if(PRM.BOX2_IS2.compare("NULL")){
			box.is.push_back(PRM.BOX2_IS2);
			box.in.push_back(PRM.BOX2_N2);
			
			box.TOT_CHARGE+=PRM.BOX2_N2*BOX2_ION2_VALENCE;
		}
		if(PRM.BOX2_IS3.compare("NULL")){
			box.is.push_back(PRM.BOX2_IS3);
			box.in.push_back(PRM.BOX2_N3);
			
			box.TOT_CHARGE+=PRM.BOX2_N3*BOX2_ION3_VALENCE;
		}
		if(PRM.BOX2_IS4.compare("NULL")){
			box.is.push_back(PRM.BOX2_IS4);
			box.in.push_back(PRM.BOX2_N4);
			
			box.TOT_CHARGE+=PRM.BOX2_N4*BOX2_ION4_VALENCE;
		}
		if(PRM.BOX2_IS5.compare("NULL")){
			box.is.push_back(PRM.BOX2_IS5);
			box.in.push_back(PRM.BOX2_N5);
			
			box.TOT_CHARGE+=PRM.BOX2_N5*BOX2_ION5_VALENCE;
		}
		if(PRM.BOX2_IS6.compare("NULL")){
			box.is.push_back(PRM.BOX2_IS6);
			box.in.push_back(PRM.BOX2_N6);
			
			box.TOT_CHARGE+=PRM.BOX2_N6*BOX2_ION6_VALENCE;
		}
		ion_boxes.push_back(box);
	}
	if(PRM.BOX3_MAX_Z-PRM.BOX3_MIN_Z>0){
		box.reset_ion_box();
		box.index=3;
		
		//~ box.MIN_X=PRM.BOX3_MIN_X;
		//~ box.MAX_X=PRM.BOX3_MAX_X;
		//~ box.MIN_Y=PRM.BOX3_MIN_Y;
		//~ box.MAX_Y=PRM.BOX3_MAX_Y;
		box.MIN_Z=PRM.BOX3_MIN_Z;
		box.MAX_Z=PRM.BOX3_MAX_Z;
		
		box.TOT_CHARGE=0.0;
		
		if(PRM.BOX3_IS1.compare("NULL")){
			box.is.push_back(PRM.BOX3_IS1);
			box.in.push_back(PRM.BOX3_N1);
			
			box.TOT_CHARGE+=PRM.BOX3_N1*BOX3_ION1_VALENCE;
		}
		if(PRM.BOX3_IS2.compare("NULL")){
			box.is.push_back(PRM.BOX3_IS2);
			box.in.push_back(PRM.BOX3_N2);
			
			box.TOT_CHARGE+=PRM.BOX3_N2*BOX3_ION2_VALENCE;
		}
		if(PRM.BOX3_IS3.compare("NULL")){
			box.is.push_back(PRM.BOX3_IS3);
			box.in.push_back(PRM.BOX3_N3);
			
			box.TOT_CHARGE+=PRM.BOX3_N3*BOX3_ION3_VALENCE;
		}
		if(PRM.BOX3_IS4.compare("NULL")){
			box.is.push_back(PRM.BOX3_IS4);
			box.in.push_back(PRM.BOX3_N4);
			
			box.TOT_CHARGE+=PRM.BOX3_N4*BOX3_ION4_VALENCE;
		}
		if(PRM.BOX3_IS5.compare("NULL")){
			box.is.push_back(PRM.BOX3_IS5);
			box.in.push_back(PRM.BOX3_N5);
			
			box.TOT_CHARGE+=PRM.BOX3_N5*BOX3_ION5_VALENCE;
		}
		if(PRM.BOX3_IS6.compare("NULL")){
			box.is.push_back(PRM.BOX3_IS6);
			box.in.push_back(PRM.BOX3_N6);
			
			box.TOT_CHARGE+=PRM.BOX3_N6*BOX3_ION6_VALENCE;
		}
		ion_boxes.push_back(box);
	}
	if(PRM.BOX4_MAX_Z-PRM.BOX4_MIN_Z>0){
		box.reset_ion_box();
		box.index=4;
		
		//~ box.MIN_X=PRM.BOX4_MIN_X;
		//~ box.MAX_X=PRM.BOX4_MAX_X;
		//~ box.MIN_Y=PRM.BOX4_MIN_Y;
		//~ box.MAX_Y=PRM.BOX4_MAX_Y;
		box.MIN_Z=PRM.BOX4_MIN_Z;
		box.MAX_Z=PRM.BOX4_MAX_Z;
		
		box.TOT_CHARGE=0.0;
		
		if(PRM.BOX4_IS1.compare("NULL")){
			box.is.push_back(PRM.BOX4_IS1);
			box.in.push_back(PRM.BOX4_N1);
			
			box.TOT_CHARGE+=PRM.BOX4_N1*BOX4_ION1_VALENCE;
		}
		if(PRM.BOX4_IS2.compare("NULL")){
			box.is.push_back(PRM.BOX4_IS2);
			box.in.push_back(PRM.BOX4_N2);
			
			box.TOT_CHARGE+=PRM.BOX4_N2*BOX4_ION2_VALENCE;
		}
		if(PRM.BOX4_IS3.compare("NULL")){
			box.is.push_back(PRM.BOX4_IS3);
			box.in.push_back(PRM.BOX4_N3);
			
			box.TOT_CHARGE+=PRM.BOX4_N3*BOX4_ION3_VALENCE;
		}
		if(PRM.BOX4_IS4.compare("NULL")){
			box.is.push_back(PRM.BOX4_IS4);
			box.in.push_back(PRM.BOX4_N4);
			
			box.TOT_CHARGE+=PRM.BOX4_N4*BOX4_ION4_VALENCE;
		}
		if(PRM.BOX4_IS5.compare("NULL")){
			box.is.push_back(PRM.BOX4_IS5);
			box.in.push_back(PRM.BOX4_N5);
			
			box.TOT_CHARGE+=PRM.BOX4_N5*BOX4_ION5_VALENCE;
		}
		if(PRM.BOX4_IS6.compare("NULL")){
			box.is.push_back(PRM.BOX4_IS6);
			box.in.push_back(PRM.BOX4_N6);
			
			box.TOT_CHARGE+=PRM.BOX4_N6*BOX4_ION6_VALENCE;
		}
		ion_boxes.push_back(box);
	}
	if(PRM.BOX5_MAX_Z-PRM.BOX5_MIN_Z>0){
		box.reset_ion_box();
		box.index=5;
		
		//~ box.MIN_X=PRM.BOX5_MIN_X;
		//~ box.MAX_X=PRM.BOX5_MAX_X;
		//~ box.MIN_Y=PRM.BOX5_MIN_Y;
		//~ box.MAX_Y=PRM.BOX5_MAX_Y;
		box.MIN_Z=PRM.BOX5_MIN_Z;
		box.MAX_Z=PRM.BOX5_MAX_Z;
		
		box.TOT_CHARGE=0.0;
		
		if(PRM.BOX5_IS1.compare("NULL")){
			box.is.push_back(PRM.BOX5_IS1);
			box.in.push_back(PRM.BOX5_N1);
			
			box.TOT_CHARGE+=PRM.BOX5_N1*BOX5_ION1_VALENCE;
		}
		if(PRM.BOX5_IS2.compare("NULL")){
			box.is.push_back(PRM.BOX5_IS2);
			box.in.push_back(PRM.BOX5_N2);
			
			box.TOT_CHARGE+=PRM.BOX5_N2*BOX5_ION2_VALENCE;
		}
		if(PRM.BOX5_IS3.compare("NULL")){
			box.is.push_back(PRM.BOX5_IS3);
			box.in.push_back(PRM.BOX5_N3);
			
			box.TOT_CHARGE+=PRM.BOX5_N3*BOX5_ION3_VALENCE;
		}
		if(PRM.BOX5_IS4.compare("NULL")){
			box.is.push_back(PRM.BOX5_IS4);
			box.in.push_back(PRM.BOX5_N4);
			
			box.TOT_CHARGE+=PRM.BOX5_N4*BOX5_ION4_VALENCE;
		}
		if(PRM.BOX5_IS5.compare("NULL")){
			box.is.push_back(PRM.BOX5_IS5);
			box.in.push_back(PRM.BOX5_N5);
			
			box.TOT_CHARGE+=PRM.BOX5_N5*BOX5_ION5_VALENCE;
		}
		if(PRM.BOX5_IS6.compare("NULL")){
			box.is.push_back(PRM.BOX5_IS6);
			box.in.push_back(PRM.BOX5_N6);
			
			box.TOT_CHARGE+=PRM.BOX5_N6*BOX5_ION6_VALENCE;
		}
		ion_boxes.push_back(box);
	}
	if(PRM.BOX6_MAX_Z-PRM.BOX6_MIN_Z>0){
		box.reset_ion_box();
		box.index=6;
		
		//~ box.MIN_X=PRM.BOX6_MIN_X;
		//~ box.MAX_X=PRM.BOX6_MAX_X;
		//~ box.MIN_Y=PRM.BOX6_MIN_Y;
		//~ box.MAX_Y=PRM.BOX6_MAX_Y;
		box.MIN_Z=PRM.BOX6_MIN_Z;
		box.MAX_Z=PRM.BOX6_MAX_Z;
		
		box.TOT_CHARGE=0.0;
		
		if(PRM.BOX6_IS1.compare("NULL")){
			box.is.push_back(PRM.BOX6_IS1);
			box.in.push_back(PRM.BOX6_N1);
			
			box.TOT_CHARGE+=PRM.BOX6_N1*BOX6_ION1_VALENCE;
		}
		if(PRM.BOX6_IS2.compare("NULL")){
			box.is.push_back(PRM.BOX6_IS2);
			box.in.push_back(PRM.BOX6_N2);
			
			box.TOT_CHARGE+=PRM.BOX6_N2*BOX6_ION2_VALENCE;
		}
		if(PRM.BOX6_IS3.compare("NULL")){
			box.is.push_back(PRM.BOX6_IS3);
			box.in.push_back(PRM.BOX6_N3);
			
			box.TOT_CHARGE+=PRM.BOX6_N3*BOX6_ION3_VALENCE;
		}
		if(PRM.BOX6_IS4.compare("NULL")){
			box.is.push_back(PRM.BOX6_IS4);
			box.in.push_back(PRM.BOX6_N4);
			
			box.TOT_CHARGE+=PRM.BOX6_N4*BOX6_ION4_VALENCE;
		}
		if(PRM.BOX6_IS5.compare("NULL")){
			box.is.push_back(PRM.BOX6_IS5);
			box.in.push_back(PRM.BOX6_N5);
			
			box.TOT_CHARGE+=PRM.BOX6_N5*BOX6_ION5_VALENCE;
		}
		if(PRM.BOX6_IS6.compare("NULL")){
			box.is.push_back(PRM.BOX6_IS6);
			box.in.push_back(PRM.BOX6_N6);
			
			box.TOT_CHARGE+=PRM.BOX6_N6*BOX6_ION6_VALENCE;
		}
		ion_boxes.push_back(box);
	}
	
	return;
}


void create_membrane_charges(){
	
	membrane_charges.clear();

	cout << "a" << endl;
	
	Charge c;
	
	double angle;
	
	for(int iring=0; iring < PRM.charge_ring_z.size(); iring++){
		angle=double(2.00)*M_PI/double(PRM.charge_ring_n.at(iring));
		for(int ichg=0; ichg<PRM.charge_ring_n.at(iring); ichg++){	
			c.reset_charge();
			c.x=PRM.charge_ring_r.at(iring)*cos(double(ichg)*angle);
			c.y=PRM.charge_ring_r.at(iring)*sin(double(ichg)*angle);
			c.z=PRM.charge_ring_z.at(iring);
			c.charge=Q*PRM.charge_ring_q.at(iring);
			c.valence=PRM.charge_ring_q.at(iring);
			c.DW_charge=c.charge/PRM.EPS_MEM;
			c.DW_valence=c.valence/PRM.EPS_MEM;
			if(fabs(c.x)<1e-3)c.x=0.00;
			if(fabs(c.y)<1e-3)c.y=0.00;
			if(fabs(c.z)<1e-3)c.z=0.00;		
			c.x=1e-12*c.x;
			c.y=1e-12*c.y;
			c.z=1e-12*c.z;
			
			c.radius=1.5e-10;
			
			membrane_charges.push_back(c);
		}

	}
	
	cout << "b" << endl;
	
	if(fileExists(PRM.FIXED_CHARGES_FILE)){

		string buffer="";
		
		ifstream fin; 
		char *tfn = new char[PRM.FIXED_CHARGES_FILE.length()+1];
		strcpy(tfn, PRM.FIXED_CHARGES_FILE.c_str());     
		fin.open(tfn, ifstream::in);  

		if(!fin){
			cout <<"ERROR: it is not possible to open the file: " << PRM.FIXED_CHARGES_FILE << endl;
			exit(4);
		}
		
		while(!fin.eof()){
			if(getline(fin, buffer)){
				if(buffer.substr(0,1).compare("#")!=0 && countTokens(buffer, " ")==4){
					c.reset_charge();
					
					c.x=1e-12*atof(getTokenbyNumber(buffer, " ", 1).c_str());
					c.y=1e-12*atof(getTokenbyNumber(buffer, " ", 2).c_str());
					c.z=1e-12*atof(getTokenbyNumber(buffer, " ", 3).c_str());
					c.charge=Q*atof(getTokenbyNumber(buffer, " ", 4).c_str());
					c.valence=atof(getTokenbyNumber(buffer, " ", 4).c_str());
					c.DW_charge=c.charge/PRM.EPS_MEM;
					c.DW_valence=c.valence/PRM.EPS_MEM;
					
					c.radius=0;
					
					membrane_charges.push_back(c);
				}
			}
		}
	
		cout << "NUMBER OF FIXED CHARGES: " << membrane_charges.size() << endl; 
		
		double total_pos_charge=0;
		double total_neg_charge=0;
		double total_charge=0;
		for(int i=0; i<membrane_charges.size(); i++){
			if(membrane_charges.at(i).valence>0){
				total_pos_charge+=membrane_charges.at(i).valence;
			}
			else if(membrane_charges.at(i).valence<0){
				total_neg_charge+=membrane_charges.at(i).valence;
			}
		}
		
		total_charge=total_pos_charge+total_neg_charge;
		
		cout << "TOTAL POS CHARGE (fixed): " << total_pos_charge << endl; 
		cout << "TOTAL NEG CHARGE (fixed): " << total_neg_charge << endl; 
		cout << "TOTAL CHARGE (fixed): " << total_charge << endl; 
		
		fin.close();
		
	}
	
	cout << "c" << endl;
		
	return;
}



















