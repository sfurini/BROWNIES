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
#include "physics_functions.h"
#include "ions_functions.h"
#include "PACO.h"
#include "classes.h"
#include "sim_structures.h"



bool verifyDistance(double ax, double ay, double az, double ar, double bx, double by, double bz, double br){
	double d_m=ar+br;
	double d=sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz));
	if(d>=d_m){
		return true;
	}
	else{
		return false;
	}  
}

void compute_c_base(){
	
	for(int i=0; i<NUM_OF_SURFACES; i++){
		C_BASE[i]=0.00;
	}
	
	
	for(int surface_index=0; surface_index<NUM_OF_SURFACES; surface_index++){

// external field		
		double elem_c=1e12*SURFACES[surface_index].area*SURFACES[surface_index].normal[2]*PRM.dielectrics_weight*PRM.EPS_W*PRM.APPLIED_FIELD;		// V/m
						
// membrane charges		
		if(!membrane_charges.empty()){
			for(int charge_index=0; charge_index<membrane_charges.size(); charge_index++){
				double distance=1e12*get_distance(SURFACES[surface_index].center_x, SURFACES[surface_index].center_y, SURFACES[surface_index].center_z, membrane_charges.at(charge_index).x, membrane_charges.at(charge_index).y, membrane_charges.at(charge_index).z);		
				if(distance<20000){
					double temp1=PRM.dielectrics_weight*membrane_charges.at(charge_index).DW_valence;	
					double field=temp1*FIELD_C[int(distance)]/COULOMB_K;
			
					double cs_x=(SURFACES[surface_index].center_x-membrane_charges.at(charge_index).x)/(1e-12*distance);
					double cs_y=(SURFACES[surface_index].center_y-membrane_charges.at(charge_index).y)/(1e-12*distance);
					double cs_z=(SURFACES[surface_index].center_z-membrane_charges.at(charge_index).z)/(1e-12*distance);
					
					double scalar=cs_x*SURFACES[surface_index].normal[0]+cs_y*SURFACES[surface_index].normal[1]+cs_z*SURFACES[surface_index].normal[2];
					elem_c+=1e24*SURFACES[surface_index].area*field*scalar;
				}
				else{
					//cut-off
				}
			}
		}
		
		C_BASE[surface_index]=elem_c;	
	}
	
	return;
}


void compute_force_on_ions(){
	
	//~ cout << "PRM.EPS_W:\t" << PRM.EPS_W <<endl;
	//~ cout << "PRM.EPS_MEM:\t" << PRM.EPS_MEM <<endl;
	//~ cout << "NUM_OF_SURFACES:\t" << NUM_OF_SURFACES <<endl;
	
	// compute induced charges		
	if(PRM.EPS_W!=PRM.EPS_MEM){
		
		for(int surface_index=0; surface_index<NUM_OF_SURFACES; surface_index++){
			
			// external field and membrane charges			
			double elem_c=C_BASE[surface_index];
			double elem_c_MD_MAP=0;
		
			if(PRM.MD_MAP_MIN_Z<PRM.MD_MAP_MAX_Z){
				double elem_c_MD_MAP=1e12*SURFACES[surface_index].area*SURFACES[surface_index].normal[2]*PRM.dielectrics_weight*PRM.EPS_W*PRM.APPLIED_FIELD;
			}
	
			// mobile ions	
			for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){
				double distance=1e12*get_distance(SURFACES[surface_index].center_x, SURFACES[surface_index].center_y, SURFACES[surface_index].center_z, IONS[INDEX_LAST_STEP][ion_index].x, IONS[INDEX_LAST_STEP][ion_index].y, IONS[INDEX_LAST_STEP][ion_index].z);		
				if(distance<20000){
					double temp1=PRM.dielectrics_weight*IONS[INDEX_LAST_STEP][ion_index].DW_valence;	
					double field=temp1*FIELD_C[int(distance)]/COULOMB_K;
			
					double cs_x=(SURFACES[surface_index].center_x-IONS[INDEX_LAST_STEP][ion_index].x)/(1e-12*distance);
					double cs_y=(SURFACES[surface_index].center_y-IONS[INDEX_LAST_STEP][ion_index].y)/(1e-12*distance);
					double cs_z=(SURFACES[surface_index].center_z-IONS[INDEX_LAST_STEP][ion_index].z)/(1e-12*distance);
					
					double scalar=cs_x*SURFACES[surface_index].normal[0]+cs_y*SURFACES[surface_index].normal[1]+cs_z*SURFACES[surface_index].normal[2];
					
					elem_c+=1e24*SURFACES[surface_index].area*field*scalar;
					
					if(PRM.MD_MAP_MIN_Z<PRM.MD_MAP_MAX_Z){
						if(1e12*IONS[INDEX_LAST_STEP][ion_index].z<PRM.MD_MAP_MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_index].z>PRM.MD_MAP_MAX_Z){
							elem_c_MD_MAP+=1e24*SURFACES[surface_index].area*field*scalar;
						}
					}
				}					
			}

			vector_c[surface_index]=elem_c;
			vector_c_MD_map[surface_index]=elem_c_MD_MAP;
		} 

		if(NUM_OF_SURFACES>0){
			// solve the ICC equation 					
			dgemv_(&_N, &matrix_size, &matrix_size, &_double_one, inverse_a, &matrix_size, vector_c, &_int_one, &_double_zero, vector_h, &_int_one);	
			if(PRM.MD_MAP_MIN_Z<PRM.MD_MAP_MAX_Z){
				dgemv_(&_N, &matrix_size, &matrix_size, &_double_one, inverse_a, &matrix_size, vector_c_MD_map, &_int_one, &_double_zero, vector_h_MD_map, &_int_one);	
			}
		}		
	}
	//~ cout << "PRM.EPS_W:\t" << PRM.EPS_W <<endl;
	//~ cout << "PRM.EPS_MEM:\t" << PRM.EPS_MEM <<endl;
	//~ cout << "NUM_OF_SURFACES:\t" << NUM_OF_SURFACES <<endl;
	
	
	//simone te devi modificare questa parte
	for(int ion_index_1=0; ion_index_1<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index_1++){

		// external field
		// all ions feel the external field, right?
		IONS[INDEX_LAST_STEP][ion_index_1].force[0]=0.00;
		IONS[INDEX_LAST_STEP][ion_index_1].force[1]=0.00;
		IONS[INDEX_LAST_STEP][ion_index_1].force[2]=Q*IONS[INDEX_LAST_STEP][ion_index_1].valence*PRM.APPLIED_FIELD;
			
		// 1- only ions outside of the map feel the protein charges, right?
		// 2- only ions outside of the map feel the induced charges, right?
		// 3- only ions outside of the map feel the SR boundary repulsion, right?
		// 4- only ions outside of the map feel the SR boxes repulsion, right? actually it's not, since the boxes are in the MD map region, but...
		
		if(1e12*IONS[INDEX_LAST_STEP][ion_index_1].z<PRM.MD_MAP_MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_index_1].z>=PRM.MD_MAP_MAX_Z){
			
			// 1- membrane charges
			if(!membrane_charges.empty()){
				for(int charge_index=0; charge_index<membrane_charges.size(); charge_index++){
					double charge_x=membrane_charges.at(charge_index).x;
					double charge_y=membrane_charges.at(charge_index).y;
					double charge_z=membrane_charges.at(charge_index).z;			
					apply_periodic_boundary(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, charge_x, charge_y, charge_z);
					double distance=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, charge_x, charge_y, charge_z);			
					
					if(distance<20000){
						double force_C=0.00;
						double force_SR=0.00;
						double dx=0.00;
						double dy=0.00;
						double dz=0.00;
						
						dx=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].x-charge_x)/distance;
						dy=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].y-charge_y)/distance;
						dz=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].z-charge_z)/distance;
						
						force_C=IONS[INDEX_LAST_STEP][ion_index_1].valence*membrane_charges.at(charge_index).DW_valence*FORCE_C[int(distance)];
						
						if(distance<4000){
							force_SR=FORCE_SR.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).at(NUM_OF_IONIC_SPECIES).at(int(distance));
						}
						
						IONS[INDEX_LAST_STEP][ion_index_1].force[0]+=(force_C+force_SR)*dx;
						IONS[INDEX_LAST_STEP][ion_index_1].force[1]+=(force_C+force_SR)*dy;
						IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=(force_C+force_SR)*dz;
					}	
					else{
						// cut-off
					}		
				}		
			}
			
			// 2- induced charges
			if(PRM.EPS_W!=PRM.EPS_MEM){
				for(int surface_index=0; surface_index<NUM_OF_SURFACES; surface_index++){		
					double charge_x=SURFACES[surface_index].center_x;
					double charge_y=SURFACES[surface_index].center_y;
					double charge_z=SURFACES[surface_index].center_z;				
					apply_periodic_boundary(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, charge_x, charge_y, charge_z);
					double distance=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, charge_x, charge_y, charge_z);	

					if(distance<20000){
						double force_C=0.00;
						double force_SR=0.00;
						double dx=0.00;
						double dy=0.00;
						double dz=0.00;
						
						dx=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].x-charge_x)/distance;
						dy=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].y-charge_y)/distance;
						dz=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].z-charge_z)/distance;
						
						double surf_charge_valence=(vector_h[surface_index]*surfaces.at(surface_index).area)/Q;
						force_C=IONS[INDEX_LAST_STEP][ion_index_1].valence*surf_charge_valence*FORCE_C[int(distance)];
						
						IONS[INDEX_LAST_STEP][ion_index_1].force[0]+=force_C*dx;
						IONS[INDEX_LAST_STEP][ion_index_1].force[1]+=force_C*dy;
						IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=force_C*dz;
					}	
					else{
						// cut-off
					}			
				}
			}
	
			// 3- SR boundary repulsion		
			if(NUM_OF_SURFACES!=0){
				int firstIndex=(int)(1e12*IONS[INDEX_LAST_STEP][ion_index_1].z-(PRM.Z_BEGIN_SR_MAP))/PRM.DELTA_SR_MAP;
				if(firstIndex>=0 && firstIndex<SR_dielectric_boundary.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).size()){	

					double dx=1e12*IONS[INDEX_LAST_STEP][ion_index_1].x;
					double dy=1e12*IONS[INDEX_LAST_STEP][ion_index_1].y;
					double dista=sqrt(dx*dx+dy*dy);
					int secondIndex=(int)(dista/PRM.DELTA_SR_MAP);
					
					if(secondIndex>=0 && secondIndex<SR_dielectric_boundary.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).at(firstIndex).size()){
						IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=SR_dielectric_boundary.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).at(firstIndex).at(secondIndex).at(2);
						
						if(fabs(dista)>1e-8){
							double angle = atan2(dy,dx);
							double p0 = SR_dielectric_boundary.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).at(firstIndex).at(secondIndex).at(0);
							
							double n0 = p0*cos(angle);
							double n1 = p0*sin(angle);
							
							IONS[INDEX_LAST_STEP][ion_index_1].force[0]+=n0;
							IONS[INDEX_LAST_STEP][ion_index_1].force[1]+=n1;
						}
					}
				}
			}
	
		
			// 4- ion boxes SR repulsion		
			if(!SR_ion_boxes.empty()){
				if(1e12*IONS[INDEX_LAST_STEP][ion_index_1].z>=PRM.Z_MOUTH_LEFT && 1e12*IONS[INDEX_LAST_STEP][ion_index_1].z<=PRM.Z_MOUTH_RIGHT){
					int index=1e12*IONS[INDEX_LAST_STEP][ion_index_1].z-PRM.Z_MOUTH_LEFT;
					IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=fabs(IONS[INDEX_LAST_STEP][ion_index_1].DW_valence)*SR_ion_boxes.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).at(index);
					if(IONS[INDEX_LAST_STEP][ion_index_1].DW_valence==0){
						IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=SR_ion_boxes.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).at(index);
					}
					
				}
				else{
					
				}
			}		
		}
		else{

			//IONS IN THE MD POTENTIAL MAP
			// 1- induced charges from external field, membrane charges
			if(NUM_OF_SURFACES>0){
				if(PRM.EPS_W!=PRM.EPS_MEM){
					for(int surface_index=0; surface_index<NUM_OF_SURFACES; surface_index++){		
						double charge_x=SURFACES[surface_index].center_x;
						double charge_y=SURFACES[surface_index].center_y;
						double charge_z=SURFACES[surface_index].center_z;				
						apply_periodic_boundary(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, charge_x, charge_y, charge_z);
						double distance=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, charge_x, charge_y, charge_z);	

						if(distance<20000){
							double force_C=0.00;
							double force_SR=0.00;
							double dx=0.00;
							double dy=0.00;
							double dz=0.00;
							
							dx=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].x-charge_x)/distance;
							dy=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].y-charge_y)/distance;
							dz=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].z-charge_z)/distance;
							
							double surf_charge_valence=(vector_h_MD_map[surface_index]*surfaces.at(surface_index).area)/Q;
							force_C=IONS[INDEX_LAST_STEP][ion_index_1].valence*surf_charge_valence*FORCE_C[int(distance)];
							
							IONS[INDEX_LAST_STEP][ion_index_1].force[0]+=force_C*dx;
							IONS[INDEX_LAST_STEP][ion_index_1].force[1]+=force_C*dy;
							IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=force_C*dz;
						}	
						else{
							// cut-off
						}			
					}
				}
			}
			//################################################################################################################################################
			//################################################################################################################################################
			//################################################################################################################################################
			//
			//
			//		  qui ci va la forza generata dalla mappa MD sullo ione IONS[INDEX_LAST_STEP][ion_index_1]
			//
			//		le interazioni ione-ione sono considerate nelle righe sotto
			//
			
			//################################################################################################################################################
			//################################################################################################################################################
			//################################################################################################################################################
		}
	}
			
		
	// ion-ion interaction
	for(int ion_index_1=0; ion_index_1<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]-1; ion_index_1++){
		for(int ion_index_2=ion_index_1+1; ion_index_2<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index_2++){
			
			//--------------------------------------------------------------------------------------------------
			// simone questa controllamela bene. l'unica volta che due ioni NON si sentono Coulombianamente e' solo quando sono entrambi nella mappa MD, giusto?
			// quindi se ce n'e' almeno uno fuori dalla mappa la forza la sentono tutti e due
			if(1e12*IONS[INDEX_LAST_STEP][ion_index_1].z<PRM.MD_MAP_MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_index_1].z>PRM.MD_MAP_MAX_Z || 1e12*IONS[INDEX_LAST_STEP][ion_index_2].z<PRM.MD_MAP_MIN_Z || 1e12*IONS[INDEX_LAST_STEP][ion_index_2].z>PRM.MD_MAP_MAX_Z){
			
				double ion2_x=IONS[INDEX_LAST_STEP][ion_index_2].x;
				double ion2_y=IONS[INDEX_LAST_STEP][ion_index_2].y;
				double ion2_z=IONS[INDEX_LAST_STEP][ion_index_2].z;			
				apply_periodic_boundary(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, ion2_x, ion2_y, ion2_z);
				double distance=1e12*get_distance(IONS[INDEX_LAST_STEP][ion_index_1].x, IONS[INDEX_LAST_STEP][ion_index_1].y, IONS[INDEX_LAST_STEP][ion_index_1].z, ion2_x, ion2_y, ion2_z);			
				
				if(distance<20000){
					double force_C=0.00;
					double force_SR=0.00;
					double force_O=0.00;
					double dx=0.00;
					double dy=0.00;
					double dz=0.00;
					
					dx=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].x-ion2_x)/distance;
					dy=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].y-ion2_y)/distance;
					dz=1e12*(IONS[INDEX_LAST_STEP][ion_index_1].z-ion2_z)/distance;
					
					force_C=IONS[INDEX_LAST_STEP][ion_index_1].valence*IONS[INDEX_LAST_STEP][ion_index_2].DW_valence*FORCE_C[int(distance)];
					
					if(distance<4000){
						force_SR=FORCE_SR.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).at(IONS[INDEX_LAST_STEP][ion_index_2].kind).at(int(distance));
						force_O=FORCE_OTHER.at(IONS[INDEX_LAST_STEP][ion_index_1].kind).at(IONS[INDEX_LAST_STEP][ion_index_2].kind).at(int(distance));
					}
					
					IONS[INDEX_LAST_STEP][ion_index_1].force[0]+=(force_C+force_SR+force_O)*dx;
					IONS[INDEX_LAST_STEP][ion_index_1].force[1]+=(force_C+force_SR+force_O)*dy;
					IONS[INDEX_LAST_STEP][ion_index_1].force[2]+=(force_C+force_SR+force_O)*dz;
					
					IONS[INDEX_LAST_STEP][ion_index_2].force[0]-=(force_C+force_SR+force_O)*dx;
					IONS[INDEX_LAST_STEP][ion_index_2].force[1]-=(force_C+force_SR+force_O)*dy;
					IONS[INDEX_LAST_STEP][ion_index_2].force[2]-=(force_C+force_SR+force_O)*dz;
				}	
				else{
					// cut-off
				}	
			}					
		}
		
	}
	
	
	
	
	//compute the potential along the axis
	if(PRM.potential!=0){
		
		if(STEPS[INDEX_LAST_STEP]>=0){
		
			for(int i=0; i<2001; i++){
				
				double ics=0;
				double ipsilon=0;
				double zeta=POINTS_ON_AXIS[i];
				
				//external field
				POTENTIALS_ON_AXIS[INDEX_LAST_STEP][i]=PRM.APPLIED_POTENTIAL*(PRM.MAX_Z-1e12*POINTS_ON_AXIS[i])/(PRM.SIM_DOMAIN_WIDTH_Z);
				
				//membrane charges
				if(!membrane_charges.empty()){
					for(int charge_index=0; charge_index<membrane_charges.size(); charge_index++){
						
						double charge_x=membrane_charges.at(charge_index).x;
						double charge_y=membrane_charges.at(charge_index).y;
						double charge_z=membrane_charges.at(charge_index).z;				
						apply_periodic_boundary(ics, ipsilon, zeta, charge_x, charge_y, charge_z);
						double distance=1e12*get_distance(ics, ipsilon, zeta, charge_x, charge_y, charge_z);	

						if(distance<20000){
							POTENTIALS_ON_AXIS[INDEX_LAST_STEP][i]+=membrane_charges.at(charge_index).DW_valence*POTENTIAL_C[int(distance)];
						}	
					}		
				}
				
				// induced charges
				if(PRM.EPS_W!=PRM.EPS_MEM){
					for(int surface_index=0; surface_index<NUM_OF_SURFACES; surface_index++){
						double charge_x=SURFACES[surface_index].center_x;
						double charge_y=SURFACES[surface_index].center_y;
						double charge_z=SURFACES[surface_index].center_z;				
						apply_periodic_boundary(ics, ipsilon, zeta, charge_x, charge_y, charge_z);
						double distance=1e12*get_distance(ics, ipsilon, zeta, charge_x, charge_y, charge_z);

						if(distance<20000){
							double surf_charge_valence=(vector_h[surface_index]*surfaces.at(surface_index).area)/Q;
							POTENTIALS_ON_AXIS[INDEX_LAST_STEP][i]+=surf_charge_valence*POTENTIAL_C[int(distance)];
						}	
					}
				}
				
				// ions
				for(int ion_index=0; ion_index<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; ion_index++){
					double ion_x=IONS[INDEX_LAST_STEP][ion_index].x;
					double ion_y=IONS[INDEX_LAST_STEP][ion_index].y;
					double ion_z=IONS[INDEX_LAST_STEP][ion_index].z;			
					apply_periodic_boundary(ics, ipsilon, zeta, ion_x, ion_y, ion_z);
					double distance=1e12*get_distance(ics, ipsilon, zeta, ion_x, ion_y, ion_z);			
					
					if(distance<20000){
						POTENTIALS_ON_AXIS[INDEX_LAST_STEP][i]+=IONS[INDEX_LAST_STEP][ion_index].DW_valence*POTENTIAL_C[int(distance)];
					}					
				}
				
				
				// use microVolts for accumulation in computing the average
				POTENTIALS_ON_AXIS[INDEX_LAST_STEP][i]=1e-6*POTENTIALS_ON_AXIS[INDEX_LAST_STEP][i];
			}
		}
	}
	
	
	
	
	return;
	
}

//#################################################
// returns the energy divided by the valence of the
// source charge and the target charge
// no dielectric - COULOMBIC
//#################################################
double energy_C_at_distance(double distance){	
	
	double temp1=0.5*8.9876*1.602176565*1.602176565;
	double temp2=temp1/distance;
	double energy=1e-17*temp2;
	
	return energy;
}

//#################################################
// returns the potential divided by the valence 
// of the source charge 
// no dielectric - COULOMBIC
//#################################################
double potential_C_at_distance(double distance){
	
	double temp1=8.9876*1.602176565;
	double temp2=temp1/distance;
	double potential=1e2*temp2;
	
	return potential;
}

//#################################################
// returns the electric field divided by the 
// valence of the source charge 
// no dielectric - COULOMBIC
//#################################################
double field_C_at_distance(double distance){
	
	double temp1=8.9876*1.602176565;
	double temp2=distance*distance;
	double temp3=temp1/temp2;
	double field=1e14*temp3;
	
	return field;
}

//#################################################
// returns the force divided by the valence of 
// the source charge and the target charge
// no dielectric - COULOMBIC
//#################################################
double force_C_at_distance(double distance){
	// = e^2 / (4*pi*eps_0*d^2) [N]
	
	double temp1=8.9876*1.602176565*1.602176565;
	double temp2=distance*distance;
	double temp3=temp1/temp2;
	double force=1e-5*temp3;
	
	return force;
}

// -BEGIN
double energy_LJ_at_distance(double distance, Ion i1, Ion i2){
	// V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
	// epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
	// Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
	double rmin_ij = (i1.charmm_half_radius + i2.charmm_half_radius)*1e12; // [pm]
	double eps_ij = sqrt(i1.charmm_eps * i2.charmm_eps); // [J]
	double repulsion_ij = pow(rmin_ij/distance,12.0);
	double london_ij = pow(rmin_ij/distance,6.0);
	return eps_ij * ( repulsion_ij - 2.0*london_ij ); // [J]
}
double force_LJ_at_distance(double distance, Ion i1, Ion i2){
	double rmin_ij = (i1.charmm_half_radius + i2.charmm_half_radius)*1e12; // [pm]
	double eps_ij = sqrt(i1.charmm_eps * i2.charmm_eps); // [J]
	double repulsion_ij = 12.0 * pow(rmin_ij/distance,12.0) / distance; // [pm^-1]
	double london_ij = 6.0 * pow(rmin_ij/distance,6.0) / distance; // [pm^-1]
	return eps_ij * ( repulsion_ij - 2.0*london_ij ) * 1e12; // [N]
}
// -END

//#################################################
// returns the energy divided by the valence of the
// source charge and the target charge
// no dielectric - SHORT-RANGE
//#################################################
double energy_SR_at_distance(double distance, Ion i1, Ion i2){
	
	double energy=0.5e-12*distance*force_SR_at_distance(distance,i1,i2); // [N/C] ;
	
	//~ if(PRM.SR_METHOD.compare("EXPONENTIAL")==0){
		//~ double F0=fabs(i1.valence*i2.valence)*force_C_at_distance(int(i1.pm_radius+i2.pm_radius))/PRM.EPS_W;
		//~ double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
		//~ energy=0.5e-12*F0*K*distance/(PRM.SHORT_RANGE_EXP-double(1.00));
	//~ // -BEGIN
	//~ } else if(PRM.SR_METHOD.compare("LJ")==0){
		//~ if(i1.kind==92 || i2.kind==92){
			//~ double F0=fabs(i1.valence*i2.valence)*force_C_at_distance(int(i1.pm_radius+i2.pm_radius))/PRM.EPS_W;
			//~ double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
			//~ energy=0.5e-12*F0*K*distance/(PRM.SHORT_RANGE_EXP-double(1.00));
		//~ }
		//~ else{
			//~ energy = 0.5 * energy_LJ_at_distance(distance,i1,i2); // [J]
		//~ }
	//~ } else if(PRM.SR_METHOD.compare("PMF")==0){
		//~ // V = PRM(1) * exp((PRM(2) - R)./PRM(3)) * cos(PRM(4)*(PRM(2)-R)*pi) 
		//~ // 		+ PRM(5)*(PRM(2)./R).^9;
		//~ vector<double> SR_PMF_PRM;
		//~ get_parameters_SR_PMF(i1,i2,SR_PMF_PRM);
		//~ if (SR_PMF_PRM.size() == 5) {
			//~ // 	SR_PMF_PRM[0] [J]
			//~ // 	SR_PMF_PRM[1] [pm]
			//~ // 	SR_PMF_PMR[2] [pm]
			//~ // 	SR_PMF_PRM[3] [pm-1]
			//~ // 	SR_PMF_PRM[4] [J]
			//~ energy = SR_PMF_PRM[0] * exp((SR_PMF_PRM[1] - distance)/SR_PMF_PRM[2])
					//~ * cos(SR_PMF_PRM[3]*(SR_PMF_PRM[1]-distance)*M_PI) 
		 		//~ + SR_PMF_PRM[4]*pow((SR_PMF_PRM[1]/distance),9.0);
			//~ energy = 0.5 * energy * PRM.EPS_W / fabs(i1.valence * i2.valence); // [J]
		//~ } else {
			//~ //cout << "WARNING: ion combination " << i1.name << " - " << i2.name << " not allowed for PMF short-range " << endl;
			//~ //cout << "WARNING: \tUsing Lennard-Jones potential instead" << endl;
			//~ energy = 0.5 * energy_LJ_at_distance(distance,i1,i2) // [J]
		//~ }
	//~ }
	//~ // -END

	return energy;
}

//#################################################
// returns the potential divided by the valence 
// of the source charge 
// no dielectric - SHORT-RANGE
//#################################################
double potential_SR_at_distance(double distance, Ion i1, Ion i2){

	double potential=0;
	
	if(PRM.SR_METHOD.compare("EXPONENTIAL")==0){
		double E0=field_C_at_distance(int(i1.pm_radius+i2.pm_radius));
		double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
		potential=1e-12*E0*K*distance/(PRM.SHORT_RANGE_EXP-double(1.00));
	// -BEGIN
	}
	else if(PRM.SR_METHOD.compare("LJ")==0){
		if(i1.kind==92 || i2.kind==92){
			double E0=field_C_at_distance(int(i1.pm_radius+i2.pm_radius));
			double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
			potential=1e-12*E0*K*distance/(PRM.SHORT_RANGE_EXP-double(1.00));
		}
		else{
			potential = 2.0 * energy_SR_at_distance(distance,i1,i2) / i2.charge; // [J/C]
		}
	}
	else if(PRM.SR_METHOD.compare("CS-LJ")==0){		
		
		
	}
	else if(PRM.SR_METHOD.compare("PMF")==0){
		
	}
	
	return potential;
}

//#################################################
// returns the electric field divided by the 
// valence of the source charge 
// no dielectric - SHORT-RANGE
//#################################################
double field_SR_at_distance(double distance, Ion i1, Ion i2){
	
	double field = force_SR_at_distance(distance,i1,i2) / i2.charge; // [N/C] 
	
	//~ if(PRM.SR_METHOD.compare("EXPONENTIAL")==0){
		//~ double E0=field_C_at_distance(int(i1.pm_radius+i2.pm_radius));
		//~ double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
		//~ field=E0*K;
	//~ // -BEGIN
	//~ } else {
		//~ field = force_SR_at_distance(PRM,distance,i1,i2) / i2.charge; // [N/C] 
	//~ }
	//~ // -END
	
	return field;
}

//#################################################
// returns the force divided by the valence of 
// the source charge and the target charge
// no dielectric - SHORT-RANGE
//#################################################
double force_SR_at_distance(double distance, Ion i1, Ion i2){
	
	double force=0;
	
	if(PRM.SR_METHOD.compare("EXPONENTIAL")==0){
		if(i1.valence==0){
			double F0=fabs(i2.valence)*force_C_at_distance(int(i1.pm_radius+i2.pm_radius));
			double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
			force=F0*K/PRM.EPS_W;
		}
		else if(i2.valence==0){
			double F0=fabs(i1.valence)*force_C_at_distance(int(i1.pm_radius+i2.pm_radius));
			double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
			force=F0*K/PRM.EPS_W;
		}
		else if(i1.valence==0 && i2.valence==0){
			double F0=force_C_at_distance(int(i1.pm_radius+i2.pm_radius));
			double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
			force=F0*K/PRM.EPS_W;
		}
		else{			
			double F0=fabs(i1.valence*i2.valence)*force_C_at_distance(int(i1.pm_radius+i2.pm_radius));
			double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
			force=F0*K/PRM.EPS_W;  // [N]
		}
	} 
	else if(PRM.SR_METHOD.compare("LJ")==0){
		if(i1.kind==92 || i2.kind==92){
			double F0=fabs(i1.valence*i2.valence)*force_C_at_distance(int(i1.pm_radius+i2.pm_radius));
			double K=pow((i1.pm_radius+i2.pm_radius)/distance, PRM.SHORT_RANGE_EXP);
			force=F0*K/PRM.EPS_W;	// [N]
		}
		else{
			//~ force = force_LJ_at_distance(distance,i1,i2) * PRM.EPS_W / fabs(i1.valence * i2.valence);
			force = force_LJ_at_distance(distance,i1,i2);
		}
	}
	else if(PRM.SR_METHOD.compare("CS-LJ")==0){		
		if(i1.kind==92){
			if(distance<pow(double(2)/double(5), double(0.17))*1e12*i2.radius-1){
				double k=200;
				double a=double(2.00)*M_PI*PRM.kT/(double(3.00)*1e-12*i2.radius);
				double b=double(18.00)*pow((i2.radius/(1e-12*distance)), 9)/double(15.00);
				double c=double(3.00)*pow((i2.radius/(1e-12*distance)), 3);
				
				force=k*a*(b-c);
			}
			else{
				force=0;
			}
		}
		else if(i2.kind==92){
			if(distance<pow(double(2)/double(5), double(0.17))*1e12*i1.radius-1){
				double k=200;
				double a=double(2.00)*M_PI*PRM.kT/(double(3.00)*1e-12*distance);
				double b=double(18.00)*pow((i1.radius/(1e-12*distance)), 9)/double(15.00);
				double c=double(3.00)*pow((i1.radius/(1e-12*distance)), 3);
				
				force=k*a*(b-c);
			}
			else{
				force=0;
			}
		}
		else{
			double sigma=i1.radius+i2.radius;
			double k=0.1;
			if(distance<pow(double(2), double(0.17))*1e12*sigma-1){
				double a=k*double(24.00)*PRM.kT/(1e-12*distance);
				double b=2*pow(sigma/(1e-12*distance), 12);
				double c=pow(sigma/(1e-12*distance), 6);
				
				force=a*(b-c);
			}
			else{
				force=0;
			}
		}
	}

	else if(PRM.SR_METHOD.compare("PMF")==0){
		vector<double> SR_PMF_PRM;
		get_parameters_SR_PMF(i1,i2,SR_PMF_PRM);
		if (SR_PMF_PRM.size() == 5) {
			// 	SR_PMF_PRM[0] [J]
			// 	SR_PMF_PRM[1] [pm]
			// 	SR_PMF_PMR[2] [pm]
			// 	SR_PMF_PRM[3] [pm-1]
			// 	SR_PMF_PRM[4] [J]
			force = SR_PMF_PRM[0] * exp((SR_PMF_PRM[1] - distance)/SR_PMF_PRM[2]) 
				* ( +(1.0/SR_PMF_PRM[2])
					* cos(SR_PMF_PRM[3]*(SR_PMF_PRM[1]-distance)*M_PI) 
				     -(M_PI*SR_PMF_PRM[3])
					* sin(SR_PMF_PRM[3]*(SR_PMF_PRM[1]-distance)*M_PI) ) 
				+ 9*(SR_PMF_PRM[4]/SR_PMF_PRM[1])*pow((SR_PMF_PRM[1]/distance),10.0);
			force = 1e12 * force * PRM.EPS_W / fabs(i1.valence * i2.valence); // [N]
		} else {
			//cout << "WARNING: ion combination " << i1.name << " - " << i2.name << " not allowed for PMF short-range " << endl;
			//cout << "WARNING: \tUsing Lennard-Jones potential instead" << endl;
			force = force_LJ_at_distance(distance,i1,i2) * PRM.EPS_W / fabs(i1.valence * i2.valence);
		}
	}
	//cout<<"DEBUG "<<i1.name<<"\t"<<i2.name<<"\t"<<distance<<"\t"<<force<<endl;
	// -END
	
	return force;
}


// EMMA DISEGUALE 
double force_other_at_distance(double distance, Ion i1, Ion i2){
	
	double force=0;
	
	return force;
}

// EMMA DISEGUALE 
double energy_other_at_distance(double distance, Ion i1, Ion i2){
	
	double energy=0.5e-12*distance*force_other_at_distance(distance,i1,i2); // [N/C] ;
		
	return energy;
}


void create_EPFFs(){
	for(int i=0; i<20000; i++){
		ENERGY_C[i]=energy_C_at_distance(double(i));
		POTENTIAL_C[i]=potential_C_at_distance(double(i));
		FIELD_C[i]=field_C_at_distance(double(i));
		FORCE_C[i]=force_C_at_distance(double(i));
	}
/*===========================================	
  kind  name
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
	
// J, Q and X are for positive, negative and neutral ions of the ion boxes respectively
	20		J00		0.25e	0.25A	
	21		J01		0.25e	0.50A	
	22		J02		0.25e	0.75A
	23		J03		0.25e	1.00A	
	24		J04		0.25e	1.25A	
	25		J05		0.25e	1.50A	
	
	26		J06		0.50e	0.25A	
	27		J07		0.50e	0.50A	
	28		J08		0.50e	0.75A	
	29		J09		0.50e	1.00A	
	30		J10		0.50e	1.25A	
	31		J11		0.50e	1.50A	
	
	32		J12		0.75e	0.25A
	33		J13		0.75e	0.50A	
	34		J14		0.75e	0.75A	
	35		J15		0.75e	1.00A		
	36		J16		0.75e	1.25A		
	37		J17		0.75e	1.50A		

	38		J18		1.00e	0.25A
	39		J19		1.00e	0.50A	
	40		J20		1.00e	0.75A	
	41		J21		1.00e	1.00A		
	42		J22		1.00e	1.25A	
	43		J23		1.00e	1.50A	

	44		J24		0.50e	1.40A	//custom 1
	45		J25		1.00e	1.20A	//custom 2
	46		J26		1.00e	0.80A	//custom 3
	47		J27		1.00e	1.15A	//custom 4
	48		J28		1.00e	1.35A	//custom 5
	49		J29		1.00e	1.88A	//custom 6


	50		Q00		-0.25e	0.25A	
	51		Q01		-0.25e	0.50A	
	52		Q02		-0.25e	0.75A
	53		Q03		-0.25e	1.00A	
	54		Q04		-0.25e	1.25A	
	55		Q05		-0.25e	1.50A	
	
	56		Q06		-0.50e	0.25A	
	57		Q07		-0.50e	0.50A	
	58		Q08		-0.50e	0.75A	
	59		Q09		-0.50e	1.00A	
	60		Q10		-0.50e	1.25A	
	61		Q11		-0.50e	1.50A	
	
	62		Q12		-0.75e	0.25A
	63		Q13		-0.75e	0.50A	
	64		Q14		-0.75e	0.75A	
	65		Q15		-0.75e	1.00A		
	66		Q16		-0.75e	1.25A		
	67		Q17		-0.75e	1.50A		

	68		Q18		-1.00e	0.25A
	69		Q19		-1.00e	0.50A	
	70		Q20		-1.00e	0.75A	
	71		Q21		-1.00e	1.00A		
	72		Q22		-1.00e	1.25A	
	73		Q23		-1.00e	1.50A	

	74		Q24		-0.50e	1.40A	//custom 1
	75		Q25		-1.00e	1.20A	//custom 2
	76		Q26		-1.00e	0.80A	//custom 3
	77		Q27		-1.00e	1.15A	//custom 4
	78		Q28		-1.00e	1.35A	//custom 5
	79		Q29		-1.00e	1.88A	//custom 6
	
	
	80		X00		0.00e	0.25A	
	81		X01		0.00e	0.50A	
	82		X02		0.00e	0.75A
	83		X03		0.00e	1.00A	
	84		X04		0.00e	1.25A	
	85		X05		0.00e	1.50A	
	
	86		X06		0.00e	1.40A	//custom 1
	87		X07		0.00e	1.20A	//custom 2
	88		X08		0.00e	0.80A	//custom 3
	89		X09		0.00e	1.15A	//custom 4
	90		X10		0.00e	1.35A	//custom 5
	91		X11		0.00e	1.88A	//custom 6
	
	92		membrane wall
														*/
//===========================================	

//==================================== short-range - begin
	vector < vector < double > > aux_vec1;
	vector < double > aux_vec2;
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		Ion ion_i(i);
		aux_vec1.clear();
		if(PRM.ions_to_simulate.at(i)){
			for(int j=0; j<NUM_OF_IONIC_SPECIES; j++){
				Ion ion_j(j);
				aux_vec2.clear();
				if(PRM.ions_to_simulate.at(j)){
					
					for(int k=0; k<4000; k++){
						aux_vec2.push_back(force_SR_at_distance(double(k), ion_i, ion_j));
					}
				}
				aux_vec1.push_back(aux_vec2);
				aux_vec2.clear();
			}
			Ion ion_membrane(92);
			aux_vec2.clear();
			for(int k=0; k<4000; k++){
				aux_vec2.push_back(force_SR_at_distance(double(k), ion_i, ion_membrane));
			}
			aux_vec1.push_back(aux_vec2);
			aux_vec2.clear();
		}
		
		FORCE_SR.push_back(aux_vec1);
		aux_vec1.clear();
	}
	
	aux_vec1.clear();
	aux_vec2.clear();
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		Ion ion_i(i);
		aux_vec1.clear();
		if(PRM.ions_to_simulate.at(i)){
			for(int j=0; j<NUM_OF_IONIC_SPECIES; j++){
				Ion ion_j(j);
				aux_vec2.clear();
				if(PRM.ions_to_simulate.at(j)){
					for(int k=0; k<4000; k++){
						aux_vec2.push_back(energy_SR_at_distance(double(k), ion_i, ion_j));
					}
				}
				aux_vec1.push_back(aux_vec2);
				aux_vec2.clear();
			}
			Ion ion_membrane(92);
			aux_vec2.clear();
			for(int k=0; k<4000; k++){
				aux_vec2.push_back(energy_SR_at_distance(double(k), ion_i, ion_membrane));
			}
			aux_vec1.push_back(aux_vec2);
			aux_vec2.clear();
		}
		
		ENERGY_SR.push_back(aux_vec1);
		aux_vec1.clear();
	}
//==================================== short-range - end

//==================================== other - begin
	aux_vec1.clear();
	aux_vec2.clear();
	
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		Ion ion_i(i);
		aux_vec1.clear();
		if(PRM.ions_to_simulate.at(i)){
			for(int j=0; j<NUM_OF_IONIC_SPECIES; j++){
				Ion ion_j(j);
				aux_vec2.clear();
				if(PRM.ions_to_simulate.at(j)){
					for(int k=0; k<4000; k++){
						// EMMA DISEGUALE 
						aux_vec2.push_back(force_other_at_distance(double(k), ion_i, ion_j));
					}
				}
				aux_vec1.push_back(aux_vec2);
				aux_vec2.clear();
			}
			Ion ion_membrane(92);
			aux_vec2.clear();
			for(int k=0; k<4000; k++){
				// EMMA DISEGUALE 
				aux_vec2.push_back(force_other_at_distance(double(k), ion_i, ion_membrane));
			}
			aux_vec1.push_back(aux_vec2);
			aux_vec2.clear();
		}
		
		FORCE_OTHER.push_back(aux_vec1);
		aux_vec1.clear();
	}
	
	aux_vec1.clear();
	aux_vec2.clear();
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		Ion ion_i(i);
		aux_vec1.clear();
		if(PRM.ions_to_simulate.at(i)){
			for(int j=0; j<NUM_OF_IONIC_SPECIES; j++){
				Ion ion_j(j);
				aux_vec2.clear();
				if(PRM.ions_to_simulate.at(j)){
					for(int k=0; k<4000; k++){
						// EMMA DISEGUALE 
						aux_vec2.push_back(energy_other_at_distance(double(k), ion_i, ion_j));
					}
				}
				aux_vec1.push_back(aux_vec2);
				aux_vec2.clear();
			}
			Ion ion_membrane(92);
			aux_vec2.clear();
			for(int k=0; k<4000; k++){
				// EMMA DISEGUALE 
				aux_vec2.push_back(energy_other_at_distance(double(k), ion_i, ion_membrane));
			}
			aux_vec1.push_back(aux_vec2);
			aux_vec2.clear();
		}
		
		ENERGY_OTHER.push_back(aux_vec1);
		aux_vec1.clear();
	}
//==================================== other - end	
	return;
}

void create_RDF_vector(vector < vector < vector < double > > >& RDF, vector < double >& RDF_samples){

	
//===========================================	
														/*
  kind  name
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
	
// J, Q and X are for positive, negative and neutral ions of the ion boxes respectively
	20		J00		0.25e	0.25A	
	21		J01		0.25e	0.50A	
	22		J02		0.25e	0.75A
	23		J03		0.25e	1.00A	
	24		J04		0.25e	1.25A	
	25		J05		0.25e	1.50A	
	
	26		J06		0.50e	0.25A	
	27		J07		0.50e	0.50A	
	28		J08		0.50e	0.75A	
	29		J09		0.50e	1.00A	
	30		J10		0.50e	1.25A	
	31		J11		0.50e	1.50A	
	
	32		J12		0.75e	0.25A
	33		J13		0.75e	0.50A	
	34		J14		0.75e	0.75A	
	35		J15		0.75e	1.00A		
	36		J16		0.75e	1.25A		
	37		J17		0.75e	1.50A		

	38		J18		1.00e	0.25A
	39		J19		1.00e	0.50A	
	40		J20		1.00e	0.75A	
	41		J21		1.00e	1.00A		
	42		J22		1.00e	1.25A	
	43		J23		1.00e	1.50A	

	44		J24		0.50e	1.40A	//custom 1
	45		J25		1.00e	1.20A	//custom 2
	46		J26		1.00e	0.80A	//custom 3
	47		J27		1.00e	1.15A	//custom 4
	48		J28		1.00e	1.35A	//custom 5
	49		J29		1.00e	1.88A	//custom 6


	50		Q00		-0.25e	0.25A	
	51		Q01		-0.25e	0.50A	
	52		Q02		-0.25e	0.75A
	53		Q03		-0.25e	1.00A	
	54		Q04		-0.25e	1.25A	
	55		Q05		-0.25e	1.50A	
	
	56		Q06		-0.50e	0.25A	
	57		Q07		-0.50e	0.50A	
	58		Q08		-0.50e	0.75A	
	59		Q09		-0.50e	1.00A	
	60		Q10		-0.50e	1.25A	
	61		Q11		-0.50e	1.50A	
	
	62		Q12		-0.75e	0.25A
	63		Q13		-0.75e	0.50A	
	64		Q14		-0.75e	0.75A	
	65		Q15		-0.75e	1.00A		
	66		Q16		-0.75e	1.25A		
	67		Q17		-0.75e	1.50A		

	68		Q18		-1.00e	0.25A
	69		Q19		-1.00e	0.50A	
	70		Q20		-1.00e	0.75A	
	71		Q21		-1.00e	1.00A		
	72		Q22		-1.00e	1.25A	
	73		Q23		-1.00e	1.50A	

	74		Q24		-0.50e	1.40A	//custom 1
	75		Q25		-1.00e	1.20A	//custom 2
	76		Q26		-1.00e	0.80A	//custom 3
	77		Q27		-1.00e	1.15A	//custom 4
	78		Q28		-1.00e	1.35A	//custom 5
	79		Q29		-1.00e	1.88A	//custom 6
	
	
	80		X00		0.00e	0.25A	
	81		X01		0.00e	0.50A	
	82		X02		0.00e	0.75A
	83		X03		0.00e	1.00A	
	84		X04		0.00e	1.25A	
	85		X05		0.00e	1.50A	
	
	86		X06		0.00e	1.40A	//custom 1
	87		X07		0.00e	1.20A	//custom 2
	88		X08		0.00e	0.80A	//custom 3
	89		X09		0.00e	1.15A	//custom 4
	90		X10		0.00e	1.35A	//custom 5
	91		X11		0.00e	1.88A	//custom 6
	
														*/
//===========================================	


	vector < vector < double > > aux_vec1;
	vector < double > aux_vec2;
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		Ion ion_i(i);
		aux_vec1.clear();
		if(PRM.ions_to_simulate.at(i)){
			for(int j=0; j<NUM_OF_IONIC_SPECIES; j++){
				Ion ion_j(j);
				aux_vec2.clear();
				if(PRM.ions_to_simulate.at(j)){
					for(int k=0; k<1000; k++){
						aux_vec2.push_back(0.00);
					}
				}
				aux_vec1.push_back(aux_vec2);
				aux_vec2.clear();
			}
			aux_vec2.clear();
		}
		
		RDF.push_back(aux_vec1);
		aux_vec1.clear();
	}
	
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		RDF_samples.push_back(0.00);
	}
	
	return;
}

void update_RDF_vector(vector < vector < vector < double > > >& RDF, vector < double >& RDF_samples){
	
	for(int ion_index_1=0; ion_index_1<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; ion_index_1++){
		RDF_samples.at(IONS[INDEX_STAT_STEP][ion_index_1].kind)+=double(1.00);
		for(int ion_index_2=0; ion_index_2<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; ion_index_2++){
			if(ion_index_1!=ion_index_2){
				double ion2_x=IONS[INDEX_STAT_STEP][ion_index_2].x_prev;
				double ion2_y=IONS[INDEX_STAT_STEP][ion_index_2].y_prev;
				double ion2_z=IONS[INDEX_STAT_STEP][ion_index_2].z_prev;			
				apply_periodic_boundary(IONS[INDEX_STAT_STEP][ion_index_1].x_prev, IONS[INDEX_STAT_STEP][ion_index_1].y_prev, IONS[INDEX_STAT_STEP][ion_index_1].z_prev, ion2_x, ion2_y, ion2_z);
				double distance=1e12*get_distance(IONS[INDEX_STAT_STEP][ion_index_1].x_prev, IONS[INDEX_STAT_STEP][ion_index_1].y_prev, IONS[INDEX_STAT_STEP][ion_index_1].z_prev, ion2_x, ion2_y, ion2_z);	
				int index=int(distance);
				if(index<1000){
					RDF.at(IONS[INDEX_STAT_STEP][ion_index_1].kind).at(IONS[INDEX_STAT_STEP][ion_index_2].kind).at(index)+=double(1.00);
				}
			}
		}
	}

	//~ for(int ion_index_1=0; ion_index_1<ions.front().size(); ion_index_1++){
		//~ RDF_samples.at(ions.front().at(ion_index_1).kind)+=double(1.00);
		//~ for(int ion_index_2=0; ion_index_2<ions.front().size(); ion_index_2++){
			//~ if(ion_index_1!=ion_index_2){
				//~ double ion2_x=ions.front().at(ion_index_2).x_prev;
				//~ double ion2_y=ions.front().at(ion_index_2).y_prev;
				//~ double ion2_z=ions.front().at(ion_index_2).z_prev;			
				//~ apply_periodic_boundary(ions.front().at(ion_index_1).x_prev, ions.front().at(ion_index_1).y_prev, ions.front().at(ion_index_1).z_prev, ion2_x, ion2_y, ion2_z);
				//~ double distance=1e12*get_distance(ions.front().at(ion_index_1).x_prev, ions.front().at(ion_index_1).y_prev, ions.front().at(ion_index_1).z_prev, ion2_x, ion2_y, ion2_z);	
				//~ int index=int(distance);
				//~ if(index<1000){
					//~ RDF.at(ions.front().at(ion_index_1).kind).at(ions.front().at(ion_index_2).kind).at(index)+=double(1.00);
				//~ }
			//~ }
		//~ }
	//~ }
	
	return;
}

void compute_RDF_vector(vector < vector < vector < double > > >& RDF){
	
	for(int ion_index_1=0; ion_index_1<NUM_OF_IONIC_SPECIES; ion_index_1++){
		if(!RDF.at(ion_index_1).empty()){
			for(int ion_index_2=0; ion_index_2<NUM_OF_IONIC_SPECIES; ion_index_2++){
				if(!RDF.at(ion_index_2).empty()){
					for(int k=0; k<1000; k++){
						RDF.at(ion_index_1).at(ion_index_2).at(k)=RDF.at(ion_index_1).at(ion_index_2).at(k)/double(PRM.SIM_STEPS);
					}
				}
			}
		}
	}
	
	return;
}

void compute_SR_dielectric_boundary(){
	
	cout << "computing SR_dielectric_boundary........"<<flush;

	bool to_compute=true;
	
	if(to_compute){

		PRM.Z_BEGIN_SR_MAP=-PRM.MEMBRANE_WIDTH/double(2.00)-double(400.00);
		PRM.Z_END_SR_MAP=PRM.MEMBRANE_WIDTH/double(2.00)+double(400.00);
		double MAX_RADIUS_SR_MAP=PRM.MAX_X*double(1.41421356);
		
		vector < vector < vector < double > > >  aux_vec_1;
		vector < vector < double > > aux_vec_2;
		vector < double > aux_vec_3;

		//~ cout << "PRM.Z_BEGIN_SR_MAP: " << PRM.Z_BEGIN_SR_MAP << endl;
		//~ cout << "PRM.Z_END_SR_MAP: " << PRM.Z_END_SR_MAP << endl;
		//~ cout << "MAX_RADIUS_SR_MAP: " << MAX_RADIUS_SR_MAP << endl;
		//~ cout << "PRM.MEMBRANE_WIDTH: " << PRM.MEMBRANE_WIDTH << endl;
		//~ cout << "limits.size(): " << limits.size() << endl;
		
		RadialPoint rp;

		vector <RadialPoint> profile;
		RadialPoint pp;
		for(int i=ceil(MAX_RADIUS_SR_MAP); i>round(limits.at(0)); i--){
			pp.z=-PRM.MEMBRANE_WIDTH/double(2.00);
			pp.r=double(i);
			profile.push_back(pp);
		}
		for(int i=0; i<limits.size(); i++){
			pp.z=-PRM.MEMBRANE_WIDTH/double(2.00)+double(i);
			pp.r=limits.at(i);
			profile.push_back(pp);
		}
		for(int i=round(limits.at(limits.size()-1)); i<=ceil(MAX_RADIUS_SR_MAP); i++){
			pp.z=PRM.MEMBRANE_WIDTH/double(2.00);
			pp.r=double(i);
			profile.push_back(pp);
		}
	
		for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
			Ion ion_i(i);
			aux_vec_1.clear();
			if(PRM.ions_to_simulate.at(i)){
				
				
				for(int z=PRM.Z_BEGIN_SR_MAP; z<=PRM.Z_END_SR_MAP; z=z+PRM.DELTA_SR_MAP){
					
					aux_vec_2.clear();
					for(int r=0; r<=MAX_RADIUS_SR_MAP; r=r+PRM.DELTA_SR_MAP){
						
						aux_vec_3.clear();
						
						rp.r=double(r);
						rp.z=double(z);
						
						double min_dist=1e40;
						int ind_on_lim=-1;
						
						for(int kk=0; kk<profile.size(); kk++){
							double distance=get_distance(rp, profile.at(kk));	
							
							if(distance<min_dist){
								min_dist=distance;
								ind_on_lim=kk;
							}
						}
						
						double force_SR=0.00;
						double dz=0.00;
						double dr=0.00;
						
						ion_i.force[0]=0.00;
						ion_i.force[1]=0.00;
						ion_i.force[2]=0.00;
						

						if(min_dist<4000){
							dz=(rp.z-profile.at(ind_on_lim).z)/min_dist;
							dr=(rp.r-profile.at(ind_on_lim).r)/min_dist;
							
							force_SR=FORCE_SR.at(ion_i.kind).at(NUM_OF_IONIC_SPECIES).at(int(min_dist))*PRM.EPS_W;
							
							
							ion_i.force[0]=force_SR*dr;
							ion_i.force[1]=0.00;
							ion_i.force[2]=force_SR*dz;
						}
						else{
						
						}
							
						aux_vec_3.push_back(ion_i.force[0]);
						aux_vec_3.push_back(ion_i.force[1]);
						aux_vec_3.push_back(ion_i.force[2]);
						
						aux_vec_2.push_back(aux_vec_3);
					}
					
					aux_vec_1.push_back(aux_vec_2);	
				}
				
			}
			SR_dielectric_boundary.push_back(aux_vec_1);	
		}
	}
	
	cout << "OK!"<<endl;

	return;
}

void create_SR_ion_boxes(){
	
	if(!ion_boxes.empty()){
		
		Ion ion_membrane(92);
	
		vector <double> aux_vec1;
		aux_vec1.clear();

		for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
			Ion ion_i(i);
			aux_vec1.clear();
			if(PRM.ions_to_simulate.at(i)){

				for(int z=-PRM.MEMBRANE_WIDTH/double(2.00); z<=PRM.MEMBRANE_WIDTH/double(2.00); z++){
					double force=0.00;
					
					for(int ib=0; ib<ion_boxes.size(); ib++){
						if(ion_boxes.at(ib).index!=0){
							for(int iv=0; iv<ion_boxes.at(ib).is.size(); iv++){
								if(ion_boxes.at(ib).is.at(iv).compare(ion_i.name)==0 && (ion_i.name.substr(0, 1).compare("J")==0 || ion_i.name.substr(0, 1).compare("Q")==0 || ion_i.name.substr(0, 1).compare("X")==0)){
									double dist1=double(z)-ion_boxes.at(ib).MIN_Z;
									if(dist1>0 && dist1<4000 && z<ion_boxes.at(ib).MAX_Z){
										force=force+force_SR_at_distance(dist1, ion_i, ion_membrane);
									}
								
									double dist2=ion_boxes.at(ib).MAX_Z-double(z);
									if(dist2>0 && dist2<4000 && z>ion_boxes.at(ib).MIN_Z){
										force=force-force_SR_at_distance(dist2, ion_i, ion_membrane);
									}
								}
							}
						}
					}
					
					aux_vec1.push_back(force);
				}
				
				aux_vec1.at(0)=aux_vec1.at(1);
				aux_vec1.at(aux_vec1.size()-1)=aux_vec1.at(aux_vec1.size()-2);
				
			}
			
			SR_ion_boxes.push_back(aux_vec1);
		}
	}		
	
	return;
}

void compute_potential_on_trajectories(){
	
	double *vector_c_T1, *vector_c_T2, *vector_c_T3, *vector_h_T1, *vector_h_T2, *vector_h_T3;
	
	posix_memalign((void **)(&vector_c_T1), 128, sizeof(double)*surfaces.size());
	posix_memalign((void **)(&vector_h_T1), 128, sizeof(double)*surfaces.size());
	posix_memalign((void **)(&vector_c_T2), 128, sizeof(double)*surfaces.size());
	posix_memalign((void **)(&vector_h_T2), 128, sizeof(double)*surfaces.size());
	posix_memalign((void **)(&vector_c_T3), 128, sizeof(double)*surfaces.size());
	posix_memalign((void **)(&vector_h_T3), 128, sizeof(double)*surfaces.size());
	
	if((!vector_c_T1) || (!vector_h_T1) || (!vector_c_T2) || (!vector_h_T2) || (!vector_c_T3) || (!vector_h_T3)){
		cerr<<"ERROR: Allocation of vector_c/h _T1/_T2/_T3 failed"<<endl;
		exit(1);
	}
	
	
	//trajectories parallel to channel axis
	for(int i=0; i<2000; i=i+100){				
		double dfa=double(0.01)*double(i);

		string s;
		stringstream ss;
		ss << dfa;
		ss >> s;
		
		string traj_file= PRM.PREFIX + "_A_" + s + ".dat";
		
		cout << "writing " <<traj_file<<"..........";
		
		createFile(traj_file);
		char *tfn1 = new char[traj_file.length()+1];
		strcpy(tfn1, traj_file.c_str());     
		ofstream fout1(tfn1, ios::app);
		
		fout1<<"#Trajectory data file"<<endl;
		fout1<<"#column 1: z (A)"<<endl;
		fout1<<"#column 2: potential from protein charges, external field and induced charges (V)"<<endl;
		fout1<<"#column 3: reaction potential for the ion (V)"<<endl;
		fout1<<"#column 4: potential from protein charges, external field and induced charges and reaction potential (V)"<<endl;
		fout1<<"#column 5: electric field along z for col. 2"<<endl;
		fout1<<"#column 6: electric field along z for col. 3"<<endl;
		fout1<<"#column 7: electric field along z for col. 4"<<endl;
		
		Ion charge(4);
		charge.x=1e-12*double(i);
		charge.y=0.00;
		charge.z=0.00;
		
		charge.charge=double(1.00)*Q;
		charge.valence=1.00;
		charge.EPS_ION=PRM.EPS_W;
		charge.DW_charge=charge.charge/charge.EPS_ION;
		charge.DW_valence=charge.valence/charge.EPS_ION;
		
		//~ string traj_file2= PRM.PREFIX + "_A_" + s + ".dat4";
		//~ createFile(traj_file2);
		//~ char *tfn2 = new char[traj_file2.length()+1];
		//~ strcpy(tfn2, traj_file2.c_str());     
		//~ ofstream fout2(tfn2, ios::app);
		
		for(int z=PRM.MIN_Z; z<=PRM.MAX_Z; z=z+10){
			
			
			
			charge.z=1e-12*double(z);
			
		// compute induced charges		
			if(NUM_OF_SURFACES>0){
				
				if(PRM.EPS_W!=PRM.EPS_MEM){
					
					
					
					for(int index_on_h=0; index_on_h<NUM_OF_SURFACES; index_on_h++){
						vector_c_T1[index_on_h]=_double_zero;
						vector_h_T1[index_on_h]=_double_zero;
						vector_c_T2[index_on_h]=_double_zero;
						vector_h_T2[index_on_h]=_double_zero;
						vector_c_T3[index_on_h]=_double_zero;
						vector_h_T3[index_on_h]=_double_zero;
					}
		
					double elem_c_T1=0.00;
					double elem_c_T2=0.00;
					double elem_c_T3=0.00;
					
					for(int surface_index=0; surface_index<NUM_OF_SURFACES; surface_index++){

				
				// external field and membrane charges		
						elem_c_T1=C_BASE[surface_index];
						elem_c_T2=0;
						elem_c_T3=C_BASE[surface_index];
						
						//~ vector_c_T1[surface_index]=c_base.at(surface_index);
						//~ vector_c_T3[surface_index]=c_base.at(surface_index);
						
						
						//source charge
						double distance=1e12*get_distance(SURFACES[surface_index].center_x, SURFACES[surface_index].center_y, SURFACES[surface_index].center_z, charge.x, charge.y, charge.z);		
						
						if(distance<20000){
							
							double temp1=PRM.dielectrics_weight*charge.DW_valence;	
							double field=temp1*FIELD_C[int(distance)]/COULOMB_K;
							
							double cs_x=(SURFACES[surface_index].center_x-charge.x)/(1e-12*distance);
							double cs_y=(SURFACES[surface_index].center_y-charge.y)/(1e-12*distance);
							double cs_z=(SURFACES[surface_index].center_z-charge.z)/(1e-12*distance);
							
							double scalar=cs_x*SURFACES[surface_index].normal[0]+cs_y*SURFACES[surface_index].normal[1]+cs_z*SURFACES[surface_index].normal[2];
							
								//~ vector_c_T2[surface_index]+=1e24*SURFACES[surface_index].area*field*scalar;
								//~ vector_c_T3[surface_index]+=1e24*SURFACES[surface_index].area*field*scalar;
								
								elem_c_T2+=1e24*SURFACES[surface_index].area*field*scalar;
								elem_c_T3+=1e24*SURFACES[surface_index].area*field*scalar;
						}
						else{
							//cut-off
						}
						
						vector_c_T1[surface_index]=elem_c_T1;
						vector_c_T2[surface_index]=elem_c_T2;
						vector_c_T3[surface_index]=elem_c_T3;
					}
					
					
					
		// solve the ICC equation 		
					dgemv_(&_N, &matrix_size, &matrix_size, &_double_one, inverse_a, &matrix_size, vector_c_T1, &_int_one, &_double_zero, vector_h_T1, &_int_one);	
					dgemv_(&_N, &matrix_size, &matrix_size, &_double_one, inverse_a, &matrix_size, vector_c_T2, &_int_one, &_double_zero, vector_h_T2, &_int_one);	
					dgemv_(&_N, &matrix_size, &matrix_size, &_double_one, inverse_a, &matrix_size, vector_c_T3, &_int_one, &_double_zero, vector_h_T3, &_int_one);	
								
//===========================================================================================	
//	SUM RULE CHECK....
	
					//~ double K=-(PRM.EPS_MEM*PRM.EPS_W)/(PRM.EPS_W-PRM.EPS_MEM);					
					//~ double total_charge_T1=0;
					//~ double total_charge_T2=0;
					//~ double total_charge_T3=0;
					//~ for(int surface_index=0; surface_index<surfaces.size(); surface_index++){		
						
						//~ double surf_charge_T1=(vector_h_T1[surface_index]*SURFACES[surface_index].area)/Q;
						//~ double surf_charge_T2=(vector_h_T2[surface_index]*SURFACES[surface_index].area)/Q;
						//~ double surf_charge_T3=(vector_h_T3[surface_index]*SURFACES[surface_index].area)/Q;

						//~ total_charge_T1+=surf_charge_T1;
						//~ total_charge_T2+=surf_charge_T2;
						//~ total_charge_T3+=surf_charge_T3;
					//~ }
					
					//~ total_charge_T1=K*total_charge_T1;
					//~ total_charge_T2=K*total_charge_T2;
					//~ total_charge_T3=K*total_charge_T3;
					
					//~ fout2<< 1e10*charge.z << " " << total_charge_T1 << " " << total_charge_T2 << " " << total_charge_T3<<endl;

//===========================================================================================					
										
					
					
				}
			}
			
			
			
			double tempAA=0.00;
			double tempBB=0.00;	
			tempAA=(PRM.MAX_Z-1e12*charge.z)/(PRM.SIM_DOMAIN_WIDTH_Z);
			tempBB=PRM.APPLIED_POTENTIAL;	
			
			double pot_T1=tempAA*tempBB;
			double pot_T2=0.00;
			double pot_T3=tempAA*tempBB;
			
			double energy_T1=pot_T1*charge.charge;
			double energy_T2=pot_T2*charge.charge;
			double energy_T3=pot_T3*charge.charge;
			
			vector <double> force_T1;
			vector <double> force_T2;
			vector <double> force_T3;
			
			for(int iff=0; iff<3; iff++){
				force_T1.push_back(0.00);
				force_T2.push_back(0.00);
				force_T3.push_back(0.00);
			}
			
			// external field
			force_T1.at(2)=0.00;
			force_T2.at(2)=0.00;
			force_T3.at(2)=0.00;
			
			force_T1.at(2)=Q*charge.valence*PRM.APPLIED_FIELD;
			force_T2.at(2)=0.00;
			force_T3.at(2)=Q*charge.valence*PRM.APPLIED_FIELD;
		
				
				
			
			
			//membrane charges
			if(!membrane_charges.empty()){
				for(int charge_index=0; charge_index<membrane_charges.size(); charge_index++){
					double charge_x=membrane_charges.at(charge_index).x;
					double charge_y=membrane_charges.at(charge_index).y;
					double charge_z=membrane_charges.at(charge_index).z;			
					
					double distance=1e12*get_distance(charge.x, charge.y, charge.z, charge_x, charge_y, charge_z);			
					
					if(distance<20000){
						double force_C=0.00;
						double force_SR=0.00;
						double dx=0.00;
						double dy=0.00;
						double dz=0.00;
						
						dx=1e12*(charge.x-charge_x)/distance;
						dy=1e12*(charge.y-charge_y)/distance;
						dz=1e12*(charge.z-charge_z)/distance;
						
						force_C=charge.valence*membrane_charges.at(charge_index).DW_valence*FORCE_C[int(distance)];
						
						force_T1.at(0)+=(force_C)*dx;
						force_T1.at(1)+=(force_C)*dy;
						force_T1.at(2)+=(force_C)*dz;
						force_T3.at(0)+=(force_C)*dx;
						force_T3.at(1)+=(force_C)*dy;
						force_T3.at(2)+=(force_C)*dz;
						
						energy_T1+=charge.valence*membrane_charges.at(charge_index).DW_valence*ENERGY_C[int(distance)];
						energy_T3+=charge.valence*membrane_charges.at(charge_index).DW_valence*ENERGY_C[int(distance)];
						
						
						pot_T1+=membrane_charges.at(charge_index).DW_valence*POTENTIAL_C[int(distance)];
						pot_T3+=membrane_charges.at(charge_index).DW_valence*POTENTIAL_C[int(distance)];
					}	
					else{
						// cut-off
					}		
				}		
			}
			
			
		
			// boundary repulsion and induced charges
			if(NUM_OF_SURFACES>0){
				if(PRM.EPS_W!=PRM.EPS_MEM){
					for(int surface_index=0; surface_index<NUM_OF_SURFACES; surface_index++){		
						double charge_x=SURFACES[surface_index].center_x;
						double charge_y=SURFACES[surface_index].center_y;
						double charge_z=SURFACES[surface_index].center_z;				
						double distance=1e12*get_distance(charge.x, charge.y, charge.z, charge_x, charge_y, charge_z);	

						if(distance<20000){
							double force_C=0.00;
							double force_SR=0.00;
							double dx=0.00;
							double dy=0.00;
							double dz=0.00;
							
							dx=1e12*(charge.x-charge_x)/distance;
							dy=1e12*(charge.y-charge_y)/distance;
							dz=1e12*(charge.z-charge_z)/distance;
							
							double surf_charge_valence_T1=(vector_h_T1[surface_index]*SURFACES[surface_index].area)/Q;
							double force_C_T1=charge.valence*surf_charge_valence_T1*FORCE_C[int(distance)];
							double surf_charge_valence_T2=(vector_h_T2[surface_index]*SURFACES[surface_index].area)/Q;
							double force_C_T2=charge.valence*surf_charge_valence_T2*FORCE_C[int(distance)];
							double surf_charge_valence_T3=(vector_h_T3[surface_index]*SURFACES[surface_index].area)/Q;
							double force_C_T3=charge.valence*surf_charge_valence_T3*FORCE_C[int(distance)];
							
							
							force_T1.at(0)+=(force_C_T1)*dx;
							force_T1.at(1)+=(force_C_T1)*dy;
							force_T1.at(2)+=(force_C_T1)*dz;
							
							force_T2.at(0)+=(force_C_T2)*dx;
							force_T2.at(1)+=(force_C_T2)*dy;
							force_T2.at(2)+=(force_C_T2)*dz;
							
							force_T3.at(0)+=(force_C_T3)*dx;
							force_T3.at(1)+=(force_C_T3)*dy;
							force_T3.at(2)+=(force_C_T3)*dz;
							
							energy_T1+=charge.valence*surf_charge_valence_T1*ENERGY_C[int(distance)];
							energy_T2+=charge.valence*surf_charge_valence_T2*ENERGY_C[int(distance)];
							energy_T3+=charge.valence*surf_charge_valence_T3*ENERGY_C[int(distance)];
							
							pot_T1+=surf_charge_valence_T1*POTENTIAL_C[int(distance)];
							pot_T2+=surf_charge_valence_T2*POTENTIAL_C[int(distance)];
							pot_T3+=surf_charge_valence_T3*POTENTIAL_C[int(distance)];
							
						}	
						else{
							// cut-off
						}			
					}
				}
			}

			
			
			
			 fout1 << 1e10*charge.z << " " 
			<< pot_T1 << " " << pot_T2 << " " << pot_T3<< " " 
			<< energy_T1 << " " << energy_T2 << " " << energy_T3<< " " 
			<< force_T1.at(0) << " " << force_T1.at(1) << " " << force_T1.at(2) << " " 
			<< force_T2.at(0) << " " << force_T2.at(1) << " " << force_T2.at(2) << " " 
			<< force_T3.at(0) << " " << force_T3.at(1) << " " << force_T3.at(2) << " " 
			<<endl;
			
			
		}
		
		fout1.close();
		
		cout << "OK!"<<endl;
		
	}
		
	
	free(vector_c_T1);
	free(vector_h_T1);
	free(vector_c_T2);
	free(vector_h_T2);
	free(vector_c_T3);
	free(vector_h_T3);
	
	
	return;
}


void compute_potential_on_2D_map(){
	
	string map_file= PRM.PREFIX + "_pot_map.dat";
	
	cout << "writing " <<map_file<<".........."<<flush;
	
	createFile(map_file);
	char *tfn1 = new char[map_file.length()+1];
	strcpy(tfn1, map_file.c_str());     
	ofstream fout1(tfn1, ios::app);


	double *vector_c_T1, *vector_h_T1;
	
						

// compute induced charges		
	if(NUM_OF_SURFACES>0){
		if(PRM.EPS_W!=PRM.EPS_MEM){
			
			posix_memalign((void **)(&vector_c_T1), 128, sizeof(double)*surfaces.size());
			posix_memalign((void **)(&vector_h_T1), 128, sizeof(double)*surfaces.size());
			
			if((!vector_c_T1) || (!vector_h_T1)){
				cerr<<"ERROR: Allocation of vector_c/h _T1 failed"<<endl;
				exit(1);
			}
			
			for(int index_on_h=0; index_on_h<surfaces.size(); index_on_h++){
				vector_c_T1[index_on_h]=_double_zero;
				vector_h_T1[index_on_h]=_double_zero;
			}

			double elem_c_T1=0.00;
			
			for(int surface_index=0; surface_index<surfaces.size(); surface_index++){

		
		// external field and membrane charges		
				
				vector_c_T1[surface_index]=C_BASE[surface_index];
			}
			
			
			
// solve the ICC equation 		
			dgemv_(&_N, &matrix_size, &matrix_size, &_double_one, inverse_a, &matrix_size, vector_c_T1, &_int_one, &_double_zero, vector_h_T1, &_int_one);	
		}
	}
			
			
	
	Point p;
	
	for(int z=PRM.MIN_Z; z<=PRM.MAX_Z; z=z+50){
		Ion charge(4);
		charge.x=0.00;
		charge.y=0.00;
		charge.z=1e-12*double(z);
		
		for(int r=0; r<=PRM.SIM_DOMAIN_WIDTH_X/2; r=r+50){
			charge.x=1e-12*double(r);
			
			double tempAA=0.00;
			double tempBB=0.00;	
			tempAA=(PRM.MAX_Z-1e12*charge.z)/(PRM.SIM_DOMAIN_WIDTH_Z);
			tempBB=PRM.APPLIED_POTENTIAL;	
			
			double pot_T1=tempAA*tempBB;
			
			double energy_T1=pot_T1*charge.charge;
			
			vector <double> force_T1;
			
			for(int iff=0; iff<3; iff++){
				force_T1.push_back(0.00);
			}
			
			// external field
			force_T1.at(2)=0.00;
			force_T1.at(2)=Q*charge.valence*PRM.APPLIED_FIELD;
			
			
			//membrane charges
			if(!membrane_charges.empty()){
				for(int charge_index=0; charge_index<membrane_charges.size(); charge_index++){
					double charge_x=membrane_charges.at(charge_index).x;
					double charge_y=membrane_charges.at(charge_index).y;
					double charge_z=membrane_charges.at(charge_index).z;			
					double distance=1e12*get_distance(charge.x, charge.y, charge.z, charge_x, charge_y, charge_z);			
					
					if(distance<20000){
						double force_C=0.00;
						double force_SR=0.00;
						double dx=0.00;
						double dy=0.00;
						double dz=0.00;
						
						dx=1e12*(charge.x-charge_x)/distance;
						dy=1e12*(charge.y-charge_y)/distance;
						dz=1e12*(charge.z-charge_z)/distance;
						
						force_C=charge.valence*membrane_charges.at(charge_index).DW_valence*FORCE_C[int(distance)];
						
						force_T1.at(0)+=(force_C)*dx;
						force_T1.at(1)+=(force_C)*dy;
						force_T1.at(2)+=(force_C)*dz;
						
						energy_T1+=charge.valence*membrane_charges.at(charge_index).DW_valence*ENERGY_C[int(distance)];
						pot_T1+=membrane_charges.at(charge_index).DW_valence*POTENTIAL_C[int(distance)];
					}	
					else{
						// cut-off
					}		
				}		
			}
		
			// boundary repulsion and induced charges
			if(NUM_OF_SURFACES>0){
				if(PRM.EPS_W!=PRM.EPS_MEM){
					for(int surface_index=0; surface_index<surfaces.size(); surface_index++){		
						double charge_x=surfaces.at(surface_index).center.x;
						double charge_y=surfaces.at(surface_index).center.y;
						double charge_z=surfaces.at(surface_index).center.z;				
						double distance=1e12*get_distance(charge.x, charge.y, charge.z, charge_x, charge_y, charge_z);	

						if(distance<20000){
							double force_C=0.00;
							double force_SR=0.00;
							double dx=0.00;
							double dy=0.00;
							double dz=0.00;
							
							dx=1e12*(charge.x-charge_x)/distance;
							dy=1e12*(charge.y-charge_y)/distance;
							dz=1e12*(charge.z-charge_z)/distance;
							
							double surf_charge_valence_T1=(vector_h_T1[surface_index]*surfaces.at(surface_index).area)/Q;
							double force_C_T1=charge.valence*surf_charge_valence_T1*FORCE_C[int(distance)];
							
							force_T1.at(0)+=(force_C_T1)*dx;
							force_T1.at(1)+=(force_C_T1)*dy;
							force_T1.at(2)+=(force_C_T1)*dz;
							
							energy_T1+=charge.valence*surf_charge_valence_T1*ENERGY_C[int(distance)];
							pot_T1+=surf_charge_valence_T1*POTENTIAL_C[int(distance)];
						}	
						else{
							// cut-off
						}			
					}
				}
			}

			fout1 << 1e10*charge.z << " " << 1e10*charge.x << " "
			<< pot_T1 << " "
			<< energy_T1 << " " 
			<< force_T1.at(0) << " " << force_T1.at(1) << " " << force_T1.at(2) << " " 
			<<endl;
			
			
		}
		
		fout1 << endl;
	}
		
	fout1.close();
	
	if(NUM_OF_SURFACES>0){
		if(PRM.EPS_W!=PRM.EPS_MEM){
			free(vector_c_T1);
			free(vector_h_T1);
		}
	}
	cout << "OK!"<<endl;
	return;
}
