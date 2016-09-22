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
#include "ions_properties.h"
#include "ions_functions.h"
#include "utils.h"
#include "classes.h"
#include "physics_functions.h"
#include "file_functions.h"
#include "classes.h"
#include "physics_functions.h"
#include "sim_structures.h"




//########################################
// Class Parameters
//########################################
Parameters::Parameters(){
	reset_parameters();
}

Parameters::~Parameters(){
	
}

void Parameters::reset_parameters(){
	HISTORY_SIZE=50;
	//~ HISTORY_SIZE=4;
	SEED=0;

	SIM_TYPE="NULL";
	PREFIX="NULL";	
	
	SR_METHOD="NULL";
	
	FIXED_CHARGES_FILE="NULL";

	MD_MAP_MIN_Z=0;
	MD_MAP_MAX_Z=0;

	PBC="NULL";
	ION_RECYCLING="NULL";
	SHORT_RANGE_EXP=12.00;
	APPLIED_POTENTIAL=0.00;
	APPLIED_FIELD=0.00;
	
	CONTROL_CELL_WIDTH=0.00;
	BATH_WIDTH=0.00;
	OUTER_REGION_WIDTH=0.00;
	MEMBRANE_WIDTH=0.00;
	SIM_DOMAIN_WIDTH_X=0.00;
	SIM_DOMAIN_WIDTH_Y=0.00;
	SIM_DOMAIN_WIDTH_Z=0.00;
	MIN_X=0.00;
	MIN_Y=0.00;
	MIN_Z=0.00;
	MAX_X=0.00;
	MAX_Y=0.00;
	MAX_Z=0.00;
	SIM_DOMAIN_VOLUME=0.00;
	SIM_DOMAIN_VOLUME=0.00;
	LEFT_CELL_MIN_Z=0.00;
	LEFT_CELL_MAX_Z=0.00;
	RIGHT_CELL_MIN_Z=0.00;
	RIGHT_CELL_MAX_Z=0.00;
	CONTROL_CELL_VOLUME=0.00;
	
	Z_BEGIN_SR_MAP=0.00;
	Z_END_SR_MAP=0.00;
	DELTA_SR_MAP=25;
	
	EPS_W=0.00;
	EPS_MEM=0.00;
	
	CONSTANT_DIFF_COEFF=true;
	DIFF_COEFF_IN_CHANNEL=1;
	
	LEFT_VESTIBULE_CURVATURE_RADIUS=0.00;
	LEFT_VESTIBULE_MIN_CHANNEL_RADIUS=0.00;
	//~ LEFT_VESTIBULE_MAX_CHANNEL_RADIUS=0.00;
	
	RIGHT_VESTIBULE_CURVATURE_RADIUS=0.00;
	RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS=0.00;
	//~ RIGHT_VESTIBULE_MAX_CHANNEL_RADIUS=0.00;
	
	TILES_PER_RING=10;
	NUM_OF_DIV=10;
	NUM_OF_SUB_DIV=100;
	MAX_TILE_WIDTH=300;
	
	delta_epsilon=0.00;
	mean_epsilon=0.00;		
	dielectrics_weight=0.00;
	
	c_0=1e-9;
	CONC_LEFT_KCL=0.00;
	CONC_LEFT_NACL=0.00;
	CONC_LEFT_CACL2=0.00;
	CONC_LEFT_MGCL2=0.00;
	CONC_RIGHT_KCL=0.00;	
	CONC_RIGHT_NACL=0.00;
	CONC_RIGHT_CACL2=0.00;
	CONC_RIGHT_MGCL2=0.00;
	
	PREP_STEPS=0;
	SIM_STEPS=0;
	FIRST_STEP=0;
	step=0;
	
	DELTA_T=0.00; //[fs]
	TEMPERATURE=0.00;
	
// ion boxes
	BOX1_MIN_X=0.00;
	BOX1_MAX_X=0.00;
	BOX1_MIN_Y=0.00;
	BOX1_MAX_Y=0.00;
	BOX1_MIN_Z=0.00;
	BOX1_MAX_Z=0.00;
	BOX1_IS1="NULL";
	BOX1_IS2="NULL";
	BOX1_IS3="NULL";
	BOX1_IS4="NULL";
	BOX1_IS5="NULL";
	BOX1_IS6="NULL";
	BOX1_N1=0;
	BOX1_N2=0;
	BOX1_N3=0;
	BOX1_N4=0;
	BOX1_N5=0;
	BOX1_N6=0;
	
	BOX2_MIN_X=0.00;
	BOX2_MAX_X=0.00;
	BOX2_MIN_Y=0.00;
	BOX2_MAX_Y=0.00;
	BOX2_MIN_Z=0.00;
	BOX2_MAX_Z=0.00;
	BOX2_IS1="NULL";
	BOX2_IS2="NULL";
	BOX2_IS3="NULL";
	BOX2_IS4="NULL";
	BOX2_IS5="NULL";
	BOX2_IS6="NULL";
	BOX2_N1=0;
	BOX2_N2=0;
	BOX2_N3=0;
	BOX2_N4=0;
	BOX2_N5=0;
	BOX2_N6=0;

	BOX3_MIN_X=0.00;
	BOX3_MAX_X=0.00;
	BOX3_MIN_Y=0.00;
	BOX3_MAX_Y=0.00;
	BOX3_MIN_Z=0.00;
	BOX3_MAX_Z=0.00;
	BOX3_IS1="NULL";
	BOX3_IS2="NULL";
	BOX3_IS3="NULL";
	BOX3_IS4="NULL";
	BOX3_IS5="NULL";
	BOX3_IS6="NULL";
	BOX3_N1=0;
	BOX3_N2=0;
	BOX3_N3=0;
	BOX3_N4=0;
	BOX3_N5=0;
	BOX3_N6=0;
	
	BOX4_MIN_X=0.00;
	BOX4_MAX_X=0.00;
	BOX4_MIN_Y=0.00;
	BOX4_MAX_Y=0.00;
	BOX4_MIN_Z=0.00;
	BOX4_MAX_Z=0.00;
	BOX4_IS1="NULL";
	BOX4_IS2="NULL";
	BOX4_IS3="NULL";
	BOX4_IS4="NULL";
	BOX4_IS5="NULL";
	BOX4_IS6="NULL";
	BOX4_N1=0;
	BOX4_N2=0;
	BOX4_N3=0;
	BOX4_N4=0;
	BOX4_N5=0;
	BOX4_N6=0;
	
	BOX5_MIN_X=0.00;
	BOX5_MAX_X=0.00;
	BOX5_MIN_Y=0.00;
	BOX5_MAX_Y=0.00;
	BOX5_MIN_Z=0.00;
	BOX5_MAX_Z=0.00;
	BOX5_IS1="NULL";
	BOX5_IS2="NULL";
	BOX5_IS3="NULL";
	BOX5_IS4="NULL";
	BOX5_IS5="NULL";
	BOX5_IS6="NULL";
	BOX5_N1=0;
	BOX5_N2=0;
	BOX5_N3=0;
	BOX5_N4=0;
	BOX5_N5=0;
	BOX5_N6=0;
	
	BOX6_MIN_X=0.00;
	BOX6_MAX_X=0.00;
	BOX6_MIN_Y=0.00;
	BOX6_MAX_Y=0.00;
	BOX6_MIN_Z=0.00;
	BOX6_MAX_Z=0.00;
	BOX6_IS1="NULL";
	BOX6_IS2="NULL";
	BOX6_IS3="NULL";
	BOX6_IS4="NULL";
	BOX6_IS5="NULL";
	BOX6_IS6="NULL";
	BOX6_N1=0;
	BOX6_N2=0;
	BOX6_N3=0;
	BOX6_N4=0;
	BOX6_N5=0;
	BOX6_N6=0;
	
//protein charges	
	charge_ring_z.clear();
	charge_ring_r.clear();
	charge_ring_n.clear();
	charge_ring_q.clear();
	
//statistics...	
	STATS_DZ=10;
	STATS_OUT_FREQ=1;

	channel_pdb_files=false;
	trajectory=false;
	flux=false;
	rdf=false;
	vel_distribution=false;
	induced_charge=false;
	potential=0;
	concentrations=0;
	
	channel_configuration=false;
	filter_configuration=false;

	currents_ZT=false;
	
	long_flight=0;
	into_membrane=0;
	out_of_SF=0;
	out_of_domain=0;
	
	ions_to_simulate.clear();
	return;
}
	
ostream& operator<<(ostream& stream, Parameters& PRM){
	
	stream << endl << "########################################" <<endl;
	stream << "    SIMULATION PARAMETERS";
	stream << endl << "########################################" <<endl;
	
//============== PAGE 2 - SIM. DOMAIN	
	stream << "\nSIMULATION TYPE: " << "\t" << PRM.SIM_TYPE << endl;  
	stream << "PREFIX: " << "\t" << PRM.PREFIX << endl <<endl;  
	
	stream << "SIMULATION DOMAIN WIDTH - X [pm]: " << "\t" << PRM.SIM_DOMAIN_WIDTH_X << endl;  
	stream << "SIMULATION DOMAIN WIDTH - Y [pm]: " << "\t" << PRM.SIM_DOMAIN_WIDTH_Y << endl;  
	stream << "CONTROL CELL WIDTH [pm]: " << "\t\t" << PRM.CONTROL_CELL_WIDTH << endl;  
	stream << "BATH WIDTH [pm]: " << "\t\t" << PRM.BATH_WIDTH << endl;  
	if(PRM.SIM_TYPE.compare("BULK")!=0){
		stream << "MEMBRANE WIDTH [pm]: " << "\t\t\t" << PRM.MEMBRANE_WIDTH << endl;  
	}
	if(PRM.SIM_TYPE.compare("PORE")==0){
		stream << "LEFT VEST. CURV. RADIUS [pm]: " << "\t" << PRM.LEFT_VESTIBULE_CURVATURE_RADIUS << endl;  
		stream << "LEFT VEST. MIN. CHANNEL RADIUS [pm]: " << "\t" << PRM.LEFT_VESTIBULE_MIN_CHANNEL_RADIUS << endl;  
		stream << "RIGHT VEST. CURV. RADIUS [pm]: " << "\t" << PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS << endl;  
		stream << "RIGHT VEST. MIN. CHANNEL RADIUS [pm]: " << "\t" << PRM.RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS << endl;  
	}
	
//============== PAGE 3 - MOTION	
	stream << endl << "DELTA T [s]:\t\t"<< PRM.DELTA_T << endl;
	stream << "SIMULATION STEPS:\t\t"<< PRM.SIM_STEPS << endl;
	stream << "SIMULATION TIME [s]:\t\t"<< PRM.DELTA_T*PRM.SIM_STEPS << endl;
	stream << "TRANSIENT STEPS:\t\t"<< PRM.PREP_STEPS << endl;
	stream << "TRANSIENT TIME [s]:\t\t"<< PRM.DELTA_T*PRM.PREP_STEPS << endl;
	if(PRM.SIM_TYPE.compare("PORE")==0){
		stream << "DIFFUSION COEFF IN THE CHANNEL:\t\t"<< PRM.DIFF_COEFF_IN_CHANNEL << endl;
	}
	stream << "SR METHOD:\t\t"<< PRM.SR_METHOD << endl;
	if(PRM.SR_METHOD.compare("EXPONENTIAL")==0){
		stream << "SR EXPONENT:\t\t"<< PRM.SHORT_RANGE_EXP << endl;
	}
	stream << "TEMPERATURE:\t\t"<< PRM.TEMPERATURE << endl;
	stream << "RANDOM NUMBER GEN. SEED:\t\t"<< PRM.SEED << endl;
	
	
	
//============== PAGE 4 - ELECTROSTATICS	
	stream << endl << "APPLIED POTENTIAL [V]:\t\t"<< PRM.APPLIED_POTENTIAL << endl;
	stream << "EPS_W:\t\t"<< PRM.EPS_W << endl;
	if(PRM.SIM_TYPE.compare("BULK")!=0 && PRM.EPS_W!=PRM.EPS_MEM){
		stream << "EPS_MEM:\t\t"<< PRM.EPS_MEM << endl;
		
		stream << "ICC ACCURACY - RADIAL:\t\t"<< PRM.TILES_PER_RING << endl;
		if(PRM.SIM_TYPE.compare("PORE")==0){
			stream << "ICC ACCURACY - AXIAL:\t\t"<< PRM.NUM_OF_DIV << endl;
		}
		stream << "ICC - SUBTILES PER TILE:\t\t"<< PRM.NUM_OF_SUB_DIV << endl;
	}
	
	

//============== PAGE 5 - ION CONCENTRATIONS	
	stream << endl << "LEFT KCL CONC [M]:\t\t"<< PRM.CONC_LEFT_KCL << endl;	
	stream << "LEFT NACL CONC [M]:\t\t"<< PRM.CONC_LEFT_NACL << endl;	
	stream << "LEFT CACL2 CONC [M]:\t\t"<< PRM.CONC_LEFT_CACL2 << endl;	
	stream << "LEFT MGCL2 CONC [M]:\t\t"<< PRM.CONC_LEFT_MGCL2 << endl;	
	stream << "RIGHT KCL CONC [M]:\t\t"<< PRM.CONC_RIGHT_KCL << endl;	
	stream << "RIGHT NACL CONC [M]:\t\t"<< PRM.CONC_RIGHT_NACL << endl;	
	stream << "RIGHT CACL2 CONC [M]:\t\t"<< PRM.CONC_RIGHT_CACL2 << endl;	
	stream << "RIGHT MGCL2 CONC [M]:\t\t"<< PRM.CONC_RIGHT_MGCL2 << endl;	
	
	stream << "PERIODIC BOUNDARY CONDITIONS:\t"<< PRM.PBC << endl;
	stream << "ION RECYCLING:\t"<< PRM.ION_RECYCLING << endl;
	
//============== PAGE 6 - IONS	
	stream << endl;
	for(int i=0; i<31; i++){
		if(PRM.ions_to_simulate.at(i)){
			Ion ion(i);
			stream << "--- ion ----------------------------------------------------------- ion ---" << endl;
			stream << "ION NAME - KIND:\t"<< ion.name << " - " << ion.kind << endl;	
			stream << "ION VALENCE - CHARGE:\t"<< ion.valence << " - " << ion.charge << endl;	
			stream << "ION MASS - RADIUS:\t"<< ion.mass << " - " << ion.radius << endl;	
			stream << "ION DIFFUSION COEFF.:\t"<< ion.diffusion_coeff << endl <<endl;	
		}
	}
	
	
	
//============== PAGE 7 - ION BOXES
	if(PRM.SIM_TYPE.compare("PORE")==0){
		stream << endl;
		if(PRM.BOX1_MIN_Z<PRM.BOX1_MAX_Z){
			stream << "BOX1 - MIN Z [pm]:\t"<< PRM.BOX1_MIN_Z << endl;
			stream << "BOX1 - MAX Z [pm]:\t"<< PRM.BOX1_MAX_Z << endl;
			if(PRM.BOX1_N1>0){				
				stream << "BOX1 - ION1:\tJ"<< "11" << endl;
				stream << "BOX1 - NUMBER OF ION J11:\t"<< PRM.BOX1_N1 << endl;
			}
			if(PRM.BOX1_N2>0){				
				stream << "BOX1 - ION2:\tJ"<< "12" << endl;
				stream << "BOX1 - NUMBER OF ION J12:\t"<< PRM.BOX1_N2 << endl;
			}
			if(PRM.BOX1_N3>0){				
				stream << "BOX1 - ION3:\tJ"<< "13" << endl;
				stream << "BOX1 - NUMBER OF ION J13:\t"<< PRM.BOX1_N3 << endl;
			}
			if(PRM.BOX1_N4>0){				
				stream << "BOX1 - ION4:\tJ"<< "14" << endl;
				stream << "BOX1 - NUMBER OF ION J14:\t"<< PRM.BOX1_N4 << endl;
			}
			if(PRM.BOX1_N5>0){				
				stream << "BOX1 - ION5:\tJ"<< "15" << endl;
				stream << "BOX1 - NUMBER OF ION J15:\t"<< PRM.BOX1_N5 << endl;
			}
			if(PRM.BOX1_N6>0){				
				stream << "BOX1 - ION6:\tJ"<< "16" << endl;
				stream << "BOX1 - NUMBER OF ION J16:\t"<< PRM.BOX1_N6 << endl;
			}
		}
		if(PRM.BOX2_MIN_Z<PRM.BOX2_MAX_Z){
			stream << endl;
			stream << "BOX2 - MIN Z [pm]:\t"<< PRM.BOX2_MIN_Z << endl;
			stream << "BOX2 - MAX Z [pm]:\t"<< PRM.BOX2_MAX_Z << endl;
			if(PRM.BOX2_N1>0){				
				stream << "BOX2 - ION1:\tJ"<< "21" << endl;
				stream << "BOX2 - NUMBER OF ION J21:\t"<< PRM.BOX2_N1 << endl;
			}
			if(PRM.BOX2_N2>0){				
				stream << "BOX2 - ION2:\tJ"<< "22" << endl;
				stream << "BOX2 - NUMBER OF ION J22:\t"<< PRM.BOX2_N2 << endl;
			}
			if(PRM.BOX2_N3>0){				
				stream << "BOX2 - ION3:\tJ"<< "23" << endl;
				stream << "BOX2 - NUMBER OF ION J23:\t"<< PRM.BOX2_N3 << endl;
			}
			if(PRM.BOX2_N4>0){				
				stream << "BOX2 - ION4:\tJ"<< "24" << endl;
				stream << "BOX2 - NUMBER OF ION J24:\t"<< PRM.BOX2_N4 << endl;
			}
			if(PRM.BOX2_N5>0){				
				stream << "BOX2 - ION5:\tJ"<< "25" << endl;
				stream << "BOX2 - NUMBER OF ION J25:\t"<< PRM.BOX2_N5 << endl;
			}
			if(PRM.BOX2_N6>0){				
				stream << "BOX2 - ION6:\tJ"<< "26" << endl;
				stream << "BOX2 - NUMBER OF ION J26:\t"<< PRM.BOX2_N6 << endl;
			}
		}
		if(PRM.BOX3_MIN_Z<PRM.BOX3_MAX_Z){
			stream << endl;
			stream << "BOX3 - MIN Z [pm]:\t"<< PRM.BOX3_MIN_Z << endl;
			stream << "BOX3 - MAX Z [pm]:\t"<< PRM.BOX3_MAX_Z << endl;
			if(PRM.BOX3_N1>0){				
				stream << "BOX3 - ION1:\tJ"<< "31" << endl;
				stream << "BOX3 - NUMBER OF ION J31:\t"<< PRM.BOX3_N1 << endl;
			}
			if(PRM.BOX3_N2>0){				
				stream << "BOX3 - ION2:\tJ"<< "32" << endl;
				stream << "BOX3 - NUMBER OF ION J32:\t"<< PRM.BOX3_N2 << endl;
			}
			if(PRM.BOX3_N3>0){				
				stream << "BOX3 - ION3:\tJ"<< "33" << endl;
				stream << "BOX3 - NUMBER OF ION J33:\t"<< PRM.BOX3_N3 << endl;
			}
			if(PRM.BOX3_N4>0){				
				stream << "BOX3 - ION4:\tJ"<< "34" << endl;
				stream << "BOX3 - NUMBER OF ION J34:\t"<< PRM.BOX3_N4 << endl;
			}
			if(PRM.BOX3_N5>0){				
				stream << "BOX3 - ION5:\tJ"<< "35" << endl;
				stream << "BOX3 - NUMBER OF ION J35:\t"<< PRM.BOX3_N5 << endl;
			}
			if(PRM.BOX3_N6>0){				
				stream << "BOX3 - ION6:\tJ"<< "36" << endl;
				stream << "BOX3 - NUMBER OF ION J36:\t"<< PRM.BOX3_N6 << endl;
			}
		}
		if(PRM.BOX4_MIN_Z<PRM.BOX4_MAX_Z){
			stream << endl;
			stream << "BOX4 - MIN Z [pm]:\t"<< PRM.BOX4_MIN_Z << endl;
			stream << "BOX4 - MAX Z [pm]:\t"<< PRM.BOX4_MAX_Z << endl;
			if(PRM.BOX4_N1>0){				
				stream << "BOX4 - ION1:\tJ"<< "41" << endl;
				stream << "BOX4 - NUMBER OF ION J41:\t"<< PRM.BOX4_N1 << endl;
			}
			if(PRM.BOX4_N2>0){				
				stream << "BOX4 - ION2:\tJ"<< "42" << endl;
				stream << "BOX4 - NUMBER OF ION J42:\t"<< PRM.BOX4_N2 << endl;
			}
			if(PRM.BOX4_N3>0){				
				stream << "BOX4 - ION3:\tJ"<< "43" << endl;
				stream << "BOX4 - NUMBER OF ION J43:\t"<< PRM.BOX4_N3 << endl;
			}
			if(PRM.BOX4_N4>0){				
				stream << "BOX4 - ION4:\tJ"<< "44" << endl;
				stream << "BOX4 - NUMBER OF ION J44:\t"<< PRM.BOX4_N4 << endl;
			}
			if(PRM.BOX4_N5>0){				
				stream << "BOX4 - ION5:\tJ"<< "45" << endl;
				stream << "BOX4 - NUMBER OF ION J45:\t"<< PRM.BOX4_N5 << endl;
			}
			if(PRM.BOX4_N6>0){				
				stream << "BOX4 - ION6:\tJ"<< "46" << endl;
				stream << "BOX4 - NUMBER OF ION J46:\t"<< PRM.BOX4_N6 << endl;
			}
		}
		if(PRM.BOX5_MIN_Z<PRM.BOX5_MAX_Z){
			stream << endl;
			stream << "BOX5 - MIN Z [pm]:\t"<< PRM.BOX5_MIN_Z << endl;
			stream << "BOX5 - MAX Z [pm]:\t"<< PRM.BOX5_MAX_Z << endl;
			if(PRM.BOX5_N1>0){				
				stream << "BOX5 - ION1:\tJ"<< "51" << endl;
				stream << "BOX5 - NUMBER OF ION J51:\t"<< PRM.BOX5_N1 << endl;
			}
			if(PRM.BOX5_N2>0){				
				stream << "BOX5 - ION2:\tJ"<< "52" << endl;
				stream << "BOX5 - NUMBER OF ION J52:\t"<< PRM.BOX5_N2 << endl;
			}
			if(PRM.BOX5_N3>0){				
				stream << "BOX5 - ION3:\tJ"<< "53" << endl;
				stream << "BOX5 - NUMBER OF ION J53:\t"<< PRM.BOX5_N3 << endl;
			}
			if(PRM.BOX5_N4>0){				
				stream << "BOX5 - ION4:\tJ"<< "54" << endl;
				stream << "BOX5 - NUMBER OF ION J54:\t"<< PRM.BOX5_N4 << endl;
			}
			if(PRM.BOX5_N5>0){				
				stream << "BOX5 - ION5:\tJ"<< "55" << endl;
				stream << "BOX5 - NUMBER OF ION J55:\t"<< PRM.BOX5_N5 << endl;
			}
			if(PRM.BOX5_N6>0){				
				stream << "BOX5 - ION6:\tJ"<< "56" << endl;
				stream << "BOX5 - NUMBER OF ION J56:\t"<< PRM.BOX5_N6 << endl;
			}
		}
		if(PRM.BOX6_MIN_Z<PRM.BOX6_MAX_Z){
			stream << endl;
			stream << "BOX6 - MIN Z [pm]:\t"<< PRM.BOX6_MIN_Z << endl;
			stream << "BOX6 - MAX Z [pm]:\t"<< PRM.BOX6_MAX_Z << endl;
			if(PRM.BOX6_N1>0){				
				stream << "BOX6 - ION1:\tJ"<< "61" << endl;
				stream << "BOX6 - NUMBER OF ION J61:\t"<< PRM.BOX6_N1 << endl;
			}
			if(PRM.BOX6_N2>0){				
				stream << "BOX6 - ION2:\tJ"<< "62" << endl;
				stream << "BOX6 - NUMBER OF ION J62:\t"<< PRM.BOX6_N2 << endl;
			}
			if(PRM.BOX6_N3>0){				
				stream << "BOX6 - ION3:\tJ"<< "63" << endl;
				stream << "BOX6 - NUMBER OF ION J63:\t"<< PRM.BOX6_N3 << endl;
			}
			if(PRM.BOX6_N4>0){				
				stream << "BOX6 - ION4:\tJ"<< "64" << endl;
				stream << "BOX6 - NUMBER OF ION J64:\t"<< PRM.BOX6_N4 << endl;
			}
			if(PRM.BOX6_N5>0){				
				stream << "BOX6 - ION5:\tJ"<< "65" << endl;
				stream << "BOX6 - NUMBER OF ION J65:\t"<< PRM.BOX6_N5 << endl;
			}
			if(PRM.BOX6_N6>0){				
				stream << "BOX6 - ION6:\tJ"<< "66" << endl;
				stream << "BOX6 - NUMBER OF ION J66:\t"<< PRM.BOX6_N6 << endl;
			}
		}

		stream << endl;
		for(int i=31; i<92; i++){
			if(PRM.ions_to_simulate.at(i)){
				Ion ion(i);
				stream << "--- ion ----------------------------------------------------------- ion ---" << endl;
				stream << "ION NAME - KIND:\t"<< ion.name << " - " << ion.kind << endl;	
				stream << "ION VALENCE - CHARGE:\t"<< ion.valence << " - " << ion.charge << endl;	
				stream << "ION MASS - RADIUS:\t"<< ion.mass << " - " << ion.radius << endl;	
				stream << "ION DIFFUSION COEFF.:\t"<< ion.diffusion_coeff << endl <<endl;	
			}
		}
	}
	
	
//============== PAGE 8 - FIXED CHARGES
	if(PRM.SIM_TYPE.compare("PORE")==0){
		stream << endl;
		for(int i=0; i<PRM.charge_ring_z.size(); i++){
			stream << "--- charge ring ----------------------------------------- charge ring  ---" << endl;
			stream << "CHARGE RING - z [pm]:\t"<< PRM.charge_ring_z.at(i) << endl;	
			stream << "CHARGE RING - r [pm]:\t"<< PRM.charge_ring_r.at(i) << endl;	
			stream << "NUMBER OF CHARGES:\t"<< PRM.charge_ring_n.at(i) << endl;	
			stream << "VALENCE:\t"<< PRM.charge_ring_q.at(i) << endl <<endl;
		}
	}
	
	
	
	
	stream << endl << endl;
	
	
	return stream;
}


//########################################
// Class Statistics
//########################################
Statistics::Statistics(){
	return;
}

Statistics::~Statistics(){
}

void Statistics::reset_statistics(){
	ERROR_nan=0;
	ERROR_long_jump=0;
	ERROR_into_membrane=0;
	ERROR_ion_boxes=0;
	
	MAX_STEPS=0;
	msds_index_1=1;
	msds_index_2=0;
	
	ions_this_step.clear();
	
	concentrations_along_z.clear();

	RDF.clear();
	RDF_samples.clear();
	velocities.clear();
	msds.clear();
	msds_start.clear();
	currents_RS.clear();
	instant_currents_RS.clear();
	currents_ZT.clear();
	instant_currents_ZT.clear();

	A_total_charge=0;
	C_total_charge.clear();
	C_total_charge_type.clear();
	concs_3D.clear();
	radial_concs.clear();
	potentials_3D.clear();
	radial_potentials.clear();
	potentials_on_axis.clear();
	
	for(int a=0; a<3; a++){
		for(int b=0; b<3; b++){
			for(int c=0; c<3; c++){
				for(int d=0; d<3; d++){
					for(int e=0; e<3; e++){
						CHANNEL_CONFIGURATIONS[a][b][c][d][e]=0;
						FILTER_CONFIGURATIONS[a][b][c][d][e]=0;
						for(int f=0; f<100; f++){
							CHANNEL_CONFIGURATIONS_FREQUENCIES[a][b][c][d][e][f]=0;
							FILTER_CONFIGURATIONS_FREQUENCIES[a][b][c][d][e][f]=0;
						}
					}
				}
			}
		}
		last_output_channel_configuration[a]=0;
		last_output_filter_configuration[a]=0;
	}

	stat_file_1="";
	stat_file_2="";
	stat_file_3="";
	stat_file_4="";
	stat_file_5="";
	stat_file_6="";
	stat_file_7="";
	stat_file_8="";
	stat_file_9="";
	stat_file_B="";
	stat_file_C="";
	
	gnuplot_file_1="";
	gnuplot_file_2="";
	gnuplot_file_3="";
	gnuplot_file_4="";
	gnuplot_file_5="";
	gnuplot_file_6="";
	gnuplot_file_7="";
	gnuplot_file_8="";
	gnuplot_file_9="";
	gnuplot_file_B="";
	gnuplot_file_C="";
	
	save_status_file="";
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
	20		J00		
	21		J01		
	22		J02		
	23		J03		
	24		J04		
	25		J05		
	
	26		J06		
	27		J07		
	28		J08		
	29		J09		
	30		J10		
	31		J11		
	
	32		J12		
	33		J13		
	34		J14		
	35		J15		
	36		J16		
	37		J17		

	38		J18		
	39		J19		
	40		J20		
	41		J21		
	42		J22		
	43		J23		

	44		J24		
	45		J25		
	46		J26		
	47		J27		
	48		J28		
	49		J29		


	50		Q00		
	51		Q01		
	52		Q02		
	53		Q03		
	54		Q04		
	55		Q05		
	
	56		Q06		
	57		Q07		
	58		Q08		
	59		Q09		
	60		Q10		
	61		Q11		
	
	62		Q12		
	63		Q13		
	64		Q14		
	65		Q15		
	66		Q16		
	67		Q17		

	68		Q18		
	69		Q19		
	70		Q20		
	71		Q21		
	72		Q22		
	73		Q23		

	74		Q24		
	75		Q25		
	76		Q26		
	77		Q27		
	78		Q28		
	79		Q29		
	
	
	80		X00		
	81		X01		
	82		X02		
	83		X03		
	84		X04		
	85		X05		
	
	86		X06		
	87		X07		
	88		X08		
	89		X09		
	90		X10		
	91		X11		
*/
//===========================================	
  
	MAX_STEPS=PRM.SIM_STEPS;
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		ions_this_step.push_back(0);
	}
	vector <double> aux_vec_1;
	vector <double> aux_vec_2;

// ions trajectories
	if(PRM.trajectory){
		string trajectory_file=PRM.PREFIX + ".trajectory.dat";
		char *tfn1 = new char[trajectory_file.length()+1];
		strcpy(tfn1, trajectory_file.c_str());
		ofstream fout(tfn1, ios::out);
		fout.close();
	}
	
// flux
	if(PRM.flux){
		num_of_dz=int(PRM.SIM_DOMAIN_WIDTH_Z/PRM.STATS_DZ);
		DELTA_Z=PRM.SIM_DOMAIN_WIDTH_Z/double(num_of_dz);
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			aux_vec_1.clear();
			aux_vec_2.clear();
			if(PRM.ions_to_simulate.at(is)){
				for(int idz=0; idz<num_of_dz; idz++){
					aux_vec_1.push_back(0.00);
					aux_vec_2.push_back(0.00);
				}
			}
			concentrations_along_z.push_back(aux_vec_1);			
			aux_vec_1.clear();
			aux_vec_2.clear();
		}
	}
	
// radial distribution function
	if(PRM.rdf){
		create_RDF_vector(RDF, RDF_samples);
	}
	
// velocity distribution computation	
	if(PRM.vel_distribution){
		velocities_samples=0;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			aux_vec_1.clear();
			if(PRM.ions_to_simulate.at(is)){
				for(int iv=0; iv<150; iv++){
					aux_vec_1.push_back(0.00);
				}
			}
		
			velocities.push_back(aux_vec_1);
			
			aux_vec_1.clear();
		}
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			
			aux_vec_1.clear();
			
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				
				double t1=ion_curr.mass/PRM.kT;
				double t2=t1*t1*t1;
				double t3=double(2.00)*t2/M_PI;
				double t4=sqrt(t3);
				double t5=-t1/double(2.00);
				
				for(int iv=0; iv<150; iv++){
					double ivv=double(iv)*double(10.00);
					double t6=double(ivv)*double(ivv);
					double t7=t4*t6;
					double t8=exp(t5*t6);
					double f_v=t7*t8;
					aux_vec_1.push_back(f_v);
				}
			}
		
			velocities_theory.push_back(aux_vec_1);
			
			aux_vec_1.clear();
		}
	}
	
// mean square displacement computation	
	if(PRM.mean_square_displ){
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			aux_vec_1.clear();
			if(PRM.ions_to_simulate.at(is)){
				for(int iv=0; iv<PRM.STATS_OUT_FREQ; iv++){
					aux_vec_1.push_back(0.00);
				}
			}
			msds.push_back(aux_vec_1);
			aux_vec_1.clear();
			msds_num.push_back(0.00);
		}
	}	

// induced charge computation
	
// potential
	if(PRM.potential == 1){
		initialize_output_potential();
	}	
	
// concentrations computation
	if(PRM.concentrations!=0){
		
		double double_zero=0.00;
		vector < double > aux_vec1;  
		vector < vector < double > > aux_vec2;  
		vector < vector < vector < double > > > aux_vec3;  
		
		
		if(PRM.concentrations==3){
			int num_of_div_x=PRM.SIM_DOMAIN_WIDTH_X/double(100.00);
			int num_of_div_y=PRM.SIM_DOMAIN_WIDTH_Y/double(100.00);
			int num_of_div_z=PRM.SIM_DOMAIN_WIDTH_Z/double(100.00);	
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				aux_vec3.clear();
				if(PRM.ions_to_simulate.at(is)){
					for(int ix=0; ix<num_of_div_x; ix++){
						aux_vec2.clear();
						for(int iy=0; iy<num_of_div_y; iy++){
							aux_vec1.clear();
							for(int iz=0; iz<num_of_div_z; iz++){
								aux_vec1.push_back(double_zero);
							}
							aux_vec2.push_back(aux_vec1);
							aux_vec1.clear();
						}
						aux_vec3.push_back(aux_vec2);
						aux_vec2.clear();
					}
				}
				concs_3D.push_back(aux_vec3);
				aux_vec3.clear();
			}
		}
		else{
			int num_of_div_on_z_stat=PRM.SIM_DOMAIN_WIDTH_Z/double(100.00);
			int num_of_div_on_r_stat=PRM.SIM_DOMAIN_WIDTH_X/double(100.00);
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				aux_vec2.clear();
				if(PRM.ions_to_simulate.at(is)){
					for(int iz=0; iz<num_of_div_on_z_stat; iz++){
						aux_vec1.clear();
						for(int ir=0; ir<num_of_div_on_r_stat; ir++){
							aux_vec1.push_back(double_zero);	
						}
						aux_vec2.push_back(aux_vec1);
						aux_vec1.clear();
					}
				}
				radial_concs.push_back(aux_vec2);
				aux_vec2.clear();
			}
		}
	}
	
// currents computation	(ZERO THRESHOLD, one threshold at z=0)
	if(PRM.currents_ZT){	
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			currents_ZT.push_back(0.00);
			instant_currents_ZT.push_back(0.00);
		}
	}
	
// statistics files creation	
	stat_file_1=PRM.PREFIX + ".ion_pos";
	stat_file_2=PRM.PREFIX + ".flux";
	stat_file_3=PRM.PREFIX + ".rdf";
	stat_file_4=PRM.PREFIX + ".vel";

	stat_file_6=PRM.PREFIX + ".channel";
	stat_file_7=PRM.PREFIX + ".concs";
	stat_file_8=PRM.PREFIX + ".filter";
	stat_file_9=PRM.PREFIX + ".curr_ZT";
	
	stat_file_A=PRM.PREFIX + ".channel.freq";
	stat_file_B=PRM.PREFIX + ".msd";
	stat_file_C=PRM.PREFIX + ".pot";
	stat_file_D=PRM.PREFIX + ".filter.freq";
	
	gnuplot_file_1=stat_file_1 + ".plot";
	gnuplot_file_2=stat_file_2 + ".plot";
	gnuplot_file_3=stat_file_3 + ".plot";
	gnuplot_file_4=stat_file_4 + ".plot";
	gnuplot_file_5=stat_file_5 + ".plot";
	gnuplot_file_6=stat_file_6 + ".plot";
	gnuplot_file_7=stat_file_7 + ".plot";
	gnuplot_file_7=stat_file_8 + ".plot";
	gnuplot_file_9=stat_file_9 + ".plot";
	gnuplot_file_A=stat_file_A + ".plot";
	gnuplot_file_B=stat_file_B + ".plot";
	gnuplot_file_C=stat_file_C + ".plot";
	gnuplot_file_D=stat_file_D + ".plot";
	
	
	deleteFile(stat_file_1);
	deleteFile(stat_file_2);
	deleteFile(stat_file_3);
	deleteFile(stat_file_4);
	deleteFile(stat_file_5);
	deleteFile(stat_file_6);
	deleteFile(stat_file_7);
	deleteFile(stat_file_8);
	deleteFile(stat_file_9);
	deleteFile(stat_file_A);
	deleteFile(stat_file_B);
	deleteFile(stat_file_C);
	deleteFile(stat_file_D);
	
	deleteFile(gnuplot_file_1);
	deleteFile(gnuplot_file_2);
	deleteFile(gnuplot_file_3);
	deleteFile(gnuplot_file_4);
	deleteFile(gnuplot_file_5);
	deleteFile(gnuplot_file_6);
	deleteFile(gnuplot_file_7);
	deleteFile(gnuplot_file_8);
	deleteFile(gnuplot_file_9);
	deleteFile(gnuplot_file_A);
	deleteFile(gnuplot_file_B);
	deleteFile(gnuplot_file_C);
	deleteFile(gnuplot_file_D);
	
	int last=0;
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		if(PRM.ions_to_simulate.at(is)){
			last=is;
		}
	}
	
	if(PRM.flux){
		createFile(stat_file_2);
		createFile(gnuplot_file_2);
		char *tfn22 = new char[gnuplot_file_2.length()+1];
		strcpy(tfn22, gnuplot_file_2.c_str());     
		ofstream fout22(tfn22, ios::app);
		fout22 << "set encoding iso_8859_1" << endl;
		fout22 << "set term wxt 1"<<endl;
		fout22 << "set title 'Ion flux'" <<endl;
		fout22 << "set xlabel 'z [\305]'" <<endl;
		fout22 << "set ylabel 'flux [# ions/s]'" <<endl;
		fout22 << "set grid"<<endl;
		fout22 << "set xtics 10"<<endl;
		fout22 << "plot ";
		int column=2;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				if(is==last){
					fout22 << "'" << stat_file_2 << "' u 1:" << column << " w l t '" << ion_curr.name << " flux'";
				}
				else{
					fout22 << "'" << stat_file_2 << "' u 1:" << column << " w l t '" << ion_curr.name << " flux', ";
				}
				column+=3;
			}
		}
		fout22 << endl;
		fout22 << "set term wxt 2"<<endl;
		fout22 << "set title 'Ion concentration'" <<endl;
		fout22 << "set ylabel 'ion concentration [# ions/\305]'" <<endl;
		fout22 << "set grid"<<endl;
		fout22 << "set xtics 10"<<endl;
		fout22 << "plot ";
		column=3;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				if(is==last){
					fout22 << "'" << stat_file_2 << "' u 1:" << column << " w l t '" << ion_curr.name << " concentration'";
				}
				else{
					fout22 << "'" << stat_file_2 << "' u 1:" << column << " w l t '" << ion_curr.name << " concentration', ";
				}
				column+=3;
			}
		}
		fout22 << endl;
		double K_Jstar=1e14/(PRM.SIM_DOMAIN_WIDTH_X*PRM.SIM_DOMAIN_WIDTH_Y);
		double K_J=K_Jstar*1e17/AVOGADRO;
		double K_I=K_J*1e-21*PRM.SIM_DOMAIN_WIDTH_X*PRM.SIM_DOMAIN_WIDTH_Y*FARADAY_K;
		double K_flux=K_Jstar;
		string units="ions  A^-2  s^-1";
		if(PRM.SIM_TYPE.compare("BULK")==0){
			K_flux=K_J;
			units="mol  m^-2  s^-1";
		}
		fout22 << "set term wxt 3"<<endl;
		fout22 << "set title 'Ion velocity'" <<endl;
		fout22 << "set ylabel 'ion velocity [m/s]'" <<endl;
		fout22 << "set grid"<<endl;
		fout22 << "set xtics 10"<<endl;
		fout22 << "plot ";
		column=4;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				if(is==last){
					fout22 << "'" << stat_file_2 << "' u 1:" << column << " w l t '" << ion_curr.name << " velocity'";
				}
				else{
					fout22 << "'" << stat_file_2 << "' u 1:" << column << " w l t '" << ion_curr.name << " velocity', ";
				}
				column+=3;
			}
		}
		fout22 << endl;
		//~ fout22 <<"pause -1" <<endl;
		fout22 << "set term wxt 4"<<endl;
		fout22 << "set title 'Ion current'" <<endl;
		fout22 << "set ylabel 'ion current [A]'" <<endl;
		fout22 << "set grid"<<endl;
		fout22 << "set xtics 10"<<endl;
		fout22 << "plot ";
		column=2;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				if(is==last){
					fout22 << "'" << stat_file_2 << "' u 1:("<< ion_curr.charge <<"*($" << column << ")) w l t '" << ion_curr.name << " current'";
				}
				else{
					fout22 << "'" << stat_file_2 << "' u 1:("<< ion_curr.charge <<"*($" << column << ")) w l t '" << ion_curr.name << " current', ";
				}
				column+=3;
			}
		}
		fout22 << endl;
		fout22 <<"pause -1" <<endl;
		fout22.close();
	}
	
	if(PRM.rdf){
		createFile(stat_file_3);
		createFile(gnuplot_file_3);
		char *tfn33 = new char[gnuplot_file_3.length()+1];
		strcpy(tfn33, gnuplot_file_3.c_str());     
		ofstream fout33(tfn33, ios::app);
		fout33 << "set encoding iso_8859_1" << endl;
		fout33 << "set title 'Radial distribution function'" <<endl;
		fout33 << "set grid" <<endl;
		fout33 << "set xlabel 'ion-ion distance [\305]'" <<endl;
		fout33 << "set ylabel 'radial distribution function'" <<endl;
		fout33 << "set xrange [0:10]" <<endl;
		fout33 << "plot ";
		int column=2;
		for(int is1=0; is1<NUM_OF_IONIC_SPECIES; is1++){
			if(PRM.ions_to_simulate.at(is1)){
				Ion ion_curr1(is1);
				for(int is2=is1; is2<NUM_OF_IONIC_SPECIES; is2++){
					if(PRM.ions_to_simulate.at(is2)){
						Ion ion_curr2(is2);
						if(is1==last && is2==last){
							fout33 << "'" << stat_file_3 << "' u 1:" << column << " w l t '" << ion_curr1.name << "-" << ion_curr2.name<< " RDF'";
						}
						else{
							fout33 << "'" << stat_file_3 << "' u 1:" << column << " w l t '" << ion_curr1.name << "-" << ion_curr2.name<< " RDF', ";
						}
						column++;
					}
				}
			}
		}
		fout33 << endl;
		fout33 <<"pause -1" <<endl;
		fout33.close();
	}
	
	if(PRM.vel_distribution){
		createFile(stat_file_4);
		createFile(gnuplot_file_4);
		char *tfn44 = new char[gnuplot_file_4.length()+1];
		strcpy(tfn44, gnuplot_file_4.c_str());     
		ofstream fout44(tfn44, ios::app);
		fout44 << "set xrange [0:1500]"<<endl;
		fout44 << "set title 'Ion velocity distribution'" <<endl;
		fout44 << "set grid"<<endl;
		fout44 << "set xlabel 'ion velocity [m/s]'" <<endl;
		fout44 << "set ylabel 'velocity density function'" <<endl;
		fout44 << "plot ";
		int column=2;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				if(is==last){
					fout44 << "'" << stat_file_4 << "' u 1:" << column << " w l t '" << ion_curr.name << " f(v) - theory', '" << stat_file_4 << "' u 1:" << column+1 << " w l t '" << ion_curr.name << " f(v) - simulation'";
				}
				else{
					fout44 << "'" << stat_file_4 << "' u 1:" << column << " w l t '" << ion_curr.name << " f(v) - theory', '" << stat_file_4 << "' u 1:" << column+1 << " w l t '" << ion_curr.name << " f(v) - simulation', ";
				}
				column+=2;
			}
		}
		fout44 << endl;
		fout44 <<"pause -1" <<endl;
		fout44.close();
	}
	
	if(PRM.mean_square_displ){		
		createFile(stat_file_B);
		createFile(gnuplot_file_B);
		char *tfnBB = new char[gnuplot_file_B.length()+1];
		strcpy(tfnBB, gnuplot_file_B.c_str());     
		ofstream foutBB(tfnBB, ios::app);
		foutBB << "set encoding iso_8859_1" << endl;
		foutBB << "set title 'Mean square displacement'" <<endl;
		foutBB << "set key left"<<endl;
		foutBB << "set grid"<<endl;
		foutBB << "set xlabel 'time [ps]'" <<endl;
		foutBB << "set ylabel 'MSD [\305^2]'" <<endl;
		foutBB << "plot ";
		int column=2;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				if(is==last){
					foutBB << "'" << stat_file_B << "' u 1:" << column << " w l t '" << ion_curr.name << " MSD - theory', '" << stat_file_B << "' u 1:" << column+1 << " w l t '" << ion_curr.name << " MSD - simulation'";
				}
				else{
					foutBB << "'" << stat_file_B << "' u 1:" << column << " w l t '" << ion_curr.name << " MSD - theory', '" << stat_file_B << "' u 1:" << column+1 << " w l t '" << ion_curr.name << " MSD - simulation', ";
				}
				column+=2;
			}
		}
		foutBB << endl;
		foutBB <<"pause -1" <<endl;
		foutBB.close();
	}
	
	if(PRM.induced_charge){
		createFile(stat_file_6);
		createFile(gnuplot_file_6);
		char *tfn66 = new char[gnuplot_file_6.length()+1];
		strcpy(tfn66, gnuplot_file_6.c_str());     
		ofstream fout66(tfn66, ios::app);
		fout66 << "set grid"<<endl;
		fout66 << "set title 'Sum rule for dielectrics'" <<endl;
		fout66 << "set xlabel 'Poisson solution index'" <<endl;
		fout66 << "set ylabel 'Total charge on the boundary [e]'" <<endl;
		fout66 << "set key left" <<endl;
		fout66 << "plot '" << stat_file_6 << "' u 0:1 w l t 'Analytic', '" << stat_file_6 << "' u 0:2 w l t 'Computed', '" << stat_file_6 << "' u 0:3 w l t 'Type'" << endl;
		fout66 <<"pause -1" <<endl;
		fout66.close();
	}

	if(PRM.concentrations!=0){
		createFile(stat_file_7);
		createFile(gnuplot_file_7);
	}

	if(PRM.potential == 1){
		createFile(stat_file_C);
		createFile(gnuplot_file_C);
		char *tfnCC = new char[gnuplot_file_C.length()+1];
		strcpy(tfnCC, gnuplot_file_C.c_str());
		ofstream foutCC(tfnCC, ios::app);
		foutCC << "set encoding iso_8859_1" << endl;
		foutCC << "set term wxt 1"<<endl;
		foutCC << "set title 'AVERAGE POTENTIAL ALONG AXIS'" <<endl;
		foutCC << "set xlabel 'z [\305]'" <<endl;
		foutCC << "set ylabel 'Electric potential [mV]'" <<endl;
		foutCC << "set grid"<<endl;
		foutCC << "set xtics 10"<<endl;
		foutCC << "plot '" << stat_file_C << "' u 1:2 w l t 'potential'"<<endl;
		foutCC <<"pause -1" <<endl;
		foutCC.close();
	}

	if(PRM.currents_ZT){
		createFile(stat_file_9);
		createFile(gnuplot_file_9);
		char *tfn99 = new char[gnuplot_file_9.length()+1];
		strcpy(tfn99, gnuplot_file_9.c_str());     
		ofstream fout99(tfn99, ios::app);
		fout99 << "set term wxt 1"<<endl;
		fout99 << "set grid"<<endl;
		fout99 << "set title 'Ion current'" <<endl;
		fout99 << "set xlabel 'simulation time [s]'" <<endl;
		fout99 << "set ylabel 'Average current [A]'" <<endl;
		fout99 << "plot ";
		int column=3;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				fout99 << "'" << stat_file_9 << "' u 1:" << column << " w l t '" << ion_curr.name << " current', ";
				column+=2;
			}
		}
		fout99 << "'" << stat_file_9 << "' u 1:" << column << " w l t 'total current'";
		fout99 << endl;
		fout99 << "set term wxt 2"<<endl;
		fout99 << "set ylabel 'Instant current [A]'" <<endl;
		fout99 << "plot ";
		column=2;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				Ion ion_curr(is);
				fout99 << "'" << stat_file_9 << "' u 1:" << column << " w l t '" << ion_curr.name << " current', ";
				column+=2;
			}
		}
		fout99 << "'" << stat_file_9 << "' u 1:" << column << " w l t 'total current'";
		fout99 << endl;
		fout99 <<"pause -1" <<endl;
		fout99.close();
	}

	if(PRM.channel_configuration){
		createFile(stat_file_6);
	}

	if(PRM.filter_configuration){
		createFile(stat_file_8);
	}

	save_status_file=PRM.PREFIX + ".save";
	return;
}

void Statistics::update_statistics(){
	int num_of_div_on_z_stat=PRM.SIM_DOMAIN_WIDTH_Z/PRM.STATS_DZ;

	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		ions_this_step.at(i)=0;
	}
	for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){
		ions_this_step.at(IONS[INDEX_STAT_STEP][i].kind)++;
	}
	
// flux computation	
	if(PRM.flux){
		for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){
			double zzz=1e12*IONS[INDEX_STAT_STEP][i].z;
			if(zzz<PRM.MIN_Z){
				zzz=PRM.SIM_DOMAIN_WIDTH_Z+zzz;
			}
			if(zzz>PRM.MAX_Z){
				zzz=zzz-PRM.SIM_DOMAIN_WIDTH_Z;
			}
			int index=int((-PRM.MIN_Z+zzz)/DELTA_Z);
			concentrations_along_z.at(IONS[INDEX_STAT_STEP][i].kind).at(index)+=1.00;
		}
	}

// radial distribution function computation	
	if(PRM.rdf){
		update_RDF_vector(RDF, RDF_samples);
	}
	
// velocity distribution computation	
	if(PRM.vel_distribution){
		
		velocities_samples=velocities_samples+1.00;
		
		for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){
			double velocity=get_ion_velocity(IONS[INDEX_STAT_STEP][i]);
				
			int indexv=int(double(velocity/double(10.00)));
			if(isinf(velocity) || isnan(velocity)){
				cout << "velocity: "<<velocity<<endl;
			}
			else{
				if(indexv>=0 && indexv<150){
					velocities.at(IONS[INDEX_STAT_STEP][i].kind).at(indexv)+=double(1.00)/double(ions_this_step.at(IONS[INDEX_STAT_STEP][i].kind));
				}
			}
		}
	}

// mean square displacement computation	
	if(PRM.mean_square_displ){
		for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){
			double start_x=1e12*msds_start.at(i).x;
			double start_y=1e12*msds_start.at(i).y;
			double start_z=1e12*msds_start.at(i).z;
			double charge_x=1e12*IONS[INDEX_STAT_STEP][i].x_prev;
			double charge_y=1e12*IONS[INDEX_STAT_STEP][i].y_prev;
			double charge_z=1e12*IONS[INDEX_STAT_STEP][i].z_prev;
			apply_periodic_boundary(start_x, start_y, start_z, charge_x, charge_y, charge_z);
			double distance=1e-2*get_distance(start_x, start_y, start_z, charge_x, charge_y, charge_z);	//MSD in A^2
			double distance_2=distance*distance;
			msds.at(IONS[INDEX_STAT_STEP][i].kind).at(msds_index_2)+=distance_2;
		}
		msds_index_2++;
	}	
	
// concentrations computation	
	if(PRM.concentrations!=0){	
		if(PRM.concentrations==3){ // average on y	
			int num_of_div_x=PRM.SIM_DOMAIN_WIDTH_X/double(100.00);
			int num_of_div_y=PRM.SIM_DOMAIN_WIDTH_Y/double(100.00);
			int num_of_div_z=PRM.SIM_DOMAIN_WIDTH_Z/double(100.00);
			for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){
				int ix=(1e12*IONS[INDEX_STAT_STEP][i].x-PRM.MIN_X)/double(100.00);
				int iy=(1e12*IONS[INDEX_STAT_STEP][i].y-PRM.MIN_Y)/double(100.00);
				int iz=(1e12*IONS[INDEX_STAT_STEP][i].z-PRM.MIN_Z)/double(100.00);
				if((ix>=0 && ix<num_of_div_x) && (iy>=0 && iy<num_of_div_y) && (iz>=0 && iz<num_of_div_z)){
					concs_3D.at(IONS[INDEX_STAT_STEP][i].kind).at(ix).at(iy).at(iz)+=1.00;
				}
			}	
		}
		else{ // rotational symmetry
			int num_of_div_on_z_stat=PRM.SIM_DOMAIN_WIDTH_Z/double(100.00);
			int num_of_div_on_r_stat=PRM.SIM_DOMAIN_WIDTH_X/double(100.00);
			for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){
				double dx=(1e12*IONS[INDEX_STAT_STEP][i].x);
				double dy=(1e12*IONS[INDEX_STAT_STEP][i].y);
				double dist=sqrt(dx*dx+dy*dy);
				int ir=(dist)/double(100.00);
				int iz=(1e12*IONS[INDEX_STAT_STEP][i].z-PRM.MIN_Z)/double(100.00);
				if((iz>=0 && iz<num_of_div_on_z_stat) && (ir>=0 && ir<num_of_div_on_r_stat)){
					radial_concs.at(IONS[INDEX_STAT_STEP][i].kind).at(iz).at(ir)+=1.00;
				}
			}
		}
	}
	
	
// potential computation	
	if(PRM.potential == 1){	
		for(int i=0; i<2001; i++){
			if(isnormal(POTENTIALS_ON_AXIS[INDEX_STAT_STEP][i]) || POTENTIALS_ON_AXIS[INDEX_STAT_STEP][i]==0){
				AVERAGE_POTENTIALS_ON_AXIS[i]+=(POTENTIALS_ON_AXIS[INDEX_STAT_STEP][i]-POTENTIALS_ON_AXIS[INDEX_STAT_STEP][2000]);
			}
		}
	}	
	
// currents computation	(ZERO THRESHOLD, one threshold at z=0)
	if(PRM.currents_ZT){	
		compute_currents_ZT();
	}
	
	if(PRM.channel_configuration){
		int counter_Na=0;	// 3
		int counter_K=0;	// 4
		int counter_Mg=0;	// 8
		int counter_Ca=0;	// 9
		int counter_Cl=0;	// 17
		
		for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){
			if(1e12*IONS[INDEX_STAT_STEP][i].z>PRM.Z_MOUTH_LEFT && 1e12*IONS[INDEX_STAT_STEP][i].z<PRM.Z_MOUTH_RIGHT){
				if(IONS[INDEX_STAT_STEP][i].kind==3){
					counter_Na++;
				}
				if(IONS[INDEX_STAT_STEP][i].kind==4){
					counter_K++;
				}
				if(IONS[INDEX_STAT_STEP][i].kind==8){
					counter_Mg++;
				}
				if(IONS[INDEX_STAT_STEP][i].kind==9){
					counter_Ca++;
				}
				if(IONS[INDEX_STAT_STEP][i].kind==17){
					counter_Cl++;
				}
			}
		}
	
		CHANNEL_CONFIGURATIONS[counter_Na][counter_K][counter_Mg][counter_Ca][counter_Cl]=CHANNEL_CONFIGURATIONS[counter_Na][counter_K][counter_Mg][counter_Ca][counter_Cl]+1.00;
		
		if(last_output_channel_configuration[0]==counter_Na && last_output_channel_configuration[1]==counter_K && last_output_channel_configuration[2]==counter_Mg && last_output_channel_configuration[3]==counter_Ca &&  last_output_channel_configuration[4]==counter_Cl){
			channel_same_configuration_counter++;
		}
		else{
		
			string s_counter_Na;
			stringstream ss_counter_Na;
			ss_counter_Na << counter_Na;
			ss_counter_Na >> s_counter_Na;
			
			string s_counter_K;
			stringstream ss_counter_K;
			ss_counter_K << counter_K;
			ss_counter_K >> s_counter_K;
			
			string s_counter_Mg;
			stringstream ss_counter_Mg;
			ss_counter_Mg << counter_Mg;
			ss_counter_Mg >> s_counter_Mg;
			
			string s_counter_Ca;
			stringstream ss_counter_Ca;
			ss_counter_Ca << counter_Ca;
			ss_counter_Ca >> s_counter_Ca;
			
			string s_counter_Cl;
			stringstream ss_counter_Cl;
			ss_counter_Cl << counter_Cl;
			ss_counter_Cl >> s_counter_Cl;
			
			stat_file_A=PRM.PREFIX + "_"+ s_counter_Na+ s_counter_K+ s_counter_Mg+ s_counter_Ca+ s_counter_Cl +".channel.freq";
			if(!fileExists(stat_file_A)){
				createFile(stat_file_A);
			}
			
			char *tfnA = new char[stat_file_A.length()+1];
			strcpy(tfnA, stat_file_A.c_str());     
			ofstream foutA(tfnA, ios::app);
			
			foutA << channel_same_configuration_counter << endl;
		
			foutA.close();
			
			last_output_channel_configuration[0]=counter_Na;
			last_output_channel_configuration[1]=counter_K;
			last_output_channel_configuration[2]=counter_Mg;
			last_output_channel_configuration[3]=counter_Ca;
			last_output_channel_configuration[4]=counter_Cl;
			
			channel_same_configuration_counter=0;
		}
			
		
	}		
	if(PRM.filter_configuration){
		int counter_Na=0;	// 3
		int counter_K=0;	// 4
		int counter_Mg=0;	// 8
		int counter_Ca=0;	// 9
		int counter_Cl=0;	// 17
		
		for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){
			if(1e12*IONS[INDEX_STAT_STEP][i].z>-500 && 1e12*IONS[INDEX_STAT_STEP][i].z<500){
				if(IONS[INDEX_STAT_STEP][i].kind==3){
					counter_Na++;
				}
				if(IONS[INDEX_STAT_STEP][i].kind==4){
					counter_K++;
				}
				if(IONS[INDEX_STAT_STEP][i].kind==8){
					counter_Mg++;
				}
				if(IONS[INDEX_STAT_STEP][i].kind==9){
					counter_Ca++;
				}
				if(IONS[INDEX_STAT_STEP][i].kind==17){
					counter_Cl++;
				}
			}
		}
	
		FILTER_CONFIGURATIONS[counter_Na][counter_K][counter_Mg][counter_Ca][counter_Cl]=FILTER_CONFIGURATIONS[counter_Na][counter_K][counter_Mg][counter_Ca][counter_Cl]+1.00;
		
		
		
		if(last_output_filter_configuration[0]==counter_Na && last_output_filter_configuration[1]==counter_K && last_output_filter_configuration[2]==counter_Mg && last_output_filter_configuration[3]==counter_Ca &&  last_output_filter_configuration[4]==counter_Cl){
			filter_same_configuration_counter++;
		}
		else{
		
			string s_counter_Na;
			stringstream ss_counter_Na;
			ss_counter_Na << counter_Na;
			ss_counter_Na >> s_counter_Na;
			
			string s_counter_K;
			stringstream ss_counter_K;
			ss_counter_K << counter_K;
			ss_counter_K >> s_counter_K;
			
			string s_counter_Mg;
			stringstream ss_counter_Mg;
			ss_counter_Mg << counter_Mg;
			ss_counter_Mg >> s_counter_Mg;
			
			string s_counter_Ca;
			stringstream ss_counter_Ca;
			ss_counter_Ca << counter_Ca;
			ss_counter_Ca >> s_counter_Ca;
			
			string s_counter_Cl;
			stringstream ss_counter_Cl;
			ss_counter_Cl << counter_Cl;
			ss_counter_Cl >> s_counter_Cl;
			
			stat_file_D=PRM.PREFIX + "_"+ s_counter_Na+ s_counter_K+ s_counter_Mg+ s_counter_Ca+ s_counter_Cl +".filter.freq";
			if(!fileExists(stat_file_D)){
				createFile(stat_file_D);
			}
			
			char *tfnD = new char[stat_file_D.length()+1];
			strcpy(tfnD, stat_file_D.c_str());     
			ofstream foutD(tfnD, ios::app);
			
			foutD << filter_same_configuration_counter << endl;
		
			foutD.close();
			
			last_output_filter_configuration[0]=counter_Na;
			last_output_filter_configuration[1]=counter_K;
			last_output_filter_configuration[2]=counter_Mg;
			last_output_filter_configuration[3]=counter_Ca;
			last_output_filter_configuration[4]=counter_Cl;
			
			filter_same_configuration_counter=0;
		}
	}
	return;
}

void Statistics::compute_currents_RS(){
	
	vector <double> this_step_currents_RS;
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		this_step_currents_RS.push_back(0.00);
	}
	
	for(int ion_index_1=0; ion_index_1<ions.front().size(); ion_index_1++){
		if(1e12*ions.front().at(ion_index_1).z>PRM.LEFT_CELL_MAX_Z && 1e12*ions.front().at(ion_index_1).z<PRM.RIGHT_CELL_MIN_Z){

			double velocity=get_ion_velocity(ions.front().at(ion_index_1));
			if(velocity>=0.00 && velocity<1.5e4){
				
			
	// external field
				double Ex=0.00;
				double Ey=0.00;
				double Ez=1e12/PRM.SIM_DOMAIN_WIDTH_Z;
								
				
	// boundary repulsion and induced charges
				if(!surfaces.empty() && PRM.EPS_W!=PRM.EPS_MEM){
					for(int surface_index=0; surface_index<surfaces.size(); surface_index++){		
						double charge_x=surfaces.at(surface_index).center.x;
						double charge_y=surfaces.at(surface_index).center.y;
						double charge_z=surfaces.at(surface_index).center.z;				
						double distance=1e12*get_distance(ions.front().at(ion_index_1).x, ions.front().at(ion_index_1).y, ions.front().at(ion_index_1).z, charge_x, charge_y, charge_z);	

						if(distance<20000){
							double field=0.00;
							double dx=0.00;
							double dy=0.00;
							double dz=0.00;
							
							dx=1e12*(ions.front().at(ion_index_1).x-charge_x)/distance;
							dy=1e12*(ions.front().at(ion_index_1).y-charge_y)/distance;
							dz=1e12*(ions.front().at(ion_index_1).z-charge_z)/distance;
							
							double surf_charge_valence=(vector_h_RAMO[surface_index]*surfaces.at(surface_index).area)/Q;
							field=surf_charge_valence*FIELD_C[int(distance)];
							
							Ex+=field*dx;
							Ey+=field*dy;
							Ez+=field*dz;
						}	
						else{
							// cut-off
						}			
					}
				}

			
	// compute current
				double ion_current_RS=ions.front().at(ion_index_1).charge*(ions.front().at(ion_index_1).velocity[0]*Ex+ions.front().at(ion_index_1).velocity[1]*Ey+ions.front().at(ion_index_1).velocity[2]*Ez);
				
				this_step_currents_RS.at(ions.front().at(ion_index_1).kind)+=ion_current_RS;
			}
		}
	}
	
	for(int i=0; i<NUM_OF_IONIC_SPECIES; i++){
		currents_RS.at(i)+=this_step_currents_RS.at(i);
		instant_currents_RS.at(i)=this_step_currents_RS.at(i);
	}
	
	return;
}

void Statistics::compute_currents_ZT(){

	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		instant_currents_ZT.at(is)=0;
	}
	
	for(int ion_index_1=0; ion_index_1<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; ion_index_1++){
		if(1e12*IONS[INDEX_STAT_STEP][ion_index_1].z_prev>PRM.LEFT_CELL_MAX_Z && 1e12*IONS[INDEX_STAT_STEP][ion_index_1].z_prev<PRM.RIGHT_CELL_MIN_Z && 		//to avoid counting ions flowing through the external circuit...
			1e12*IONS[INDEX_STAT_STEP][ion_index_1].z>PRM.LEFT_CELL_MAX_Z && 1e12*IONS[INDEX_STAT_STEP][ion_index_1].z<PRM.RIGHT_CELL_MIN_Z){
		
			//CHECK LEFT TO RIGHT
			if(IONS[INDEX_STAT_STEP][ion_index_1].z_prev<0 && IONS[INDEX_STAT_STEP][ion_index_1].z>=0){
				currents_ZT.at(IONS[INDEX_STAT_STEP][ion_index_1].kind)=currents_ZT.at(IONS[INDEX_STAT_STEP][ion_index_1].kind)+IONS[INDEX_STAT_STEP][ion_index_1].charge/PRM.DELTA_T;
				instant_currents_ZT.at(IONS[INDEX_STAT_STEP][ion_index_1].kind)+=IONS[INDEX_STAT_STEP][ion_index_1].charge/PRM.DELTA_T;
			}
			else if(IONS[INDEX_STAT_STEP][ion_index_1].z_prev>0 && IONS[INDEX_STAT_STEP][ion_index_1].z<=0){
				currents_ZT.at(IONS[INDEX_STAT_STEP][ion_index_1].kind)=currents_ZT.at(IONS[INDEX_STAT_STEP][ion_index_1].kind)-IONS[INDEX_STAT_STEP][ion_index_1].charge/PRM.DELTA_T;
				instant_currents_ZT.at(IONS[INDEX_STAT_STEP][ion_index_1].kind)+=(-IONS[INDEX_STAT_STEP][ion_index_1].charge/PRM.DELTA_T);
			}
			else{
				
			}
		}
	}
	
	return;
}

void Statistics::print_statistics(){
	
	//~ double time_computed=double(step_window.at(0))*PRM.DELTA_T;
	double time_computed=double(STEPS[INDEX_STAT_STEP])*PRM.DELTA_T;
	
// ions position
	if(PRM.channel_pdb_files){
		string trajectory_file=PRM.PREFIX + ".pdb";
		char *tfn1 = new char[trajectory_file.length()+1];
		strcpy(tfn1, trajectory_file.c_str());
		ofstream fout(tfn1, ios::out);
		for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_LAST_STEP]; i++){
			fout<<"ATOM  "<<setw(5)<<i+1<<" "
				<<setw(4)<<IONS[INDEX_LAST_STEP][i].name<<" "
				<<setw(3)<<IONS[INDEX_LAST_STEP][i].name<<" "
				<<setw(5)<<i+1<<"    "			
				<<setw(8)<<setprecision(3)<<IONS[INDEX_LAST_STEP][i].x*1e10 //atom X coordinate
				<<setw(8)<<setprecision(3)<<IONS[INDEX_LAST_STEP][i].y*1e10 //atom Y coordinate
				<<setw(8)<<setprecision(3)<<IONS[INDEX_LAST_STEP][i].z*1e10 //atom Z coordinate
				<<setw(6)<<"0.0"		//occupancy
				<<setw(6)<<"0.0"		//tempFactor
				<<endl;
		}
		fout.close();
	}
	
// flux
	if(PRM.flux){
		createFile(stat_file_2);
		char *tfn2 = new char[stat_file_2.length()+1];
		strcpy(tfn2, stat_file_2.c_str());     
		ofstream fout(tfn2, ios::app);
		float z_pm,section_area_pm2, volume_slice_A3, ionsA3_2_mol, concentration_ions_slice;
		ionsA3_2_mol = 1e27/AVOGADRO;
		for(int i=0; i<num_of_dz; i++){
			z_pm = (PRM.MIN_Z+double(i)*DELTA_Z+double(0.50)*DELTA_Z);
			if(PRM.SIM_TYPE.compare("BULK") == 0){
				section_area_pm2 = PRM.SIM_DOMAIN_WIDTH_X*PRM.SIM_DOMAIN_WIDTH_Y;
			} else {
				if (((z_pm+double(0.50)*DELTA_Z) < -0.5*PRM.MEMBRANE_WIDTH) || ((z_pm-double(0.50)*DELTA_Z) > 0.5*PRM.MEMBRANE_WIDTH)) {
					section_area_pm2 = PRM.SIM_DOMAIN_WIDTH_X*PRM.SIM_DOMAIN_WIDTH_Y;
				} else {
					int ind_in_limits = round(z_pm + 0.5*PRM.MEMBRANE_WIDTH) - 1;
					if (ind_in_limits < 0) {
						ind_in_limits = 0;
					} else if (ind_in_limits >= limits.size()) {
						ind_in_limits =  limits.size() - 1;
					}
					section_area_pm2 = M_PI*limits.at(ind_in_limits)*limits.at(ind_in_limits);
				}
			}
			fout << z_pm*1e-2 << "\t" << section_area_pm2*1e-4 << "\t"; // output is in A
			volume_slice_A3 = DELTA_Z*section_area_pm2*1e-6;
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				if(PRM.ions_to_simulate.at(is)){
					concentration_ions_slice = concentrations_along_z.at(is).at(i)/double(STEPS[INDEX_STAT_STEP]);
					fout << currents_ZT.at(is)/double(STEPS[INDEX_STAT_STEP]+1)/sample_ions[is].charge << " "
					     << concentration_ions_slice << " " 
					     << ionsA3_2_mol*concentration_ions_slice/volume_slice_A3 << " " 
					     << (currents_ZT.at(is)/double(STEPS[INDEX_STAT_STEP]+1)/sample_ions[is].charge)/concentrations_along_z.at(is).at(i)/(1e-2*DELTA_Z*double(STEPS[INDEX_STAT_STEP])) << " " ;
				}
			}
			fout<<endl;
		}
		fout.close();
	}

	
// radial distribution function computation	
	if(PRM.rdf){
		createFile(stat_file_3);
		char *tfn3 = new char[stat_file_3.length()+1];
		strcpy(tfn3, stat_file_3.c_str());     
		ofstream fout3(tfn3, ios::app);
		
		//~ double total_volume=1.33333e-36*M_PI*pow(double(1000.00), double(3.00));
		//~ double total_volume=1.33333e-36*PRM.SIM_DOMAIN_WIDTH_X*PRM.SIM_DOMAIN_WIDTH_Y*PRM.SIM_DOMAIN_WIDTH_Z;
		
		for(int i=0; i<1000; i++){
			fout3 << double(i)/double(100.00)+double(0.005) << " ";
			double volume=1.33333e-6*M_PI*(pow(double(i+1), double(3.00))-pow(double(i), double(3.00)));	//volume is in Angstroms^3
			for(int is1=0; is1<NUM_OF_IONIC_SPECIES; is1++){
				if(PRM.ions_to_simulate.at(is1)){
					for(int is2=is1; is2<NUM_OF_IONIC_SPECIES; is2++){
						if(PRM.ions_to_simulate.at(is2)){
							//~ fout3 << (RDF.at(is1).at(is2).at(i)*total_volume)/(volume*RDF_samples.at(is1))<<" ";
							//~ fout3 << (RDF.at(is1).at(is2).at(i)*total_volume)/(volume*RDF_samples.at(is1))<<" ";
							
							fout3 << (RDF.at(is1).at(is2).at(i))/(volume*RDF_samples.at(is1))<<" ";
							//~ fout3 << RDF.at(is1).at(is2).at(i)<<" ";
							
							
							
							
							//~ fout3 << (RDF.at(is1).at(is2).at(i))<<" ";
							
							//~ fout3 << (RDF.at(is1).at(is2).at(i))/(AVOGADRO*volume*concs.at(is2)*RDF_samples.at(is1).at(is2))<<" ";
							//~ fout3 << RDF.at(is1).at(is2).at(i)/(volume)<<" ";
							//~ cout  << "RDF_samples.at(is1).at(is2): " <<RDF_samples.at(is1).at(is2)<<endl; 
							//~ cout  << "RDF.at(is1).at(is2).at(i): " <<RDF.at(is1).at(is2).at(i)<<endl; 
							//~ cout  << "total_volume: " <<total_volume<<endl; 
							//~ cout  << "volume: " <<volume<<endl; 
							//~ sleep(1);
						}
					}
				}
			}
			fout3<<" " << volume << endl;
		}
		
		fout3.close();
	}				

// velocity distribution computation	
	if(PRM.vel_distribution){
		
		createFile(stat_file_4);
		char *tfn4 = new char[stat_file_4.length()+1];
		strcpy(tfn4, stat_file_4.c_str());     
		ofstream fout4(tfn4, ios::app);
		
		for(int iv=0; iv<150; iv++){
			fout4 << iv*10 << " ";
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				if(PRM.ions_to_simulate.at(is)){
					Ion ion_curr(is);
					fout4 << velocities_theory.at(is).at(iv) << " " << velocities.at(is).at(iv)/(velocities_samples*double(10.00)) << " ";
				}
			}
			fout4<<endl;
		}
		
		fout4.close();
	}
	

// mean square displacement computation	
	if(PRM.mean_square_displ){
		if(step_window.at(0)<=0){
			if(step_window.at(0)==0){
				for(int i=0; i<ions.front().size(); i++){
					msds_num.at(ions.front().at(i).kind)=msds_num.at(ions.front().at(i).kind)+double(1.00);
				}
			}
			
			msds_index_1=1;
			for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
				if(PRM.ions_to_simulate.at(is)){
					for(int iv=0; iv<PRM.STATS_OUT_FREQ; iv++){
						msds.at(is).at(iv)=0.00;
					}
				}
			}
		}
		else{
			createFile(stat_file_B);
			char *tfnB = new char[stat_file_B.length()+1];
			strcpy(tfnB, stat_file_B.c_str());     
			ofstream foutB(tfnB, ios::app);
			
			for(int iv=0; iv<PRM.STATS_OUT_FREQ; iv++){
				foutB << 1e12*iv*PRM.DELTA_T << " ";
				for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
					if(PRM.ions_to_simulate.at(is)){
						Ion ion_curr(is);
						foutB << 1e20*double(6.00)*ion_curr.diffusion_coeff*iv*PRM.DELTA_T << " " << (msds.at(is).at(iv))/(double(msds_index_1)*msds_num.at(is)) << " ";
					}
				}
				foutB<<endl;
			}
			
			foutB.close();
		}
		
		msds_start.clear();
		Point p;
		for(int i=0; i<ions.front().size(); i++){
			p.x=ions.front().at(i).x;
			p.y=ions.front().at(i).y;
			p.z=ions.front().at(i).z;
			
			//~ cout << endl<<ions.front().at(i).kind <<"\t"<<1e12*p.x<<"\t"<<1e12*p.y<<"\t"<<1e12*p.z;
			msds_start.push_back(p);
		}
		//~ cout <<endl<<endl;
		//~ sleep(3);
		
		msds_index_1++;
		msds_index_2=1;
	}	

	
	
	if(PRM.induced_charge){
		char *tfn6 = new char[stat_file_6.length()+1];
		strcpy(tfn6, stat_file_6.c_str());     
		ofstream fout6(tfn6, ios::app);
		
		for(int i=0; i<C_total_charge.size(); i++){
			fout6 << A_total_charge << " " << C_total_charge.at(i)  << " " << C_total_charge_type.at(i) <<endl;
		}
		fout6.close();
		C_total_charge.clear();
		C_total_charge_type.clear();
	}

// potential computation	
	if(PRM.potential == 1){	
		createFile(stat_file_C);
		char *tfnC = new char[stat_file_C.length()+1];
		strcpy(tfnC, stat_file_C.c_str());     
		ofstream foutC(tfnC, ios::app);
		for(int i=0; i<2001; i++){
			// convert from microVolts for accumulation in computing the average
			foutC << 1e10*POINTS_ON_AXIS[i] << "\t" << 1e9*AVERAGE_POTENTIALS_ON_AXIS[i]/double(STEPS[INDEX_STAT_STEP]+1) << endl;
		}
		foutC.close();
	}
	
// concentrations computation	
	if(PRM.concentrations!=0){	
		createFile(stat_file_7);
		char *tfn7 = new char[stat_file_7.length()+1];
		strcpy(tfn7, stat_file_7.c_str());     
		ofstream fout7(tfn7, ios::app);
			
		if(PRM.concentrations==3){ //average on y (?)
			int num_of_div_x=PRM.SIM_DOMAIN_WIDTH_X/double(100.00);
			int num_of_div_y=PRM.SIM_DOMAIN_WIDTH_Y/double(100.00);
			int num_of_div_z=PRM.SIM_DOMAIN_WIDTH_Z/double(100.00);
			
			double cell_volume_for_M=1e-33*AVOGADRO*double(100.00)*double(100.00)*double(100.00);
			
			for(int iz=0; iz<num_of_div_z; iz++){
			
				for(int ix=0; ix<num_of_div_x; ix++){
					
					for(int iy=0; iy<num_of_div_y; iy++){
					
						fout7 	<< ((PRM.MIN_Z)+double(50.0)+double(iz)*double(100.00))/double(100.00) << " "
						<< ((PRM.MIN_X)+double(50.0)+double(ix)*double(100.00))/double(100.00) << " "
						<< ((PRM.MIN_Y)+double(50.0)+double(iy)*double(100.00))/double(100.00) << " ";
						
						for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
							if(PRM.ions_to_simulate.at(is)){
								fout7 	<< concs_3D.at(is).at(ix).at(iy).at(iz)/(cell_volume_for_M*double(step_window.at(0))) << " ";
							}
						}
						fout7 << endl;
					}
					fout7<<endl;
				}
				fout7<<endl<<endl<<endl;
			}
		} else{ // rotational
			int num_of_div_on_z_stat=PRM.SIM_DOMAIN_WIDTH_Z/double(100.00);
			int num_of_div_on_r_stat=PRM.SIM_DOMAIN_WIDTH_X/double(100.00);
			
			for(int iz=0; iz<num_of_div_on_z_stat; iz++){
			
				for(int ir=0; ir<num_of_div_on_r_stat; ir++){
					double max_rad=double(ir+1)*double(100.00);
					double max_vol=M_PI*max_rad*max_rad;
					double min_rad=double(ir)*double(100.00);
					double min_vol=M_PI*min_rad*min_rad;
					double cell_volume_for_M=1e-33*AVOGADRO*double(100.00)*(max_vol-min_vol);
					
					vector <double> this_concs;
					for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
						this_concs.push_back(0.00);
					}
					
					for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
						if(PRM.ions_to_simulate.at(is)){
							this_concs.at(is)=radial_concs.at(is).at(iz).at(ir)/(cell_volume_for_M*double(step_window.at(0)));
						}
					}
					
					fout7 	<< ((PRM.MIN_Z)+double(50.0)+double(iz)*double(100.00))/double(100.00) << " "
						<< (double(50.0)+double(ir)*double(100.00))/double(100.00) << " ";
					
					for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
						if(PRM.ions_to_simulate.at(is)){
							fout7 	<< this_concs.at(is) << " ";
						}
					}
					
					
					fout7 << endl;
					this_concs.clear();
				}
				
				fout7<<endl;
			}
		}
		fout7.close();
	}
	
// currents computation	(ZERO THRESHOLD, one threshold at z=0)
	if(PRM.currents_ZT){	
		
		char *tfn9 = new char[stat_file_9.length()+1];
		strcpy(tfn9, stat_file_9.c_str());     
		ofstream fout9(tfn9, ios::app);
		
		fout9 << time_computed << " ";
			
		double total_current=0.00;
		double instant_total_current=0.00;
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			if(PRM.ions_to_simulate.at(is)){
				fout9 << instant_currents_ZT.at(is) << " " << currents_ZT.at(is)/double(STEPS[INDEX_STAT_STEP]+1) << " ";
				total_current+=currents_ZT.at(is)/double(STEPS[INDEX_STAT_STEP]+1);
				instant_total_current+=instant_currents_ZT.at(is);
			}
		}
			
		fout9 << instant_total_current << " " << total_current <<endl;
		
		fout9.close();
	}	
	
	if(PRM.channel_configuration){
		
		createFile(stat_file_6);
		char *tfn6 = new char[stat_file_6.length()+1];
		strcpy(tfn6, stat_file_6.c_str());     
		ofstream fout6(tfn6, ios::app);
		
		for(int a=0; a<3; a++){
			for(int b=0; b<3; b++){
				for(int c=0; c<3; c++){
					for(int d=0; d<3; d++){
						for(int e=0; e<3; e++){
							fout6 << CHANNEL_CONFIGURATIONS[a][b][c][d][e]/double(STEPS[INDEX_STAT_STEP]) << "\n";
						}
					}
				}
			}
		}
		
		fout6.close();
		
		
		//~ createFile(stat_file_A);
		//~ char *tfnA = new char[stat_file_A.length()+1];
		//~ strcpy(tfnA, stat_file_A.c_str());     
		//~ ofstream foutA(tfnA, ios::app);
		
		//~ for(int f=0; f<100; f++){
			//~ for(int a=0; a<3; a++){
				//~ for(int b=0; b<3; b++){
					//~ for(int c=0; c<3; c++){
						//~ for(int d=0; d<3; d++){
							//~ for(int e=0; e<3; e++){
								//~ foutA << CHANNEL_CONFIGURATIONS_FREQUENCIES[a][b][c][d][e][f] << "\t";
							//~ }
						//~ }
					//~ }
				//~ }
			//~ }
			//~ foutA <<endl;
		//~ }
		
		//~ foutA.close();
		
	}
		
			
	if(PRM.filter_configuration){
		createFile(stat_file_8);
		char *tfn8 = new char[stat_file_8.length()+1];
		strcpy(tfn8, stat_file_8.c_str());     
		ofstream fout8(tfn8, ios::app);
		
		for(int a=0; a<3; a++){
			for(int b=0; b<3; b++){
				for(int c=0; c<3; c++){
					for(int d=0; d<3; d++){
						for(int e=0; e<3; e++){
							fout8 << FILTER_CONFIGURATIONS[a][b][c][d][e]/double(STEPS[INDEX_STAT_STEP]) << "\n";
						}
					}
				}
			}
		}
		
		fout8.close();
		
		//~ createFile(stat_file_D);
		//~ char *tfnD = new char[stat_file_D.length()+1];
		//~ strcpy(tfnD, stat_file_D.c_str());     
		//~ ofstream foutD(tfnD, ios::app);
		
		//~ for(int f=0; f<100; f++){
			//~ for(int a=0; a<3; a++){
				//~ for(int b=0; b<3; b++){
					//~ for(int c=0; c<3; c++){
						//~ for(int d=0; d<3; d++){
							//~ for(int e=0; e<3; e++){
								//~ foutD << FILTER_CONFIGURATIONS_FREQUENCIES[a][b][c][d][e][f] << "\t";
							//~ }
						//~ }
					//~ }
				//~ }
			//~ }
			//~ foutD <<endl;
		//~ }
		
		//~ foutD.close();
	}
	
	return;
}

void Statistics::reset_after_restarting_from_the_beginning(){
	
	if(PRM.potential == 1){	
		for(int i=0; i<50; i++){
			for(int j=0; j<2001; j++){
				POTENTIALS_ON_AXIS[i][j]=0.00;
				AVERAGE_POTENTIALS_ON_AXIS[j]=0.00;
			}
		}
	}
	if(PRM.flux){
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			for(int j=0; j<concentrations_along_z[is].size(); j++){
				concentrations_along_z[is][j]=0.00;
			}
		}
	}
	if(PRM.currents_ZT){	
		for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
			currents_ZT.at(is)=0.00;
			instant_currents_ZT.at(is)=0.00;
		}
	}
	
	
	if(PRM.channel_configuration){
		for(int a=0; a<5; a++){
			for(int b=0; b<5; b++){
				for(int c=0; c<5; c++){
					for(int d=0; d<5; d++){
						for(int e=0; e<5; e++){
							CHANNEL_CONFIGURATIONS[a][b][c][d][e]=0;
						}
					}
				}
			}
		}
	}
	
	if(PRM.filter_configuration){
		for(int a=0; a<5; a++){
			for(int b=0; b<5; b++){
				for(int c=0; c<5; c++){
					for(int d=0; d<5; d++){
						for(int e=0; e<5; e++){
							FILTER_CONFIGURATIONS[a][b][c][d][e]=0;
						}
					}
				}
			}
		}
	}
	
	return;
}

ostream& operator<<(ostream& stream, Statistics& STAT){
	
	stream << endl << "########################################" <<endl;
	stream << "Statistics" <<endl;
	
	return stream;
}

Control_cell::Control_cell(){
	
}

Control_cell::~Control_cell(){
	
}

void Control_cell::reset_control_cell(){
	
	SIDE=-1;
	
	eta.clear();
	SumEta.clear();
	SumNi.clear();
	Xi.clear();
	
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){
		eta.push_back(0.00);
		SumEta.push_back(0.00);
		SumNi.push_back(0.00);
		Xi.push_back(0.00);
		in_control_cell.push_back(false);
	}

	return;
}

ostream& operator<<(ostream& stream, Control_cell& cell){
	
	stream << endl << "########################################" <<endl;
	stream << "Control_cell" <<endl;
	
	return stream;
}

//########################################
// Class Point
//########################################
Point::Point(){
	reset_point();
}

Point::~Point(){
	
}

void Point::reset_point(){
	
	x=0.00;
	y=0.00;
	z=0.00;

	potential=0.00;
	energy=0.00;
	
	field[0]=0.00;
	field[1]=0.00;
	field[2]=0.00;
	
}

ostream& operator<<(ostream& stream, Point& p){
	
	stream << endl << "########################################" <<endl;
	stream << "Point" <<endl;
	stream << "position:\t"<<p.x << "\t" << p.y << "\t" << p.z << endl;
	stream << "pot/energy:\t"<<p.potential << "\t" << p.energy << endl;
	stream << "field:\t\t"<<p.field[0] << "\t" << p.field[1] << "\t" << p.field[2] << endl;
	
	return stream;
}

//########################################
// Class RadialPoint
//########################################
RadialPoint::RadialPoint(){
	reset_radialpoint();
}

RadialPoint::~RadialPoint(){
	
}

void RadialPoint::reset_radialpoint(){
	
	z=0.00;
	r=0.00;

	potential=0.00;
	energy=0.00;
	
	field[0]=0.00;
	field[1]=0.00;
	
}

ostream& operator<<(ostream& stream, RadialPoint& rp){
	
	stream << endl << "########################################" <<endl;
	stream << "RadialPoint" <<endl;
	stream << "position:\t"<<rp.z << "\t" << rp.r << endl;
	stream << "pot/energy:\t"<<rp.potential << "\t" << rp.energy << endl;
	stream << "field:\t\t"<<rp.field[0] << "\t" << rp.field[1] << endl;
	
	return stream;
}

//########################################
// Class Charge
//########################################
Charge::Charge(){
	reset_charge();
}

Charge::~Charge(){
	
}

void Charge::reset_charge(){
	
	x=0.00;
	y=0.00;
	z=0.00;
	
	charge=0.00;
	radius=0.00;
	
}

ostream& operator<<(ostream& stream, Charge& c){

	stream << endl << "########################################" <<endl;
	stream << "Charge" <<endl;
	stream << "charge:\t"<<c.charge<<endl;
	stream << "radius:\t"<<c.charge<<endl;
	stream << "position:\t"<<c.x << "\t" << c.y << "\t" << c.z << endl;

	return stream;
}

//########################################
// Class Ion
//########################################
Ion::Ion(){}

Ion::Ion(int type){
	

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
	
// 
	20		P00		
	21		P01		
	22		P02		
	23		P03		
	24		P04		
	25		P05		
	26		P06		
	27		P07		
	28		P08		
	29		P09		
	30		P10		
	31		P11		
	32		P12		
	33		P13		
	34		P14		
	35		P15		
	36		P16		
	37		P17		
	38		P18		
	39		P19		
	40		N00		
	41		N01		
	42		N02		
	43		N03		
	44		N04		
	45		N05		
	46		N06		
	47		N07		
	48		N08		
	49		N09		
	50		N10		
	51		N11		
	52		N12		
	53		N13		
	54		N14		
	55		N15		
	56		N16		
	57		N17		
	58		N18		
	59		N19		

	60		Q10		
	61		Q11		
	
	62		Q12		
	63		Q13		
	64		Q14		
	65		Q15		
	66		Q16		
	67		Q17		

	68		Q18		
	69		Q19		
	70		Q20		
	71		Q21		
	72		Q22		
	73		Q23		

	74		Q24		
	75		Q25		
	76		Q26		
	77		Q27		
	78		Q28		
	79		Q29		
	
	80		X00		
	81		X01		
	82		X02		
	83		X03		
	84		X04		
	85		X05		
	
	86		X06		
	87		X07		
	88		X08		
	89		X09		
	90		X10		
	91		X11		

	92		MEMBRANE
	
	93		FIXED CHARGE
	
														*/
//===========================================	

	reset_ion();
	
	kind=type;
	
	
	if(type==0){
		name="O";
		charge=-double(0.50)*Q;
		valence=-0.50;
		mass=26.552e-27;
		radius=1.40e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; //[J]
		diffusion_coeff=2.2e-9;
	}
	else if(type==1){
		name="H";
		charge=Q;
		valence=1.00;
		mass=1.67262158e-27;
		radius=0.8775e-15;
		charmm_half_radius = 0.2245e-10; // [m]
		charmm_eps = -0.046 * 4186.8 / AVOGADRO; //[J]
		diffusion_coeff=9e-9;
	}
	else if(type==2){
		name="LI";
		charge=Q;
		valence=1.00;
		mass=0.00;
		radius=0.60e-10;
		diffusion_coeff=0.00;
	}
	else if(type==3){
		name="NA";
		charge=Q;
		valence=1.00;
		mass=38.4532e-27;
		radius=0.95e-10;
		charmm_half_radius = 1.36375e-10; // [m]
		charmm_eps = -0.0469 * 4186.8 / AVOGADRO; //[J]
		diffusion_coeff=1.33e-9;
//-------------------------------------------------------------- new GUI
		valence=NA_VALENCE;
		charge=NA_VALENCE*Q;
		radius=1e-10*NA_RADIUS;
		diffusion_coeff=1e-9*NA_DIFF_COEFF;
		mass=1e-27*NA_MASS;	
//-------------------------------------------------------------- new GUI
	}
	else if(type==4){
		name="K";
		charge=Q*K_VALENCE;
		valence=K_VALENCE;
		mass=65.3967e-27;
		radius=1.33e-10;
		charmm_half_radius = 1.76375e-10; // [m]
		charmm_eps = -0.0870 * 4186.8 / AVOGADRO; //[J]
		diffusion_coeff=1.96e-9;
//-------------------------------------------------------------- new GUI
		valence=K_VALENCE;
		charge=K_VALENCE*Q;
		radius=1e-10*K_RADIUS;
		diffusion_coeff=1e-9*K_DIFF_COEFF;
		mass=1e-27*K_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==5){
		name="RB";
		charge=Q;
		valence=1.00;
		mass=0.00;
		radius=1.48e-10;
		diffusion_coeff=0.00;
	}
	else if(type==6){
		name="CS";
		charge=Q;
		valence=1.00;
		mass=0.00;
		radius=1.69e-10;
		diffusion_coeff=0.00;
	}
	else if(type==7){
		name="TL";
		charge=Q;
		valence=1.00;
		mass=0.00;
		radius=1.40e-10;
		diffusion_coeff=0.00;
	}
	else if(type==8){
		name="MG";
		charge=double(2.00)*Q;
		valence=2.00;
		mass=0.00;
		radius=0.65e-10;
		charmm_half_radius = 1.185e-10; // [m]
		charmm_eps = -0.0150 * 4186.8 / AVOGADRO; //[J]
		diffusion_coeff=0.00;
//-------------------------------------------------------------- new GUI
		valence=MG_VALENCE;
		charge=MG_VALENCE*Q;
		radius=1e-10*MG_RADIUS;
		diffusion_coeff=1e-9*MG_DIFF_COEFF;
		mass=1e-27*MG_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==9){
		name="CA";
		charge=double(2.00)*Q;
		valence=2.00;
		mass=67.0353e-27;
		radius=0.99e-10;
		charmm_half_radius = 1.367e-10; // [m]
		charmm_eps = -0.120 * 4186.8 / AVOGADRO; //[J]
		diffusion_coeff=0.79e-9;
//-------------------------------------------------------------- new GUI
		valence=CA_VALENCE;
		charge=CA_VALENCE*Q;
		radius=1e-10*CA_RADIUS;
		diffusion_coeff=1e-9*CA_DIFF_COEFF;
		mass=1e-27*CA_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==10){
		name="SR";
		charge=double(2.00)*Q;
		valence=2.00;
		mass=0.00;
		radius=1.13e-10;
		diffusion_coeff=0.00;
	}
	else if(type==11){
		name="BA";
		charge=double(2.00)*Q;
		valence=2.00;
		mass=0.00;
		radius=1.35e-10;
		diffusion_coeff=0.00;
	}
	else if(type==12){
		name="MN";
		charge=double(2.00)*Q;
		valence=2.00;
		mass=0.00;
		radius=0.80e-10;
		diffusion_coeff=0.00;
	}
	else if(type==13){
		name="CO";
		charge=double(2.00)*Q;
		valence=2.00;
		mass=0.00;
		radius=0.74e-10;
		diffusion_coeff=0.00;
	}
	else if(type==14){
		name="NI";
		charge=double(2.00)*Q;
		valence=2.00;
		mass=0.00;
		radius=0.72e-10;
		diffusion_coeff=0.00;
	}
	else if(type==15){
		name="ZN";
		charge=double(2.00)*Q;
		valence=2.00;
		mass=0.00;
		radius=0.74e-10;
		diffusion_coeff=0.00;
	}
	else if(type==16){
		name="F";
		charge=-Q;
		valence=-1.00;
		mass=0.00;
		radius=1.36e-10;
		diffusion_coeff=0.00;
	}
	else if(type==17){
		name="CL";
		charge=-Q;
		valence=-1.00;
		mass=59.2995e-27;
		radius=1.81e-10;
		charmm_half_radius = 2.27e-10; // [m]
		charmm_eps = -0.150 * 4186.8 / AVOGADRO; //[J]
		diffusion_coeff=2.03e-9;
//-------------------------------------------------------------- new GUI
		valence=CL_VALENCE;
		charge=CL_VALENCE*Q;
		radius=1e-10*CL_RADIUS;
		diffusion_coeff=1e-9*CL_DIFF_COEFF;
		mass=1e-27*CL_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==18){
		name="BR";
		charge=-Q;
		valence=-1.00;
		mass=0.00;
		radius=1.95e-10;
		diffusion_coeff=0.00;
	}
	else if(type==19){
		name="I";
		charge=-Q;
		valence=-1.00;
		mass=0.00;
		radius=2.16e-10;
		diffusion_coeff=0.00;
	}
	else if(type==31){
		name="J11";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX1_ION1_VALENCE;
		charge=BOX1_ION1_VALENCE*Q;
		radius=1e-10*BOX1_ION1_RADIUS;
		diffusion_coeff=1e-9*BOX1_ION1_DIFF_COEFF;
		mass=1e-27*BOX1_ION1_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==32){
		name="J12";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX1_ION2_VALENCE;
		charge=BOX1_ION2_VALENCE*Q;
		radius=1e-10*BOX1_ION2_RADIUS;
		diffusion_coeff=1e-9*BOX1_ION2_DIFF_COEFF;
		mass=1e-27*BOX1_ION2_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==33){
		name="J13";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX1_ION3_VALENCE;
		charge=BOX1_ION3_VALENCE*Q;
		radius=1e-10*BOX1_ION3_RADIUS;
		diffusion_coeff=1e-9*BOX1_ION3_DIFF_COEFF;
		mass=1e-27*BOX1_ION3_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==34){
		name="J14";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX1_ION4_VALENCE;
		charge=BOX1_ION4_VALENCE*Q;
		radius=1e-10*BOX1_ION4_RADIUS;
		diffusion_coeff=1e-9*BOX1_ION4_DIFF_COEFF;
		mass=1e-27*BOX1_ION4_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==35){
		name="J15";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX1_ION5_VALENCE;
		charge=BOX1_ION5_VALENCE*Q;
		radius=1e-10*BOX1_ION5_RADIUS;
		diffusion_coeff=1e-9*BOX1_ION5_DIFF_COEFF;
		mass=1e-27*BOX1_ION5_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==36){
		name="J16";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX1_ION6_VALENCE;
		charge=BOX1_ION6_VALENCE*Q;
		radius=1e-10*BOX1_ION6_RADIUS;
		diffusion_coeff=1e-9*BOX1_ION6_DIFF_COEFF;
		mass=1e-27*BOX1_ION6_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==41){
		name="J21";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX2_ION1_VALENCE;
		charge=BOX2_ION1_VALENCE*Q;
		radius=1e-10*BOX2_ION1_RADIUS;
		diffusion_coeff=1e-9*BOX2_ION1_DIFF_COEFF;
		mass=1e-27*BOX2_ION1_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==42){
		name="J22";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX2_ION2_VALENCE;
		charge=BOX2_ION2_VALENCE*Q;
		radius=1e-10*BOX2_ION2_RADIUS;
		diffusion_coeff=1e-9*BOX2_ION2_DIFF_COEFF;
		mass=1e-27*BOX2_ION2_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==43){
		name="J23";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX2_ION3_VALENCE;
		charge=BOX2_ION3_VALENCE*Q;
		radius=1e-10*BOX2_ION3_RADIUS;
		diffusion_coeff=1e-9*BOX2_ION3_DIFF_COEFF;
		mass=1e-27*BOX2_ION3_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==44){
		name="J24";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX2_ION4_VALENCE;
		charge=BOX2_ION4_VALENCE*Q;
		radius=1e-10*BOX2_ION4_RADIUS;
		diffusion_coeff=1e-9*BOX2_ION4_DIFF_COEFF;
		mass=1e-27*BOX2_ION4_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==45){
		name="J25";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX2_ION5_VALENCE;
		charge=BOX2_ION5_VALENCE*Q;
		radius=1e-10*BOX2_ION5_RADIUS;
		diffusion_coeff=1e-9*BOX2_ION5_DIFF_COEFF;
		mass=1e-27*BOX2_ION5_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==46){
		name="J26";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX2_ION6_VALENCE;
		charge=BOX2_ION6_VALENCE*Q;
		radius=1e-10*BOX2_ION6_RADIUS;
		diffusion_coeff=1e-9*BOX2_ION6_DIFF_COEFF;
		mass=1e-27*BOX2_ION6_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==51){
		name="J31";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX3_ION1_VALENCE;
		charge=BOX3_ION1_VALENCE*Q;
		radius=1e-10*BOX3_ION1_RADIUS;
		diffusion_coeff=1e-9*BOX3_ION1_DIFF_COEFF;
		mass=1e-27*BOX3_ION1_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==52){
		name="J32";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX3_ION2_VALENCE;
		charge=BOX3_ION2_VALENCE*Q;
		radius=1e-10*BOX3_ION2_RADIUS;
		diffusion_coeff=1e-9*BOX3_ION2_DIFF_COEFF;
		mass=1e-27*BOX3_ION2_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==53){
		name="J33";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX3_ION3_VALENCE;
		charge=BOX3_ION3_VALENCE*Q;
		radius=1e-10*BOX3_ION3_RADIUS;
		diffusion_coeff=1e-9*BOX3_ION3_DIFF_COEFF;
		mass=1e-27*BOX3_ION3_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==54){
		name="J34";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX3_ION4_VALENCE;
		charge=BOX3_ION4_VALENCE*Q;
		radius=1e-10*BOX3_ION4_RADIUS;
		diffusion_coeff=1e-9*BOX3_ION4_DIFF_COEFF;
		mass=1e-27*BOX3_ION4_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==55){
		name="J35";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX3_ION5_VALENCE;
		charge=BOX3_ION5_VALENCE*Q;
		radius=1e-10*BOX3_ION5_RADIUS;
		diffusion_coeff=1e-9*BOX3_ION5_DIFF_COEFF;
		mass=1e-27*BOX3_ION5_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==56){
		name="J36";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX3_ION6_VALENCE;
		charge=BOX3_ION6_VALENCE*Q;
		radius=1e-10*BOX3_ION6_RADIUS;
		diffusion_coeff=1e-9*BOX3_ION6_DIFF_COEFF;
		mass=1e-27*BOX3_ION6_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==61){
		name="J41";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX4_ION1_VALENCE;
		charge=BOX4_ION1_VALENCE*Q;
		radius=1e-10*BOX4_ION1_RADIUS;
		diffusion_coeff=1e-9*BOX4_ION1_DIFF_COEFF;
		mass=1e-27*BOX4_ION1_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==62){
		name="J42";
		charge=double(+0.51)*Q;
		valence=0.51;
		mass=12.011e-27;
		radius=2.00e-10;
		charmm_half_radius = 2.0e-10; // [m]
		charmm_eps = -0.11 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX4_ION2_VALENCE;
		charge=BOX4_ION2_VALENCE*Q;
		radius=1e-10*BOX4_ION2_RADIUS;
		diffusion_coeff=1e-9*BOX4_ION2_DIFF_COEFF;
		mass=1e-27*BOX4_ION2_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==63){
		name="J43";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX4_ION3_VALENCE;
		charge=BOX4_ION3_VALENCE*Q;
		radius=1e-10*BOX4_ION3_RADIUS;
		diffusion_coeff=1e-9*BOX4_ION3_DIFF_COEFF;
		mass=1e-27*BOX4_ION3_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==64){
		name="J44";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX4_ION4_VALENCE;
		charge=BOX4_ION4_VALENCE*Q;
		radius=1e-10*BOX4_ION4_RADIUS;
		diffusion_coeff=1e-9*BOX4_ION4_DIFF_COEFF;
		mass=1e-27*BOX4_ION4_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==65){
		name="J45";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX4_ION5_VALENCE;
		charge=BOX4_ION5_VALENCE*Q;
		radius=1e-10*BOX4_ION5_RADIUS;
		diffusion_coeff=1e-9*BOX4_ION5_DIFF_COEFF;
		mass=1e-27*BOX4_ION5_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==66){
		name="J46";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX4_ION6_VALENCE;
		charge=BOX4_ION6_VALENCE*Q;
		radius=1e-10*BOX4_ION6_RADIUS;
		diffusion_coeff=1e-9*BOX4_ION6_DIFF_COEFF;
		mass=1e-27*BOX4_ION6_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==71){
		name="J51";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX5_ION1_VALENCE;
		charge=BOX5_ION1_VALENCE*Q;
		radius=1e-10*BOX5_ION1_RADIUS;
		diffusion_coeff=1e-9*BOX5_ION1_DIFF_COEFF;
		mass=1e-27*BOX5_ION1_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==72){
		name="J52";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX5_ION2_VALENCE;
		charge=BOX5_ION2_VALENCE*Q;
		radius=1e-10*BOX5_ION2_RADIUS;
		diffusion_coeff=1e-9*BOX5_ION2_DIFF_COEFF;
		mass=1e-27*BOX5_ION2_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==73){
		name="J53";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX5_ION3_VALENCE;
		charge=BOX5_ION3_VALENCE*Q;
		radius=1e-10*BOX5_ION3_RADIUS;
		diffusion_coeff=1e-9*BOX5_ION3_DIFF_COEFF;
		mass=1e-27*BOX5_ION3_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==74){
		name="J54";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX5_ION4_VALENCE;
		charge=BOX5_ION4_VALENCE*Q;
		radius=1e-10*BOX5_ION4_RADIUS;
		diffusion_coeff=1e-9*BOX5_ION4_DIFF_COEFF;
		mass=1e-27*BOX5_ION4_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==75){
		name="J55";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX5_ION5_VALENCE;
		charge=BOX5_ION5_VALENCE*Q;
		radius=1e-10*BOX5_ION5_RADIUS;
		diffusion_coeff=1e-9*BOX5_ION5_DIFF_COEFF;
		mass=1e-27*BOX5_ION5_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==76){
		name="J56";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX5_ION6_VALENCE;
		charge=BOX5_ION6_VALENCE*Q;
		radius=1e-10*BOX5_ION6_RADIUS;
		diffusion_coeff=1e-9*BOX5_ION6_DIFF_COEFF;
		mass=1e-27*BOX5_ION6_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==81){
		name="J61";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX6_ION1_VALENCE;
		charge=BOX6_ION1_VALENCE*Q;
		radius=1e-10*BOX6_ION1_RADIUS;
		diffusion_coeff=1e-9*BOX6_ION1_DIFF_COEFF;
		mass=1e-27*BOX6_ION1_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==82){
		name="J62";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX6_ION2_VALENCE;
		charge=BOX6_ION2_VALENCE*Q;
		radius=1e-10*BOX6_ION2_RADIUS;
		diffusion_coeff=1e-9*BOX6_ION2_DIFF_COEFF;
		mass=1e-27*BOX6_ION2_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==83){
		name="J63";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX6_ION3_VALENCE;
		charge=BOX6_ION3_VALENCE*Q;
		radius=1e-10*BOX6_ION3_RADIUS;
		diffusion_coeff=1e-9*BOX6_ION3_DIFF_COEFF;
		mass=1e-27*BOX6_ION3_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==84){
		name="J64";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX6_ION4_VALENCE;
		charge=BOX6_ION4_VALENCE*Q;
		radius=1e-10*BOX6_ION4_RADIUS;
		diffusion_coeff=1e-9*BOX6_ION4_DIFF_COEFF;
		mass=1e-27*BOX6_ION4_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==85){
		name="J65";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX6_ION5_VALENCE;
		charge=BOX6_ION5_VALENCE*Q;
		radius=1e-10*BOX6_ION5_RADIUS;
		diffusion_coeff=1e-9*BOX6_ION5_DIFF_COEFF;
		mass=1e-27*BOX6_ION5_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==86){
		name="J65";
		charge=double(-0.51)*Q;
		valence=-0.51;
		mass=15.999e-27;
		radius=1.70e-10;
		charmm_half_radius = 1.7e-10; // [m]
		charmm_eps = -0.12 * 4186.8 / AVOGADRO; // [J]
		diffusion_coeff=2.299e-9; // [m^2/s] Diffusion coefficient of water at 25 C
//-------------------------------------------------------------- new GUI
		valence=BOX6_ION6_VALENCE;
		charge=BOX6_ION6_VALENCE*Q;
		radius=1e-10*BOX6_ION6_RADIUS;
		diffusion_coeff=1e-9*BOX6_ION6_DIFF_COEFF;
		mass=1e-27*BOX6_ION6_MASS;	
//-------------------------------------------------------------- new GUI		
	}
	else if(type==92){
		name="MEMBRANE";
		charge=double(1.00)*Q;
		valence=1.00;
		mass=0.00;
		radius=0.00e-10;
		charmm_half_radius = 0.00e-10; // [m]
		charmm_eps = -0.07 * 4186.8 / AVOGADRO; //[J]
		diffusion_coeff=0.00;
	}
	else{
		
		name="UNKNOWN";
		charge=0;
		valence=0.00;
		mass=0.00;
		radius=0.00;
		diffusion_coeff=0.00;
	
		//~ cout << "Ion of kind " << type << " does not exist yet!" <<endl;
		//~ exit(100);
	}
	
	EPS_ION=PRM.EPS_W;
	
	DW_charge=charge/EPS_ION;
	DW_valence=valence/EPS_ION;
	pm_radius=radius*1e12;
	
	set_brownian_parameters();
	
	x=0.00;
	y=0.00;
	z=0.00;
	
	x_prev=0.00;
	y_prev=0.00;
	z_prev=0.00;
	
	x_next=0.00;
	y_next=0.00;
	z_next=0.00;

	for(int i=0; i<3; i++){
		field[i]=0.00;
		velocity[i]=0.00;
		force[i]=0.00;
		force_old[i]=0.00;
		X_n[i]=0.00;
		X_n_old[i]=0.00;
		X_n_minus[i]=0.00;
	}
	
}
/*
Ion::Ion(Parameters PRM, int type, double valence, double radius, double diffusion_coeff, double mass){

	reset_ion();
	
	kind=type;
	
	valence=valence;
	charge=valence*Q;
	radius=1e-10*radius;
	diffusion_coeff=1e-9*diffusion_coeff;
	mass=1e-27*mass;	
	
	if(type==0)name="O";
	else if(type==1)name="H";
	else if(type==2)name="LI";
	else if(type==3)name="NA";
	else if(type==4)name="K";
	else if(type==5)name="RB";
	else if(type==6)name="CS";
	else if(type==7)name="TL";
	else if(type==8)name="MG";
	else if(type==9)name="CA";
	else if(type==10)name="SR";
	else if(type==11)name="BA";
	else if(type==12)name="MN";
	else if(type==13)name="CO";
	else if(type==14)name="NI";
	else if(type==15)name="ZN";
	else if(type==16)name="F";
	else if(type==17)name="CL";
	else if(type==18)name="BR";
	else if(type==19)name="I";
	
	else if(type==31)name="J11";
	else if(type==32)name="J12";
	else if(type==33)name="J13";
	else if(type==34)name="J14";
	else if(type==35)name="J15";
	else if(type==36)name="J16";
	
	else if(type==41)name="J21";
	else if(type==42)name="J22";
	else if(type==43)name="J23";
	else if(type==44)name="J24";
	else if(type==45)name="J25";
	else if(type==46)name="J26";
	
	else if(type==51)name="J31";
	else if(type==52)name="J32";
	else if(type==53)name="J33";
	else if(type==54)name="J34";
	else if(type==55)name="J35";
	else if(type==56)name="J36";
	
	else if(type==61)name="J41";
	else if(type==62)name="J42";
	else if(type==63)name="J43";
	else if(type==64)name="J44";
	else if(type==65)name="J45";
	else if(type==66)name="J46";
	
	else if(type==71)name="J51";
	else if(type==72)name="J52";
	else if(type==73)name="J53";
	else if(type==74)name="J54";
	else if(type==75)name="J55";
	else if(type==76)name="J56";
	
	else if(type==81)name="J61";
	else if(type==82)name="J62";
	else if(type==83)name="J63";
	else if(type==84)name="J64";
	else if(type==85)name="J65";
	else if(type==86)name="J66";

	EPS_ION=PRM.EPS_W;
	
	DW_charge=charge/EPS_ION;
	DW_valence=valence/EPS_ION;
	pm_radius=radius*1e12;
	
	set_brownian_parameters(PRM);
	created=true;
	displaced=false;
	
	x=0.00;
	y=0.00;
	z=0.00;
	
	x_prev=0.00;
	y_prev=0.00;
	z_prev=0.00;
	
	x_next=0.00;
	y_next=0.00;
	z_next=0.00;

	potential=0.00;
	energy=0.00;
	
	for(int i=0; i<3; i++){
		field.at(i)=0.00;
		velocity.at(i)=0.00;
		force.at(i)=0.00;
		force_old.at(i)=0.00;
		X_n.at(i)=0.00;
		X_n_old.at(i)=0.00;
		X_n_minus.at(i)=0.00;
	}
	
}
*/
Ion::~Ion(){
	
}
void Ion::reset_ion(){
	
	kind=-1;
	
	name="";
	charge=0.00;
	valence=0.00;
	mass=0.00;
	radius=0.00;
	charmm_half_radius=0.00; // 
	charmm_eps=0.00; // 
	diffusion_coeff=0.00;
	
	x=0.00;
	y=0.00;
	z=0.00;
	
	x_prev=0.00;
	y_prev=0.00;
	z_prev=0.00;
	
	x_next=0.00;
	y_next=0.00;
	z_next=0.00;
	
	EPS_ION=0.00;

	for(int i=0; i<3; i++){
		field[i]=0.00;
		velocity[i]=0.00;
		force[i]=0.00;
		force_old[i]=0.00;
		X_n[i]=0.00;
		X_n_old[i]=0.00;
		X_n_minus[i]=0.00;
	}
	
	side=-1;

	return;
}

void Ion::set_brownian_parameters(){
	
	double gamma=PRM.kT/(mass*diffusion_coeff);
	double gamma_dt=gamma*PRM.DELTA_T;
	
	double exp_gamma_dt=exp(gamma_dt);
	double exp_minus_gamma_dt=exp(-gamma_dt);
	double exp_2gamma_dt=exp(double(2.00)*gamma_dt);
	double exp_minus_2gamma_dt=exp(-double(2.00)*gamma_dt);
	
	double temp1=PRM.DELTA_T/(mass*gamma);
	double temp2=PRM.DELTA_T/(mass*gamma*gamma);
	double temp3=(double(1.00))/(mass*gamma*gamma);
	double temp4=(double(1.00))/(mass*gamma*gamma*gamma);
	
	double temp11=double(16.00)*(exp_gamma_dt+exp_minus_gamma_dt)-double(4.00)*(exp_2gamma_dt+exp_minus_2gamma_dt);
	double temp12=-double(24.00)-double(4.00)*gamma_dt*(exp_gamma_dt-exp_minus_gamma_dt);
	double temp13=double(2.00)*gamma_dt*(exp_2gamma_dt-exp_minus_2gamma_dt);
	E=temp11+temp12+temp13;
	C=double(2.00)*gamma_dt-double(3.00)+double(4.00)*exp_minus_gamma_dt-exp_minus_2gamma_dt;
	G=exp_gamma_dt-double(2.00)*gamma_dt-exp_minus_gamma_dt;
	H=gamma_dt/(exp_gamma_dt-exp_minus_gamma_dt);
	
	K_pos_1=double(1.00)+exp_minus_gamma_dt;
	K_pos_2=-exp_minus_gamma_dt;
	K_pos_3=temp1*(double(1.00)-exp_minus_gamma_dt);
	K_pos_4=temp2*(double(0.5)*gamma_dt*(double(1.00)+exp_minus_gamma_dt)-(double(1.00)-exp_minus_gamma_dt));
	K_pos_5=exp_minus_gamma_dt;

	K_vel_1=temp3*G;
	K_vel_2=-temp4*G;
	K_vel_3=H/PRM.DELTA_T;
	
	K_Y=sqrt((PRM.kT*temp3*E)/C);
	K_X=sqrt(PRM.kT*temp3*C);
	
	K_step0_1=(double(1.00)-exp_minus_gamma_dt)/gamma;
	K_step0_2=temp3*(gamma_dt-(double(1.00)-exp_minus_gamma_dt));
	
	bro_KF=(diffusion_coeff*PRM.DELTA_T)/(PRM.kT);
	bro_KR=sqrt(double(2.00)*diffusion_coeff*PRM.DELTA_T);
	
	//~ cout << "gamma: " << gamma <<endl;
	//~ cout << "gamma_dt: " << gamma_dt <<endl;
	//~ cout << "diffusion_coeff: " << diffusion_coeff <<endl;
	//~ cout << "mass: " << mass <<endl;
	//~ cout << "PRM.kT: " << PRM.kT <<endl<<endl;
	//~ sleep(1);
	
	
	return;
}

ostream& operator<<(ostream& stream, Ion& i){

	stream << endl << "########################################" <<endl;
	stream << "Ion" <<endl;
	stream << "kind/name:\t"<<i.kind<<"\t" << i.name<<endl;
	stream << "charge:\t"<<i.charge<<"\t"<<i.DW_charge<<endl;
	stream << "valence:\t"<<i.valence<<"\t"<<i.DW_valence<<endl;
	stream << "mass/radius:\t"<<i.mass<<"\t"<<i.radius<<endl;
	stream << "diffusion_coeff/EPS_ION:\t\t"<<i.diffusion_coeff<<"\t"<<i.EPS_ION<<endl;
	stream << "prev pos (pm):\t"<<1e12*i.x_prev << "\t" << 1e12*i.y_prev << "\t" << 1e12*i.z_prev << endl;
	stream << "position (pm):\t"<<1e12*i.x << "\t" << 1e12*i.y << "\t" << 1e12*i.z << endl;
	stream << "next pos (pm):\t"<<1e12*i.x_next << "\t" << 1e12*i.y_next << "\t" << 1e12*i.z_next << endl;
	stream << "velocity:\t"<<i.velocity[0] << "\t" << i.velocity[1] << "\t" << i.velocity[2] << endl;
	stream << "field:\t\t"<<i.field[0] << "\t" << i.field[1] << "\t" << i.field[2] << endl;
	stream << "force:\t"<<i.force[0] << "\t" << i.force[1] << "\t" << i.force[2] << endl;
	stream << "force_old:\t"<<i.force_old[0] << "\t" << i.force_old[1] << "\t" << i.force_old[2] << endl;
	stream << "X_n:\t"<<i.X_n[0] << "\t" << i.X_n[1] << "\t" << i.X_n[2] << endl;
	stream << "X_n_old:\t"<<i.X_n_old[0] << "\t" << i.X_n_old[1] << "\t" << i.X_n_old[2] << endl;
	stream << "X_n_minus:\t"<<i.X_n_minus[0] << "\t" << i.X_n_minus[1] << "\t" << i.X_n_minus[2] << endl;
	stream << "K_pos:\t"<<i.K_pos_1 << "\t" << i.K_pos_2 << "\t" << i.K_pos_3 << "\t" << i.K_pos_4 << "\t" << i.K_pos_5 << endl;
	stream << "K_vel:\t"<<i.K_vel_1 << "\t" << i.K_vel_2 << "\t" << i.K_vel_3 << endl;
	stream << "K_Y/K_X:\t"<<i.K_Y << "\t" << i.K_X << endl;
	stream << "K_step0_1/02:\t"<<i.K_step0_1 << "\t" << i.K_step0_2 << endl;
	stream << "C/G/E/H:\t"<<i.C << "\t" << i.G << "\t" << i.E << "\t" << i.H << endl;
	return stream;
}

//########################################
// Class Surface
//########################################
Surface::Surface(){
	reset_surface();
}

Surface::~Surface(){
	
}

void Surface::reset_surface(){
	
	pA.reset_point();
	pB.reset_point();
	pC.reset_point();
	pD.reset_point();
	center.reset_point();
	
	normal[0]=0.00;
	normal[1]=0.00;
	normal[2]=0.00;
	
	area=0.00;
	
	return;
}

ostream& operator<<(ostream& stream, Surface& s){

	stream << endl << "########################################" <<endl;
	stream << "Surface" <<endl;
	
	return stream;
}

//########################################
// Class SurfaceLight
//########################################
SurfaceLight::SurfaceLight(){
	reset_surfaceLight();
}

SurfaceLight::~SurfaceLight(){
	
}

void SurfaceLight::reset_surfaceLight(){
	
	center_x=0.00;
	center_y=0.00;
	center_z=0.00;
	
	normal[0]=0.00;
	normal[1]=0.00;
	normal[2]=0.00;
	
	area=0.00;
	
	return;
}

ostream& operator<<(ostream& stream, SurfaceLight& s){

	stream << endl << "########################################" <<endl;
	stream << "SurfaceLight" <<endl;
	
	return stream;
}


//########################################
// Class Tile
//########################################
Tile::Tile(){
	reset_Tile();
}

Tile::~Tile(){
	
}

void Tile::reset_Tile(){
	
	pA.reset_point();
	pB.reset_point();
	pC.reset_point();
	pD.reset_point();
	center.reset_point();
	
	normal[0]=0.00;
	normal[1]=0.00;
	normal[2]=0.00;
	
	area=0.00;
	
	return;
}

ostream& operator<<(ostream& stream, Tile& t){

	stream << endl << "########################################" <<endl;
	stream << "Tile" <<endl;
	
	return stream;
}

//########################################
// Class Ion_box
//########################################
Ion_box::Ion_box(){
	reset_ion_box();
}

Ion_box::~Ion_box(){
	
}

void Ion_box::reset_ion_box(){
	
	index=0;
		
	//~ MIN_X=0.00;
	//~ MIN_Y=0.00;
	MIN_Z=0.00;
	//~ MAX_X=0.00;
	//~ MAX_Y=0.00;
	MAX_Z=0.00;
	
	TOT_CHARGE=0.00;
	
	is.clear();
	in.clear();
	ion_indexes.clear();
	
}


//########################################
// Class Box
//########################################
Box::Box(){
	reset_box();
}

Box::~Box(){
	
}

void Box::reset_box(){
	
	min_x=0.00;
	min_y=0.00;
	min_z=0.00;
	max_x=0.00;
	max_y=0.00;
	max_z=0.00;
	volume_dimension_x=max_x-min_x;
	volume_dimension_y=max_y-min_y;
	volume_dimension_z=max_z-min_z;
	permittivity=1.00;
	
}

ostream& operator<<(ostream& stream, Box& box){

	stream << endl << "########################################" <<endl;
	stream << "Box" <<endl;
	stream << "permittivity:\t"<<box.permittivity<<endl;
	stream << "x:\t"<< box.min_x << "\t" << box.max_x << "\t" << box.volume_dimension_x << endl;
	stream << "y:\t"<< box.min_y << "\t" << box.max_y << "\t" << box.volume_dimension_y << endl;
	stream << "z:\t"<< box.min_z << "\t" << box.max_z << "\t" << box.volume_dimension_z << endl;

	return stream;
}
