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
#include <iomanip>

#include "constants.h"
#include "ions_properties.h"
#include "utils.h"
#include "classes.h"
#include "file_functions.h"
#include "input_output.h"
#include "sim_structures.h"


void retrieve_parameters(string conf_file){
	if(fileExists(conf_file)){
		if(verifyExtension(conf_file, "conf")){
			string buffer="";
			ifstream fin; 
			char *tfn = new char[conf_file.length()+1];
			strcpy(tfn, conf_file.c_str());     
			fin.open(tfn, ifstream::in);  
			if(!fin){
				cout <<"ERROR: it is not possible to open the file: " << conf_file << endl;
				exit(2);
			}
			while(!fin.eof()){
				if(getline(fin, buffer)){
					if(buffer.substr(0,1).compare("#")!=0){
						size_t found=buffer.find("=");
						if(found!=string::npos){
							removeAllWhite(buffer);
							if(getTokenbyNumber(buffer, "=", 1).compare("SR_METHOD")==0){
								PRM.SR_METHOD=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("MD_MAP_MIN_Z")==0){
								PRM.MD_MAP_MIN_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("MD_MAP_MAX_Z")==0){
								PRM.MD_MAP_MAX_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("SIM_TYPE")==0){
								PRM.SIM_TYPE=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("PREFIX")==0){
								PRM.PREFIX=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONTROL_CELL_WIDTH")==0){
								PRM.CONTROL_CELL_WIDTH=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("OUTER_REGION_WIDTH")==0){
								PRM.OUTER_REGION_WIDTH=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BATH_WIDTH")==0){
								PRM.BATH_WIDTH=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("MEMBRANE_WIDTH")==0){
								PRM.MEMBRANE_WIDTH=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("SIM_DOMAIN_WIDTH_X")==0){
								PRM.SIM_DOMAIN_WIDTH_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("SIM_DOMAIN_WIDTH_Y")==0){
								PRM.SIM_DOMAIN_WIDTH_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("APPLIED_POTENTIAL")==0){
								PRM.APPLIED_POTENTIAL=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("PBC")==0){
								PRM.PBC=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("ION_RECYCLING")==0){
								PRM.ION_RECYCLING=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONC_LEFT_KCL")==0){
								PRM.CONC_LEFT_KCL=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONC_LEFT_NACL")==0){
								PRM.CONC_LEFT_NACL=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONC_LEFT_CACL2")==0){
								PRM.CONC_LEFT_CACL2=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONC_LEFT_MGCL2")==0){
								PRM.CONC_LEFT_MGCL2=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONC_RIGHT_KCL")==0){
								PRM.CONC_RIGHT_KCL=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONC_RIGHT_NACL")==0){
								PRM.CONC_RIGHT_NACL=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONC_RIGHT_CACL2")==0){
								PRM.CONC_RIGHT_CACL2=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CONC_RIGHT_MGCL2")==0){
								PRM.CONC_RIGHT_MGCL2=1e-3*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("SHORT_RANGE_EXP")==0){
								PRM.SHORT_RANGE_EXP=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("EPS_W")==0){
								PRM.EPS_W=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("EPS_MEM")==0){
								PRM.EPS_MEM=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("DIFF_COEFF_IN_CHANNEL")==0){
								PRM.DIFF_COEFF_IN_CHANNEL=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("DELTA_T")==0){
								PRM.DELTA_T=1e-15*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("TEMPERATURE")==0){
								PRM.TEMPERATURE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("HISTORY_SIZE")==0){
								PRM.HISTORY_SIZE=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("SEED")==0){
								PRM.SEED=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("PREP_STEPS")==0){
								PRM.PREP_STEPS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("SIM_STEPS")==0){
								PRM.SIM_STEPS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("STATS_OUT_FREQ")==0){
								PRM.STATS_OUT_FREQ=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("STATS_DZ")==0){
								PRM.STATS_DZ=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("channel_pdb_files")==0){
								PRM.channel_pdb_files=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("trajectory")==0){
								PRM.trajectory=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("flux")==0){
								PRM.flux=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("rdf")==0){
								PRM.rdf=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("vel_distribution")==0){
								PRM.vel_distribution=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("mean_square_displ")==0){
								PRM.mean_square_displ=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("induced_charge")==0){
								PRM.induced_charge=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("potential")==0){
								PRM.potential=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("concentrations")==0){
								PRM.concentrations=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("currents_ZT")==0){
								PRM.currents_ZT=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("LEFT_VESTIBULE_CURVATURE_RADIUS")==0){
								PRM.LEFT_VESTIBULE_CURVATURE_RADIUS=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("LEFT_VESTIBULE_MIN_CHANNEL_RADIUS")==0){
								PRM.LEFT_VESTIBULE_MIN_CHANNEL_RADIUS=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("RIGHT_VESTIBULE_CURVATURE_RADIUS")==0){
								PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS")==0){
								PRM.RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CHANNEL_PROFILE_POINT")==0){
								RadialPoint rp;
								rp.z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
								rp.r=double(100.00)*atof(getTokenbyNumber(buffer, "=", 3).c_str());
								
								if(PRM.channel_profile_points.empty()){
									PRM.channel_profile_points.push_back(rp);
								}
								else{
									bool canInsert=true;
									for(int i=0; i<PRM.channel_profile_points.size(); i++){
										if(int(PRM.channel_profile_points.at(i).z)==int(rp.z)){
											canInsert=false;
										}
									}
									if(canInsert){
										PRM.channel_profile_points.push_back(rp);
									}
								}
								
								if(PRM.channel_profile_points.size()>1){
									int flag = 1;    // set flag to 1 to start first pass
									RadialPoint tempRp;             // holding variable

									for(int i=1; (i<=PRM.channel_profile_points.size()) && flag; i++){
										flag=0;
										for(int j=0;j<(PRM.channel_profile_points.size()-1);j++){
											if (PRM.channel_profile_points.at(j+1).z<PRM.channel_profile_points.at(j).z){ 
												tempRp.z=PRM.channel_profile_points.at(j).z;
												tempRp.r=PRM.channel_profile_points.at(j).r;
												PRM.channel_profile_points.at(j).z=PRM.channel_profile_points.at(j+1).z;
												PRM.channel_profile_points.at(j).r=PRM.channel_profile_points.at(j+1).r;
												PRM.channel_profile_points.at(j+1).z=tempRp.z;
												PRM.channel_profile_points.at(j+1).r=tempRp.r;
												flag=1;               // indicates that a swap occurred.
											}
										}
									}
								}
								
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("TILES_PER_RING")==0){
								PRM.TILES_PER_RING=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("NUM_OF_DIV")==0){
								PRM.NUM_OF_DIV=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("NUM_OF_SUB_DIV")==0){
								PRM.NUM_OF_SUB_DIV=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("MAX_TILE_WIDTH")==0){
								PRM.MAX_TILE_WIDTH=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CHARGE_RING")==0){
								PRM.charge_ring_z.push_back(double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str()));
								PRM.charge_ring_r.push_back(double(100.00)*atof(getTokenbyNumber(buffer, "=", 3).c_str()));
								PRM.charge_ring_n.push_back(atoi(getTokenbyNumber(buffer, "=", 4).c_str()));
								PRM.charge_ring_q.push_back(atof(getTokenbyNumber(buffer, "=", 5).c_str()));
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("FIXED_CHARGES_FILE")==0){
								PRM.FIXED_CHARGES_FILE=getTokenbyNumber(buffer, "=", 2);
								cout << "PRM.FIXED_CHARGES_FILE:\t" << PRM.FIXED_CHARGES_FILE << endl;
							}
//=============================================================================== ION BOX1 begin							
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_MIN_X")==0){
								PRM.BOX1_MIN_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_MAX_X")==0){
								PRM.BOX1_MAX_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_MIN_Y")==0){
								PRM.BOX1_MIN_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_MAX_Y")==0){
								PRM.BOX1_MAX_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_MIN_Z")==0){
								PRM.BOX1_MIN_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_MAX_Z")==0){
								PRM.BOX1_MAX_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_IS1")==0){
								PRM.BOX1_IS1=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_N1")==0){
								PRM.BOX1_N1=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_IS2")==0){
								PRM.BOX1_IS2=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_N2")==0){
								PRM.BOX1_N2=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_IS3")==0){
								PRM.BOX1_IS3=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_N3")==0){
								PRM.BOX1_N3=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_IS4")==0){
								PRM.BOX1_IS4=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_N4")==0){
								PRM.BOX1_N4=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_IS5")==0){
								PRM.BOX1_IS5=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_N5")==0){
								PRM.BOX1_N5=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_IS6")==0){
								PRM.BOX1_IS6=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_N6")==0){
								PRM.BOX1_N6=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
//=============================================================================== ION BOX1 end
//=============================================================================== ION BOX2 begin							
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_MIN_X")==0){
								PRM.BOX2_MIN_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_MAX_X")==0){
								PRM.BOX2_MAX_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_MIN_Y")==0){
								PRM.BOX2_MIN_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_MAX_Y")==0){
								PRM.BOX2_MAX_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_MIN_Z")==0){
								PRM.BOX2_MIN_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_MAX_Z")==0){
								PRM.BOX2_MAX_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_IS1")==0){
								PRM.BOX2_IS1=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_N1")==0){
								PRM.BOX2_N1=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_IS2")==0){
								PRM.BOX2_IS2=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_N2")==0){
								PRM.BOX2_N2=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_IS3")==0){
								PRM.BOX2_IS3=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_N3")==0){
								PRM.BOX2_N3=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_IS4")==0){
								PRM.BOX2_IS4=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_N4")==0){
								PRM.BOX2_N4=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_IS5")==0){
								PRM.BOX2_IS5=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_N5")==0){
								PRM.BOX2_N5=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_IS6")==0){
								PRM.BOX2_IS6=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_N6")==0){
								PRM.BOX2_N6=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
//=============================================================================== ION BOX2 end
//=============================================================================== ION BOX3 begin							
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_MIN_X")==0){
								PRM.BOX3_MIN_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_MAX_X")==0){
								PRM.BOX3_MAX_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_MIN_Y")==0){
								PRM.BOX3_MIN_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_MAX_Y")==0){
								PRM.BOX3_MAX_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_MIN_Z")==0){
								PRM.BOX3_MIN_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_MAX_Z")==0){
								PRM.BOX3_MAX_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_IS1")==0){
								PRM.BOX3_IS1=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_N1")==0){
								PRM.BOX3_N1=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_IS2")==0){
								PRM.BOX3_IS2=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_N2")==0){
								PRM.BOX3_N2=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_IS3")==0){
								PRM.BOX3_IS3=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_N3")==0){
								PRM.BOX3_N3=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_IS4")==0){
								PRM.BOX3_IS4=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_N4")==0){
								PRM.BOX3_N4=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_IS5")==0){
								PRM.BOX3_IS5=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_N5")==0){
								PRM.BOX3_N5=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_IS6")==0){
								PRM.BOX3_IS6=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_N6")==0){
								PRM.BOX3_N6=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
//=============================================================================== ION BOX3 end
//=============================================================================== ION BOX4 begin							
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_MIN_X")==0){
								PRM.BOX4_MIN_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_MAX_X")==0){
								PRM.BOX4_MAX_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_MIN_Y")==0){
								PRM.BOX4_MIN_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_MAX_Y")==0){
								PRM.BOX4_MAX_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_MIN_Z")==0){
								PRM.BOX4_MIN_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_MAX_Z")==0){
								PRM.BOX4_MAX_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_IS1")==0){
								PRM.BOX4_IS1=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_N1")==0){
								PRM.BOX4_N1=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_IS2")==0){
								PRM.BOX4_IS2=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_N2")==0){
								PRM.BOX4_N2=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_IS3")==0){
								PRM.BOX4_IS3=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_N3")==0){
								PRM.BOX4_N3=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_IS4")==0){
								PRM.BOX4_IS4=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_N4")==0){
								PRM.BOX4_N4=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_IS5")==0){
								PRM.BOX4_IS5=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_N5")==0){
								PRM.BOX4_N5=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_IS6")==0){
								PRM.BOX4_IS6=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_N6")==0){
								PRM.BOX4_N6=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
//=============================================================================== ION BOX4 end
//=============================================================================== ION BOX5 begin							
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_MIN_X")==0){
								PRM.BOX5_MIN_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_MAX_X")==0){
								PRM.BOX5_MAX_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_MIN_Y")==0){
								PRM.BOX5_MIN_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_MAX_Y")==0){
								PRM.BOX5_MAX_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_MIN_Z")==0){
								PRM.BOX5_MIN_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_MAX_Z")==0){
								PRM.BOX5_MAX_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_IS1")==0){
								PRM.BOX5_IS1=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_N1")==0){
								PRM.BOX5_N1=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_IS2")==0){
								PRM.BOX5_IS2=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_N2")==0){
								PRM.BOX5_N2=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_IS3")==0){
								PRM.BOX5_IS3=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_N3")==0){
								PRM.BOX5_N3=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_IS4")==0){
								PRM.BOX5_IS4=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_N4")==0){
								PRM.BOX5_N4=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_IS5")==0){
								PRM.BOX5_IS5=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_N5")==0){
								PRM.BOX5_N5=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_IS6")==0){
								PRM.BOX5_IS6=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_N6")==0){
								PRM.BOX5_N6=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
//=============================================================================== ION BOX5 end			
//=============================================================================== ION BOX6 begin							
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_MIN_X")==0){
								PRM.BOX6_MIN_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_MAX_X")==0){
								PRM.BOX6_MAX_X=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_MIN_Y")==0){
								PRM.BOX6_MIN_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_MAX_Y")==0){
								PRM.BOX6_MAX_Y=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_MIN_Z")==0){
								PRM.BOX6_MIN_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_MAX_Z")==0){
								PRM.BOX6_MAX_Z=double(100.00)*atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_IS1")==0){
								PRM.BOX6_IS1=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_N1")==0){
								PRM.BOX6_N1=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_IS2")==0){
								PRM.BOX6_IS2=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_N2")==0){
								PRM.BOX6_N2=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_IS3")==0){
								PRM.BOX6_IS3=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_N3")==0){
								PRM.BOX6_N3=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_IS4")==0){
								PRM.BOX6_IS4=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_N4")==0){
								PRM.BOX6_N4=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_IS5")==0){
								PRM.BOX6_IS5=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_N5")==0){
								PRM.BOX6_N5=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_IS6")==0){
								PRM.BOX6_IS6=getTokenbyNumber(buffer, "=", 2);
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_N6")==0){
								PRM.BOX6_N6=atoi(getTokenbyNumber(buffer, "=", 2).c_str());
							}
//=============================================================================== ION BOX6 end						
//=============================================================================== ION PARAMETERS begin
							else if(getTokenbyNumber(buffer, "=", 1).compare("K_VALENCE")==0){
								K_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("K_RADIUS")==0){
								K_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("K_DIFF_COEFF")==0){
								K_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("K_MASS")==0){
								K_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("NA_VALENCE")==0){
								NA_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("NA_RADIUS")==0){
								NA_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("NA_DIFF_COEFF")==0){
								NA_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("NA_MASS")==0){
								NA_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CA_VALENCE")==0){
								CA_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CA_RADIUS")==0){
								CA_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CA_DIFF_COEFF")==0){
								CA_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CA_MASS")==0){
								CA_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("MG_VALENCE")==0){
								MG_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("MG_RADIUS")==0){
								MG_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("MG_DIFF_COEFF")==0){
								MG_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("MG_MASS")==0){
								MG_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CL_VALENCE")==0){
								CL_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CL_RADIUS")==0){
								CL_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CL_DIFF_COEFF")==0){
								CL_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("CL_MASS")==0){
								CL_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION1_VALENCE")==0){
								BOX1_ION1_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION1_RADIUS")==0){
								BOX1_ION1_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION1_DIFF_COEFF")==0){
								BOX1_ION1_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION1_MASS")==0){
								BOX1_ION1_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION2_VALENCE")==0){
								BOX1_ION2_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION2_RADIUS")==0){
								BOX1_ION2_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION2_DIFF_COEFF")==0){
								BOX1_ION2_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION2_MASS")==0){
								BOX1_ION2_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION3_VALENCE")==0){
								BOX1_ION3_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION3_RADIUS")==0){
								BOX1_ION3_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION3_DIFF_COEFF")==0){
								BOX1_ION3_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION3_MASS")==0){
								BOX1_ION3_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION4_VALENCE")==0){
								BOX1_ION4_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION4_RADIUS")==0){
								BOX1_ION4_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION4_DIFF_COEFF")==0){
								BOX1_ION4_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION4_MASS")==0){
								BOX1_ION4_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION5_VALENCE")==0){
								BOX1_ION5_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION5_RADIUS")==0){
								BOX1_ION5_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION5_DIFF_COEFF")==0){
								BOX1_ION5_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION5_MASS")==0){
								BOX1_ION5_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION6_VALENCE")==0){
								BOX1_ION6_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION6_RADIUS")==0){
								BOX1_ION6_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION6_DIFF_COEFF")==0){
								BOX1_ION6_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX1_ION6_MASS")==0){
								BOX1_ION6_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION1_VALENCE")==0){
								BOX2_ION1_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION1_RADIUS")==0){
								BOX2_ION1_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION1_DIFF_COEFF")==0){
								BOX2_ION1_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION1_MASS")==0){
								BOX2_ION1_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION2_VALENCE")==0){
								BOX2_ION2_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION2_RADIUS")==0){
								BOX2_ION2_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION2_DIFF_COEFF")==0){
								BOX2_ION2_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION2_MASS")==0){
								BOX2_ION2_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION3_VALENCE")==0){
								BOX2_ION3_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION3_RADIUS")==0){
								BOX2_ION3_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION3_DIFF_COEFF")==0){
								BOX2_ION3_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION3_MASS")==0){
								BOX2_ION3_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION4_VALENCE")==0){
								BOX2_ION4_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION4_RADIUS")==0){
								BOX2_ION4_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION4_DIFF_COEFF")==0){
								BOX2_ION4_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION4_MASS")==0){
								BOX2_ION4_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION5_VALENCE")==0){
								BOX2_ION5_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION5_RADIUS")==0){
								BOX2_ION5_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION5_DIFF_COEFF")==0){
								BOX2_ION5_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION5_MASS")==0){
								BOX2_ION5_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION6_VALENCE")==0){
								BOX2_ION6_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION6_RADIUS")==0){
								BOX2_ION6_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION6_DIFF_COEFF")==0){
								BOX2_ION6_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX2_ION6_MASS")==0){
								BOX2_ION6_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION1_VALENCE")==0){
								BOX3_ION1_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION1_RADIUS")==0){
								BOX3_ION1_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION1_DIFF_COEFF")==0){
								BOX3_ION1_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION1_MASS")==0){
								BOX3_ION1_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION2_VALENCE")==0){
								BOX3_ION2_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION2_RADIUS")==0){
								BOX3_ION2_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION2_DIFF_COEFF")==0){
								BOX3_ION2_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION2_MASS")==0){
								BOX3_ION2_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION3_VALENCE")==0){
								BOX3_ION3_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION3_RADIUS")==0){
								BOX3_ION3_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION3_DIFF_COEFF")==0){
								BOX3_ION3_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION3_MASS")==0){
								BOX3_ION3_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION4_VALENCE")==0){
								BOX3_ION4_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION4_RADIUS")==0){
								BOX3_ION4_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION4_DIFF_COEFF")==0){
								BOX3_ION4_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION4_MASS")==0){
								BOX3_ION4_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION5_VALENCE")==0){
								BOX3_ION5_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION5_RADIUS")==0){
								BOX3_ION5_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION5_DIFF_COEFF")==0){
								BOX3_ION5_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION5_MASS")==0){
								BOX3_ION5_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION6_VALENCE")==0){
								BOX3_ION6_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION6_RADIUS")==0){
								BOX3_ION6_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION6_DIFF_COEFF")==0){
								BOX3_ION6_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX3_ION6_MASS")==0){
								BOX3_ION6_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION1_VALENCE")==0){
								BOX4_ION1_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION1_RADIUS")==0){
								BOX4_ION1_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION1_DIFF_COEFF")==0){
								BOX4_ION1_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION1_MASS")==0){
								BOX4_ION1_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION2_VALENCE")==0){
								BOX4_ION2_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION2_RADIUS")==0){
								BOX4_ION2_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION2_DIFF_COEFF")==0){
								BOX4_ION2_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION2_MASS")==0){
								BOX4_ION2_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION3_VALENCE")==0){
								BOX4_ION3_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION3_RADIUS")==0){
								BOX4_ION3_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION3_DIFF_COEFF")==0){
								BOX4_ION3_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION3_MASS")==0){
								BOX4_ION3_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION4_VALENCE")==0){
								BOX4_ION4_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION4_RADIUS")==0){
								BOX4_ION4_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION4_DIFF_COEFF")==0){
								BOX4_ION4_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION4_MASS")==0){
								BOX4_ION4_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION5_VALENCE")==0){
								BOX4_ION5_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION5_RADIUS")==0){
								BOX4_ION5_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION5_DIFF_COEFF")==0){
								BOX4_ION5_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION5_MASS")==0){
								BOX4_ION5_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION6_VALENCE")==0){
								BOX4_ION6_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION6_RADIUS")==0){
								BOX4_ION6_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION6_DIFF_COEFF")==0){
								BOX4_ION6_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX4_ION6_MASS")==0){
								BOX4_ION6_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION1_VALENCE")==0){
								BOX5_ION1_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION1_RADIUS")==0){
								BOX5_ION1_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION1_DIFF_COEFF")==0){
								BOX5_ION1_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION1_MASS")==0){
								BOX5_ION1_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION2_VALENCE")==0){
								BOX5_ION2_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION2_RADIUS")==0){
								BOX5_ION2_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION2_DIFF_COEFF")==0){
								BOX5_ION2_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION2_MASS")==0){
								BOX5_ION2_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION3_VALENCE")==0){
								BOX5_ION3_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION3_RADIUS")==0){
								BOX5_ION3_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION3_DIFF_COEFF")==0){
								BOX5_ION3_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION3_MASS")==0){
								BOX5_ION3_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION4_VALENCE")==0){
								BOX5_ION4_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION4_RADIUS")==0){
								BOX5_ION4_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION4_DIFF_COEFF")==0){
								BOX5_ION4_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION4_MASS")==0){
								BOX5_ION4_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION5_VALENCE")==0){
								BOX5_ION5_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION5_RADIUS")==0){
								BOX5_ION5_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION5_DIFF_COEFF")==0){
								BOX5_ION5_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION5_MASS")==0){
								BOX5_ION5_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION6_VALENCE")==0){
								BOX5_ION6_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION6_RADIUS")==0){
								BOX5_ION6_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION6_DIFF_COEFF")==0){
								BOX5_ION6_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX5_ION6_MASS")==0){
								BOX5_ION6_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION1_VALENCE")==0){
								BOX6_ION1_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION1_RADIUS")==0){
								BOX6_ION1_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION1_DIFF_COEFF")==0){
								BOX6_ION1_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION1_MASS")==0){
								BOX6_ION1_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION2_VALENCE")==0){
								BOX6_ION2_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION2_RADIUS")==0){
								BOX6_ION2_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION2_DIFF_COEFF")==0){
								BOX6_ION2_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION2_MASS")==0){
								BOX6_ION2_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION3_VALENCE")==0){
								BOX6_ION3_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION3_RADIUS")==0){
								BOX6_ION3_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION3_DIFF_COEFF")==0){
								BOX6_ION3_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION3_MASS")==0){
								BOX6_ION3_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION4_VALENCE")==0){
								BOX6_ION4_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION4_RADIUS")==0){
								BOX6_ION4_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION4_DIFF_COEFF")==0){
								BOX6_ION4_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION4_MASS")==0){
								BOX6_ION4_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION5_VALENCE")==0){
								BOX6_ION5_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION5_RADIUS")==0){
								BOX6_ION5_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION5_DIFF_COEFF")==0){
								BOX6_ION5_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION5_MASS")==0){
								BOX6_ION5_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION6_VALENCE")==0){
								BOX6_ION6_VALENCE=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION6_RADIUS")==0){
								BOX6_ION6_RADIUS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION6_DIFF_COEFF")==0){
								BOX6_ION6_DIFF_COEFF=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
							else if(getTokenbyNumber(buffer, "=", 1).compare("BOX6_ION6_MASS")==0){
								BOX6_ION6_MASS=atof(getTokenbyNumber(buffer, "=", 2).c_str());
							}
//=============================================================================== ION PARAMETERS end
							else {
								cerr << "The format of the file " << conf_file << " is not correct" << endl;
								cerr << "Line: " << buffer << endl;
								exit(2);
							}
						}
					}
				}					
			}
			fin.close();
		}
		else{
			cerr << "The file " << conf_file << " is not a regular .conf file. Please give a correct path to .conf file" << endl;
			exit(2);
		}
	}
	else{
		cerr << "The file " << conf_file << " does not exist. Please give a correct path to .conf file" <<endl;
		exit(2);
	}
	PRM.SIM_STEPS=PRM.SIM_STEPS+PRM.HISTORY_SIZE;
	if(PRM.SIM_TYPE.compare("BULK")==0){
		PRM.Z_MOUTH_LEFT=-double(100.00);
		PRM.Z_MOUTH_RIGHT=double(100.00);
	}
    	else if(PRM.SIM_TYPE.compare("MEMBRANE")==0){
		PRM.Z_MOUTH_LEFT=-PRM.MEMBRANE_WIDTH/double(2.00);
		PRM.Z_MOUTH_RIGHT=PRM.MEMBRANE_WIDTH/double(2.00);
	}
	else if(PRM.SIM_TYPE.compare("PORE")==0){
		PRM.Z_MOUTH_LEFT=-PRM.MEMBRANE_WIDTH/double(2.00);
		PRM.Z_MOUTH_RIGHT=PRM.MEMBRANE_WIDTH/double(2.00);
	}
	return;
}

void check_parameters(){
	if(PRM.SIM_TYPE.compare("BULK")!=0 && PRM.SIM_TYPE.compare("MEMBRANE")!=0 && PRM.SIM_TYPE.compare("PORE")!=0){
		cout << "WRONG PARAMETER!" <<endl;
		cout << "Parameter SIM_TYPE is set to: " <<PRM.SIM_TYPE<<endl;
		cout << PRM.SIM_TYPE << " is not a keyword (BULK, MEMBRANE, PORE)" <<endl;
		exit(5);
	}
	if(PRM.SIM_TYPE.compare("BULK")==0) {
		if(PRM.induced_charge != 0) {
			cout << "WRONG PARAMETER! induced_charge must be zero for BULK simulation" <<endl;
			exit(5);
		}
		if(PRM.mean_square_displ != 0) {
			cout << "WRONG PARAMETER! mean_square_displ must be zero for BULK simulation" <<endl;
			exit(5);
		}
	}
	if(PRM.PREFIX.compare("NULL")==0){
		cout << "WRONG PARAMETER! PREFIX must be set! e.g. PREFIX = sim_prefix" <<endl;
		exit(5);
	}
	if(PRM.CONTROL_CELL_WIDTH<0){
		cout << "WRONG PARAMETER! CONTROL_CELL_WIDTH must be >0 (picometers)." <<endl;
		cout << "CONTROL CELLS MANAGE ION CREATION. PROVIDE A CONTROL CELL LARGE ENOUGH!"<<endl;
		exit(5);
	}
	if(PRM.BATH_WIDTH<0){
		cout << "WRONG PARAMETER! BATH_WIDTH must be >=0 (picometers)." <<endl;
		exit(5);
	}
	if(PRM.SIM_TYPE.compare("BULK")==0 && PRM.MEMBRANE_WIDTH!=0.00){
		cout << "WRONG PARAMETER!" <<endl;
		cout << "Parameter SIM_TYPE is set to: " <<PRM.SIM_TYPE<<". You are attempting to simulate a bulk solution."<<endl;
		cout << "Keep SIM_TYPE=BULK and remove or set to 0 MEMBRANE_WIDTH parameter in the .conf file to simulate a bulk solution."<<endl;  
		cout << "Or keep MEMBRANE_WIDTH (picometers) parameter and change SIM_TYPE parameter to simulate an ion channel"<<endl;
		exit(5);
	}
	if(PRM.SIM_DOMAIN_WIDTH_X<=0){
		cout << "WRONG PARAMETER! SIM_DOMAIN_WIDTH_X must be >0 (picometers)." <<endl;
		exit(5);
	}
	if(PRM.SIM_DOMAIN_WIDTH_Y<=0){
		cout << "WRONG PARAMETER! SIM_DOMAIN_WIDTH_Y must be >0 (picometers)." <<endl;
		exit(5);
	}
	if(PRM.CONC_LEFT_KCL<0){
		cout << "WRONG PARAMETER! CONC_LEFT_KCL must be >=0 (Moles)." << endl;
		exit(5);
	}
	if(PRM.CONC_LEFT_NACL<0){
		cout << "WRONG PARAMETER! CONC_LEFT_NACL must be >=0 (Moles)." <<endl;
		exit(5);
	}
	if(PRM.CONC_LEFT_CACL2<0){
		cout << "WRONG PARAMETER! CONC_LEFT_CACL2 must be >=0 (Moles)." <<endl;
		exit(5);
	}
	if(PRM.CONC_LEFT_MGCL2<0){
		cout << "WRONG PARAMETER! CONC_LEFT_MGCL2 must be >=0 (Moles)." <<endl;
		exit(5);
	}
	if(PRM.CONC_RIGHT_KCL<0){
		cout << "WRONG PARAMETER! CONC_RIGHT_KCL must be >=0 (Moles)." <<endl;
		exit(5);
	}
	if(PRM.CONC_RIGHT_NACL<0){
		cout << "WRONG PARAMETER! CONC_RIGHT_NACL must be >=0 (Moles)." <<endl;
		exit(5);
	}
	if(PRM.CONC_RIGHT_CACL2<0){
		cout << "WRONG PARAMETER! CONC_RIGHT_CACL2 must be >=0 (Moles)." <<endl;
		exit(5);
	}
	if(PRM.CONC_RIGHT_MGCL2<0){
		cout << "WRONG PARAMETER! CONC_RIGHT_MGCL2 must be >=0 (Moles)." <<endl;
		exit(5);
	}
	if(PRM.SHORT_RANGE_EXP<0){
		cout << "WRONG PARAMETER! SHORT_RANGE_EXP must be >=0 (12 default for L-J)." <<endl;
		exit(5);
	}
	if(!(PRM.potential==0 || PRM.potential==1 || PRM.potential==2)){
		cout << "WRONG PARAMETER! STATISTICS - potential must be chosen in:" <<endl;
		cout << "\t0 - no potential statistics collected"<<endl;
		cout << "\t1 - 1-D potential statistics collected along channel axis"<<endl;
		cout << "\t2 - (1) + 2-D potential statistics collected (rotational symmetry)"<<endl;
		exit(5);
	}
	if(!(PRM.concentrations==0 || PRM.concentrations==2 || PRM.concentrations==3)){
		cout << "WRONG PARAMETER! STATISTICS - concentrations must be chosen in:" <<endl;
		cout << "0 - no concentrations statistics collected"<<endl;
		cout << "2 - 2-D concentrations statistics collected (rotational symmetry)"<<endl;
		cout << "3 - 3-D concentrations statistics collected"<<endl;
		exit(5);
	}
	if(PRM.PBC.compare("000")!=0 && PRM.PBC.compare("001")!=0 && PRM.PBC.compare("110")!=0 && PRM.PBC.compare("111")!=0){
		cout << "WRONG PARAMETER! PBC must be chosen in [000, 001, 110, 111]. The digits are for x, y and z respectively" <<endl;
		cout << "0: isolated box, no periodicity"<<endl;
		cout << "1: periodic boundary conditions, closest image convention"<<endl;
		exit(5);
	}
	if(PRM.ION_RECYCLING.compare("000")!=0 && PRM.ION_RECYCLING.compare("001")!=0 && PRM.ION_RECYCLING.compare("002")!=0 && PRM.ION_RECYCLING.compare("110")!=0 && PRM.ION_RECYCLING.compare("111")!=0 && PRM.ION_RECYCLING.compare("112")!=0 && PRM.ION_RECYCLING.compare("220")!=0 && PRM.ION_RECYCLING.compare("221")!=0 && PRM.ION_RECYCLING.compare("222")!=0){
		cout << "WRONG PARAMETER! ION_RECYCLING must be chosen in [000, 001, 002, 110, 111, 112, 220, 221, 222]. The digits are for x, y and z respectively" <<endl;
		cout << "0: reflective boundary, NO ion recycling"<<endl;
		cout << "1: open boundary, ion recycling"<<endl;
		cout << "2: open boundary, ion deletion"<<endl;
		exit(5);
	}
	
	if(PRM.SR_METHOD.compare("EXPONENTIAL") && PRM.SR_METHOD.compare("LJ") && PRM.SR_METHOD.compare("PMF") && PRM.SR_METHOD.compare("CS-LJ")){
		cout << "WRONG PARAMETER!" <<endl;
		cout << "Parameter SR_METHOD is set to: " <<PRM.SR_METHOD<<endl;
		cout << "Parameter SR_METHOD must be chosen in EXPONENTIAL, LJ, PMF or CS-LJ" <<endl;
		exit(5);
	}
	
	if(PRM.EPS_W<1){
		cout << "WRONG PARAMETER! EPS_W must be >=1. Use EPS_W=1 for vacuum." <<endl;
		exit(5);
	}
	if(PRM.EPS_MEM<1 && PRM.SIM_TYPE.compare("BULK")!=0){
		cout << "WRONG PARAMETER! EPS_MEM must be >=1. Use EPS_MEM=1 for vacuum." <<endl;
		exit(5);
	}
	if(PRM.DIFF_COEFF_IN_CHANNEL!=1){
		PRM.CONSTANT_DIFF_COEFF=false;
	}
	if(PRM.DELTA_T==0){
		cout << "WRONG PARAMETER! DELTA_T must have a positive value (femtoseconds)." <<endl;
		exit(5);
	}
	if(PRM.TEMPERATURE==0){
		cout << "WRONG PARAMETER! TEMPERATURE must have a positive value (Kelvin)." <<endl;
		exit(5);
	}
	if(PRM.PREP_STEPS<0){
		cout << "WRONG PARAMETER! PREP_STEPS must have a positive value." <<endl;
		exit(5);
	}
	if(PRM.SIM_STEPS<0){
		cout << "WRONG PARAMETER! SIM_STEPS must have a positive value." <<endl;
		exit(5);
	}
	if(PRM.STATS_DZ==0){
		cout << "WRONG PARAMETER! STATS_DZ must have a positive value (picometer)." <<endl;
		exit(5);
	}
	if(PRM.STATS_OUT_FREQ==0){
		cout << "WRONG PARAMETER! STATS_OUT_FREQ must have a positive value (number of steps)." <<endl;
		exit(5);
	}
	else{
		PRM.NUM_OF_SUB_DIV=round(sqrt(PRM.NUM_OF_SUB_DIV));
	}
	PRM.kT=BOLTZMANN_K*PRM.TEMPERATURE;
	return;
}

void print_PDB_filter_file(){
	string ions_file=PRM.PREFIX+".filter.pdb";
	char *tfn = new char[ions_file.length()+1];
	strcpy(tfn, ions_file.c_str());     
	ofstream fout(tfn, ios::app);
	fout << setprecision(2);
	int num_atom = 1;
	int ind_last_output;
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){ // all the possible ionic species
		if(PRM.ions_to_simulate.at(is)){ // any ion of this kind in the simulation ?
			Ion iti(is);
			
			if(iti.name.substr(0, 1).compare("J")){
			
				int real_ions=0;
				int fake_ions=5;			
				
				for(int i=0; i<ions.front().size(); i++){ 
					if (ions.front().at(i).kind == is) {
						//~ cout << ions.front().at(i).z << "\t" << 1e-12*PRM.Z_MOUTH_LEFT << "\t" << 1e-12*PRM.Z_MOUTH_RIGHT <<endl;
						
						if(ions.front().at(i).z >= 1e-12*PRM.Z_MOUTH_LEFT && ions.front().at(i).z <= 1e-12*PRM.Z_MOUTH_RIGHT) {
						
							fout<<"ATOM  "<<		//recname
								setw(5)<<num_atom	//serial
								<<" "			//space
								<<setw(4)<<ions.front().at(i).name	//atom
								<<" "				//altLoc
								<<setw(3)<<ions.front().at(i).name	//resName
								<<" "				//space
								<<" "				//cainID
								<<setw(5)<<num_atom	//Seqno
								<<"   "				//three spaces
								<<setw(8)<<1e10*ions.front().at(i).x	//atom X coordinate
								<<setw(8)<<1e10*ions.front().at(i).y	//atom Y coordinate
								<<setw(8)<<1e10*ions.front().at(i).z	//atom Z coordinate
								<<setw(6)<<i			//occupancy
								<<setw(6)<<i			//tempFactor
								<<endl;
							
							real_ions++;
							fake_ions--;
							
							num_atom++;
						}
					}
				}
				
				
				
				for(int i=0; i<fake_ions; i++){
					fout<<"ATOM  "<<			//recname
						setw(5)<<num_atom			//serial
						<<" "					//space
						<<setw(4)<<iti.name	//atom
						<<" "					//altLoc
						<<setw(3)<<iti.name	//resName
						<<" "					//space
						<<" "					//cainID
						<<setw(5)<<real_ions+i			//Seqno
						<<"   "					//three spaces
						<<setw(8)<<1e-2*PRM.MEMBRANE_WIDTH	//atom X coordinate
						<<setw(8)<<1e-2*PRM.MEMBRANE_WIDTH	//atom Y coordinate
						<<setw(8)<<1e-2*PRM.MEMBRANE_WIDTH	//atom Z coordinate
						<<setw(6)<<i				//occupancy
						<<setw(6)<<i				//tempFactor
						<<endl;
				}
			}
		}
	}
	
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){ // all the possible ionic species
		if(PRM.ions_to_simulate.at(is)){ // any ion of this kind in the simulation ?
			Ion iti(is);
			
			if(iti.name.substr(0, 1).compare("J")==0){
				for(int i=0; i<ions.front().size(); i++){ 
					if (ions.front().at(i).kind == is) {
						
						fout<<"ATOM  "<<		//recname
							setw(5)<<num_atom	//serial
							<<" "			//space
							<<setw(4)<<"J"	//atom
							<<" "				//altLoc
							<<setw(3)<<"J"	//resName
							<<" "				//space
							<<" "				//cainID
							<<setw(5)<<num_atom	//Seqno
							<<"   "				//three spaces
							<<setw(8)<<1e10*ions.front().at(i).x	//atom X coordinate
							<<setw(8)<<1e10*ions.front().at(i).y	//atom Y coordinate
							<<setw(8)<<1e10*ions.front().at(i).z	//atom Z coordinate
							<<setw(6)<<i			//occupancy
							<<setw(6)<<i			//tempFactor
							<<endl;
					}
				}
			}
		}
	}
	

	fout<<"END"<<endl;
	fout.close();
	return;
}

void print_ions_trajectories(){
	string ions_file=PRM.PREFIX+".trajectory.dat";
	char *tfn = new char[ions_file.length()+1];
	strcpy(tfn, ions_file.c_str());     
	ofstream fout(tfn, ios::app);
	fout << setprecision(5);
	fout << "// step = " << INDEX_LAST_STEP << endl;
	for(int i=0; i<NUM_OF_IONS_IN_STEP[INDEX_STAT_STEP]; i++){ 
		fout<<i<<"\t"
			<<IONS[INDEX_LAST_STEP][i].kind<<"\t"
			<<IONS[INDEX_LAST_STEP][i].x<<"\t"
			<<IONS[INDEX_LAST_STEP][i].y<<"\t"
			<<IONS[INDEX_LAST_STEP][i].z<<"\t"
			<<IONS[INDEX_LAST_STEP][i].x_prev<<"\t"
			<<IONS[INDEX_LAST_STEP][i].y_prev<<"\t"
			<<IONS[INDEX_LAST_STEP][i].z_prev<<"\t"
			<<endl;
	}
	fout.close();
	return;
}

void print_PDB_file(vector <Ion>& ions){
	
	string ions_file=PRM.PREFIX+".pdb";
	char *tfn = new char[ions_file.length()+1];
	strcpy(tfn, ions_file.c_str());     
	ofstream fout(tfn, ios::app);
	fout << setprecision(2);

	int num_atom = 1;
	int ind_last_output;
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){ // all the possible ionic species
		if(PRM.ions_to_simulate.at(is)){ // any ion of this kind in the simulation ?
			int num_ions_out = 0;
			for(int i=0; i<ions.size(); i++){ 
				if (ions.at(i).kind == is) {
					fout<<"ATOM  "<<		//recname
						setw(5)<<num_atom	//serial
						<<" "			//space
						<<setw(4)<<ions.at(i).name	//atom
						<<" "				//altLoc
						<<setw(3)<<ions.at(i).name	//resName
						<<" "				//space
						<<" "				//cainID
						<<setw(5)<<num_atom	//Seqno
						<<"   "				//three spaces
						<<setw(8)<<1e10*ions.at(i).x	//atom X coordinate
						<<setw(8)<<1e10*ions.at(i).y	//atom Y coordinate
						<<setw(8)<<1e10*ions.at(i).z	//atom Z coordinate
						<<setw(6)<<i			//occupancy
						<<setw(6)<<i			//tempFactor
						<<endl;
					num_ions_out += 1;
					ind_last_output = i;
					num_atom += 1;
				}
			}
			for(int i=num_ions_out; i<2000; i++){
				fout<<"ATOM  "<<			//recname
					setw(5)<<num_atom			//serial
					<<" "					//space
					<<setw(4)<<ions.at(ind_last_output).name	//atom
					<<" "					//altLoc
					<<setw(3)<<ions.at(ind_last_output).name	//resName
					<<" "					//space
					<<" "					//cainID
					<<setw(5)<<num_atom			//Seqno
					<<"   "					//three spaces
					<<setw(8)<<1e-2*PRM.SIM_DOMAIN_WIDTH_X	//atom X coordinate
					<<setw(8)<<1e-2*PRM.SIM_DOMAIN_WIDTH_Y	//atom Y coordinate
					<<setw(8)<<1e-2*PRM.SIM_DOMAIN_WIDTH_Z	//atom Z coordinate
					<<setw(6)<<i				//occupancy
					<<setw(6)<<i				//tempFactor
					<<endl;
			}
		}
	}

	fout<<"END"<<endl;
	fout.close();
	return;
}

void print_PDB_file_debug(vector <Ion>& ions, int index){
	string s;
	stringstream ss;
	ss << index;
	ss >> s;
	
	string ions_file=PRM.PREFIX+"."+s+".pdb";
	char *tfn = new char[ions_file.length()+1];
	strcpy(tfn, ions_file.c_str());     
	ofstream fout(tfn, ios::app);
	fout << setprecision(4);

	int num_atom = 1;
	int ind_last_output;
	for(int is=0; is<NUM_OF_IONIC_SPECIES; is++){ // all the possible ionic species
		if(PRM.ions_to_simulate.at(is)){ // any ion of this kind in the simulation ?
			int num_ions_out = 0;
			for(int i=0; i<ions.size(); i++){ 
				if (ions.at(i).kind == is && ions.at(i).z>-2e-9 && ions.at(i).z<2e-9) {
					fout<<"ATOM  "<<		//recname
						setw(5)<<num_atom	//serial
						<<" "			//space
						<<setw(4)<<ions.at(i).name	//atom
						<<" "				//altLoc
						<<setw(3)<<ions.at(i).name	//resName
						<<" "				//space
						<<" "				//cainID
						<<setw(5)<<num_atom	//Seqno
						<<"   "				//three spaces
						<<setw(8)<<1e10*ions.at(i).x	//atom X coordinate
						<<setw(8)<<1e10*ions.at(i).y	//atom Y coordinate
						<<setw(8)<<1e10*ions.at(i).z	//atom Z coordinate
						<<setw(6)<<i			//occupancy
						<<setw(6)<<i			//tempFactor
						<<endl;
					num_ions_out += 1;
					ind_last_output = i;
					num_atom += 1;
				}
			}
			for(int i=num_ions_out; i<10; i++){
				fout<<"ATOM  "<<			//recname
					setw(5)<<num_atom			//serial
					<<" "					//space
					<<setw(4)<<ions.at(ind_last_output).name	//atom
					<<" "					//altLoc
					<<setw(3)<<ions.at(ind_last_output).name	//resName
					<<" "					//space
					<<" "					//cainID
					<<setw(5)<<num_atom			//Seqno
					<<"   "					//three spaces
					<<setw(8)<<1e-2*PRM.SIM_DOMAIN_WIDTH_X	//atom X coordinate
					<<setw(8)<<1e-2*PRM.SIM_DOMAIN_WIDTH_Y	//atom Y coordinate
					<<setw(8)<<1e-2*PRM.SIM_DOMAIN_WIDTH_Z	//atom Z coordinate
					<<setw(6)<<i				//occupancy
					<<setw(6)<<i				//tempFactor
					<<endl;
			}
		}
	}

	fout<<"END"<<endl;
	fout.close();
	return;
}


void print_surfaces(){
	double total_area=0.00;
	
    string fl1 = PRM.PREFIX + ".surfaces";
    string fl2 = PRM.PREFIX + ".n1";
    string fl3 = PRM.PREFIX + ".n2";
    
    if(surfaces.empty()){
        cerr << "NO surfaces!!!"<<endl;
    }
    else{
    
	    string hhh="";
	    createFile(fl1);
	    char *tfn1 = new char[fl1.length()+1];
	    strcpy(tfn1, fl1.c_str());     
	    ofstream fout1(tfn1, ios::app);
	    createFile(fl2);
	    char *tfn2 = new char[fl2.length()+1];
	    strcpy(tfn2, fl2.c_str());     
	    ofstream fout2(tfn2, ios::app);
	    createFile(fl3);
	    char *tfn3 = new char[fl3.length()+1];
	    strcpy(tfn3, fl3.c_str());     
	    ofstream fout3(tfn3, ios::app);
	    
		for(int i=0; i<surfaces.size(); i++){	
			fout1 <<"# tile "<<i<<"\n";
			fout1 <<surfaces.at(i).pA.x << " "<<surfaces.at(i).pA.y << " "<<surfaces.at(i).pA.z << "\n";
			fout1 <<surfaces.at(i).pB.x << " "<<surfaces.at(i).pB.y << " "<<surfaces.at(i).pB.z << "\n";
			fout1 <<surfaces.at(i).pC.x << " "<<surfaces.at(i).pC.y << " "<<surfaces.at(i).pC.z << "\n";	
			if(surfaces.at(i).pD.x!=0 && surfaces.at(i).pD.y!=0 && surfaces.at(i).pD.z!=0){
				fout1 <<surfaces.at(i).pD.x << " "<<surfaces.at(i).pD.y << " "<<surfaces.at(i).pD.z << "\n";				
			}			
			fout1 <<surfaces.at(i).pA.x << " "<<surfaces.at(i).pA.y << " "<<surfaces.at(i).pA.z << "\n\n\n";
			
			fout2 <<"# tile "<<i<<"\n";
			fout2 <<surfaces.at(i).center.x+surfaces.at(i).normal[0]*100 << " "<<surfaces.at(i).center.y+surfaces.at(i).normal[1]*100 << " "<<surfaces.at(i).center.z+surfaces.at(i).normal[2]*100 << "\n\n\n";

			fout3 <<"# tile "<<i<<"\n";
			fout3 <<surfaces.at(i).center.x-surfaces.at(i).normal[0]*100 << " "<<surfaces.at(i).center.y-surfaces.at(i).normal[1]*100 << " "<<surfaces.at(i).center.z-surfaces.at(i).normal[2]*100 << "\n\n\n";
			
			total_area+=surfaces.at(i).area;
		}
		
		cout << "total surfaces area: " << total_area <<endl;
		
	    fout1.close();	
	    fout2.close();	
	    fout3.close();	
	}

	return;
}

void print_subSurfaces(){
	double total_area=0.00;
	
    string fl1 = PRM.PREFIX + ".sub_surfaces";
    string fl2 = PRM.PREFIX + ".sub_n1";
    string fl3 = PRM.PREFIX + ".sub_n2";

    if(subSurfaces.empty()){
        cerr << "NO subSurfaces!!!"<<endl;
    }
	else{    
    
	 
    string hhh="";
    createFile(fl1);
    char *tfn1 = new char[fl1.length()+1];
    strcpy(tfn1, fl1.c_str());     
    ofstream fout1(tfn1, ios::app);
    createFile(fl2);
    char *tfn2 = new char[fl2.length()+1];
    strcpy(tfn2, fl2.c_str());     
    ofstream fout2(tfn2, ios::app);
    createFile(fl3);
    char *tfn3 = new char[fl3.length()+1];
    strcpy(tfn3, fl3.c_str());     
    ofstream fout3(tfn3, ios::app);
    
	
	for(int i=0; i<subSurfaces.size(); i++){	
		for(int j=0; j<subSurfaces.at(i).size(); j++){	
			fout1 << "# tile "<<i<<"\n";
			fout1 << subSurfaces.at(i).at(j).pA.x << " "<<subSurfaces.at(i).at(j).pA.y << " "<<subSurfaces.at(i).at(j).pA.z << "\n";
			fout1 << subSurfaces.at(i).at(j).pB.x << " "<<subSurfaces.at(i).at(j).pB.y << " "<<subSurfaces.at(i).at(j).pB.z << "\n";
			fout1 << subSurfaces.at(i).at(j).pC.x << " "<<subSurfaces.at(i).at(j).pC.y << " "<<subSurfaces.at(i).at(j).pC.z << "\n";	
			if(subSurfaces.at(i).at(j).pD.x!=0 && subSurfaces.at(i).at(j).pD.y!=0 && subSurfaces.at(i).at(j).pD.z!=0){
				fout1 << subSurfaces.at(i).at(j).pD.x << " "<<subSurfaces.at(i).at(j).pD.y << " "<<subSurfaces.at(i).at(j).pD.z << "\n";	
			}						
			fout1 << subSurfaces.at(i).at(j).pA.x << " "<<subSurfaces.at(i).at(j).pA.y << " "<<subSurfaces.at(i).at(j).pA.z << "\n\n\n";


			fout2 << "# tile "<<i<<"\n";
			fout2 << subSurfaces.at(i).at(j).center.x+subSurfaces.at(i).at(j).normal[0]*100 << " "<<subSurfaces.at(i).at(j).center.y+subSurfaces.at(i).at(j).normal[1]*100<< " "<<subSurfaces.at(i).at(j).center.z+subSurfaces.at(i).at(j).normal[2]*100 << "\n\n\n";


			fout3 << "# tile "<<i<<"\n";
			fout3 << subSurfaces.at(i).at(j).center.x-subSurfaces.at(i).at(j).normal[0]*100 << " "<<subSurfaces.at(i).at(j).center.y-subSurfaces.at(i).at(j).normal[1]*100 << " "<<subSurfaces.at(i).at(j).center.z-subSurfaces.at(i).at(j).normal[2]*100 << "\n\n\n";

			total_area+=subSurfaces.at(i).at(j).area;
		}
	}
		
	cout << "total subsurfaces area: " << total_area <<endl;
    
    fout1.close();	
    fout2.close();	
    fout3.close();	
}
    
	return;
}


void output_vmd_sim_domain(ostream& stream){
	
	double lx,ly,lz,hx,hy,hz;
    int i;
    
    if(surfaces.empty()){
        
    }
	else{
		
		stream << "draw materials off\ndraw color 17"<< endl;
		if(!surfaces.empty()){
			
			for(i=0; i<surfaces.size(); i++){
				if(surfaces.at(i).pD.x!=0 && surfaces.at(i).pD.y!=0){
					stream << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pA.x) << " " <<double(0.01)*(surfaces.at(i).pA.y) << " " <<double(0.01)*(surfaces.at(i).pA.z) << "\" \"" <<
					double(0.01)*(surfaces.at(i).pB.x) << " " <<double(0.01)*(surfaces.at(i).pB.y) << " " <<double(0.01)*(surfaces.at(i).pB.z)<< "\" style dashed \n";
					stream << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pB.x) << " " <<double(0.01)*(surfaces.at(i).pB.y) << " " <<double(0.01)*(surfaces.at(i).pB.z) << "\" \"" <<
					double(0.01)*(surfaces.at(i).pC.x) << " " <<double(0.01)*(surfaces.at(i).pC.y) << " " <<double(0.01)*(surfaces.at(i).pC.z)<< "\" style dashed \n";
					stream << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pC.x) << " " <<double(0.01)*(surfaces.at(i).pC.y) << " " <<double(0.01)*(surfaces.at(i).pC.z) << "\" \"" <<
					double(0.01)*(surfaces.at(i).pD.x) << " " <<double(0.01)*(surfaces.at(i).pD.y) << " " <<double(0.01)*(surfaces.at(i).pD.z)<< "\" style dashed \n";
					stream << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pD.x) << " " <<double(0.01)*(surfaces.at(i).pD.y) << " " <<double(0.01)*(surfaces.at(i).pD.z) << "\" \"" <<
					double(0.01)*(surfaces.at(i).pA.x) << " " <<double(0.01)*(surfaces.at(i).pA.y) << " " <<double(0.01)*(surfaces.at(i).pA.z)<< "\" style dashed \n";
				}
				else{
					stream << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pA.x) << " " <<double(0.01)*(surfaces.at(i).pA.y) << " " <<double(0.01)*(surfaces.at(i).pA.z) << "\" \"" <<
					double(0.01)*(surfaces.at(i).pB.x) << " " <<double(0.01)*(surfaces.at(i).pB.y) << " " <<double(0.01)*(surfaces.at(i).pB.z)<< "\" style dashed \n";
					stream << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pB.x) << " " <<double(0.01)*(surfaces.at(i).pB.y) << " " <<double(0.01)*(surfaces.at(i).pB.z) << "\" \"" <<
					double(0.01)*(surfaces.at(i).pC.x) << " " <<double(0.01)*(surfaces.at(i).pC.y) << " " <<double(0.01)*(surfaces.at(i).pC.z)<< "\" style dashed \n";
					stream << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pC.x) << " " <<double(0.01)*(surfaces.at(i).pC.y) << " " <<double(0.01)*(surfaces.at(i).pC.z) << "\" \"" <<
					double(0.01)*(surfaces.at(i).pA.x) << " " <<double(0.01)*(surfaces.at(i).pA.y) << " " <<double(0.01)*(surfaces.at(i).pA.z)<< "\" style dashed \n";
				}
			}
		}
	}
		
	lx=double(0.01)*(PRM.MIN_X);
	ly=double(0.01)*(PRM.MIN_Y);
	lz=double(0.01)*(PRM.MIN_Z);
	hx=double(0.01)*(PRM.MAX_X);
	hy=double(0.01)*(PRM.MAX_Y);
	hz=double(0.01)*(PRM.MAX_Z);
	stream << "draw materials off\ndraw color 11"<< endl;
	stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	stream << "draw color 9"<< endl;
	
	lx=double(0.01)*(PRM.MIN_X);
	ly=double(0.01)*(PRM.MIN_Y);
	lz=double(0.01)*(PRM.MIN_Z+PRM.CONTROL_CELL_WIDTH);
	hx=double(0.01)*(PRM.MAX_X);
	hy=double(0.01)*(PRM.MAX_Y);
	hz=double(0.01)*(PRM.MAX_Z-PRM.CONTROL_CELL_WIDTH);

	stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
	stream << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
	stream << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";

	
	if(PRM.BOX1_MAX_Z-PRM.BOX1_MIN_Z>0){
		stream << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX1_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX1_MAX_Z);
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX2_MAX_Z-PRM.BOX2_MIN_Z>0){
		stream << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX2_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX2_MAX_Z);
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX3_MAX_Z-PRM.BOX3_MIN_Z>0){
		stream << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX3_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX3_MAX_Z);
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX4_MAX_Z-PRM.BOX4_MIN_Z>0){
		stream << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX4_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX4_MAX_Z);
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX5_MAX_Z-PRM.BOX5_MIN_Z>0){
		stream << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX5_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX5_MAX_Z);
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX6_MAX_Z-PRM.BOX6_MIN_Z>0){
		stream << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX6_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX6_MAX_Z);
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		stream << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(!membrane_charges.empty()){
		for(int i=0; i<membrane_charges.size(); i++){
			if(membrane_charges.at(i).charge>0){
				stream << "draw color 23" <<endl;
			}
			else{
				stream << "draw color 13" <<endl;
			}
			stream << "draw sphere {" << 1e10*membrane_charges.at(i).x << " " << 1e10*membrane_charges.at(i).y << " " << 1e10*membrane_charges.at(i).z << "} radius 0.5 \n";
		}
	}
	
	stream << endl << endl;
	
	
	return;
}


void print_vmd_channel_files(){
    double lx,ly,lz,hx,hy,hz;
    int i;

    string fileout2 = PRM.PREFIX + "_draw_sim_domain.tcl";
    createFile(fileout2);
    char *tfn2 = new char[fileout2.length()+1];
    strcpy(tfn2, fileout2.c_str());     
    ofstream fout2(tfn2);

    if(surfaces.empty()){
    }
    else{
	string fileout = PRM.PREFIX + "_draw_tiles.tcl";
	createFile(fileout);
	char *tfn = new char[fileout.length()+1];
	strcpy(tfn, fileout.c_str());     
	ofstream fout(tfn);
	
	
	fout << "draw materials off\ndraw color 17"<< endl;
	if(!surfaces.empty()){
		
		for(i=0; i<surfaces.size(); i++){
				fout << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pA.x) << " " <<double(0.01)*(surfaces.at(i).pA.y) << " " <<double(0.01)*(surfaces.at(i).pA.z) << "\" \"" <<
				double(0.01)*(surfaces.at(i).pB.x) << " " <<double(0.01)*(surfaces.at(i).pB.y) << " " <<double(0.01)*(surfaces.at(i).pB.z)<< "\" style dashed \n";
				fout << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pB.x) << " " <<double(0.01)*(surfaces.at(i).pB.y) << " " <<double(0.01)*(surfaces.at(i).pB.z) << "\" \"" <<
				double(0.01)*(surfaces.at(i).pC.x) << " " <<double(0.01)*(surfaces.at(i).pC.y) << " " <<double(0.01)*(surfaces.at(i).pC.z)<< "\" style dashed \n";
				fout << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pC.x) << " " <<double(0.01)*(surfaces.at(i).pC.y) << " " <<double(0.01)*(surfaces.at(i).pC.z) << "\" \"" <<
				double(0.01)*(surfaces.at(i).pD.x) << " " <<double(0.01)*(surfaces.at(i).pD.y) << " " <<double(0.01)*(surfaces.at(i).pD.z)<< "\" style dashed \n";
				fout << "draw line "<< "\"" <<double(0.01)*(surfaces.at(i).pD.x) << " " <<double(0.01)*(surfaces.at(i).pD.y) << " " <<double(0.01)*(surfaces.at(i).pD.z) << "\" \"" <<
				double(0.01)*(surfaces.at(i).pA.x) << " " <<double(0.01)*(surfaces.at(i).pA.y) << " " <<double(0.01)*(surfaces.at(i).pA.z)<< "\" style dashed \n";
		}
	}
	
	fout.close();
    }
	lx=double(0.01)*(PRM.MIN_X);
	ly=double(0.01)*(PRM.MIN_Y);
	lz=double(0.01)*(PRM.MIN_Z);
	hx=double(0.01)*(PRM.MAX_X);
	hy=double(0.01)*(PRM.MAX_Y);
	hz=double(0.01)*(PRM.MAX_Z);
	fout2 << "draw materials off\ndraw color 11"<< endl;
	fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	fout2 << "draw color 9"<< endl;
	
	lx=double(0.01)*(PRM.MIN_X);
	ly=double(0.01)*(PRM.MIN_Y);
	lz=double(0.01)*(PRM.MIN_Z+PRM.CONTROL_CELL_WIDTH);
	hx=double(0.01)*(PRM.MAX_X);
	hy=double(0.01)*(PRM.MAX_Y);
	hz=double(0.01)*(PRM.MAX_Z-PRM.CONTROL_CELL_WIDTH);

	fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
	fout2 << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";

	
	if(PRM.BOX1_MAX_Z-PRM.BOX1_MIN_Z>0){
		fout2 << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX1_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX1_MAX_Z);
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX2_MAX_Z-PRM.BOX2_MIN_Z>0){
		fout2 << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX2_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX2_MAX_Z);
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX3_MAX_Z-PRM.BOX3_MIN_Z>0){
		fout2 << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX3_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX3_MAX_Z);
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX4_MAX_Z-PRM.BOX4_MIN_Z>0){
		fout2 << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX4_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX4_MAX_Z);
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX5_MAX_Z-PRM.BOX5_MIN_Z>0){
		fout2 << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX5_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX5_MAX_Z);
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	
	if(PRM.BOX6_MAX_Z-PRM.BOX6_MIN_Z>0){
		fout2 << "draw color 25"<< endl;
		lx = 0.01*(-1000);
		ly = 0.01*(-1000);
		lz = 0.01*(PRM.BOX6_MIN_Z);
		hx = 0.01*(1000);
		hy = 0.01*(1000);
		hz = 0.01*(PRM.BOX6_MAX_Z);
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << ly << " " << lz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << lz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << lz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << lz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << lz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << ly << " " << hz << "\" \"" << hx << " " << ly << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << ly << " " << hz << "\" \"" << hx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << hx << " " << hy << " " << hz << "\" \"" << lx << " " << hy << " " << hz << "\" style solid \n";
		fout2 << "draw line \"" << lx << " " << hy << " " << hz << "\" \"" << lx << " " << ly << " " << hz << "\" style solid \n";
	}
	fout2.close();
	 
	string fileout3 = PRM.PREFIX + "_mem_charges.pdbrq";
	createFile(fileout3);
	if(!membrane_charges.empty()){
		char *tfn3 = new char[fileout3.length()+1];
		strcpy(tfn3, fileout3.c_str());     
		ofstream fout3(tfn3);
		fout3 << setprecision(4);
		for(int i=0; i<membrane_charges.size(); i++){
			//~ cout << membrane_charges.at(i).charge <<endl;
			fout3<<"ATOM  "<<				//recname
			setw(5)<<i+1					//serial
			<<" ";								//space
			if(membrane_charges.at(i).charge<0){
				fout3<<setw(4)<<"   C"		//atom
				<<" "								//altLoc
				<<setw(3)<<"  C";		//resName
			}
			else{
				fout3<<setw(4)<<"   N"		//atom
				<<" "								//altLoc
				<<setw(3)<<"  N";		//resName
			}
			fout3<<" "					//space
			<<" "						//chainID
			<<setw(5)<<i+1					//Seqno
			<<"   "							//three spaces
			<<setw(8)<<1e10*membrane_charges.at(i).x	//atom X coordinate
			<<setw(8)<<1e10*membrane_charges.at(i).y	//atom Y coordinate
			<<setw(8)<<1e10*membrane_charges.at(i).z	//atom Z coordinate
			<<setw(6)<<i				//occupancy
			<<setw(6)<<i				//tempFactor
			<<endl;
		}
		fout3<<"END"<<endl;
		fout3.close();	
	}
    return;
}



void print_conf_file(ostream& stream){	
	
	string fileout = PRM.PREFIX + ".conf";
	createFile(fileout);
	char *tfn = new char[fileout.length()+1];
	strcpy(tfn, fileout.c_str());     
	ofstream fout(tfn);
	
	stream << endl << "########################################" <<endl;
	stream << "    SIMULATION PARAMETERS";
	stream << endl << "########################################" <<endl;
	fout << endl << "########################################" <<endl;
	fout << "    SIMULATION PARAMETERS";
	fout << endl << "########################################" <<endl;
	
//============== PAGE 2 - SIM. DOMAIN	
	stream << "\nSIM_TYPE = " << " " << PRM.SIM_TYPE << endl;  
	stream << "PREFIX = " << " " << PRM.PREFIX << endl <<endl;  
	
	stream << "SIM_DOMAIN_WIDTH_X = " << " " << PRM.SIM_DOMAIN_WIDTH_X << endl;  
	stream << "SIM_DOMAIN_WIDTH_Y = " << " " << PRM.SIM_DOMAIN_WIDTH_Y << endl;  
	stream << "CONTROL_CELL_WIDTH = " << "  " << PRM.CONTROL_CELL_WIDTH << endl;  
	stream << "BATH_WIDTH = " << "  " << PRM.BATH_WIDTH << endl;  
	if(PRM.SIM_TYPE.compare("BULK")!=0){
		stream << "MEMBRANE_WIDTH = " << "   " << PRM.MEMBRANE_WIDTH << endl;  
	}
	if(PRM.SIM_TYPE.compare("PORE")==0){
		stream << "LEFT_VESTIBULE_CURVATURE_RADIUS = " << " " << PRM.LEFT_VESTIBULE_CURVATURE_RADIUS << endl;  
		stream << "LEFT_VESTIBULE_MIN_CHANNEL_RADIUS = " << " " << PRM.LEFT_VESTIBULE_MIN_CHANNEL_RADIUS << endl;  
		stream << "RIGHT_VESTIBULE_CURVATURE_RADIUS = " << " " << PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS << endl;  
		stream << "RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS = " << " " << PRM.RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS << endl;  
	}
	
	fout << "\nSIM_TYPE = " << " " << PRM.SIM_TYPE << endl;
	
	fout << "PREFIX = " << " " << PRM.PREFIX << endl <<endl;  
	
	fout << "SIM_DOMAIN_WIDTH_X = " << " " << PRM.SIM_DOMAIN_WIDTH_X << endl;  
	fout << "SIM_DOMAIN_WIDTH_Y = " << " " << PRM.SIM_DOMAIN_WIDTH_Y << endl;  
	fout << "CONTROL_CELL_WIDTH = " << "  " << PRM.CONTROL_CELL_WIDTH << endl;  
	fout << "BATH_WIDTH = " << "  " << PRM.BATH_WIDTH << endl;  
	if(PRM.SIM_TYPE.compare("BULK")!=0){
		fout << "MEMBRANE_WIDTH = " << "   " << PRM.MEMBRANE_WIDTH << endl;  
	}
	if(PRM.SIM_TYPE.compare("PORE")==0){
		
		fout << "LEFT_VESTIBULE_CURVATURE_RADIUS = " << " " << PRM.LEFT_VESTIBULE_CURVATURE_RADIUS << endl;  
		fout << "LEFT_VESTIBULE_MIN_CHANNEL_RADIUS = " << " " << PRM.LEFT_VESTIBULE_MIN_CHANNEL_RADIUS << endl;  
		fout << "RIGHT_VESTIBULE_CURVATURE_RADIUS = " << " " << PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS << endl;  
		fout << "RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS = " << " " << PRM.RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS << endl;  

	}
	
		
//============== PAGE 3 - MOTION	
	stream << endl << "DELTA_T =   "<< 1e15*PRM.DELTA_T << endl;
	stream << "SIM_STEPS =   "<< PRM.SIM_STEPS << endl;
	stream << "PREP_STEPS =   "<< PRM.PREP_STEPS << endl;
	
	if(PRM.SIM_TYPE.compare("PORE")==0){
		stream << "DIFF_COEFF_IN_CHANNEL =   "<< PRM.DIFF_COEFF_IN_CHANNEL << endl;
	}
	stream << "SR_METHOD =   "<< PRM.SR_METHOD << endl;
	if(PRM.SR_METHOD.compare("EXPONENTIAL")==0){
		stream << "SHORT_RANGE_EXP =   "<< PRM.SHORT_RANGE_EXP << endl;
	}
	stream << "TEMPERATURE =   "<< PRM.TEMPERATURE << endl;
	stream << "SEED =   "<< PRM.SEED << endl;
	
	fout << endl << "DELTA_T =   "<< 1e15*PRM.DELTA_T << endl;
	fout << "SIM_STEPS =   "<< PRM.SIM_STEPS << endl;
	fout << "PREP_STEPS =   "<< PRM.PREP_STEPS << endl;
	
	if(PRM.SIM_TYPE.compare("PORE")==0){
		fout << "DIFF_COEFF_IN_CHANNEL =   "<< PRM.DIFF_COEFF_IN_CHANNEL << endl;
	}
	fout << "SR_METHOD =   "<< PRM.SR_METHOD << endl;
	if(PRM.SR_METHOD.compare("EXPONENTIAL")==0){
		fout << "SHORT_RANGE_EXP =   "<< PRM.SHORT_RANGE_EXP << endl;
	}
	fout << "TEMPERATURE =   "<< PRM.TEMPERATURE << endl;
	fout << "SEED =   "<< PRM.SEED << endl;
	
	
	
	
	
//============== PAGE 4 - ELECTROSTATICS	
	stream << endl << "APPLIED_POTENTIAL =   "<< PRM.APPLIED_POTENTIAL << endl;
	stream << "EPS_W =   "<< PRM.EPS_W << endl;
	if(PRM.SIM_TYPE.compare("BULK")!=0){
		stream << "EPS_MEM =   "<< PRM.EPS_MEM << endl;	
		stream << "TILES_PER_RING =   "<< PRM.TILES_PER_RING << endl;
		if(PRM.SIM_TYPE.compare("PORE")==0){
			stream << "NUM_OF_DIV =   "<< PRM.NUM_OF_DIV << endl;
		}
		stream << "NUM_OF_SUB_DIV =   "<< PRM.NUM_OF_SUB_DIV << endl;
	}
	
	fout << endl << "APPLIED_POTENTIAL =   "<< PRM.APPLIED_POTENTIAL << endl;
	fout << "EPS_W =   "<< PRM.EPS_W << endl;
	if(PRM.SIM_TYPE.compare("BULK")!=0){
		fout << "EPS_MEM =   "<< PRM.EPS_MEM << endl;
		fout << "TILES_PER_RING =   "<< PRM.TILES_PER_RING << endl;
		if(PRM.SIM_TYPE.compare("PORE")==0){
			fout << "NUM_OF_DIV =   "<< PRM.NUM_OF_DIV << endl;
		}
		fout << "NUM_OF_SUB_DIV =   "<< PRM.NUM_OF_SUB_DIV << endl;
	}
		
	
	
//============== PAGE 5 - ION CONCENTRATIONS	
	stream << endl << "CONC_LEFT_KCL =   "<< PRM.CONC_LEFT_KCL << endl;	
	stream << "CONC_LEFT_NACL =   "<< PRM.CONC_LEFT_NACL << endl;	
	stream << "CONC_LEFT_CACL2 =   "<< PRM.CONC_LEFT_CACL2 << endl;	
	stream << "CONC_LEFT_MGCL2 =   "<< PRM.CONC_LEFT_MGCL2 << endl;	
	stream << "CONC_RIGHT_KCL =   "<< PRM.CONC_RIGHT_KCL << endl;	
	stream << "CONC_RIGHT_NACL =   "<< PRM.CONC_RIGHT_NACL << endl;	
	stream << "CONC_RIGHT_CACL2 =   "<< PRM.CONC_RIGHT_CACL2 << endl;	
	stream << "CONC_RIGHT_MGCL2 =   "<< PRM.CONC_RIGHT_MGCL2 << endl;	

	stream << "PERIODIC BOUNDARY CONDITIONS = "<< PRM.PBC << endl;
	stream << "ION RECYCLING = "<< PRM.ION_RECYCLING << endl;	
	
	fout << endl << "CONC_LEFT_KCL =   "<< PRM.CONC_LEFT_KCL << endl;	
	fout << "CONC_LEFT_NACL =   "<< PRM.CONC_LEFT_NACL << endl;	
	fout << "CONC_LEFT_CACL2 =   "<< PRM.CONC_LEFT_CACL2 << endl;	
	fout << "CONC_LEFT_MGCL2 =   "<< PRM.CONC_LEFT_MGCL2 << endl;	
	fout << "CONC_RIGHT_KCL =   "<< PRM.CONC_RIGHT_KCL << endl;	
	fout << "CONC_RIGHT_NACL =   "<< PRM.CONC_RIGHT_NACL << endl;	
	fout << "CONC_RIGHT_CACL2 =   "<< PRM.CONC_RIGHT_CACL2 << endl;	
	fout << "CONC_RIGHT_MGCL2 =   "<< PRM.CONC_RIGHT_MGCL2 << endl;	
	
	fout << "PBC = "<< PRM.PBC << endl;
	fout << "ION_RECYCLING = "<< PRM.ION_RECYCLING << endl;	
	
	
	
//============== PAGE 6 - IONS	
	stream << endl;
	if(PRM.CONC_LEFT_KCL>0 || PRM.CONC_RIGHT_KCL>0){
		stream << "K_VALENCE = " <<  K_VALENCE <<endl;
		stream << "K_RADIUS = " <<  K_RADIUS <<endl;
		stream << "K_DIFF_COEFF = " <<  K_DIFF_COEFF <<endl;
		stream << "K_MASS = " <<  K_MASS <<endl;
	}
	if(PRM.CONC_LEFT_NACL>0 || PRM.CONC_RIGHT_NACL>0){
		stream << "NA_VALENCE = " <<  NA_VALENCE <<endl;
		stream << "NA_RADIUS = " <<  NA_RADIUS <<endl;
		stream << "NA_DIFF_COEFF = " <<  NA_DIFF_COEFF <<endl;
		stream << "NA_MASS = " <<  NA_MASS <<endl;
	}
	if(PRM.CONC_LEFT_CACL2>0 || PRM.CONC_RIGHT_CACL2>0){
		stream << "CA_VALENCE = " <<  CA_VALENCE <<endl;
		stream << "CA_RADIUS = " <<  CA_RADIUS <<endl;
		stream << "CA_DIFF_COEFF = " <<  CA_DIFF_COEFF <<endl;
		stream << "CA_MASS = " <<  CA_MASS <<endl;
	}
	if(PRM.CONC_LEFT_MGCL2>0 || PRM.CONC_RIGHT_MGCL2>0){
		stream << "MG_VALENCE = " <<  MG_VALENCE <<endl;
		stream << "MG_RADIUS = " <<  MG_RADIUS <<endl;
		stream << "MG_DIFF_COEFF = " <<  MG_DIFF_COEFF <<endl;
		stream << "MG_MASS = " <<  MG_MASS <<endl;
	}
	if(PRM.CONC_LEFT_KCL>0 || PRM.CONC_RIGHT_KCL>0 || PRM.CONC_LEFT_NACL>0 || PRM.CONC_RIGHT_NACL>0 || PRM.CONC_LEFT_CACL2>0 || PRM.CONC_RIGHT_CACL2>0 || PRM.CONC_LEFT_MGCL2>0 || PRM.CONC_RIGHT_MGCL2>0){
		stream << "CL_VALENCE = " <<  CL_VALENCE <<endl;
		stream << "CL_RADIUS = " <<  CL_RADIUS <<endl;
		stream << "CL_DIFF_COEFF = " <<  CL_DIFF_COEFF <<endl;
		stream << "CL_MASS = " <<  CL_MASS <<endl;
	}
	
	fout << endl;
	if(PRM.CONC_LEFT_KCL>0 || PRM.CONC_RIGHT_KCL>0){
		fout << "K_VALENCE = " <<  K_VALENCE <<endl;
		fout << "K_RADIUS = " <<  K_RADIUS <<endl;
		fout << "K_DIFF_COEFF = " <<  K_DIFF_COEFF <<endl;
		fout << "K_MASS = " <<  K_MASS <<endl;
	}
	if(PRM.CONC_LEFT_NACL>0 || PRM.CONC_RIGHT_NACL>0){
		fout << "NA_VALENCE = " <<  NA_VALENCE <<endl;
		fout << "NA_RADIUS = " <<  NA_RADIUS <<endl;
		fout << "NA_DIFF_COEFF = " <<  NA_DIFF_COEFF <<endl;
		fout << "NA_MASS = " <<  NA_MASS <<endl;
	}
	if(PRM.CONC_LEFT_CACL2>0 || PRM.CONC_RIGHT_CACL2>0){
		fout << "CA_VALENCE = " <<  CA_VALENCE <<endl;
		fout << "CA_RADIUS = " <<  CA_RADIUS <<endl;
		fout << "CA_DIFF_COEFF = " <<  CA_DIFF_COEFF <<endl;
		fout << "CA_MASS = " <<  CA_MASS <<endl;
	}
	if(PRM.CONC_LEFT_MGCL2>0 || PRM.CONC_RIGHT_MGCL2>0){
		fout << "MG_VALENCE = " <<  MG_VALENCE <<endl;
		fout << "MG_RADIUS = " <<  MG_RADIUS <<endl;
		fout << "MG_DIFF_COEFF = " <<  MG_DIFF_COEFF <<endl;
		fout << "MG_MASS = " <<  MG_MASS <<endl;
	}
	if(PRM.CONC_LEFT_KCL>0 || PRM.CONC_RIGHT_KCL>0 || PRM.CONC_LEFT_NACL>0 || PRM.CONC_RIGHT_NACL>0 || PRM.CONC_LEFT_CACL2>0 || PRM.CONC_RIGHT_CACL2>0 || PRM.CONC_LEFT_MGCL2>0 || PRM.CONC_RIGHT_MGCL2>0){
		fout << "CL_VALENCE = " <<  CL_VALENCE <<endl;
		fout << "CL_RADIUS = " <<  CL_RADIUS <<endl;
		fout << "CL_DIFF_COEFF = " <<  CL_DIFF_COEFF <<endl;
		fout << "CL_MASS = " <<  CL_MASS <<endl;
	}
	
	
	
//============== PAGE 7 - ION BOXES
	if(PRM.SIM_TYPE.compare("PORE")==0){
		if(PRM.BOX1_MIN_Z<PRM.BOX1_MAX_Z){
			stream << endl;
			stream << "BOX1_MIN_Z =  "<< PRM.BOX1_MIN_Z << endl;
			stream << "BOX1_MAX_Z =  "<< PRM.BOX1_MAX_Z << endl;
			if(PRM.BOX1_N1>0){	
				stream << "BOX1_IS1 = J11" << endl;
				stream << "BOX1_N1 = "<< PRM.BOX1_N1 << endl;
			}
			if(PRM.BOX1_N2>0){	
				stream << "BOX1_IS2 = J12" << endl;
				stream << "BOX1_N2 = "<< PRM.BOX1_N2 << endl;
			}
			if(PRM.BOX1_N3>0){	
				stream << "BOX1_IS3 = J13" << endl;
				stream << "BOX1_N3 = "<< PRM.BOX1_N3 << endl;
			}
			if(PRM.BOX1_N4>0){	
				stream << "BOX1_IS4 = J14" << endl;
				stream << "BOX1_N4 = "<< PRM.BOX1_N4 << endl;
			}
			if(PRM.BOX1_N5>0){	
				stream << "BOX1_IS5 = J15" << endl;
				stream << "BOX1_N5 = "<< PRM.BOX1_N5 << endl;
			}
			if(PRM.BOX1_N6>0){	
				stream << "BOX1_IS6 = J16" << endl;
				stream << "BOX1_N6 = "<< PRM.BOX1_N6 << endl;
			}
		}	
		if(PRM.BOX2_MIN_Z<PRM.BOX2_MAX_Z){
			stream << endl;
			stream << "BOX2_MIN_Z =  "<< PRM.BOX2_MIN_Z << endl;
			stream << "BOX2_MAX_Z =  "<< PRM.BOX2_MAX_Z << endl;
			if(PRM.BOX2_N1>0){	
				stream << "BOX2_IS1 = J21" << endl;
				stream << "BOX2_N1 = "<< PRM.BOX2_N1 << endl;
			}
			if(PRM.BOX2_N2>0){	
				stream << "BOX2_IS2 = J22" << endl;
				stream << "BOX2_N2 = "<< PRM.BOX2_N2 << endl;
			}
			if(PRM.BOX2_N3>0){	
				stream << "BOX2_IS3 = J23" << endl;
				stream << "BOX2_N3 = "<< PRM.BOX2_N3 << endl;
			}
			if(PRM.BOX2_N4>0){	
				stream << "BOX2_IS4 = J24" << endl;
				stream << "BOX2_N4 = "<< PRM.BOX2_N4 << endl;
			}
			if(PRM.BOX2_N5>0){	
				stream << "BOX2_IS5 = J25" << endl;
				stream << "BOX2_N5 = "<< PRM.BOX2_N5 << endl;
			}
			if(PRM.BOX2_N6>0){	
				stream << "BOX2_IS6 = J26" << endl;
				stream << "BOX2_N6 = "<< PRM.BOX2_N6 << endl;
			}
		}
		if(PRM.BOX3_MIN_Z<PRM.BOX3_MAX_Z){
			stream << endl;
			stream << "BOX3_MIN_Z =  "<< PRM.BOX3_MIN_Z << endl;
			stream << "BOX3_MAX_Z =  "<< PRM.BOX3_MAX_Z << endl;
			if(PRM.BOX3_N1>0){	
				stream << "BOX3_IS1 = J31" << endl;
				stream << "BOX3_N1 = "<< PRM.BOX3_N1 << endl;
			}
			if(PRM.BOX3_N2>0){	
				stream << "BOX3_IS2 = J32" << endl;
				stream << "BOX3_N2 = "<< PRM.BOX3_N2 << endl;
			}
			if(PRM.BOX3_N3>0){	
				stream << "BOX3_IS3 = J33" << endl;
				stream << "BOX3_N3 = "<< PRM.BOX3_N3 << endl;
			}
			if(PRM.BOX3_N4>0){	
				stream << "BOX3_IS4 = J34" << endl;
				stream << "BOX3_N4 = "<< PRM.BOX3_N4 << endl;
			}
			if(PRM.BOX3_N5>0){	
				stream << "BOX3_IS5 = J35" << endl;
				stream << "BOX3_N5 = "<< PRM.BOX3_N5 << endl;
			}
			if(PRM.BOX3_N6>0){	
				stream << "BOX3_IS6 = J36" << endl;
				stream << "BOX3_N6 = "<< PRM.BOX3_N6 << endl;
			}
		}
		if(PRM.BOX4_MIN_Z<PRM.BOX4_MAX_Z){
			stream << endl;
			stream << "BOX4_MIN_Z =  "<< PRM.BOX4_MIN_Z << endl;
			stream << "BOX4_MAX_Z =  "<< PRM.BOX4_MAX_Z << endl;
			if(PRM.BOX4_N1>0){	
				stream << "BOX4_IS1 = J41" << endl;
				stream << "BOX4_N1 = "<< PRM.BOX4_N1 << endl;
			}
			if(PRM.BOX4_N2>0){	
				stream << "BOX4_IS2 = J42" << endl;
				stream << "BOX4_N2 = "<< PRM.BOX4_N2 << endl;
			}
			if(PRM.BOX4_N3>0){	
				stream << "BOX4_IS3 = J43" << endl;
				stream << "BOX4_N3 = "<< PRM.BOX4_N3 << endl;
			}
			if(PRM.BOX4_N4>0){	
				stream << "BOX4_IS4 = J44" << endl;
				stream << "BOX4_N4 = "<< PRM.BOX4_N4 << endl;
			}
			if(PRM.BOX4_N5>0){	
				stream << "BOX4_IS5 = J45" << endl;
				stream << "BOX4_N5 = "<< PRM.BOX4_N5 << endl;
			}
			if(PRM.BOX4_N6>0){	
				stream << "BOX4_IS6 = J46" << endl;
				stream << "BOX4_N6 = "<< PRM.BOX4_N6 << endl;
			}
		}
		if(PRM.BOX5_MIN_Z<PRM.BOX5_MAX_Z){
			stream << endl;
			stream << "BOX5_MIN_Z =  "<< PRM.BOX5_MIN_Z << endl;
			stream << "BOX5_MAX_Z =  "<< PRM.BOX5_MAX_Z << endl;
			if(PRM.BOX5_N1>0){	
				stream << "BOX5_IS1 = J51" << endl;
				stream << "BOX5_N1 = "<< PRM.BOX5_N1 << endl;
			}
			if(PRM.BOX5_N2>0){	
				stream << "BOX5_IS2 = J52" << endl;
				stream << "BOX5_N2 = "<< PRM.BOX5_N2 << endl;
			}
			if(PRM.BOX5_N3>0){	
				stream << "BOX5_IS3 = J53" << endl;
				stream << "BOX5_N3 = "<< PRM.BOX5_N3 << endl;
			}
			if(PRM.BOX5_N4>0){	
				stream << "BOX5_IS4 = J54" << endl;
				stream << "BOX5_N4 = "<< PRM.BOX5_N4 << endl;
			}
			if(PRM.BOX5_N5>0){	
				stream << "BOX5_IS5 = J55" << endl;
				stream << "BOX5_N5 = "<< PRM.BOX5_N5 << endl;
			}
			if(PRM.BOX5_N6>0){	
				stream << "BOX5_IS6 = J56" << endl;
				stream << "BOX5_N6 = "<< PRM.BOX5_N6 << endl;
			}
		}
		if(PRM.BOX6_MIN_Z<PRM.BOX6_MAX_Z){
			stream << endl;
			stream << "BOX6_MIN_Z =  "<< PRM.BOX6_MIN_Z << endl;
			stream << "BOX6_MAX_Z =  "<< PRM.BOX6_MAX_Z << endl;
			if(PRM.BOX6_N1>0){	
				stream << "BOX6_IS1 = J61" << endl;
				stream << "BOX6_N1 = "<< PRM.BOX6_N1 << endl;
			}
			if(PRM.BOX6_N2>0){	
				stream << "BOX6_IS2 = J62" << endl;
				stream << "BOX6_N2 = "<< PRM.BOX6_N2 << endl;
			}
			if(PRM.BOX6_N3>0){	
				stream << "BOX6_IS3 = J63" << endl;
				stream << "BOX6_N3 = "<< PRM.BOX6_N3 << endl;
			}
			if(PRM.BOX6_N4>0){	
				stream << "BOX6_IS4 = J64" << endl;
				stream << "BOX6_N4 = "<< PRM.BOX6_N4 << endl;
			}
			if(PRM.BOX6_N5>0){	
				stream << "BOX6_IS5 = J65" << endl;
				stream << "BOX6_N5 = "<< PRM.BOX6_N5 << endl;
			}
			if(PRM.BOX6_N6>0){	
				stream << "BOX6_IS6 = J66" << endl;
				stream << "BOX6_N6 = "<< PRM.BOX6_N6 << endl;
			}
		}
	}
	
	if(PRM.SIM_TYPE.compare("PORE")==0){
		if(PRM.BOX1_MIN_Z<PRM.BOX1_MAX_Z){
			fout << endl;
			fout << "BOX1_MIN_Z =  "<< PRM.BOX1_MIN_Z << endl;
			fout << "BOX1_MAX_Z =  "<< PRM.BOX1_MAX_Z << endl;
			if(PRM.BOX1_N1>0){	
				fout << "BOX1_IS1 = J11" << endl;
				fout << "BOX1_N1 = "<< PRM.BOX1_N1 << endl;
			}
			if(PRM.BOX1_N2>0){	
				fout << "BOX1_IS2 = J12" << endl;
				fout << "BOX1_N2 = "<< PRM.BOX1_N2 << endl;
			}
			if(PRM.BOX1_N3>0){	
				fout << "BOX1_IS3 = J13" << endl;
				fout << "BOX1_N3 = "<< PRM.BOX1_N3 << endl;
			}
			if(PRM.BOX1_N4>0){	
				fout << "BOX1_IS4 = J14" << endl;
				fout << "BOX1_N4 = "<< PRM.BOX1_N4 << endl;
			}
			if(PRM.BOX1_N5>0){	
				fout << "BOX1_IS5 = J15" << endl;
				fout << "BOX1_N5 = "<< PRM.BOX1_N5 << endl;
			}
			if(PRM.BOX1_N6>0){	
				fout << "BOX1_IS6 = J16" << endl;
				fout << "BOX1_N6 = "<< PRM.BOX1_N6 << endl;
			}
		}	
		if(PRM.BOX2_MIN_Z<PRM.BOX2_MAX_Z){
			fout << endl;
			fout << "BOX2_MIN_Z =  "<< PRM.BOX2_MIN_Z << endl;
			fout << "BOX2_MAX_Z =  "<< PRM.BOX2_MAX_Z << endl;
			if(PRM.BOX2_N1>0){	
				fout << "BOX2_IS1 = J21" << endl;
				fout << "BOX2_N1 = "<< PRM.BOX2_N1 << endl;
			}
			if(PRM.BOX2_N2>0){	
				fout << "BOX2_IS2 = J22" << endl;
				fout << "BOX2_N2 = "<< PRM.BOX2_N2 << endl;
			}
			if(PRM.BOX2_N3>0){	
				fout << "BOX2_IS3 = J23" << endl;
				fout << "BOX2_N3 = "<< PRM.BOX2_N3 << endl;
			}
			if(PRM.BOX2_N4>0){	
				fout << "BOX2_IS4 = J24" << endl;
				fout << "BOX2_N4 = "<< PRM.BOX2_N4 << endl;
			}
			if(PRM.BOX2_N5>0){	
				fout << "BOX2_IS5 = J25" << endl;
				fout << "BOX2_N5 = "<< PRM.BOX2_N5 << endl;
			}
			if(PRM.BOX2_N6>0){	
				fout << "BOX2_IS6 = J26" << endl;
				fout << "BOX2_N6 = "<< PRM.BOX2_N6 << endl;
			}
		}
		if(PRM.BOX3_MIN_Z<PRM.BOX3_MAX_Z){
			fout << endl;
			fout << "BOX3_MIN_Z =  "<< PRM.BOX3_MIN_Z << endl;
			fout << "BOX3_MAX_Z =  "<< PRM.BOX3_MAX_Z << endl;
			if(PRM.BOX3_N1>0){	
				fout << "BOX3_IS1 = J31" << endl;
				fout << "BOX3_N1 = "<< PRM.BOX3_N1 << endl;
			}
			if(PRM.BOX3_N2>0){	
				fout << "BOX3_IS2 = J32" << endl;
				fout << "BOX3_N2 = "<< PRM.BOX3_N2 << endl;
			}
			if(PRM.BOX3_N3>0){	
				fout << "BOX3_IS3 = J33" << endl;
				fout << "BOX3_N3 = "<< PRM.BOX3_N3 << endl;
			}
			if(PRM.BOX3_N4>0){	
				fout << "BOX3_IS4 = J34" << endl;
				fout << "BOX3_N4 = "<< PRM.BOX3_N4 << endl;
			}
			if(PRM.BOX3_N5>0){	
				fout << "BOX3_IS5 = J35" << endl;
				fout << "BOX3_N5 = "<< PRM.BOX3_N5 << endl;
			}
			if(PRM.BOX3_N6>0){	
				fout << "BOX3_IS6 = J36" << endl;
				fout << "BOX3_N6 = "<< PRM.BOX3_N6 << endl;
			}
		}
		if(PRM.BOX4_MIN_Z<PRM.BOX4_MAX_Z){
			fout << endl;
			fout << "BOX4_MIN_Z =  "<< PRM.BOX4_MIN_Z << endl;
			fout << "BOX4_MAX_Z =  "<< PRM.BOX4_MAX_Z << endl;
			if(PRM.BOX4_N1>0){	
				fout << "BOX4_IS1 = J41" << endl;
				fout << "BOX4_N1 = "<< PRM.BOX4_N1 << endl;
			}
			if(PRM.BOX4_N2>0){	
				fout << "BOX4_IS2 = J42" << endl;
				fout << "BOX4_N2 = "<< PRM.BOX4_N2 << endl;
			}
			if(PRM.BOX4_N3>0){	
				fout << "BOX4_IS3 = J43" << endl;
				fout << "BOX4_N3 = "<< PRM.BOX4_N3 << endl;
			}
			if(PRM.BOX4_N4>0){	
				fout << "BOX4_IS4 = J44" << endl;
				fout << "BOX4_N4 = "<< PRM.BOX4_N4 << endl;
			}
			if(PRM.BOX4_N5>0){	
				fout << "BOX4_IS5 = J45" << endl;
				fout << "BOX4_N5 = "<< PRM.BOX4_N5 << endl;
			}
			if(PRM.BOX4_N6>0){	
				fout << "BOX4_IS6 = J46" << endl;
				fout << "BOX4_N6 = "<< PRM.BOX4_N6 << endl;
			}
		}
		if(PRM.BOX5_MIN_Z<PRM.BOX5_MAX_Z){
			fout << endl;
			fout << "BOX5_MIN_Z =  "<< PRM.BOX5_MIN_Z << endl;
			fout << "BOX5_MAX_Z =  "<< PRM.BOX5_MAX_Z << endl;
			if(PRM.BOX5_N1>0){	
				fout << "BOX5_IS1 = J51" << endl;
				fout << "BOX5_N1 = "<< PRM.BOX5_N1 << endl;
			}
			if(PRM.BOX5_N2>0){	
				fout << "BOX5_IS2 = J52" << endl;
				fout << "BOX5_N2 = "<< PRM.BOX5_N2 << endl;
			}
			if(PRM.BOX5_N3>0){	
				fout << "BOX5_IS3 = J53" << endl;
				fout << "BOX5_N3 = "<< PRM.BOX5_N3 << endl;
			}
			if(PRM.BOX5_N4>0){	
				fout << "BOX5_IS4 = J54" << endl;
				fout << "BOX5_N4 = "<< PRM.BOX5_N4 << endl;
			}
			if(PRM.BOX5_N5>0){	
				fout << "BOX5_IS5 = J55" << endl;
				fout << "BOX5_N5 = "<< PRM.BOX5_N5 << endl;
			}
			if(PRM.BOX5_N6>0){	
				fout << "BOX5_IS6 = J56" << endl;
				fout << "BOX5_N6 = "<< PRM.BOX5_N6 << endl;
			}
		}
		if(PRM.BOX6_MIN_Z<PRM.BOX6_MAX_Z){
			fout << endl;
			fout << "BOX6_MIN_Z =  "<< PRM.BOX6_MIN_Z << endl;
			fout << "BOX6_MAX_Z =  "<< PRM.BOX6_MAX_Z << endl;
			if(PRM.BOX6_N1>0){	
				fout << "BOX6_IS1 = J61" << endl;
				fout << "BOX6_N1 = "<< PRM.BOX6_N1 << endl;
			}
			if(PRM.BOX6_N2>0){	
				fout << "BOX6_IS2 = J62" << endl;
				fout << "BOX6_N2 = "<< PRM.BOX6_N2 << endl;
			}
			if(PRM.BOX6_N3>0){	
				fout << "BOX6_IS3 = J63" << endl;
				fout << "BOX6_N3 = "<< PRM.BOX6_N3 << endl;
			}
			if(PRM.BOX6_N4>0){	
				fout << "BOX6_IS4 = J64" << endl;
				fout << "BOX6_N4 = "<< PRM.BOX6_N4 << endl;
			}
			if(PRM.BOX6_N5>0){	
				fout << "BOX6_IS5 = J65" << endl;
				fout << "BOX6_N5 = "<< PRM.BOX6_N5 << endl;
			}
			if(PRM.BOX6_N6>0){	
				fout << "BOX6_IS6 = J66" << endl;
				fout << "BOX6_N6 = "<< PRM.BOX6_N6 << endl;
			}
		}
	}
	
	if(PRM.SIM_TYPE.compare("PORE")==0){
		if(PRM.BOX1_MIN_Z<PRM.BOX1_MAX_Z){
			stream <<endl;
			if(PRM.BOX1_N1>0){
				stream << "BOX1_ION1_VALENCE = " <<  BOX1_ION1_VALENCE <<endl;
				stream << "BOX1_ION1_RADIUS = " <<  BOX1_ION1_RADIUS <<endl;
				stream << "BOX1_ION1_DIFF_COEFF = " <<  BOX1_ION1_DIFF_COEFF <<endl;
				stream << "BOX1_ION1_MASS = " <<  BOX1_ION1_MASS <<endl;
			}
			if(PRM.BOX1_N2>0){
				stream << "BOX1_ION2_VALENCE = " <<  BOX1_ION2_VALENCE <<endl;
				stream << "BOX1_ION2_RADIUS = " <<  BOX1_ION2_RADIUS <<endl;
				stream << "BOX1_ION2_DIFF_COEFF = " <<  BOX1_ION2_DIFF_COEFF <<endl;
				stream << "BOX1_ION2_MASS = " <<  BOX1_ION2_MASS <<endl;
			}
			if(PRM.BOX1_N3>0){
				stream << "BOX1_ION3_VALENCE = " <<  BOX1_ION3_VALENCE <<endl;
				stream << "BOX1_ION3_RADIUS = " <<  BOX1_ION3_RADIUS <<endl;
				stream << "BOX1_ION3_DIFF_COEFF = " <<  BOX1_ION3_DIFF_COEFF <<endl;
				stream << "BOX1_ION3_MASS = " <<  BOX1_ION3_MASS <<endl;
			}
			if(PRM.BOX1_N4>0){
				stream << "BOX1_ION4_VALENCE = " <<  BOX1_ION4_VALENCE <<endl;
				stream << "BOX1_ION4_RADIUS = " <<  BOX1_ION4_RADIUS <<endl;
				stream << "BOX1_ION4_DIFF_COEFF = " <<  BOX1_ION4_DIFF_COEFF <<endl;
				stream << "BOX1_ION4_MASS = " <<  BOX1_ION4_MASS <<endl;
			}
			if(PRM.BOX1_N5>0){
				stream << "BOX1_ION5_VALENCE = " <<  BOX1_ION5_VALENCE <<endl;
				stream << "BOX1_ION5_RADIUS = " <<  BOX1_ION5_RADIUS <<endl;
				stream << "BOX1_ION5_DIFF_COEFF = " <<  BOX1_ION5_DIFF_COEFF <<endl;
				stream << "BOX1_ION5_MASS = " <<  BOX1_ION5_MASS <<endl;
			}
			if(PRM.BOX1_N6>0){
				stream << "BOX1_ION6_VALENCE = " <<  BOX1_ION6_VALENCE <<endl;
				stream << "BOX1_ION6_RADIUS = " <<  BOX1_ION6_RADIUS <<endl;
				stream << "BOX1_ION6_DIFF_COEFF = " <<  BOX1_ION6_DIFF_COEFF <<endl;
				stream << "BOX1_ION6_MASS = " <<  BOX1_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX2_MIN_Z<PRM.BOX2_MAX_Z){
			stream <<endl;
			if(PRM.BOX2_N1>0){
				stream << "BOX2_ION1_VALENCE = " <<  BOX2_ION1_VALENCE <<endl;
				stream << "BOX2_ION1_RADIUS = " <<  BOX2_ION1_RADIUS <<endl;
				stream << "BOX2_ION1_DIFF_COEFF = " <<  BOX2_ION1_DIFF_COEFF <<endl;
				stream << "BOX2_ION1_MASS = " <<  BOX2_ION1_MASS <<endl;
			}
			if(PRM.BOX2_N2>0){
				stream << "BOX2_ION2_VALENCE = " <<  BOX2_ION2_VALENCE <<endl;
				stream << "BOX2_ION2_RADIUS = " <<  BOX2_ION2_RADIUS <<endl;
				stream << "BOX2_ION2_DIFF_COEFF = " <<  BOX2_ION2_DIFF_COEFF <<endl;
				stream << "BOX2_ION2_MASS = " <<  BOX2_ION2_MASS <<endl;
			}
			if(PRM.BOX2_N3>0){
				stream << "BOX2_ION3_VALENCE = " <<  BOX2_ION3_VALENCE <<endl;
				stream << "BOX2_ION3_RADIUS = " <<  BOX2_ION3_RADIUS <<endl;
				stream << "BOX2_ION3_DIFF_COEFF = " <<  BOX2_ION3_DIFF_COEFF <<endl;
				stream << "BOX2_ION3_MASS = " <<  BOX2_ION3_MASS <<endl;
			}
			if(PRM.BOX2_N4>0){
				stream << "BOX2_ION4_VALENCE = " <<  BOX2_ION4_VALENCE <<endl;
				stream << "BOX2_ION4_RADIUS = " <<  BOX2_ION4_RADIUS <<endl;
				stream << "BOX2_ION4_DIFF_COEFF = " <<  BOX2_ION4_DIFF_COEFF <<endl;
				stream << "BOX2_ION4_MASS = " <<  BOX2_ION4_MASS <<endl;
			}
			if(PRM.BOX2_N5>0){
				stream << "BOX2_ION5_VALENCE = " <<  BOX2_ION5_VALENCE <<endl;
				stream << "BOX2_ION5_RADIUS = " <<  BOX2_ION5_RADIUS <<endl;
				stream << "BOX2_ION5_DIFF_COEFF = " <<  BOX2_ION5_DIFF_COEFF <<endl;
				stream << "BOX2_ION5_MASS = " <<  BOX2_ION5_MASS <<endl;
			}
			if(PRM.BOX2_N6>0){
				stream << "BOX2_ION6_VALENCE = " <<  BOX2_ION6_VALENCE <<endl;
				stream << "BOX2_ION6_RADIUS = " <<  BOX2_ION6_RADIUS <<endl;
				stream << "BOX2_ION6_DIFF_COEFF = " <<  BOX2_ION6_DIFF_COEFF <<endl;
				stream << "BOX2_ION6_MASS = " <<  BOX2_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX3_MIN_Z<PRM.BOX3_MAX_Z){
			stream <<endl;
			if(PRM.BOX3_N1>0){
				stream << "BOX3_ION1_VALENCE = " <<  BOX3_ION1_VALENCE <<endl;
				stream << "BOX3_ION1_RADIUS = " <<  BOX3_ION1_RADIUS <<endl;
				stream << "BOX3_ION1_DIFF_COEFF = " <<  BOX3_ION1_DIFF_COEFF <<endl;
				stream << "BOX3_ION1_MASS = " <<  BOX3_ION1_MASS <<endl;
			}
			if(PRM.BOX3_N2>0){
				stream << "BOX3_ION2_VALENCE = " <<  BOX3_ION2_VALENCE <<endl;
				stream << "BOX3_ION2_RADIUS = " <<  BOX3_ION2_RADIUS <<endl;
				stream << "BOX3_ION2_DIFF_COEFF = " <<  BOX3_ION2_DIFF_COEFF <<endl;
				stream << "BOX3_ION2_MASS = " <<  BOX3_ION2_MASS <<endl;
			}
			if(PRM.BOX3_N3>0){
				stream << "BOX3_ION3_VALENCE = " <<  BOX3_ION3_VALENCE <<endl;
				stream << "BOX3_ION3_RADIUS = " <<  BOX3_ION3_RADIUS <<endl;
				stream << "BOX3_ION3_DIFF_COEFF = " <<  BOX3_ION3_DIFF_COEFF <<endl;
				stream << "BOX3_ION3_MASS = " <<  BOX3_ION3_MASS <<endl;
			}
			if(PRM.BOX3_N4>0){
				stream << "BOX3_ION4_VALENCE = " <<  BOX3_ION4_VALENCE <<endl;
				stream << "BOX3_ION4_RADIUS = " <<  BOX3_ION4_RADIUS <<endl;
				stream << "BOX3_ION4_DIFF_COEFF = " <<  BOX3_ION4_DIFF_COEFF <<endl;
				stream << "BOX3_ION4_MASS = " <<  BOX3_ION4_MASS <<endl;
			}
			if(PRM.BOX3_N5>0){
				stream << "BOX3_ION5_VALENCE = " <<  BOX3_ION5_VALENCE <<endl;
				stream << "BOX3_ION5_RADIUS = " <<  BOX3_ION5_RADIUS <<endl;
				stream << "BOX3_ION5_DIFF_COEFF = " <<  BOX3_ION5_DIFF_COEFF <<endl;
				stream << "BOX3_ION5_MASS = " <<  BOX3_ION5_MASS <<endl;
			}
			if(PRM.BOX3_N6>0){
				stream << "BOX3_ION6_VALENCE = " <<  BOX3_ION6_VALENCE <<endl;
				stream << "BOX3_ION6_RADIUS = " <<  BOX3_ION6_RADIUS <<endl;
				stream << "BOX3_ION6_DIFF_COEFF = " <<  BOX3_ION6_DIFF_COEFF <<endl;
				stream << "BOX3_ION6_MASS = " <<  BOX3_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX4_MIN_Z<PRM.BOX4_MAX_Z){
			stream <<endl;
			if(PRM.BOX4_N1>0){
				stream << "BOX4_ION1_VALENCE = " <<  BOX4_ION1_VALENCE <<endl;
				stream << "BOX4_ION1_RADIUS = " <<  BOX4_ION1_RADIUS <<endl;
				stream << "BOX4_ION1_DIFF_COEFF = " <<  BOX4_ION1_DIFF_COEFF <<endl;
				stream << "BOX4_ION1_MASS = " <<  BOX4_ION1_MASS <<endl;
			}
			if(PRM.BOX4_N2>0){
				stream << "BOX4_ION2_VALENCE = " <<  BOX4_ION2_VALENCE <<endl;
				stream << "BOX4_ION2_RADIUS = " <<  BOX4_ION2_RADIUS <<endl;
				stream << "BOX4_ION2_DIFF_COEFF = " <<  BOX4_ION2_DIFF_COEFF <<endl;
				stream << "BOX4_ION2_MASS = " <<  BOX4_ION2_MASS <<endl;
			}
			if(PRM.BOX4_N3>0){
				stream << "BOX4_ION3_VALENCE = " <<  BOX4_ION3_VALENCE <<endl;
				stream << "BOX4_ION3_RADIUS = " <<  BOX4_ION3_RADIUS <<endl;
				stream << "BOX4_ION3_DIFF_COEFF = " <<  BOX4_ION3_DIFF_COEFF <<endl;
				stream << "BOX4_ION3_MASS = " <<  BOX4_ION3_MASS <<endl;
			}
			if(PRM.BOX4_N4>0){
				stream << "BOX4_ION4_VALENCE = " <<  BOX4_ION4_VALENCE <<endl;
				stream << "BOX4_ION4_RADIUS = " <<  BOX4_ION4_RADIUS <<endl;
				stream << "BOX4_ION4_DIFF_COEFF = " <<  BOX4_ION4_DIFF_COEFF <<endl;
				stream << "BOX4_ION4_MASS = " <<  BOX4_ION4_MASS <<endl;
			}
			if(PRM.BOX4_N5>0){
				stream << "BOX4_ION5_VALENCE = " <<  BOX4_ION5_VALENCE <<endl;
				stream << "BOX4_ION5_RADIUS = " <<  BOX4_ION5_RADIUS <<endl;
				stream << "BOX4_ION5_DIFF_COEFF = " <<  BOX4_ION5_DIFF_COEFF <<endl;
				stream << "BOX4_ION5_MASS = " <<  BOX4_ION5_MASS <<endl;
			}
			if(PRM.BOX4_N6>0){
				stream << "BOX4_ION6_VALENCE = " <<  BOX4_ION6_VALENCE <<endl;
				stream << "BOX4_ION6_RADIUS = " <<  BOX4_ION6_RADIUS <<endl;
				stream << "BOX4_ION6_DIFF_COEFF = " <<  BOX4_ION6_DIFF_COEFF <<endl;
				stream << "BOX4_ION6_MASS = " <<  BOX4_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX5_MIN_Z<PRM.BOX5_MAX_Z){
			stream <<endl;
			if(PRM.BOX5_N1>0){
				stream << "BOX5_ION1_VALENCE = " <<  BOX5_ION1_VALENCE <<endl;
				stream << "BOX5_ION1_RADIUS = " <<  BOX5_ION1_RADIUS <<endl;
				stream << "BOX5_ION1_DIFF_COEFF = " <<  BOX5_ION1_DIFF_COEFF <<endl;
				stream << "BOX5_ION1_MASS = " <<  BOX5_ION1_MASS <<endl;
			}
			if(PRM.BOX5_N2>0){
				stream << "BOX5_ION2_VALENCE = " <<  BOX5_ION2_VALENCE <<endl;
				stream << "BOX5_ION2_RADIUS = " <<  BOX5_ION2_RADIUS <<endl;
				stream << "BOX5_ION2_DIFF_COEFF = " <<  BOX5_ION2_DIFF_COEFF <<endl;
				stream << "BOX5_ION2_MASS = " <<  BOX5_ION2_MASS <<endl;
			}
			if(PRM.BOX5_N3>0){
				stream << "BOX5_ION3_VALENCE = " <<  BOX5_ION3_VALENCE <<endl;
				stream << "BOX5_ION3_RADIUS = " <<  BOX5_ION3_RADIUS <<endl;
				stream << "BOX5_ION3_DIFF_COEFF = " <<  BOX5_ION3_DIFF_COEFF <<endl;
				stream << "BOX5_ION3_MASS = " <<  BOX5_ION3_MASS <<endl;
			}
			if(PRM.BOX5_N4>0){
				stream << "BOX5_ION4_VALENCE = " <<  BOX5_ION4_VALENCE <<endl;
				stream << "BOX5_ION4_RADIUS = " <<  BOX5_ION4_RADIUS <<endl;
				stream << "BOX5_ION4_DIFF_COEFF = " <<  BOX5_ION4_DIFF_COEFF <<endl;
				stream << "BOX5_ION4_MASS = " <<  BOX5_ION4_MASS <<endl;
			}
			if(PRM.BOX5_N5>0){
				stream << "BOX5_ION5_VALENCE = " <<  BOX5_ION5_VALENCE <<endl;
				stream << "BOX5_ION5_RADIUS = " <<  BOX5_ION5_RADIUS <<endl;
				stream << "BOX5_ION5_DIFF_COEFF = " <<  BOX5_ION5_DIFF_COEFF <<endl;
				stream << "BOX5_ION5_MASS = " <<  BOX5_ION5_MASS <<endl;
			}
			if(PRM.BOX5_N6>0){
				stream << "BOX5_ION6_VALENCE = " <<  BOX5_ION6_VALENCE <<endl;
				stream << "BOX5_ION6_RADIUS = " <<  BOX5_ION6_RADIUS <<endl;
				stream << "BOX5_ION6_DIFF_COEFF = " <<  BOX5_ION6_DIFF_COEFF <<endl;
				stream << "BOX5_ION6_MASS = " <<  BOX5_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX6_MIN_Z<PRM.BOX6_MAX_Z){
			stream <<endl;
			if(PRM.BOX6_N1>0){
				stream << "BOX6_ION1_VALENCE = " <<  BOX6_ION1_VALENCE <<endl;
				stream << "BOX6_ION1_RADIUS = " <<  BOX6_ION1_RADIUS <<endl;
				stream << "BOX6_ION1_DIFF_COEFF = " <<  BOX6_ION1_DIFF_COEFF <<endl;
				stream << "BOX6_ION1_MASS = " <<  BOX6_ION1_MASS <<endl;
			}
			if(PRM.BOX6_N2>0){
				stream << "BOX6_ION2_VALENCE = " <<  BOX6_ION2_VALENCE <<endl;
				stream << "BOX6_ION2_RADIUS = " <<  BOX6_ION2_RADIUS <<endl;
				stream << "BOX6_ION2_DIFF_COEFF = " <<  BOX6_ION2_DIFF_COEFF <<endl;
				stream << "BOX6_ION2_MASS = " <<  BOX6_ION2_MASS <<endl;
			}
			if(PRM.BOX6_N3>0){
				stream << "BOX6_ION3_VALENCE = " <<  BOX6_ION3_VALENCE <<endl;
				stream << "BOX6_ION3_RADIUS = " <<  BOX6_ION3_RADIUS <<endl;
				stream << "BOX6_ION3_DIFF_COEFF = " <<  BOX6_ION3_DIFF_COEFF <<endl;
				stream << "BOX6_ION3_MASS = " <<  BOX6_ION3_MASS <<endl;
			}
			if(PRM.BOX6_N4>0){
				stream << "BOX6_ION4_VALENCE = " <<  BOX6_ION4_VALENCE <<endl;
				stream << "BOX6_ION4_RADIUS = " <<  BOX6_ION4_RADIUS <<endl;
				stream << "BOX6_ION4_DIFF_COEFF = " <<  BOX6_ION4_DIFF_COEFF <<endl;
				stream << "BOX6_ION4_MASS = " <<  BOX6_ION4_MASS <<endl;
			}
			if(PRM.BOX6_N5>0){
				stream << "BOX6_ION5_VALENCE = " <<  BOX6_ION5_VALENCE <<endl;
				stream << "BOX6_ION5_RADIUS = " <<  BOX6_ION5_RADIUS <<endl;
				stream << "BOX6_ION5_DIFF_COEFF = " <<  BOX6_ION5_DIFF_COEFF <<endl;
				stream << "BOX6_ION5_MASS = " <<  BOX6_ION5_MASS <<endl;
			}
			if(PRM.BOX6_N6>0){
				stream << "BOX6_ION6_VALENCE = " <<  BOX6_ION6_VALENCE <<endl;
				stream << "BOX6_ION6_RADIUS = " <<  BOX6_ION6_RADIUS <<endl;
				stream << "BOX6_ION6_DIFF_COEFF = " <<  BOX6_ION6_DIFF_COEFF <<endl;
				stream << "BOX6_ION6_MASS = " <<  BOX6_ION6_MASS <<endl;
			}
		}
	}
	
	if(PRM.SIM_TYPE.compare("PORE")==0){	
		if(PRM.BOX1_MIN_Z<PRM.BOX1_MAX_Z){
			fout <<endl;
			if(PRM.BOX1_N1>0){
				fout << "BOX1_ION1_VALENCE = " <<  BOX1_ION1_VALENCE <<endl;
				fout << "BOX1_ION1_RADIUS = " <<  BOX1_ION1_RADIUS <<endl;
				fout << "BOX1_ION1_DIFF_COEFF = " <<  BOX1_ION1_DIFF_COEFF <<endl;
				fout << "BOX1_ION1_MASS = " <<  BOX1_ION1_MASS <<endl;
			}
			if(PRM.BOX1_N2>0){
				fout << "BOX1_ION2_VALENCE = " <<  BOX1_ION2_VALENCE <<endl;
				fout << "BOX1_ION2_RADIUS = " <<  BOX1_ION2_RADIUS <<endl;
				fout << "BOX1_ION2_DIFF_COEFF = " <<  BOX1_ION2_DIFF_COEFF <<endl;
				fout << "BOX1_ION2_MASS = " <<  BOX1_ION2_MASS <<endl;
			}
			if(PRM.BOX1_N3>0){
				fout << "BOX1_ION3_VALENCE = " <<  BOX1_ION3_VALENCE <<endl;
				fout << "BOX1_ION3_RADIUS = " <<  BOX1_ION3_RADIUS <<endl;
				fout << "BOX1_ION3_DIFF_COEFF = " <<  BOX1_ION3_DIFF_COEFF <<endl;
				fout << "BOX1_ION3_MASS = " <<  BOX1_ION3_MASS <<endl;
			}
			if(PRM.BOX1_N4>0){
				fout << "BOX1_ION4_VALENCE = " <<  BOX1_ION4_VALENCE <<endl;
				fout << "BOX1_ION4_RADIUS = " <<  BOX1_ION4_RADIUS <<endl;
				fout << "BOX1_ION4_DIFF_COEFF = " <<  BOX1_ION4_DIFF_COEFF <<endl;
				fout << "BOX1_ION4_MASS = " <<  BOX1_ION4_MASS <<endl;
			}
			if(PRM.BOX1_N5>0){
				fout << "BOX1_ION5_VALENCE = " <<  BOX1_ION5_VALENCE <<endl;
				fout << "BOX1_ION5_RADIUS = " <<  BOX1_ION5_RADIUS <<endl;
				fout << "BOX1_ION5_DIFF_COEFF = " <<  BOX1_ION5_DIFF_COEFF <<endl;
				fout << "BOX1_ION5_MASS = " <<  BOX1_ION5_MASS <<endl;
			}
			if(PRM.BOX1_N6>0){
				fout << "BOX1_ION6_VALENCE = " <<  BOX1_ION6_VALENCE <<endl;
				fout << "BOX1_ION6_RADIUS = " <<  BOX1_ION6_RADIUS <<endl;
				fout << "BOX1_ION6_DIFF_COEFF = " <<  BOX1_ION6_DIFF_COEFF <<endl;
				fout << "BOX1_ION6_MASS = " <<  BOX1_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX2_MIN_Z<PRM.BOX2_MAX_Z){
			fout <<endl;
			if(PRM.BOX2_N1>0){
				fout << "BOX2_ION1_VALENCE = " <<  BOX2_ION1_VALENCE <<endl;
				fout << "BOX2_ION1_RADIUS = " <<  BOX2_ION1_RADIUS <<endl;
				fout << "BOX2_ION1_DIFF_COEFF = " <<  BOX2_ION1_DIFF_COEFF <<endl;
				fout << "BOX2_ION1_MASS = " <<  BOX2_ION1_MASS <<endl;
			}
			if(PRM.BOX2_N2>0){
				fout << "BOX2_ION2_VALENCE = " <<  BOX2_ION2_VALENCE <<endl;
				fout << "BOX2_ION2_RADIUS = " <<  BOX2_ION2_RADIUS <<endl;
				fout << "BOX2_ION2_DIFF_COEFF = " <<  BOX2_ION2_DIFF_COEFF <<endl;
				fout << "BOX2_ION2_MASS = " <<  BOX2_ION2_MASS <<endl;
			}
			if(PRM.BOX2_N3>0){
				fout << "BOX2_ION3_VALENCE = " <<  BOX2_ION3_VALENCE <<endl;
				fout << "BOX2_ION3_RADIUS = " <<  BOX2_ION3_RADIUS <<endl;
				fout << "BOX2_ION3_DIFF_COEFF = " <<  BOX2_ION3_DIFF_COEFF <<endl;
				fout << "BOX2_ION3_MASS = " <<  BOX2_ION3_MASS <<endl;
			}
			if(PRM.BOX2_N4>0){
				fout << "BOX2_ION4_VALENCE = " <<  BOX2_ION4_VALENCE <<endl;
				fout << "BOX2_ION4_RADIUS = " <<  BOX2_ION4_RADIUS <<endl;
				fout << "BOX2_ION4_DIFF_COEFF = " <<  BOX2_ION4_DIFF_COEFF <<endl;
				fout << "BOX2_ION4_MASS = " <<  BOX2_ION4_MASS <<endl;
			}
			if(PRM.BOX2_N5>0){
				fout << "BOX2_ION5_VALENCE = " <<  BOX2_ION5_VALENCE <<endl;
				fout << "BOX2_ION5_RADIUS = " <<  BOX2_ION5_RADIUS <<endl;
				fout << "BOX2_ION5_DIFF_COEFF = " <<  BOX2_ION5_DIFF_COEFF <<endl;
				fout << "BOX2_ION5_MASS = " <<  BOX2_ION5_MASS <<endl;
			}
			if(PRM.BOX2_N6>0){
				fout << "BOX2_ION6_VALENCE = " <<  BOX2_ION6_VALENCE <<endl;
				fout << "BOX2_ION6_RADIUS = " <<  BOX2_ION6_RADIUS <<endl;
				fout << "BOX2_ION6_DIFF_COEFF = " <<  BOX2_ION6_DIFF_COEFF <<endl;
				fout << "BOX2_ION6_MASS = " <<  BOX2_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX3_MIN_Z<PRM.BOX3_MAX_Z){
			fout <<endl;
			if(PRM.BOX3_N1>0){
				fout << "BOX3_ION1_VALENCE = " <<  BOX3_ION1_VALENCE <<endl;
				fout << "BOX3_ION1_RADIUS = " <<  BOX3_ION1_RADIUS <<endl;
				fout << "BOX3_ION1_DIFF_COEFF = " <<  BOX3_ION1_DIFF_COEFF <<endl;
				fout << "BOX3_ION1_MASS = " <<  BOX3_ION1_MASS <<endl;
			}
			if(PRM.BOX3_N2>0){
				fout << "BOX3_ION2_VALENCE = " <<  BOX3_ION2_VALENCE <<endl;
				fout << "BOX3_ION2_RADIUS = " <<  BOX3_ION2_RADIUS <<endl;
				fout << "BOX3_ION2_DIFF_COEFF = " <<  BOX3_ION2_DIFF_COEFF <<endl;
				fout << "BOX3_ION2_MASS = " <<  BOX3_ION2_MASS <<endl;
			}
			if(PRM.BOX3_N3>0){
				fout << "BOX3_ION3_VALENCE = " <<  BOX3_ION3_VALENCE <<endl;
				fout << "BOX3_ION3_RADIUS = " <<  BOX3_ION3_RADIUS <<endl;
				fout << "BOX3_ION3_DIFF_COEFF = " <<  BOX3_ION3_DIFF_COEFF <<endl;
				fout << "BOX3_ION3_MASS = " <<  BOX3_ION3_MASS <<endl;
			}
			if(PRM.BOX3_N4>0){
				fout << "BOX3_ION4_VALENCE = " <<  BOX3_ION4_VALENCE <<endl;
				fout << "BOX3_ION4_RADIUS = " <<  BOX3_ION4_RADIUS <<endl;
				fout << "BOX3_ION4_DIFF_COEFF = " <<  BOX3_ION4_DIFF_COEFF <<endl;
				fout << "BOX3_ION4_MASS = " <<  BOX3_ION4_MASS <<endl;
			}
			if(PRM.BOX3_N5>0){
				fout << "BOX3_ION5_VALENCE = " <<  BOX3_ION5_VALENCE <<endl;
				fout << "BOX3_ION5_RADIUS = " <<  BOX3_ION5_RADIUS <<endl;
				fout << "BOX3_ION5_DIFF_COEFF = " <<  BOX3_ION5_DIFF_COEFF <<endl;
				fout << "BOX3_ION5_MASS = " <<  BOX3_ION5_MASS <<endl;
			}
			if(PRM.BOX3_N6>0){
				fout << "BOX3_ION6_VALENCE = " <<  BOX3_ION6_VALENCE <<endl;
				fout << "BOX3_ION6_RADIUS = " <<  BOX3_ION6_RADIUS <<endl;
				fout << "BOX3_ION6_DIFF_COEFF = " <<  BOX3_ION6_DIFF_COEFF <<endl;
				fout << "BOX3_ION6_MASS = " <<  BOX3_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX4_MIN_Z<PRM.BOX4_MAX_Z){
			if(PRM.BOX4_N1>0){
				fout << "BOX4_ION1_VALENCE = " <<  BOX4_ION1_VALENCE <<endl;
				fout << "BOX4_ION1_RADIUS = " <<  BOX4_ION1_RADIUS <<endl;
				fout << "BOX4_ION1_DIFF_COEFF = " <<  BOX4_ION1_DIFF_COEFF <<endl;
				fout << "BOX4_ION1_MASS = " <<  BOX4_ION1_MASS <<endl;
			}
			if(PRM.BOX4_N2>0){
				fout << "BOX4_ION2_VALENCE = " <<  BOX4_ION2_VALENCE <<endl;
				fout << "BOX4_ION2_RADIUS = " <<  BOX4_ION2_RADIUS <<endl;
				fout << "BOX4_ION2_DIFF_COEFF = " <<  BOX4_ION2_DIFF_COEFF <<endl;
				fout << "BOX4_ION2_MASS = " <<  BOX4_ION2_MASS <<endl;
			}
			if(PRM.BOX4_N3>0){
				fout << "BOX4_ION3_VALENCE = " <<  BOX4_ION3_VALENCE <<endl;
				fout << "BOX4_ION3_RADIUS = " <<  BOX4_ION3_RADIUS <<endl;
				fout << "BOX4_ION3_DIFF_COEFF = " <<  BOX4_ION3_DIFF_COEFF <<endl;
				fout << "BOX4_ION3_MASS = " <<  BOX4_ION3_MASS <<endl;
			}
			if(PRM.BOX4_N4>0){
				fout << "BOX4_ION4_VALENCE = " <<  BOX4_ION4_VALENCE <<endl;
				fout << "BOX4_ION4_RADIUS = " <<  BOX4_ION4_RADIUS <<endl;
				fout << "BOX4_ION4_DIFF_COEFF = " <<  BOX4_ION4_DIFF_COEFF <<endl;
				fout << "BOX4_ION4_MASS = " <<  BOX4_ION4_MASS <<endl;
			}
			if(PRM.BOX4_N5>0){
				fout << "BOX4_ION5_VALENCE = " <<  BOX4_ION5_VALENCE <<endl;
				fout << "BOX4_ION5_RADIUS = " <<  BOX4_ION5_RADIUS <<endl;
				fout << "BOX4_ION5_DIFF_COEFF = " <<  BOX4_ION5_DIFF_COEFF <<endl;
				fout << "BOX4_ION5_MASS = " <<  BOX4_ION5_MASS <<endl;
			}
			if(PRM.BOX4_N6>0){
				fout << "BOX4_ION6_VALENCE = " <<  BOX4_ION6_VALENCE <<endl;
				fout << "BOX4_ION6_RADIUS = " <<  BOX4_ION6_RADIUS <<endl;
				fout << "BOX4_ION6_DIFF_COEFF = " <<  BOX4_ION6_DIFF_COEFF <<endl;
				fout << "BOX4_ION6_MASS = " <<  BOX4_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX5_MIN_Z<PRM.BOX5_MAX_Z){
			if(PRM.BOX5_N1>0){
				fout << "BOX5_ION1_VALENCE = " <<  BOX5_ION1_VALENCE <<endl;
				fout << "BOX5_ION1_RADIUS = " <<  BOX5_ION1_RADIUS <<endl;
				fout << "BOX5_ION1_DIFF_COEFF = " <<  BOX5_ION1_DIFF_COEFF <<endl;
				fout << "BOX5_ION1_MASS = " <<  BOX5_ION1_MASS <<endl;
			}
			if(PRM.BOX5_N2>0){
				fout << "BOX5_ION2_VALENCE = " <<  BOX5_ION2_VALENCE <<endl;
				fout << "BOX5_ION2_RADIUS = " <<  BOX5_ION2_RADIUS <<endl;
				fout << "BOX5_ION2_DIFF_COEFF = " <<  BOX5_ION2_DIFF_COEFF <<endl;
				fout << "BOX5_ION2_MASS = " <<  BOX5_ION2_MASS <<endl;
			}
			if(PRM.BOX5_N3>0){
				fout << "BOX5_ION3_VALENCE = " <<  BOX5_ION3_VALENCE <<endl;
				fout << "BOX5_ION3_RADIUS = " <<  BOX5_ION3_RADIUS <<endl;
				fout << "BOX5_ION3_DIFF_COEFF = " <<  BOX5_ION3_DIFF_COEFF <<endl;
				fout << "BOX5_ION3_MASS = " <<  BOX5_ION3_MASS <<endl;
			}
			if(PRM.BOX5_N4>0){
				fout << "BOX5_ION4_VALENCE = " <<  BOX5_ION4_VALENCE <<endl;
				fout << "BOX5_ION4_RADIUS = " <<  BOX5_ION4_RADIUS <<endl;
				fout << "BOX5_ION4_DIFF_COEFF = " <<  BOX5_ION4_DIFF_COEFF <<endl;
				fout << "BOX5_ION4_MASS = " <<  BOX5_ION4_MASS <<endl;
			}
			if(PRM.BOX5_N5>0){
				fout << "BOX5_ION5_VALENCE = " <<  BOX5_ION5_VALENCE <<endl;
				fout << "BOX5_ION5_RADIUS = " <<  BOX5_ION5_RADIUS <<endl;
				fout << "BOX5_ION5_DIFF_COEFF = " <<  BOX5_ION5_DIFF_COEFF <<endl;
				fout << "BOX5_ION5_MASS = " <<  BOX5_ION5_MASS <<endl;
			}
			if(PRM.BOX5_N6>0){
				fout << "BOX5_ION6_VALENCE = " <<  BOX5_ION6_VALENCE <<endl;
				fout << "BOX5_ION6_RADIUS = " <<  BOX5_ION6_RADIUS <<endl;
				fout << "BOX5_ION6_DIFF_COEFF = " <<  BOX5_ION6_DIFF_COEFF <<endl;
				fout << "BOX5_ION6_MASS = " <<  BOX5_ION6_MASS <<endl;
			}
		}
		if(PRM.BOX6_MIN_Z<PRM.BOX6_MAX_Z){
			if(PRM.BOX6_N1>0){
				fout << "BOX6_ION1_VALENCE = " <<  BOX6_ION1_VALENCE <<endl;
				fout << "BOX6_ION1_RADIUS = " <<  BOX6_ION1_RADIUS <<endl;
				fout << "BOX6_ION1_DIFF_COEFF = " <<  BOX6_ION1_DIFF_COEFF <<endl;
				fout << "BOX6_ION1_MASS = " <<  BOX6_ION1_MASS <<endl;
			}
			if(PRM.BOX6_N2>0){
				fout << "BOX6_ION2_VALENCE = " <<  BOX6_ION2_VALENCE <<endl;
				fout << "BOX6_ION2_RADIUS = " <<  BOX6_ION2_RADIUS <<endl;
				fout << "BOX6_ION2_DIFF_COEFF = " <<  BOX6_ION2_DIFF_COEFF <<endl;
				fout << "BOX6_ION2_MASS = " <<  BOX6_ION2_MASS <<endl;
			}
			if(PRM.BOX6_N3>0){
				fout << "BOX6_ION3_VALENCE = " <<  BOX6_ION3_VALENCE <<endl;
				fout << "BOX6_ION3_RADIUS = " <<  BOX6_ION3_RADIUS <<endl;
				fout << "BOX6_ION3_DIFF_COEFF = " <<  BOX6_ION3_DIFF_COEFF <<endl;
				fout << "BOX6_ION3_MASS = " <<  BOX6_ION3_MASS <<endl;
			}
			if(PRM.BOX6_N4>0){
				fout << "BOX6_ION4_VALENCE = " <<  BOX6_ION4_VALENCE <<endl;
				fout << "BOX6_ION4_RADIUS = " <<  BOX6_ION4_RADIUS <<endl;
				fout << "BOX6_ION4_DIFF_COEFF = " <<  BOX6_ION4_DIFF_COEFF <<endl;
				fout << "BOX6_ION4_MASS = " <<  BOX6_ION4_MASS <<endl;
			}
			if(PRM.BOX6_N5>0){
				fout << "BOX6_ION5_VALENCE = " <<  BOX6_ION5_VALENCE <<endl;
				fout << "BOX6_ION5_RADIUS = " <<  BOX6_ION5_RADIUS <<endl;
				fout << "BOX6_ION5_DIFF_COEFF = " <<  BOX6_ION5_DIFF_COEFF <<endl;
				fout << "BOX6_ION5_MASS = " <<  BOX6_ION5_MASS <<endl;
			}
			if(PRM.BOX6_N6>0){
				fout << "BOX6_ION6_VALENCE = " <<  BOX6_ION6_VALENCE <<endl;
				fout << "BOX6_ION6_RADIUS = " <<  BOX6_ION6_RADIUS <<endl;
				fout << "BOX6_ION6_DIFF_COEFF = " <<  BOX6_ION6_DIFF_COEFF <<endl;
				fout << "BOX6_ION6_MASS = " <<  BOX6_ION6_MASS <<endl;
			}
		}
	}
	
//============== PAGE 8 - FIXED CHARGES
	if(PRM.SIM_TYPE.compare("PORE")==0){
		stream << endl;
		for(int i=0; i<PRM.charge_ring_z.size(); i++){
			stream << "CHARGE_RING = " << PRM.charge_ring_z.at(i) << " = " << PRM.charge_ring_r.at(i) << " = " << PRM.charge_ring_n.at(i) << " = " << PRM.charge_ring_q.at(i) << endl;
		}
	}	
	
	if(PRM.SIM_TYPE.compare("PORE")==0){
		fout << endl;
		for(int i=0; i<PRM.charge_ring_z.size(); i++){
			fout << "CHARGE_RING = " << PRM.charge_ring_z.at(i) << " = " << PRM.charge_ring_r.at(i) << " = " << PRM.charge_ring_n.at(i) << " = " << PRM.charge_ring_q.at(i) << endl;
		}
	}
	
	fout.close();
	
	return;
}
