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
#include "sim_domain.h"
#include "icc.h"
#include "sim_structures.h"



void Square_Matrix_Mult(vector < vector < double > >& A, vector < vector < double > >& B, vector < vector < double > >& C){
	
	int ORDER_OF_MATRIX_A_1=A.size();
	int ORDER_OF_MATRIX_A_2=A.at(0).size();
	int ORDER_OF_MATRIX_B_1=B.size();
	int ORDER_OF_MATRIX_B_2=B.at(0).size();
	int ORDER_OF_MATRIX_C_1=C.size();
	int ORDER_OF_MATRIX_C_2=C.at(0).size();
	
	if(ORDER_OF_MATRIX_A_1!=ORDER_OF_MATRIX_A_2){
		cout << "Matrix A is not square!!" << endl;
		exit(901);
	}
	if(ORDER_OF_MATRIX_B_1!=ORDER_OF_MATRIX_B_2){
		cout << "Matrix B is not square!!" << endl;
		exit(902);
	}
	if(ORDER_OF_MATRIX_C_1!=ORDER_OF_MATRIX_C_2){
		cout << "Matrix C is not square!!" << endl;
		exit(903);
	}
	if(ORDER_OF_MATRIX_A_1!=ORDER_OF_MATRIX_B_1){
		cout << "Matrix A and B have different sizes!!" << endl;
		exit(904);
	}
	if(ORDER_OF_MATRIX_A_1!=ORDER_OF_MATRIX_C_1){
		cout << "Matrix A and C have different sizes!!" << endl;
		exit(905);
	}
	if(ORDER_OF_MATRIX_B_1!=ORDER_OF_MATRIX_C_1){
		cout << "Matrix B and C have different sizes!!" << endl;
		exit(906);
	}
	
	int ORDER_OF_MATRICES=A.size();
	
	for(int row=0; row<ORDER_OF_MATRICES; row++){
		for(int col=0; col<ORDER_OF_MATRICES; col++){
			C.at(row).at(col)=double(0.00);			
			for(int k=0; k<ORDER_OF_MATRICES; k++){
				C.at(row).at(col)+=A.at(row).at(k)*B.at(k).at(col);
			}
		}
	}
			
	return;
}

void create_surfaces(){
	if(PRM.SIM_TYPE.compare("PORE")==0){
		double PORE_RADIUS=0;
		// interpolation of the profile with NUM_OF_DIV parameter
		cout << "create_channel_surfaces............"<<endl;
		surfaces.clear();
		vector <double> zetas; 
		vector <double> erres; 
		//define initial limits guess
		for(int i=PRM.Z_MOUTH_LEFT; i<=PRM.Z_MOUTH_RIGHT; i++){
			zetas.push_back(i);
			erres.push_back(16000);
		}
		double LEFT_VESTIBULE_MAX_CHANNEL_RADIUS=PRM.LEFT_VESTIBULE_MIN_CHANNEL_RADIUS+PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;
		double RIGHT_VESTIBULE_MAX_CHANNEL_RADIUS=PRM.RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS+PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;
		double center_of_circle_z=PRM.Z_MOUTH_LEFT+PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;
		double center_of_circle_r=PRM.LEFT_VESTIBULE_MIN_CHANNEL_RADIUS+PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;
		for(int i=0; i<=int(PRM.LEFT_VESTIBULE_CURVATURE_RADIUS); i++){
			double coseno=PRM.Z_MOUTH_LEFT+double(i)-center_of_circle_z;
			coseno=coseno/PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;
			double angle=acos(coseno);
			double seno=sin(angle);
			erres.at(i)=center_of_circle_r-seno*PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;
		}
		//adjust limits at right mouth
		center_of_circle_z=PRM.Z_MOUTH_RIGHT-PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;
		center_of_circle_r=PRM.RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS+PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;
		for(int i=erres.size()-int(PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS); i<=erres.size()-1; i++){
			double coseno=PRM.Z_MOUTH_LEFT+double(i)-center_of_circle_z;
			coseno=coseno/PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;
			double angle=acos(coseno);
			double seno=sin(angle);
			erres.at(i)=center_of_circle_r-seno*PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;
		}
//==========================================================================================
//==========================================================================================
//		working on channel_profile_points
		for(int i=PRM.channel_profile_points.size()-1; i>=0; i--){
			if(PRM.channel_profile_points.at(i).z<PRM.Z_MOUTH_LEFT+PRM.LEFT_VESTIBULE_CURVATURE_RADIUS || PRM.channel_profile_points.at(i).z>PRM.Z_MOUTH_RIGHT-PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS){
				PRM.channel_profile_points.erase(PRM.channel_profile_points.begin()+i);
			}
		}
		
		
		RadialPoint rp;
		rp.z=PRM.Z_MOUTH_LEFT+PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;
		rp.r=PRM.LEFT_VESTIBULE_MIN_CHANNEL_RADIUS;
		
		if(PRM.channel_profile_points.empty()){
			PRM.channel_profile_points.push_back(rp);
		}
		else{
			for(int i=PRM.channel_profile_points.size()-1; i>=0; i--){
				if(int(PRM.channel_profile_points.at(i).z)==int(rp.z)){
					PRM.channel_profile_points.erase(PRM.channel_profile_points.begin()+i);
				}
			}
			
			PRM.channel_profile_points.push_back(rp);
		}
		
		rp.z=PRM.Z_MOUTH_RIGHT-PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;
		rp.r=PRM.RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS;
		cout << "=========\t" << rp.z << "\t" << rp.r <<endl;
		
		if(PRM.channel_profile_points.empty()){
			PRM.channel_profile_points.push_back(rp);
		}
		else{
			for(int i=PRM.channel_profile_points.size()-1; i>=0; i--){
				if(int(PRM.channel_profile_points.at(i).z)==int(rp.z)){
					PRM.channel_profile_points.erase(PRM.channel_profile_points.begin()+i);
				}
			}
			
			PRM.channel_profile_points.push_back(rp);
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
		
	
		for(int i=0; i<PRM.channel_profile_points.size()-1; i++){
			double z1=PRM.channel_profile_points.at(i).z;
			double r1=PRM.channel_profile_points.at(i).r;
			
			double z2=PRM.channel_profile_points.at(i+1).z;
			double r2=PRM.channel_profile_points.at(i+1).r;
			
			for(int i=int(z1-PRM.Z_MOUTH_LEFT); i<=int(z2-PRM.Z_MOUTH_LEFT); i++){
				erres.at(i)=r1+(r2-r1)*(double(i+PRM.Z_MOUTH_LEFT)-z1)/(z2-z1);
			}
		}
//=================================================
		for(int i=0; i<erres.size(); i++){
			rp.z=zetas.at(i);
			rp.r=erres.at(i);
			channelProfile.push_back(rp);
		}	
		
		limits.clear();
		for(int i=0; i<channelProfile.size(); i++){
			limits.push_back(channelProfile.at(i).r);
		}
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
		int NUM_OF_DIV=PRM.NUM_OF_DIV;
		if(NUM_OF_DIV<=0){
			NUM_OF_DIV=1;
		}
		int NUM_OF_SUB_DIV=PRM.NUM_OF_SUB_DIV;

		double PHI=double(2.00)*M_PI/double(PRM.TILES_PER_RING);
		double SUB_PHI=PHI/double(NUM_OF_SUB_DIV);
		
		double SRL=PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;
		double SRR=PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;
		
		//CENTRAL SECTION
		
		for(int i=0; i<PRM.channel_profile_points.size()-1; i++){
			
			if(PRM.channel_profile_points.at(i).r==PRM.channel_profile_points.at(i+1).r){
				//cylindrical section
				PORE_RADIUS=PRM.channel_profile_points.at(i).r;
				
				NUM_OF_DIV=PRM.NUM_OF_DIV;
				int NUM_OF_DIV2=round((PRM.channel_profile_points.at(i+1).z-PRM.channel_profile_points.at(i).z)/PRM.MAX_TILE_WIDTH);
				
				if(NUM_OF_DIV2<PRM.NUM_OF_DIV){
					NUM_OF_DIV=NUM_OF_DIV2;
				}
				if(NUM_OF_DIV<1){
					NUM_OF_DIV=1;
				}					
				
				
				double cyl_start_z=PRM.channel_profile_points.at(i).z;
				double cyl_end_z=PRM.channel_profile_points.at(i+1).z;
				double cyl_length=cyl_end_z-cyl_start_z;
				double DELTA_Z=cyl_length/double(NUM_OF_DIV);
				double SUB_DELTA_Z=DELTA_Z/double(NUM_OF_SUB_DIV);
				
				
				if(DELTA_Z>0){
	
					for(int i=0; i<NUM_OF_DIV; i++){
						double z_start=cyl_start_z+double(i)*DELTA_Z;
						double z_center=z_start+double(0.5)*DELTA_Z;
						double z_end=z_start+DELTA_Z;
						
						for(int j=0; j<PRM.TILES_PER_RING; j++){
							double phi_start=double(j)*PHI;
							double phi_center=phi_start+double(0.5)*PHI;
							double phi_end=phi_start+PHI;
							
							Surface surface;
							
							//===================== surface vertexes
							surface.pA.x=PORE_RADIUS*cos(phi_start);
							surface.pA.y=PORE_RADIUS*sin(phi_start);
							surface.pA.z=z_start;
							
							surface.pB.x=PORE_RADIUS*cos(phi_start);
							surface.pB.y=PORE_RADIUS*sin(phi_start);
							surface.pB.z=z_end;
							
							surface.pC.x=PORE_RADIUS*cos(phi_end);
							surface.pC.y=PORE_RADIUS*sin(phi_end);
							surface.pC.z=z_end;
							
							surface.pD.x=PORE_RADIUS*cos(phi_end);
							surface.pD.y=PORE_RADIUS*sin(phi_end);
							surface.pD.z=z_start;
						
							surface.center.x=PORE_RADIUS*cos(phi_center);
							surface.center.y=PORE_RADIUS*sin(phi_center);
							surface.center.z=z_center;
							
							
							//===================== surface area
							surface.area=PHI*DELTA_Z*PORE_RADIUS;
							
							
							//===================== surface normal
							surface.normal[0]=cos(phi_center);
							surface.normal[1]=sin(phi_center);
							surface.normal[2]=0.00;
							
							
							//===================== subSurfaces
							vector <Surface> subSurfaces_of_current_surface;
							subSurfaces_of_current_surface.clear();		
							
							for(int ii=0; ii<NUM_OF_SUB_DIV; ii++){
								double sub_z_start=z_start+double(ii)*SUB_DELTA_Z;
								double sub_z_center=sub_z_start+double(0.5)*SUB_DELTA_Z;
								double sub_z_end=sub_z_start+SUB_DELTA_Z;
								
								for(int jj=0; jj<NUM_OF_SUB_DIV; jj++){
									double sub_phi_start=phi_start+double(jj)*SUB_PHI;
									double sub_phi_center=sub_phi_start+double(0.5)*SUB_PHI;
									double sub_phi_end=sub_phi_start+SUB_PHI;
							
									Surface subsurface;
									
									//===================== subsurface vertexes
									subsurface.pA.x=PORE_RADIUS*cos(sub_phi_start);
									subsurface.pA.y=PORE_RADIUS*sin(sub_phi_start);
									subsurface.pA.z=sub_z_start;
									
									subsurface.pB.x=PORE_RADIUS*cos(sub_phi_start);
									subsurface.pB.y=PORE_RADIUS*sin(sub_phi_start);
									subsurface.pB.z=sub_z_end;
									
									subsurface.pC.x=PORE_RADIUS*cos(sub_phi_end);
									subsurface.pC.y=PORE_RADIUS*sin(sub_phi_end);
									subsurface.pC.z=sub_z_end;
									
									subsurface.pD.x=PORE_RADIUS*cos(sub_phi_end);
									subsurface.pD.y=PORE_RADIUS*sin(sub_phi_end);
									subsurface.pD.z=sub_z_start;
								
									subsurface.center.x=PORE_RADIUS*cos(sub_phi_center);
									subsurface.center.y=PORE_RADIUS*sin(sub_phi_center);
									subsurface.center.z=sub_z_center;
									
									
									//===================== subsubsurface area
									subsurface.area=SUB_PHI*SUB_DELTA_Z*PORE_RADIUS;
									
									
									//===================== subsubsurface normal
									subsurface.normal[0]=cos(sub_phi_center);
									subsurface.normal[1]=sin(sub_phi_center);
									subsurface.normal[2]=0.00;
									
									subSurfaces_of_current_surface.push_back(subsurface);       
								}
							}
									
								
							surfaces.push_back(surface);
							subSurfaces.push_back(subSurfaces_of_current_surface);
							subSurfaces_of_current_surface.clear();
							
						}
					}	
				}
				
				
			}
			else{
				//conic section
			
				double cyl_start_z=PRM.channel_profile_points.at(i).z;
				double cyl_end_z=PRM.channel_profile_points.at(i+1).z;
				
				
				//~ cout <<  "cyl_start_z:\t" << cyl_start_z<<endl; 
				//~ cout <<  "cyl_end_z:\t" << cyl_end_z<<endl; 
				
				//~ double cyl_start_z=PRM.Z_MOUTH_LEFT;
				//~ double cyl_end_z=PRM.Z_MOUTH_RIGHT;
				
				//~ double cyl_length=cyl_end_z-cyl_start_z;
				//~ double cyl_length=PRM.Z_MOUTH_LEFT-PRM.Z_MOUTH_RIGHT;

//========================================================================================			
				NUM_OF_DIV=PRM.NUM_OF_DIV;
				if(NUM_OF_DIV<=0){
					NUM_OF_DIV=1;
				}
							
				double DELTA_Z=(cyl_end_z-cyl_start_z)/double(NUM_OF_DIV);
				double SUB_DELTA_Z=DELTA_Z/double(NUM_OF_SUB_DIV);
				
				//~ cout << "DELTA_Z:\t" << DELTA_Z << endl;
				
				int number_of_rows=PRM.NUM_OF_SUB_DIV;    
				
				vector <RadialPoint> nodes;
				RadialPoint node;
				if(DELTA_Z>0){
					
					bool gotIt=true;
					NUM_OF_DIV=0;
					
					do{
						gotIt=true;
						nodes.clear();
						NUM_OF_DIV++;
											
						DELTA_Z=(cyl_end_z-cyl_start_z)/double(NUM_OF_DIV);
						double SUB_DELTA_Z=DELTA_Z/double(NUM_OF_SUB_DIV);
						
						for(int i=0; i<=NUM_OF_DIV; i++){
							node.z=cyl_start_z+double(i)*DELTA_Z;
							int index=node.z-PRM.Z_MOUTH_LEFT;
							node.r=channelProfile.at(index).r;
							nodes.push_back(node);		
						}

						if(NUM_OF_DIV<PRM.NUM_OF_DIV){
							
							//check distance
							for(int i=0; i<nodes.size()-1; i++){
								double distance=get_distance(nodes.at(i), nodes.at(i+1));
								if(distance>PRM.MAX_TILE_WIDTH){
									gotIt=false;
								}
							}
						}
						
					}
					while(!gotIt);
				}	
				
				if(DELTA_Z>0){

					for(int i=0; i<nodes.size()-1; i++){
						double z_start=nodes.at(i).z;
						double z_end=nodes.at(i+1).z;
						double z_center=(z_start+z_end)/double(2.00);
						double r_start=nodes.at(i).r;
						double r_end=nodes.at(i+1).r;
						double r_center=(r_start+r_end)/double(2.00);
						
						int indp1=0, ind0=0;

						int kp1=0, k0=0;
						if(PRM.TILES_PER_RING%2!=0){
							cout << "Parameter TILES_PER_RING must be even!" <<endl;
							exit(99);
						}
						
						double phi_start;
						double phi_center;
						double phi_end;
						
						for(int j=0; j<PRM.TILES_PER_RING*2; j++){
							Surface surface;
							
							
						
							
							//===================== surface vertexes
							if(j%2==0){						
								phi_start=double(j/2)*PHI;
								phi_center=phi_start+double(0.5)*PHI;
								phi_end=phi_start+PHI;
								
								k0=j/2;
								kp1=j/2+1;
								
								surface.pA.x=r_start*cos(phi_start);
								surface.pA.y=r_start*sin(phi_start);
								surface.pA.z=z_start;
								
								surface.pB.x=r_end*cos(phi_start);
								surface.pB.y=r_end*sin(phi_start);
								surface.pB.z=z_end;
								
								surface.pC.x=r_start*cos(phi_end);
								surface.pC.y=r_start*sin(phi_end);
								surface.pC.z=z_start;
								
								surface.pD.x=r_start*cos(phi_end);
								surface.pD.y=r_start*sin(phi_end);
								surface.pD.z=z_start;
							}
							else{						
								k0=j/2;
								kp1=j/2+1;
								surface.pB.x=r_end*cos(phi_start);
								surface.pB.y=r_end*sin(phi_start);
								surface.pB.z=z_end;
							
								surface.pA.x=r_start*cos(phi_end);
								surface.pA.y=r_start*sin(phi_end);
								surface.pA.z=z_start;
								
								surface.pC.x=r_end*cos(phi_end);
								surface.pC.y=r_end*sin(phi_end);
								surface.pC.z=z_end;
								
								surface.pD.x=r_end*cos(phi_end);
								surface.pD.y=r_end*sin(phi_end);
								surface.pD.z=z_end;
							}
							
							surface.center.x=(surface.pA.x+surface.pB.x+surface.pC.x)/3;
							surface.center.y=(surface.pA.y+surface.pB.y+surface.pC.y)/3;
							surface.center.z=(surface.pA.z+surface.pB.z+surface.pC.z)/3;
							
							//===================== surface area
							double dAB, dAC, dCB, semiP;
							dAB=get_distance(surface.pA, surface.pB);
							dAC=get_distance(surface.pA, surface.pC);
							dCB=get_distance(surface.pC, surface.pB);

							semiP=(dAB+dAC+dCB)/2;
							surface.area=sqrt(semiP*(semiP-dAB)*(semiP-dAC)*(semiP-dCB));
							
							
							//===================== surface normal
							double vAB[3], vAC[3], vCB[3];
							dAB=get_distance(surface.pA, surface.pB);
							dAC=get_distance(surface.pA, surface.pC);
							
							vAB[0]=surface.pB.x-surface.pA.x;
							vAB[1]=surface.pB.y-surface.pA.y;
							vAB[2]=surface.pB.z-surface.pA.z;
							
							vAC[0]=surface.pC.x-surface.pA.x;
							vAC[1]=surface.pC.y-surface.pA.y;
							vAC[2]=surface.pC.z-surface.pA.z;
							
							vCB[0]=surface.pB.x-surface.pC.x;
							vCB[1]=surface.pB.y-surface.pC.y;
							vCB[2]=surface.pB.z-surface.pC.z;
							
							for(int iD=0; iD<3; iD++){
								vAB[iD]=vAB[iD]/dAB;
								vAC[iD]=vAC[iD]/dAC;
								vCB[iD]=vCB[iD]/dCB;
							}
							
							double n[3];
							
							n[0]=vAB[1]*vAC[2]-vAB[2]*vAC[1];
							n[1]=vAB[2]*vAC[0]-vAB[0]*vAC[2];
							n[2]=vAB[0]*vAC[1]-vAB[1]*vAC[0];
							
							double normaN=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
							for(int iR=0; iR<3; iR++){
								n[iR]=n[iR]/normaN;
							}
							
							surface.normal[0]=-n[0];
							surface.normal[1]=-n[1];
							surface.normal[2]=-n[2];		

							
							//===================== subSurfaces				
							vector <Surface> subSurfaces_of_current_surface;
							subSurfaces_of_current_surface.clear();	

							//~ int number_of_rows=round(sqrt(double(PRM.NUM_OF_SUB_DIV))); 
							int number_of_rows=PRM.NUM_OF_SUB_DIV; 
							double length_on_AB=dAB/double(number_of_rows);
							double length_on_AC=dAC/double(number_of_rows);
							double length_on_CB=dCB/double(number_of_rows);


							for(int j=0; j<number_of_rows; j++){
								Point begin_row_on_AB;
								begin_row_on_AB.x=surface.pA.x+double(j)*length_on_AB*vAB[0];
								begin_row_on_AB.y=surface.pA.y+double(j)*length_on_AB*vAB[1];
								begin_row_on_AB.z=surface.pA.z+double(j)*length_on_AB*vAB[2];
								
								Point end_row_on_AB;
								end_row_on_AB.x=surface.pA.x+double(j+1)*length_on_AB*vAB[0];
								end_row_on_AB.y=surface.pA.y+double(j+1)*length_on_AB*vAB[1];
								end_row_on_AB.z=surface.pA.z+double(j+1)*length_on_AB*vAB[2];
								
								Point begin_row_on_AC;
								begin_row_on_AC.x=surface.pA.x+double(j)*length_on_AC*vAC[0];
								begin_row_on_AC.y=surface.pA.y+double(j)*length_on_AC*vAC[1];
								begin_row_on_AC.z=surface.pA.z+double(j)*length_on_AC*vAC[2];
								
								Point end_row_on_AC;
								end_row_on_AC.x=surface.pA.x+double(j+1)*length_on_AC*vAC[0];
								end_row_on_AC.y=surface.pA.y+double(j+1)*length_on_AC*vAC[1];
								end_row_on_AC.z=surface.pA.z+double(j+1)*length_on_AC*vAC[2];
								
								
								vCB[0]=end_row_on_AC.x-end_row_on_AB.x;
								vCB[1]=end_row_on_AC.y-end_row_on_AB.y;
								vCB[2]=end_row_on_AC.z-end_row_on_AB.z;
								dCB=get_distance(end_row_on_AB, end_row_on_AC);
								
								for(int iD=0; iD<3; iD++){
									vCB[iD]=vCB[iD]/dCB;
								}
								
								int num_of_tria_per_row=(j+1)*2-1;
								int num_of_sub_low=num_of_tria_per_row/2;
								int num_of_sub_high=num_of_sub_low+1;
								
								for(int h=0; h<num_of_tria_per_row; h++){
									Surface subsurface;

									double length_on_CB_begin=get_distance(begin_row_on_AB, begin_row_on_AC)/double(num_of_sub_low);
									double length_on_CB_end=get_distance(end_row_on_AB, end_row_on_AC)/double(num_of_sub_high);
									
									if(j==0){
										length_on_CB_begin=0.00;
									}
									
									int uu=0, vv=0;
									if(h%2==0){
										uu=h/2;
										vv=uu+1;
										
										subsurface.pA.x=begin_row_on_AB.x+double(uu)*length_on_CB_begin*vCB[0];
										subsurface.pA.y=begin_row_on_AB.y+double(uu)*length_on_CB_begin*vCB[1];
										subsurface.pA.z=begin_row_on_AB.z+double(uu)*length_on_CB_begin*vCB[2];
										
										subsurface.pB.x=end_row_on_AB.x+double(uu)*length_on_CB_end*vCB[0];
										subsurface.pB.y=end_row_on_AB.y+double(uu)*length_on_CB_end*vCB[1];
										subsurface.pB.z=end_row_on_AB.z+double(uu)*length_on_CB_end*vCB[2];
										
										subsurface.pC.x=end_row_on_AB.x+double(vv)*length_on_CB_end*vCB[0];
										subsurface.pC.y=end_row_on_AB.y+double(vv)*length_on_CB_end*vCB[1];
										subsurface.pC.z=end_row_on_AB.z+double(vv)*length_on_CB_end*vCB[2];
										
										subsurface.pD.x=end_row_on_AB.x+double(vv)*length_on_CB_end*vCB[0];
										subsurface.pD.y=end_row_on_AB.y+double(vv)*length_on_CB_end*vCB[1];
										subsurface.pD.z=end_row_on_AB.z+double(vv)*length_on_CB_end*vCB[2];
									}
									else{
										uu=h/2;
										vv=uu+1;
										
										subsurface.pA.x=end_row_on_AB.x+double(vv)*length_on_CB_end*vCB[0];
										subsurface.pA.y=end_row_on_AB.y+double(vv)*length_on_CB_end*vCB[1];
										subsurface.pA.z=end_row_on_AB.z+double(vv)*length_on_CB_end*vCB[2];
										
										subsurface.pB.x=begin_row_on_AB.x+double(vv)*length_on_CB_begin*vCB[0];
										subsurface.pB.y=begin_row_on_AB.y+double(vv)*length_on_CB_begin*vCB[1];
										subsurface.pB.z=begin_row_on_AB.z+double(vv)*length_on_CB_begin*vCB[2];
										
										subsurface.pC.x=begin_row_on_AB.x+double(uu)*length_on_CB_begin*vCB[0];
										subsurface.pC.y=begin_row_on_AB.y+double(uu)*length_on_CB_begin*vCB[1];
										subsurface.pC.z=begin_row_on_AB.z+double(uu)*length_on_CB_begin*vCB[2];
										
										subsurface.pD.x=begin_row_on_AB.x+double(uu)*length_on_CB_begin*vCB[0];
										subsurface.pD.y=begin_row_on_AB.y+double(uu)*length_on_CB_begin*vCB[1];
										subsurface.pD.z=begin_row_on_AB.z+double(uu)*length_on_CB_begin*vCB[2];
									}
									
									subsurface.center.x=(subsurface.pA.x+subsurface.pB.x+subsurface.pC.x)/3;
									subsurface.center.y=(subsurface.pA.y+subsurface.pB.y+subsurface.pC.y)/3;
									subsurface.center.z=(subsurface.pA.z+subsurface.pB.z+subsurface.pC.z)/3;
									
									double sub_vAB[3], sub_vAC[3], sub_vCB[3];
									double sub_dAB, sub_dAC, sub_dCB, sub_semiP;
									
									sub_dAB=get_distance(subsurface.pA, subsurface.pB);
									sub_dAC=get_distance(subsurface.pA, subsurface.pC);
									
									sub_vAB[0]=subsurface.pB.x-subsurface.pA.x;
									sub_vAB[1]=subsurface.pB.y-subsurface.pA.y;
									sub_vAB[2]=subsurface.pB.z-subsurface.pA.z;
									
									sub_vAC[0]=subsurface.pC.x-subsurface.pA.x;
									sub_vAC[1]=subsurface.pC.y-subsurface.pA.y;
									sub_vAC[2]=subsurface.pC.z-subsurface.pA.z;
									
									sub_vCB[0]=subsurface.pB.x-subsurface.pC.x;
									sub_vCB[1]=subsurface.pB.y-subsurface.pC.y;
									sub_vCB[2]=subsurface.pB.z-subsurface.pC.z;
									
									for(int iD=0; iD<3; iD++){
										sub_vAB[iD]=sub_vAB[iD]/sub_dAB;
										sub_vAC[iD]=sub_vAC[iD]/sub_dAC;
									}
									
									double sub_n[3];
									
									sub_n[0]=sub_vAB[1]*sub_vAC[2]-sub_vAB[2]*sub_vAC[1];
									sub_n[1]=sub_vAB[2]*sub_vAC[0]-sub_vAB[0]*sub_vAC[2];
									sub_n[2]=sub_vAB[0]*sub_vAC[1]-sub_vAB[1]*sub_vAC[0];
									
									double sub_normaN=sqrt(sub_n[0]*sub_n[0]+sub_n[1]*sub_n[1]+sub_n[2]*sub_n[2]);
									for(int iR=0; iR<3; iR++){
										//change verse and normalize
										sub_n[iR]=-sub_n[iR]/sub_normaN;
									}
									
									subsurface.normal[0]=sub_n[0];
									subsurface.normal[1]=sub_n[1];
									subsurface.normal[2]=sub_n[2];		
									
									sub_dAB=get_distance(subsurface.pA, subsurface.pB);
									sub_dAC=get_distance(subsurface.pA, subsurface.pC);
									sub_dCB=get_distance(subsurface.pC, subsurface.pB);

									sub_semiP=(sub_dAB+sub_dAC+sub_dCB)/2;
									subsurface.area=sqrt(sub_semiP*(sub_semiP-sub_dAB)*(sub_semiP-sub_dAC)*(sub_semiP-sub_dCB));			
									
									
									subSurfaces_of_current_surface.push_back(subsurface);
								}
							
							}		
							
							surfaces.push_back(surface);
							subSurfaces.push_back(subSurfaces_of_current_surface);
							subSurfaces_of_current_surface.clear();
						}
					}			
				}		
			}
			//~ cout << "qui?" <<endl;
		}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
		
		NUM_OF_DIV=PRM.NUM_OF_DIV;		
		NUM_OF_DIV=double(1.5)*PRM.LEFT_VESTIBULE_CURVATURE_RADIUS/PRM.MAX_TILE_WIDTH;
		if(NUM_OF_DIV<=0){
			NUM_OF_DIV=1;
		}
		//~ cout << "NUM_OF_DIV:\t" << NUM_OF_DIV <<endl;
		
		//MOUTHS	
		double cyl_start_z=PRM.Z_MOUTH_LEFT+PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;	
		double cyl_end_z=PRM.Z_MOUTH_RIGHT-PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;	
		double cyl_length=cyl_end_z-cyl_start_z;
		double DELTA_Z=cyl_length/double(NUM_OF_DIV);
		double SUB_DELTA_Z=DELTA_Z/double(NUM_OF_SUB_DIV);
		
		
		double THETA=M_PI/(double(2.00)*double(NUM_OF_DIV));
		
		//~ double circ_seg=THETA*SRL;
		//~ THETA=DELTA_Z/SRL;
		NUM_OF_DIV=round(M_PI/(double(2.00)*THETA));
		THETA=M_PI/(double(2.00)*double(NUM_OF_DIV));

		double SUB_THETA=THETA/double(NUM_OF_SUB_DIV);
		
	//LEFT_MOUTH
		double OZ=PRM.Z_MOUTH_LEFT+SRL;			
		double OR=channelProfile.at(0).r;
		
		//~ cout << "NUM_OF_DIV:\t" << NUM_OF_DIV <<endl;
		
		for(int i=0; i<NUM_OF_DIV; i++){
			double theta_start=double(i)*THETA;
			double theta_center=theta_start+double(0.5)*THETA;
			double theta_end=theta_start+THETA;
			
			double z_start=OZ-cos(theta_start)*SRL;
			double z_center=OZ-cos(theta_center)*SRL;
			double z_end=OZ-cos(theta_end)*SRL;
			
			double R_start=OR-sin(theta_start)*SRL;
			double R_center=OR-sin(theta_center)*SRL;
			double R_end=OR-sin(theta_end)*SRL;
			
			
			for(int j=0; j<PRM.TILES_PER_RING; j++){
				double phi_start=double(j)*PHI;
				double phi_center=phi_start+double(0.5)*PHI;
				double phi_end=phi_start+PHI;
				
				Surface surface;
				
				//===================== surface vertexes
				surface.pA.x=R_start*cos(phi_start);
				surface.pA.y=R_start*sin(phi_start);
				surface.pA.z=z_start;
				
				surface.pB.x=R_end*cos(phi_start);
				surface.pB.y=R_end*sin(phi_start);
				surface.pB.z=z_end;
				
				surface.pC.x=R_end*cos(phi_end);
				surface.pC.y=R_end*sin(phi_end);
				surface.pC.z=z_end;
				
				surface.pD.x=R_start*cos(phi_end);
				surface.pD.y=R_start*sin(phi_end);
				surface.pD.z=z_start;
			
				surface.center.x=R_center*cos(phi_center);
				surface.center.y=R_center*sin(phi_center);
				surface.center.z=z_center;
				
				
				//===================== surface area
				double k1=PHI*SRL;
				double t1=OR*THETA;
				double t2=SRL*(cos(theta_end)-cos(theta_start));
				surface.area=k1*(t1+t2);
				
				
				//===================== surface normal
				surface.normal[0]=sin(theta_center)*cos(phi_center);
				surface.normal[1]=sin(theta_center)*sin(phi_center);
				surface.normal[2]=cos(theta_center);
				
				
				//===================== subSurfaces
				vector <Surface> subSurfaces_of_current_surface;
				subSurfaces_of_current_surface.clear();		
				
				for(int ii=0; ii<NUM_OF_SUB_DIV; ii++){
					double sub_theta_start=theta_start+double(ii)*SUB_THETA;
					double sub_theta_center=sub_theta_start+double(0.5)*SUB_THETA;
					double sub_theta_end=sub_theta_start+SUB_THETA;
					
					double sub_z_start=OZ-cos(sub_theta_start)*SRL;
					double sub_z_center=OZ-cos(sub_theta_center)*SRL;
					double sub_z_end=OZ-cos(sub_theta_end)*SRL;
					
					double sub_R_start=OR-sin(sub_theta_start)*SRL;
					double sub_R_center=OR-sin(sub_theta_center)*SRL;
					double sub_R_end=OR-sin(sub_theta_end)*SRL;
					
					
					for(int jj=0; jj<NUM_OF_SUB_DIV; jj++){
						double sub_phi_start=phi_start+double(jj)*SUB_PHI;
						double sub_phi_center=sub_phi_start+double(0.5)*SUB_PHI;
						double sub_phi_end=sub_phi_start+SUB_PHI;
				
						Surface subsurface;
						
						//===================== subsurface vertexes
						subsurface.pA.x=sub_R_start*cos(sub_phi_start);
						subsurface.pA.y=sub_R_start*sin(sub_phi_start);
						subsurface.pA.z=sub_z_start;
						
						subsurface.pB.x=sub_R_end*cos(sub_phi_start);
						subsurface.pB.y=sub_R_end*sin(sub_phi_start);
						subsurface.pB.z=sub_z_end;
						
						subsurface.pC.x=sub_R_end*cos(sub_phi_end);
						subsurface.pC.y=sub_R_end*sin(sub_phi_end);
						subsurface.pC.z=sub_z_end;
						
						subsurface.pD.x=sub_R_start*cos(sub_phi_end);
						subsurface.pD.y=sub_R_start*sin(sub_phi_end);
						subsurface.pD.z=sub_z_start;
					
						subsurface.center.x=sub_R_center*cos(sub_phi_center);
						subsurface.center.y=sub_R_center*sin(sub_phi_center);
						subsurface.center.z=sub_z_center;
						
						
						//===================== subsubsurface area
						double sub_k1=SUB_PHI*SRL;
						double sub_t1=OR*SUB_THETA;
						double sub_t2=SRL*(cos(sub_theta_end)-cos(sub_theta_start));
						subsurface.area=sub_k1*(sub_t1+sub_t2);
						
						
						//===================== subsubsurface normal
						subsurface.normal[0]=sin(sub_theta_center)*cos(sub_phi_center);
						subsurface.normal[1]=sin(sub_theta_center)*sin(sub_phi_center);
						subsurface.normal[2]=cos(sub_theta_center);
						
						
						
						subSurfaces_of_current_surface.push_back(subsurface);       
					}
				}
						
					
				surfaces.push_back(surface);
				subSurfaces.push_back(subSurfaces_of_current_surface);
				subSurfaces_of_current_surface.clear();
			}
		}
		
		
//RIGHT_MOUTH
		
		NUM_OF_DIV=PRM.NUM_OF_DIV;		
		NUM_OF_DIV=double(1.5)*PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS/PRM.MAX_TILE_WIDTH;
		if(NUM_OF_DIV<=0){
			NUM_OF_DIV=1;
		}
		//~ cout << "NUM_OF_DIV:\t" << NUM_OF_DIV <<endl;
		
		//MOUTHS	
		//~ double cyl_start_z=PRM.Z_MOUTH_LEFT+PRM.LEFT_VESTIBULE_CURVATURE_RADIUS;	
		//~ double cyl_end_z=PRM.Z_MOUTH_RIGHT-PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;	
		//~ double cyl_length=cyl_end_z-cyl_start_z;
		//~ double DELTA_Z=cyl_length/double(NUM_OF_DIV);
		//~ double SUB_DELTA_Z=DELTA_Z/double(NUM_OF_SUB_DIV);
		
		
		THETA=M_PI/(double(2.00)*double(NUM_OF_DIV));
		
		//~ double circ_seg=THETA*SRL;
		//~ THETA=DELTA_Z/SRL;
		NUM_OF_DIV=round(M_PI/(double(2.00)*THETA));
		THETA=M_PI/(double(2.00)*double(NUM_OF_DIV));

		SUB_THETA=THETA/double(NUM_OF_SUB_DIV);
		//~ cout << "NUM_OF_DIV:\t" << NUM_OF_DIV <<endl;
		
	
		SRR=PRM.RIGHT_VESTIBULE_CURVATURE_RADIUS;
		OZ=PRM.Z_MOUTH_RIGHT-SRR;			
		OR=channelProfile.at(channelProfile.size()-1).r;
		
		for(int i=0; i<NUM_OF_DIV; i++){
			double theta_start=double(i)*THETA;
			double theta_center=theta_start+double(0.5)*THETA;
			double theta_end=theta_start+THETA;
			
			double z_start=OZ+sin(theta_start)*SRR;
			double z_center=OZ+sin(theta_center)*SRR;
			double z_end=OZ+sin(theta_end)*SRR;
			
			double R_start=OR-cos(theta_start)*SRR;
			double R_center=OR-cos(theta_center)*SRR;
			double R_end=OR-cos(theta_end)*SRR;
			
			
			for(int j=0; j<PRM.TILES_PER_RING; j++){
				double phi_start=double(j)*PHI;
				double phi_center=phi_start+double(0.5)*PHI;
				double phi_end=phi_start+PHI;
				
				Surface surface;
				
				//===================== surface vertexes
				surface.pA.x=R_start*cos(phi_start);
				surface.pA.y=R_start*sin(phi_start);
				surface.pA.z=z_start;
				
				surface.pB.x=R_end*cos(phi_start);
				surface.pB.y=R_end*sin(phi_start);
				surface.pB.z=z_end;
				
				surface.pC.x=R_end*cos(phi_end);
				surface.pC.y=R_end*sin(phi_end);
				surface.pC.z=z_end;
				
				surface.pD.x=R_start*cos(phi_end);
				surface.pD.y=R_start*sin(phi_end);
				surface.pD.z=z_start;
			
				surface.center.x=R_center*cos(phi_center);
				surface.center.y=R_center*sin(phi_center);
				surface.center.z=z_center;
				
				
				//===================== surface area
				double k1=PHI*SRR;
				double t1=OR*THETA;
				double t2=SRR*(sin(theta_start)-sin(theta_end));
				surface.area=k1*(t1+t2);
				
				
				//===================== surface normal
				surface.normal[0]=cos(theta_center)*cos(phi_center);
				surface.normal[1]=cos(theta_center)*sin(phi_center);
				surface.normal[2]=-sin(theta_center);
				
				
				//===================== subSurfaces
				vector <Surface> subSurfaces_of_current_surface;
				subSurfaces_of_current_surface.clear();		
				
				for(int ii=0; ii<NUM_OF_SUB_DIV; ii++){
					double sub_theta_start=theta_start+double(ii)*SUB_THETA;
					double sub_theta_center=sub_theta_start+double(0.5)*SUB_THETA;
					double sub_theta_end=sub_theta_start+SUB_THETA;
					
					double sub_z_start=OZ+sin(sub_theta_start)*SRR;
					double sub_z_center=OZ+sin(sub_theta_center)*SRR;
					double sub_z_end=OZ+sin(sub_theta_end)*SRR;
					
					double sub_R_start=OR-cos(sub_theta_start)*SRR;
					double sub_R_center=OR-cos(sub_theta_center)*SRR;
					double sub_R_end=OR-cos(sub_theta_end)*SRR;
					
					for(int jj=0; jj<NUM_OF_SUB_DIV; jj++){
						double sub_phi_start=phi_start+double(jj)*SUB_PHI;
						double sub_phi_center=sub_phi_start+double(0.5)*SUB_PHI;
						double sub_phi_end=sub_phi_start+SUB_PHI;
				
						Surface subsurface;
						
						//===================== subsurface vertexes
						subsurface.pA.x=sub_R_start*cos(sub_phi_start);
						subsurface.pA.y=sub_R_start*sin(sub_phi_start);
						subsurface.pA.z=sub_z_start;
						
						subsurface.pB.x=sub_R_end*cos(sub_phi_start);
						subsurface.pB.y=sub_R_end*sin(sub_phi_start);
						subsurface.pB.z=sub_z_end;
						
						subsurface.pC.x=sub_R_end*cos(sub_phi_end);
						subsurface.pC.y=sub_R_end*sin(sub_phi_end);
						subsurface.pC.z=sub_z_end;
						
						subsurface.pD.x=sub_R_start*cos(sub_phi_end);
						subsurface.pD.y=sub_R_start*sin(sub_phi_end);
						subsurface.pD.z=sub_z_start;
					
						subsurface.center.x=sub_R_center*cos(sub_phi_center);
						subsurface.center.y=sub_R_center*sin(sub_phi_center);
						subsurface.center.z=sub_z_center;
						
						
						//===================== subsubsurface area
						double sub_k1=SUB_PHI*SRR;
						double sub_t1=OR*SUB_THETA;
						double sub_t2=SRR*(sin(sub_theta_start)-sin(sub_theta_end));
						subsurface.area=sub_k1*(sub_t1+sub_t2);
						
						
						//===================== subsubsurface normal
						subsurface.normal[0]=cos(sub_theta_center)*cos(sub_phi_center);
						subsurface.normal[1]=cos(sub_theta_center)*sin(sub_phi_center);
						subsurface.normal[2]=-sin(sub_theta_center);
						
						
						
						subSurfaces_of_current_surface.push_back(subsurface);       
					}
				}
						
					
				surfaces.push_back(surface);
				subSurfaces.push_back(subSurfaces_of_current_surface);
				subSurfaces_of_current_surface.clear();
				
			}
		}	
				
//MEMBRANES
	
		int TILES_PER_RING_MEMBRANE=round(double(1.50)*double(PRM.TILES_PER_RING));
		PHI=double(2.00)*M_PI/(double(TILES_PER_RING_MEMBRANE));
		SUB_PHI=PHI/double(NUM_OF_SUB_DIV);
		

		double mem_start_R=channelProfile.at(0).r;
		double mem_end_R=PRM.SIM_DOMAIN_WIDTH_X*double(0.50)*double(1.44);	
		
		double mem_length=mem_start_R;
		int mem_NUM_OF_DIV=NUM_OF_DIV*4;
		
		double RING_GROWING_FACTOR=1.33;

		
		DELTA_Z=PRM.MAX_TILE_WIDTH;
		double surface_length=PRM.MAX_TILE_WIDTH;
		mem_NUM_OF_DIV=1;
		for(int i=0; i<1000; i++){
			if(i!=0){
				surface_length=RING_GROWING_FACTOR*surface_length;
			}
			mem_length=mem_length+surface_length;
			if(mem_length<mem_end_R){
				mem_NUM_OF_DIV++;
			}
		}
		
		if(mem_NUM_OF_DIV>0){

//LEFT MEMBRANE
			double DELTA_R=DELTA_Z;
			double SUB_DELTA_R=DELTA_R/double(NUM_OF_SUB_DIV);
			double R_start=mem_start_R;
			
			for(int i=0; i<mem_NUM_OF_DIV; i++){
				double R_center=R_start+double(0.5)*DELTA_R;
				double R_end=R_start+DELTA_R;
				SUB_DELTA_R=DELTA_R/double(NUM_OF_SUB_DIV);
				
				for(int j=0; j<TILES_PER_RING_MEMBRANE; j++){
					double phi_start=double(j)*PHI;
					double phi_center=phi_start+double(0.5)*PHI;
					double phi_end=phi_start+PHI;
					
					Surface surface;
					
					//===================== surface vertexes
					surface.pA.x=R_start*cos(phi_start);
					surface.pA.y=R_start*sin(phi_start);
					surface.pA.z=PRM.Z_MOUTH_LEFT;
					
					surface.pB.x=R_end*cos(phi_start);
					surface.pB.y=R_end*sin(phi_start);
					surface.pB.z=PRM.Z_MOUTH_LEFT;
					
					surface.pC.x=R_end*cos(phi_end);
					surface.pC.y=R_end*sin(phi_end);
					surface.pC.z=PRM.Z_MOUTH_LEFT;
					
					surface.pD.x=R_start*cos(phi_end);
					surface.pD.y=R_start*sin(phi_end);
					surface.pD.z=PRM.Z_MOUTH_LEFT;
				
					surface.center.x=R_center*cos(phi_center);
					surface.center.y=R_center*sin(phi_center);
					surface.center.z=PRM.Z_MOUTH_LEFT;
					
					
					//===================== surface area
					double R_start_2=R_start*R_start;
					double R_end_2=R_end*R_end;
					surface.area=double(0.5)*PHI*(R_end_2-R_start_2);
					
					
					//===================== surface normal
					surface.normal[0]=0.00;
					surface.normal[1]=0.00;
					surface.normal[2]=1.00;
					
					
					//===================== subSurfaces
					vector <Surface> subSurfaces_of_current_surface;
					subSurfaces_of_current_surface.clear();		
					
					for(int ii=0; ii<NUM_OF_SUB_DIV; ii++){
						double sub_R_start=R_start+double(ii)*SUB_DELTA_R;
						double sub_R_center=sub_R_start+double(0.5)*SUB_DELTA_R;
						double sub_R_end=sub_R_start+SUB_DELTA_R;
						
						for(int jj=0; jj<NUM_OF_SUB_DIV; jj++){
							double sub_phi_start=phi_start+double(jj)*SUB_PHI;
							double sub_phi_center=sub_phi_start+double(0.5)*SUB_PHI;
							double sub_phi_end=sub_phi_start+SUB_PHI;
					
							Surface subsurface;
							
							//===================== subsurface vertexes
							subsurface.pA.x=sub_R_start*cos(sub_phi_start);
							subsurface.pA.y=sub_R_start*sin(sub_phi_start);
							subsurface.pA.z=PRM.Z_MOUTH_LEFT;
							
							subsurface.pB.x=sub_R_end*cos(sub_phi_start);
							subsurface.pB.y=sub_R_end*sin(sub_phi_start);
							subsurface.pB.z=PRM.Z_MOUTH_LEFT;
							
							subsurface.pC.x=sub_R_end*cos(sub_phi_end);
							subsurface.pC.y=sub_R_end*sin(sub_phi_end);
							subsurface.pC.z=PRM.Z_MOUTH_LEFT;
							
							subsurface.pD.x=sub_R_start*cos(sub_phi_end);
							subsurface.pD.y=sub_R_start*sin(sub_phi_end);
							subsurface.pD.z=PRM.Z_MOUTH_LEFT;
						
							subsurface.center.x=sub_R_center*cos(sub_phi_center);
							subsurface.center.y=sub_R_center*sin(sub_phi_center);
							subsurface.center.z=PRM.Z_MOUTH_LEFT;
							
							
							//===================== subsubsurface area
							subsurface.area=SUB_PHI*SUB_DELTA_Z*PORE_RADIUS;
							
							double sub_R_start_2=sub_R_start*sub_R_start;
							double sub_R_end_2=sub_R_end*sub_R_end;
							subsurface.area=double(0.5)*SUB_PHI*(sub_R_end_2-sub_R_start_2);
									
							
							//===================== subsubsurface normal
							subsurface.normal[0]=0.00;
							subsurface.normal[1]=0.00;
							subsurface.normal[2]=1.00;
							
							
							
							subSurfaces_of_current_surface.push_back(subsurface);       
						}
					}
							
						
					surfaces.push_back(surface);
					subSurfaces.push_back(subSurfaces_of_current_surface);
					subSurfaces_of_current_surface.clear();
					
				}
				
				R_start=R_start+DELTA_R;
				DELTA_R=RING_GROWING_FACTOR*DELTA_R;
			}	
		
//RIGHT MEMBRANE
			//~ DELTA_R=DELTA_Z;
			DELTA_R=PRM.MAX_TILE_WIDTH;
			SUB_DELTA_R=DELTA_R/double(NUM_OF_SUB_DIV);
			mem_start_R=channelProfile.at(channelProfile.size()-1).r;
			R_start=mem_start_R;
			
			
			for(int i=0; i<mem_NUM_OF_DIV; i++){
				double R_center=R_start+double(0.5)*DELTA_R;
				double R_end=R_start+DELTA_R;
				SUB_DELTA_R=DELTA_R/double(NUM_OF_SUB_DIV);
				
				for(int j=0; j<TILES_PER_RING_MEMBRANE; j++){
					double phi_start=double(j)*PHI;
					double phi_center=phi_start+double(0.5)*PHI;
					double phi_end=phi_start+PHI;
					
					Surface surface;
					
					//===================== surface vertexes
					surface.pA.x=R_start*cos(phi_start);
					surface.pA.y=R_start*sin(phi_start);
					surface.pA.z=PRM.Z_MOUTH_RIGHT;
					
					surface.pB.x=R_end*cos(phi_start);
					surface.pB.y=R_end*sin(phi_start);
					surface.pB.z=PRM.Z_MOUTH_RIGHT;
					
					surface.pC.x=R_end*cos(phi_end);
					surface.pC.y=R_end*sin(phi_end);
					surface.pC.z=PRM.Z_MOUTH_RIGHT;
					
					surface.pD.x=R_start*cos(phi_end);
					surface.pD.y=R_start*sin(phi_end);
					surface.pD.z=PRM.Z_MOUTH_RIGHT;
				
					surface.center.x=R_center*cos(phi_center);
					surface.center.y=R_center*sin(phi_center);
					surface.center.z=PRM.Z_MOUTH_RIGHT;
					
					
					//===================== surface area
					double R_start_2=R_start*R_start;
					double R_end_2=R_end*R_end;
					surface.area=double(0.5)*PHI*(R_end_2-R_start_2);
					
					
					//===================== surface normal
					surface.normal[0]=0.00;
					surface.normal[1]=0.00;
					surface.normal[2]=-1.00;
					
					
					//===================== subSurfaces
					vector <Surface> subSurfaces_of_current_surface;
					subSurfaces_of_current_surface.clear();		
					
					for(int ii=0; ii<NUM_OF_SUB_DIV; ii++){
						double sub_R_start=R_start+double(ii)*SUB_DELTA_R;
						double sub_R_center=sub_R_start+double(0.5)*SUB_DELTA_R;
						double sub_R_end=sub_R_start+SUB_DELTA_R;
						
						for(int jj=0; jj<NUM_OF_SUB_DIV; jj++){
							double sub_phi_start=phi_start+double(jj)*SUB_PHI;
							double sub_phi_center=sub_phi_start+double(0.5)*SUB_PHI;
							double sub_phi_end=sub_phi_start+SUB_PHI;
					
							Surface subsurface;
							
							//===================== subsurface vertexes
							subsurface.pA.x=sub_R_start*cos(sub_phi_start);
							subsurface.pA.y=sub_R_start*sin(sub_phi_start);
							subsurface.pA.z=PRM.Z_MOUTH_RIGHT;
							
							subsurface.pB.x=sub_R_end*cos(sub_phi_start);
							subsurface.pB.y=sub_R_end*sin(sub_phi_start);
							subsurface.pB.z=PRM.Z_MOUTH_RIGHT;
							
							subsurface.pC.x=sub_R_end*cos(sub_phi_end);
							subsurface.pC.y=sub_R_end*sin(sub_phi_end);
							subsurface.pC.z=PRM.Z_MOUTH_RIGHT;
							
							subsurface.pD.x=sub_R_start*cos(sub_phi_end);
							subsurface.pD.y=sub_R_start*sin(sub_phi_end);
							subsurface.pD.z=PRM.Z_MOUTH_RIGHT;
						
							subsurface.center.x=sub_R_center*cos(sub_phi_center);
							subsurface.center.y=sub_R_center*sin(sub_phi_center);
							subsurface.center.z=PRM.Z_MOUTH_RIGHT;
							
							
							//===================== subsubsurface area
							subsurface.area=SUB_PHI*SUB_DELTA_Z*PORE_RADIUS;
							
							double sub_R_start_2=sub_R_start*sub_R_start;
							double sub_R_end_2=sub_R_end*sub_R_end;
							subsurface.area=double(0.5)*SUB_PHI*(sub_R_end_2-sub_R_start_2);
									
							
							//===================== subsubsurface normal
							subsurface.normal[0]=0.00;
							subsurface.normal[1]=0.00;
							subsurface.normal[2]=-1.00;
							
							
							
							subSurfaces_of_current_surface.push_back(subsurface);       
						}
					}
							
						
					surfaces.push_back(surface);
					subSurfaces.push_back(subSurfaces_of_current_surface);
					subSurfaces_of_current_surface.clear();
					
				}
				
				R_start=R_start+DELTA_R;
				DELTA_R=RING_GROWING_FACTOR*DELTA_R;
			}		
		}
		
		
		
	}
		

	else if(PRM.SIM_TYPE.compare("MEMBRANE")==0){
		cout << "PRM.SIM_TYPE=MEMBRANE" <<endl;
		create_membrane_tiles();
	}
		
	PRM.NUMBER_OF_TILES=surfaces.size();
	PRM.NUMBER_OF_SUBTILES_PER_TILE=subSurfaces.at(0).size();
	cout << 	"NUMBER_OF_TILES: " << PRM.NUMBER_OF_TILES<<endl;
	cout << 	"NUMBER_OF_SUBTILES_PER_TILE: " << PRM.NUMBER_OF_SUBTILES_PER_TILE<<endl;
	if(surfaces.empty()){
		cerr << "There are no surfaces!!!"<<endl;
		exit(1);
	}
	else{
		//~ print_surfaces(); // write the discretizion tiles to output files
		//~ print_subSurfaces();		
	}
	double delta_epsilon=PRM.EPS_MEM-PRM.EPS_W;
	double mean_epsilon=(PRM.EPS_MEM+PRM.EPS_W)/double(2.00);
	double den_for_K=double(4.00)*M_PI*mean_epsilon;
	double K_for_integral=delta_epsilon/den_for_K;
	vector<double> aux_double_vec;
	double aux_double=0.00;	
	for(int row=0; row<surfaces.size(); row++){
		aux_double_vec.clear();
		for(int col=0; col<surfaces.size(); col++){
			aux_double_vec.push_back(aux_double);
		}
		MATRIX_A.push_back(aux_double_vec);
		aux_double_vec.clear();
	}	 



	for(int row=0; row<surfaces.size(); row++){
		
		for(int col=0; col<surfaces.size(); col++){

			double delta_kronecker=0.00;
			if(row==col){
				delta_kronecker=1.00;
			}

			double I=0.00;

			for(int st=0; st<subSurfaces.at(row).size(); st++){
				double cc_0=subSurfaces.at(row).at(st).center.x-surfaces.at(col).center.x;
				double cc_1=subSurfaces.at(row).at(st).center.y-surfaces.at(col).center.y;
				double cc_2=subSurfaces.at(row).at(st).center.z-surfaces.at(col).center.z;
				double cc_0_square=cc_0*cc_0;
				double cc_1_square=cc_1*cc_1;
				double cc_2_square=cc_2*cc_2;
				double temp_cc=cc_0_square+cc_1_square+cc_2_square;
				double ccLen=sqrt(temp_cc);
				if(ccLen>1e-5){
					double ccLen_square=ccLen*ccLen;
					double cc_00=cc_0/ccLen;
					double cc_11=cc_1/ccLen;
					double cc_22=cc_2/ccLen;

					double sc_0=cc_00*subSurfaces.at(row).at(st).normal[0];
					double sc_1=cc_11*subSurfaces.at(row).at(st).normal[1];
					double sc_2=cc_22*subSurfaces.at(row).at(st).normal[2];
					double scalar=sc_0+sc_1+sc_2;
					double temp1=scalar*subSurfaces.at(row).at(st).area;
					double temp2=temp1/ccLen_square;

					I+=temp2;	
				}
			}

			double elem=delta_kronecker*surfaces.at(row).area+K_for_integral*surfaces.at(col).area*I;

			MATRIX_A.at(row).at(col)=elem;
		}
	}
	
	

	for(int i=0; i<subSurfaces.size(); i++){
		surfaces.at(i).area=surfaces.at(i).area*1e-24;
		for(int j=0; j<subSurfaces.at(i).size(); j++){
			subSurfaces.at(i).at(j).area=subSurfaces.at(i).at(j).area*1e-24;
		}
	}

	vector < vector < double > > INVERSE_MATRIX_A;
	aux_double_vec.clear();
	
	
	ORDER_OF_MATRIX = MATRIX_A.size();
	matrix_size = MATRIX_A.size();
	
	_int_one=1;
	_double_one=1.00;
	_double_zero=0.00;
	_N='N';
	
	
	
	
	int *ipiv;
	double *work,*a,*inverse_a_tmp;
	double workspace;
	double numero=0;
	int lwork, tmp = -1, info = 0;
	int i,j;

	cout<<"ORDER OF MATRIX FOR ICC: "<<ORDER_OF_MATRIX<<endl;
	for(i=1; i<=ORDER_OF_MATRIX; i++){
		aux_double_vec.clear();
		for(j=1; j<=ORDER_OF_MATRIX; j++){
			numero=0.00;
			aux_double_vec.push_back(numero);
		}
		INVERSE_MATRIX_A.push_back(aux_double_vec);
		aux_double_vec.clear();
	}
	
	posix_memalign((void **)(&ipiv), 128, sizeof(int)*ORDER_OF_MATRIX);
	posix_memalign((void **)(&a), 128, sizeof(double)*ORDER_OF_MATRIX*ORDER_OF_MATRIX);
	
	if( (a == NULL) || (ipiv == NULL) ){
		cerr<<"a/ipiv malloc error"<<endl;
		exit(1);
	}
	for(i=0; i<ORDER_OF_MATRIX; i++){
		for(j=0; j<ORDER_OF_MATRIX; j++){
			a[i+j*ORDER_OF_MATRIX]=MATRIX_A[i][j];
		}
	}
	
	
	//---------Call Cell LAPACK library---------
	dgetrf_(&ORDER_OF_MATRIX, &ORDER_OF_MATRIX, a, &ORDER_OF_MATRIX, ipiv, &info);
	if( info != 0 ){
		cout << "ORDER_OF_MATRIX: " << ORDER_OF_MATRIX << endl;
		cerr<<"Call dgetrf error"<<endl;
		exit(1);
	}
	
	//---------Query workspace-------
	dgetri_(&ORDER_OF_MATRIX, a, &ORDER_OF_MATRIX, ipiv, &workspace, &tmp, &info);
	lwork = (int)workspace;
	work = (double*)malloc(sizeof(double)*lwork);
	if(work == NULL){
		cerr<<"work malloc error"<<endl;
		exit(2);
	}
	
	/*---------Call Cell LAPACK library---------*/
	dgetri_(&ORDER_OF_MATRIX, a, &ORDER_OF_MATRIX, ipiv, work, &lwork, &info);
	if( info != 0 ){
		cerr<<"Call dgetri error"<<endl;
		exit(3);
	}
	for(i=0;i<ORDER_OF_MATRIX;i++){
		for(j=0;j<ORDER_OF_MATRIX;j++){
			INVERSE_MATRIX_A[i][j]=a[i+j*ORDER_OF_MATRIX];
		}
	}
	cout<<"Matrix inversion completed !"<<endl;
	free(ipiv);
	free(a);
	posix_memalign((void **)(&inverse_a_tmp), 128, sizeof(double)*ORDER_OF_MATRIX*ORDER_OF_MATRIX);
	if (inverse_a_tmp == NULL) {
		cerr<<"ERROR: malloc error"<<endl;
		exit(1);
	}
	for(i=0;i<ORDER_OF_MATRIX;i++){
		for(j=0;j<ORDER_OF_MATRIX;j++){
			inverse_a_tmp[i+j*ORDER_OF_MATRIX]=INVERSE_MATRIX_A[i][j];
		}
	}
	(inverse_a) = inverse_a_tmp;
	
	
	

	//~ vector < vector < double > > IDENTITY_MATRIX_A; 	


	//~ aux_double_vec.clear();
	//~ aux_double=0.00;	
	//~ for(int row=0; row<ORDER_OF_MATRIX; row++){
		//~ aux_double_vec.clear();
		//~ for(int col=0; col<ORDER_OF_MATRIX; col++){
			//~ aux_double_vec.push_back(aux_double);
		//~ }
		//~ IDENTITY_MATRIX_A.push_back(aux_double_vec);
		//~ aux_double_vec.clear();
	//~ }
	//~ Square_Matrix_Mult(MATRIX_A, INVERSE_MATRIX_A, IDENTITY_MATRIX_A);
	//~ for(int row=0; row<ORDER_OF_MATRIX; row++){

		//~ for(int col=0; col<ORDER_OF_MATRIX; col++){
			//~ if(row==col){
				//~ if(fabs(IDENTITY_MATRIX_A.at(row).at(col)-double(1.00))>1e-13){
					//~ cout << "IDENTITY_MATRIX_A["<<row<<"]["<<col<<"] is not 1!!! --- IDENTITY_MATRIX_A["<<row<<"]["<<col<<"]: "<<IDENTITY_MATRIX_A.at(row).at(col) << endl;
					//~ exit(1001);
				//~ }
			//~ }
			//~ else{
				//~ if(fabs(IDENTITY_MATRIX_A.at(row).at(col))>1e-13){
					//~ cout << "IDENTITY_MATRIX_A["<<row<<"]["<<col<<"] is not 0!!! --- IDENTITY_MATRIX_A["<<row<<"]["<<col<<"]: "<<IDENTITY_MATRIX_A.at(row).at(col) << endl;
					//~ exit(1001);
				//~ }
			//~ }
		//~ }
	//~ }
	//~ for(int row=0; row<ORDER_OF_MATRIX; row++){
		//~ for(int col=0; col<ORDER_OF_MATRIX; col++){
			//~ if(INVERSE_MATRIX_A.at(row).at(col)!=INVERSE_MATRIX_A.at(row).at(col) || isinf(INVERSE_MATRIX_A.at(row).at(col)) || isnan(INVERSE_MATRIX_A.at(row).at(col))){
				//~ cout << "cacchio!" << endl;
			//~ }
		//~ }
	//~ }
	
	return;
}

void adjust_surfaces(){
	
	for(int i=0; i<subSurfaces.size(); i++){
		for(int j=0; j<surfaces.size(); j++){
			double distance=get_distance(surfaces.at(i).center.x, surfaces.at(i).center.y, surfaces.at(i).center.z, surfaces.at(j).center.x, surfaces.at(j).center.y, surfaces.at(j).center.z);	
			if(i!=j && distance<1){
				cout << i <<" " <<j <<" " <<distance<<endl;
			}
		}
		
		
		for(int j=0; j<subSurfaces.at(i).size(); j++){
			subSurfaces.at(i).at(j).center.x=1e-12*subSurfaces.at(i).at(j).center.x;
			subSurfaces.at(i).at(j).center.y=1e-12*subSurfaces.at(i).at(j).center.y;
			subSurfaces.at(i).at(j).center.z=1e-12*subSurfaces.at(i).at(j).center.z;
		}
		
		surfaces.at(i).center.x=1e-12*surfaces.at(i).center.x;
		surfaces.at(i).center.y=1e-12*surfaces.at(i).center.y;
		surfaces.at(i).center.z=1e-12*surfaces.at(i).center.z;
	}
	
	return;
}

void create_membrane_tiles(){
	
	cout << "create_membrane_surfaces............"<<endl;
	surfaces.clear();
	
	
	double tcl=PRM.MEMBRANE_WIDTH;
	for(int i=0; i<=tcl; i++){
		limits.push_back(0.00);
	}
	
	int sqrtNoD=sqrt(int(PRM.NUM_OF_DIV));
	PRM.NUM_OF_DIV=sqrtNoD*sqrtNoD;
	
	int sqrtNoDSub=int(PRM.NUM_OF_SUB_DIV);
	PRM.NUM_OF_SUB_DIV=sqrtNoDSub*sqrtNoDSub;
	
	cout << "sqrtNoD: " << sqrtNoD <<endl; 
	cout << "PRM.NUM_OF_DIV: " << PRM.NUM_OF_DIV <<endl; 
	
	cout << "sqrtNoDSub: " << sqrtNoDSub <<endl; 
	cout << "PRM.NUM_OF_SUB_DIV: " << PRM.NUM_OF_SUB_DIV <<endl; 
	
	
	
	double dx=PRM.SIM_DOMAIN_WIDTH_X/double(sqrtNoD);
	double dy=PRM.SIM_DOMAIN_WIDTH_Y/double(sqrtNoD);
	
	double dxSub=dx/double(sqrtNoDSub);
	double dySub=dy/double(sqrtNoDSub);
	
	
	
	//LEFT MEMBRANE
	for(int ix=0; ix<sqrtNoD; ix++){
		
		for(int iy=0; iy<sqrtNoD; iy++){
			
			Surface surface;
				
			surface.center.x=double(0.5)*dx+PRM.MIN_X+double(ix)*dx;
			surface.center.y=double(0.5)*dy+PRM.MIN_Y+double(iy)*dy;
			surface.center.z=PRM.Z_MOUTH_LEFT;
			
			//===================== surface vertexes
			surface.pA.x=surface.center.x-double(0.5)*dx;
			surface.pA.y=surface.center.y-double(0.5)*dy;
			surface.pA.z=PRM.Z_MOUTH_LEFT;
			
			surface.pB.x=surface.center.x+double(0.5)*dx;
			surface.pB.y=surface.center.y-double(0.5)*dy;
			surface.pB.z=PRM.Z_MOUTH_LEFT;
			
			surface.pC.x=surface.center.x+double(0.5)*dx;
			surface.pC.y=surface.center.y+double(0.5)*dy;
			surface.pC.z=PRM.Z_MOUTH_LEFT;
			
			surface.pD.x=surface.center.x-double(0.5)*dx;
			surface.pD.y=surface.center.y+double(0.5)*dy;
			surface.pD.z=PRM.Z_MOUTH_LEFT;
		
			
				
			
			//===================== surface area
			surface.area=dx*dy;
				
				
			//===================== surface normal
			surface.normal[0]=0.00;
			surface.normal[1]=0.00;
			surface.normal[2]=1.00;
		
			//===================== subSurfaces
			vector <Surface> subSurfaces_of_current_surface;
			subSurfaces_of_current_surface.clear();	
			
			for(int ixSub=0; ixSub<sqrtNoDSub; ixSub++){
		
				for(int iySub=0; iySub<sqrtNoDSub; iySub++){
					Surface subsurface;
					
					
					subsurface.center.x=double(0.5)*dxSub+surface.pA.x+double(ix)*dxSub;
					subsurface.center.y=double(0.5)*dySub+surface.pA.y+double(iy)*dySub;
					subsurface.center.z=PRM.Z_MOUTH_LEFT;
					
					//===================== subsurface vertexes
					subsurface.pA.x=subsurface.center.x-double(0.5)*dxSub;
					subsurface.pA.y=subsurface.center.y-double(0.5)*dySub;
					subsurface.pA.z=PRM.Z_MOUTH_LEFT;
					
					subsurface.pB.x=subsurface.center.x+double(0.5)*dxSub;
					subsurface.pB.y=subsurface.center.y-double(0.5)*dySub;
					subsurface.pB.z=PRM.Z_MOUTH_LEFT;
					
					subsurface.pC.x=subsurface.center.x+double(0.5)*dxSub;
					subsurface.pC.y=subsurface.center.y+double(0.5)*dySub;
					subsurface.pC.z=PRM.Z_MOUTH_LEFT;
					
					subsurface.pD.x=subsurface.center.x-double(0.5)*dxSub;
					subsurface.pD.y=subsurface.center.y+double(0.5)*dySub;
					subsurface.pD.z=PRM.Z_MOUTH_LEFT;
					
					
					//===================== subsubsurface area
					subsurface.area=dxSub*dySub;
					
					//===================== subsubsurface normal
					subsurface.normal[0]=0.00;
					subsurface.normal[1]=0.00;
					subsurface.normal[2]=1.00;
					
					subSurfaces_of_current_surface.push_back(subsurface);    
					
				}
			}
			
			surfaces.push_back(surface);
			subSurfaces.push_back(subSurfaces_of_current_surface);
			subSurfaces_of_current_surface.clear();
		}
	}
	
	
	//RIGHT MEMBRANE
	for(int ix=0; ix<sqrtNoD; ix++){
		
		for(int iy=0; iy<sqrtNoD; iy++){
			
			Surface surface;
				
			surface.center.x=double(0.5)*dx+PRM.MIN_X+double(ix)*dx;
			surface.center.y=double(0.5)*dy+PRM.MIN_Y+double(iy)*dy;
			surface.center.z=PRM.Z_MOUTH_RIGHT;
			
			//===================== surface vertexes
			surface.pA.x=surface.center.x-double(0.5)*dx;
			surface.pA.y=surface.center.y-double(0.5)*dy;
			surface.pA.z=PRM.Z_MOUTH_RIGHT;
			
			surface.pB.x=surface.center.x+double(0.5)*dx;
			surface.pB.y=surface.center.y-double(0.5)*dy;
			surface.pB.z=PRM.Z_MOUTH_RIGHT;
			
			surface.pC.x=surface.center.x+double(0.5)*dx;
			surface.pC.y=surface.center.y+double(0.5)*dy;
			surface.pC.z=PRM.Z_MOUTH_RIGHT;
			
			surface.pD.x=surface.center.x-double(0.5)*dx;
			surface.pD.y=surface.center.y+double(0.5)*dy;
			surface.pD.z=PRM.Z_MOUTH_RIGHT;
		
			
				
			
			//===================== surface area
			surface.area=dx*dy;
				
				
			//===================== surface normal
			surface.normal[0]=0.00;
			surface.normal[1]=0.00;
			surface.normal[2]=-1.00;
		
			//===================== subSurfaces
			vector <Surface> subSurfaces_of_current_surface;
			subSurfaces_of_current_surface.clear();	
			
			for(int ixSub=0; ixSub<sqrtNoDSub; ixSub++){
		
				for(int iySub=0; iySub<sqrtNoDSub; iySub++){
					Surface subsurface;
					
					
					subsurface.center.x=double(0.5)*dxSub+surface.pA.x+double(ix)*dxSub;
					subsurface.center.y=double(0.5)*dySub+surface.pA.y+double(iy)*dySub;
					subsurface.center.z=PRM.Z_MOUTH_RIGHT;
					
					//===================== subsurface vertexes
					subsurface.pA.x=subsurface.center.x-double(0.5)*dxSub;
					subsurface.pA.y=subsurface.center.y-double(0.5)*dySub;
					subsurface.pA.z=PRM.Z_MOUTH_RIGHT;
					
					subsurface.pB.x=subsurface.center.x+double(0.5)*dxSub;
					subsurface.pB.y=subsurface.center.y-double(0.5)*dySub;
					subsurface.pB.z=PRM.Z_MOUTH_RIGHT;
					
					subsurface.pC.x=subsurface.center.x+double(0.5)*dxSub;
					subsurface.pC.y=subsurface.center.y+double(0.5)*dySub;
					subsurface.pC.z=PRM.Z_MOUTH_RIGHT;
					
					subsurface.pD.x=subsurface.center.x-double(0.5)*dxSub;
					subsurface.pD.y=subsurface.center.y+double(0.5)*dySub;
					subsurface.pD.z=PRM.Z_MOUTH_RIGHT;
					
					
					//===================== subsubsurface area
					subsurface.area=dxSub*dySub;
					
					//===================== subsubsurface normal
					subsurface.normal[0]=0.00;
					subsurface.normal[1]=0.00;
					subsurface.normal[2]=-1.00;
					
					subSurfaces_of_current_surface.push_back(subsurface);    
					
				}
			}
			
			surfaces.push_back(surface);
			subSurfaces.push_back(subSurfaces_of_current_surface);
			subSurfaces_of_current_surface.clear();
		}
	}
	
	return;
}
