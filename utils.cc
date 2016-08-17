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
#include "sim_structures.h"

string convertInt(int number){
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

string convertDouble(double number){
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}


//copy a source Point into a destination Point
void copy_point(Point ps, Point& pd){
	
	pd.reset_point();
	
	pd.x=ps.x;
	pd.y=ps.y;
	pd.z=ps.z;
	
	pd.potential=ps.potential;
	pd.energy=ps.energy;
	
	pd.field[0]=ps.field[0];
	pd.field[1]=ps.field[1];
	pd.field[2]=ps.field[2];
	
}

//copy a source RadialPoint into a destination RadialPoint
void copy_radialpoint(RadialPoint rps, RadialPoint& rpd){
	rpd.reset_radialpoint();
	
	rpd.z=rps.z;
	rpd.r=rps.r;
	
	rpd.potential=rps.potential;
	rpd.energy=rps.energy;
	
	rpd.field[0]=rps.field[0];
	rpd.field[1]=rps.field[1];
	
}

//copy a source Charge into a destination Charge
void copy_charge(Charge cs, Charge& cd){
	cd.reset_charge();
	
	cd.x=cs.x;
	cd.y=cs.y;
	cd.z=cs.z;
	
	cd.charge=cs.charge;
	
}
	
void copy_ion(Ion is, Ion& id){
	
	id.kind=is.kind;
	
	id.name=is.name;
	id.charge=is.charge;
	id.valence=is.valence;
	id.DW_charge=is.DW_charge;
	id.DW_valence=is.DW_valence;
	
	
	id.mass=is.mass;
	id.radius=is.radius;
	id.pm_radius=is.pm_radius;
	id.charmm_half_radius=is.charmm_half_radius; // 
	id.charmm_eps=is.charmm_eps; // 
	id.diffusion_coeff=is.diffusion_coeff;
	id.EPS_ION=is.EPS_ION;
	
	id.x=is.x;
	id.y=is.y;
	id.z=is.z;
	
	id.x_prev=is.x_prev;
	id.y_prev=is.y_prev;
	id.z_prev=is.z_prev;
	
	id.x_next=is.x_next;
	id.y_next=is.y_next;
	id.z_next=is.z_next;

	id.field[0]=is.field[0];
	id.field[1]=is.field[1];
	id.field[2]=is.field[2];
	
	id.velocity[0]=is.velocity[0];
	id.velocity[1]=is.velocity[1];
	id.velocity[2]=is.velocity[2];
	
	id.force[0]=is.force[0];
	id.force[1]=is.force[1];
	id.force[2]=is.force[2];
	
	id.force_old[0]=is.force_old[0];
	id.force_old[1]=is.force_old[1];
	id.force_old[2]=is.force_old[2];
	
	id.X_n[0]=is.X_n[0];
	id.X_n[1]=is.X_n[1];
	id.X_n[2]=is.X_n[2];
	
	id.X_n_old[0]=is.X_n_old[0];
	id.X_n_old[1]=is.X_n_old[1];
	id.X_n_old[2]=is.X_n_old[2];
	
	id.X_n_minus[0]=is.X_n_minus[0];
	id.X_n_minus[1]=is.X_n_minus[1];
	id.X_n_minus[2]=is.X_n_minus[2];
	
	id.K_pos_1=is.K_pos_1;
	id.K_pos_2=is.K_pos_2;
	id.K_pos_3=is.K_pos_3;
	id.K_pos_4=is.K_pos_4;
	id.K_pos_5=is.K_pos_5;
	
	id.K_vel_1=is.K_vel_1;
	id.K_vel_2=is.K_vel_2;
	id.K_vel_3=is.K_vel_3;
	
	id.K_Y=is.K_Y;
	id.K_step0_1=is.K_step0_1;
	id.K_step0_2=is.K_step0_2;
	id.K_X=is.K_X;
	
	id.C=is.C;
	id.G=is.G;
	id.E=is.E;
	id.H=is.H;

	id.side=is.side;
	
	return;
}


void apply_periodic_boundary(double& x_fix, double& y_fix, double& z_fix, double& x_mv, double& y_mv, double& z_mv){

	double dist_x=x_fix-x_mv;
	double dist_y=y_fix-y_mv;
	double dist_z=z_fix-z_mv;
	
	if(PRM.PBC.substr(0, 2).compare("11")==0){
		if(dist_x>PRM.MAX_X){
			x_mv+=PRM.SIM_DOMAIN_WIDTH_X;
		}
		else if(dist_x<PRM.MIN_X){
			x_mv-=PRM.SIM_DOMAIN_WIDTH_X;
		}
		if(dist_y>PRM.MAX_Y){
			y_mv+=PRM.SIM_DOMAIN_WIDTH_Y;
		}
		else if(dist_y<PRM.MIN_Y){
			y_mv-=PRM.SIM_DOMAIN_WIDTH_Y;
		}
	}
	if(PRM.PBC.substr(2, 1).compare("1")==0){
		if(dist_z>PRM.MAX_Z){
			z_mv+=PRM.SIM_DOMAIN_WIDTH_Z;
		}
		else if(dist_z<PRM.MIN_Z){
			z_mv-=PRM.SIM_DOMAIN_WIDTH_Z;
		}
	}
	
	return;
}
void apply_periodic_boundary_left_cell(double& x_fix, double& y_fix, double& z_fix, double& x_mv, double& y_mv, double& z_mv){
	
	double dist_x=x_fix-x_mv;
	double dist_y=y_fix-y_mv;
	double dist_z=z_fix-z_mv;
	
	if(dist_x>0.5e-12*PRM.SIM_DOMAIN_WIDTH_X){
		x_mv+=1e-12*PRM.SIM_DOMAIN_WIDTH_X;
	}
	else if(dist_x<-0.5e-12*PRM.SIM_DOMAIN_WIDTH_X){
		x_mv-=1e-12*PRM.SIM_DOMAIN_WIDTH_X;
	}
	if(dist_y>0.5e-12*PRM.SIM_DOMAIN_WIDTH_Y){
		y_mv+=1e-12*PRM.SIM_DOMAIN_WIDTH_Y;
	}
	else if(dist_y<-0.5e-12*PRM.SIM_DOMAIN_WIDTH_Y){
		y_mv-=1e-12*PRM.SIM_DOMAIN_WIDTH_Y;
	}
	if(dist_z>0.5e-12*PRM.CONTROL_CELL_WIDTH){
		z_mv+=1e-12*PRM.CONTROL_CELL_WIDTH;
	}
	else if(dist_z<-0.5e-12*PRM.CONTROL_CELL_WIDTH){
		z_mv-=1e-12*PRM.CONTROL_CELL_WIDTH;
	}
    
	return;
}

void apply_periodic_boundary_right_cell(double& x_fix, double& y_fix, double& z_fix, double& x_mv, double& y_mv, double& z_mv){
	
	double dist_x=x_fix-x_mv;
	double dist_y=y_fix-y_mv;
	double dist_z=z_fix-z_mv;
	
	if(dist_x>0.5e-12*PRM.SIM_DOMAIN_WIDTH_X){
		x_mv+=1e-12*PRM.SIM_DOMAIN_WIDTH_X;
	}
	else if(dist_x<-0.5e-12*PRM.SIM_DOMAIN_WIDTH_X){
		x_mv-=1e-12*PRM.SIM_DOMAIN_WIDTH_X;
	}
	if(dist_y>0.5e-12*PRM.SIM_DOMAIN_WIDTH_Y){
		y_mv+=1e-12*PRM.SIM_DOMAIN_WIDTH_Y;
	}
	else if(dist_y<-0.5e-12*PRM.SIM_DOMAIN_WIDTH_Y){
		y_mv-=1e-12*PRM.SIM_DOMAIN_WIDTH_Y;
	}
	if(dist_z>0.5e-12*PRM.CONTROL_CELL_WIDTH){
		z_mv+=1e-12*PRM.CONTROL_CELL_WIDTH;
	}
	else if(dist_z<-0.5e-12*PRM.CONTROL_CELL_WIDTH){
		z_mv-=1e-12*PRM.CONTROL_CELL_WIDTH;
	}
    
	return;
}

string getNumberWithScale(double numero){
	string numeroConScala="";
	double max=1000;
	
	double ab=fabs(numero);
	ab=ab/max;
	
	ostringstream ss;
	
	if(ab>=double(1e21) && ab<double(1e24)){
		ss << numero*double(1e-24);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " Y";
	}	
	else if(ab>=double(1e18) && ab<double(1e21)){
		ss << numero*double(1e-21);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " Z";
	}	
	else if(ab>=double(1e15) && ab<double(1e18)){
		ss << numero*double(1e-18);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " E";
	}	
	else if(ab>=double(1e12) && ab<double(1e15)){
		ss << numero*double(1e-15);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " P";
	}	
	else if(ab>=double(1e9) && ab<double(1e12)){
		ss << numero*double(1e-12);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " T";
	}	
	else if(ab>=double(1e6) && ab<double(1e9)){
		ss << numero*double(1e-9);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " G";
	}	
	else if(ab>=double(1e3) && ab<double(1e6)){
		ss << numero*double(1e-6);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " M";
	}	
	else if(ab>=double(1) && ab<double(1e3)){
		ss << numero*double(1e-3);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " k";
	}	
	else if((ab>=double(1e-3) && ab<double(1)) || numero==0){
		ss << numero;
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " ";
	}	
	else if(ab>=double(1e-6) && ab<double(1e-3)){
		ss << numero*double(1e3);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " m";
	}
	else if(ab>=double(1e-9) && ab<double(1e-6)){
		ss << numero*double(1e6);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " u";
	} 
	else if(ab>=double(1e-12) && ab<double(1e-9)){
		ss << numero*double(1e9);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " n";
	} 
	else if(ab>=double(1e-15) && ab<double(1e-12)){
		ss << numero*double(1e12);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " p";
	} 
	else if(ab>=double(1e-18) && ab<double(1e-15)){
		ss << numero*double(1e15);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " f";
	} 
	else if(ab>=double(1e-21) && ab<double(1e-18)){
		ss << numero*double(1e18);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " a";
	} 
	else if(ab>=double(1e-24) && ab<double(1e-21)){
		ss << numero*double(1e21);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " z";
	} 
	else if(ab>=double(1e-27) && ab<double(1e-24)){
		ss << numero*double(1e24);
		numeroConScala=ss.str();
		numeroConScala=numeroConScala + " y";
	} 
	
	return numeroConScala;
}
