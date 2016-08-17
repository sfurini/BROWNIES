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

#ifndef UTILS_H
#define UTILS_H


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
#include <boost/random/linear_congruential.hpp> 
#include <boost/random/lagged_fibonacci.hpp>	
#include <boost/random/uniform_int.hpp>		
#include <boost/random/uniform_real.hpp>	
#include <boost/random/normal_distribution.hpp>	
#include <boost/random/variate_generator.hpp>	

#include "constants.h"
#include "classes.h"
#include "physics_functions.h"


typedef boost::lagged_fibonacci19937 base_generator_type; // RANDOM NUMBER GENERATOR
typedef boost::normal_distribution<> base_distribution_type; // DISTRIBUTION RANDOM NUMBER
base_generator_type generator((unsigned)time(NULL)); // Definizione generatore di numeri casuali
base_distribution_type dist(0,1); // Definizione distribuzione
boost::variate_generator<base_generator_type&, base_distribution_type> uni(generator, dist);   

using namespace std;

string convertInt(int number);

string convertDouble(double number);






//gets the distance between two 3-D points
inline double get_distance(Point p1, Point p2){
	
	double distance;
	
	double dx=p1.x-p2.x;
	double dy=p1.y-p2.y;
	double dz=p1.z-p2.z;
	
	double dx_2=dx*dx;
	double dy_2=dy*dy;
	double dz_2=dz*dz;
	
	double distance_2=dx_2+dy_2+dz_2;
	
	distance=sqrt(distance_2);
	
	return distance;
}

//gets the distance between two 2-D points
inline double get_distance(RadialPoint rp1, RadialPoint rp2){
	
	double distance;
	
	double dz=rp1.z-rp2.z;
	double dr=rp1.r-rp2.r;
	
	double dz_2=dz*dz;
	double dr_2=dr*dr;
	
	double distance_2=dz_2+dr_2;
	
	distance=sqrt(distance_2);
	
	return distance;
}

//gets the distance between two charges
inline double get_distance(Charge c1, Charge c2){
	
	double distance;
	
	double dx=c1.x-c2.x;
	double dy=c1.y-c2.y;
	double dz=c1.z-c2.z;
	
	double dx_2=dx*dx;
	double dy_2=dy*dy;
	double dz_2=dz*dz;
	
	double distance_2=dx_2+dy_2+dz_2;
	
	distance=sqrt(distance_2);
	
	return distance;
}

//gets the distance between two ions
inline double get_distance(Ion i1, Ion i2){
	
	double distance;
	
	double dx=i1.x-i2.x;
	double dy=i1.y-i2.y;
	double dz=i1.z-i2.z;
	
	double dx_2=dx*dx;
	double dy_2=dy*dy;
	double dz_2=dz*dz;
	
	double distance_2=dx_2+dy_2+dz_2;
	
	distance=sqrt(distance_2);
	
	return distance;
}

//gets the distance between a point and an ion
inline double get_distance(Point p1, Ion i2){
	
	double distance;
	
	double dx=p1.x-i2.x;
	double dy=p1.y-i2.y;
	double dz=p1.z-i2.z;
	
	double dx_2=dx*dx;
	double dy_2=dy*dy;
	double dz_2=dz*dz;
	
	double distance_2=dx_2+dy_2+dz_2;
	
	distance=sqrt(distance_2);
	
	return distance;
}

//gets the distance between a point and a charge
inline double get_distance(Point p1, Charge c2){
	
	double distance;
	
	double dx=p1.x-c2.x;
	double dy=p1.y-c2.y;
	double dz=p1.z-c2.z;
	
	double dx_2=dx*dx;
	double dy_2=dy*dy;
	double dz_2=dz*dz;
	
	double distance_2=dx_2+dy_2+dz_2;
	
	distance=sqrt(distance_2);
	
	return distance;
}

//gets the distance between two entities
inline double get_distance(double x1, double y1, double z1, double x2, double y2, double z2){
	
	double distance;
	
	double dx=x1-x2;
	double dy=y1-y2;
	double dz=z1-z2;
	
	double dx_2=dx*dx;
	double dy_2=dy*dy;
	double dz_2=dz*dz;
	
	double distance_2=dx_2+dy_2+dz_2;
	
	distance=sqrt(distance_2);
	
	return distance;
}







//copy an entity to another (source --> destination)
void copy_point(Point ps, Point& pd);
void copy_radialpoint(RadialPoint rps, RadialPoint& rpd);
void copy_charge(Charge cs, Charge& cd);
void copy_ion(Ion is, Ion& id);

inline double gauss_cut(){
	double tmp=uni();
	if(tmp>double(3.00)){ 
		return double(3.00);
	}
	else if(tmp<double(-3.00)){ 
		return double(-3.00);
	}
	else{
		return tmp;
	}
}

void apply_periodic_boundary(double& x_fix, double& y_fix, double& z_fix, double& x_mv, double& y_mv, double& z_mv);



void apply_periodic_boundary_left_cell(double& x_fix, double& y_fix, double& z_fix, double& x_mv, double& y_mv, double& z_mv);
void apply_periodic_boundary_right_cell(double& x_fix, double& y_fix, double& z_fix, double& x_mv, double& y_mv, double& z_mv);



inline double sign(double value){
	
	return value/fabs(value);
}



string getNumberWithScale(double numero);

#endif
