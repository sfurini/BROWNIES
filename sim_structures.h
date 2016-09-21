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

#ifndef SIM_STRUCTURES_H
#define SIM_STRUCTURES_H


Parameters PRM;

vector <Surface> surfaces;
vector < vector < Surface > > subSurfaces;

vector<RadialPoint> channelProfile;

vector <Charge> membrane_charges;



vector < vector < vector < double > > > FORCE_SR;
vector < vector < vector < double > > > ENERGY_SR;	
vector < vector < vector < double > > > FORCE_OTHER;
vector < vector < vector < double > > > ENERGY_OTHER;



double ENERGY_C[20000];
double POTENTIAL_C[20000];
double FIELD_C[20000];
double FORCE_C[20000];




vector < vector < vector < vector < double  > > > > SR_dielectric_boundary;

vector <double> limits;

vector < vector < double > > SR_ion_boxes;

vector <Ion_box> ion_boxes;	

vector < vector < double > > MATRIX_A;

double *vector_c, *vector_h;
double *vector_c_MD_map, *vector_h_MD_map;
double *inverse_a;

int ORDER_OF_MATRIX;
int matrix_size;

int _int_one;
double _double_one;
double _double_zero;
char _N;

Statistics STAT;

vector < vector < double > > diffusion_coefficients;

vector < vector < Ion > > ions;

Ion sample_ions[100];

vector <Control_cell> left_cell;
vector <Control_cell> right_cell;


double step;
vector <double> step_window;

int left_LFother=0;
int left_LBother=0;
int left_RFother=0;
int left_RBother=0;
int left_LFthis=0;
int left_LBthis=0;
int left_RFthis=0;
int left_RBthis=0;
int left_LFtotal=0;
int left_LBtotal=0;
int left_RFtotal=0;
int left_RBtotal=0;

int right_LFother=0;
int right_LBother=0;
int right_RFother=0;
int right_RBother=0;
int right_LFthis=0;
int right_LBthis=0;
int right_RFthis=0;
int right_RBthis=0;
int right_LFtotal=0;
int right_LBtotal=0;
int right_RFtotal=0;
int right_RBtotal=0;



int retrace_window;
	
double last_step_statistics;
double last_step_output;

double step_cap;

double local_step_cap1;
double local_step_cap2;
double local_step_cap3;
double local_step_cap4;
double local_step_cap5;

double attempts;		


Ion IONS[50][2000];

int INDEX_LAST_STEP;
int INDEX_STAT_STEP;

double STEPS[50];
int NUM_OF_IONS_IN_STEP[50];


SurfaceLight SURFACES[1200];
int NUM_OF_SURFACES=0;
double C_BASE[1200];


Control_cell LEFT_CELLS[50];
Control_cell RIGHT_CELLS[50];

int indexesGap;

int* test_arrayInt;
int* test_arrayArrayInt;

double POINTS_ON_AXIS[2001];
double POTENTIALS_ON_AXIS[50][2001];
double AVERAGE_POTENTIALS_ON_AXIS[2001];

double CHANNEL_CONFIGURATIONS[3][3][3][3][3];
double FILTER_CONFIGURATIONS[3][3][3][3][3];

double CHANNEL_CONFIGURATIONS_FREQUENCIES[3][3][3][3][3][100];
double FILTER_CONFIGURATIONS_FREQUENCIES[3][3][3][3][3][100];

int last_output_channel_configuration[5];
int last_output_filter_configuration[5];

int channel_same_configuration_counter=0;
int filter_same_configuration_counter=0;

int pdb_file_index=0;

#endif
