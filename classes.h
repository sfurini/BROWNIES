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

#ifndef CLASSES_H
#define CLASSES_H


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
#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw

//~ #include <include/constants.h>
//~ #include <include/utils.h>
//~ #include <include/classes.h>

//~ #include "physics_functions.h"
//~ #include "physics_functions.cc"

using namespace std;

//########################################
// Class Point
//########################################
// point in a 3-D space
// with physical properties
// position is in pm
// potential is in Volt
// field is in V/m
// energy is in J
class Point{
	public:
		
	double x;
	double y;
	double z;

	double potential;
	double energy;
	
	double field[3];
	
	Point();
	~Point();
	void reset_point();
	
	friend ostream& operator<<(ostream& stream, Point& p);
	
};		


//########################################
// Class RadialPoint
//########################################
// point in a 2-D space
// with physical properties
// position is in pm
// potential is in Volt
// field is in V/m
// energy is in J
class RadialPoint{
	public:
		
	double z;
	double r;

	double potential;
	double energy;
	
	double field[2];
	
	RadialPoint();
	~RadialPoint();
	void reset_radialpoint();
	
	friend ostream& operator<<(ostream& stream, RadialPoint& rp);
};	

//########################################
// Class Charge
//########################################
// point in a 3-D space
// charge is in Coulomb

class Charge{
	public:

	double x;
	double y;
	double z;
	double charge;
	double valence;
	double DW_charge;
	double DW_valence;
	double radius;
	
	Charge();
	~Charge();
	void reset_charge();
	
	friend ostream& operator<<(ostream& stream, Charge& c);
};


//########################################
// Class Parameters
//########################################
class Parameters{
	public:
		
	int HISTORY_SIZE;
	
	unsigned int SEED;
		
	string SIM_TYPE;

	string PREFIX;
	
	string SR_METHOD;
	
	double MD_MAP_MIN_Z;
	double MD_MAP_MAX_Z;
	
	string PBC;
	string ION_RECYCLING;
	
	string FIXED_CHARGES_FILE;
	
	double CONTROL_CELL_WIDTH;
	double BATH_WIDTH;
	double OUTER_REGION_WIDTH;
	double MEMBRANE_WIDTH;
	double SIM_DOMAIN_WIDTH_X;
	double SIM_DOMAIN_WIDTH_Y;
	double SIM_DOMAIN_WIDTH_Z;
	double MIN_X;
	double MIN_Y;
	double MIN_Z;
	double MAX_X;
	double MAX_Y;
	double MAX_Z;
	
	double SIM_DOMAIN_VOLUME;
	double LEFT_CELL_MIN_Z;
	double LEFT_CELL_MAX_Z;
	double RIGHT_CELL_MIN_Z;
	double RIGHT_CELL_MAX_Z;
	double CONTROL_CELL_VOLUME;
	
	double Z_BEGIN_SR_MAP;
	double Z_END_SR_MAP;
	double DELTA_SR_MAP;

	double APPLIED_POTENTIAL;
	double APPLIED_FIELD;
	
	double c_0;	
	double CONC_LEFT_KCL;
	double CONC_LEFT_NACL;
	double CONC_LEFT_CACL2;
	double CONC_LEFT_MGCL2;
	double CONC_RIGHT_KCL;	
	double CONC_RIGHT_NACL;
	double CONC_RIGHT_CACL2; 
	double CONC_RIGHT_MGCL2; 
	
	double SHORT_RANGE_EXP;
	
	double EPS_W;
	double EPS_MEM;
	
	bool CONSTANT_DIFF_COEFF;
	double DIFF_COEFF_IN_CHANNEL;
	
	double delta_epsilon;
	double mean_epsilon;
	double dielectrics_weight;
	
	double LEFT_VESTIBULE_CURVATURE_RADIUS;
	double LEFT_VESTIBULE_MIN_CHANNEL_RADIUS;
	
	double RIGHT_VESTIBULE_CURVATURE_RADIUS;
	double RIGHT_VESTIBULE_MIN_CHANNEL_RADIUS;
	
	int TILES_PER_RING;
	int NUM_OF_DIV;
	int NUM_OF_SUB_DIV;
	double MAX_TILE_WIDTH;
	
	int NUMBER_OF_TILES;
	int NUMBER_OF_SUBTILES_PER_TILE;
	
	double Z_MOUTH_LEFT;
	double Z_MOUTH_RIGHT;
	
	double DELTA_T;
	double TEMPERATURE;
	double kT;
	
	double MAX_VEL;
	double MAX_FLIGHT;
	
	double PREP_STEPS;
	double SIM_STEPS;
	double FIRST_STEP;
	
	int STATS_OUT_FREQ;
	double STATS_DZ;
	
	bool channel_pdb_files;
	bool trajectory;
	bool flux;
	bool rdf;
	bool vel_distribution;
	bool mean_square_displ;
	bool induced_charge;	
	int potential;
	int concentrations;
	bool channel_configuration;
	bool filter_configuration;
	
	bool currents_ZT;

	int long_flight;
	int into_membrane;
	int out_of_SF;
	int out_of_domain;
	
	double BOX1_MIN_X;
	double BOX1_MAX_X;
	double BOX1_MIN_Y;
	double BOX1_MAX_Y;
	double BOX1_MIN_Z;
	double BOX1_MAX_Z;
	string BOX1_IS1;
	string BOX1_IS2;
	string BOX1_IS3;
	string BOX1_IS4;
	string BOX1_IS5;
	string BOX1_IS6;
	int BOX1_N1;
	int BOX1_N2;
	int BOX1_N3;
	int BOX1_N4;
	int BOX1_N5;
	int BOX1_N6;
							
	double BOX2_MIN_X;
	double BOX2_MAX_X;
	double BOX2_MIN_Y;
	double BOX2_MAX_Y;
	double BOX2_MIN_Z;
	double BOX2_MAX_Z;
	string BOX2_IS1;
	string BOX2_IS2;
	string BOX2_IS3;
	string BOX2_IS4;
	string BOX2_IS5;
	string BOX2_IS6;
	int BOX2_N1;
	int BOX2_N2;
	int BOX2_N3;
	int BOX2_N4;
	int BOX2_N5;
	int BOX2_N6;
						
	double BOX3_MIN_X;
	double BOX3_MAX_X;
	double BOX3_MIN_Y;
	double BOX3_MAX_Y;
	double BOX3_MIN_Z;
	double BOX3_MAX_Z;
	string BOX3_IS1;
	string BOX3_IS2;
	string BOX3_IS3;
	string BOX3_IS4;
	string BOX3_IS5;
	string BOX3_IS6;
	int BOX3_N1;
	int BOX3_N2;
	int BOX3_N3;
	int BOX3_N4;
	int BOX3_N5;
	int BOX3_N6;
						
	double BOX4_MIN_X;
	double BOX4_MAX_X;
	double BOX4_MIN_Y;
	double BOX4_MAX_Y;
	double BOX4_MIN_Z;
	double BOX4_MAX_Z;
	string BOX4_IS1;
	string BOX4_IS2;
	string BOX4_IS3;
	string BOX4_IS4;
	string BOX4_IS5;
	string BOX4_IS6;
	int BOX4_N1;
	int BOX4_N2;
	int BOX4_N3;
	int BOX4_N4;
	int BOX4_N5;
	int BOX4_N6;
						
	double BOX5_MIN_X;
	double BOX5_MAX_X;
	double BOX5_MIN_Y;
	double BOX5_MAX_Y;
	double BOX5_MIN_Z;
	double BOX5_MAX_Z;
	string BOX5_IS1;
	string BOX5_IS2;
	string BOX5_IS3;
	string BOX5_IS4;
	string BOX5_IS5;
	string BOX5_IS6;
	int BOX5_N1;
	int BOX5_N2;
	int BOX5_N3;
	int BOX5_N4;
	int BOX5_N5;
	int BOX5_N6;
					
	double BOX6_MIN_X;
	double BOX6_MAX_X;
	double BOX6_MIN_Y;
	double BOX6_MAX_Y;
	double BOX6_MIN_Z;
	double BOX6_MAX_Z;
	string BOX6_IS1;
	string BOX6_IS2;
	string BOX6_IS3;
	string BOX6_IS4;
	string BOX6_IS5;
	string BOX6_IS6;
	int BOX6_N1;
	int BOX6_N2;
	int BOX6_N3;
	int BOX6_N4;
	int BOX6_N5;
	int BOX6_N6;
	
	vector <RadialPoint> channel_profile_points;
	
	vector <double> charge_ring_z;
	vector <double> charge_ring_r;
	vector <int> charge_ring_n;
	vector <double> charge_ring_q;
	
	vector <bool> ions_to_simulate;
	
	Parameters();
	~Parameters();
	void reset_parameters();
	
	friend ostream& operator<<(ostream& stream, Parameters& PRM);
	
};


//########################################
// Class Ion
//########################################
// charge
// physical properties

class Ion{
	
	public:
		
		static int IDENTIFIER;
//===========================================	
	int kind;	
	string name;									/*

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

	double charge;
	double valence;
	double DW_charge;
	double DW_valence;
	
	double mass;
	double radius;
	double pm_radius;
	double charmm_half_radius; // 
	double charmm_eps; // 
	double diffusion_coeff;

	double EPS_ION;

	double x;
	double y;
	double z;
	
	double x_prev;
	double y_prev;
	double z_prev;
	
	double x_next;
	double y_next;
	double z_next;
	
	double field[3];
	double velocity[3];
	double force[3];
	double force_old[3];
	double X_n[3];
	double X_n_old[3];
	double X_n_minus[3];
	
	double K_pos_1;
	double K_pos_2;
	double K_pos_3;
	double K_pos_4;
	double K_pos_5;
	double K_vel_1;
	double K_vel_2;
	double K_vel_3;
	double K_Y;
	double K_step0_1;
	double K_step0_2;
	double K_X;
	
	double C;
	double G;
	double E;
	double H;

	double bro_KF;
	double bro_KR;
	
	int side;
	
	Ion();
	Ion(int type);
	~Ion();
	void reset_ion();
	void set_brownian_parameters();
	
	friend ostream& operator<<(ostream& stream, Ion& i);
};

//########################################
// Class Surface
//########################################
class Surface{
	
	public:
	Point pA;
	Point pB;
	Point pC;
	Point pD;
	Point center;
	double normal[3];
	double area;

	Surface();
	~Surface();
	void reset_surface();

	friend ostream& operator<<(ostream& stream, Surface& s);
};

//########################################
// Class SurfaceLight
//########################################
class SurfaceLight{
	
	public:
	double center_x;
	double center_y;
	double center_z;
	double normal[3];
	double area;

	SurfaceLight();
	~SurfaceLight();
	void reset_surfaceLight();

	friend ostream& operator<<(ostream& stream, SurfaceLight& s);
};
//########################################
// Class Tile
//########################################
class Tile{
	
	public:
	Point pA;
	Point pB;
	Point pC;
	Point pD;
	Point center;
	double normal[3];
	double area;
	double P_B;

	Tile();
	~Tile();
	void reset_Tile();
	
	friend ostream& operator<<(ostream& stream, Tile& t);
};


//########################################
// Class Statistics
//########################################
class Statistics{
	public:
	//~ int step;
	int MAX_STEPS;
	int num_of_dz;
	double DELTA_Z;
	vector <int> ions_this_step;
	
// flux computation	
	vector < vector < double > > concentrations_along_z;
	
// radial distribution function computation
	vector < vector < vector < double > > > RDF;
	vector < double > RDF_samples;
	
// velocity distribution computation	
	vector < vector < double > > velocities;
	vector < vector < double > > velocities_theory;
	double velocities_samples;
	
// mean square displacement computation	
	int msds_index_1;
	int msds_index_2;
	vector < vector < double > > msds;	
	vector < double > msds_num;	
	vector <Point> msds_start;	
	
// induced charge computation
	double A_total_charge;
	vector <double> C_total_charge;
	vector <double> C_total_charge_type;
		
// concentrations computation
	vector < vector < vector < vector < double > > > > concs_3D;
	vector < vector < vector <  double > > > radial_concs;
			
// potential computation
	vector < double > average_potentials_on_axis;
	vector < vector < double > > potentials_on_axis;
	vector < vector < double > > average_radial_potentials;
	vector < vector < vector <  double > > > radial_potentials;
	
// currents computation	(RAMO-SCHOCKLEY THEOREM)	
	vector <double> currents_RS;	
	vector <double> instant_currents_RS;	
	double *vector_h_RAMO;	
	
// currents computation	(ZERO THRESHOLD, one threshold at z=0)
	vector <double> currents_ZT;
	vector <double> instant_currents_ZT;


	string stat_file_1;
	string stat_file_2;
	string stat_file_3;
	string stat_file_4;
	string stat_file_5;
	string stat_file_6;
	string stat_file_7;
	string stat_file_8;
	string stat_file_9;
	string stat_file_A;
	string stat_file_B;
	string stat_file_C;
	string stat_file_D;
	
	string gnuplot_file_1;
	string gnuplot_file_2;
	string gnuplot_file_3;
	string gnuplot_file_4;
	string gnuplot_file_5;
	string gnuplot_file_6;
	string gnuplot_file_7;
	string gnuplot_file_8;
	string gnuplot_file_9;
	string gnuplot_file_A;
	string gnuplot_file_B;
	string gnuplot_file_C;
	string gnuplot_file_D;
	
	string save_status_file;
	
	
	double ERROR_nan;
	double ERROR_long_jump;
	double ERROR_into_membrane;
	double ERROR_ion_boxes;


	Statistics();
	~Statistics();
	void reset_statistics();
	void update_statistics();
	void compute_currents_RS();
	void compute_currents_ZT();
	
	void print_statistics();
	
	void reset_after_restarting_from_the_beginning();
	
	friend ostream& operator<<(ostream& stream, Statistics& STAT);
	
};

//########################################
// Class Control_cell
//########################################
class Control_cell{
	public:
		
	int SIDE;

	vector <double> eta;
	vector <double> SumEta;	
	vector <double> SumNi;
	vector <double> Xi;
	
	vector <bool> in_control_cell;
	
	Control_cell();
	~Control_cell();
	void reset_control_cell();
	
	friend ostream& operator<<(ostream& stream, Control_cell& cell);
	
};
	


	



//########################################
// Class Box
//########################################
// simulation box
// limits, dimensions, permittivity
class Box{
	public:
		
	double min_x;
	double min_y;
	double min_z;
	double max_x;
	double max_y;
	double max_z;
	double volume_dimension_x;
	double volume_dimension_y;
	double volume_dimension_z;
	
	double permittivity;


	Box();
	~Box();
	void reset_box();
	
	friend ostream& operator<<(ostream& stream, Box& box);
	
};



//########################################
// Class Ion_Box
//########################################
// ion box
// limits, dimensions, 
class Ion_box{
	public:
		
	int index;
		
	//~ double MIN_X;
	//~ double MIN_Y;
	double MIN_Z;
	//~ double MAX_X;
	//~ double MAX_Y;
	double MAX_Z;
	
	double TOT_CHARGE;
	
	vector <string> is;
	vector <int> in;
	vector <int> ion_indexes;

	Ion_box();
	~Ion_box();
	void reset_ion_box();
	
	friend ostream& operator<<(ostream& stream, Ion_box& ion_box);
	
};





#endif
