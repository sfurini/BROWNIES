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


#ifndef CONSTANTS_H
#define CONSTANTS_H

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

#include "constants.h"
#include "utils.h"
#include "classes.h"
#include "physics_functions.h"

using namespace std;

const double zero = 1e-20;
const double INF = 1e37;
const double degToRad=2.0*M_PI/360.0;
const double radToDeg=360.0/(2.0*M_PI);


const double Q = 1.602176487e-19;			// [C] elementary charge
const double EPS_0=8.854187817e-12;			// [F/m] vacuum permitivity
const double COULOMB_K=8.987551787e9;		// [N*m^2/C^2] Coulomb'constant 1/(4*PI*EPS_0)
const double AVOGADRO=6.0221415e23;      	// [#] Avogadro' number
const double BOLTZMANN_K=1.3806504e-23;	// [J/K] Boltzmann's constant
const double FARADAY_K=9.64853399e4;	// [C/mol] Faraday's constant

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

const int NUM_OF_IONIC_SPECIES=92;




#endif
