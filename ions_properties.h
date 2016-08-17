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

#ifndef IONS_PROPERTIES_H
#define IONS_PROPERTIES_H


double K_VALENCE;  /* for input.phase(page_6).number(K_VALENCE) */
double K_RADIUS;  /* for input.phase(page_6).number(K_RADIUS) */
double K_DIFF_COEFF;  /* for input.phase(page_6).number(K_DIFF_COEFF) */
double K_MASS;  /* for input.phase(page_6).number(K_MASS) */
double NA_VALENCE;  /* for input.phase(page_6).number(NA_VALENCE) */
double NA_RADIUS;  /* for input.phase(page_6).number(NA_RADIUS) */
double NA_DIFF_COEFF;  /* for input.phase(page_6).number(NA_DIFF_COEFF) */
double NA_MASS;  /* for input.phase(page_6).number(NA_MASS) */
double CA_VALENCE;  /* for input.phase(page_6).number(CA_VALENCE) */
double CA_RADIUS;  /* for input.phase(page_6).number(CA_RADIUS) */
double CA_DIFF_COEFF;  /* for input.phase(page_6).number(CA_DIFF_COEFF) */
double CA_MASS;  /* for input.phase(page_6).number(CA_MASS) */
double MG_VALENCE;  /* for input.phase(page_6).number(MG_VALENCE) */
double MG_RADIUS;  /* for input.phase(page_6).number(MG_RADIUS) */
double MG_DIFF_COEFF;  /* for input.phase(page_6).number(MG_DIFF_COEFF) */
double MG_MASS;  /* for input.phase(page_6).number(MG_MASS) */
double CL_VALENCE;  /* for input.phase(page_6).number(CL_VALENCE) */
double CL_RADIUS;  /* for input.phase(page_6).number(CL_RADIUS) */
double CL_DIFF_COEFF;  /* for input.phase(page_6).number(CL_DIFF_COEFF) */
double CL_MASS;  /* for input.phase(page_6).number(CL_MASS) */
    


double BOX1_ION1_VALENCE;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_1).number(BOX1_ION1_VALENCE) */
double BOX1_ION1_RADIUS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_1).number(BOX1_ION1_RADIUS) */
double BOX1_ION1_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_1).number(BOX1_ION1_DIFF_COEFF) */
double BOX1_ION1_MASS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_1).number(BOX1_ION1_MASS) */
double BOX1_ION2_VALENCE;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_2).number(BOX1_ION2_VALENCE) */
double BOX1_ION2_RADIUS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_2).number(BOX1_ION2_RADIUS) */
double BOX1_ION2_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_2).number(BOX1_ION2_DIFF_COEFF) */
double BOX1_ION2_MASS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_2).number(BOX1_ION2_MASS) */
double BOX1_ION3_VALENCE;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_3).number(BOX1_ION3_VALENCE) */
double BOX1_ION3_RADIUS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_3).number(BOX1_ION3_RADIUS) */
double BOX1_ION3_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_3).number(BOX1_ION3_DIFF_COEFF) */
double BOX1_ION3_MASS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_3).number(BOX1_ION3_MASS) */
double BOX1_ION4_VALENCE;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_4).number(BOX1_ION4_VALENCE) */
double BOX1_ION4_RADIUS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_4).number(BOX1_ION4_RADIUS) */
double BOX1_ION4_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_4).number(BOX1_ION4_DIFF_COEFF) */
double BOX1_ION4_MASS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_4).number(BOX1_ION4_MASS) */
double BOX1_ION5_VALENCE;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_5).number(BOX1_ION5_VALENCE) */
double BOX1_ION5_RADIUS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_5).number(BOX1_ION5_RADIUS) */
double BOX1_ION5_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_5).number(BOX1_ION5_DIFF_COEFF) */
double BOX1_ION5_MASS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_5).number(BOX1_ION5_MASS) */
double BOX1_ION6_VALENCE;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_6).number(BOX1_ION6_VALENCE) */
double BOX1_ION6_RADIUS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_6).number(BOX1_ION6_RADIUS) */
double BOX1_ION6_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_6).number(BOX1_ION6_DIFF_COEFF) */
double BOX1_ION6_MASS;  /* for input.phase(page_7).group(ion_box_1).group(ion_box_1_ion_6).number(BOX1_ION6_MASS) */
double BOX2_ION1_VALENCE;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_1).number(BOX2_ION1_VALENCE) */
double BOX2_ION1_RADIUS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_1).number(BOX2_ION1_RADIUS) */
double BOX2_ION1_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_1).number(BOX2_ION1_DIFF_COEFF) */
double BOX2_ION1_MASS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_1).number(BOX2_ION1_MASS) */
double BOX2_ION2_VALENCE;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_2).number(BOX2_ION2_VALENCE) */
double BOX2_ION2_RADIUS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_2).number(BOX2_ION2_RADIUS) */
double BOX2_ION2_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_2).number(BOX2_ION2_DIFF_COEFF) */
double BOX2_ION2_MASS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_2).number(BOX2_ION2_MASS) */
double BOX2_ION3_VALENCE;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_3).number(BOX2_ION3_VALENCE) */
double BOX2_ION3_RADIUS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_3).number(BOX2_ION3_RADIUS) */
double BOX2_ION3_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_3).number(BOX2_ION3_DIFF_COEFF) */
double BOX2_ION3_MASS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_3).number(BOX2_ION3_MASS) */
double BOX2_ION4_VALENCE;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_4).number(BOX2_ION4_VALENCE) */
double BOX2_ION4_RADIUS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_4).number(BOX2_ION4_RADIUS) */
double BOX2_ION4_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_4).number(BOX2_ION4_DIFF_COEFF) */
double BOX2_ION4_MASS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_4).number(BOX2_ION4_MASS) */
double BOX2_ION5_VALENCE;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_5).number(BOX2_ION5_VALENCE) */
double BOX2_ION5_RADIUS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_5).number(BOX2_ION5_RADIUS) */
double BOX2_ION5_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_5).number(BOX2_ION5_DIFF_COEFF) */
double BOX2_ION5_MASS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_5).number(BOX2_ION5_MASS) */
double BOX2_ION6_VALENCE;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_6).number(BOX2_ION6_VALENCE) */
double BOX2_ION6_RADIUS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_6).number(BOX2_ION6_RADIUS) */
double BOX2_ION6_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_6).number(BOX2_ION6_DIFF_COEFF) */
double BOX2_ION6_MASS;  /* for input.phase(page_7).group(ion_box_2).group(ion_box_2_ion_6).number(BOX2_ION6_MASS) */
double BOX3_ION1_VALENCE;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_1).number(BOX3_ION1_VALENCE) */
double BOX3_ION1_RADIUS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_1).number(BOX3_ION1_RADIUS) */
double BOX3_ION1_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_1).number(BOX3_ION1_DIFF_COEFF) */
double BOX3_ION1_MASS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_1).number(BOX3_ION1_MASS) */
double BOX3_ION2_VALENCE;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_2).number(BOX3_ION2_VALENCE) */
double BOX3_ION2_RADIUS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_2).number(BOX3_ION2_RADIUS) */
double BOX3_ION2_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_2).number(BOX3_ION2_DIFF_COEFF) */
double BOX3_ION2_MASS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_2).number(BOX3_ION2_MASS) */
double BOX3_ION3_VALENCE;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_3).number(BOX3_ION3_VALENCE) */
double BOX3_ION3_RADIUS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_3).number(BOX3_ION3_RADIUS) */
double BOX3_ION3_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_3).number(BOX3_ION3_DIFF_COEFF) */
double BOX3_ION3_MASS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_3).number(BOX3_ION3_MASS) */
double BOX3_ION4_VALENCE;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_4).number(BOX3_ION4_VALENCE) */
double BOX3_ION4_RADIUS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_4).number(BOX3_ION4_RADIUS) */
double BOX3_ION4_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_4).number(BOX3_ION4_DIFF_COEFF) */
double BOX3_ION4_MASS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_4).number(BOX3_ION4_MASS) */
double BOX3_ION5_VALENCE;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_5).number(BOX3_ION5_VALENCE) */
double BOX3_ION5_RADIUS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_5).number(BOX3_ION5_RADIUS) */
double BOX3_ION5_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_5).number(BOX3_ION5_DIFF_COEFF) */
double BOX3_ION5_MASS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_5).number(BOX3_ION5_MASS) */
double BOX3_ION6_VALENCE;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_6).number(BOX3_ION6_VALENCE) */
double BOX3_ION6_RADIUS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_6).number(BOX3_ION6_RADIUS) */
double BOX3_ION6_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_6).number(BOX3_ION6_DIFF_COEFF) */
double BOX3_ION6_MASS;  /* for input.phase(page_7).group(ion_box_3).group(ion_box_3_ion_6).number(BOX3_ION6_MASS) */
double BOX4_ION1_VALENCE;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_1).number(BOX4_ION1_VALENCE) */
double BOX4_ION1_RADIUS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_1).number(BOX4_ION1_RADIUS) */
double BOX4_ION1_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_1).number(BOX4_ION1_DIFF_COEFF) */
double BOX4_ION1_MASS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_1).number(BOX4_ION1_MASS) */
double BOX4_ION2_VALENCE;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_2).number(BOX4_ION2_VALENCE) */
double BOX4_ION2_RADIUS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_2).number(BOX4_ION2_RADIUS) */
double BOX4_ION2_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_2).number(BOX4_ION2_DIFF_COEFF) */
double BOX4_ION2_MASS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_2).number(BOX4_ION2_MASS) */
double BOX4_ION3_VALENCE;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_3).number(BOX4_ION3_VALENCE) */
double BOX4_ION3_RADIUS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_3).number(BOX4_ION3_RADIUS) */
double BOX4_ION3_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_3).number(BOX4_ION3_DIFF_COEFF) */
double BOX4_ION3_MASS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_3).number(BOX4_ION3_MASS) */
double BOX4_ION4_VALENCE;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_4).number(BOX4_ION4_VALENCE) */
double BOX4_ION4_RADIUS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_4).number(BOX4_ION4_RADIUS) */
double BOX4_ION4_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_4).number(BOX4_ION4_DIFF_COEFF) */
double BOX4_ION4_MASS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_4).number(BOX4_ION4_MASS) */
double BOX4_ION5_VALENCE;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_5).number(BOX4_ION5_VALENCE) */
double BOX4_ION5_RADIUS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_5).number(BOX4_ION5_RADIUS) */
double BOX4_ION5_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_5).number(BOX4_ION5_DIFF_COEFF) */
double BOX4_ION5_MASS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_5).number(BOX4_ION5_MASS) */
double BOX4_ION6_VALENCE;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_6).number(BOX4_ION6_VALENCE) */
double BOX4_ION6_RADIUS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_6).number(BOX4_ION6_RADIUS) */
double BOX4_ION6_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_6).number(BOX4_ION6_DIFF_COEFF) */
double BOX4_ION6_MASS;  /* for input.phase(page_7).group(ion_box_4).group(ion_box_4_ion_6).number(BOX4_ION6_MASS) */
double BOX5_ION1_VALENCE;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_1).number(BOX5_ION1_VALENCE) */
double BOX5_ION1_RADIUS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_1).number(BOX5_ION1_RADIUS) */
double BOX5_ION1_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_1).number(BOX5_ION1_DIFF_COEFF) */
double BOX5_ION1_MASS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_1).number(BOX5_ION1_MASS) */
double BOX5_ION2_VALENCE;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_2).number(BOX5_ION2_VALENCE) */
double BOX5_ION2_RADIUS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_2).number(BOX5_ION2_RADIUS) */
double BOX5_ION2_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_2).number(BOX5_ION2_DIFF_COEFF) */
double BOX5_ION2_MASS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_2).number(BOX5_ION2_MASS) */
double BOX5_ION3_VALENCE;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_3).number(BOX5_ION3_VALENCE) */
double BOX5_ION3_RADIUS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_3).number(BOX5_ION3_RADIUS) */
double BOX5_ION3_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_3).number(BOX5_ION3_DIFF_COEFF) */
double BOX5_ION3_MASS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_3).number(BOX5_ION3_MASS) */
double BOX5_ION4_VALENCE;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_4).number(BOX5_ION4_VALENCE) */
double BOX5_ION4_RADIUS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_4).number(BOX5_ION4_RADIUS) */
double BOX5_ION4_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_4).number(BOX5_ION4_DIFF_COEFF) */
double BOX5_ION4_MASS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_4).number(BOX5_ION4_MASS) */
double BOX5_ION5_VALENCE;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_5).number(BOX5_ION5_VALENCE) */
double BOX5_ION5_RADIUS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_5).number(BOX5_ION5_RADIUS) */
double BOX5_ION5_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_5).number(BOX5_ION5_DIFF_COEFF) */
double BOX5_ION5_MASS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_5).number(BOX5_ION5_MASS) */
double BOX5_ION6_VALENCE;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_6).number(BOX5_ION6_VALENCE) */
double BOX5_ION6_RADIUS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_6).number(BOX5_ION6_RADIUS) */
double BOX5_ION6_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_6).number(BOX5_ION6_DIFF_COEFF) */
double BOX5_ION6_MASS;  /* for input.phase(page_7).group(ion_box_5).group(ion_box_5_ion_6).number(BOX5_ION6_MASS) */
double BOX6_ION1_VALENCE;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_1).number(BOX6_ION1_VALENCE) */
double BOX6_ION1_RADIUS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_1).number(BOX6_ION1_RADIUS) */
double BOX6_ION1_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_1).number(BOX6_ION1_DIFF_COEFF) */
double BOX6_ION1_MASS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_1).number(BOX6_ION1_MASS) */
double BOX6_ION2_VALENCE;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_2).number(BOX6_ION2_VALENCE) */
double BOX6_ION2_RADIUS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_2).number(BOX6_ION2_RADIUS) */
double BOX6_ION2_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_2).number(BOX6_ION2_DIFF_COEFF) */
double BOX6_ION2_MASS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_2).number(BOX6_ION2_MASS) */
double BOX6_ION3_VALENCE;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_3).number(BOX6_ION3_VALENCE) */
double BOX6_ION3_RADIUS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_3).number(BOX6_ION3_RADIUS) */
double BOX6_ION3_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_3).number(BOX6_ION3_DIFF_COEFF) */
double BOX6_ION3_MASS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_3).number(BOX6_ION3_MASS) */
double BOX6_ION4_VALENCE;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_4).number(BOX6_ION4_VALENCE) */
double BOX6_ION4_RADIUS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_4).number(BOX6_ION4_RADIUS) */
double BOX6_ION4_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_4).number(BOX6_ION4_DIFF_COEFF) */
double BOX6_ION4_MASS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_4).number(BOX6_ION4_MASS) */
double BOX6_ION5_VALENCE;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_5).number(BOX6_ION5_VALENCE) */
double BOX6_ION5_RADIUS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_5).number(BOX6_ION5_RADIUS) */
double BOX6_ION5_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_5).number(BOX6_ION5_DIFF_COEFF) */
double BOX6_ION5_MASS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_5).number(BOX6_ION5_MASS) */
double BOX6_ION6_VALENCE;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_6).number(BOX6_ION6_VALENCE) */
double BOX6_ION6_RADIUS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_6).number(BOX6_ION6_RADIUS) */
double BOX6_ION6_DIFF_COEFF;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_6).number(BOX6_ION6_DIFF_COEFF) */
double BOX6_ION6_MASS;  /* for input.phase(page_7).group(ion_box_6).group(ion_box_6_ion_6).number(BOX6_ION6_MASS) */









#endif