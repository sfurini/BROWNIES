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


bool fileExists(string fileName){
	bool result=false;
	char *fN = new char[fileName.length()+1];
	strcpy(fN, fileName.c_str()); 
	ifstream myFile(fN);
 	if(myFile.is_open()){
		result=true;	
		myFile.close();
	}
	return result;
}

bool directoryExists(string dirName){
	bool result=false;
	
	char *dName = new char[dirName.length()+1];
    strcpy(dName, dirName.c_str());
	
	DIR *pdir=opendir(dName);
	
	if(pdir!=NULL){
		result=true;
	}
	
	return result;
}

string deleteDirectory(string dirName){
	string result="";
	int i=0;
	
	char *dName = new char[dirName.length()+1];
   strcpy(dName, dirName.c_str());
	
	string cmd="rm -rf ";
	cmd=cmd+dirName;
	char *command = new char[cmd.length()+1];
	strcpy(command, cmd.c_str()); 	
	
	ifstream inputStream;
 	inputStream.open(dName);
	
 	if(!inputStream){
		result="ok";
 	}    
	else{
		//cout << "Deleting directory " << dirName << endl;	
		i=system(command);
		
		if(i==0){
			result="ok";	
		}
		else{
			result="NO!!!";	
		}	
	}
	inputStream.close();
	
	return result;	
}

string createDirectory(string dirName){
	string result="", temp="";	
	int i=0;
	temp=deleteDirectory(dirName);	
	if(ok(temp)){
		string cmd="mkdir -p ";
		cmd=cmd+dirName;
		char *command = new char[cmd.length()+1];
		strcpy(command, cmd.c_str()); 
		i=system(command);
	}
	else{
		i=1000;	
	}	
	if(i==0){
		result="ok";	
	}
	else{
		result="NO!!!";	
	}
	return result;
}

string deleteFile(string fileName){
	string result="";
	int i=0;
	
	char *fName = new char[fileName.length()+1];
	strcpy(fName, fileName.c_str());
	
	string cmd="rm ";
	cmd=cmd+fileName;
	char *command = new char[cmd.length()+1];
   strcpy(command, cmd.c_str()); 	
	
	ifstream inputStream;
 	inputStream.open(fName);
	
 	if(!inputStream){
		result="ok";
 	}    
	else{
		//cout << "Deleting file " << fileName << endl;	
		i=system(command);
		
		if(i==0){
			result="ok";	
		}
		else{
			result="NO!!!";	
		}	
	}
	inputStream.close();
	
	return result;	
}


string createFile(string fileName){
	string result="", temp="";	
	int i=0;
	temp=deleteFile(fileName);	
	if(ok(temp)){
		string cmd="touch ";
		cmd=cmd+fileName;
		char *command = new char[cmd.length()+1];
		strcpy(command, cmd.c_str()); 
		i=system(command);
	} else {
		i=1000;	
	}	
	if(i==0){
		result="ok";	
	} else {
		result="NO!!!";	
	}
	return result;
}


int countLinesInFile(string fileIn){
	int totalLines=0;
	string buffer="";	
		
	char *tfn = new char[fileIn.length()+1];
	strcpy(tfn, fileIn.c_str()); 
    
	ifstream fin;        
	fin.open(tfn, ifstream::in);   
	while(!fin.eof()){
		if(getline(fin, buffer)){
			totalLines++;  
		}        	
	}
	fin.close();
    
	return totalLines;		
}

bool ok(string test){
	bool result=false;
	if(!test.compare("ok")){
		result=true;	
	}
	return result;
}

bool verifyExtension(string fileName, string extension){
	bool result=false;
	string name="", ext="";
	int n=fileName.find_last_of("/");
	if(n==string::npos){
		name=fileName;
	} else { 
		name=fileName.substr(n+1);
	}
	int i=name.find_last_of(".");
	if(i==string::npos){
		ext="NONE";
	} else {
		ext=name.substr(i+1);
	}
	if(!ext.compare(extension)){
		result=true;	
	} 
	return result;	
}

int countTokens(string myString, string separator){
	int count=0;
	int firstSep=0, firstNotSep=0;

	firstSep=myString.find_first_of(separator);
	firstNotSep=myString.find_first_not_of(separator);

	if(firstSep==string::npos && firstNotSep==string::npos){
		//stringa vuota
		count=0;
	}
		
	else if(firstSep==string::npos){
		//stringa senza separatori
		count=1;
	}
		
	else if(firstNotSep==string::npos){
		//stringa solo separatori
		count=0;
	}
	else{
		count=1+countSeparators(myString, separator);
	}

	return count;	
}

string getTokenbyNumber(string myString, string separator, int tokenNumber){
	
	int numberOfTokens=1+countSeparators(myString, separator);
	stringstream oss;
    oss << "La stringa " << myString << " contiene solo " << numberOfTokens + " token.";
	string result=(oss.str());
	string buffer=myString;
	int i=0, begin=0, end=0, offset=0;
	
	int stringLength=myString.length();	
	
	if(tokenNumber>numberOfTokens){		
		
	}
	else if(tokenNumber<=0){
		result="Occorre inserire un numero intero positivo per identificare i token.";
	}
	
	else{
		begin=buffer.find_first_not_of(separator);
		end=buffer.find_first_of(separator);
		if(end==string::npos && begin!=string::npos){
			result=1;
		}
		else{
			do{
				begin=buffer.find_first_not_of(separator);
				buffer=myString.substr(offset+begin);
				end=buffer.find_first_of(separator);
				offset=offset+begin+end;
				result=buffer.substr(0, end);
				buffer=myString.substr(offset);  
				i++;			             
			}
			while(i<tokenNumber);	
		}
	}	
	return result;
}

int countSeparators(string myString, string separator){
	int firstSep=0, firstNotSep=0;
	bool exit=false;
	int count=0;
	int stringLength=myString.length();
	
	string buffer=myString;
	string temp="";

	do{	
		firstSep=buffer.find_first_of(separator);
		firstNotSep=buffer.find_first_not_of(separator);

		if(firstSep==string::npos && firstNotSep==string::npos){
			exit=true;
		}
		
		else if(firstSep==string::npos){
			exit=true;
		}
		
		else if(firstNotSep==string::npos){
			exit=true;
		}
		else if(firstSep==0){
			temp=buffer.substr(1, buffer.length()-1);
			buffer=temp;			
		}
		else if(firstSep==buffer.length()-1){
			temp=buffer.substr(0, buffer.length()-1);
			buffer=temp;		
		} 
		else{	
			temp=buffer.substr(firstSep+1);
			firstNotSep=temp.find_first_not_of(separator);			
			if(firstNotSep==string::npos){
				exit=true;		
			}
			else{
				buffer=temp;
				count++;
			}
		}
		
	}
	while(!exit);
	
	return count;  

}

void removeAllWhite(string &str){
	string temp;
	for(unsigned int i=0; i<str.length(); i++){
		if(str[i] != ' ' && str[i] != '\t'){
			temp+=str[i];
		}
	}
	str=temp;

	return;
}





