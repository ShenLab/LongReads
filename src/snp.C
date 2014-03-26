#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>

//#define DEBUG

using namespace std;

#include "snp.h"

//SNP::SNP(string chr_s, long snp_pos, char snp_ref, char snp_alt, vector<string> gl3, int known_par, double qualscore)
SNP::SNP(string chr_s, long snp_pos, char snp_ref, char snp_alt, int known_par, double qualscore)
{
	strcpy(chr,chr_s.c_str());
	known = known_par;
	ref = snp_ref;
	alt = snp_alt;
	refcount = altcount = errcount = 0;
	count = -1;
	overlap_count = 0;
	known_overlap_count = 0;
	position = snp_pos;
	qual = qualscore;
	likelihood_ratio = -1;
	reads = new READ*[200];
	//gl = new double[3];
	//posterior = new double[3];
	somatic_posterior = new double[3];

	//gl[0] = pow(10.0, -(atof(gl3[0].c_str()))/10.0);
	//gl[1] = pow(10.0, -(atof(gl3[1].c_str()))/10.0);
	//gl[2] = pow(10.0, -(atof(gl3[2].c_str()))/10.0);
	//double norm = gl[0]+gl[1]+gl[2];

	//posterior[0] = gl[0]/norm;
	//posterior[1] = gl[1]/norm;
	//posterior[2] = gl[2]/norm;

	somatic_posterior[0] = 1.0;
	somatic_posterior[1] = 1.0;
	somatic_posterior[2] = 1.0;

	emission[1] = 0.0;
	emission[2] = 0.0;
}

SNP::~SNP()
{
	delete [] reads;
	//delete [] gl;
	//delete [] posterior;
	delete [] somatic_posterior;
}

char* SNP::GetChr()
{
	return chr;
}

long SNP::GetPos()
{
	return position;
}

void SNP::append(int type, READ *read)
{
	add_read(type, read);
}

void SNP::add_read(int type, READ *read)
{
	if(type==0)
		refcount++;
	else if(type==1)
		altcount++;
	else
		errcount++;
	reads[++count] = read;
}

void SNP::add_somatic_posteriors(double post[3])
{
	for(int i=0; i<3; i++) {
		somatic_posterior[i] = post[i];
	}
}

void SNP::addEmission(int haplotype, double probability)
{
	emission[haplotype] = probability;
}

void SNP::IncrOverlapCount()
{
	overlap_count++; 
}

void SNP::IncrKnownOverlapCount()
{
	known_overlap_count++; 
}

READ* SNP::GetRead(int pos)
{
	if(pos==-1)
		return NULL;
	else
		return reads[pos];
}

char SNP::GetRef()
{
	return ref;
}

char SNP::GetAlt()
{
	return alt;
}

int SNP::GetReadCount()
{
	return count + 1;
}

int SNP::GetRefCount()
{
	return refcount;
}

int SNP::GetAltCount()
{
	return altcount;
}

int SNP::GetErrCount()
{
	return errcount;
}

int SNP::GetKnown()
{
	return known;
}

double SNP::GetQualScore()
{
	return qual;
}

double* SNP::GetSomaticPosteriors()
{
	return somatic_posterior;
}

int SNP::GetEmissionState()
{
	if(emission[1]>emission[2])
		return 1;
	else if(emission[2]>emission[1])
		return 2;
	else {
		cout << "Error for snp " << position << " with equal emissions requested bearing value " << emission[1] << " " << emission[2] << endl;
	}
}

int SNP::GetOverlapCount()
{
	return overlap_count;
}

int SNP::GetKnownOverlapCount()
{
	return known_overlap_count;
}

/*
void SNP::assign_genotype(int gt, double genp)
{
	genotype = gt;
	genprob = genp;
}
*/

/*
double* SNP::GetPosteriors()
{
	return posterior;
}
*/
/*
double* SNP::GetGenLik()
{
	return gl;
}
*/
/*
void SNP::add_posteriors(double post[3])
{
	for(int i=0; i<3; i++) {
		posterior[i] = post[i];
	}
}
*/
