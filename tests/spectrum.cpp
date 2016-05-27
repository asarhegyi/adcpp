/************************************************************
* real_FFT -- unit test the functions of the spectrum class.*
*		The program reads in a data set from a file	and does*
*		an FFT on it.										*
*															*
* Author: Attila Sarhegyi									*
*															*
* Version: 1.0												*
*															*
* Purpose: Do an initial estimate of the sine-wave			*
*		parameters before the least-suqares fitting.		*
*															*
* Usage:													*
*		make -f ../Makefile TARGET=spectrum					*
*		./spectrum											*
************************************************************/
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_math.h>

#include <iostream>		// std::cin, std::cout
#include <fstream>		// std::ifstream, std::ofstream
#include <vector>		// std::vector
#include "../src/utils.h"	// dbl_array2file, dbl_file2array
#include "../src/spectrum_class.h"



int main (void)
{
	/************************************************************
	* Open the input file (tes.int) and load the content into	*
	* the data vector. Return the number of elements in the		*
	* vector.													*
	************************************************************/
	std::vector<double> data;
	int N = dbl_file2array ("../charData/test.int", data);


	/************************************************************
	* Set-up the gsl-fft workspace, calculate the FFT and find	*
	* the maximum of the magnitude spectrum.					*
	************************************************************/
	std::cout << "Calculating FFT..." << std::endl << std::endl;

	class spectrum myfft(data);
	dbl_array2file("magnitude.dat", myfft.magnitude, N/2, "");

	std::cout << "Post-processing FFT results..." << std::endl << std::endl;

	double *maxMag;
	maxMag = myfft.findMaxMagnitude();

	double fSample;
	fSample	= 5.0*pow(10.0,5.0);

	std::cout << "	maxIndex           : " << maxMag[1] << std::endl;
	std::cout << "	magnitude[maxIndex]: " << maxMag[0] << std::endl;
	std::cout << "	maxFrequency       : " << maxMag[1]*(fSample/N) << std::endl;
	std::cout << "	magnitude[0]       : " << myfft.magnitude[0] << std::endl;


	double lambda;
	lambda = myfft.IpFFT();

	//std::cout << std::setprecision(11) << std::scientific << std::endl;
	std::cout << "	Ip. Frequency      : " << lambda*(fSample/N) << std::endl;

	return(0);
}





