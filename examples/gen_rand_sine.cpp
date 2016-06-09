/************************************************************
* gen_rand_sine -- generate random parameter sine-wave		*
*		for Monte Carlo simulations.						*
*															*
* Author: Attila Sarhegyi									*
*															*
* Version: 1.0												*
*															*
* Purpose: Experiment various Least-Squares algorithms and	*
*		compare it with MATLAB								*
*															*
* Usage:													*
*		sfit3_linear -help									*
************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_rng.h>	  // random number generator
#include <gsl/gsl_randist.h>
#include "../src/stimulus.h"


/************************************************************
*****    *****    *********   *********  ****    *****  *****
*****  *  ***  *  ********  *  ********  ****  *  ****  *****
*****  **  *  **  *******  ***  *******  ****  **  ***  *****
*****  ***   ***  ******         ******  ****  ***  **  *****
*****  *********  *****  *******  *****  ****  ****  *  *****
*****  *********  ****  *********  ****  ****  *****    *****
************************************************************/

int main ()
{
	int nofSamples = 32768;      // number of samples to generate
	double fs = 60000.52;        // sampling frequency

	// parameters of the quantizer
	double vmax = 2.1;           // maximum of the input voltage range
	double vmin = -2.1;          // minimum of the input voltage range
	int qbits = 12;              // number of quantizer bits

	// parameters of the sine-wave
	double amplitude = vmax/1.5; // Amplitude of the sine wave
	double frequency = 1000.23;  // Frequency of the sine wave
	double phase = M_PIl;        // Phase of the sine wave
	double dc = 0.0;             // dc of the sine wave

	// additive Gaussian noise parameters for the signal
	double mean_s = 0.0;
	double sigma_s = 0.03;
	

	// initialize the random number generator
	gsl_rng * r;
	gsl_rng_env_setup();
	r = gsl_rng_alloc (gsl_rng_default);

	// read in the random number state from file
	// into the random number generator if file exits
	if( access ("gsl_rng_uniform_state", F_OK) != -1 )
   	{
		FILE *FRIN;
		FRIN = fopen("gsl_rng_uniform_state", "r");
		gsl_rng_fread (FRIN, r);
		fclose(FRIN);
	}
	
	// perturb parameters
	double u = gsl_rng_uniform(r);
	frequency = (u+0.5)*frequency;

	u = gsl_rng_uniform(r);
	phase = (u+0.5)*phase;

	u = gsl_rng_uniform(r);
	dc = vmin+(0.25+0.5*u)*(vmax-vmin);

	do
	{
		u = gsl_rng_uniform(r);
		amplitude = fmin(vmax-dc,dc-vmin)*(1-6*sigma_s)*u;
	} 
	while (amplitude < 0.15);

	u = gsl_rng_uniform(r);
	fs = (u+0.5)*fs;

	u = gsl_rng_uniform(r);
	nofSamples = (u+0.5)*nofSamples;

	u = gsl_rng_uniform(r);
	qbits = qbits+(u-0.5)*8;  //qbits+/-4

	// save the random number state of the random number generator
	FILE *FROUT;
	FROUT = fopen("gsl_rng_uniform_state", "w");
	gsl_rng_fwrite (FROUT, r);
	fclose(FROUT);

	gsl_rng_free(r);


	// create a stimulus class object
	class stimulus mystimulus;

	// generate noise free sine-wave
	mystimulus.cosine(amplitude, frequency, phase, dc, nofSamples, fs);

	// generate ideal code transition levels
	mystimulus.get_tl_ideal_quantizer(vmax, vmin, qbits);

	// add Gauss noise to the signal stored in &mystimulus.data
	mystimulus.add_noise_signal(mean_s, sigma_s, &mystimulus.data[0]);

	// add Gauss noise to CTL stored in &mystimulus.Tl
	mystimulus.add_noise_Tl(&mystimulus.Tl[0]);

	// quantize noise contaminated signal with a non-ideal quantizer
	//mystimulus.quantize(&mystimulus.noisyTl[0], &mystimulus.noisyData[0]);
	mystimulus.quantize(&mystimulus.Tl[0], &mystimulus.noisyData[0]);


	// save quantized data to file
	FILE *FQOUT;
	FQOUT = fopen("quantizedData.dat", "w");
	for(int index = 0; index < mystimulus.nofSamples; index++)
	{
		fprintf(FQOUT, "%d\n", mystimulus.quantizedData[index]);
	}
	fclose(FQOUT);


	// save parameters to file
	FILE *FPOUT;
	FPOUT = fopen("parameters.dat", "w");
	fprintf(FPOUT, "%.20f %.20f %.20f %.20f %.20f %d %d %.20f %.6f %.6f\n", amplitude, frequency, phase, dc, fs, nofSamples, qbits, mystimulus.Q, vmax, vmin);
	fclose(FPOUT);

	return 0;
}
