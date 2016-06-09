/************************************************************
* debug_resolution -- test the sensitivity of the sine fit	*
*		algorithm with respect to resolution for Monte Carlo*
*		simulations.										*
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
#include <gsl/gsl_rng.h>	 // random number generator
#include <gsl/gsl_randist.h>
#include "../src/stimulus.h"

double vmaxIn;
char *program_name;         // name of the program (for errors)



/************************************************************
* usage -- tell the user how to use this program and exit	*
************************************************************/
void usage()
{
	std::cerr << "Usage is " << program_name << " [options]\n";
	std::cerr << "Options\n";
	std::cerr << "  -v<number>  set Vmax\n";
	std::cerr << "  -h[elp]    or\n";
	std::cerr << "  -u[sage]    Print this help message\n";
	exit(8);
}


/************************************************************
*****    *****    *********   *********  ****    *****  *****
*****  *  ***  *  ********  *  ********  ****  *  ****  *****
*****  **  *  **  *******  ***  *******  ****  **  ***  *****
*****  ***   ***  ******         ******  ****  ***  **  *****
*****  *********  *****  *******  *****  ****  ****  *  *****
*****  *********  ****  *********  ****  ****  *****    *****
************************************************************/

int main (int argc, char *argv[])
{
	// Save the program name for future use
	program_name = argv[0];

	// Loop for each option. Stop if we run out of arguments or
	// we get an argument without a dash.
	while ((argc > 1) && (argv[1][0] == '-'))
   	{
		// argv[1][1] is the actual option character
		switch (argv[1][1])
	   	{
			// -v<number> set vmax
			case 'v':
				vmaxIn = strtod(&argv[1][2], NULL);
				break;
			// -h[elp] or -u[sage] prints help
			case 'h':
			case 'u':
				usage();
				break;
			default:
				std::cerr << "Bad option " << argv[1] << "\n\n";
				usage();
		}

		// move the argument list up one, and move the count down one
		++argv;
		--argc;
	}


	// parameters of the quantizer
	double vmax = vmaxIn;                         // maximum of the input voltage range
	double vmin = -vmax;                          // minimum of the input voltage range
	int qbits = 15;                               // number of quantizer bits

	// parameters of the sine-wave
	double amplitude = 1.27651245586844375168;    // Amplitude of the sine wave
	double frequency = 1500.08668950892047178058; // Frequency of the sine wave
	double phase = 2.08259279452065149130;        // Phase of the sine wave
	double dc = -0.45650260888505728163;          // dc of the sine wave

	double fs = 43899.77302618323301430792;       // sampling frequency
	int nofSamples = 32275;                       // number of samples to generate

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
