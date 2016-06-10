/************************************************************
* gen_stimulus -- generate a noise free sine-wave, and		*
*		ideal code transition levels. Add Gaussian noise to	*
*		both of them, and quantize the resulting signal		*
*		with the non-ideal quantizer.						*
*															*
* Author: Attila Sarhegyi									*
*															*
* Date : 11:29 Dec. 22 2015									*
*															*
* Version: 1.0												*
*															*
* Purpose: Experiment various Least-Squares algorithms and	*
*		compare it with MATLAB								*
*															*
* Usage:													*
*		test_stimulus										*
************************************************************/
#include <stdio.h>
#include <iostream>
#include "../src/stimulus.h"    // generate input signal for simulation

double nofBits;              // number of quantizer bits
char *program_name;          // name of the program (for errors)


/************************************************************
* usage -- tell the user how to use this program and exit	*
************************************************************/
void usage()
{
	std::cerr << "Usage is " << program_name << " [options]\n";
	std::cerr << "Options\n";
	std::cerr << "  -n<number>  Number of quantizer bits\n";
	std::cerr << "  -h[elp]    or\n";
	std::cerr << "  -u[sage]    Print this help message\n";
	exit(8);
}


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
			// -n<number> set number of quantizer bits
			case 'n':
				nofBits = atoi(&argv[1][2]);
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


	// create a stimulus class object
	class stimulus mystimulus;


	/************************************************************
	 * Generate noise free sine-wave or ramp signals
	 ************************************************************/
	//define parameters
	int nofSamples = 512;        // number of samples to generate
	double fs = 60000.56;        // sampling frequency

	double amplitude = 2.1;	     // Amplitude of the sine wave
	//double frequency = 1000.23;  // Frequency of the sine wave
	double frequency = 151.23;  // Frequency of the sine wave
	double phase = M_PIl;        // Phase of the sine wave
	double dc = 0.0;             // dc of the sine wave

	// generate noise free sine wave
	mystimulus.cosine(amplitude, frequency, phase, dc, nofSamples, fs);

	// generate noise free ramp (overwrites the sine-wave!)
	//mystimulus.ramp(amplitude, dc, nofSamples, fs);

	// save signal to file
	FILE *FOUT;
	FOUT = fopen("data.dat", "w");
	for(int index = 0; index < mystimulus.nofSamples; index++)
	{
		fprintf(FOUT, "%.14f\n", mystimulus.data[index]);
	}
	fclose(FOUT);


	/************************************************************
	 * Add noise to the generated noise-free signal
	 ************************************************************/
	double mean = 0.0;
	double sigma = 0.03;

	// add Gauss noise to the signal stored in &mystimulus.data
	mystimulus.add_noise_signal(mean, sigma, &mystimulus.data[0]);

	// save noisy signal to file
	FILE *FNOUT;
	FNOUT = fopen("noisyData.dat", "w");
	for(int index = 0; index < mystimulus.nofSamples; index++)
	{
		fprintf(FNOUT, "%.14f\n", mystimulus.noisyData[index]);
	}
	fclose(FNOUT);


	/************************************************************
	 * Generate ideal code transition levels for the quantizer
	 ************************************************************/
	double vmax = 2.15;
	double vmin = -2.15;
	double qbits = nofBits;

	// define code transition levels
	mystimulus.get_tl_ideal_quantizer(vmax, vmin, qbits);

	// save code transition levels to file
	FILE *FTOUT;
	FTOUT = fopen("Tl.dat", "w");
	for(int index = 0; index < mystimulus.nofCodes-1; index++)
	{
		fprintf(FTOUT, "%.14f\n", mystimulus.Tl[index]);
	}
	fclose(FTOUT);


	/************************************************************
	 * quantize the noise contaminated signal with an ideal quantizer 
	 ************************************************************/
	mystimulus.quantize(&mystimulus.Tl[0], &mystimulus.noisyData[0]);

	// save quantized signal to file
	FILE *FQOUT;
	FQOUT = fopen("quantizedData.dat", "w");
	for(int index = 0; index < mystimulus.nofSamples; index++)
	{
		fprintf(FQOUT, "%d\n", mystimulus.quantizedData[index]);
	}
	fclose(FQOUT);


	/************************************************************
	 * Add noise to CTLs
	 ************************************************************/
	// add Gauss noise to CTL stored in &mystimulus.Tl
	mystimulus.add_noise_Tl(&mystimulus.Tl[0]);

	// save noisy CTL to file
	FILE *FNTOUT;
	FNTOUT = fopen("noisyTl.dat", "w");
	for(int index = 0; index < mystimulus.nofCodes-1; index++)
	{
		fprintf(FNTOUT, "%.14f\n", mystimulus.noisyTl[index]);
	}
	fclose(FNTOUT);


	/************************************************************
	 * quantize noise-free signal with an ideal quantizer
	 ************************************************************/
	mystimulus.quantize(&mystimulus.Tl[0], &mystimulus.data[0]);

	// save quantized signal to file
	FILE *FQIOUT;
	FQIOUT = fopen("quantizedData_idealTl.dat", "w");
	for(int index = 0; index < mystimulus.nofSamples; index++)
	{
		fprintf(FQIOUT, "%d\n", mystimulus.quantizedData[index]);
	}
	fclose(FQIOUT);


	/************************************************************
	 * quantize noise-free with a non-ideal quantizer
	 * NOTE that this will overwrite the content of
	 * mystimulus.quantizedData
	 ************************************************************/
	mystimulus.quantize(&mystimulus.noisyTl[0], &mystimulus.data[0]);

	// save quantized signal to file
	FILE *FQNOUT;
	FQNOUT = fopen("quantizedData_noisyTl.dat", "w");
	for(int index = 0; index < mystimulus.nofSamples; index++)
	{
		fprintf(FQNOUT, "%d\n", mystimulus.quantizedData[index]);
	}
	fclose(FQNOUT);


	return 0;
}
