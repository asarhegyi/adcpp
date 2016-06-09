#ifndef GENERATOR_CLASS
#define GENERATOR_CLASS

#include <math.h>
#include <unistd.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include "utils.h"		// save array into file

/************************************************************
* stimulus class											*
*		generate a sine-wave stimulus with noise and		*
*		distortion for the sine-fit simulations.			*
*															*
* Parameters:												*
*		See parameters below.								*
*															*
* Member functions											*
*		cosine -- generate x=A*cos(2*pi*f*index/fs*phi)+dc	*
*		ramp   -- generate x=-A+index*(2*A)/(N-1)+dc		*
*		add_noise_signal -- add Gaussian noise to the signal*
*		get_tl_ideal_quantizer -- generate the Code			*
*			 Transition Levels (CTL) for an ideal quantizer	*
*		add_noise_Tl -- add uniformly distributed noise to	*
*			 CTLs											*
*		quantize -- quantize the input signal using the CTLs*
************************************************************/
class stimulus {
	public:
		int verbose;         // verbosity level
		// parameters of the sine wave and ramp
		double amplitude;    // Amplitude of the sine wave
		double frequency;    // Frequency of the sine wave
		double omega;        // Angular frequency of the sine wave
		double phase;        // Phase of the sine wave
		double dc;           // dc of the sine wave
		// parameters of the data sampling
		int nofSamples;      // number of samples to generate
		double fs;           // sampling frequency
		double * data;       // samples of the generated sine wave (pointer)
		// parameters of the addition Gaussian noise added to the signal
		double mean_s;       // mean of the noise added to the signal
		double sigma_s;      // standard deviation if the noise added to the signal
		double * noisyData;  // the generated signal + noise (pointer)
		// parameters of the quantizer
		double Vmax;         // maximum input voltage of the quantizer
		double Vmin;         // minimum input voltage of the quantizer
		int qBits;           // number of bits of the quantizer
		int nofCodes;        // number of output digital code 2^qBits
		double Q;            // step size of the ideal quantizer 
		double * Tl;         // code transition levels (pointer)
		double * noisyTl;     // the generated signal + noise (pointer)

		int * quantizedData;  // quantized data

		// Constructor
		stimulus()
		{
			verbose = 0;
		}

		//stimulus(const stimulus& old_stimulus);
			// Use automatically generated copy constructor

		//stimulus operator = (const stimulus& old_stimulus);
			// Use automatically generated assign constructor

		//~stimulus();
			// Use automatically generated destructor

		// generate cosine x_index = A*cos(2*pi*f*index/fs+phi)+dc
		void cosine(double ampIn, double freqIn, double  phiIn, double dcIn, int N, double fsIn);

		// generate ramp x_index = -A+index*(2*A)/(N-1)+d
		void ramp(double ampIn, double dcIn, int N, double fsIn);

		// add Gaussian noise to the signal y_i = x_i + e_i
		void add_noise_signal(double meanIn, double sigmaIn, double dataIn[]);

		// generate the code transition levels of an ideal quantizer 
		void get_tl_ideal_quantizer(double vmaxIn, double vminIn, int qbitsIn);

		// add uniformly distributed (+/- Q/2) noise to CTLs
		void add_noise_Tl(double TlIn[]);

		// quantize data
		void quantize(double TlIn[], double dataIn[]);

};


/************************************************************
* stimulus::cosine -- generate a cosine waveform and		*
*		store it in the class array variable called data.	*
*	x_index = ampIn*cos(2*pi*freqIn*index/fsIn+phiIn)+dcIn	*
*															*
* Parameters												*
*		ampIn  -- amplitude of the ramp						*
*		freqIn -- frequency of the sine wave				*
*		phiIn  -- phase of the sine wave					*
*		dcIn   -- dc offset of the ramp						*
*		N      -- number of samples to generate				*
*		fsIn   -- sampling frequency						*
*															*
* Return													*
*		waveform is stored in the class array variable data	*
************************************************************/
inline void stimulus::cosine(double ampIn, double freqIn, double phiIn, double dcIn, int N, double fsIn)
{
	nofSamples = N;
	fs = fsIn;
	amplitude = ampIn;
	frequency = freqIn;
	phase = phiIn;
	dc = dcIn;

	omega = 2*M_PIl*frequency;
	data = new double[nofSamples];

	for (int index = 0; index < nofSamples; index++)
		data[index] = amplitude*cos(omega*index/fs+phase)+dc;
}


/************************************************************
* stimulus::ramp -- generate a ramp waveform and store		*
*		it in the class array variable called data.			*
*		x_index = -ampIn+index*(2*ampIn)/(N-1)+dcIn			*
*															*
* Parameters												*
*		ampIn -- amplitude of the sine wave					*
*		dcIn  -- dc of the sine wave						*
*		N     -- number of samples to generate				*
*		fsIn  -- sampling frequency							*
*															*
* Return													*
*		waveform is stored in the class array variable data	*
************************************************************/
inline void stimulus::ramp(double ampIn, double dcIn, int N, double fsIn)
{
	nofSamples = N;
	fs = fsIn;
	amplitude = ampIn;
	dc = dcIn;

	data = new double[nofSamples];

	for (int index = 0; index < nofSamples; index++)
		data[index] = -amplitude+index*(2*amplitude)/(nofSamples-1)+dc;
}


/************************************************************
* stimulus::add_noise_signal -- add Gaussian noise to		*
*		the signal that gets passed in.						*
*															*
* Parameters												*
*		meanIn  -- mean value of the additive Gauss noise	*
*		sigmaIn -- standard deviation of the additive noise	* 
*		dataIn  -- input signal								*
*															*
* Return													*
*		the contaminated waveform is stored in noisyData	*
************************************************************/
inline void stimulus::add_noise_signal(double meanIn, double sigmaIn, double dataIn[])
{
	// initialize the random number generator
	gsl_rng * rg;
	gsl_rng_env_setup();
	rg = gsl_rng_alloc (gsl_rng_default);

	// read in the random number state from file
	// into the random number generator if file exits
	if( access ("gsl_rng_gaussian_state", F_OK) != -1 )
   	{
		FILE *FRGIN;
		FRGIN = fopen("gsl_rng_gaussian_state", "r");
		gsl_rng_fread (FRGIN, rg);
		fclose(FRGIN);
	}

	mean_s = meanIn;
	sigma_s = sigmaIn;
	// declare the size of the array
	noisyData = new double[nofSamples];

	for (int index = 0; index < nofSamples; index++)
		noisyData[index] = dataIn[index] + mean_s + gsl_ran_gaussian (rg, sigma_s);

	// save the random number state of the random number generator
	FILE *FRGOUT;
	FRGOUT = fopen("gsl_rng_gaussian_state", "w");
	gsl_rng_fwrite (FRGOUT, rg);
	fclose(FRGOUT);

	gsl_rng_free(rg);
}


/************************************************************
* stimulus::get_tl_ideal_quantizer -- generate the code		*
*		transition levels of an ideal quantizer.			*
*															*
* Parameters												*
*		vmaxIn -- maximum input voltage of the quantizer	*
*		vminIn -- minimum input voltage of the quantizer	*
*		qBitsIn -- number of bits of the quantizer			*
*															*
* Return													*
*		the code transition levels are stored in the array Tl*
************************************************************/
inline void stimulus::get_tl_ideal_quantizer(double vmaxIn, double vminIn, int qBitsIn)
{
	Vmax = vmaxIn;
	Vmin = vminIn;
	qBits = qBitsIn;
	nofCodes = pow(2.0, qBits);

	Q = (Vmax-Vmin)/nofCodes;          // step sise of the ideal quantizer
	Tl = new double[nofCodes-1];       // code transition levels Tl[0] <= T1

	for (int index = 0; index < (nofCodes-1); index++)
		Tl[index] = Vmin+Q/2+index*Q;
}


/************************************************************
* stimulus::add_noise_Tl -- add uniformly distributed		*
*		noise to the Code Transition Levels (CTL) that gets	*
*		passed in. Tl falls in the range of Tl+/-Q/2 when	*
*		Q is the step size of the ideal quantizer.			*
*															*
* Parameters												*
*		TlIn    -- input Code Transition Levels				*
*															*
* Return													*
*		the contaminated CTLs are stored in noisyTl			*
************************************************************/
inline void stimulus::add_noise_Tl(double TlIn[])
{
	// initialize the random number generator
	gsl_rng * ru;
	gsl_rng_env_setup();
	ru = gsl_rng_alloc (gsl_rng_default);

	// read in the random number state from file
	// into the random number generator if file exits
	if( access ("gsl_rng_uniform_state", F_OK) != -1 )
   	{
		FILE *FRUIN;
		FRUIN = fopen("gsl_rng_uniform_state", "r");
		gsl_rng_fread (FRUIN, ru);
		fclose(FRUIN);
	}

	// declare the size of the array
	noisyTl = new double[nofCodes-1];       // code transition levels Tl[0] <= T1

	for (int index = 0; index < (nofCodes-1); index++)
	{
		noisyTl[index] = TlIn[index] + (gsl_rng_uniform_pos(ru)-0.5)*Q;

		// check if noisyTl is monotonic
		if ((index > 0) && (noisyTl[index-1] > noisyTl[index]))
		{
			fprintf(stderr, "### ERROR: noisyTl is non-monotonic at %d - %d\nExiting...\n", index-1, index);
			dbl_array2file("noisyTl.dat", &noisyTl[0], index+1, "");
			exit(1);
		}
	}

	// save the random number state of the random number generator
	FILE *FRUOUT;
	FRUOUT = fopen("gsl_rng_uniform_state", "w");
	gsl_rng_fwrite (FRUOUT, ru);
	fclose(FRUOUT);

	gsl_rng_free(ru);
}


/************************************************************
* stimulus::quantizer -- quantize the input signal			*
*		using the code transition levels passed in.			*
*															*
* Parameters												*
*		dataIn  -- input signal								*
*		TlIn    -- input Code Transition Levels				*
*															*
* Return													*
*		the result of the quantization is stored in			*
*		quantizedData.										*
************************************************************/
inline void stimulus::quantize(double TlIn[], double dataIn[])
{
	int dataIndex;
	int TlIndex;

	quantizedData = new int[nofSamples];    // quantized signal

	for (dataIndex = 0; dataIndex < nofSamples; dataIndex++)
	{
		/*
		 * qBits number of bits decodes 2^qBits digital code
		 * from 0 to 2^qBits-1. There are one less code
		 * transition levels which means Tl goes from
		 * 1 to 2^qBits-1. T1 (the first code transition level)
		 * is stored in TlIn[0] hence the TlIndex goes
		 * from 0 to 2^qBits-2
		 */
		for (TlIndex = 0; TlIndex < (nofCodes-1); TlIndex++)
		{
			if (verbose >= 2 && nofSamples <= 1024)     // keep it for debug purposes
				fprintf(stdout, "dataIn[%d] = %f, TlIn[%d] = %f\n", dataIndex, dataIn[dataIndex], TlIndex, TlIn[TlIndex]);

			if (dataIn[dataIndex] <= TlIn[TlIndex])
			{
				quantizedData[dataIndex] = (TlIndex);
				break;
			}
		}

		if (TlIndex == (nofCodes-1))
			quantizedData[dataIndex] = (nofCodes-1);

		if (verbose >= 2 && nofSamples <= 1024)     // keep it for debug purposes
			fprintf(stdout, "quantizedData[%d] = %d\n---------\n", dataIndex, quantizedData[dataIndex]);
	}
}


#endif
