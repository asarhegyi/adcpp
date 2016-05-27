#ifndef SPECTRUM_CLASS
#define SPECTRUM_CLASS

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_math.h>
/************************************************************
* spectrum class											*
*		calculates the FFT of the data vector that gets		*
*		passed in. The member functions return with the		*
*		magnitude spectrum, the maximum magnitude and the	*
*		corresponding index of the magnitude array.			*
*		The interpolated FFT (IpFFT) is also calculated for	*
*		the fundamental frequency.							*
*															*
* Parameters:												*
*		data -- the time domain samples; it is expected to	*
*				be a real data-set							*
*															*
* Member functions											*
*		findMaxMagnitude -- returns the maximum magnitude of*
*					the spectrum and the corresponding index*
*					of the magnitude array					*
************************************************************/
class spectrum {
	public:
		int verbose;              // verbosity level
		std::vector<double> DFT;		// the half-complex spectrum of the real data sequence
		std::vector<double> magnitude;  // the magnitude spectrum 2*sqrt(REAL^2 + IMAG^2)/numSamples,
		                                // except magnitude[0]=REAL[0], calculate only 0 - (numSamples/2)-1
		int numSamples;           // number of data samples
		int prevpow2;             // the previous power of 2 from the size of the data vector
		                          // f.i. numSamples = 24512 -> prevpow2 = 16384					
		int nextpow2;             // the previous power of 2 from the size of the data vector
		                          // f.i. numSamples = 24512 -> nextpow2 = 32768
		int maxIndex;		      // the index of the maximum element of the magnitude spectrum
		double maxMagnitude[2];   // store the maxMagnitude and its index in this array


		// Constructor
		spectrum(std::vector<double>& data);

		//spectrum(const spectrum& old_spectrum);
			// Use automatically generated copy constructor

		//spectrum operator = (const spectrum& old_spectrum);
			// Use automatically generated assign constructor

		//~spectrum();
			// Use automatically generated destructor
		
		// return the maximum magnitude and the corresponding index of the spectrum
		double* findMaxMagnitude();

		// calculate the Interpolated FFT (IpFFT)
		double IpFFT();
};


/************************************************************
* spectrum::spectrum -- initialize the spectrum class by	*
* calculating the half-complex DFT of the real data-set,	*
* the magnitude, and determining whether the number of		*
* samples falls on a radix2 bin.							*
************************************************************/
inline spectrum::spectrum(std::vector<double>& data)
{
	verbose = 0;
	DFT = data;
	numSamples = DFT.size();
	magnitude.resize(numSamples/2);

	prevpow2 = pow(2, (int)log2 (numSamples));
	nextpow2 = pow(2, ((int)log2 (numSamples)+1));

	if (verbose >= 1)
	{
		std::cout << std::endl;
		std::cout << "Size of the data vector: " << numSamples << std::endl;
		std::cout << "The previous power of 2 from the size of the data vector: " << prevpow2 << std::endl;
		std::cout << "The next power of 2 from the size of the data vector: " << nextpow2 << std::endl;
		std::cout << std::endl;
	} 

	if ( numSamples == prevpow2 )
	{
		if (verbose >= 1)
			std::cout << "Number of Samples: " << numSamples << " => using Radix-2 FFT routines" << std::endl;

		gsl_fft_real_radix2_transform (&DFT[0], 1, numSamples);

		// Calculate the magnitude of the spectrum: 2*sqrt(REAL^2 + IMAG^2)/numSamples
		int index = 0;
		while ( index < numSamples/2 )
		{
			if ((index == 0) || (index == numSamples/2))
				magnitude[index] = DFT[index]/numSamples;
			else
				magnitude[index] = 2*gsl_hypot(DFT[index], DFT[(numSamples-index)])/numSamples;

			++index;
		}
	}
	else
	{
		if (verbose >= 1)
			std::cout << "Number of Samples: " << numSamples << " => using Mixed-radix FFT routines" << std::endl;

		gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc (numSamples);
		gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc (numSamples); 

		gsl_fft_real_transform (&DFT[0], 1, numSamples, real, work);
	
		gsl_fft_real_wavetable_free (real);
		gsl_fft_real_workspace_free (work);

		// Calculate the magnitude of the spectrum: 2*sqrt(REAL^2 + IMAG^2)/numSamples
		int index = 0;
		while ( index < numSamples/2 )
		{
			if ( index == 0 )
				magnitude[index] = DFT[index]/numSamples;
			else
				magnitude[index] = 2*gsl_hypot(DFT[(2*index)-1], DFT[(2*index)])/numSamples;

			++index;
		}
	}
}


/************************************************************
* spectrum::findMaxMagnitude --returns with a pointer to the*
* maxMagnitude array which stores the maximum magnitude of	*
* the spectrum and the corresponding index of the magnitude	*
* array.													*
************************************************************/
inline double* spectrum::findMaxMagnitude()
{
	double *maxMagnitude_ptr;
	maxMagnitude_ptr = maxMagnitude;

	maxIndex = 1;
	int index = 0;
	while ( index < (int)magnitude.size() )
	{
		// find the maximum and ignore the DC component
		if (index != 0)
		{
			if ( magnitude[index] > magnitude[maxIndex] )
			{
				maxIndex = index;
			}
		}
		++index;
	}

	maxMagnitude[0] = magnitude[maxIndex];
	maxMagnitude[1] = maxIndex;
	return maxMagnitude_ptr;
}


/************************************************************
* spectrum::IpFFT -- calculate the interpolated FFT	based on*
* the following article:									*
* H. Renders, J. Schoukens, and G. Vilain, “High-Accuracy	*
* Spectrum Analysis of Sampled Discrete Frequency Signals	*
* by Analytical Leakage Compensation” in IEEE Transaction on*
* Instrumentation and Measurement, VOL. IM-33, NO. 4,		*
* December 1984												*
************************************************************/
inline double spectrum::IpFFT()
{
	int index;      // spectral line index
	double U[2];	// the REAL term of the two larges spectral lines
	double V[2];	// the IMAGINARY term of the two larges spectral lines
	int spectLine;	// spectral line counter
	double n, Kopt, Z1, Z2, lambda;  // variables for calculating the IpFFT

	findMaxMagnitude();

	// choose the largest two spectral line  
	if ( magnitude[maxIndex-1] > magnitude[maxIndex+1] )
		index = maxIndex-1;
	else
		index = maxIndex;

	if (verbose >= 1)
	{
		std::cout << std::endl;
		std::cout << "	maxIndex     : " << maxMagnitude[1] << std::endl;
		std::cout << "	magnitude[" << maxIndex << "]: " << maxMagnitude[0] << std::endl;
		std::cout << "	magnitude[" << maxIndex-1 << "]: " << magnitude[maxIndex-1] << std::endl;
		std::cout << "	magnitude[" << maxIndex+1 << "]: " << magnitude[maxIndex+1] << std::endl;
		std::cout << "	IpFFT index  : " << index << std::endl;
		std::cout << std::endl;
	}

	// load the real and imaginary part of the two largest spectral lines (index; index+1) into U[], and V[]
	if ( numSamples == prevpow2 )  // using Radix-2 FFT routines
	{
		for ( spectLine = 0; spectLine < 2; spectLine++ )
		{
			if ( (index == 0) || (index == numSamples/2) )
			{
				U[spectLine] = DFT[index];
				V[spectLine] = 0;
			}
			else
			{
				U[spectLine] = DFT[index];
				V[spectLine] = DFT[(numSamples-index)];
			}
			++index;
		}
	}
	else    // using Mixed-radix FFT routines
	{
		for ( spectLine = 0; spectLine < 2; spectLine++ )
		{
			if ( index == 0 )
			{
				U[spectLine] = DFT[index];
				V[spectLine] = 0;
			}
			else
			{
				U[spectLine] = DFT[(2*index)-1];
				V[spectLine] = DFT[(2*index)];
			}
			++index;
		}
	}

	index = index-2;  // restore index after loading U[] and V[]

	// calculate the Interpolated FFT
	n = 2*M_PIl/numSamples;
	Kopt = ( sin(n*index) * (V[1]-V[0]) + cos(n*index) * (U[1]-U[0]) ) / (U[1]-U[0]);
	Z1 = V[0] * ( Kopt-cos(n*index) ) / sin(n*index) + U[0];
	Z2 = V[1] * ( Kopt-cos(n*(index+1)) ) / sin(n*(index+1)) + U[1];
	lambda = acos( (Z2*cos(n*(index+1)) - Z1*cos(n*index)) / (Z2-Z1) )/n;

	return(lambda);
}

#endif
