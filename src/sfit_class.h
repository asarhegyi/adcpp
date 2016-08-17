#ifndef SFIT_CLASS
#define SFIT_CLASS

#include <cmath>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_statistics_double.h>

#include "utils.h"
#include "spectrum_class.h"
/************************************************************
* sfit class												*
*	This class contains the three and four parameter		*
*	sine-wave fitting algorithms which minimizes			*
*	chisq = ||y - Xc||^2 where y is the observation vector	*
*	and Xc is the output of the fitted model.				*
*															*
*	For the 3-parameter sine-wave fitting algorithms each	*
*	sine-wave can be described by 2 coefficients + DC (the	*
*	frequency is known in this case) which results in		*
*	nofCoeffs=3+2*nofHarmonics coefficients to estimate.	*
*	For the 4-parameter sine-wave fitting algorithms the	*
*	fundamental frequency is unknown but all harmonics are	*
*	tied to this frequency which results in					*
*	nofCoeffs=(3+1)+2*nofHarmonics coefficients to estimate.*
*															*
*	Note that he indices of the sine-waves are shifted by 1	*
*	since the first harmonic is considered to be the		*
*	fundamental frequency "f" and the second harmonics is	*
*	at the frequency 2*f and so on.							*
*															*
* Parameters:												*
*		See parameters below.								*
*															*
* Member functions											*
*	sfit3_linear -- 3-parameter sine-wave fitting algorithm	*
*				for processing a small set of data (<32768).*
*															*
*	sfit3_large -- 3-parameter sine-wave fitting algorithm	*
*				for processing a large set of data.			*
************************************************************/
class sfit {
	public:
		std::string method;       // the sine-fit algorithm to run
		std::string mode;         // used to distinguish normal and test mode
		int verbose;              // verbosity level
		std::string inputFile;    // input file name
		std::vector<double> time; // the sample timing information
		std::vector<double> data; // the input samples are stored in this vector
		int nofSamples;           // the number of samples stored in the data vector
		int nofHarmonics;         // number of sine-wave harmonics
		double frequency;         // the fundamental frequency of the sine-wave
		double fs;                // sampling frequency
		double omega;             // angular frequency of the sine wave
		int nofCoeffs;            // number of coefficients to estimate

		double chisq;             // the weighted sum of squared residuals, the object of minimization 
		gsl_matrix *X;            // matrix of predictor variables
		gsl_vector_view y;        // observation vector is a gsl_vector_view of data
		gsl_vector *c;            // best fit coefficients; results of the LS-fit 
		gsl_matrix *cov;          // covariance matrix of the coefficients
		int blockSize;            // the size of the block to exhaust at a time for large set of data
		double rnorm, snorm;      // residual and solution norms for regularized LS (sfit3_large)

		double *amplitude;        // the estimated amplitude of the sine-wave in digital code
		int ampSize;              // number of element of the amplitude[] array
		double *phase;            // the estimated phase of the sine-wave in digital code
		double dc;                // the estimated dc of the sine-wave in digital code

		gsl_vector *residuals;    // vector of residuals : y-Xc
		double *yf, *yf_err;      // the fitted sine-wave and its standard deviation
		double tss;               // total sum of suqares tss = sum(x_i-mean)^2
		double Rsq;               // coefficient of determination : Rsq = 1-chisq/tss;
		double erms;              // erms = sqrt(chisq/nofSamples)

		// Use constructor to initialize macros and variables
		sfit();

		//sfit(const sfit& old_sfit);
			// Use automatically generated copy constructor

		//sfit operator = (const sfit& old_sfit);
			// Use automatically generated assign constructor

		// Use destructor to free-up the memory 
		~sfit();

		// 3-paramter linear sine-wave fit
		void sfit3_linear(std::string fileNameIn, double freqIn, int harmonicsIn, double fsIn);

		// 3-paramter linear sine-wave fit for large sets
		void sfit3_large(std::string fileNameIn, double freqIn, int harmonicsIn, double fsIn, int blockSizeIn);

		// convert the A,B,C coefficients to Amplitude, Pahse, and DC 
		void ABC2AmpPhaseDC();

		// calculate the vector of residuals = y-Xc for sfit3_linear
		void residual_linear();

		// calculate the fitted function values and its standard deviation for sfit3_linear
		void fitted_linear();

		// Use amplitude, phase, and dc to calculate the fitted sine-wave instead of A, B, and C 
		void fitted2_linear();

		// print coefficients to the standard output
		void print_coeffs();

		// print the parameters of the sine-wave
		void print_parameters();

		// print the covariance matrix from sfit3_linear
		void print_cov();

		// print chisq, degree of freedom, and coefficient of determination.
		void print_chisqRsq();

		// calculate the vector of residuals = y-Xc for sfit3_large
		void residual_large();

		// 4-paramter sine-wave fit
		void sfit4_linear(std::string fileNameIn, double freqIn, int harmonicsIn, double fsIn, std::string FFT, double freqError, int maxCycle );

		// 4-paramter sine-wave fit
		void sfit4_large(std::string fileNameIn, double freqIn, int harmonicsIn, double fsIn, int blockSizeIn, std::string FFT, double freqError, int maxCycle );

};


/************************************************************
* sfit::sfit3 -- constructor to initialize macros & variables*
************************************************************/
inline sfit::sfit()
{
	#define C(i) (gsl_vector_get(c,(i)))
	#define CV(i) (gsl_vector_get(&cv.vector,(i)))

	verbose = 0;
	nofSamples = 0;
	blockSize = 32768;
	chisq = 0;
	tss = 0;

	X = gsl_matrix_alloc (1,1);
	c = gsl_vector_alloc (1);
	cov = gsl_matrix_alloc (1,1);
	residuals = gsl_vector_alloc (1);
}


/************************************************************
* sfit::~sfit3 -- destructor to free up memory				*
************************************************************/
inline sfit::~sfit()
{
	//std::vector<double> data; supposed to be deallocated by its class destructor 
	gsl_matrix_free (X);
	gsl_vector_free (c);
	gsl_matrix_free (cov);
	gsl_vector_free (residuals);
}


/************************************************************
* sfit::sfit3_linear -- 3-parameter sine-wave fitting		*
*		algorithm using the gsl_multifit_linear workspace.	*
*		This method is	suitable for processing a small set	*
*		of data (<32768) since it loads the whole set into	*
*		the memory at the same time.						*
*		Note that this is non-iterative	least-squares method*
*		There is no need to do any initial estimate on the	*
*		coefficients.										*
*															*
* Parameters												*
*		fileNameIn -- the input file stored the observations*
*		harmonicsIn -- number of harmonics to estimate		*
*		freqIn -- fundamental frequancy of the sine-wave	*
*		fsIn   -- sampling frequency						*
*															*
* Return													*
*		the estimented coefficinets are stored in the		*
*		gsl_vector *c										*
************************************************************/
inline void sfit::sfit3_linear(std::string fileNameIn, double freqIn, int harmonicsIn, double fsIn)
{
	inputFile = fileNameIn;
	frequency = freqIn;
	nofHarmonics = harmonicsIn;
	fs = fsIn;
	method = "linear";
	double sine_arg;            // the argument of the base functions of X

	// Load the content of the input file into the data vector
	// and return with the number of elements in the vector.
	// use member functions from fileio class
	class fileio sfit3_fileio;
	sfit3_fileio.verbose = verbose;
	nofSamples = sfit3_fileio.load_data (inputFile, time, data, fs);

	if (verbose >= 2)
	{
		// Save off data for debug purposes
		dbl_array2file("input_data.dat", data, nofSamples, "");
		dbl_array2file("input_time.dat", data, nofSamples, "");
	}

	omega = 2*M_PIl*frequency;
	nofCoeffs = 3 + 2*nofHarmonics;

	X = gsl_matrix_alloc (nofSamples, nofCoeffs);  // matrix of predictor variables
	c = gsl_vector_alloc (nofCoeffs);              // best fit coefficients
	cov = gsl_matrix_alloc (nofCoeffs, nofCoeffs); // covariance matrix

	// Compose the matrix of predictor variables X
	for (int index = 0; index < nofSamples; index++)
	{
		// DC component
		gsl_matrix_set (X, index, 0, 1);

		// Sine-wave components
		for (int index2 = 0; index2 <= nofHarmonics; index2++)
		{
			//sine_arg = (index2+1) * omega * index/fs;
			sine_arg = (index2+1) * omega * time[index];
			gsl_matrix_set (X, index, 2*index2+1, cos(sine_arg));
			gsl_matrix_set (X, index, 2*index2+2, sin(sine_arg));
		}
	}

	// Create the observation vector y
	y = gsl_vector_view_array (&data[0], nofSamples);

	// Save off the predictor matrix and observation vector in verbose mode for debug
	if (verbose >= 2)
	{
		sfit3_fileio.save_matrix (X, "m_X.dat", "%.6e", "w");
		sfit3_fileio.save_vector (&y.vector, "v_y.dat", "%g", "w");
	}

	// Do the linear least squares fitting: c = inv(X^T*X) * X^T*y
	{
		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (nofSamples, nofCoeffs);
		gsl_multifit_linear (X, &y.vector, c, cov, &chisq, work);
		gsl_multifit_linear_free (work);
	}
}


/************************************************************
* sfit::sfit3_large -- 3 parameter sine-wave fitting		*
*		algorithm using the	gsl_multilarge_linear workspace.*
*		This method was	designed to exhaust a larger set of	*
*		data by reading	it in recursively. This algorithm	*
*		still solves the same LINEAR least-squares equation	*
*		as sfit3_linear.									*
*		Note that this is non-iterative	least-squares method*
*		There is no need to do any initial estimate on the	*
*		coefficients.										*
*
* Parameters												*
*		fileNameIn -- the input file stored the observations*
*		harmonicsIn -- number of harmonics to estimate		*
*		freqIn -- fundamental frequancy of the sine-wave	*
*		fsIn   -- sampling frequency						*
*		blockSizeIn -- the size of the data block to exhaust*
*					at a time								*
*															*
* Return													*
*		the estimented coefficinets are stored in the		*
*		gsl_vector *c										*
************************************************************/
inline void sfit::sfit3_large(std::string fileNameIn, double freqIn, int harmonicsIn, double fsIn, int blockSizeIn = 32768)
{
	inputFile = fileNameIn;
	frequency = freqIn;
	nofHarmonics = harmonicsIn;
	fs = fsIn;
	blockSize = blockSizeIn;

	// use member functions from fileio class
	class fileio sfit3_fileio;
	sfit3_fileio.verbose = verbose;

	method = "large";
	omega = 2*M_PIl*frequency;
	nofCoeffs = 3 + 2*nofHarmonics;

	gsl_multilarge_linear_workspace * work = gsl_multilarge_linear_alloc(gsl_multilarge_linear_normal, nofCoeffs);
	X = gsl_matrix_alloc (blockSize, nofCoeffs);   // matrix of predictor variables
	c = gsl_vector_alloc (nofCoeffs);              // best fit coefficients
	cov = gsl_matrix_alloc (nofCoeffs, nofCoeffs); // dummy covariance matrix

	int loadSize = 0;           // Size of data loaded from file; the last block can be smaller than blockSize
	bool openFlag = false;      // Indicate the state of the input file (already opened or closed)
	bool eofFlag = false;       // End of file flag for the input file
	std::ifstream inFile;       // Initialize an IO stream object for the input file

	size_t rowIndex = 0;        // used to iterate through the whole observation vector
	size_t rowLeft = 0;         // number of rows left to accumulate
	size_t nofRows = 0;         // number of rows in the current block
	size_t blockRowIndex = 0;   // row index in the current block
	double sine_arg;            // the argument of the base functions of X
	char fmode;                 // I/O mode of the output files (write and append modes utilized)
	fmode = 'w';                // write mode

	// block loop
	while (true)
	{
		// Load blockSize chunks of data intermittently from the input file 
		loadSize = sfit3_fileio.stream_data (openFlag, eofFlag, inFile, inputFile, time, data, fs, blockSize, nofSamples);
		nofSamples += loadSize;

		rowLeft = nofSamples - rowIndex;
		nofRows = GSL_MIN(blockSize, rowLeft);

		// Filter out corner case when blockSize == number of lines in the input file
		if (nofRows == 0 )
			break;

		// Create matrix and vector views of X and data for the block -> Xv and y
		y = gsl_vector_view_array (&data[0], nofRows);
		gsl_matrix_view Xv = gsl_matrix_submatrix (X, 0, 0, nofRows, nofCoeffs);

		// Build predictor matrix (X) blocks with nofRows
		for (blockRowIndex = 0; blockRowIndex < nofRows; blockRowIndex++)
		{
			// DC component of Xv
			gsl_matrix_set (&Xv.matrix, blockRowIndex, 0, 1);

			// Sine-wave components of Xv
			for (int index2 = 0; index2 <= nofHarmonics; index2++)
			{
				//sine_arg = (index2+1) * omega * (rowIndex+blockRowIndex)/fs;
				sine_arg = (index2+1) * omega * time[blockRowIndex];
				gsl_matrix_set (&Xv.matrix, blockRowIndex, 2*index2+1, cos(sine_arg));
				gsl_matrix_set (&Xv.matrix, blockRowIndex, 2*index2+2, sin(sine_arg));
			}
		}

		// accumulate (X,y) block into LS system
		gsl_multilarge_linear_accumulate(&Xv.matrix, &y.vector, work);

		if (verbose >= 1)    // print data for debug purposes
			fprintf (stdout, "Done with accumulating a %d line block.\n", nofRows);

		// Save off the predictor matrix and observation vector in verbose mode for debug
		if (verbose >= 2)
		{
			sfit3_fileio.save_matrix (&Xv.matrix, "m_X.dat", "%.6e", &fmode);
			sfit3_fileio.save_vector (&y.vector, "v_y.dat", "%g", &fmode);
			fmode = 'a';  // set the fmode to append after first time through
		}

		rowIndex += nofRows;

		// Exit while loop after the last block is processed
		if (eofFlag)
			break;
	}

	// solve large LS system and store solution in c
	double lambda = 0.0;          // regularization parameter
	gsl_multilarge_linear_solve(lambda, c, &rnorm, &snorm, work);
}



/************************************************************
* sfit::ABC2AmpPhaseDC -- convert the A,B,C parameters of	*
* the sine fit to Amplitude, Phase, and DC.					*
************************************************************/
inline void sfit::ABC2AmpPhaseDC()
{
	int cSize = static_cast<int> (c->size);
	ampSize = (cSize -1)/2;				// fundamental + harmonics
	amplitude = new double[ampSize];
	phase = new double[ampSize];

	for (int index = 0; index < ampSize; index++)
	{
		gsl_vector_view cv = gsl_vector_subvector(c, (2*index+1), 2);

		// convert amplitude
		amplitude[index] = gsl_blas_dnrm2 (&cv.vector);

		// convert phase
		if (CV(0) > 0)
			phase[index] = atan(-CV(1)/CV(0));
		else if (CV(0) < 0)
			phase[index] = atan(-CV(1)/CV(0))+M_PIl;
		else  // CV(0) == 0
		{
			if (CV(1) < 0)
				phase[index] = M_PIl/2;
			else
				phase[index] = 3*M_PIl/2;
		}
		//phase[index] *= (180/M_PIl);
	}

	// convert dc
	dc = C(0);
}


/************************************************************
* sfit::residual_linear --  compute the vector of residuals	*
* (residuals = y-Xc) for sfit3_linear.						*
* The whole X matrix and y vector are stored in memory.		*
* A standard GSL function will do the work.					*
************************************************************/
inline void sfit::residual_linear()
{
	if (method.compare("linear") != 0)
	{
		fprintf (stderr, "Use sfit::residual_linear() with sfit::sfit3_linear\nExiting...\n");
		exit(1);
	}
		
	residuals = gsl_vector_alloc (nofSamples);
	gsl_multifit_linear_residuals (X, &y.vector, c, residuals);

	if (mode.compare("test_time") != 0)
	{
		dbl_array2file("residuals.dat", residuals->data, nofSamples, "");
	}
}


/************************************************************
* sfit::fitted_linear -- compute the fitted function values	*
* yf, its standard deviation yf_err, and save the waveforms	*
* to a file.												*
************************************************************/
inline void sfit::fitted_linear()
{
	if (method.compare("linear") != 0)
	{
		fprintf (stderr, "Use sfit::fitted_linear() with sfit::sfit3_linear\nExiting...\n");
		exit(1);
	}

	int index;					// used to loop through the dataset
	gsl_vector_view X_row;		// vector view of the rows of matrix X
	yf = new double[nofSamples];
	yf_err= new double[nofSamples];

	index = 0;
	while (true)
	{
		// Utilize the row view of predictor matrix X
		X_row = gsl_matrix_row (X, index);
		
		// Calculate the fitted value yf and its standard deviation yf_err
		gsl_multifit_linear_est (&X_row.vector, c, cov, &yf[index], &yf_err[index]);

		// Exit condition
		if (index >= (nofSamples-1))
			break;

		++index;
	}

	// Save the fitted waveform to file
	if (mode.compare("test_time") != 0)
	{
		dbl_array2file("fitted.dat", yf, nofSamples, "");
		//dbl_array2file("fitted_err.dat", yf_err, nofSamples, "");
	}
}


/************************************************************
* sfit::fitted2_linear -- Use amplitude, phase, and dc to	*
* calculate the fitted sine-wave instead of A, B, and C to	*
* verify the ABC2AmpPhaseDC conversion.						*
************************************************************/
inline void sfit::fitted2_linear()
{
	double yfm[nofSamples];

	for(int index = 0; index < nofSamples; index++)
	{
		yfm[index] = dc;

		for (int harmonicsIndx = 0; harmonicsIndx < ampSize; harmonicsIndx++)
			yfm[index] += amplitude[harmonicsIndx] * cos((harmonicsIndx+1)*omega*time[index] + phase[harmonicsIndx]);
	}

	dbl_array2file("fitted2.dat", yfm, nofSamples, "");
}


/************************************************************
* sfit::print_coeffs -- print the coefficients estimated	*
* by Least-Squares to the standard output					*
************************************************************/
inline void sfit::print_coeffs()
{
	int cSize = static_cast<int> (c->size);
	bool isDeltaOmega = false; // set true if 4-param. sfit is active

	if ((cSize%2) == 0)        // 4-parameter sine-fit
	{
		isDeltaOmega = true;
		--cSize;
	}

	fprintf(stdout, "Best fit results:\n");
	for (int index2 = 0; index2 < cSize; index2++)
	{
		if (index2 == 0)
			fprintf (stdout, "\tC : %.20e\n", C(index2));
		else if ((index2%2) == 0)
			fprintf (stdout, "\tB%d: %.20e\n", (index2-1)/2, C(index2));
		else
			fprintf (stdout, "\tA%d: %.20e\n", (index2-1)/2, C(index2));
	}

	if ( isDeltaOmega )
		fprintf (stdout, "\tdelta_omega: %.20e\n", C(cSize));

	fprintf (stdout, "\n");
}


/************************************************************
* sfit::print_parameters -- print the parameters of the		*
* sine-wave such as amplitude, frequency, phase, and dc.	*
************************************************************/
inline void sfit::print_parameters()
{
	fprintf(stdout, "Parameters of the sine-wave:\n");
	// print amplitude and phase of the fundamental and harmonics
	for (int index = 0; index < ampSize; index++)
	{
		fprintf(stdout, "\tAmp%d  : %.20e\n", index, amplitude[index]);
		fprintf(stdout, "\tPhase%d: %.20e\n", index, phase[index]);
	}

	// Plot frequency
	fprintf(stdout, "\tFreq  : %.20e\n", frequency);

	// print dc
	fprintf(stdout, "\tdc    : %.20e\n", dc);
	fprintf (stdout, "\n");
}


/************************************************************
* sfit::print_cov -- print the covariance matrix from		*
* sfit3_linear												*
************************************************************/
inline void sfit::print_cov()
{
	if (method.compare("linear") != 0)
	{
		fprintf (stderr, "Use sfit::print_cov() with sfit::sfit3_linear\nExiting...\n");
		exit(1);
	}

	#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

	fprintf (stdout, "Covariance matrix:\n");
	fprintf (stdout, "[ %+.5e, %+.5e, %+.5e\n",   COV(0,0), COV(0,1), COV(0,2));
	fprintf (stdout, "  %+.5e, %+.5e, %+.5e\n",   COV(1,0), COV(1,1), COV(1,2));
	fprintf (stdout, "  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));
	fprintf (stdout, "\n");
}


/************************************************************
* sfit::print_chisqRsq -- print chisq, degree of freedom,	*
* and coefficient of determination.							*
************************************************************/
inline void sfit::print_chisqRsq()
{
	if (method.compare("linear") == 0)
	{
		tss = gsl_stats_tss(&y.vector.data[0], 1, *(&y.vector.size));
	}
fprintf(stdout, "\ttss: %.20e\n", tss);

	Rsq = 1.0 - chisq / tss;
	erms = sqrt(chisq/nofSamples);
	
	fprintf(stdout, "Misc. parameters:\n");
	fprintf(stdout, "\tchisq: %.20e\n", chisq);
	fprintf(stdout, "\tdof: %d\n", nofSamples - nofCoeffs);
	fprintf(stdout, "\tchisq/dof: %.20e\n", chisq / (nofSamples - nofCoeffs));
	fprintf(stdout, "\tRsq: %.20e\n", Rsq);
	fprintf(stdout, "\terms: %.20e\n", erms);
} 


/************************************************************
* sfit::residual_large --  compute the vector of residuals	*
* (residuals = y-Xc) for sfit3_large.						*
* The observations vector (y) is read in recursively due to	*
* its size. The whole residuals vector is saved in			*
* residuals.dat.											*
* The assumption was made that most of the parameters,		*
* vectors, and matrices were already defined.				*
************************************************************/
inline void sfit::residual_large()
{
	if (method.compare("large") != 0)
	{
		fprintf (stderr, "Use residual_large with sfit::sfit3_large\nExiting...\n");
		exit(1);
	}

	int loadSize = 0;           // Size of data loaded from file; the last block can be smaller than blockSize
	bool openFlag = false;      // Indicate the state of the input file (already opened or closed)
	bool eofFlag = false;       // End of file flag for the input file
	std::ifstream inFile;       // Initialize an IO stream object for the input file

	size_t rowIndex = 0;        // used to iterate through the whole observation vector
	size_t rowLeft = 0;         // number of rows left to accumulate
	size_t nofRows = 0;         // number of rows in the current block
	size_t blockRowIndex = 0;   // row index in the current block
	double sine_arg;            // the argument of the base functions of X
	char fmode;                 // I/O mode of the output files (write and append modes utilized)
	fmode = 'w';                // set fmode

	yf = new double[blockSize];
	residuals = gsl_vector_alloc (blockSize);

	int nofSamples_backup = nofSamples; // save off nofSamples before post-processing the sfit results
	nofSamples = 0;                     // make sure nofSamples_backup == nofSamples at the end of post-processing

	// use member functions from fileio class
	class fileio sfit_fileio;
	sfit_fileio.verbose = verbose;

	// block loop
	while (true)
	{
		// Load blockSize chunks of data intermittently from the input file 
		loadSize = sfit_fileio.stream_data (openFlag, eofFlag, inFile, inputFile, time, data, fs, blockSize, nofSamples);
		nofSamples += loadSize;

		rowLeft = nofSamples - rowIndex;
		nofRows = GSL_MIN(blockSize, rowLeft);

		// Filter out corner case when blockSize == number of lines in the input file
		if (nofRows == 0 )
			break;

		// Create matrix and vector views of X and data for the block => Xv and y
		gsl_matrix_view Xv = gsl_matrix_submatrix (X, 0, 0, nofRows, nofCoeffs);
		y = gsl_vector_view_array (&data[0], nofRows);

		// Build predictor matrix (X) blocks with nofRows
		for (blockRowIndex = 0; blockRowIndex < nofRows; blockRowIndex++)
		{
			// DC component of Xv
			gsl_matrix_set (&Xv.matrix, blockRowIndex, 0, 1);

			// Sine-wave components of Xv
			for (int index2 = 0; index2 <= nofHarmonics; index2++)
			{
				//sine_arg = (index2+1) * omega * (rowIndex+blockRowIndex)/fs;
				sine_arg = (index2+1) * omega * time[blockRowIndex];
				gsl_matrix_set (&Xv.matrix, blockRowIndex, 2*index2+1, cos(sine_arg));
				gsl_matrix_set (&Xv.matrix, blockRowIndex, 2*index2+2, sin(sine_arg));
			}
		}

		// the fitted sine-wave yf = X*c
		gsl_vector_view yfv = gsl_vector_view_array (yf, nofRows);
		gsl_vector_set_zero(&yfv.vector);  // clear yfv
		float alpha = 1.0;
		float beta = 1.0;
		// yfv = alpha*X*c + beta*yfv (yfv cleared above)
		gsl_blas_dgemv(CblasNoTrans, alpha, &Xv.matrix, c, beta, &yfv.vector);

		if (mode.compare("test_time") != 0)
		{
			// Save off the fitted sine-wave
			FILE *yfd;
			yfd = fopen ("fitted.dat", &fmode);
			gsl_vector_fprintf (yfd, &yfv.vector, "%g");
			fclose(yfd);
		}

		// residuals = y-X*c
		gsl_vector_view resv = gsl_vector_subvector (residuals, 0, nofRows);
		gsl_vector_memcpy (&resv.vector, &y.vector);  // y -> resv
		gsl_vector_sub (&resv.vector, &yfv.vector);  // resv = resv -yfv (yfv=Xc) => residuals = y-Xc

		if (mode.compare("test_time") != 0)
		{
			// Save off residuals
			FILE *resfd;
			resfd = fopen ("residuals.dat", &fmode);
			gsl_vector_fprintf (resfd, &resv.vector, "%g");
			fclose(resfd);
		}

		// resq = resv^T*resv; the sum of resq for all blocks equals to
		// chisq which should be equal to rnorm^2
		double resq;
		gsl_blas_ddot(&resv.vector, &resv.vector, &resq);
		chisq += resq;

		// calculate tss for Rsq 
		tss += gsl_stats_tss(&y.vector.data[0], 1, *(&y.vector.size));

		if (verbose >= 1)    // print data for debug purposes
			fprintf (stdout, "Done with processing a %d line block.\n", nofRows);

		rowIndex += nofRows;

		// First time through; set the fmode to append
		fmode = 'a';

		// Exit while loop after the last block is processed
		if (eofFlag)
			break;
	}

	if (nofSamples_backup != nofSamples)
	{
		fprintf (stderr, "nofSamples_backup != nofSamples");
		exit(1);
	}
}


/************************************************************
* sfit::sfit4_linear -- 4-parameter sine-wave fitting		*
*		algorithm using the gsl_multifit_linear workspace.	*
*		This method is	suitable for processing a small set	*
*		of data (<32768) since it loads the whole set into	*
*		the memory at the same time.						*
*		Note that this is an iterative least-squares method	*
*		which requires to do an initial estimate on the		*
*		coefficients.										*
*															*
* Parameters												*
*		fileNameIn -- the input file stored the observations*
*		harmonicsIn -- number of harmonics to estimate		*
*		freqIn -- fundamental frequency of the sine-wave	*
*		fsIn   -- sampling frequency						*
*		FFT -- defines the method to estimate the initial	*
*			value of the fundamental frequency.				*
*			Options: "FFT", "IpFFT", "input"				*
*		freqError -- the threshold of the frequency error.	*
*			If the frequency change between two consecutive	*
*			steps is smaller than this threshold, the		*
*			iteration terminates. Set to 10^-6 by default.	*
*		maxCycle -- specifies the maximum number of			*
*			iterations. It terminates the iteration even if	*
*			the frequency error is higher than freqError.	*
*															*
* Return													*
*		the estimated coefficients are stored in the		*
*		gsl_vector *c										*
************************************************************/
inline void sfit::sfit4_linear(std::string fileNameIn, double freqIn, int harmonicsIn, double fsIn, std::string FFT = "IpFFT", double freqError = pow(10,-6), int maxCycle = 30 )
{
	inputFile = fileNameIn;
	nofHarmonics = harmonicsIn;
	fs = fsIn;
	method = "linear";
	double sine_arg;            // the argument of the base functions of X
	double omegaError;			// stopping criteria 1: the smallest steps size in omega before exiting the iteration
	struct estimate				// store the results of the sfit4 iteration
	{
		std::vector<double> A_iter;
		std::vector<double> B_iter;
		std::vector<double> C_iter;
		std::vector<double> f_iter;
	} est;


	// Load the content of the input file into the data vector
	// and return with the number of elements in the vector.
	// use member functions from fileio class
	class fileio sfit4_fileio;
	sfit4_fileio.verbose = verbose;
	nofSamples = sfit4_fileio.load_data (inputFile, time, data, fs);

	if (verbose >= 1)
		std::cout << "Starting sfit4_linear with calculating the FFT." << std::endl;

	// instantiate a spectrum class to run FFT
	// for the initial estimate of the frequency
	class spectrum sfit4_fft(data);
	sfit4_fft.verbose = verbose;

	if ( FFT.compare("FFT") == 0 )
	{
		double *maxMag;
		maxMag = sfit4_fft.findMaxMagnitude();
		frequency = maxMag[1]*(fs/nofSamples);

		if (verbose >= 1)
		{
			std::cout << std::endl;
			std::cout << "	maxIndex     : " << maxMag[1] << std::endl;
			std::cout << "	magnitude[" << maxMag[1] << "]: " << maxMag[0] << std::endl;
			std::cout << "	maxFrequency : " << frequency << std::endl;
			std::cout << "	magnitude[0] : " << sfit4_fft.magnitude[0] << std::endl;
			std::cout << std::endl;
		}
	}
	else if ( FFT.compare("IpFFT") == 0 )
	{
		double lambda;
		lambda = sfit4_fft.IpFFT();
		frequency = lambda*(fs/nofSamples);
	}
	else
		frequency = freqIn;

	// run 3-paramter fit and get the initial estimate of A, B, and C
	sfit3_linear(inputFile, frequency, nofHarmonics, fs);

	est.C_iter.push_back (C(0));
	est.A_iter.push_back (C(1));
	est.B_iter.push_back (C(2));
	est.f_iter.push_back (frequency);

	if (verbose >= 1)
	{
		fprintf (stdout, "Cyc#	 frequency[Hz]	    Change of frequency  	     A		     B		     C     \n");
		fprintf (stdout, "----	 -------------	  -----------------------	 ----------	 ----------	 ----------\n");
		fprintf (stdout, " %02d\t%14.14g\t%+25.25g\t%+10.10g\t%+10.10g\t%+10.10g\n",\
			   	0, est.f_iter[0], est.f_iter[0], est.A_iter[0], est.B_iter[0], est.C_iter[0]);
	}


	omega = 2*M_PIl*frequency;
	nofCoeffs = 4 + 2*nofHarmonics;      // number of parameter: A0, B0, DC, frequency = 4+2*harmonics 
	omegaError = (2*M_PIl) * freqError;  // stopping criteria 1: smallest steps size before exiting the iteration

	X = gsl_matrix_alloc (nofSamples, nofCoeffs);  // matrix of predictor variables
	c = gsl_vector_alloc (nofCoeffs);              // best fit coefficients
	cov = gsl_matrix_alloc (nofCoeffs, nofCoeffs); // covariance matrix
	
	int iCycle = 0;
	
	do
	{
		// Compose the matrix of predictor variables X
		for (int index = 0; index < nofSamples; index++)
		{
			// DC component
			gsl_matrix_set (X, index, 0, 1);

			// Sine-wave components
			for (int index2 = 0; index2 <= nofHarmonics; index2++)
			{
				sine_arg = (index2+1) * omega * time[index];
				gsl_matrix_set (X, index, 2*index2+1, cos(sine_arg));
				gsl_matrix_set (X, index, 2*index2+2, sin(sine_arg));
			}
	
			// the delta_omega component
			double grad;
			grad=-est.A_iter[iCycle]*time[index]*sin(omega*time[index])+\
				  est.B_iter[iCycle]*time[index]*cos(omega*time[index]);
			gsl_matrix_set (X, index, 3+2*nofHarmonics, grad);
		}

		// Create the observation vector y
		y = gsl_vector_view_array (&data[0], nofSamples);

		// Save off the predictor matrix and observation vector in verbose mode for debug
		if (verbose >= 2)
		{
			sfit4_fileio.save_matrix (X, "m_X.dat", "%.6e", "w");
			sfit4_fileio.save_vector (&y.vector, "v_y.dat", "%g", "w");
		}

		// Do the linear least squares fitting: c = inv(X^T*X) * X^T*y
		{
			gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (nofSamples, nofCoeffs);
			gsl_multifit_linear (X, &y.vector, c, cov, &chisq, work);
			gsl_multifit_linear_free (work);
		}

		// adjust the step size of omega to minimize oscillation around the minimum
		if ( std::abs( C(3+2*nofHarmonics) ) < 0.05*omega )
			omega = omega + C(3+2*nofHarmonics);
		else
			omega = omega + 0.05*C(3+2*nofHarmonics);

		// update frequency
		frequency = omega/(2*M_PIl);

		est.C_iter.push_back (C(0));
		est.A_iter.push_back (C(1));
		est.B_iter.push_back (C(2));
		est.f_iter.push_back (frequency);

		++iCycle;

		if (verbose >= 1)
		{
			double delta_f = C(3+2*nofHarmonics)/(2*M_PIl);
			fprintf (stdout, " %02d\t%14.14g\t%+21.21g\t%+10.10g\t%+10.10g\t%+10.10g\n",\
					iCycle, est.f_iter[iCycle], delta_f, est.A_iter[iCycle], est.B_iter[iCycle], est.C_iter[iCycle]);
		}

	} while ( ( iCycle < maxCycle ) && ( std::abs( C(3+2*nofHarmonics) ) > omegaError ) );

	fprintf (stdout, "\n");
}

/************************************************************
* sfit::sfit4_large -- 4-parameter sine-wave fitting		*
*		algorithm using the	gsl_multilarge_linear workspace.*
*		This method was	designed to exhaust a larger set of	*
*		data by reading	it in recursively. This algorithm	*
*		still solves the same least-squares equation as		*
*		sfit4_linear.										*
*		Note that this is an iterative least-squares method	*
*		which requires to do an initial estimate on the		*
*		coefficients.										*
*															*
* Parameters												*
*		fileNameIn -- the input file stored the observations*
*		harmonicsIn -- number of harmonics to estimate		*
*		freqIn -- fundamental frequency of the sine-wave	*
*		fsIn   -- sampling frequency						*
*		blockSizeIn -- the size of the data block to exhaust*
*					at a time								*
*		FFT -- defines the method to estimate the initial	*
*			value of the fundamental frequency.				*
*			Options: "FFT", "IpFFT", "input". Set to IpFFT	*
*			by default.										*
*		freqError -- the threshold of the frequency error.	*
*			If the frequency change between two consecutive	*
*			steps is smaller than this threshold, the		*
*			iteration terminates. Set to 10^-6 by default.	*
*		maxCycle -- specifies the maximum number of			*
*			iterations. It terminates the iteration even if	*
*			the frequency error is higher than freqError.	*
*															*
* Return													*
*		the estimented coefficinets are stored in the		*
*		gsl_vector *c										*
************************************************************/
inline void sfit::sfit4_large(std::string fileNameIn, double freqIn, int harmonicsIn, double fsIn, int blockSizeIn, std::string FFT = "IpFFT", double freqError = pow(10,-6), int maxCycle = 30 )
{
	inputFile = fileNameIn;
	frequency = freqIn;
	nofHarmonics = harmonicsIn;
	fs = fsIn;
	blockSize = blockSizeIn;

	method = "large";
	double omegaError;			// stopping criteria 1: the smallest steps size in omega before exiting the iteration
	struct estimate				// store the results of the sfit4 iteration
	{
		std::vector<double> A_iter;
		std::vector<double> B_iter;
		std::vector<double> C_iter;
		std::vector<double> f_iter;
	} est;

	// use member functions from fileio class
	class fileio sfit4_fileio;
	sfit4_fileio.verbose = verbose;

	// count number of lines in the input file
	int numLines = sfit4_fileio.count_lines(inputFile);
	std::cout << "Number of Lines : " << numLines << std::endl;

	// calculate decimation ratio
	// set maximum size for FFT to 32768
	int halfFFTsize = 16384;
	int decRatio = numLines/halfFFTsize;
	std::cout << "Minimum (half) FFT size: "<< halfFFTsize << std::endl;
	std::cout << "Decimation ratio: " << decRatio << std::endl;

	// load data through recursive running sum decimation filter
	time.resize(2*halfFFTsize);
	data.resize(2*halfFFTsize);
	nofSamples = sfit4_fileio.load_N_decimate(inputFile, time, data, fs, decRatio);
	std::cout << "Number of decimated samples:  " << nofSamples << std::endl;

	if (verbose >= 2)
	{
		// Save off data for debug purposes
		dbl_array2file("rsum_out_data.dat", data, nofSamples, "");
		dbl_array2file("rsum_out_time.dat", time, nofSamples, "");
	}

	// instantiate a spectrum class to run FFT
	// for the initial estimate of the frequency
	class spectrum sfit4_fft(data);
	sfit4_fft.verbose = verbose;

	if ( FFT.compare("FFT") == 0 )
	{
		double *maxMag;
		maxMag = sfit4_fft.findMaxMagnitude();
		frequency = maxMag[1]*(fs/nofSamples/decRatio);

		if (verbose >= 1)
		{
			std::cout << std::endl;
			std::cout << "	maxIndex     : " << maxMag[1] << std::endl;
			std::cout << "	magnitude[" << maxMag[1] << "]: " << maxMag[0] << std::endl;
			std::cout << "	maxFrequency : " << frequency << std::endl;
			std::cout << "	magnitude[0] : " << sfit4_fft.magnitude[0] << std::endl;
			std::cout << std::endl;
		}
	}
	else if ( FFT.compare("IpFFT") == 0 )
	{
		double lambda;
		lambda = sfit4_fft.IpFFT();
		frequency = lambda*(fs/nofSamples/decRatio);
	}
	else
		frequency = freqIn;

	// run 3-paramter fit and get the initial estimate of A, B, and C
	sfit3_large(inputFile, frequency, nofHarmonics, fs);

	est.C_iter.push_back (C(0));
	est.A_iter.push_back (C(1));
	est.B_iter.push_back (C(2));
	est.f_iter.push_back (frequency);

	if (verbose >= 1)
	{
		fprintf (stdout, "Cyc#	 frequency[Hz]	    Change of frequency  	     A		     B		     C     \n");
		fprintf (stdout, "----	 -------------	  -----------------------	 ----------	 ----------	 ----------\n");
		fprintf (stdout, " %02d\t%14.14g\t%+25.25g\t%+10.10g\t%+10.10g\t%+10.10g\n",\
			   	0, est.f_iter[0], est.f_iter[0], est.A_iter[0], est.B_iter[0], est.C_iter[0]);
	}

	omega = 2*M_PIl*frequency;
	nofCoeffs = 4 + 2*nofHarmonics;      // number of parameter: A0, B0, DC, frequency = 4+2*harmonics 
	omegaError = (2*M_PIl) * freqError;  // stopping criteria 1: smallest steps size before exiting the iteration

	gsl_multilarge_linear_workspace * work = gsl_multilarge_linear_alloc(gsl_multilarge_linear_normal, nofCoeffs);
	X = gsl_matrix_alloc (blockSize, nofCoeffs);   // matrix of predictor variables
	c = gsl_vector_alloc (nofCoeffs);              // best fit coefficients
	cov = gsl_matrix_alloc (nofCoeffs, nofCoeffs); // dummy covariance matrix

	int loadSize = 0;           // Size of data loaded from file; the last block can be smaller than blockSize
	bool openFlag = false;      // Indicate the state of the input file (already opened or closed)
	bool eofFlag = false;       // End of file flag for the input file
	std::ifstream inFile;       // Initialize an IO stream object for the input file

	size_t rowIndex = 0;        // used to iterate through the whole observation vector
	size_t rowLeft = 0;         // number of rows left to accumulate
	size_t nofRows = 0;         // number of rows in the current block
	size_t blockRowIndex = 0;   // row index in the current block
	double sine_arg;            // the argument of the base functions of X
	char fmode;                 // I/O mode of the output files (write and append modes utilized)
	fmode = 'w';                // write mode


	int iCycle = 0;

	do
	{
		nofSamples =0;   // reset number of samples after using it during the decimation
		loadSize = 0;
		openFlag = false;
		eofFlag = false;
		rowIndex = 0;
		rowLeft = 0;
		nofRows = 0;
		blockRowIndex = 0;
		fmode = 'w';

		// block loop
		// read in large data file recursively
		while (true)
		{
			// Load blockSize chunks of data intermittently from the input file 
			loadSize = sfit4_fileio.stream_data (openFlag, eofFlag, inFile, inputFile, time, data, fs, blockSize, nofSamples);
			nofSamples += loadSize;

			rowLeft = nofSamples - rowIndex;
			nofRows = GSL_MIN(blockSize, rowLeft);

			// Filter out corner case when blockSize == number of lines in the input file
			if (nofRows == 0 )
				break;

			// Create matrix and vector views of X and data for the block -> Xv and y
			y = gsl_vector_view_array (&data[0], nofRows);
			gsl_matrix_view Xv = gsl_matrix_submatrix (X, 0, 0, nofRows, nofCoeffs);

			// Build predictor matrix (X) blocks with nofRows
			for (blockRowIndex = 0; blockRowIndex < nofRows; blockRowIndex++)
			{
				// DC component of Xv
				gsl_matrix_set (&Xv.matrix, blockRowIndex, 0, 1);

				// Sine-wave components of Xv
				for (int index2 = 0; index2 <= nofHarmonics; index2++)
				{
					sine_arg = (index2+1) * omega * time[blockRowIndex];
					gsl_matrix_set (&Xv.matrix, blockRowIndex, 2*index2+1, cos(sine_arg));
					gsl_matrix_set (&Xv.matrix, blockRowIndex, 2*index2+2, sin(sine_arg));
				}

				// the delta_omega component of Xv
				double grad;
				grad=-est.A_iter[iCycle]*time[blockRowIndex]*sin(omega*time[blockRowIndex])+\
					  est.B_iter[iCycle]*time[blockRowIndex]*cos(omega*time[blockRowIndex]);
				gsl_matrix_set (&Xv.matrix, blockRowIndex, 3+2*nofHarmonics, grad);
			}

			// accumulate (X,y) block into LS system
			gsl_multilarge_linear_accumulate(&Xv.matrix, &y.vector, work);

			if (verbose >= 1)    // print data for debug purposes
				fprintf (stdout, "Done with accumulating a %d line block.\n", nofRows);

			// Save off the predictor matrix and observation vector in verbose mode for debug
			if (verbose >= 2)
			{
				sfit4_fileio.save_matrix (&Xv.matrix, "m_X.dat", "%.6e", &fmode);
				sfit4_fileio.save_vector (&y.vector, "v_y.dat", "%g", &fmode);
				fmode = 'a';  // set the fmode to append after first time through
			}

			rowIndex += nofRows;

			// Exit while loop after the last block is processed
			if (eofFlag)
				break;
		}

		// solve large LS system and store solution in c
		double lambda = 0.0;          // regularization parameter
		gsl_multilarge_linear_solve(lambda, c, &rnorm, &snorm, work);


		// adjust the step size of omega to minimize oscillation around the minimum
		if ( std::abs( C(3+2*nofHarmonics) ) < 0.05*omega )
			omega = omega + C(3+2*nofHarmonics);
		else
			omega = omega + 0.05*C(3+2*nofHarmonics);

		// update frequency
		frequency = omega/(2*M_PIl);

		est.C_iter.push_back (C(0));
		est.A_iter.push_back (C(1));
		est.B_iter.push_back (C(2));
		est.f_iter.push_back (frequency);

		++iCycle;

		if (verbose >= 1)
		{
			double delta_f = C(3+2*nofHarmonics)/(2*M_PIl);
			fprintf (stdout, " %02d\t%14.14g\t%+21.21g\t%+10.10g\t%+10.10g\t%+10.10g\n",\
					iCycle, est.f_iter[iCycle], delta_f, est.A_iter[iCycle], est.B_iter[iCycle], est.C_iter[iCycle]);
		}
		
	} while ( ( iCycle < maxCycle ) && ( std::abs( C(3+2*nofHarmonics) ) > omegaError ) );

	fprintf (stdout, "\n");

}


#endif
