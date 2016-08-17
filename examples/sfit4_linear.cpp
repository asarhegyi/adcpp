/************************************************************
* sfit3_linear -- program to do a 4 parameter sine-wave fit *
*        using the gsl_multifit_linear function. Note that  *
*        this is an iterative algorithm. Initial estimate of*
*        the coefficients are needed.                       *
*                                                           *
* Author: Attila Sarhegyi                                   *
*                                                           *
* Version: 1.0                                              *
*                                                           *
* Purpose: Experiment various Least-Squares algorithms and  *
*        compare it with MATLAB                             *
*                                                           *
* Usage:                                                    *
*        sfit4_linear -help                                 *
*                                                           *
* Example:                                                  *
*        sfit4_linear -f1000.0445712 -s500000 -r0 \         *
*        ../charData/test.int                               *
*                                                           *
************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctime>
#include "../src/sfit_class.h"

int verbose = 0;               // verbose mode (default = false)
int debug = 0;                 // debug mode (default = false)
std::string mode = "user";     // used to distinguish user and test time mode
char *program_name;            // name of the program (for errors)
int nofHarmonics;              // number of sine-wave harmonics
double frequency;              // the fundamental frequency of the sine-wave
double fSample;                // sampling frequency
std::string inputFile;         // input file name


/************************************************************
* usage -- tell the user how to use this program and exit   *
************************************************************/
void usage()
{
    std::cerr << "Usage is " << program_name << " [options] [file-name]\n";
    std::cerr << "Options\n";
    std::cerr << "  -v          Verbose mode. Overwrites verbosity level if used with debug mode!\n";
    std::cerr << "  -d          Debug mode. It also set the verbosity level to max\n";
    std::cerr << "  -t          Test time mode for the Monte Carlo simulations\n";
    std::cerr << "  -f<number>  Fundamental frequency of the sine-wave\n";
    std::cerr << "  -s<number>  Sampling frequency\n";
    std::cerr << "  -r<number>  Number of harmonics\n";
    std::cerr << "  -h[elp]     or\n";
    std::cerr << "  -u[sage]    Print this help message\n";
    std::cerr << "Example: ./sfit4_linear -f1000.0445712 -s500000 -r2 ../charData/test.int\n";
    exit(8);
}


int main (int argc, char *argv[])
{


    /************************************************************
    *                                                           *
    * Parse command line options                                *
    *                                                           *
    ************************************************************/

    // Save the program name for future use
    program_name = argv[0];

    // Loop for each option. Stop if we run out of arguments or
    // we get an argument without a dash.
    while ((argc > 1) && (argv[1][0] == '-'))
       {
        // argv[1][1] is the actual option character
        switch (argv[1][1])
           {
            // -v verbose 
            case 'v':
                verbose = atoi(&argv[1][2]);
                break;
            // -d debug: set maximum verbosity level
            case 'd':
                debug = 1;
                verbose = 100;
                break;
            // -t test time mode for the Monte Carlo simulations
            case 't':
                mode = "test_time";
                break;
            // -f<number> set fundamental frequency
            case 'f':
                frequency = strtod(&argv[1][2], NULL);
                break;
            // -s<number> set sampling frequency
            case 's':
                fSample = strtod(&argv[1][2],NULL);
                break;
            // -r<number> set number of harmonics
            case 'r':
                nofHarmonics = atoi(&argv[1][2]);
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

    /************************************************************
    * At this point all the options have been processed. Check  *
    * to see if we have no files in the list and if so, we need *
    * to list the usage and exit the program.                   *
    ************************************************************/
    if (argc == 1)
    {
        std::cerr << "There is no input file-name specified" << std::endl;
        usage();
        exit(1);
    }
    else
        inputFile = argv[1];

    if (verbose >= 1)
    {
        std::cout << "Command line arguments parsed:" << std::endl;
        std::cout << "\tVerbose: " << verbose << std::endl;
        std::cout << "\tMode: " << mode << std::endl;
        std::cout << "\tFundamental frequency: " << frequency << std::endl;
        std::cout << "\tSampling frequency: " << fSample << std::endl;
        std::cout << "\tNumber of harmonics: " << nofHarmonics << std::endl;
        std::cout << "\tInput file: " << inputFile << std::endl;
        std::cout << std::endl;
    }


    class sfit mysfit;
    mysfit.verbose = verbose;
    mysfit.mode = mode;

    clock_t start = clock();   // Start timer

    mysfit.sfit4_linear(inputFile, frequency, nofHarmonics, fSample, "FFT");
    mysfit.ABC2AmpPhaseDC();
    mysfit.residual_linear();

    clock_t end = clock();     // Stop timer
    double dif = difftime(end, start)/CLOCKS_PER_SEC;

    mysfit.print_coeffs();
    mysfit.print_parameters();

    fprintf (stdout, "fs: %.20f\n", mysfit.fs);
    fprintf (stdout, "Execution time: %.6f seconds\n", dif);
    fprintf (stdout, "\n");

    mysfit.fitted_linear();
    mysfit.print_chisqRsq();

    return 0;
}
