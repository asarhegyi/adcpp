/************************************************************
* fileio -- unit test the file interface:					*
*		read data from file, load it into a vector and save *
*		off the content of the vector into a different file.*
*		Compare the input and output file manually and		*
*		verify behaviour.									*
*															*
* Author: Attila Sarhegyi									*
*															*
* Version: 1.0												*
*															*
* Purpose: Standardize the read and write of a data file	*
*															*
* Usage:													*
*		make -f ../Makefile TARGET=fileio					*
*		./fileio											*
************************************************************/
#include <iostream>		// std::cin, std::cout
#include <fstream>		// std::ifstream, std::ofstream
#include <vector>		// std::vector
#include "../src/utils.h"


int verbose;
std::string inputFile;
double fs;
int N;


int main (void)
{
	std::string test = "stream";
	std::vector<double> time;
	std::vector<double> data;
	fs = 36000;

	//inputFile = "../charData/test.int";
	inputFile = "../charData/dacout_sine_short.txt";


	if ( test.compare("load") == 0 )
	{
		// read in two columns which are time and data or one
		// column which is data and generate the time vector from index/fs
		// use member functions from fileio class
		class fileio myfileio;
		myfileio.verbose = verbose;
		N = myfileio.load_data (inputFile, time, data, fs);

		dbl_array2file("check_time.dat", time, N, "");
		dbl_array2file("check_data.dat", data, N, "");
	}
	else if ( test.compare("stream") == 0 )
	{
		int nofSamples = 0;
		int blockSize;
		bool openFlag = false;
		bool eofFlag = false;
		std::ifstream inFile;
		std::string fmode = "";

		blockSize = 5000;

		// use member functions from fileio class
		class fileio sfit_fileio;
		sfit_fileio.verbose = verbose;

		while (true)
		{
			// Load blockSize chunks of data intermittently from the input file 
			N = sfit_fileio.stream_data (openFlag, eofFlag, inFile, inputFile, time, data, fs, blockSize, nofSamples);
			nofSamples += N;
			
			// Filter out corner case when blockSize == number of lines in the input file
			if ( N == 0 )
				break;

			fprintf (stdout, "Done with buffering a %d line block.\n", N);

			dbl_array2file("check_time.dat", time, N, fmode);
			dbl_array2file("check_data.dat", data, N, fmode);

			fmode = "Append";

			// Exit while loop after the last block is processed
			if (eofFlag)
				break;
		}
	}

	return(0);
}

