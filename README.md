# Construction-of-CSA-Compressed-Suffix-Array-
A repository about construction of compressed suffix array (bio-info technology)

Usage manual:
	modify main.c to change input file path, output file path, output file header. 

details: 
	/* ##############################

	char* FILEPATH = "testdata/testdata_1000.fna";   // file path
	int ARRAYLENGTH = 0; // length of T ~ n
	int PARTLENGTH = 0; // part length of T ~ l
	int PARTNUM = 0; // number of parts ~ ceil(n/l)

	char* T = NULL; // DNA sequence of (A,C,G,T) plus a '$'
	int* SA = NULL; // SA of T
	int* SA_inverse = NULL; // inverse of SA
	int* Psi = NULL; // Psi of T - the compressed suffix array

	char* BWT = NULL; // BWT of T - Burrows-Wheeler Transform

	char* BWTFILEPATH = "outputdata/testdata_1000.bwt";
	char* BWTFILEHEADER = ">gi|110640213|ref|NC_008253.1| Escherichia coli 536, bwt array";
	int LINELENGTH = 70;

	############################## */

FILEPATH ~ input file path
BWTFILEPATH ~ output file path
BWTFILEHEADER ~ output file header

TODO
	maybe will add function using format like "./xxx -arg ... ... ..."
