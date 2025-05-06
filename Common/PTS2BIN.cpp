// Copyright (c) 2025, Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski
// of Imperial College London and Dr Marc Seidel of Siemens Gamesa Renewable
// Energy.

// Copyright under a BSD 3-Clause License, see
// https://github.com/KathirkamanathanL/AxisymmetricShellGDTProtocolPublic.git

// Last modified at 22.01 on 28/04/2025


// A simple C++ program to transform the .pts text-based file containing tower coordinates
// into a .bin file with same but in a much more efficient (for storage and reading) binary format

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		std::cout << "Usage: PTS2BIN pts_file_root_no_extension" << std::endl;
		std::cout << "e.g.: PTS2BIN test // (for test.pts file)" << std::endl << std::endl;
		return 1;
	}

	std::string sptsFileName = (std::string)argv[1];
	std::fstream ptsFile, BINFile;

	ptsFile.open(sptsFileName + ".pts", std::ios::in);
	BINFile.open(sptsFileName + ".bin", std::ios::out | std::ios::binary);
	if (ptsFile.is_open() && BINFile.is_open())
	{
		std::cout << "File " << sptsFileName + ".pts" << " opened." << std::endl;
		uint64_t unCounter = 0;
		BINFile.write((char*)&unCounter, sizeof(uint64_t));

		double x, y, z;
		std::string sLine;
		while (std::getline(ptsFile, sLine))
		{
			std::istringstream isLine(sLine);
			isLine >> x;
			isLine >> y;
			isLine >> z;
			BINFile.write((char*)&x, sizeof(double));
			BINFile.write((char*)&y, sizeof(double));
			BINFile.write((char*)&z, sizeof(double));
			unCounter += 1;
			if (!(unCounter%1000000)) std::cout << "Read " << unCounter/1e6 << " million pts rows." << std::endl;
		}
		ptsFile.close();
		BINFile.seekp(0);
		BINFile.write((char*)&unCounter, sizeof(uint64_t));
		BINFile.close();
		std::cout << "File " << sptsFileName + ".bin" << " created." << std::endl;
		std::cout << "Total " << unCounter << " data points." << std::endl << std::endl;
	}
	else
	{
		std::cout << "File " << sptsFileName + ".pts" << " not found." << std::endl;
	}

	return 0;
}