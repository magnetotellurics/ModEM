#include <string>
#include <vector>
#include <iostream>
//
#include "Period.h"
#include "TextFile.h"
//
#ifndef DATA_FILE_H_
#define DATA_FILE_H_
//
class DataFile : public TextFile {
public:
	DataFile( std::string, std::string );
	~DataFile( void );
	//
	void loadValues( void );
	void createHtml( void );
	int existPeriodByIndex( int );
	int existPeriodByValue( double );
	//
	std::string getMathBoxBundleContent( void );
private:
	std::string math_box_path;
	std::vector< Period* > _periods;
	//
	double lx, ly, lz;
	double hx, hy, hz;
	//
	int max_polarization_index;
};
//
#endif
