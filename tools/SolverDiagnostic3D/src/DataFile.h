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
	DataFile( std::string );
	~DataFile( void );
	//
	void loadValues( void );
	void createHtml( void );
	int existPeriod( double );
private:
	std::vector< Period* > _periods;
	//
	double lx, ly, lz;
	double hx, hy, hz;
};
//
#endif
