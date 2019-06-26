#include <string>
#include <vector>
#include <iostream>

//
//
#ifndef DATA_FILE_H_
#define DATA_FILE_H_
//
class DataFile {
public:
	DataFile( std::string );
	~DataFile( void );
	//
	void load( void );
	double norm1( DataFile* );
	double euclidianNorm( DataFile* );
	//
	// GETTERS & SETTERS
	size_t getSize( void ) const;
	std::vector< std::string > getLines( void );
	double** getValues( void ) const;
	std::string getPath( void ) const;
	void setPath( std::string );

private:
	std::string path;
	std::vector< std::string > _lines;
	double ** __values;
};

#endif

//
// SEPARA OS ESPACOS
std :: vector< std::string > split( std::string );
std :: vector< std::string > split( std::string, char );
