#include <string>
#include <vector>
#include <iostream>
//
#define TOL6 0.000001
//
#ifndef DATA_BASE_H_
#define DATA_BASE_H_
//
class DataBase {
public:
	DataBase( std::string );
	~DataBase( void );
	//
	void loadLines( void );
	void loadInfo( void );
	bool existType( std::string );
	//
	// GETTERS & SETTERS
	std::string getPath( void ) const;
	void setPath( std::string );
	//
	std::vector< std::string > getLines( void );
private:
	std::string path;
	int header_counter;
	std::vector< std::string > 		_types;
	std::vector< std::string > 		_lines;
};
//
#endif
//
// SEPARA OS ESPACOS
std :: vector< std::string > split( std::string );
std :: vector< std::string > split( std::string, char );
std::string getValidStringType( std::string );
