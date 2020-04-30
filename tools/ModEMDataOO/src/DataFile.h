#include <string>
#include <vector>
#include <iostream>
//
#include "DataEntryMT.h"
#include "DataEntryCSEM.h"
#include "StationEntries.h"
//
#define TOL6 0.000001
//
#ifndef DATA_FILE_H_
#define DATA_FILE_H_
//
class DataFile {
public:
	DataFile( std::string );
	~DataFile( void );
	//
	void loadLines( void );
	void loadEntries( void );
	bool existType( std::string );
	int getStationIndexByPeriod( std::string, double );
	void loadStations( void );
	//
	// GETTERS & SETTERS
	std::string getPath( void ) const;
	void setPath( std::string );
	//
	std::vector< std::string > getLines( void );
private:
	std::string path, actual_type;
	int header_counter;
	std::vector< std::string > 		_types;
	std::vector< std::string > 		_lines;
	std::vector< DataEntryBase* > 	_data_entries;
	std::vector< StationEntries* > 	_station_entries;

};
//
#endif
//
// SEPARA OS ESPACOS
std :: vector< std::string > split( std::string );
std :: vector< std::string > split( std::string, char );
