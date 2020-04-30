#include <string>
#include <vector>
#include <iostream>
//
#include "DataEntryMT.h"
#include "DataEntryCSEM.h"
//
#ifndef STATION_ENTRIES_H_
#define STATION_ENTRIES_H_
//
class StationEntries {
public:
	StationEntries( int );
	~StationEntries( void );
	//
	// GETTERS & SETTERS
	int getPolarizationCounter( void ) const;
	void setPolarizationCounter( int );
	//
	std::vector< DataEntryBase* > getDataEntries( void ) const;
	void setDataEntries( std::vector< DataEntryBase* > );
	void addDataEntry( DataEntryBase* );

protected:
	int polarization_counter;
	std::vector< DataEntryBase* > _data_entries;
};

#endif
