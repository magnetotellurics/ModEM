#include <stddef.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "StationEntries.h"
//
StationEntries :: StationEntries( int polarization_counter )
{
	this->polarization_counter = polarization_counter;
}
StationEntries :: ~StationEntries( void ){}
//
// GETTERS & SETTERS
int StationEntries :: getPolarizationCounter( void ) const
{
	return polarization_counter;
}
void StationEntries :: setPolarizationCounter( int polarization_counter )
{
	this->polarization_counter = polarization_counter;
}
//
std::vector< DataEntryBase* > StationEntries :: getDataEntries( void ) const
{
	return _data_entries;
}
void StationEntries :: setDataEntries( std::vector< DataEntryBase* > _data_entries )
{
	this->_data_entries = _data_entries;
}
void StationEntries :: addDataEntry( DataEntryBase* data_entry )
{
	_data_entries.push_back( data_entry );
}
