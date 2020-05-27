#include "DataFile.h"

#include <stddef.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

DataFile :: DataFile( std::string path )
{
	this->path = path;
	//
	std::cout << "Loading [" << path << "]:" << std::endl;
	//
	loadLines();
	//
	std::cout << "Defining entries:" << std::endl;
	//
	loadEntries();
	//
	std::cout << "Grouping stations by Type and Period:" << std::endl;
	//
	loadStations();
}
DataFile :: ~DataFile( void ){}
//
void DataFile :: loadLines( void )
{
	// READ FILE LINES TO VECTOR _LINES
	std::ifstream in( path.c_str() );
	//
	if ( in.is_open() )
	{
		std::string file_content( ( std :: istreambuf_iterator<char>( in ) ), std :: istreambuf_iterator<char>() );
		_lines = split( file_content, '\n' );
		//
		std::cout << "\t" << _lines.size() << " lines!" << std::endl;
	}
	else
		std::cout << "\tCan't open [" << path << "]" << std::endl;
}
bool DataFile :: existType( std::string type )
{
	for( size_t tid = 0; tid < _types.size(); tid++ )
		if( _types[ tid ] == type )
			return true;
	return false;
}
void DataFile :: loadEntries( void )
{
	int header_line_counter = 0;
	//
	for( size_t lid = 0; lid < _lines.size(); lid++ )
	{
		std::vector< std::string > _line_parts = split( _lines[ lid ] );
		//
		if( _line_parts[ 0 ][ 0 ] != '#' &&
				_line_parts[ 0 ][ 0 ] != '>' )
		{
			if( header_line_counter )
			{
				header_line_counter = 0;
				if( !existType( actual_type ) )
					_types.push_back( actual_type );
			}
			//
			if( actual_type == "Full_Impedance" || actual_type == "Full_Vertical_Components" )
			{
					_data_entries.push_back( new DataEntryMT( actual_type,
							atof( _line_parts[ 0 ].c_str() ),
							_line_parts[ 1 ],
							atof( _line_parts[ 2 ].c_str() ),
							atof( _line_parts[ 3 ].c_str() ),
							atof( _line_parts[ 4 ].c_str() ),
							atof( _line_parts[ 5 ].c_str() ),
							atof( _line_parts[ 6 ].c_str() ),
							_line_parts[ 7 ],
							atof( _line_parts[ 8 ].c_str() ),
							atof( _line_parts[ 9 ].c_str() ),
							atof( _line_parts[ 10 ].c_str() ) ) );
			}
			else if( actual_type == "Ex_Field" || actual_type == "Ey_Field" || actual_type == "Bx_Field" ||
					actual_type == "By_Field" || actual_type == "Bz_Field" )
			{
				//# Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(x) Tx_Z(m) Code X(m) Y(x) Z(m) Component Real Imag, Error
				std::string dipole 		= _line_parts[ 0 ];
				double period 			= atof( _line_parts[ 1 ].c_str() );
				double moment 			= atof( _line_parts[ 2 ].c_str() );
				double azimuth 			= atof( _line_parts[ 3 ].c_str() );
				double dip 			= atof( _line_parts[ 4 ].c_str() );
				double tx_x 			= atof( _line_parts[ 5 ].c_str() );
				double tx_y 			= atof( _line_parts[ 6 ].c_str() );
				double tx_z 			= atof( _line_parts[ 7 ].c_str() );
				std::string code 		= _line_parts[ 8 ];
				double x 				= atof( _line_parts[ 9 ].c_str() );
				double y 				= atof( _line_parts[ 10 ].c_str() );
				double z 				= atof( _line_parts[ 11 ].c_str() );
				std::string component 	= _line_parts[ 12 ];
				double real 			= atof( _line_parts[ 13 ].c_str() );
				double imaginary 		= atof( _line_parts[ 14 ].c_str() );
				double error 			= atof( _line_parts[ 15 ].c_str() );
				//
				_data_entries.push_back( new DataEntryCSEM( actual_type, dipole, period, moment, azimuth, dip,
						tx_x, tx_y, tx_z, code, x, y, z, component, real, imaginary, error ) );
			}
			else
				std::cout << "Not implemented for type: [" << actual_type << "]" << std::endl;
		}
		else
		{
			header_line_counter++;
			if( header_line_counter == 3 ) // TYPE
			{
				actual_type = getValidStringType( _lines[ lid ] );
			}
		}
	}
	//
	std::cout << "\t" << _data_entries.size() << " entries!" << std::endl;
	std::cout << "\t" << _types.size() << " header types:" << std::endl;
	for( size_t tid = 0; tid < _types.size(); tid++ )
		std::cout << "\t\t" << _types[ tid ] << std::endl;
}
int DataFile :: getStationIndexByPeriod( std::string type, double period )
{
	if( _station_entries.size() )
	{
		for( size_t seid = 0; seid < _station_entries.size(); seid++ )
		{
			for( size_t deid = 0; deid < _station_entries[ seid ]->getDataEntries().size(); deid++ )
			{
				if( _station_entries[ seid ]->getDataEntries()[ deid ]->getType() == type &&
						fabs( _station_entries[ seid ]->getDataEntries()[ deid ]->getPeriod() - period ) < TOL6 )
				{
					return seid;
				}
			}
		}
		return -1;
	}
	else
		return -1;
}
void DataFile :: loadStations( void )
{
	for( size_t deid = 0; deid < _data_entries.size(); deid++ )
	{
		int founded_station_index = getStationIndexByPeriod( _data_entries[ deid ]->getType(), _data_entries[ deid ]->getPeriod() );
		if( founded_station_index == -1 )
		{
			StationEntries* station_entries = new StationEntries( 1 );
			station_entries->addDataEntry( _data_entries[ deid ] );
			_station_entries.push_back( station_entries );
		}
		else
			_station_entries[ founded_station_index ]->addDataEntry( _data_entries[ deid ] );
	}
	//
	std::cout << "\t" << _station_entries.size() << " stations!" << std::endl;
	//
	std::cout << "\t\tSize\tType\tPeriod:" << std::endl;
	for( size_t seid = 0; seid < _station_entries.size(); seid++ )
		std::cout << "\t\t(" << _station_entries[ seid ]->getDataEntries().size() << ")\t" << _station_entries[ seid ]->getDataEntries()[0]->getType() << "\t" << _station_entries[ seid ]->getDataEntries()[0]->getPeriod() << std::endl;
}
//
//
std::vector< std::string > split( std::string target )
{
	std::stringstream ss( target );
	std::istream_iterator< std::string > begin( ss ), end;
	std::vector< std::string > vstrings( begin, end );
	//
	return vstrings;
}
std::vector< std::string > split( std::string target, char del )
{
	std::vector< std::string > internal;
	std::stringstream ss( target );
	std::string tok;
	//
	while( getline( ss, tok, del ) )
		internal.push_back( tok );
	//
	return internal;
}
