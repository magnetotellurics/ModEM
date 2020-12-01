#include <stddef.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "DataBase.h"

DataBase :: DataBase( std::string path )
{
	this->path = path;
	//
	loadLines();
	//
	loadInfo();
}
DataBase :: ~DataBase( void ){}
//
void DataBase :: loadLines( void )
{
	// READ FILE LINES TO VECTOR _LINES
	std::ifstream in( path.c_str() );
	//
	if ( in.is_open() )
	{
		std::string file_content( ( std :: istreambuf_iterator<char>( in ) ), std :: istreambuf_iterator<char>() );
		_lines = split( file_content, '\n' );
	}
	else
		std::cout << "\tCan't open [" << path << "]" << std::endl;
}
bool DataBase :: existType( std::string type )
{
	for( size_t tid = 0; tid < _types.size(); tid++ )
		if( _types[ tid ] == type )
			return true;
	return false;
}
void DataBase :: loadInfo( void )
{
	int header_line_counter = 0;
	std::string actual_type;
	int actual_periods, nprocs = 0;
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
			// MT
			if( actual_type == "Full_Impedance" || actual_type == "Full_Vertical_Components" )
			{
				if( actual_periods * 2 + 1 > nprocs )
					nprocs = actual_periods * 2 + 1;
			}
			// CSEM
			else if( actual_type == "Ex_Field" || actual_type == "Ey_Field" || actual_type == "Bx_Field" ||
					actual_type == "By_Field" || actual_type == "Bz_Field" )
			{
				if( actual_periods + 1 > nprocs )
					nprocs = actual_periods + 1;
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
			if( header_line_counter == 8 ) // TYPE
			{
				std::vector< std::string > _line_parts = split( _lines[ lid ] );
				actual_periods = atoi( _line_parts[ 1 ].c_str() );
			}
		}
	}
	//
	int ppn = nprocs / 100 ? nprocs / 100 + 7 : 7;
	//
	std::stringstream prm_output;
	prm_output << "nodes=1:ppn=1+";
	nprocs--;
	//
	if( nprocs <= ppn )
	{
		prm_output << "1:ppn=" << nprocs;
	}
	else
	{
		int nodes = 0;
		//
		while( ( nodes + 1 ) * ppn  <= nprocs )
			nodes++;
		//
		prm_output << nodes << ":ppn=" << ppn;
		//
		if( ( nprocs - nodes * ppn ) > 0 )
		{
			//
			prm_output << "+1:ppn=" << nprocs - nodes * ppn;
			//
		}
	}
	//
	std::cout << prm_output.str() << std::endl;
}
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
//
//
std::string getValidStringType( std::string type )
{
	if( type.find( "Full_Impedance" ) != std::string::npos )
		return "Full_Impedance";
	else if( type.find( "Full_Vertical_Components" ) != std::string::npos )
		return "Full_Vertical_Components";
	else if( type.find( "Ex_Field" ) != std::string::npos )
		return "Ex_Field";
	else if( type.find( "Ey_Field" ) != std::string::npos )
		return "Ey_Field";
	else if( type.find( "Bx_Field" ) != std::string::npos )
		return "Bx_Field";
	else if( type.find( "By_Field" ) != std::string::npos )
		return "By_Field";
	else if( type.find( "Bz_Field" ) != std::string::npos )
		return "Bz_Field";
	else
	{
		std::cout << "Invalid type [" << type << "]" << std::endl;
		return "Invalid_Type";
	}
}
