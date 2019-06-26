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
	load();
}
DataFile :: ~DataFile( void )
{
	for ( size_t i = 0; i < _lines.size(); i++ )
		delete [] __values[ i ];
	delete [] __values;
}
void DataFile :: load( void )
{
	// READ FILE LINES TO VECTOR _LINES
	std::ifstream in( path.c_str() );
	std::string file_content( ( std :: istreambuf_iterator<char>( in ) ), std :: istreambuf_iterator<char>() );
	std::vector< std::string > _file_content_lines = split( file_content, '\n' );
	//
	size_t lid = 0;
	//
	while( ++lid < _file_content_lines.size() )
		if( _file_content_lines[ lid ][ 0 ] != '#' &&
				_file_content_lines[ lid ][ 0 ] != '>' )
			_lines.push_back( _file_content_lines[ lid ] );
	//
	// BUILD VALUE MATRIX (N, 2)
	__values = new double*[ _lines.size() ];
	for ( size_t i = 0; i < _lines.size(); i++ )
		__values[ i ] = new double[ 2 ];
	//
	for( lid = 0; lid < _lines.size(); lid++ )
	{
		std::vector< std::string > _line_parts = split( _lines[ lid ] );
		__values[ lid ][ 0 ] = atof( _line_parts[ 8 ].c_str() );
		__values[ lid ][ 1 ] = atof( _line_parts[ 9 ].c_str() );
	}
}
double DataFile :: euclidianNorm( DataFile* other )
{
	double diff = 0.0;
	//
	for( int lid = 0; lid < ( ( getSize() < other->getSize() ) ? getSize() : other->getSize() ); lid++ )
		diff += sqrt( ( ( __values[ lid ][ 0 ] - other->getValues()[ lid ][ 0 ] ) * ( __values[ lid ][ 0 ] - other->getValues()[ lid ][ 0 ] ) ) + ( ( __values[ lid ][ 1 ] - other->getValues()[ lid ][ 1 ] ) * ( __values[ lid ][ 1 ] - other->getValues()[ lid ][ 1 ] ) ) );
	//
	return diff;
}
double DataFile :: norm1( DataFile* other )
{
	double diff = 0.0;
	//
	for( int lid = 0; lid < ( ( getSize() < other->getSize() ) ? getSize() : other->getSize() ); lid++ )
		diff += fabs( __values[ lid ][ 0 ] - other->getValues()[ lid ][ 0 ] ) + fabs( __values[ lid ][ 1 ] - other->getValues()[ lid ][ 1 ] );
	//
	return diff;
}
//
// GETTERS & SETTERS

double** DataFile :: getValues( void ) const
{
	return __values;
}
std::vector< std::string > DataFile :: getLines( void )
{
	return _lines;
}
size_t DataFile :: getSize( void ) const
{
	return _lines.size();
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


