#include <stddef.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "TextFile.h"
//
TextFile :: TextFile( std::string path )
{
	this->path = path;
	//
	std::vector< std::string > _path_parts = split( path, '/' );
	this->name = _path_parts[ _path_parts.size() - 1 ];
	//
	load();
}
//
TextFile :: ~TextFile( void ){}
//
void TextFile :: load( void )
{
	// READ FILE LINES TO VECTOR _LINES
	std::ifstream in( path.c_str() );
	//
	if ( in.is_open() )
	{
		std::string file_content( ( std :: istreambuf_iterator<char>( in ) ), std :: istreambuf_iterator<char>() );
		std::vector< std::string > _file_content_lines = split( file_content, '\n' );
		//
		size_t lid = -1;
		//
		while( ++lid < _file_content_lines.size() )
			if( trim( _file_content_lines[ lid ] )[ 0 ] != '#' &&
					trim( _file_content_lines[ lid ] )[ 0 ] != '>' )
			{
				_lines.push_back( _file_content_lines[ lid ] );
			}
	}
	else
		std::cout << "Can't open [" << path << "]" << std::endl;
}
//
// GETTERS & SETTERS
std::vector< std::string > TextFile :: getLines( void )
{
	return _lines;
}
//
size_t TextFile :: getSize( void ) const
{
	return _lines.size();
}
//
//
std::string trim( std::string target )
{
	std::string trimed = target;
	//
	while( isspace( trimed[0] ) || trimed[0] == ' ' || trimed[0] == '\t' )
		trimed.erase( trimed.begin() );
	//
	while( isspace( trimed[ trimed.size() - 1 ] ) || trimed[ trimed.size() - 1 ] == ' ' || trimed[ trimed.size() - 1 ] == '\t' )
		trimed.erase( trimed.begin() + trimed.size() - 1  );
	//
	return trimed;
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


