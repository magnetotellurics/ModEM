#include <string>
#include <vector>
#include <iostream>
//
//
#ifndef TEXT_FILE_H_
#define TEXT_FILE_H_
//
class TextFile {
public:
	TextFile( std::string );
	~TextFile( void );
	//
	void load( void );
	//
	// GETTERS & SETTERS
	size_t getSize( void ) const;
	std::vector< std::string > getLines( void );
	std::string getPath( void ) const;
	void setPath( std::string );

protected:
	std::string path, name;
	std::vector< std::string > _lines;
};

#endif

//
// SEPARA OS ESPACOS
std::string trim( std::string );
std :: vector< std::string > split( std::string );
std :: vector< std::string > split( std::string, char );
