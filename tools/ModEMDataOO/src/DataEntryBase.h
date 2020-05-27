#include <string>
#include <vector>
#include <iostream>
//
#ifndef DATA_ENTRY_BASE_H_
#define DATA_ENTRY_BASE_H_
//
class DataEntryBase
{
public:
	// TYPE, PERIOD, CODE, X, Y, Z, COMPONENT, REAL, IMAGINARY, ERROR
	DataEntryBase( std::string, double, std::string, double, double, double, std::string, double, double, double );
	~DataEntryBase( void );
	//
	void show( void );
	int typeStringToInt( std::string );
	//
	// GETTERS & SETTERS
	//
	std::string getType( void ) const;
	void setType( std::string );
	//
	double getPeriod( void ) const;
	void setPeriod( double );
	//
	std::string getCode( void ) const;
	void setCode( std::string );
	//
	double getLat( void ) const;
	void setLat( double );
	//
	double getLong( void ) const;
	void setLong( double );
	//
	double getX( void ) const;
	void setX( double );
	//
	double getY( void ) const;
	void setY( double );
	//
	double getZ( void ) const;
	void setZ( double );
	//
	std::string getComp( void ) const;
	void setComp( std::string );
	//
	double getReal( void ) const;
	void setReal( double );
	//
	double getImag( void ) const;
	void setImag( double );
	//
	double getError( void ) const;
	void setError( double );

protected:
	std::string type;
	double period, latitude, longitude, x, y, z, real, imaginary, error;
	std::string code, component;
};

#endif

std::string getValidStringType( std::string );
