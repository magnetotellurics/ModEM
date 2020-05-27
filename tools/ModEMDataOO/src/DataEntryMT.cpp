#include "DataEntryMT.h"

#include <stddef.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
//
// Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error
DataEntryMT :: DataEntryMT( std::string type, double period, std::string code, double latitude,
		double longitude, double x, double y, double z, std::string component,
		double real, double imaginary, double error )
:DataEntryBase( type, period, code, x, y, z, component, real, imaginary, error )
{
	this->latitude = latitude;
	this->longitude = longitude;
}
DataEntryMT :: ~DataEntryMT( void ){}
//
void DataEntryMT :: show( void )
{
	std:: cout << "MT DATA ENTRY:" <<
	"\n\tperiod:\t" << period <<
	"\n\tcode:\t" << code <<
	"\n\tlatitude:\t" << latitude <<
	"\n\tlongitude:\t" << longitude <<
	"\n\t(x,y,z):\t(" << x << ", " << y << ", " << z << ")" <<
	"\n\tcomponent:\t" << component <<
	"\n\treal:\t" << real <<
	"\n\timaginary:\t" << imaginary <<
	"\n\terror:\t" << error  << std::endl;
}
//
// GETTERS & SETTERS
//
double DataEntryMT :: getLat( void ) const
{
	return latitude;
}
void DataEntryMT :: setLat( double latitude )
{
	this->latitude = latitude;
}
//
double DataEntryMT :: getLong( void ) const
{
	return longitude;
}
void DataEntryMT :: setLong( double longitude )
{
	this->longitude = longitude;
}
