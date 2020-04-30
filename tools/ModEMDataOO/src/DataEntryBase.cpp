#include <stddef.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "DataEntryBase.h"
//
//Period(s) Code X(m) Y(m) Z(m) Component Real Imag Error
DataEntryBase :: DataEntryBase( std::string type, double period, std::string code,
		double x, double y, double z, std::string component,
		double real, double imaginary, double error )
{
	this->type = type;
	this->period = period;
	this->code = code;
	this->latitude = latitude;
	this->longitude = longitude;
	this->x = x;
	this->y = y;
	this->z = z;
	this->component = component;
	this->real = real;
	this->imaginary = imaginary;
	this->error = error;

}
DataEntryBase :: ~DataEntryBase( void ){}
//
void DataEntryBase :: show( void )
{
	std:: cout << "BASE DATA ENTRY:" <<
	"\n\tperiod:\t" << period <<
	"\ncode:\t" << code <<
	"\nlatitude:\t" << latitude <<
	"\nlongitude:\t" << longitude <<
	"\n(x,y,z):\t(" << x << ", " << y << ", " << z << ")" <<
	"\ncomponent:\t" << component <<
	"\nreal:\t" << real <<
	"\nimaginary:\t" << imaginary <<
	"\nerror:\t" << error  << std::endl;
}
//
// GETTERS & SETTERS
//
std::string DataEntryBase :: getType( void ) const
{
	return type;
}
void DataEntryBase :: setType( std::string type )
{
	this->type = type;
}
//
double DataEntryBase :: getPeriod( void ) const
{
	return period;
}
void DataEntryBase :: setPeriod( double period )
{
	this->period =period;
}
//
std::string DataEntryBase :: getCode( void ) const
{
	return code;
}
void DataEntryBase :: setCode( std::string code )
{
	this->code = code;
}
//
double DataEntryBase :: getLat( void ) const
{
	return latitude;
}
void DataEntryBase :: setLat( double latitude )
{
	this->latitude = latitude;
}
//
double DataEntryBase :: getLong( void ) const
{
	return longitude;
}
void DataEntryBase :: setLong( double longitude )
{
	this->longitude = longitude;
}
//
double DataEntryBase :: getX( void ) const
{
	return x;
}
void DataEntryBase :: setX( double x )
{
	this->x = x;
}
//
double DataEntryBase :: getY( void ) const
{
	return y;
}
void DataEntryBase :: setY( double y )
{
	this->y = y;
}
//
double DataEntryBase :: getZ( void ) const
{
	return z;
}
void DataEntryBase :: setZ( double z )
{
	this->z = z;
}
//
std::string DataEntryBase :: getComp( void ) const
{
	return component;
}
void DataEntryBase :: setComp( std::string component )
{
	this->component = component;
}
//
double DataEntryBase :: getReal( void ) const
{
	return real;
}
void DataEntryBase :: setReal( double real )
{
	this->real =real;
}
//
double DataEntryBase :: getImag( void ) const
{
	return imaginary;
}
void DataEntryBase :: setImag( double imaginary )
{
	this->imaginary =imaginary;
}
//
double DataEntryBase :: getError( void ) const
{
	return error;
}
void DataEntryBase :: setError( double error )
{
	this->error = error;
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
