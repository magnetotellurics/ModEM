#include "DataFile.h"
//
//
int main( int argc, char ** argv )
{
	if( argc == 2 )
		DataFile* data_file_1 = new DataFile( argv[ 1 ] );
	else
		std::cout << "Please input a file" << std::endl;
}
