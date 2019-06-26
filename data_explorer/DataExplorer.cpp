#include "DataFile.h"
//
//
int main( int argc, char ** argv )
{
	if( argc == 3 )
	{
		DataFile * data_file_1 = new DataFile( argv[ 1 ] );
		//
		DataFile * data_file_2 = new DataFile( argv[ 2 ] );
		//
		if( data_file_1->getLines().size() == 0 ||
			data_file_2->getLines().size() == 0 ||
				data_file_1->euclidianNorm( data_file_2 ) )
		{
			std::cout << argv[ 1 ] << "\t" << argv[ 2 ] << "\t" << data_file_1->euclidianNorm( data_file_2 ) << "\tFAILED" << std::endl;
			return 1;
		}
		else
		{
			std::cout << argv[ 1 ] << "\t" << argv[ 2 ] << "\t" << data_file_1->euclidianNorm( data_file_2 ) << "\tSUCESS" << std::endl;
			return 0;
		}
		return 0;
	}
	else
	{
		std::cout << "MODEL\tDATA_TYPE\tERROR\tSUCESS" << std::endl;
		return 13;
	}
}
