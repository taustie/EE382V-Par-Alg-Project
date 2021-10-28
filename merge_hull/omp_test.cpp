#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>

int main( int argc, char* argv[] )
{
	omp_set_num_threads( 32 );	

	double pi = acos(-1.0);

	int size = 1024 / sizeof(double); // 1 KB of doubles 
	size *= 1024; // 1 MB
	size *= 1024; // 1 GB 
	size *= 4; // 4 GB 
	
	std::cout << "Allocating " << size * sizeof(double) / (1024*1024) << " MB of memory ..." << std::endl;
	std::vector<double> my_vector( size , 0.0 ); // 1280000000
	std::cout << "Done!" << std::endl << std::endl;
	
	std::cout << "Entering main loop ..." << std::endl;
	
	#pragma omp parallel for
	for( int i=0; i < my_vector.size(); i++ )
	{
		my_vector[i] = exp( -sin( i*i + pi*log(i+1.0) ) );
		my_vector[i] += -exp( -cos( i*i*i + pi*log10(i+10.0) ) );
		my_vector[i] -= exp( -tan( i*i*i*i + pi*log(i+99.0) ) );
		my_vector[i] += exp( -sin( i*i*i*i*i + pi*log(i+2.0) ) );
		my_vector[i] -= -exp( -cos( i*i*i*i*i*i + pi*log10(i+20.0) ) );
		my_vector[i] += exp( -tan( i*i*i*i*i*i*i + pi*log(i+199.0) ) );		
		
		// add this if your test is too fast!! :-)
	/*
		my_vector[i] = 0.0; 
		my_vector[i] += -exp( -cos( i*i*i + pi*log10(i+40.0) ) );
		my_vector[i] -= exp( -tan( i*i*i*i + pi*log(i+19.0) ) );
		my_vector[i] += exp( -sin( i*i*i*i*i + pi*log(i+3.0) ) );
		my_vector[i] -= -exp( -cos( i*i*i*i*i*i + pi*log10(i+17.0) ) );
		my_vector[i] += exp( -tan( i*i*i*i*i*i*i + pi*log(i+499.0) ) );		
	*/		
	}
	std::cout << "Done!" << std::endl;
	
	return 0;
}
