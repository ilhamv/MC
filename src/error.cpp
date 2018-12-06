#include "error.h"

void error(const std::string message){
	std::cout<<"[ERROR] "+message+"\n\n\n";
	std::exit(EXIT_FAILURE);
}
