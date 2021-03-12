#include <iostream>

#include <stdlib.h>

#include <boost/graph/adjacency_list.hpp>

using namespace boost;

int main(int argc, char* argv[])

{

	adjacency_list<> mygraph;

	add_edge(1, 2, mygraph);

	add_edge(1, 3, mygraph);

	add_edge(1, 4, mygraph);

	add_edge(2, 4, mygraph);

	return 0;

}