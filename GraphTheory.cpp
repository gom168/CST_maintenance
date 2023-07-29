#include "MSTruss.h"
#include "core.h"
#define DIRNAME "graphviz"

using namespace std;



int main()
{
	core *graph = new core("4.txt",16); 
	
	graph->read_graph(false,0); // read graph,  0 for number begin from 1, 0 for number begin from 1 
	graph->Dynamic_maintainences("update.txt", 1, 0, 0); // files; 0 for dynamic 1 for staric; 0 for delete 1 for insert; 0 for number begin form 1, 1 for number begin from 0
	graph->Dynamic_kcores("update.txt", 0, 1, 0, 2, 30); // files; 0 for dynamic 1 for staric; 0 for delete 1 for insert; 0 for number begin form 1, 1 for number begin from 0
}


