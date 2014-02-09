#include <stdio.h>

void plotdata(char fname[])
{
	//char fname[] = "fourtest1.txt";
	FILE *pipe = popen("gnuplot -persist","w");
//	fprintf(pipe, "set data s lines\n");
	fprintf(pipe,"plot '%s' using 1:2 title '%s' with lines\n",fname,fname);
	close(pipe);


}


