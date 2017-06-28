#include <stdio.h>
#include <math.h>
int main(void)
{
    FILE *fp;
    char *name = "y0.dat";
    double dp, dy, dyp, dtof;
    fp = fopen(name, "r");
    while(fscanf(fp, "  %lf %lf %lf %lf\n", &dp, &dy, &dyp, &dtof)
	  !=EOF){
	if(dtof<-8.e-5){
	    if(dtof>-2.e-5){
		printf(" %12.9lf %12.9lf %12.9lf %12.9lf\n",
		       dp, dy, dyp, dtof);
	    }
	}
    }
    return 0;
}
