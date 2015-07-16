#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

#define NUM_ELEMENTS 10000
#define MAX 1000.0

struct bin{
	double *B;
	int place;
	long int count,index;
};

double rand_normal(double mean, double stddev);
int isSorted(double a[],long int size);
void bucket_sort(double *A,long int num_elem,int buckets);
double get_wall_time();
static int compare (const void * a, const void * b);
static int compare_desc (const void * a, const void * b);

int main(int argc, char *argv[]){

	long int num_elem,i;
	int buckets,x;
	double *A;
	double start_time,end_time;

	num_elem = NUM_ELEMENTS;

	if (argc == 2){
		// First argument is the size of the list
		num_elem = atol(argv[1]);
	}else{
		printf("\nNo arguments given, continuing with default value: NUM_ELEMENTS 10000\n");
	}

	A = malloc(num_elem * sizeof(double));

	// Choose one seed - either srand48 for drand48 or srand for rand_normal
	//srand48((unsigned int)time(NULL));
	srand(time(NULL));

	// Initialize the list with random values - drand48 for uniform distribution and rand_normal for normal distribution
	for (i=0;i<num_elem;i++){
		//A[i] = drand48() * MAX;
		A[i] = rand_normal(MAX / 2.0, MAX / 5.0);
	}

	// By default openmp assigns as many threads as the available cores
	x = omp_get_max_threads();
	printf("\nNumber of threads: %d\n", x);

	buckets = x;

	start_time = get_wall_time();
	bucket_sort(A,num_elem,buckets);
	end_time = get_wall_time();

	// Check to see if the list was sorted correctly
	if (!isSorted(A, num_elem)){
		printf("\nList did not get sorted dummy!\n");
	}else{
		printf("\nEverything went great, the list is sorted!\n");
	}

	// Print the processing time
	printf("\nProcessing (Wall) time: %f s\n\n", (end_time-start_time) );

	free(A);
	return 0;
}

double rand_normal(double mean, double stddev){

	static double n2 = 0.0;
	static int n2_cached = 0;

	if (!n2_cached){
		double x, y, r;

		do{
			x = 2.0*rand()/RAND_MAX - 1;
			y = 2.0*rand()/RAND_MAX - 1;
			r = x*x + y*y;
		}while (r == 0.0 || r > 1.0);

		double d = sqrt(-2.0*log(r)/r);
		double n1 = x*d;
		n2 = y*d;
		double result = n1*stddev + mean;
		n2_cached = 1;
		return result;
	}else{
		n2_cached = 0;
		return n2*stddev + mean;
	}
}

int isSorted(double a[],long int size){

	long int i;
	for (i = 1;i < size;i ++){
		if (a[i] < a[i-1]){
			printf("\nError: at loc %li, %lf < %lf \n", i, a[i], a[i-1]);
			return 0;
		}
	}
	return 1;
}

double get_wall_time(){

	struct timeval time;

	if (gettimeofday(&time,NULL)){
		return 0;
	}

	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void bucket_sort(double *A,long int num_elem,int buckets){

	long int i,n,index;
	int y,z,buckettobe;
	double bucketsize=MAX/buckets;
	struct bin bins[buckets];

	for(y=0;y<buckets;y++){
		bins[y].count = 0;
		bins[y].index = 0;
	}

	// First pass through the array, find out the size of each bin
	for(i=0;i<num_elem;i++){
		buckettobe = (int)(floor(A[i])/bucketsize);
		// If the number is negative (normal distribution), put it in the first bucket
		if (buckettobe < 0){
			buckettobe = 0;
		// If the number is very high (normal distribution), put it in the last bucket
		}else if (buckettobe >= buckets){
			buckettobe=buckets-1;
		}
		bins[buckettobe].count++;
	}

	// Allocate space for each bin
	for(y=0;y<buckets;y++){
		bins[y].B = malloc(bins[y].count * sizeof(double));
		bins[y].place = y;
	}

	// Second pass through the array, assign each element to a bin
	for(i=0;i<num_elem;i++){
		buckettobe = (int)(floor(A[i])/bucketsize);
		// If the number is negative (normal distribution), put it in the first bucket
		if (buckettobe < 0){
			buckettobe = 0;
		// If the number is very high (normal distribution), put it in the last bucket
		}else if (buckettobe >= buckets){
			buckettobe=buckets-1;
		}
		bins[buckettobe].B[bins[buckettobe].index] = A[i];
		bins[buckettobe].index++;
	}

	// Sort the bins in descending order according to their size
	qsort(bins,buckets,sizeof(struct bin),compare_desc);

	// Send the bins for sorting - highest sized bins are sent first, each thread gets a bin in cyclic order.
	#pragma omp parallel
	{
		#pragma omp for nowait
		for(y=0;y<buckets;y++){
			qsort(bins[y].B,bins[y].count,sizeof(bins[y].B[0]),compare);
		}
	}

	index = 0;

	// Join all the sorted bins into the master table
	for(y=0;y<buckets;y++){
		for(z=0;z<buckets;z++){
			if (bins[z].place == y){
				for(n=0;n<bins[z].count;n++){
					A[index] = bins[z].B[n];
					index++;
				}
				break;
			}
		}
	}

	for(y=0;y<buckets;y++){
		free(bins[y].B);
	}

}

static int compare(const void *a, const void *b){

	if (*(double*)a > *(double*)b)
		return 1;
	else if (*(double*)a < *(double*)b)
		return -1;
	else return 0;
}

static int compare_desc(const void *a, const void *b){

	struct bin *b1 = (struct bin *)a;
	struct bin *b2 = (struct bin *)b;

	if (b1->count > b2->count)
		return -1;
	else if (b1->count < b2->count)
		return 1;
	else return 0;
}
