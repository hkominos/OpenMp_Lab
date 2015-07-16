#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#define NUM_ELEMENTS 1000
#define MAX 1000

struct number{
	double numb;
	struct number *next;
};

void copy2array(struct number *ptr_table[],double* matrixes[],int matrixes_size[],int index);
void putinbucket(struct number* new,int buckettobe,struct number *ptr_table[],int buckets);
static int compare (const void * a, const void * b);
void bucket_sort(double *A,int num_elem,int buckets);


double get_wall_time(){

	struct timeval time;

	if (gettimeofday(&time,NULL)){
		return 0;
	}

	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int isSorted(double a[], int size)
{
	int i;
	for (i = 1;i < size;i ++){
		if (a[i] < a[i-1]){
			printf("\nat loc %d, %e < %e \n", i, a[i], a[i-1]);
			return 0;
		}
	}
	return 1;
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

int main(int argc, char *argv[]){

	int i,num_elem,j,buckets,x;
	double *A;
	double start_time,end_time;

	num_elem = NUM_ELEMENTS;

	if (argc == 2){
		// First argument is the size of the list
		num_elem = atoi(argv[1]);
	}else{
		printf("\nNo arguments given, continuing with default value: NUM_ELEMENTS 1000\n");
	}

	A = malloc(num_elem * sizeof(double));

	// Choose one seed - either srand48 for drand48 or srand for rand_normal

	//srand48((unsigned int)time(NULL));
	srand(time(NULL));

	// Initialize the list with random values - drand48 for uniform distribution and rand_normal for normal distribution
	for (i=0;i<num_elem;i++){
		//A[i] = drand48() * MAX;
		A[i] = rand_normal(500.0, 250.0);
	}

	// By default openmp assigns as many threads as the available cores
	x = omp_get_max_threads();
	// Assign 5 more threads than the available cores
	omp_set_num_threads(x+5);
	printf("\nNumber of threads: %d\n", x+5);
	// Assign the same number of buckets as threads
	buckets = x+5;

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
	printf("Processing (Wall) time: %f s\n\n", (end_time-start_time) );

}


void bucket_sort(double *A,int num_elem,int buckets){

	int bucketsize=MAX/buckets;
	int i,j;
	double* matrixes[buckets];
	struct number* ptr_table[buckets];
	int matrixes_size[buckets];

	// Creating pointer table
	for(i=0;i<buckets;i++){
		ptr_table[i]=NULL;
	}

	/* Put every element of the input array into a bucket. To decide which bucket to put it into,
	 * we divide the element by (MAX/buckets). We have a MAX of 1000 and so our elements are of 3
	 * digits (in the 0-1000 range). So for example, if the element 345 comes up and we have 2 buckets,
	 * 345/500 = 0.69. (int)0.69 = 0 so it goes into the first bucket. For normal distribution, we just
	 * have to take care to catch the wild numbers on both ends - put all the really small elements
	 * into the first bucket and all the very large elements into the last bucket.*/
	for(i=0;i<num_elem;i++){
		struct number* new=malloc(sizeof(struct number));
		new->numb=A[i];
		new->next=NULL;
		putinbucket(new,(int)(new->numb/bucketsize),ptr_table,buckets);
	}

	/* We use the schedule(auto) command to let the run-time system and/or compiler decide automatically
	 * and handle the work load for us.*/
	#pragma omp parallel
	{
		//#pragma omp for nowait schedule(runtime)
		#pragma omp for nowait schedule(auto)
		for(i=0;i<buckets;i++){
			// Every bucket is a linked list up until now. Convert them into arrays one by one and sort them with qsort
			copy2array(ptr_table,matrixes,matrixes_size,i);
			qsort(matrixes[i],matrixes_size[i],sizeof(double),compare);
		}
	}

	// Join all the sorted buckets into the master table
	int index=0;

	for(i=0;i<buckets;i++){
		for(j=0;j<matrixes_size[i];j++){
			A[index++]=*(matrixes[i]+j);
			}
	}
}

void copy2array(struct number *ptr_table[],double* matrixes[],int matrixes_size[],int index){

	int newsize=1;
	int j=0;

	double *table=(double*)malloc(1*sizeof(double));
	table[j++]=ptr_table[index]->numb;
	struct number *temp=ptr_table[index]->next;

	// Copy the linked list to an array
	do{
		newsize+=1;
		table=(double*)realloc(table,newsize*sizeof(double));
		table[j++]=temp->numb;
		temp=temp->next;
	}while(temp!=NULL);

	matrixes_size[index]=j;
	matrixes[index]=table;
	struct number *temp2=ptr_table[index];

	// Free the linked list - we don't need it anymore
	do{
		temp=temp2;
		temp2=temp2->next;
		free(temp);
	}while(temp2!=NULL);
}





void putinbucket(struct number* new,int buckettobe,struct number *ptr_table[],int buckets){

	// If the number is negative (normal distribution), put it in the first bucket
	if (buckettobe < 0){
		buckettobe = 0;
	}

	// If the number is very high (normal distribution), put it in the last bucket
	if (buckettobe >= buckets){
		buckettobe=buckets-1;
	}

	if (ptr_table[buckettobe]==NULL){
		ptr_table[buckettobe]=new;
	}else{
		struct number *temp=ptr_table[buckettobe];
		ptr_table[buckettobe]=new;
		new->next=temp;
	}
}

static int compare(const void *a, const void *b){

	if (*(double*)a > *(double*)b)
		return 1;
	else if (*(double*)a < *(double*)b)
		return -1;
	else return 0;
}

