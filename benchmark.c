//
//  main.c
//  PrimeNumbersParallelOpenMP
//
//  Created by Luís Felipe Rabello Taveira on 17/03/15.
//  Copyright (c) 2015 Luís Felipe Rabello Taveira. All rights reserved.
//

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

int IsPrime(unsigned int number);
double process (int n,int*out, char* out_type, int nthreads, omp_sched_t kind, int chunk_size);
void parse_params (int argc, char *argv[], int* n, char* out_type);
void benchmark();

int main (int argc, char *argv[])
{
    benchmark();
    return 0;
    
}
int IsPrime(unsigned int number) {
    if (number <= 1) return 0; // zero and one are not prime
    if(number != 2 && number % 2 == 0) return 0; //there is no need to check even numbers
    unsigned int i;
    for (i=3; i*i<=number; i+=2) {
        if (number % i == 0) return 0;
    }
    return 1;
}
double process (int n, int * out, char* out_type, int nthreads, omp_sched_t kind, int chunk_size){
    
    omp_set_schedule(kind,chunk_size);
    
    double begin, end;
    double time_spent;
    begin = omp_get_wtime();
    //Put number 2 in the list
    if (n > 2) {
        out[0] = 2;
        out[(int)ceil(n/2)] = 0;
    }
    //Check which numbers greater than 2 are prime (no need to check even numbers)
#pragma omp parallel num_threads (nthreads)
    {
#pragma omp for schedule (runtime)
        for(int i=3 ; i<=n; i+=2)
        {
            if(IsPrime(i)){
                out[ (int)floor(i/2)] = i;
            }
            else {
                out[(int)floor(i/2)]=0;
            }
        }
    }
    end = omp_get_wtime();
    time_spent = (end - begin);
    
    return time_spent;
}

//benchmark routine
void benchmark(){
    
    int* out = (int*) malloc ( 100000000 * (sizeof(int)));
    if (out == NULL){
        printf("Problems! Memory allocation related. Quitting program...\n");
        exit(2);
    }
    char benchmark[] = "benchmark.output";
    omp_sched_t kind;
    double time;
    
    FILE* ofp;
    ofp = fopen (benchmark, "w");
    if (ofp == NULL) {
        printf("Problems! related to write permission to directory. Quitting program...\n");
        exit(2);
    }
    //Ranging the size of n
    for(int N=pow(10,6) ; N<= pow(10,8); N*=10){
        fprintf(ofp,"n = %d\n",N);
        printf("Running with n=%d\n", N);
        //Ranging the chunk_size
        for (int CHUNK_SIZE=2500;CHUNK_SIZE <= 10000;CHUNK_SIZE= CHUNK_SIZE*2){
            fprintf(ofp, "\tchunk_size = %d\n",CHUNK_SIZE);
            printf("\tRunning with Chunk_Size = %d\n", CHUNK_SIZE);
            //Ranging the scheduler kind
            for (int k=0;k<=2;k++){
                if (k==0){kind = omp_sched_static;
                    fprintf(ofp, "\t\tscheduler = Static\n");}
                else if (k==1){ kind = omp_sched_dynamic;
                    fprintf(ofp, "\t\tscheduler = Dynamic\n");}
                else {kind = omp_sched_guided;
                    fprintf(ofp, "\t\tscheduler = Guided\n");}
                //Ranging the number of threads
                for(int nthreads=1; nthreads<=8; nthreads*=2){
                    fprintf(ofp,"\t\t\tnthreads = %d\n",nthreads);
                    time = 0;
                    //Repeat each test 3 times
                    for(int j=1; j<=3; j++)time += process (N, out, "t", nthreads, kind, CHUNK_SIZE);
                    time = time/3;
                    fprintf(ofp,"\t\t\t\texec_time = %lf\n",(time));
                }
            }
        }
    }
    fclose(ofp);
    
}
//
//  main.c
//  PrimeNumbersParallelOpenMP
//
//  Created by Luís Felipe Rabello Taveira on 17/03/15.
//  Copyright (c) 2015 Luís Felipe Rabello Taveira. All rights reserved.
//
