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
    int n;
    char out_type [10];
    parse_params(argc, argv, &n, out_type);
    //Allocate half the size of the input (no need to check even numbers)
    int* out = (int*) malloc ( ceil(n/2)* (sizeof(int)));
    if (out == NULL){
        printf("Problems! Memory allocation related. How big is your memory? Quitting program...\n");
        exit(2);
    }
    
    int nthreads = 8;
    omp_sched_t kind = omp_sched_guided;
    int chunk_size = 64;
    process (n, out, out_type, nthreads, kind, chunk_size);
    //benchmark();
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
    if (n >= 2) {
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
    
    //Print to stdout
    if(out_type[0] == 'l' || out_type[0] == 'a'){
        for(int k=0; k<=floor(n/2); k++){
            if(out[k]!=0) printf("%d ",out[k]);
        }
        printf("\n");
    }
    
    double end_print = omp_get_wtime();
    double time_spent_after_print;
    time_spent_after_print =(end_print - begin);
    
    //Print to stdout
    if(out_type[0] == 'a' || out_type[0] == 't')printf("%lf\n",time_spent_after_print);
    
    
    return time_spent;
}

//benchmark routine
void benchmark(){
    
    int* out = (int*) malloc ( ceil(100000000/2) * (sizeof(int)));
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

void parse_params (int argc, char *argv[], int* n, char* out_type){
    
    if(argc != 3) {
        printf ("\n Usage: %s NUMBER a \n or\n Usage: %s NUMBER l \n or \n Usage: %s NUMBER t\n\n", argv[0],argv[0],argv[0]);
        printf("Quitting application...\n");
        exit (1);
    }
    else{
        *n = atoi(argv[1]);
        strcpy(out_type,argv[2]);
        if(out_type[0] != 'a' && out_type[0] != 'l' && out_type[0] != 't'){
            printf("Bad params... \n");	
            printf ("\n Usage: %s NUMBER a \n or\n Usage: %s NUMBER l \n or \n Usage: %s NUMBER t\n\n", argv[0],argv[0],argv[0]);
            printf("Quitting application...\n");
            exit (1);
            
        }		
    }
}
//
//  main.c
//  PrimeNumbersParallelOpenMP
//
//  Created by Luís Felipe Rabello Taveira on 17/03/15.
//  Copyright (c) 2015 Luís Felipe Rabello Taveira. All rights reserved.
//
