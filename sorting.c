/**
 * MPI Program that performs simple sorting
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

// mpi functions
int MPI_Direct_Merge_sort(int n, double * array, int root, MPI_Comm comm); 
int MPI_Exchange( int n, double * array, int rank1, int rank2, MPI_Comm comm );
int MPI_Sort( int n, double * array, int root, MPI_Comm comm );
int MPI_Sort_ranking(int n, double * a, int root, MPI_Comm comm);
int MPI_Sort_bucket(int n, double * a, double max, int root, MPI_Comm comm);

double * merge( int n, double * array, int m, double * b );
void merge_array(int n, double * a, int m, double * b, double * c);

void     merge_sort(int n, double * a);
void     swap (double * a, double * b);

int main (int argc, char *argv[])
{

	 int size, rank;

	int n = 100000, i, j, k, x, q, l, shell, pair, *nr;
	double m = 10.0;
	double * scattered_array, * arrayA, * arrayB, * arrayC, * arrayD, * arrayE;
	

	// Init + rank + size
	MPI_Init(&argc, &argv);
   	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_Comm_size( MPI_COMM_WORLD, &size );

	scattered_array=(double *)calloc(n*size, sizeof(double));
	arrayA = (double *) calloc( n , sizeof(double) );

	if( rank == 0 )
	{

	   //initialise the array with random values, then scatter to all processors
           
	   srand( ( ( unsigned ) time( NULL ) + rank) );

	   for( i = 0; i < n; i++ )
	   {
	      arrayA[i] =( (double) rand() / RAND_MAX) * m;

	      //printf("%f\n",array[i]);
	   }
	   
	}
	
	// The different sorting algorithms
        MPI_Direct_Merge_sort(n, arrayA, 0, MPI_COMM_WORLD);
	
	MPI_Sort_ranking(n, arrayA, 0, MPI_COMM_WORLD);
	
	MPI_Odd_Even_Sort(n, arrayA, 0, MPI_COMM_WORLD);
	
	MPI_Sort_bucket(n, arrayA, m, 0, MPI_COMM_WORLD);
	
	MPI_Finalize();

}

/**
 * This method determines whether the array is sorted or not.
 */
int is_sorted(int n, double * a, int root, MPI_Comm comm) { 
  
  double * first, * last; 
  int i, rank, size, s, sorted = 1; 

  MPI_Comm_size(comm, &size); MPI_Comm_rank(comm, &rank); 
  // Number of elements in each array 
  s = n / size; 
  
  first = (double *) calloc( size , sizeof(double) ); 
  last = (double *) calloc( size, sizeof(double) ); 
  
  // Gather the first elements
  MPI_Gather(&a[0], 1, MPI_DOUBLE, &first[0], 1, MPI_DOUBLE, root, comm); 
  // Gather the last elements
  MPI_Gather(&a[s - 1], 1, MPI_DOUBLE, &last[0], 1, MPI_DOUBLE, root, comm); 
  
  if (rank == 0) { 
    // Loop through the elements and detemine whether the any of the first elements
    // are actually less than the last elements. If a first element of a later 
    // array is less than a last element of previous array, then we know the lists
    // aren't sorted
    for(i = 0; i < size; i++) { 
      if (first[i + 1] < last[i]) { 
	 sorted = 0;
	 break;
      }   
    }
  }
  
  // Broadcast whether the array is sorted or not.
  MPI_Bcast(&sorted, 1, MPI_INT, root, comm);
  
  return sorted;
  
}

/**
 * BUCKET SORT
 */
int MPI_Sort_bucket(int n, double * a, double max, int root, MPI_Comm comm) {
 
    // dclarations
    int size, rank, i, j, count = 0, * overallCount, * displs, result;
    double * bucket;
    double comTime = 0, commTime = 0, time1;
    
    // find rank and service
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
  
    // allocate space for the buckets (n elements)
   bucket = (double *) calloc (n, sizeof(double));
   overallCount = (int *) calloc (n, sizeof(int));
   displs = (int *) calloc(n, sizeof(int));
  
    // bcast elements of a onto processors
    time1 = MPI_Wtime();
    MPI_Bcast(&a[0], n, MPI_DOUBLE, root, comm);
    
    time1 = MPI_Wtime() - time1;
    commTime += time1;
   
    // traverse and fill bucket for P rank
    time1 = MPI_Wtime();
    for (i = 0; i < n; i++) {
      if ((rank * max /size <= a[i]) && (a[i] < (rank + 1) * max / size)) {
	bucket[count++] = a[i];
      }
    }
    
    // sort the elements in the bucket
    merge_sort(count, bucket);
    //MPI_Odd_Even_Sort(count, bucket, rank, comm);
    
    // Calculate the time to sort
    time1 = MPI_Wtime() - time1;
    comTime += time1;
    
    time1 = MPI_Wtime();
    // gather count for all buckets 
      result = MPI_Gather(&count, 1, MPI_INT, &overallCount[0], 1, MPI_INT, root, comm);
      
      time1 = MPI_Wtime() - time1;
      commTime += time1;
     
      if (result != MPI_SUCCESS) return result;
      
      // Get the displacement for each bucket
      if (rank == root) {	
	displs[0] = 0;
      
	for (i = 1; i < size; i++) {
	    displs[i] = displs[i - 1] + overallCount[i - 1];      
	}
      }
      
    
      // gatherv the buckets
      time1 = MPI_Wtime();
      
      result = MPI_Gatherv(&bucket[0], count, MPI_DOUBLE, &a[0], &overallCount[0], displs, MPI_DOUBLE, root, comm);
      
      time1 = MPI_Wtime() - time1;
      commTime += time1;
      
     // if (rank == 0) {
	//for(i=0; i < n; i++) printf("value: %f\n",a[i]);
      //}

      
      if (result != MPI_SUCCESS) return result;
      
     printf("BUCKET_SORT:Processor %d COMMUNICATED for %lf and COMPUTED for %lf\n", rank, commTime, comTime);
    
    // print execution times - maybe?
    return MPI_SUCCESS;
  
}

/**
 * Direct merge sort.
 */
int MPI_Direct_Merge_sort(int n, double * array, int root, MPI_Comm comm ) {
  
  double * scattered_array;
  int rank, size, error, i;
  double comTime = 0, commTime = 0, time1;
  
  // get rank and size of comm
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  
  // Set the size of the scattered array
  scattered_array = (double *) calloc(n/size, sizeof(double));
  
  
  time1 = MPI_Wtime();
  
  // scatter array to local array (scattered_array)
  error = MPI_Scatter(
    array, 
    n/size, 
    MPI_DOUBLE,
    scattered_array,
    n/size,
    MPI_DOUBLE,
    root,
    comm);
  
  // Calculate the communication time to scatter the elements
  time1 = MPI_Wtime() - time1;
  commTime += time1;
  
  // If the scatter fails we exit the sort
  if (error != MPI_SUCCESS) return error;
  
  time1 = MPI_Wtime();
  
  // Merge sort the scattered elements
  merge_sort(n/size, scattered_array);
  
  // Calculate the time to sort the scattered elements
  time1 = MPI_Wtime() - time1;
  comTime += time1;
  
  time1 = MPI_Wtime();
  // Gather the elements back to the root
  error = MPI_Gather (
    scattered_array,
    n/size, 
    MPI_DOUBLE, 
    array, 
    n/size, 
    MPI_DOUBLE, 
    root, 
    comm);
  
  time1 = MPI_Wtime() - time1;
  commTime += time1;
  
  time1 = MPI_Wtime();
  
  if (rank == root){
    // Merge the elements the arrays togethar
    for (i = 1; i < size; i++)
    {
      merge_array(i * n/size, array, n/size, array + i * n/size, array);
    }
    
  }
  
  //if (rank == 0) {
  //    for(i=0; i<n; i++) printf("%f\n",array[i]);
  //}
  
  // Calculate the time to compute the merges
  time1 = MPI_Wtime() - time1;
  comTime += time1;
  
  printf("DIRECT_SORT: Processor %d COMMUNICATED for %lf and COMPUTED for %lf\n", rank, commTime, comTime);
  
  if (error != MPI_SUCCESS) return error;
  
  return MPI_SUCCESS; 
}

/**
 * Odd even sort.
 */
int MPI_Odd_Even_Sort( int n, double * array, int root, MPI_Comm comm ) {
      
      int rank, size, step, result, i;
      
      double * local_a;
      double comTime = 0, commTime = 0, time1;
      
      // get rank and size of comm
      MPI_Comm_size(comm, &size);
      MPI_Comm_rank(comm, &rank);
      
      //allocate space for numElements/numProcessors amount of doubles
      local_a = (double *) calloc(n/size, sizeof(double));   
      
      //scatter the elements of array "a" and Calculate time communication time.
      time1 = MPI_Wtime();
      result = MPI_Scatter(&array[0], n/size, MPI_DOUBLE, &local_a[0], n/size, MPI_DOUBLE, root, comm);
      time1 = MPI_Wtime() - time1;
      commTime += time1;
      
      
      // Result must be succesful
      if (result != MPI_SUCCESS) {
	  return result;
      }
       
      //sort local array using mergeSort calculate the time it takes
      time1 = MPI_Wtime();
      // Merge sort the arrays and 
      merge_sort(n/size, local_a);
      time1 = MPI_Wtime() - time1;
      comTime += time1;
      
      time1 = MPI_Wtime();
      //odd-even iterations
      for (step = 0; step < size; step ++) {
	  // check if even
	  if ((rank + step) % 2 == 0) {
	      // exchange between one processor and the next
	      if (rank < size - 1) {
		  MPI_Exchange(n/size, local_a, rank, rank + 1, comm);
	      }
	  } else {
              // Exchange between one processor and the previous
	      if (rank > 0) {
		  MPI_Exchange(n/size, local_a, rank - 1, rank, comm);
	      }
	  }
	  
	  MPI_Barrier(comm);
	  
	  // do MPI_Is_sorted() test
	  is_sorted(n, local_a, root, comm);
      }
      
      // Calculate the time to make the exchanges
      time1 = MPI_Wtime() - time1;
      commTime += time1;
      
      //gather local_a
      time1 = MPI_Wtime();
      result = MPI_Gather(&local_a[0], n/size, MPI_DOUBLE, &array[0], n/size, MPI_DOUBLE, root, comm);
      time1 = MPI_Wtime() - time1;
      
      comTime += time1;
      
      // print array
     // if (rank == root) {

	//  for(i=0; i<n; i++) printf("%f\n",array[i]);
	//}
	
      // Result must be succesful
      if (result != MPI_SUCCESS) {
	  return result;
      }
      
      printf("OddEven: Processor %d COMMUNICATED for %lf and COMPUTED for %lf\n", rank, commTime, comTime);
      
      return MPI_SUCCESS;
}

/**
 * Ranking sort.
 * 
 */
int MPI_Sort_ranking(int n, double * a, int root, MPI_Comm comm)
{
	
	int rank, size, i, j, *ranking, * overallRanking, result;
	double * b;
	double comTime = 0, commTime = 0, time1;
	
	// find rank and size
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	
	// allocate the extra memory / arrays needed
	ranking = (int *) calloc(n/size, sizeof(int));
	overallRanking = (int *) calloc(n, sizeof(int));
	b = (double *) calloc(n, sizeof(double));
	
	// Brodcast the array to all processo in the communicator
	time1 = MPI_Wtime();
	MPI_Bcast(&a[0], n, MPI_DOUBLE, root, comm);
	time1 = MPI_Wtime() - time1;
	commTime += time1;
	
	time1 = MPI_Wtime();
        // 
	for (i = 0; i < n/size; i++) {
	    ranking[i] = 0;
	    
	    for (j = 0; j < n; j++) {
		if (a[j] < a[i + rank * n/size]) ranking[i]++;
	    }
	  
	}
	
	time1 = MPI_Wtime() - time1;
	comTime += time1;
	
	// Gather the array ranking to finalRanking
	time1 = MPI_Wtime();
	result = MPI_Gather(&ranking[0], n/size, MPI_INT, &overallRanking[0], n/size, MPI_INT, root, comm);
	time1 = MPI_Wtime() - time1;
	commTime += time1;
        
        if (result != MPI_SUCCESS) {
	  return result;
        }
	
	// if processor 0 then restore the order in the array b and move b back to a
	time1 = MPI_Wtime();
	if (rank == root) {
	    // restore order in b
	    for (i = 0; i < n; i++) {
		b[overallRanking[i]] = a[i]; 
	    }
	  
	    for (i = 0; i < n; i++) {
		a[i] = b[i]; 
	    }
	 // for(i=0; i<n; i++) printf("%f\n",a[i]);
	}
	
	
	
	time1 = MPI_Wtime() - time1;
	comTime += time1;
	
	printf("RANKING_SORT: Processor %d COMMUNICATED for %lf and COMPUTED for %lf\n", rank, commTime, comTime);
	
	return MPI_SUCCESS;
	
}

/**
 * 
 */
int MPI_Exchange( int n, double * array, int rank1, int rank2, MPI_Comm comm ) {
      int rank, size, result, i, tag1 = 0, tag2 = 1;
      double * b = ( double * ) calloc( n, sizeof( double ) );
      double * c;
       
      MPI_Status status;
      MPI_Comm_rank( comm, &rank );
      MPI_Comm_size( comm, &size );
     
      //L8.6
      if( rank == rank1 )
      {
        result = MPI_Send( &array[ 0 ], n, MPI_DOUBLE, rank2, tag1, comm );
        result = MPI_Recv( &b[ 0 ], n, MPI_DOUBLE, rank2, tag2, comm, &status );
        c = merge( n, array, n, b );
        for( i = 0; i < n; i++ )
        {
          array[ i ] = c[ i ];
        }
      }
      else if( rank == rank2 )
      {
        result = MPI_Recv( &b[ 0 ], n, MPI_DOUBLE, rank1, tag1, comm, &status );
        result = MPI_Send( &array[ 0 ], n, MPI_DOUBLE, rank1, tag2, comm) ;
        c = merge( n, array, n, b );
        for( i =0; i < n; i++ )
        {
          array[ i ] = c[ i + n ];
        }
      }
      return MPI_SUCCESS;
}

// function to merge sort the array a with n elements

void merge_sort(int n, double * a){

   double * c;
   int i;

   if (n<=1) return;

   if(n==2) {

      if(a[0]>a[1])swap(&a[0],&a[1]);
      return;
   }

   merge_sort(n/2,a);
   merge_sort(n-n/2,a+n/2);

   c = merge(n/2,a,n-n/2,a+n/2);

   for(i=0;i<n;i++)a[i]=c[i];

   return;
}
    
//notes
double * merge( int n, double * a, int m, double * b ) {
       int i, j, k;
       double * c = ( double * ) calloc( n + m, sizeof( double ) );
     
       for( i=j=k=0; ( i < n ) && ( j < m ); )
       {
          if( a[ i ] <= b[ j ] )
          {
            c[ k++ ] = a[ i++ ];
          }
          else
          {
            c[ k++ ] = b[ j++ ];
          }
       }
      if( i == n )
      {
        for( ; j < m; )
        {
          c[ k++ ] = b[ j++ ];
        }
      }
      else
      {
        for( ; i < n; )
        {
          c[ k++ ] = a[ i++ ];
        }
      }
      return c;
}

void merge_array(int n, double * a, int m, double * b, double * c) {
   // Initialize iterators.
   int i = 0, j = 0, k = 0;

   // Initialize temp merged_array.
   double * local_c = (double *) calloc(n+m, sizeof(double));

   while((i < n) && (j < m)) {
      if(a[i] <= b[j]) {
          local_c[k++] = a[i++];
      }
      else {
          local_c[k++] = b[j++];
      }
   }

   if(i==n) {
       while(j<m) {
           local_c[k++] = b[j++];
       }
   }

   else {
       while(i<n) {
           local_c[k++] = a[i++];
       }
   }

   for (i = 0; i < n+m; i++) {
       c[i] = local_c[i];
   }
}





// swap two doubles
void swap (double * a, double * b){

   double temp;

   temp=*a;*a=*b;*b=temp;

}