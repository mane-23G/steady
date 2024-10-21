/*
  Title          : steady_state.c
  Author         : Mane Galstyan
  Created on     : Apri 17, 2024
  Description    : parralize a random walk function
  Purpose        : To use MPI broadcast MPI scatter
  Usage          : mpirun --use-hwthread-cpus steady_state <file> <rows> <columns>
  Build with     : mpicc -Wall -o steady_state steady_state.c
  Modifications  : comments
 
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define SUCCESS             0
#define MALLOC_ERROR        1

#define CONVERGENCE_THRESHOLD  0.05
#define ROOT 0

#define NORTH   1
#define EAST    2
#define SOUTH   3
#define WEST    4

typedef struct point2d_tag {
    int x;
    int y;
} point2d;

const point2d East  = {1, 0};
const point2d West  = {-1,0};
const point2d North = {0, 1};
const point2d South = {0,-1};

/* prints error*/
void fatal_error(int errornum, const char *msg) {
    fprintf(stderr, "%s\n", msg);
    printf("%s\n",msg);
    MPI_Abort(MPI_COMM_WORLD, 1); // Use MPI_Abort to terminate all processes.
}

/* prints usage error*/
void usage_error(const char *msg) {
    fprintf(stderr, "Usage: %s", msg);
    printf("%s\n",msg);
    MPI_Abort(MPI_COMM_WORLD, 1); // Use MPI_Abort to terminate all processes.
}

/* function to allocate a matrix*/
void alloc_matrix(
        int     nrows,          /* number of rows in matrix                   */
        int     ncols,          /* number of columns in matrix                */
        size_t  element_size,   /* number of bytes per matrix element         */
        void  **matrix_storage, /* address of linear storage array for matrix */
        void ***matrix,         /* address of start of matrix                 */
        int    *errvalue)       /* return code for error, if any              */
{
    int   i;
    void *ptr_to_row_in_storage;
    void **matrix_row_start;
    size_t total_bytes;
    total_bytes = nrows * ncols * element_size;

    *matrix_storage = malloc(total_bytes);
    if ( NULL == *matrix_storage ) {
        *errvalue = MALLOC_ERROR;
        return;
    }

    memset(*matrix_storage, 0, total_bytes );

    *matrix = malloc (nrows * sizeof(void*));
    if ( NULL == *matrix ) {
        *errvalue = MALLOC_ERROR;
        return;
    }

    matrix_row_start = (void*) &(*matrix[0]);

    ptr_to_row_in_storage = (void*) *matrix_storage;

    for ( i = 0; i < nrows; i++ ) {
        *matrix_row_start = (void*) ptr_to_row_in_storage;
        matrix_row_start++;
        ptr_to_row_in_storage +=  ncols * element_size;
    }
    *errvalue = SUCCESS;
}

/* function to calculate where each processor will get from the array*/
void displace_matrix(int displs[], int distribute[], int cols, int p) {
    displs[0] = 0;
    int c = 0;
    for(int i = 1; i < p; i++) {
        c += distribute[i-1];
        displs[i] = displs[i-1] + distribute[i-1] * cols;
    }
}

/* function to allocate amount of elements in the array and how many rows each processor gets*/
void distribute_matrix(int rows, int cols, int distribute[], int amount[], int p) {
    for(int i = 0; i < p; i++) {
        distribute[i] = ((((i+1)*rows)/p) - ((i * rows)/p));
        amount[i] = distribute[i] * cols;
    }
}

/* function for each processor to generate random numbers based on id*/
void  init_random( int id ) {
    unsigned int seed = (unsigned int) time(NULL) + id;
    srandom(seed);
}

/* function called to get random double*/
double uniform_random() {
    return (double) (random()) / RAND_MAX;
}

/* calculate randome value and choose direction based on the value*/
point2d next_dir() {
    double u;
    u = uniform_random();
    if ( u < 0.25 )
        return North;
    else if ( u < 0.50 )
        return East;
    else if ( u < 0.75 )
        return South;
    else
        return West;
}

/* return the location used to calculate avg temp*/
int on_boundary(point2d point, int cols, int rows) {
    if ( 0 == point.x )
        return WEST;
    else if ( cols -1 == point.x )
        return EAST;
    else if ( 0 == point.y )
        return NORTH;
    else if ( rows - 1 == point.y )
        return SOUTH;
    else
        return 0;
}

/* adds two points together and return the new point*/
point2d next_point(point2d oldpoint, point2d direction) {
    point2d temp;
    temp.x = oldpoint.x + direction.x;
    temp.y = oldpoint.y + direction.y;
    return temp;
}

/* function used to convert command line arg coords to valid number*/
/* returns number or call fatal error to quit if arg is not valid*/
int arg_to_num(char *arg, int interval) {
    errno = 0;
    int arg_num = strtol(arg, NULL, 10);
    if(errno == ERANGE || arg_num < 0 || arg_num >= interval) {
        fatal_error(errno, "Coordinate is invalid");
    }
    return arg_num;
}


int main(int argc, char *argv[]) {
    /* initilazie MPI */
    /* rank, size*/
    int id, p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    char usage_msg[512];
    
    
    /* if less than 4 args returns usage error*/
    if (argc != 4 && id == ROOT) {
        sprintf(usage_msg, "%s <file> <row> <column>", argv[0]);
        usage_error(usage_msg);
    }
    
    int rows, columns; /* total number or rows and cols*/
    int expected_input = 6; /* expected number of file arguments*/
    double boundary_temp[4]; /* array to hold north,east,south,west temps*/
    int xcord,ycord; /* coord provided by user*/
    double elapsed_time; /* used to keep track of elapsed time*/
    
    /* used to allocate matrix for ROOT and the partial matrixs*/
    double **plate, **local_plate;
    double *plate_storage, **local_plate_storage;
    int error = 0;
    
    int distribute[p]; /* the amount of rows each processor gets            */
    int displs[p];     /* the place from the array each processor gets      */
    int amount[p];     /* the number of element arrays each processor gets  */
    
    
    /* parse the command line and set plate up for distrubution with the ROOT*/
    if (ROOT == id) {
        /* open the file and read contents into varibles*/
        errno = 0;
        FILE *file = fopen(argv[1], "r");
        if(file == NULL) {
            fatal_error(errno, "Error opening file");
        }
        if(expected_input != fscanf(file, "%d %d %lf %lf %lf %lf", &rows, &columns, &boundary_temp[0], &boundary_temp[1], &boundary_temp[2], &boundary_temp[3])) {
            fatal_error(errno, "File is in incorrect format: <rows> <columns> <north_temp> <east_temp> <south_temp> <west_temp>");
        }
        if(errno != 0) {
            fatal_error(errno, "Error reading file");
        }
        fclose(file);
        
        /* if there is not enough matrix to calculate temp distribution*/
        if(rows < 3 || columns < 3) {
            fatal_error(errno, "Grid is too small");
        }
        
        /* convert user given coordinate values into numbers*/
        ycord = arg_to_num(argv[2], rows);
        xcord = arg_to_num(argv[3], columns);
        
        /* allocate plate matrix using alloc_matrix*/
        alloc_matrix(rows, columns, sizeof(double), (void**) &plate_storage, (void***) &plate, &error);
        if(error != 0) {
            fatal_error(error, "Error while creating array");
        }
        
        /* set the corners of the matrix to average of its respective sides*/
        plate[0][0] = (boundary_temp[0] + boundary_temp[3])/2; /* avg of north + west*/
        plate[0][columns-1] = (boundary_temp[0] + boundary_temp[1])/2; /* avg of north + east*/
        plate[rows-1][0] = (boundary_temp[3] + boundary_temp[2])/2; /* avg of west + south*/
        plate[rows-1][columns-1] = (boundary_temp[1] + boundary_temp[2])/2; /* avg of south + east*/

        
        /* set the top and bottom borders to north and south*/
        for(int i = 1; i < columns-1; i++) {
            plate[0][i] = boundary_temp[0];
            plate[rows-1][i] = boundary_temp[2];
        }

        
        /* set the side borders to west and east*/
        for(int i = 1; i < rows-1; i++) {
            plate[i][0] = boundary_temp[3];
            plate[i][columns-1] = boundary_temp[1];
        }

        /* set the inside (excluding border and corners) to zero*/
        for (int i = 1; i < rows-1 ; i++ ) {
            for (int j = 1; j < columns -1 ; j++) {
                plate[i][j] = 0.0;
            }
        }
        
        /* set amount, distrbute, and displs arrays*/
        distribute_matrix(rows, columns, distribute, amount, p);
        displace_matrix(displs, distribute, columns, p);
    }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = MPI_Wtime();
    
    /* broadcast the distribute, displs, amount, boundary_temp, rows and cols to the rest of the processors*/
    MPI_Bcast(distribute, p, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(displs, p, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(amount, p, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(boundary_temp, 4, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&rows, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&columns, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    init_random(id); /* call for the processors to be able to produce random numbers*/
    
    int local_rows = distribute[id]; /* the number of rows the processor is in charge of*/
    
    /* allocate space in size of local rows and columns in prepartion for receving parts of the plate*/
    alloc_matrix(distribute[id], columns, sizeof(double), (void**) &local_plate_storage, (void***) &local_plate, &error);
    if(error != 0) {
        fatal_error(error, "error allocating memeory");
    }
    
    /* send each processor their portion of the array*/
    MPI_Scatterv(plate_storage, amount, displs, MPI_DOUBLE, *local_plate, amount[id], MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    
    /* calculate the row start and end since two process will have to de mindful of the north and south border */
    int row_start, row_end;
    if(id == ROOT) {
        row_start = 1;
        row_end = local_rows;
    }
    else if(id == p - 1) {
        row_start = 0;
        row_end = local_rows-1;
    }
    else {
        row_start = 0;
        row_end = local_rows;
    }
    /* start the random walk and adjust each plate cord so*/
    int count = 0;
    int location;
    point2d current, next;
    double old_value, diff, maxdiff;
    while(count < 10000) {
        maxdiff = 0;
        for(int i = row_start; i < row_end; i++) {
            for(int j = 1; j < columns - 1; j++) {
                current.x = j;
                current.y = i + distribute[i];
                while(0 == (location = on_boundary(current,columns,rows))) {
                    next = next_point(current, next_dir(id));
                    current = next;
                }
                old_value = local_plate[i][j];
                local_plate[i][j] = (old_value*count + boundary_temp[location-1]) / (count+1);
                diff = fabs(local_plate[i][j] - old_value);
                if(diff > maxdiff) {
                    maxdiff = diff;
                }
            }
        }
        if(maxdiff < CONVERGENCE_THRESHOLD) {
            break;
        }
        else count++;
    }
    
    /* gather the local_plates into the main plate in the root*/
    MPI_Gatherv(local_plate_storage, amount[id], MPI_DOUBLE, plate_storage, amount, displs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    
    /* free the memory*/
    free(local_plate);
    free(local_plate_storage);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();

    /* ROOT prints the coordinate wanted by the user than frees its own memory*/
    if(id == ROOT) {
        printf("%.2lf\n", plate[ycord][xcord]);
//        printf("%.6lf sec\n", elapsed_time);
        free(plate);
        free(plate_storage);
    }
    
    /* finialize mpi*/
    MPI_Finalize();
    return 0;
}
