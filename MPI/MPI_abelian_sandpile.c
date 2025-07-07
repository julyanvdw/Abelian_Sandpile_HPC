#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>

int ROWS;
int COLS;
int INIT_VALUE;

//Per process global variables
int** local_chunk = NULL; 
int rows_per_chunk,  cols_per_chunk;    

//Printing functions
void print_2d_array_to_file(int** array, int rows, int cols) {
    FILE* file = fopen("printed.txt", "w");
    if (file == NULL) {
        perror("Failed to open printed.txt");
        return;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%4d ", array[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

//Reference: ChatGPT-4.1 VS Copilot
void save_as_ppm(const char *filename, int** array) { 
    FILE *fp = fopen(filename, "w");
   
    // PPM header
    fprintf(fp, "P3\n%d %d\n255\n", COLS, ROWS); // Exclude sink border

    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            int val = array[i][j];
            int r=0, g=0, b=0;
            if (val == 1) { g = 255; }
            else if (val == 2) { b = 255; }
            else if (val == 3) { r = 255; }
            // else (val == 0) stays black
            fprintf(fp, "%d %d %d ", r, g, b);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

//Setup functions
int** allocate_2d_array(int rows, int cols) {
    int** array = (int**)malloc(rows * sizeof(int*));
    
    for (int i = 0; i < rows; i++) {
        array[i] = (int*)malloc(cols * sizeof(int));
    }
    
    return array;
}

void setup(int** array) {

    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            array[i][j] = INIT_VALUE;
        }
    }
}

// Partition array using MPI_Dims_create and MPI_Cart_create, then scatter to all ranks
MPI_Comm partition_array(int** array, int rank, int size) {

    //MPI TOPOLOGY SET UP
    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);

    // Check if grid can be evenly divided among processes (assert that the input is valid)
    if (ROWS % dims[0] != 0 || COLS % dims[1] != 0) {
        if (rank == 0) {
            printf("Error: Grid size (%d x %d) is not divisible by process grid (%d x %d).\n",
                    ROWS, COLS, dims[0], dims[1]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int periods[2] = {0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);

    //get the coords of this rank
    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    //calculate the local chunk dimensions
    rows_per_chunk = ROWS / dims[0];
    cols_per_chunk = COLS / dims[1];

    // Each processor will get a chunk of size rows_per_proc x cols_per_proc
    int chunk_size = rows_per_chunk * cols_per_chunk;
    int* local_chunk_buff = (int*)malloc(chunk_size * sizeof(int));
    int* sendbuf = NULL;

    if (rank == 0) {
        sendbuf = (int*)malloc(size * chunk_size * sizeof(int));
        int processor_coords[2];
        int start_row, start_col;

        //Flatten the entire global array into a 1D array to be scattered
        for (int p = 0; p < size; p++) {
            MPI_Cart_coords(cart_comm, p, 2, processor_coords);
            start_row = processor_coords[0] * rows_per_chunk;
            start_col = processor_coords[1] * cols_per_chunk;

            for (int i = 0; i < rows_per_chunk; i++) {
                for (int j = 0; j < cols_per_chunk; j++) {
                    sendbuf[p * chunk_size + i * cols_per_chunk + j] = array[start_row + i][start_col + j];
                }
            }
        }
    }

    // Scatter the chunks to all processors
    MPI_Scatter(sendbuf, chunk_size, MPI_INT, local_chunk_buff, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate global local_chunk pointer then reconstruct the 2D array
    local_chunk = (int**)malloc(rows_per_chunk * sizeof(int*));

    for (int i = 0; i < rows_per_chunk; i++) {
        local_chunk[i] = (int*)malloc(cols_per_chunk * sizeof(int));
        for (int j = 0; j < cols_per_chunk; j++) {
            local_chunk[i][j] = local_chunk_buff[i * cols_per_chunk + j];
        }
    }

    // Free buffers after we're done
    if (rank == 0) {
        free(sendbuf);
    }

    free(local_chunk_buff);

    //return the MPI Com to be used in execute
    return cart_comm;
}

int probe_and_receive(int rank, int source_rank, char direction, int** local_array, int chunk_height, int chunk_width, MPI_Comm cart_comm) {
    
    int flag;
    MPI_Status status;

    MPI_Iprobe(source_rank, 0, cart_comm, &flag, &status);
    if (!flag) return 0;

    int count; 
    MPI_Get_count(&status, MPI_INT, &count);

    int* buffer = malloc(count * sizeof(int));
    MPI_Recv(buffer, count, MPI_INT, source_rank, 0, cart_comm, &status);

    //halos are packages as follows: <edge_coord>,<value>,<edge_coord>,<value> etc
    //go through the halo, get the coordinate , value pair, then apply the value at that coordinate at that edge

    // Apply received halo data
    int index, value;
    for (int i = 0; i < count; i += 2) {
        index = buffer[i];
        value = buffer[i + 1];

        if (direction == 't') {
            local_array[0][index] += value;
        } else if (direction == 'b') {
            local_array[chunk_height - 1][index] += value;
        } else if (direction == 'l') {
            local_array[index][0] += value;
        } else if (direction == 'r') {
            local_array[index][chunk_width - 1] += value;
        }
    }
    free(buffer);
    return 1;
}

void gather_and_stitch(int** global_array, MPI_Comm cart_comm, int rank, int size) {
    int chunk_size = rows_per_chunk * cols_per_chunk;

    // Allocate send buffer and flatten local_chunk
    int* sendbuf = malloc(chunk_size * sizeof(int));
    for (int i = 0; i < rows_per_chunk; i++) {
        for (int j = 0; j < cols_per_chunk; j++) {
            sendbuf[i * cols_per_chunk + j] = local_chunk[i][j];
        }
    }

    int* recvbuf = NULL;
    if (rank == 0) {
        recvbuf = malloc(size * chunk_size * sizeof(int));
    }

    // Gather all chunks flattened into recvbuf at rank 0
    MPI_Gather(sendbuf, chunk_size, MPI_INT,
               recvbuf, chunk_size, MPI_INT,
               0, cart_comm);

    free(sendbuf);

    if (rank == 0) {
        int dims[2], periods[2], coords[2];
        MPI_Cart_get(cart_comm, 2, dims, periods, coords);

        // Stitch chunks into global_array
        for (int p = 0; p < size; p++) {
            MPI_Cart_coords(cart_comm, p, 2, coords);
            int row_start = coords[0] * rows_per_chunk;
            int col_start = coords[1] * cols_per_chunk;

            for (int i = 0; i < rows_per_chunk; i++) {
                for (int j = 0; j < cols_per_chunk; j++) {
                    global_array[row_start + i][col_start + j] = 
                        recvbuf[p * chunk_size + i * cols_per_chunk + j];
                }
            }
        }
        free(recvbuf);
    }
}

void execute(MPI_Comm cart_comm) {
    //SETUP
    //define halo's (buffers to be sent)
    int* top_halo = (int*)malloc(cols_per_chunk * 2 *  sizeof(int)); //x 2 since we need one space for the int value and one for the coordinate
    int* bottom_halo = (int*)malloc(cols_per_chunk * 2 * sizeof(int));
    int* right_halo = (int*)malloc(rows_per_chunk * 2 * sizeof(int));
    int* left_halo = (int*)malloc(rows_per_chunk * 2 * sizeof(int));

    int halo_used_top, halo_used_bottom, halo_used_left, halo_used_right;

    //Sending request declaration 
    MPI_Request reqs[4] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    //Setup for global convergence check
    bool unstable;
    int total_sent = 0;
    int total_received = 0;
    int global_sent, global_received, global_unstable;
    int local_unstable;

    //Topple setup
    int split;

    //Determining neigbour ranks
    int top_rank, bottom_rank, left_rank, right_rank;
    MPI_Cart_shift(cart_comm, 0, 1, &top_rank, &bottom_rank);
    MPI_Cart_shift(cart_comm, 1, 1, &left_rank, &right_rank);

    //Determining this rank
    int rank;
    MPI_Comm_rank(cart_comm, &rank);

    //Execute the procedure
    do {
        //FLAGS
        unstable = false;
        halo_used_top = 0;
        halo_used_bottom = 0;
        halo_used_left = 0;
        halo_used_right = 0;

        //CHECK FOR MESSAGES and APPLY 
        // Receive halos from neighbors - ONLY IF there is a neighbour (ie: chunks on the border won't receive in the direction of the borders)

        if (top_rank != MPI_PROC_NULL)
            total_received += probe_and_receive(rank, top_rank, 't', local_chunk, rows_per_chunk, cols_per_chunk, cart_comm);

        if (bottom_rank != MPI_PROC_NULL)
            total_received += probe_and_receive(rank, bottom_rank, 'b', local_chunk, rows_per_chunk, cols_per_chunk, cart_comm);

        if (left_rank != MPI_PROC_NULL)
            total_received += probe_and_receive(rank, left_rank, 'l', local_chunk, rows_per_chunk, cols_per_chunk, cart_comm);

        if (right_rank != MPI_PROC_NULL)
            total_received += probe_and_receive(rank, right_rank, 'r', local_chunk, rows_per_chunk, cols_per_chunk, cart_comm);

        //loop through and topple (only if you need to)
        for (int i = 0; i < rows_per_chunk; i++) {
            for (int j = 0; j < cols_per_chunk; j++) {
                if (local_chunk[i][j] >= 4) {
                    unstable = true;

                    split = local_chunk[i][j] / 4;
                    local_chunk[i][j] -= split * 4;

                    //check toppling on the border of the cell
                    // Top
                    if (i == 0) { //checks if we're at the top edge of chunk
                        if (top_rank != MPI_PROC_NULL) { //checks if there even is a neight at the top (if this chunk is not on the border of the global grid)
                            top_halo[halo_used_top++] = j; //adds the coord
                            top_halo[halo_used_top++] = split; //adds to the halo instead of altering on the chunk
                        }
                    } else {
                        local_chunk[i - 1][j] += split; //else we're not on the edge of a chunk, then we just add to the neighbour on the chunk
                    }

                    // Bottom
                    if (i == rows_per_chunk - 1) {
                        if (bottom_rank != MPI_PROC_NULL) {
                            bottom_halo[halo_used_bottom++] = j;
                            bottom_halo[halo_used_bottom++] = split;
                        }
                    } else {
                        local_chunk[i + 1][j] += split;
                    }

                    // Left
                    if (j == 0) {
                        if (left_rank != MPI_PROC_NULL) {
                            left_halo[halo_used_left++] = i;
                            left_halo[halo_used_left++] = split;
                        }
                    } else {
                        local_chunk[i][j - 1] += split;
                    }

                    // Right
                    if (j == cols_per_chunk - 1) {
                        if (right_rank != MPI_PROC_NULL) {
                            right_halo[halo_used_right++] = i;
                            right_halo[halo_used_right++] = split;
                        }
                    } else {
                        local_chunk[i][j + 1] += split;
                    }
                }
            }
        }
        
        //SENDING HALOs WHERE APPLICABLE 
        //checking that prev sent has happened - waits if not (in a specific direction)
        for (int i = 0; i < 4; i++) {
            if (reqs[i] != MPI_REQUEST_NULL) {
                MPI_Wait(&reqs[i], MPI_STATUS_IGNORE);
                reqs[i] = MPI_REQUEST_NULL;
            }
        }

        //sending
        //check that the halo to be sent will actually go somewhere (this chunk is not on the border of anything)
        // also check that there is actually something to be sent (used_halo > 0)
        // Send top 
        if (top_rank != MPI_PROC_NULL && halo_used_top > 0) {
            MPI_Isend(top_halo, halo_used_top, MPI_INT, top_rank, 0, cart_comm, &reqs[0]);
            total_sent++;
        }

        // Send bottom
        if (bottom_rank != MPI_PROC_NULL && halo_used_bottom > 0) {
            MPI_Isend(bottom_halo, halo_used_bottom, MPI_INT, bottom_rank, 0, cart_comm, &reqs[1]);
            total_sent++;
        } 

        // Send left
        if (left_rank != MPI_PROC_NULL && halo_used_left > 0) {
            MPI_Isend(left_halo, halo_used_left, MPI_INT, left_rank, 0, cart_comm, &reqs[2]);
            total_sent++;
        }

        // Send right
        if (right_rank != MPI_PROC_NULL && halo_used_right > 0) {
            MPI_Isend(right_halo, halo_used_right, MPI_INT, right_rank, 0, cart_comm, &reqs[3]);
            total_sent++;
        }

        //CONVERGENCE CHECK
        //check that there are no unstable chunks across processors 
        //and check that all the sent messages across the entire simulation have been received (total_sent == total_received)

        local_unstable = unstable ? 1 : 0;

        MPI_Allreduce(&local_unstable, &global_unstable, 1, MPI_INT, MPI_LOR, cart_comm);
        MPI_Allreduce(&total_sent, &global_sent, 1, MPI_INT, MPI_SUM, cart_comm);
        MPI_Allreduce(&total_received, &global_received, 1, MPI_INT, MPI_SUM, cart_comm);

        // printf("total sends: %d, total received: %d, global_unstable: %d \n",total_sent, total_received, global_unstable);

    } while(global_unstable || (global_sent != global_received));

    //gather local arrays and stitch together
}

int main(int argc, char** argv) {

    //take in arguemnts
    if (argc < 4) {
        printf("Usage: %s <ROWS> <COLS> <INIT_VALUE>\n", argv[0]);
        return 1;
    }

    ROWS = atoi(argv[1]);
    COLS = atoi(argv[2]);
    INIT_VALUE = atoi(argv[3]);

    //MPI initialisation
    int rank, size;
    MPI_Init(&argc, &argv); 

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    //Create and setup the array - done by MASTER
    int** array = NULL;
    if (rank == 0) {
        array = allocate_2d_array(ROWS, COLS);
        setup(array);

        //test
        array[ROWS / 2][COLS / 2] = 526338;
    }
   
    // START TIMER HERE - includes the time to split the array, execute compute and stitch back together again
    double start_time = MPI_Wtime();

    //Split the array among processors and obtain the cartesian communicator with topology info
    MPI_Comm cart_comm = partition_array(array, rank, size);

    //have each thread execute on their local copy
    execute(cart_comm);

    //Gather local chunks and stitch them together
    if (rank == 0) {
        int** global_array = allocate_2d_array(ROWS, COLS);
        gather_and_stitch(global_array, cart_comm, rank, size);

        // END TIMER HERE
        double end_time = MPI_Wtime();
        double elapsed = end_time - start_time;
        printf("%.6f seconds\n", elapsed);

        // Now global_array contains the full grid; print or save as needed
        //print_2d_array_to_file(global_array, ROWS, COLS);

        //save_as_ppm("output.ppm", global_array);   

        // Free global_array when done
        for (int i = 0; i < ROWS; i++) {
            free(global_array[i]);
        }
        free(global_array);
    }
    else {
        // other ranks just call gather_and_stitch to send their chunk
        gather_and_stitch(NULL, cart_comm, rank, size);
    }
    
    //Free relevant memory
    if (rank == 0 && array != NULL) {
        for (int i = 0; i < ROWS; i++) free(array[i]);
        free(array);
    }

    MPI_Finalize();
    return 0;
}