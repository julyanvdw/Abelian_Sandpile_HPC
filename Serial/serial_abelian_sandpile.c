// Julyan van der Westhuizen (VWSJUL003) 
// HPC Assignment - Serial Version with Dynamic Allocation
// 18/05/25

//first we create an unoptimised version of the serial code (using brute force )

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

int ROWS;
int COLS;
int INIT_VALUE;

//Testing function
void print_2d_array(int** array, int rows, int cols) {
    FILE* file = fopen("printed.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
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
// Function to allocate 2D array dynamically
int** allocate_2d_array(int rows, int cols) {
    int** array = (int**)malloc(rows * sizeof(int*));
    
    for (int i = 0; i < rows; i++) {
        array[i] = (int*)malloc(cols * sizeof(int));
    }
    
    return array;
}

// Function to free 2D array
void free_2d_array(int** array, int rows) {
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

void setup(int** array) {
    
    // Initialise the array with zeros (including sink padding)
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            array[i][j] = INIT_VALUE;
        }
    }
    
}

void execute(int** array) {
    bool unstable;
    int split = 0;
    
    do {
        unstable = false;
        
        for (int i = 0; i < ROWS; i++) {
            for (int j = 0; j < COLS; j++) {
                if (array[i][j] > 3) {
                    unstable = true;
                    
                    //topple the cell
                    split = array[i][j] / 4;
                    array[i][j] -= 4 * split;

                    //check if the toppling is on the grid border
                    if (i > 0) {
                        array[i-1][j] += split;
                    }

                    if (i < ROWS - 1) {
                        array[i+1][j] += split;
                    }
                    
                    if (j > 0) {
                        array[i][j-1] += split;
                    }

                    if (j < COLS - 1) {
                        array[i][j+1] += split;
                    }
                }
            }
        }
    } while (unstable);

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

int main(int argc, char *argv[]) {  

    //take in arguemnts
    if (argc < 4) {
        printf("Usage: %s <ROWS> <COLS> <INIT_VALUE>\n", argv[0]);
        return 1;
    }

    ROWS = atoi(argv[1]);
    COLS = atoi(argv[2]);
    INIT_VALUE = atoi(argv[3]);
    
    int** array = allocate_2d_array(ROWS, COLS);
    setup(array);

    //test
    array[ROWS / 2][COLS / 2] = 526338;

    clock_t start = clock();  
    execute(array);
    clock_t end = clock();

    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%.6f\n", elapsed);

    //save_as_ppm("output.ppm", array);
    
    //print_2d_array(array, ROWS, COLS);
    
    free_2d_array(array, ROWS);

    return 0;
}
