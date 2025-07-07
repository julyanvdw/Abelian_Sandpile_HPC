SERIAL:
    Build the program by running:
        make

    Run a default test by running:
        make run

    Clean build artifacts by running:
        make clean

OMP:
    The OMP solution only simluates the middle toppling (INIT_VALUE set to 0 and middle set to 526338)
        Everywhere toppling (INIT_VALUE = 4) can be simulated by editing the OMP_abelian_sandpile.c:
            1. INIT_VALUE = 4 (on line 6)
            2. Removing line 125

    Build the program by running:
        make
        This compiles with -O3 and -fopenmp, producing the executable omp_sandpile

    Run a default test by running:
        make run
        This sets OMP_NUM_THREADS=4 and runs: ./omp_sandpile 1024 1024 4
            NB: Note that the last number of the argument is for the number of thread, not INIT_VALUE
    Custom runs:
        To change the run target, edit the 'run' line in the Makefile.
    
    Or run manually:
        Manually run ./omp_sandpile <rows> <columns> <num threads>

    Clean build artifacts by running:
    make clean

MPI:
    Build the program by running:
        make all

    Run a default test by running:
        make run
        This sets number of processes = 8 and runs: mpirun -np 8 ./MPI_abelian_sandpile 1024 1024 526338

    Custom runs:
        To change the run target, edit the 'run' line in the Makefile.
    
    Or run manually:
        Manually run:
             mpirun -np <num processes> ./MPI_abelian_sandpile <rows> <columns> <INIT_VALUE>

    Clean build artifacts by running:
        make clean

Cluster:
    1. Change directories by:
        cd /mnt/lustre/users/<userID>
    2. Compile the code using make
    3. Create your batchscript
    4. Submit a job by:
        qsub <batchscript>.pbs