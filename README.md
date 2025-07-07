# Parallel Abelian Sandpile Simulation

This repository contains implementations of the **Abelian Sandpile Model** using both **serial** and **MPI-based parallel** approaches. These were developed as part of a high-performance computing project to study scalability and simulation efficiency.

---

##  Contents

- `SERIAL/` – Serial CPU implementation  
- `MPI/` – Parallel implementation using MPI  

---

##  Usage
SERIAL:
    Build the program by running:
        make

    Run a default test by running:
        make run

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

##  Dependencies

`open-mpi` — MPI-enabled C compiler



