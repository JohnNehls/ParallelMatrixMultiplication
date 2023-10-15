# ParallelMatrixMultiplication
A short study of OpenMP and MPI by way of matrix multiplication, results and context reported in [report.pdf](./report.pdf).  The examples are written in C++ and make use of the Vector class from the standard library.


## Run
See the [makefile](./makefile) for the build commands for the MPI, OpenMP, and serial codes.  For running instructions, execute any of the executables (.exe) without any command line inputs.

## Dependencies
On Fedora/REHL the following commands are needed for dependencies and environment setup:
``` bash
sudo dnf install g++ openmpi-devel
source /etc/profile.d/modules.sh
module load mpi/openmpi-x86_64
```
