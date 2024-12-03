# dna2bit
dna2bit: an ultra-fast and accurate genomic distance estimation software
## **Compilation**
dna2bit is a software tool developed in C++11, leveraging the capabilities of OpenMP for parallel computing and the popcount technique for efficient bit manipulation. It has been thoroughly tested using the g++ and clang compilers on both Linux and MacOS platforms.

To install from source code, please clone this repository first.

```
git clone https://github.com/lijuzeng/dna2bit.git
cd dna2bit/src
```

For Linux users, we recommend instructions
```
cd dna2bit/src
make
```
Detailed instruction in the makefile is
```
g++ dna2bit.cpp -o dna2bit -O3 -Wall -march=native -fopenmp -lz
```

For OSX users, we recommend installing `llvm` and `lib` via Homebrew before compiling dna2bit to ensure compatibility and utilize the latest compiler features.
```
brew install llvm libomp
g++ dna2bit.cpp -o dna2bit -O3 -Wall -march=native -Xpreprocessor -fopenmp -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include -lomp -lz
```
## **Usage**
`dna2bit <command> [option] [arguments ...]`  
  
dna2bit is equipped with two primary commands designed for distinct tasks in genomic analysis:   

`sketch`: This command will generate sketches from input sequence files in the disk.  
`dist`: This command estimates the genomic distance between sketched files and will generate a distance matrix.  

### **sketch**
For single input genome file please try: 
`dna2bit sketch genome1.fna`  
This command will create a sketched bit file named `genome1.fna.k.17.l.8192.bit`, where the suffixes `.k.17` and `.l.8192` represent the default settings for k-mer size and bit array length, respectively.

For a list of genome files please try:
`dna2bit sketch -n 10 -F files.txt`  
Here, the `-n` option is used to specify the number of threads which accelerate the processing. The `-F` option directs the software to read from a file where each line contains the path to a genome file.

### **dist**
For two sketched files to calculate distance metric, please try:
`dna2bit dist -p genome1.fna.k.17.l.8192.bit genome2.fna.k.17.l.8192.bit`  

For a list of sketched files please try:
`dna2bit dist -n 10 -F bits.txt`  

The `-n` and `-F` options are used in the same way for both sketch and dist command.
For detailed usage instructions specific to each command, please use the help option like `dna2bit <command> --help`.
