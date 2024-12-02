# dna2bit
dna2bit: an ultra-fast and accurate genomic distance estimation software
## **Compilation**
```
git clone https://github.com/lijuzeng/dna2bit.git
cd dna2bit/src
make
./dna2bit --help
```
## **Usage**
`dna2bit <command> [option] [arguments ...]`  
  
there are two commands for dna2bit: sketch and dist.   

sketch: create sketch of input sequence file  
dist: estimate distance between sketched files  
  
in each command, use `dna2bit <command> --help` to see command-specific usage
### **sketch**
for one input genome file  
`dna2bit sketch genome1.fna`  
this command line will generate a sketched bit file 'genome1.fna.k.17.l.8192.bit'  
the number behind k. and l. are the default setting (k-mer and bit_array_length) of dna2bit

for a list of genome files  
`dna2bit sketch -n 10 -F files.txt`  
-n 10: 10 threads, -F files.txt: the each line of 'files.txt' is the path to genome file  

### **dist**
for two sketched file to calculate distance metric
`dna2bit dist -p genome1.fna.k.17.l.8192.bit genome2.fna.k.17.l.8192.bit`  

for a list of sketched files  
`dna2bit dist -n 10 -F bits.txt`  
-n 10: 10 threads, -F bits.txt: the each line of 'bits.txt' is the path to sketched file
