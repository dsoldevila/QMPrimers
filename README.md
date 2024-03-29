# QMPrimers
A tool to select the best primers to perform quantitative metabarcoding studies. It aims to be used by all kind of users as it comes with a graphical interface. It has two main stages:
* **MATCHING:** Given a genome file and a primer pairs file, it matches every genome sequence with every primer pair and generates a template containing the best matches found as well as additional information. A part from the template, it also generates a table containing the genomes sequences and primer pairs that didn't match and another file containing statistics. If the output directory was "output", the program would generate respectively output_positive.csv, output_negative.csv and output_stats.txt.

* **SIMULATION:** Given the template (output_positive.csv) generated by the MATCHING (or any other program as long as it has the required format), performs a simulated PCR-amplification (or _in silico_ PCR) of a random sample of genomes (or species); the effciency of the PCR of each species in the mixture is mediated by its number of template-primer mismatches. The rationale for the model is taken from [Piñol et al. (2019)](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14776). In short:
    1. We generate a random sample of S species from the pool P of species.
    2. We generate the initial (before PCR) abundance of each of the S species in the sample using a geometric distribution of parameter k (equation 3 of Piñol et al., 2019). 
    3. We estimate the relative amplification efficiency of each species in the sample using the parameter beta and the number of primer-template mismatches of each species according equation 4 of Piñol et al. (2019).
    4. We calculate the relative final (after PCR) DNA concentration of each species in the mixture using equation 2 of Piñol et al. (2019).
    5. Finally, we calculate the Pearson correlation coefficient between the initial and the final relative abundance of species (pre and post-PCR). 
    6. The above steps are repeated N times with a fresh random set of species.

    P is the numer of total genomes for which there was a positive matching in the template file (output_positive.csv).
S, k, beta, and N are given by the user.

    The results are summarized in two files. (1) raw results containing a table with the correlation coefficient for each sample and each primer pair; (2) a summary containing the basic statistics of the simulation for each primer pair.

**WARNING:** This program is still in development, so expect bugs and crashes. Its use is still not recommended!

# Installation
QMPrimers can e executed as a Python script or as a one-file executable.

## As a Python script
To run it as a Python script it is mandatory to use Python3.x, as well as to install the following dependencies.

**Dependencies:**
- [Biopython](https://biopython.org)
- [Numpy](http://www.numpy.org)
- [Pandas](https://pandas.pydata.org)

Once fullfilled the requirements, [execute the program](#run-as-a-python-script).

## As an executable
Download the program via one of the links below, dependiing on your platform:
- [Windows 10](https://github.com/dsoldevila/QMPrimers/releases/download/v0.9.0-beta/QMPrimers.exe)
- insert link to Linux executable
- insert link to MacOS executable

Once downloaded the program, [execute it](#run-as-an-executable).

# Wiki
## First steps

### Run as a Python script
- To execute the program in GUI mode, just type in a terminal:
```bash
python QMPrimers.py
#or python3 or python3.x, where X depends on your python version
```

- To execute the program in command line mode, type in a terminal:
```bash
python QMPrimers.py <paramters....>
```
### Run as an executable
- To execute the program (in GUI mode) double click the program file.

- To execute the program in command line mode, type in a terminal:
```bash
./QMPrimers --help
```
## How to use the program
### Matching
Everything related to the Matching is explained [here](https://github.com/masenar/QMPrimers/wiki/Matching).

### Simulation
Everything related to the Simulation is explained [here](https://github.com/masenar/QMPrimers/wiki/Simulation).

## Development
Feel free to colaborate and to send me a message for any doubts. If you have never contributed to any respository, I recommend to take a look on this [website](https://www.firsttimersonly.com/)

To know what is going on under the hood, you can check this yet-unfinished [page](https://github.com/dsoldevila/QMPrimers/wiki/Implementation)


