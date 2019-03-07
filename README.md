# QMPrimers
A tool to select the best primers to perform quantitative metabarcoding studies. It aims to be used by all kind of users as it comes with a graphical interface.

**WARNING:** This program is still in development, so expect bugs and crashes. Its use is still not recommended!

### Dependencies

- [Biopython](https://biopython.org)
- [Numpy](http://www.numpy.org)
- [Pandas](https://pandas.pydata.org)

### Brief tutorial

This programm has a command line mode and a GUI mode:

* Command line mode:
  ```bash
  python3.7 QMPrimers --nogui <other required parameters>
  
  #Use --help to print a brief manual 
  python3.7 QMPrimers --help
  ```
  
* GUI mode:
  ```bash
  python3.7 QMPrimers <no options>
  ```
  ![screenshot](./Misc/screenshot.png)

  * **Hanging primers**: mf = forward maximum missmatches, mr = reverse max. miss. 
  Primer pairs are allowed to match between [0-mf,len(genome)+mr] instead of just between the length of the genome
  * The **primer pairs file** must have the following header. The order does not matter:                        id;forwardPrimer;fPDNA;reversePrimer;rPDNA;ampliconMinLength;ampiconMaxLength




