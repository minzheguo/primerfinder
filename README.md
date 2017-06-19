# primerfinder

This program integrates primer3 and OligoCalc functions to pick candidate primers for staining.

## Workflow

* Take RNA sequence as input
* Call Primer3 to find hybridization probe (internal oligo) candidates with 21 bases with melting temperature 45~65 (celcius)
* Pick candidates ending with “TA” or “AG” 
* Generate DNA primer pairs (20+20) for valid candidates
* Use the OligoCalc algorithm to check the melting temperature and self-complementarity of primer pairs 
* Output valid primer pairs

## Installation and usage

* Download the entire "primerfinder" project to local machine. Unzip the "primer3.zip" inside the project folder.
* Open Windows command line
* Go to the “primerfinder” folder
* Save nucleotide sequence in a file in the input folder, such as “.\input\SCGB1A1.txt”
* Make sure Rscript.exe is in your system path
* Run the program by typing “Rscript primerfinder.R .\input\SCGB1A1.txt”
* The output will be saved in a “uidxxx” folder in output directory containing the following files
** primer3.in.txt: the input file to primer3
** primers.1.all.txt: all primer candidates picked by primer3
** primers.2.validending.txt: all primer candidates with valid ending
** primers.3.pairs.txt: all primer pairs
** primers.4.final.txt: valid primer pairs

## Notes

* Currently support Windows system.
* Currently only accept bases: A, C, T, G, U.
* The program is in the testing phase. Highly recommend double checking the properties of the final primers in OligoCalc (http://biotools.nubic.northwestern.edu/OligoCalc.html).
* The implementation of the self-complementarity checking algorithm is more stringent than the implementation in OligoCalc. 
* The Windows executables of primer3 (libprimer3 release 2.3.6) were downloaded from http://primer3.sourceforge.net/. 
* The primer3 execution results are consistent with this web version of primer3: http://primer3.ut.ee/.



