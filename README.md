# VASP-Tools
Useful auxiliary programs for  VASP calculations

# Converter POSCAR (VASP) format to .xyz

C code. Just copy (or download) the source code and compile it as usual. The program takes the first argument for input and the second is the file where you want to direct the output. The program recognize if the input file have Direct or Cartesian format and convert fractional coordinates to Cartesian (legible .xyz format) if it is necessary. It also works if  the file has selective dynamics.

Compile: 

~$ gcc -o POSCARtoXYZ POSCARtoXYZ.c -lm 

Run:

~$ ./POSCARtoXYZ File_to_convert(InputFile) File_converted(OutputFile).xyz

# Converter  .xyz to POSCAR format

Shell script. .xyz don't has information about the cell, so the program don't  write the Vectors nor the scale factor but it carefully writes the appropriate coordinates with the corresponding atomic symbol. The program also recognize and transform selective dynamics. Copy or download the source code. Save it for example in XYZtoPOSCAR.sh. 

Change file  permissions to executable mode.

~$ chmod +x XYZtoPOSCAR.sh

Run:

~$ ./XYZtoPOSCAR.sh File_to_convert(InputFile).xyz File_converted(OutputFile)



