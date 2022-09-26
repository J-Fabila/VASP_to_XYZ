# VASP-Tools

Useful auxiliary programs for  VASP calculations


### Converter POSCAR (VASP) format to .xyz

C code. Just copy (or download) the source code and compile it as usual.

##### Compile: 
```
~$ gcc -o to_xyz POSCARtoXYZ.c -lm 
```
The program takes the first argument for input and the second is the file where you want to direct the output. If you  skip the second argument the output will be shown in display. The program recognize if the input file has Direct or Cartesian format and convert fractional coordinates to Cartesian (legible .xyz format) if it is necessary. It also works if  the file has selective dynamics.

##### Run:
```
~$ ./POSCARtoXYZ File_to_convert(InputFile) File_converted(OutputFile).xyz
```
You can put the program in some directory and add it to the environment variable PATH in your .bashrc, so you can use the executable like any other bash command. For example, in your /home/$USER directory you can create the folder 'Programs', in there compile the code. Now open with any text editor (like Vi or Nano) the .bashrc file (it is in your /home/$USER/ directory). At the end of that file add the line
```
~$ PATH=$PATH:/home/$USER/Programs
```
Save the file and execute
```
~$ source .bashrc
```
After this you can run

##### Run: 
```
~$ to_xyz File_to_convert(InputFile) File_converted(OutputFile).xyz
```


### Converter  .xyz to POSCAR format

Shell script. The .xyz file has no information about the cell, so the program don't  write the Vectors nor the scale factor but it carefully writes the appropriate coordinates with the corresponding atomic symbol. The program also recognize and transform selective dynamics. Copy or download the source code. Save it for example in to_vasp 

##### Change file  permissions to executable mode:
```
~$ chmod u+x to_vasp
```
##### Run:
```
~$ ./to_vasp File_to_convert(InputFile).xyz File_converted(OutputFile)
```

You can skip the second argument and it will shown the output in display. You can also add this program to /home/$USER/Programs and use it like any other bash command.
