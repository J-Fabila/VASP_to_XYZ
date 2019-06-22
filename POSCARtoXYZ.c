#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
//#include"/home/jorge/Proyectitos/Bibliotecas/atomic.h"

struct Atom
{
   char Symbol[3];
   float x[3];
};




/********************* Function POSCARtoXYZ  ************************/

void POSCARtoXYZ(char *inputfile, char *outputfile)
{
struct Atom xyz[500];
struct Atom poscar[500];

int i, j, l, m;
char Symbol[500][3];
char command[150];
int N_Symbol[500];
float Factor;
float M[3][3];
float Mi[3][3];
float suma;
int Nat;
float s;
int Ntyp;
int Cartesian;
int Sel;
int Direct;
/*************** Selectivedynamics **********************/
strcpy(command,"cat ");
strcat(command, inputfile);
strcat(command,"  >>  aux");
system(command);
system("grep \"elective\" aux | wc -l >>sel ; rm aux ");

FILE *po=fopen("sel", "r");
fscanf(po,"%i\n",&Sel);
fclose(po);
system("rm sel");
strcpy(command, "");

if(Sel==1)//Si hay selective dynamics)
{
strcpy(command, "tail -$(($(grep -A 500 \"elective\" ");
strcat(command, inputfile);
strcat(command, " | wc -l )-2)) ");
strcat(command, inputfile);
strcat(command, " | grep .  >> selective  ");
FILE *tg=fopen("commandi","w");
fprintf(tg,"%s\n",command);
fclose(tg);
system("chmod +x commandi ");
strcpy(command,"");
strcpy(command," ./commandi  ");
system(command);
strcpy(command, "");
strcpy(command,"Nat=$(cat selective | wc -l  ) " );
system(command);
strcpy(command,"");
//strcat(command,"for((i=0;i<$(($Nat+1));i++));"); //SI da problemas el for quitarlo, igual funcioinara
system(" awk '{print $4 \" \" $5 \" \" $6}' selective | tr 'F' '0' | tr  'T' '1' >>selectivedynamicsaux ");
system( " echo \" \" >> selectivedynamics ; echo \" \" >>selectivedynamics ; cat selectivedynamicsaux >> selectivedynamics");
system("rm selective selectivedynamicsaux");
}
else
{

}
/************************************************/

strcpy(command, "grep \"irect\" ");
strcat(command,inputfile);
strcat(command," | wc -l >> Direct");
system(command);
FILE *x=fopen("Direct", "r");
fscanf(x,"%i\n",&Direct);
fclose(x);
system("rm Direct");
strcpy(command, "");

strcpy(command, "grep \"artesian\" ");
strcat(command,inputfile);
strcat(command," | wc -l >> Cartesian");
system(command);
FILE *xi=fopen("Cartesian", "r");
fscanf(xi,"%i\n",&Cartesian);
fclose(xi);
system("rm Cartesian");

if(Direct==1)
{
strcpy(command, "head -5 ");
strcat(command, inputfile);
strcat(command, " | tail -3 >> Matriz");
system(command);
strcpy(command, "");
strcpy(command, "head -2 ");
strcat(command, inputfile);
strcat(command, " | tail -1 >> Factor");
system(command);
strcpy(command, "");
strcpy(command, "tail -$(($(grep -A 500 \"irect\" ");
strcat(command, inputfile);
strcat(command, " | wc -l )-1)) ");
strcat(command, inputfile);
strcat(command, " | awk '{ print $1 \"  \" $2 \"  \" $3  }' | grep . >> Posiciones ");
//strcat(command, " >> Posiciones  ");
system(command);
strcpy(command, "");
strcpy(command, "head -6 ");
strcat(command, inputfile);
strcat(command, "  | tail -1 >> aux");
system(command);
strcpy(command, "");
strcpy(command, "head -7 ");
strcat(command, inputfile);
strcat(command, "  | tail -1 >> aux2");
system(command);
strcpy(command, "");
strcpy(command, "tr -s '[:blank:]' '\n' < aux >> Simbolos1");
system(command);
strcpy(command, "");
strcpy(command, "tr -s '[:blank:]' '\n' < aux2 >> Numeros");
system(command);
strcpy(command, "");
system("grep  . Simbolos1 >> Simbolos");
system("rm Simbolos1");
system("cat Simbolos | wc -l >>  Ntyp");//Guarda el numero de especies como archivo

FILE *f= fopen("Matriz", "r");
FILE *g= fopen("Factor", "r");
FILE *h= fopen("Posiciones", "r");
FILE *q= fopen("Simbolos", "r");
FILE *r= fopen("Numeros", "r");
FILE *p= fopen("Ntyp", "r");
FILE *o= fopen(outputfile, "w");
fscanf(p,"%i \n", &Ntyp);
for(i=0;i<Ntyp;i++)
{
fscanf(q,"%2s \n", Symbol[i]);
fscanf(r,"%i \n", &N_Symbol[i]);
suma=suma+N_Symbol[i];
}
Nat=suma;
l=0;
for(i=0;i<Ntyp;i++) //Reconoce que tipo de atomo es cada vector posicion
{
for(j=0;j<N_Symbol[i];j++)
{
strcpy( xyz[l].Symbol, Symbol[i] );
l=l+1;
}
}

for(i=0;i<3;i++) //Lee los elementos de matriz
{
fscanf(f,"%f %f %f \n", &Mi[i][0], &Mi[i][1], &Mi[i][2]);
}
/*
for(i=0;i<3;i++)
{
for(j=0;j<3;j++)
{
M[i][j]=Mi[i][j];
}
}
*/
fscanf(g,"%f \n", &Factor);  //Lee el factor de escala
for(i=0;i<3;i++)//Multiplica el factor de escala a la matriz
{
for(j=0;j<3;j++)
{
M[i][j]=Factor*M[j][i];
}
}
for(i=0;i<Nat;i++) //Extrae las coordenadas de Positions
{
fscanf(h,"%f %f %f \n", &poscar[i].x[0], &poscar[i].x[1], &poscar[i].x[2]);
}
/*
for(i=0;i<Nat;i++)//Aplica la matriz  de cambio a cada atomo
{
for(l=0;l<3;l++)
    {
     s=0;
     for(m=0;m<3;m++)
        {
        s=s+(M[l][m]*(poscar[i].x[m]));
        }
     xyz[i].x[l]=s;
    }
}*/
for(i=0;i<Nat;i++)//Aplica la matriz  de cambio a cada atomo
{
for(l=0;l<3;l++)
{
for(m=0;m<3;m++)
{
M[l][m]=poscar[i].x[l]*Mi[l][m];
//printf("%f\t", M[l][m]);
}
//printf(" \n");
}
//printf(" \n");
for(m=0;m<3;m++)
    {
     s=0;
     for(l=0;l<3;l++)
        {
        s=s+M[l][m];
        }
     xyz[i].x[m]=s;
    }
}
system("rm Matriz Factor Posiciones Ntyp Simbolos Numeros aux aux2"); //Borra archivos residuales
fprintf(o,"%i\n\n", Nat); //Imprime el numero de atomos como primera linea del xyz
for(i=0;i<Nat;i++) //Imprime las lineas de atomos
{
fprintf(o,"%s %f %f %f\n", xyz[i].Symbol, xyz[i].x[0], xyz[i].x[1], xyz[i].x[2]);

}
fclose(o);

}//Cierra if/********************************************************/
else{
if(Cartesian==1)
{
strcpy(command, "");
strcpy(command, "tail -$(($(grep -A 500 \"artesian\" ");
strcat(command, inputfile);
strcat(command, " | wc -l )-1)) ");
strcat(command, inputfile);
strcat(command, " | awk '{ print $1 \"  \" $2 \"  \" $3  }' | grep . >> Posiciones ");
//strcat(command, " >> Posiciones  ");
system(command);
system(command);
strcpy(command, "");
strcpy(command, "head -6 ");
strcat(command, inputfile);
strcat(command, "  | tail -1 >> aux");
system(command);
strcpy(command, "");
strcpy(command, "head -7 ");
strcat(command, inputfile);
strcat(command, "  | tail -1 >> aux2");
system(command);
strcpy(command, "");
strcpy(command, "tr -s '[:blank:]' '\n' < aux >> Simbolos1");
system(command);
strcpy(command, "");
strcpy(command, "tr -s '[:blank:]' '\n' < aux2 >> Numeros");
system(command);
strcpy(command, "");
system("grep  . Simbolos1 >> Simbolos");
system("rm Simbolos1");
system("cat Simbolos | wc -l >>  Ntyp");//Guarda el numero de especies como archivo

FILE *hi= fopen("Posiciones", "r");
FILE *qi= fopen("Simbolos", "r");
FILE *ri= fopen("Numeros", "r");
FILE *pi= fopen("Ntyp", "r");
FILE *oi= fopen(outputfile, "w");
fscanf(pi,"%i \n", &Ntyp);
for(i=0;i<Ntyp;i++)
{
fscanf(qi,"%2s \n", Symbol[i]);
fscanf(ri,"%i \n", &N_Symbol[i]);
suma=suma+N_Symbol[i];
}
Nat=suma;
l=0;
for(i=0;i<Ntyp;i++) //Reconoce que tipo de atomo es cada vector posicion
{
for(j=0;j<N_Symbol[i];j++)
{
strcpy( xyz[l].Symbol, Symbol[i] );
l=l+1;
}
}
for(i=0;i<Nat;i++) //Extrae las coordenadas de Positions
{
fscanf(hi,"%f %f %f \n", &poscar[i].x[0], &poscar[i].x[1], &poscar[i].x[2]);
}
system("rm  Posiciones Ntyp Simbolos Numeros aux aux2"); //Borra archivos residuales
fprintf(oi,"%i\n\n", Nat); //Imprime el numero de atomos como primera linea del xyz
for(i=0;i<Nat;i++) //Imprime las lineas de atomos
{
fprintf(oi,"%s %f %f %f\n", xyz[i].Symbol, poscar[i].x[0], poscar[i].x[1], poscar[i].x[2]);

}
fclose(oi);

}//Cierra if cartesian
}//Cierra else

if (Sel==1)
{
strcpy(command,"cat ");
strcat(command,outputfile);
strcat(command," >> aux ;  rm ");
strcat(command,outputfile);
strcat(command," ; paste aux selectivedynamics  >> ");
strcat(command,outputfile);
system(command);
system("rm  aux selectivedynamics");
}
}//Cierra funcion


int main(int argc, char **argv)
{

if(argc<3)
{
POSCARtoXYZ(argv[1],"tostd");
system("cat tostd ; rm tostd");
}
else if(argc==3)
{
POSCARtoXYZ(argv[1],argv[2]);
}
return 0;
}


