#include"atomicpp.h"
Atomic_Structure est;

int main(int argc, char **argv)
{
// void Atomic_Structure::print_VASP(string outputfile,string Titulo, float Factor, double M[3][3])
double M[3][3]={{25.0, 0.0 , 0.0},{0.0, 25.0 , 0.0},{0.0, 0.0 , 25.0}};
est.read_xyz(argv[1]);
if(argc<3)
{
  est.print_VASP("tostd","XYZ to VASP", 1.0, M);
  system("cat tostd ; rm tostd");
}
else if(argc==3)
{
  est.print_xyz(argv[2]);
}
return 0;
}
