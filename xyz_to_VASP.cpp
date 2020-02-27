#include"atomicpp.h"
Atomic_Structure est;

int main(int argc, char **argv)
{
est.read_xyz(argv[1]);
if(argc<3)
{
  est.print_VASP("tostd");
  system("cat tostd ; rm tostd");
}
else if(argc==3)
{
  est.print_xyz(argv[2]);
}
return 0;
}
