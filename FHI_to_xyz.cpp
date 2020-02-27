#include"atomicpp.h"
Atomic_Structure est;

int main(int argc, char **argv)
{
est.read_fhi(argv[1]);
if(argc<3)
{
  est.print_xyz("tostd");
  system("cat tostd ; rm tostd");
}
else if(argc==3)
{
  est.print_xyz(argv[2]);
}
return 0;
}
