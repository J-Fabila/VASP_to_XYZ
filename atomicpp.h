//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/ atomicpp.h library _/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_ For simulation of atoms and molecules _/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

#include<iostream>
#include<stdlib.h>
#include<string.h>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include  <string>
#include   <cmath>
#include   <ctime>


#include <string>
#include <utility>
#include <map>
#include <iomanip>

/********* * * * * * *  *  *   *   *   *  *  *  * * * ***************/

int Nat;
int i,j,k,ll;
using namespace std;


/********************************************************************/
/*********************** Electron Definition ************************/
/********************************************************************/

class Electron
{
   public:
     double x[3];
     double v[3];
     double a[3];
     double M;
     double Z;
     double spin;
     Electron();
     void initialize_Electron(double, double, double, double);
     void print_electron(string);
};

Electron::Electron()
{
   spin=0;
   Z=-1; //-1 e ; −1.602 176 565(35)×10−19 C​
   M=1/1822.8884; // 1822.8884845 (14)−1 uma ; 9,109 382 91(40)×10−31 kg;
   x[0]=0.0; x[1]=0.0; x[2]=0.0;
   v[0]=0.0; v[1]=0.0; v[2]=0.0;
   a[0]=0.0; a[1]=0.0; a[2]=0.0;
}


void Electron::initialize_Electron(double _x, double _y, double _z,double _spin)
{
   x[0] = _x;
   x[1] = _y;
   x[2] = _z;
   spin = _spin;
}

void Electron::print_electron(string file)
{
   ofstream out_f(file);
   out_f<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<spin<<endl;
}
/********************************************************************/
/************************* Atom Definition **************************/
/********************************************************************/

class Atom
{
   public:
     string Symbol;
     double x[3];
     double v[3];
     double a[3];
     double R;
     double M;
     int Z;
     int charge;
     Electron *electron;
     Atom();
     void read_Atom(string, float, float, float);
};

Atom::Atom()
{
   Symbol="AAA"; R=0; Z=0; M=0;
   electron=NULL; charge=0;
   x[0]=0.0; x[1]=0.0; x[2]=0.0;
   v[0]=0.0; v[1]=0.0; v[2]=0.0;
   a[0]=0.0; a[1]=0.0; a[2]=0.0;
}

void Atom::read_Atom(string _Symbol, float _x, float _y, float _z)
{
   Symbol = _Symbol;
   x[0] = _x;
   x[1] = _y;
   x[2] = _z;
}

double Atomic_Distance(Atom atomo1, Atom atomo2)
{
   double suma=0;
   for(ll=0;ll<3;ll++)
   {
      suma=suma+pow(atomo1.x[ll]-atomo2.x[ll],2);
   }
   return sqrt(suma);
}

/********************************************************************/
/******************* Simulation_Cell Definition *********************/
/********************************************************************/

class Simulation_Cell{
    public:
       double M[3][3]={ {25.00,0.00,0.00}, {0.00,25.00,0.00},{0.00,0.00,25.00} };
       bool periodicity=false;
       Simulation_Cell(double M[3][3],bool periodicity);
       Simulation_Cell();
       ~Simulation_Cell();
       void read_VASP(string file,bool periodicity);
       void print_screen();
};

Simulation_Cell default_box;

/********************************************************************/
/******************* Atomic_Structure Definition ********************/
/********************************************************************/

class Atomic_Structure{
     public:
        int Nat;
        int total_charge;
        Atom *atom;
        Atomic_Structure(string);
        Atomic_Structure();
        ~Atomic_Structure();
        void read_xyz(string file);
        void read_fhi(string file);
        void read_VASP(string file);
        void print_xyz(string file, string tag);
        void print_fhi(string file);
        void print_VASP(string, string, float, double (*)[3]);
        void move(double, double, double);
        double x_min(); double x_max();
        double y_min(); double y_max();
        double z_min(); double z_max();
        bool fit_in(float, float, float, float, float, float);
        bool fit_in(Simulation_Cell);
        void show(string visualizer);
        void molecular_dynamic(Simulation_Cell box,  int nsteps, string movie_file, float temperature, double cutoff_radii,  string force_field);
        void geometry_optimization(string movie_file,int max_steps, Simulation_Cell box, double tolerance, double cutoff_radii, int optimization_algorithm, string force_field);
        void initialize_electronic_density(int initialization, int charge);
};

/********************************************************************/
/************************ Cluster Definition ************************/
/********************************************************************/


class Cluster : public Atomic_Structure{
     public:
        string type;
        double quirality;
        Cluster(string);
        Cluster();
        ~Cluster();
        void rotate_Rad(float, float);
        void rotate_Deg(float, float);
        void kick(float);
        void kick_lennard(float);
        void swap(int);
        void srand_generator(string, int, string, int, float);
        void srand_generator_box(string, int, string, int, float);
        void rand_generator(string, int, string, int);
        void centroid();
        void simulated_annealing(float kick_s, string movie_file,int max_steps, double tolerance, double cutoff_radii, int optimization_algorithm, string force_field);
};


/********************************************************************/
/************************ Crystal Definition ************************/
/********************************************************************/


class Crystal : public Atomic_Structure{
     public:
       //Parámetros de celda
       double lattice[3][3];
       double factor=1;
        Crystal(string);
        Crystal();
        ~Crystal();
        void read_fhi(string file);
        void read_VASP(string file);
        void print_fhi(string file);
        void print_VASP(string, string, float);
};

/********************************************************************/
/*********************** Molecule Definition ************************/
/********************************************************************/

class Molecule : public Atomic_Structure{
    public:
       Molecule(string);
       Molecule();
       ~Molecule();
      void rotate_Rad(float, float);
      void rotate_Deg(float, float);
      void centroid();

};


/********************************************************************/
/********************** Emptiness Constructor ***********************/
/******************* Simulation_Cell  super_cell; *******************/
/********************************************************************/

Simulation_Cell::Simulation_Cell(){}


/********************************************************************/
/******************** Constructor from Matrix ***********************/
/************* Simulation_Cell  super_cell(M,false); ****************/
/********************************************************************/

Simulation_Cell::Simulation_Cell(double _M[3][3],bool _periodicity)
{
  M[0][0]=_M[0][0]; M[1][0]=_M[1][0]; M[2][0]=_M[2][0];
  M[0][1]=_M[0][1]; M[1][1]=_M[1][1]; M[2][1]=_M[2][1];
  M[0][2]=_M[0][2]; M[1][2]=_M[1][2]; M[2][2]=_M[2][2];
  periodicity=_periodicity;
}


/********************************************************************/
/********************* Constructor from file ************************/
/************* Simulation_Cell  super_cell(file,false); ****************/
/********************************************************************/

void Simulation_Cell::read_VASP(string file,bool _periodicity)
{
   string command="head -5 ";
          command+=file;
          command+=" | tail -3 > matriz ";
   system(command.c_str());
   ifstream fm("matriz");
   for(i=0;i<3;i++) //Lee los elementos de matriz
   {
      fm>>M[i][0]>> M[i][1]>>M[i][2];
   }
   fm.close();
   system("rm matriz");
  periodicity=_periodicity;
}


/********************************************************************/
/********************** Emptiness Constructor ***********************/
/***************** Atomic_Structure  molecule_name; *****************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Atomic_Structure::Atomic_Structure()
{
  Nat=0;
  total_charge=0;
  atom=NULL;
}

/********************************************************************/
/********************** Emptiness Constructor ***********************/
/********************* Molecule  molecule_name; *********************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Molecule::Molecule() : Atomic_Structure()
{
  Nat=0;
  atom=NULL;
}

/********************************************************************/
/********************** Emptiness Constructor ***********************/
/********************** Cluster  cluster_name; **********************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Cluster::Cluster() : Atomic_Structure()
{
  Nat=0;
  atom=NULL;
}

/********************************************************************/
/********************** Emptiness Constructor ***********************/
/********************** Crystal  crystal_name; **********************/
/************************* Molecule Cysteine; ***********************/
/********************************************************************/

Crystal::Crystal() : Atomic_Structure()
{
  Nat=0;
  atom=NULL;
}

Atomic_Structure::~Atomic_Structure(){}
Molecule::~Molecule(){}
Cluster::~Cluster(){}
Crystal::~Crystal(){}
Simulation_Cell::~Simulation_Cell(){}

/**************************** random_number *************************/
/********************* float r=random_number(0,10); *****************/
/********************************************************************/

float random_number(float min,float max)
{
   float random = min + ((float)rand())/((float)RAND_MAX ) * (max-min);
   return random;
}


/********************************************************************/
/************************ minimum_separation ************************/
/**************** Molecule  molecule_1, molecule_2; *****************/
/******* separation=minimum_separation(molecule_1,molecule_2); ******/
/********************************************************************/

double minimun_separation(Atomic_Structure Molecule1, Atomic_Structure Molecule2)
{
   double  min;
   int tot=Molecule1.Nat*Molecule2.Nat;
   int j;
   int k=0, l=0;
   double distances[tot];
   ofstream distancias("Distancias");

   for(i=0;i<Molecule1.Nat;i++)
   {
      for(j=0;j<Molecule2.Nat;j++)
      {
         distances[k]=Atomic_Distance(Molecule1.atom[i],Molecule2.atom[j])-Molecule1.atom[i].R-Molecule2.atom[j].R;
         distancias<<distances[k]<<endl;
         k++;
      }
   }
   distancias.close();
   system("cat Distancias | sort -n | head -1 >>Minimum_Separation");
   ifstream Minimum("Minimum_Separation");
   Minimum>>min;
   system("rm Minimum_Separation");
  return min;
}

/********************************************************************/
/************************ Radii Constructor *************************/
/******************** map<string, double> Radios; *******************/
/********************* Radios=radii_dictionary(); *******************/
/********************************************************************/

map<string, double> radii_dictionary()
{
  typedef pair<string, double> radio_atomico;

     map<string, double> Radii;
     Radii.insert( radio_atomico("H"  , 0.46 ) );
     Radii.insert( radio_atomico("He" , 1.22 ) );
     Radii.insert( radio_atomico("Li" , 1.57 ) );
     Radii.insert( radio_atomico("Be" , 1.12 ) );
     Radii.insert( radio_atomico("B"  , 0.81 ) );
     Radii.insert( radio_atomico("C"  , 0.77 ) );
     Radii.insert( radio_atomico("N"  , 0.74 ) );
     Radii.insert( radio_atomico("O"  , 0.74 ) );
     Radii.insert( radio_atomico("F"  , 0.72 ) );
     Radii.insert( radio_atomico("Ne" , 1.6  ) );
     Radii.insert( radio_atomico("Na" , 1.91 ) );
     Radii.insert( radio_atomico("Mg" , 1.6  ) );
     Radii.insert( radio_atomico("Al" , 1.43 ) );
     Radii.insert( radio_atomico("Si" , 1.18 ) );
     Radii.insert( radio_atomico("P"  , 1.10 ) );
     Radii.insert( radio_atomico("S"  , 1.04 ) );
     Radii.insert( radio_atomico("Cl" , 0.99 ) );
     Radii.insert( radio_atomico("Ar" , 1.92 ) );
     Radii.insert( radio_atomico("K"  , 2.35 ) );
     Radii.insert( radio_atomico("Ca" , 1.97 ) );
     Radii.insert( radio_atomico("Sc" , 1.64 ) );
     Radii.insert( radio_atomico("Ti" , 1.44 ) );
     Radii.insert( radio_atomico("Va" , 0.8  ) );
     Radii.insert( radio_atomico("Cr" , 1.29 ) );
     Radii.insert( radio_atomico("Mn" , 1.37 ) );
     Radii.insert( radio_atomico("Fe" , 1.26 ) );
     Radii.insert( radio_atomico("Co" , 1.25 ) );
     Radii.insert( radio_atomico("Ni" , 1.25 ) );
     Radii.insert( radio_atomico("Cu" , 1.28 ) );
     Radii.insert( radio_atomico("Zn" , 1.37 ) );
     Radii.insert( radio_atomico("Ga" , 1.53 ) );
     Radii.insert( radio_atomico("Ge" , 1.22 ) );
     Radii.insert( radio_atomico("As" , 1.21 ) );
     Radii.insert( radio_atomico("Se" , 1.04 ) );
     Radii.insert( radio_atomico("Br" , 1.14 ) );
     Radii.insert( radio_atomico("Kr" , 1.98 ) );
     Radii.insert( radio_atomico("Rb" , 2.5  ) );
     Radii.insert( radio_atomico("Sr" , 2.15 ) );
     Radii.insert( radio_atomico("Y"  , 1.82 ) );
     Radii.insert( radio_atomico("Zr" , 1.6  ) );
     Radii.insert( radio_atomico("Nb" , 1.47 ) );
     Radii.insert( radio_atomico("Mo" , 1.4  ) );
     Radii.insert( radio_atomico("Tc" , 1.35 ) );
     Radii.insert( radio_atomico("Ru" , 1.34 ) );
     Radii.insert( radio_atomico("Rh" , 1.34 ) );
     Radii.insert( radio_atomico("Pd" , 1.37 ) );
     Radii.insert( radio_atomico("Ag" , 1.66 ) );
     Radii.insert( radio_atomico("Cd" , 1.52 ) );
     Radii.insert( radio_atomico("In" , 1.67 ) );
     Radii.insert( radio_atomico("Sn" , 1.58 ) );
     Radii.insert( radio_atomico("Sb" , 1.41 ) );
     Radii.insert( radio_atomico("Te" , 1.37 ) );
     Radii.insert( radio_atomico("I"  , 1.33 ) );
     Radii.insert( radio_atomico("Xe" , 2.18 ) );
     Radii.insert( radio_atomico("Cs" , 2.72 ) );
     Radii.insert( radio_atomico("Ba" , 2.24 ) );
     Radii.insert( radio_atomico("La" , 1.88 ) );
     Radii.insert( radio_atomico("Ce" , 1.82 ) );
     Radii.insert( radio_atomico("Pr" , 1.82 ) );
     Radii.insert( radio_atomico("Nd" , 1.82 ) );
     Radii.insert( radio_atomico("Pm" , 1.81 ) );
     Radii.insert( radio_atomico("Sm" , 1.81 ) );
     Radii.insert( radio_atomico("Eu" , 2.06 ) );
     Radii.insert( radio_atomico("Gd" , 1.79 ) );
     Radii.insert( radio_atomico("Tb" , 1.77 ) );
     Radii.insert( radio_atomico("Dy" , 1.77 ) );
     Radii.insert( radio_atomico("Ho" , 1.76 ) );
     Radii.insert( radio_atomico("Er" , 1.75 ) );
     Radii.insert( radio_atomico("Tm" , 1.00 ) );
     Radii.insert( radio_atomico("Yb" , 1.94 ) );
     Radii.insert( radio_atomico("Lu" , 1.72 ) );
     Radii.insert( radio_atomico("Hf" , 1.59 ) );
     Radii.insert( radio_atomico("Ta" , 1.47 ) );
     Radii.insert( radio_atomico("W" ,  1.41 ) );
     Radii.insert( radio_atomico("Re" , 1.37 ) );
     Radii.insert( radio_atomico("Os" , 1.35 ) );
     Radii.insert( radio_atomico("Ir" , 1.35 ) );
     Radii.insert( radio_atomico("Pt" , 1.39 ) );
     Radii.insert( radio_atomico("Au" , 1.44 ) );
     Radii.insert( radio_atomico("Hg" , 1.55 ) );
     Radii.insert( radio_atomico("Tl" , 1.71 ) );
     Radii.insert( radio_atomico("Pb" , 1.75 ) );
     Radii.insert( radio_atomico("Bi" , 1.82 ) );
     Radii.insert( radio_atomico("Po" , 1.77 ) );
     Radii.insert( radio_atomico("At" , 1.62 ) );
     Radii.insert( radio_atomico("Rn" , 1.8  ) );
     Radii.insert( radio_atomico("Fr" , 1.00 ) );
     Radii.insert( radio_atomico("Ra" , 2.35 ) );
     Radii.insert( radio_atomico("Ac" , 2.03 ) );
     Radii.insert( radio_atomico("To" , 1.8  ) );
     Radii.insert( radio_atomico("Pa" , 1.63 ) );
     Radii.insert( radio_atomico("U" ,  1.56 ) );
     Radii.insert( radio_atomico("Np" , 1.56 ) );
     Radii.insert( radio_atomico("Pl" , 1.5  ) );
     Radii.insert( radio_atomico("Am" , 1.73 ) );
     Radii.insert( radio_atomico("Cm" , 1.8  ) );
     Radii.insert( radio_atomico("Bk" , 1.8  ) );
     Radii.insert( radio_atomico("Cf" , 1.8  ) );
     Radii.insert( radio_atomico("Es" , 1.8  ) );

    return Radii;
}


/********************************************************************/
/************************* Radii Assignment *************************/
/********** atom[i].R=assign_radii(Radios,atom[i].Symbol); **********/
/********************************************************************/

double assign_radii(map<string, double> _Radios, string _Symbol)
{
     map<string, double>::iterator p = _Radios.find(_Symbol);
     double _radii;
     if(p != _Radios.end())
     _radii= p->second;
     return _radii;
}

/********************************************************************/
/************************** Masas Constructor ***********************/
/******************** map<string, double> Masas; ********************/
/********************* Masas=masses_dictionary(); ********************/
/********************************************************************/

map<string, double> masses_dictionary()
{
  typedef pair<string, double> masa_atomica;
     // Lower Standard Atomic Weight: NIST
     // https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
     map<string, double> Masas;
     Masas.insert( masa_atomica("H"  , 1.00784     ) );
     Masas.insert( masa_atomica("He" , 4.002602    ) );
     Masas.insert( masa_atomica("Li" , 6.938       ) );
     Masas.insert( masa_atomica("Be" , 9.0121831   ) );
     Masas.insert( masa_atomica("B"  , 10.806      ) );
     Masas.insert( masa_atomica("C"  , 12.0096     ) );
     Masas.insert( masa_atomica("N"  , 14.00643    ) );
     Masas.insert( masa_atomica("O"  , 15.99903    ) );
     Masas.insert( masa_atomica("F"  , 18.9984031  ) );
     Masas.insert( masa_atomica("Ne" , 20.1797     ) );
     Masas.insert( masa_atomica("Na" , 22.989769   ) );
     Masas.insert( masa_atomica("Mg" , 24.304      ) );
     Masas.insert( masa_atomica("Al" , 26.981538   ) );
     Masas.insert( masa_atomica("Si" , 28.084      ) );
     Masas.insert( masa_atomica("P"  , 30.973761998) );
     Masas.insert( masa_atomica("S"  , 32.059      ) );
     Masas.insert( masa_atomica("Cl" , 35.446      ) );
     Masas.insert( masa_atomica("Ar" , 39.0983     ) );
     Masas.insert( masa_atomica("K"  , 39.0983     ) );
     Masas.insert( masa_atomica("Ca" , 40.0788     ) );
     Masas.insert( masa_atomica("Sc" , 44.95590    ) );
     Masas.insert( masa_atomica("Ti" , 47.867      ) );
     Masas.insert( masa_atomica("Va" , 50.9415     ) );
     Masas.insert( masa_atomica("Cr" , 51.9961     ) );
     Masas.insert( masa_atomica("Mn" , 54.938044   ) );
     Masas.insert( masa_atomica("Fe" , 55.845      ) );
     Masas.insert( masa_atomica("Co" , 58.933194   ) );
     Masas.insert( masa_atomica("Ni" , 58.6934     ) );
     Masas.insert( masa_atomica("Cu" , 63.546      ) );
     Masas.insert( masa_atomica("Zn" , 65.38       ) );
     Masas.insert( masa_atomica("Ga" , 69.723      ) );
     Masas.insert( masa_atomica("Ge" , 72.630      ) );
     Masas.insert( masa_atomica("As" , 74.921595   ) );
     Masas.insert( masa_atomica("Se" , 78.971      ) );
     Masas.insert( masa_atomica("Br" , 79.901      ) );
     Masas.insert( masa_atomica("Kr" , 83.798      ) );
     Masas.insert( masa_atomica("Rb" , 85.467      ) );
     Masas.insert( masa_atomica("Sr" , 87.62       ) );
     Masas.insert( masa_atomica("Y"  , 88.90584    ) );
     Masas.insert( masa_atomica("Zr" , 88.90584    ) );
     Masas.insert( masa_atomica("Nb" , 92.90637    ) );
     Masas.insert( masa_atomica("Mo" , 95.95       ) );
     Masas.insert( masa_atomica("Tc" , 98          ) );
     Masas.insert( masa_atomica("Ru" , 101.07      ) );
     Masas.insert( masa_atomica("Rh" , 102.90550   ) );
     Masas.insert( masa_atomica("Pd" , 106.42      ) );
     Masas.insert( masa_atomica("Ag" , 107.8682    ) );
     Masas.insert( masa_atomica("Cd" , 112.414     ) );
     Masas.insert( masa_atomica("In" , 114.818     ) );
     Masas.insert( masa_atomica("Sn" , 118.710     ) );
     Masas.insert( masa_atomica("Sb" , 121.760     ) );
     Masas.insert( masa_atomica("Te" , 127.60      ) );
     Masas.insert( masa_atomica("I"  , 126.90447   ) );
     Masas.insert( masa_atomica("Xe" , 131.293     ) );
     Masas.insert( masa_atomica("Cs" , 132.905451  ) );
     Masas.insert( masa_atomica("Ba" , 137.327     ) );
     Masas.insert( masa_atomica("La" , 138.90547   ) );
     Masas.insert( masa_atomica("Ce" , 140.116     ) );
     Masas.insert( masa_atomica("Pr" , 140.90766   ) );
     Masas.insert( masa_atomica("Nd" , 144.242     ) );
     Masas.insert( masa_atomica("Pm" , 145         ) );
     Masas.insert( masa_atomica("Sm" , 150.36      ) );
     Masas.insert( masa_atomica("Eu" , 151.964     ) );
     Masas.insert( masa_atomica("Gd" , 157.25      ) );
     Masas.insert( masa_atomica("Tb" , 158.92535   ) );
     Masas.insert( masa_atomica("Dy" , 162.500     ) );
     Masas.insert( masa_atomica("Ho" , 164.93033   ) );
     Masas.insert( masa_atomica("Er" , 167.259     ) );
     Masas.insert( masa_atomica("Tm" , 168.93422   ) );
     Masas.insert( masa_atomica("Yb" , 173.054     ) );
     Masas.insert( masa_atomica("Lu" , 174.9668    ) );
     Masas.insert( masa_atomica("Hf" , 178.49      ) );
     Masas.insert( masa_atomica("Ta" , 180.94788   ) );
     Masas.insert( masa_atomica("W"  , 183.84      ) );
     Masas.insert( masa_atomica("Re" , 183.84      ) );
     Masas.insert( masa_atomica("Os" , 190.23      ) );
     Masas.insert( masa_atomica("Ir" , 192.217     ) );
     Masas.insert( masa_atomica("Pt" , 195.084     ) );
     Masas.insert( masa_atomica("Au" , 196.966569  ) );
     Masas.insert( masa_atomica("Hg" , 200.592     ) );
     Masas.insert( masa_atomica("Tl" , 204.382     ) );
     Masas.insert( masa_atomica("Pb" , 207.2       ) );
     Masas.insert( masa_atomica("Bi" , 208.98040   ) );
     Masas.insert( masa_atomica("Po" , 209         ) );
     Masas.insert( masa_atomica("At" , 210         ) );
     Masas.insert( masa_atomica("Rn" , 222         ) );
     Masas.insert( masa_atomica("Fr" , 223         ) );
     Masas.insert( masa_atomica("Ra" , 226         ) );
     Masas.insert( masa_atomica("Ac" , 227         ) );
     Masas.insert( masa_atomica("To" , 232.0377    ) );
     Masas.insert( masa_atomica("Pa" , 231.03588   ) );
     Masas.insert( masa_atomica("U"  , 238.02891   ) );
     Masas.insert( masa_atomica("Np" , 237         ) );
     Masas.insert( masa_atomica("Pl" , 244         ) );
     Masas.insert( masa_atomica("Am" , 241.0568293 ) );
     Masas.insert( masa_atomica("Cm" , 243.0613893 ) );
     Masas.insert( masa_atomica("Bk" , 247.0703073 ) );
     Masas.insert( masa_atomica("Cf" , 249.0748539 ) );
     Masas.insert( masa_atomica("Es" , 252.082980  ) );

    return Masas;
}


/********************************************************************/
/************************ Masses Assignment *************************/
/*********** atom[i].M=assign_mass(Masas,atom[i].Symbol); **********/
/********************************************************************/

double assign_mass(map<string, double> _Masas, string _Symbol)
{
     map<string, double>::iterator p = _Masas.find(_Symbol);
     double _mass;
     if(p != _Masas.end())
     _mass= p->second;
     return _mass;
}

/********************************************************************/
/************************* Charge Constructor ***********************/
/******************** map<string, double> Cargas; *******************/
/******************* Radios=charges_dictionary(); *******************/
/********************************************************************/

map<string, int> charges_dictionary()
{
  typedef pair<string, int> carga_nucleo;

     map<string, int> Cargas;
     Cargas.insert( carga_nucleo("H"  , 1  ) );
     Cargas.insert( carga_nucleo("He" , 2  ) );
     Cargas.insert( carga_nucleo("Li" , 3  ) );
     Cargas.insert( carga_nucleo("Be" , 4  ) );
     Cargas.insert( carga_nucleo("B"  , 5  ) );
     Cargas.insert( carga_nucleo("C"  , 6  ) );
     Cargas.insert( carga_nucleo("N"  , 7  ) );
     Cargas.insert( carga_nucleo("O"  , 8  ) );
     Cargas.insert( carga_nucleo("F"  , 9  ) );
     Cargas.insert( carga_nucleo("Ne" , 10 ) );
     Cargas.insert( carga_nucleo("Na" , 11 ) );
     Cargas.insert( carga_nucleo("Mg" , 12 ) );
     Cargas.insert( carga_nucleo("Al" , 13 ) );
     Cargas.insert( carga_nucleo("Si" , 14 ) );
     Cargas.insert( carga_nucleo("P"  , 15 ) );
     Cargas.insert( carga_nucleo("S"  , 16 ) );
     Cargas.insert( carga_nucleo("Cl" , 17 ) );
     Cargas.insert( carga_nucleo("Ar" , 18 ) );
     Cargas.insert( carga_nucleo("K"  , 19 ) );
     Cargas.insert( carga_nucleo("Ca" , 20 ) );
     Cargas.insert( carga_nucleo("Sc" , 21 ) );
     Cargas.insert( carga_nucleo("Ti" , 22 ) );
     Cargas.insert( carga_nucleo("Va" , 23 ) );
     Cargas.insert( carga_nucleo("Cr" , 24 ) );
     Cargas.insert( carga_nucleo("Mn" , 25 ) );
     Cargas.insert( carga_nucleo("Fe" , 26 ) );
     Cargas.insert( carga_nucleo("Co" , 27 ) );
     Cargas.insert( carga_nucleo("Ni" , 28 ) );
     Cargas.insert( carga_nucleo("Cu" , 29 ) );
     Cargas.insert( carga_nucleo("Zn" , 30 ) );
     Cargas.insert( carga_nucleo("Ga" , 31 ) );
     Cargas.insert( carga_nucleo("Ge" , 32 ) );
     Cargas.insert( carga_nucleo("As" , 33 ) );
     Cargas.insert( carga_nucleo("Se" , 34 ) );
     Cargas.insert( carga_nucleo("Br" , 35 ) );
     Cargas.insert( carga_nucleo("Kr" , 36 ) );
     Cargas.insert( carga_nucleo("Rb" , 37 ) );
     Cargas.insert( carga_nucleo("Sr" , 38 ) );
     Cargas.insert( carga_nucleo("Y"  , 39 ) );
     Cargas.insert( carga_nucleo("Zr" , 40 ) );
     Cargas.insert( carga_nucleo("Nb" , 41 ) );
     Cargas.insert( carga_nucleo("Mo" , 42 ) );
     Cargas.insert( carga_nucleo("Tc" , 43 ) );
     Cargas.insert( carga_nucleo("Ru" , 44 ) );
     Cargas.insert( carga_nucleo("Rh" , 45 ) );
     Cargas.insert( carga_nucleo("Pd" , 46 ) );
     Cargas.insert( carga_nucleo("Ag" , 47 ) );
     Cargas.insert( carga_nucleo("Cd" , 48 ) );
     Cargas.insert( carga_nucleo("In" , 49 ) );
     Cargas.insert( carga_nucleo("Sn" , 50 ) );
     Cargas.insert( carga_nucleo("Sb" , 51 ) );
     Cargas.insert( carga_nucleo("Te" , 52 ) );
     Cargas.insert( carga_nucleo("I"  , 53 ) );
     Cargas.insert( carga_nucleo("Xe" , 54 ) );
     Cargas.insert( carga_nucleo("Cs" , 55 ) );
     Cargas.insert( carga_nucleo("Ba" , 56 ) );
     Cargas.insert( carga_nucleo("La" , 57 ) );
     Cargas.insert( carga_nucleo("Ce" , 58 ) );
     Cargas.insert( carga_nucleo("Pr" , 59 ) );
     Cargas.insert( carga_nucleo("Nd" , 60 ) );
     Cargas.insert( carga_nucleo("Pm" , 61 ) );
     Cargas.insert( carga_nucleo("Sm" , 62 ) );
     Cargas.insert( carga_nucleo("Eu" , 63 ) );
     Cargas.insert( carga_nucleo("Gd" , 64 ) );
     Cargas.insert( carga_nucleo("Tb" , 65 ) );
     Cargas.insert( carga_nucleo("Dy" , 66 ) );
     Cargas.insert( carga_nucleo("Ho" , 67 ) );
     Cargas.insert( carga_nucleo("Er" , 68 ) );
     Cargas.insert( carga_nucleo("Tm" , 69 ) );
     Cargas.insert( carga_nucleo("Yb" , 70 ) );
     Cargas.insert( carga_nucleo("Lu" , 71 ) );
     Cargas.insert( carga_nucleo("Hf" , 72 ) );
     Cargas.insert( carga_nucleo("Ta" , 73 ) );
     Cargas.insert( carga_nucleo("W" , 74 ) );
     Cargas.insert( carga_nucleo("Re" , 75 ) );
     Cargas.insert( carga_nucleo("Os" , 76 ) );
     Cargas.insert( carga_nucleo("Ir" , 77 ) );
     Cargas.insert( carga_nucleo("Pt" , 78 ) );
     Cargas.insert( carga_nucleo("Au" , 79 ) );
     Cargas.insert( carga_nucleo("Hg" , 80 ) );
     Cargas.insert( carga_nucleo("Tl" , 81 ) );
     Cargas.insert( carga_nucleo("Pb" , 82 ) );
     Cargas.insert( carga_nucleo("Bi" , 83 ) );
     Cargas.insert( carga_nucleo("Po" , 84 ) );
     Cargas.insert( carga_nucleo("At" , 85 ) );
     Cargas.insert( carga_nucleo("Rn" , 86 ) );
     Cargas.insert( carga_nucleo("Fr" , 87 ) );
     Cargas.insert( carga_nucleo("Ra" , 88 ) );
     Cargas.insert( carga_nucleo("Ac" , 89 ) );
     Cargas.insert( carga_nucleo("To" , 90 ) );
     Cargas.insert( carga_nucleo("Pa" , 91 ) );
     Cargas.insert( carga_nucleo("U"  , 92 ) );
     Cargas.insert( carga_nucleo("Np" , 93 ) );
     Cargas.insert( carga_nucleo("Pl" , 94 ) );
     Cargas.insert( carga_nucleo("Am" , 95 ) );
     Cargas.insert( carga_nucleo("Cm" , 96 ) );
     Cargas.insert( carga_nucleo("Bk" , 97 ) );
     Cargas.insert( carga_nucleo("Cf" , 98 ) );
     Cargas.insert( carga_nucleo("Es" , 99 ) );

    return Cargas;
}


/********************************************************************/
/*********************** Charges Assignment *************************/
/********** atom[i].Z=assign_charge(Cargas,atom[i].Symbol); **********/
/********************************************************************/

double assign_charge(map<string, int> _Cargas, string _Symbol)
{
     map<string, int>::iterator p = _Cargas.find(_Symbol);
     double _carga;
     if(p != _Cargas.end())
     _carga= p->second;
     return _carga;
}



/********************************************************************/
/****************** Constructor (from .xyz file) ********************/
/************** Molecule  molecule_name(string file); ***************/
/******************* Molecule Cysteine(argv[1]); ********************/
/********************************************************************/

Atomic_Structure::Atomic_Structure(string file)
{

   float x,y,z;
   string Symbol;
   string command="awk '{print $1\" \"$2\" \"$3\" \"$4}' ";
          command+=file;
          command+=" > coords.aux";
   system(command.c_str());
   system("cat coords.aux");
          command.clear();
          command="head -1 coords.aux >> Nat";
   system(command.c_str());

   ifstream Nat_file("Nat");
            Nat_file >> Nat;
            Nat_file.close();

   system("rm Nat");
          command.clear();
          command="head -";
          command+=to_string(Nat+2);
          command+=" ";
          command+=file;
          command+=" | tail -";
          command+=to_string(Nat);
          command+=" >> coordinatesAux ";
   system(command.c_str());

   ifstream coordinates_file("coordinatesAux");
   atom=new Atom[Nat+1];
   for(i=0;i<Nat;i++)
   {
      coordinates_file>>Symbol>>x>>y>>z;
      atom[i].read_Atom(Symbol,x,y,z);
   }
   coordinates_file.close();
   system("rm coordinatesAux");
   map<string, double> Radios;
   map<string, double> Masas;
   map<string, int > Cargas;
   Radios=radii_dictionary();
   Masas=masses_dictionary();
   Cargas=charges_dictionary();
   for(i=0;i<Nat;i++)
   {
      atom[i].R=assign_radii(Radios,atom[i].Symbol);
      atom[i].M=assign_mass(Masas,atom[i].Symbol);
      atom[i].Z=assign_charge(Cargas,atom[i].Symbol);
   }
}

/********************************************************************/
/******************************* z_max ******************************/
/************************** molecule.z_max() ************************/
/********************************************************************/

double Atomic_Structure::z_max()
{
   double current, maximum=atom[0].x[2]+atom[0].R;
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[2]+atom[i].R;
      if( current > maximum )
      {
         maximum=current;
      }
   }
   return maximum;
}


double Atomic_Structure::y_max()
{
   double current, maximum=atom[0].x[1]+atom[0].R;
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[1]+atom[i].R;
      if( current > maximum )
      {
         maximum=current;
      }
   }
   return maximum;
}


double Atomic_Structure::x_max()
{
   double current, maximum=atom[0].x[0]+atom[0].R;
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[0]+atom[i].R;
      if( current > maximum )
      {
         maximum=current;
      }
   }
   return maximum;
}
double Atomic_Structure::z_min()
{
   double current, minimum=atom[0].x[2]-atom[0].R;
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[2]-atom[i].R;
      if( current < minimum )
      {
         minimum=current;
      }
   }
   return minimum;
}

double Atomic_Structure::y_min()
{
   double current, minimum=atom[0].x[1]-atom[0].R;
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[1]-atom[i].R;
      if( current < minimum )
      {
         minimum=current;
      }
   }
   return minimum;
}

double Atomic_Structure::x_min()
{
   double current, minimum=atom[0].x[0]-atom[0].R;
   for (i=1;i<Nat;i++)
   {
      current=atom[i].x[0]-atom[i].R;
      if( current < minimum )
      {
         minimum=current;
      }
   }
   return minimum;
}


/******************************* fit_in *****************************/
/********** molecule.fit_in(xmin,xmax,ymin,ymax,zmin,zmax) **********/
/********************************************************************/

bool Atomic_Structure::fit_in(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax)
{
   float xrange=xmax-xmin;
   float yrange=ymax-ymin;
   float zrange=zmax-zmin;
   bool fit;
   float mol_xrange=x_max()-x_min();
   float mol_yrange=y_max()-y_min();
   float mol_zrange=z_max()-z_min();

   if(mol_xrange < xrange && mol_yrange < yrange && mol_zrange < zrange)
   {
      fit=true;
   }
   else
   {
      fit=false;
   }
   return fit;
}


/******************************* fit_in *****************************/
/********************** molecule.fit_in(super_cell) *****************/
/********************************************************************/

bool Atomic_Structure::fit_in(Simulation_Cell super_cell)
{
   float xrange=super_cell.M[0][0]+super_cell.M[1][0]+super_cell.M[2][0];
   float yrange=super_cell.M[0][1]+super_cell.M[1][1]+super_cell.M[2][1];
   float zrange=super_cell.M[0][2]+super_cell.M[1][2]+super_cell.M[2][2];
   bool fit;
   float mol_xrange=x_max()-x_min();
   float mol_yrange=y_max()-y_min();
   float mol_zrange=z_max()-z_min();

   if(mol_xrange < xrange && mol_yrange < yrange && mol_zrange < zrange)
   {
      fit=true;
   }
   else
   {
      fit=false;
   }
   return fit;
}

/***************************** read_xyz *****************************/
/********************* Molecule  molecule_name; *********************/
/***************** molecule_name.read_xyz(argv[1]); *****************/
/********************************************************************/

void Atomic_Structure::read_xyz(string file)
{

   float x,y,z;
   string Symbol;
   string command="awk '{print $1\" \"$2\" \"$3\" \"$4}' ";
          command+=file;
          command+=" > coords.aux";
   system(command.c_str());
          command.clear();
          command="head -1 coords.aux >> Nat";
   system(command.c_str());

   ifstream Nat_file("Nat");
            Nat_file >> Nat;
            Nat_file.close();

   system("rm Nat");
          command.clear();
          command="head -";
          command+=to_string(Nat+2);
          command+=" coords.aux | tail -";
          command+=to_string(Nat);
          command+=" > coordinatesAux ";
   system(command.c_str());
   system("rm coords.aux");

   ifstream coordinates_file("coordinatesAux");
   atom=new Atom[Nat+1];

   for(i=0;i<Nat;i++)
   {
      coordinates_file>>Symbol>>x>>y>>z;
      atom[i].read_Atom(Symbol,x,y,z);
   }

   coordinates_file.close();
   system("rm coordinatesAux");
   map<string, double> Radios;
   map<string, double> Masas;
   map<string, int> Cargas;
   Radios=radii_dictionary();
   Masas=masses_dictionary();
   Cargas=charges_dictionary();
   for(i=0;i<Nat;i++)
   {
      atom[i].R=assign_radii(Radios,atom[i].Symbol);
      atom[i].M=assign_mass(Masas,atom[i].Symbol);
      atom[i].Z=assign_charge(Cargas,atom[i].Symbol);
   }
}

/***************************** read_fhi *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.read_fhi("file.in"); ****************/
/********************************************************************/

void Atomic_Structure::read_fhi(string file)
{
   string com="cp ";
          com+=file;
          com+=" tmp.in";
   system(com.c_str());
   system("grep \"atom\" tmp.in | wc -l > tmp.xyz");
   system("echo \" \" >> tmp.xyz");
   system("grep \"atom\" tmp.in | awk '{print $5\" \"$2\" \"$3\" \"$4}' >> tmp.xyz");
   system("rm tmp.in");
   read_xyz("tmp.xyz");
   system("rm tmp.xyz");
}


/***************************** read_fhi *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.read_fhi("file.in"); ****************/
/********************************************************************/

void Crystal::read_fhi(string file)
{
   string com="cp ";
          com+=file;
          com+=" tmp.in";
   system(com.c_str());
   system("grep \"lattice_vector\" tmp.in | awk '{print $2\" \"$3\" \"$4}' >>vectors ");
   ifstream f("vectors");
   for(i=0;i<3;i++) //Lee los elementos de matriz
   {
      f>>lattice[i][0]>> lattice[i][1]>>lattice[i][2];
   }
   f.close();
   system("rm vectors");
   system("grep \"atom\" tmp.in | wc -l > tmp.xyz");
   system("echo \" \" >> tmp.xyz");
   system("grep \"atom\" tmp.in | awk '{print $5\" \"$2\" \"$3\" \"$4}' >> tmp.xyz");
   system("rm tmp.in");
   read_xyz("tmp.xyz");
   system("rm tmp.xyz");
}

/*************************** rand_generator *************************/
/********************** Cluster  Cluster_name; **********************/
/********** Cluster_name.rand_generator("Au",5,"Ir",3,range) ********/
/*************** Cluster_name.rand_generator("Ir",3) ****************/
/************ Cluster_name.rand_generator("Ir",3,range) *************/


void Cluster::rand_generator(string Symbol_1, int N_Symbol_1, string Symbol_2="AAA", int N_Symbol_2=0 )
{
   map<string, double> Radios;
   map<string, double> Masas;
   map<string, int> Cargas;
   Radios=radii_dictionary();
   Masas=masses_dictionary();
   Cargas=charges_dictionary();

   Nat=N_Symbol_1+N_Symbol_2;
   srand(time(NULL));
   atom=new Atom[Nat+1];
   int cont_S1=1;
   int cont_S2=0;
   int randomS;

   if(strcmp(Symbol_2.c_str(),"AAA") != 0 )
   {
      type="bimetallic";
   }
   else
   {
      type="monometallic";
   }
   double pseudo_rand[81][3];
   pseudo_rand[1][0]=-0.0195029 ; pseudo_rand[1][1]=0.0847575 ; pseudo_rand[1][2]=-0.0264053;
   pseudo_rand[2][0]=-0.029306 ; pseudo_rand[2][1]=-2.51719 ; pseudo_rand[2][2]=-0.0385384;
   pseudo_rand[3][0]=1.32679 ; pseudo_rand[3][1]=-1.18825 ; pseudo_rand[3][2]=1.79442;
   pseudo_rand[4][0]=1.24944 ; pseudo_rand[4][1]=-1.2225 ; pseudo_rand[4][2]=-1.89543;
   pseudo_rand[5][0]=-2.17871 ; pseudo_rand[5][1]=-1.22141 ; pseudo_rand[5][2]=0.640053;
   pseudo_rand[7][0]=-0.838767 ; pseudo_rand[7][1]=0.0835188 ; pseudo_rand[7][2]=2.45218;
   pseudo_rand[8][0]=2.59034 ; pseudo_rand[8][1]=0.0824258 ; pseudo_rand[8][2]=-0.0840119;
   pseudo_rand[6][0]=-1.35706 ; pseudo_rand[6][1]=-1.23125 ; pseudo_rand[6][2]=-1.83431;
   pseudo_rand[9][0]=1.32544 ; pseudo_rand[9][1]=1.38315 ; pseudo_rand[9][2]=1.79148;
   pseudo_rand[10][0]=1.24746 ; pseudo_rand[10][1]=1.38283 ; pseudo_rand[10][2]=-1.89791;
   pseudo_rand[11][0]=-2.17943 ; pseudo_rand[11][1]=1.38392 ; pseudo_rand[11][2]=0.636637;
   pseudo_rand[12][0]=-1.35706 ; pseudo_rand[12][1]=1.39101 ; pseudo_rand[12][2]=-1.83544;
   pseudo_rand[13][0]=-0.0259896 ; pseudo_rand[13][1]=2.68706 ; pseudo_rand[13][2]=-0.0362973;
   pseudo_rand[14][0]=-0.829979 ; pseudo_rand[14][1]=-2.48665 ; pseudo_rand[14][2]=2.42262;
   pseudo_rand[15][0]=2.55844 ; pseudo_rand[15][1]=-2.48773 ; pseudo_rand[15][2]=-0.0834814;
   pseudo_rand[16][0]=-0.0928186 ; pseudo_rand[16][1]=0.077041 ; pseudo_rand[16][2]=-3.65928;
   pseudo_rand[17][0]=-3.47169 ; pseudo_rand[17][1]=0.078118 ; pseudo_rand[17][2]=-1.16024;
   pseudo_rand[18][0]=-0.826892 ; pseudo_rand[18][1]=2.65276 ; pseudo_rand[18][2]=2.41527;
   pseudo_rand[19][0]=2.55262 ; pseudo_rand[19][1]=2.65168 ; pseudo_rand[19][2]=-0.0842503;
   pseudo_rand[20][0]=1.31501 ; pseudo_rand[20][1]=-3.74875 ; pseudo_rand[20][2]=1.7796;
   pseudo_rand[22][0]=-2.16839 ; pseudo_rand[22][1]=-3.82939 ; pseudo_rand[22][2]=0.623685;
   pseudo_rand[23][0]=-1.356 ; pseudo_rand[23][1]=-3.82677 ; pseudo_rand[23][2]=-1.83176;
   pseudo_rand[21][0]=1.22969 ; pseudo_rand[21][1]=-3.83047 ; pseudo_rand[21][2]=-1.88957;
   pseudo_rand[24][0]=0.524113 ; pseudo_rand[24][1]=-1.2125 ; pseudo_rand[24][2]=4.26029;
   pseudo_rand[25][0]=3.9196 ; pseudo_rand[25][1]=-1.21359 ; pseudo_rand[25][2]=1.74896;
   pseudo_rand[26][0]=-0.099016 ; pseudo_rand[26][1]=-2.54239 ; pseudo_rand[26][2]=-3.68792;
   pseudo_rand[27][0]=-3.49834 ; pseudo_rand[27][1]=-2.54131 ; pseudo_rand[27][2]=-1.17375;
   pseudo_rand[28][0]=2.66405 ; pseudo_rand[28][1]=0.0885545 ; pseudo_rand[28][2]=3.60194;
   pseudo_rand[29][0]=-2.99977 ; pseudo_rand[29][1]=-1.21263 ; pseudo_rand[29][2]=3.1135;
   pseudo_rand[30][0]=3.8549 ; pseudo_rand[30][1]=-1.21482 ; pseudo_rand[30][2]=-1.95627;
   pseudo_rand[31][0]=0.521314 ; pseudo_rand[31][1]=1.36949 ; pseudo_rand[31][2]=4.25938;
   pseudo_rand[32][0]=3.92062 ; pseudo_rand[32][1]=1.3684 ; pseudo_rand[32][2]=1.74523;
   pseudo_rand[33][0]=2.52371 ; pseudo_rand[33][1]=0.0752852 ; pseudo_rand[33][2]=-3.75257;
   pseudo_rand[34][0]=-4.32703 ; pseudo_rand[34][1]=0.0774689 ; pseudo_rand[34][2]=1.3143;
   pseudo_rand[36][0]=3.85082 ; pseudo_rand[36][1]=1.36547 ; pseudo_rand[36][2]=-1.95956;
   pseudo_rand[35][0]=-3.00066 ; pseudo_rand[35][1]=1.36766 ; pseudo_rand[35][2]=3.10785;
   pseudo_rand[37][0]=-2.69634 ; pseudo_rand[37][1]=0.0758621 ; pseudo_rand[37][2]=-3.64567;
   pseudo_rand[39][0]=-0.101921 ; pseudo_rand[39][1]=2.6949 ; pseudo_rand[39][2]=-3.69677;
   pseudo_rand[38][0]=1.31946 ; pseudo_rand[38][1]=3.98582 ; pseudo_rand[38][2]=1.78228;
   pseudo_rand[40][0]=-3.5038 ; pseudo_rand[40][1]=2.69599 ; pseudo_rand[40][2]=-1.18072;
   pseudo_rand[41][0]=-2.17698 ; pseudo_rand[41][1]=3.98769 ; pseudo_rand[41][2]=0.613141;
   pseudo_rand[42][0]=1.22535 ; pseudo_rand[42][1]=3.98661 ; pseudo_rand[42][2]=-1.90325;
   pseudo_rand[43][0]=-1.35967 ; pseudo_rand[43][1]=3.98151 ; pseudo_rand[43][2]=-1.84008;
   pseudo_rand[44][0]=-0.00349797 ; pseudo_rand[44][1]=-5.09938 ; pseudo_rand[44][2]=-0.00253142;
   pseudo_rand[45][0]=2.64629 ; pseudo_rand[45][1]=-2.50334 ; pseudo_rand[45][2]=3.57904;
   pseudo_rand[46][0]=2.48502 ; pseudo_rand[46][1]=-2.51022 ; pseudo_rand[46][2]=-3.72354;
   pseudo_rand[47][0]=-4.28901 ; pseudo_rand[47][1]=-2.50806 ; pseudo_rand[47][2]=1.28659;
   pseudo_rand[48][0]=-1.62841 ; pseudo_rand[48][1]=0.0750585 ; pseudo_rand[48][2]=4.88148;
   pseudo_rand[49][0]=5.14437 ; pseudo_rand[49][1]=0.0728996 ; pseudo_rand[49][2]=-0.127727;
   pseudo_rand[50][0]=-2.68054 ; pseudo_rand[50][1]=-2.49415 ; pseudo_rand[50][2]=-3.62319;
   pseudo_rand[51][0]=2.63614 ; pseudo_rand[51][1]=2.68213 ; pseudo_rand[51][2]=3.56308;
   pseudo_rand[52][0]=-4.29393 ; pseudo_rand[52][1]=2.66266 ; pseudo_rand[52][2]=1.28107;
   pseudo_rand[53][0]=2.48332 ; pseudo_rand[53][1]=2.6605 ; pseudo_rand[53][2]=-3.73144;
   pseudo_rand[54][0]=-2.68333 ; pseudo_rand[54][1]=2.6446 ; pseudo_rand[54][2]=-3.62919;
   pseudo_rand[55][0]=-0.0111239 ; pseudo_rand[55][1]=5.24589 ; pseudo_rand[55][2]=-0.0173006;
   pseudo_rand[57][0]=3.86127 ; pseudo_rand[57][1]=-3.78637 ; pseudo_rand[57][2]=1.71154;
   pseudo_rand[58][0]=-2.93427 ; pseudo_rand[58][1]=-3.76959 ; pseudo_rand[58][2]=3.06773;
   pseudo_rand[56][0]=0.504357 ; pseudo_rand[56][1]=-3.7853 ; pseudo_rand[56][2]=4.19434;
   pseudo_rand[59][0]=3.7909 ; pseudo_rand[59][1]=-3.77173 ; pseudo_rand[59][2]=-1.90626;
   pseudo_rand[60][0]=-0.82015 ; pseudo_rand[60][1]=-5.06436 ; pseudo_rand[60][2]=2.43582;
   pseudo_rand[62][0]=-1.61014 ; pseudo_rand[62][1]=-2.4979 ; pseudo_rand[62][2]=4.84251;
   pseudo_rand[61][0]=2.56712 ; pseudo_rand[61][1]=-5.06544 ; pseudo_rand[61][2]=-0.0694326;
   pseudo_rand[64][0]=1.14477 ; pseudo_rand[64][1]=-1.20367 ; pseudo_rand[64][2]=-5.4834;
   pseudo_rand[63][0]=5.1007 ; pseudo_rand[63][1]=-2.50004 ; pseudo_rand[63][2]=-0.120882;
   pseudo_rand[65][0]=-5.57876 ; pseudo_rand[65][1]=-1.20153 ; pseudo_rand[65][2]=-0.510616;
   pseudo_rand[66][0]=-1.40233 ; pseudo_rand[66][1]=-1.20338 ; pseudo_rand[66][2]=-5.44063;
   pseudo_rand[68][0]=-1.62222 ; pseudo_rand[68][1]=2.64533 ; pseudo_rand[68][2]=4.83796;
   pseudo_rand[69][0]=5.10201 ; pseudo_rand[69][1]=2.64319 ; pseudo_rand[69][2]=-0.135334;
   pseudo_rand[70][0]=1.14364 ; pseudo_rand[70][1]=1.35008 ; pseudo_rand[70][2]=-5.48765;
   pseudo_rand[71][0]=-5.58144 ; pseudo_rand[71][1]=1.35222 ; pseudo_rand[71][2]=-0.513724;
   pseudo_rand[67][0]=-4.79206 ; pseudo_rand[67][1]=-1.2023 ; pseudo_rand[67][2]=-2.93356;
   pseudo_rand[72][0]=3.87373 ; pseudo_rand[72][1]=3.92977 ; pseudo_rand[72][2]=1.68879;
   pseudo_rand[73][0]=0.482131 ; pseudo_rand[73][1]=3.93085 ; pseudo_rand[73][2]=4.19725;
   pseudo_rand[74][0]=-2.93769 ; pseudo_rand[74][1]=3.92401 ; pseudo_rand[74][2]=3.06034;
   pseudo_rand[75][0]=3.78801 ; pseudo_rand[75][1]=3.92187 ; pseudo_rand[75][2]=-1.91404;
   pseudo_rand[76][0]=-1.40406 ; pseudo_rand[76][1]=1.35043 ; pseudo_rand[76][2]=-5.44481;
   pseudo_rand[77][0]=-4.79449 ; pseudo_rand[77][1]=1.35151 ; pseudo_rand[77][2]=-2.93722;
   pseudo_rand[78][0]=-0.832405 ; pseudo_rand[78][1]=5.20814 ; pseudo_rand[78][2]=2.41864;
   pseudo_rand[79][0]=2.55852 ; pseudo_rand[79][1]=5.20706 ; pseudo_rand[79][2]=-0.0893115;
   pseudo_rand[80][0]=1.36354 ; pseudo_rand[80][1]=-6.22286 ; pseudo_rand[80][2]=1.84628;

   atom[0].Symbol=Symbol_1;  //o podría ser que
   atom[0].x[0]=pseudo_rand[1][0];           //se genere en un
   atom[0].x[1]=pseudo_rand[1][1];           //punto aleatorio
   atom[0].x[2]=pseudo_rand[1][2];
   atom[0].R=assign_radii(Radios, atom[0].Symbol);

   for(i=1;i<Nat;i++)
   {
      if(N_Symbol_2!=0)
      {
         randomS=rand()%2;
         if(cont_S1<N_Symbol_1 && randomS==0)
         {
            atom[i].Symbol=Symbol_1;
            cont_S1++;
         }
         else
         {
            if(cont_S2<N_Symbol_2 && randomS==1)
            {
               atom[i].Symbol=Symbol_2;
               cont_S2++;
            }
            else
            {
               if(cont_S1>=N_Symbol_1)
               {
                  atom[i].Symbol=Symbol_2;
                  cont_S2++;
               }
               else
               {
                  if(cont_S2>=N_Symbol_2)
                  {
                     atom[i].Symbol=Symbol_1;
                     cont_S1++;
                  }
               }

            }
        }
      }
      else
      {
         atom[i].Symbol=Symbol_1;
      }
      atom[i].R=assign_radii(Radios, atom[i].Symbol);

      for(j=0;j<3;j++)
      {
         atom[i].x[j]=pseudo_rand[i+1][j];
      }
   }
   for(i=0;i<Nat;i++)
   {
      atom[i].M=assign_mass(Masas,atom[i].Symbol);
      atom[i].Z=assign_charge(Cargas,atom[i].Symbol);
   }
}



/************************* srand_generator **************************/
/********************** Cluster  Cluster_name; **********************/
/********* Cluster_name.srand_generator("Au",5,"Ir",3,range) ********/
/************** Cluster_name.srand_generator("Ir",3) ****************/
/*********** Cluster_name.srand_generator("Ir",3,range) *************/


void Cluster::srand_generator(string Symbol_1, int N_Symbol_1, string Symbol_2="AAA", int N_Symbol_2=0,float compactness=1.5)
{

   map<string, double> Radios;
   Radios=radii_dictionary();
   map<string, double> Masas;
   Masas=masses_dictionary();
   map<string, int > Cargas;
   Cargas=charges_dictionary();
   Nat=N_Symbol_1+N_Symbol_2;
   int random;
   int randomS;
   double criterio;
   int accepted=0;
   int rejected=0;
   float Mx,Nx,My,Ny,Mz,Nz;
   float x,y,z;
   double Distance;
   srand(time(NULL));
   atom=new Atom[Nat+1];
   int cont_S1=1;
   int cont_S2=0;
   float _radio,phi,theta;
   int contador_radios;
   float volumen;
   float r_1,r_2, radio_max;
   int maximo_pasos_radio=Nat*10;

   if(strcmp(Symbol_2.c_str(),"AAA") != 0 )
   {
      type="bimetallic";
      r_1=assign_radii(Radios,Symbol_1);
      r_2=assign_radii(Radios,Symbol_2);
      volumen=8*(N_Symbol_1* pow(r_1,3)+N_Symbol_2*pow(r_2,3));
   }
   else
   {
      type="monometallic";
      r_1=assign_radii(Radios,Symbol_1);
      volumen=8*(N_Symbol_1* pow(r_1,3));
   }
   radio_max=compactness*pow((3*volumen/(4*3.1415926535)),0.333); //*0.75*pow(Nat,0.333);
   //////////// Con probabilidades ///////////////
   int num_enlaces;
   float probabilidad;
   ///////////////////////////////////////////////

///////////////////////////////Coloca el primer
///////////////////////////////átomo en el origen
   atom[0].Symbol=Symbol_1;  //o podría ser que
   atom[0].x[0]=0;           //se genere en un
   atom[0].x[1]=0;           //punto aleatorio
   atom[0].x[2]=0;           //o que admita como
///////////////////////////////argumento este punto
///////////////////////////////por default el origen


   atom[0].R=assign_radii(Radios, atom[0].Symbol);
   for(i=1;i<Nat;i++)
   {
      accepted=0;
      while(accepted==0)
      {
         random=(i-1)+(double)rand()/((double)RAND_MAX/(0-(i-1)+1)+1);

         if(N_Symbol_2!=0)
         {
            randomS=rand()%2;

            if(cont_S1<N_Symbol_1 && randomS==0)
            {
               atom[i].Symbol=Symbol_1;
            }
            else
            {
               if(cont_S2<N_Symbol_2 && randomS==1)
               {
                  atom[i].Symbol=Symbol_2;
               }
               else
               {
                  if(cont_S1>=N_Symbol_1)
                  {
                     atom[i].Symbol=Symbol_2;
                  }
                  else
                  {
                     if(cont_S2>=N_Symbol_2)
                     {
                        atom[i].Symbol=Symbol_1;
                     }
                  }

               }
            }
         }
         else
         {
            atom[i].Symbol=Symbol_1;
         }
         atom[i].R=assign_radii(Radios, atom[i].Symbol);

         phi=random_number(0,2*3.1415926535);
         theta=random_number(0,3.1415926535);
         //_radio=random_number(atom[random].R,epsilon*(atom[random].R+atom[i].R));
         _radio=atom[random].R+atom[i].R;

         Mx =_radio*sin(theta)*cos(phi);
         My =_radio*sin(theta)*sin(phi);
         Mz =_radio*cos(theta);

         atom[i].x[0]=atom[random].x[0]+Mx;
         atom[i].x[1]=atom[random].x[1]+My;
         atom[i].x[2]=atom[random].x[2]+Mz;

         rejected=0;

         for(j=0;j<i;j++)
         {
            if(contador_radios>maximo_pasos_radio)
            {
               radio_max=radio_max*1.01;
               contador_radios=0;
            }
            Distance=sqrt(pow(atom[i].x[0]-atom[j].x[0],2)+pow(atom[i].x[1]-atom[j].x[1],2)+pow(atom[i].x[2]-atom[j].x[2],2));
            criterio=atom[i].R+atom[j].R;
            if(Distance<criterio || Distance > radio_max)
            {
               rejected++;
               contador_radios++;
            }
         }

         if(rejected>0)
         {
            accepted=0;
         }
         else
         {
            accepted=1;
            contador_radios=0;
            if(strcmp(atom[i].Symbol.c_str(),Symbol_1.c_str()) == 0 )
            {
               cont_S1++;
            }
            else
            {
               cont_S2++;
            }

         }
         if(i<Nat)
         {
            continue;
         }
         else
         {
            break;
         }

      }

   }

   for(i=0;i<Nat;i++)
   {
      atom[i].M=assign_mass(Masas,atom[i].Symbol);
      atom[i].Z=assign_charge(Cargas,atom[i].Symbol);
   }
}


/************************* srand_generator **************************/
/********************** Cluster  Cluster_name; **********************/
/********* Cluster_name.srand_generator("Au",5,"Ir",3,range) ********/
/************** Cluster_name.srand_generator("Ir",3) ****************/
/*********** Cluster_name.srand_generator("Ir",3,range) *************/


void Cluster::srand_generator_box(string Symbol_1, int N_Symbol_1, string Symbol_2="AAA", int N_Symbol_2=0, float epsilon=2.5)
{

   map<string, double> Radios;
   Radios=radii_dictionary();
   map<string, double> Masas;
   Masas=masses_dictionary();
   map<string, int > Cargas;
   Cargas=charges_dictionary();
   Nat=N_Symbol_1+N_Symbol_2;
   int random;
   int randomS;
   double criterio;
   int accepted=0;
   int rejected=0;
   float Mx,Nx,My,Ny,Mz,Nz;
   float x,y,z;
   double Distance;
   srand(time(NULL));
   atom=new Atom[Nat+1];
   int cont_S1=1;
   int cont_S2=0;

   if(strcmp(Symbol_2.c_str(),"AAA") != 0 )
   {
      type="bimetallic";
   }
   else
   {
      type="monometallic";
   }
///////////////////////////////Coloca el primer
///////////////////////////////átomo en el origen
   atom[0].Symbol=Symbol_1;  //o podría ser que
   atom[0].x[0]=0;           //se genere en un
   atom[0].x[1]=0;           //punto aleatorio
   atom[0].x[2]=0;           //o que admita como
///////////////////////////////argumento este punto
///////////////////////////////por default el origen


   atom[0].R=assign_radii(Radios, atom[0].Symbol);
   for(i=1;i<Nat;i++)
   {
      accepted=0;
      while(accepted==0)
      {
         random=(i-1)+(double)rand()/((double)RAND_MAX/(0-(i-1)+1)+1);

         if(N_Symbol_2!=0)
         {
            randomS=rand()%2;

            if(cont_S1<N_Symbol_1 && randomS==0)
            {
               atom[i].Symbol=Symbol_1;
            }
            else
            {
               if(cont_S2<N_Symbol_2 && randomS==1)
               {
                  atom[i].Symbol=Symbol_2;
               }
               else
               {
                  if(cont_S1>=N_Symbol_1)
                  {
                     atom[i].Symbol=Symbol_2;
                  }
                  else
                  {
                     if(cont_S2>=N_Symbol_2)
                     {
                        atom[i].Symbol=Symbol_1;
                     }
                  }

               }
            }
         }
         else
         {
            atom[i].Symbol=Symbol_1;
         }
         atom[i].R=assign_radii(Radios, atom[i].Symbol);

         Mx=atom[random].x[0]+epsilon;  //Maximum range for random X
         Nx=atom[random].x[0]-epsilon;  //Minimun range for random x
         My=atom[random].x[1]+epsilon;  //Maximum range for random y
         Ny=atom[random].x[1]-epsilon;  //Minimum range for random y
         Mz=atom[random].x[2]+epsilon;  //Maximum range for random z
         Nz=atom[random].x[2]-epsilon;  //Minimum range for randon z

         atom[i].x[0]=Nx+ ((double)rand())/((double)RAND_MAX )* (Mx-Nx);
         atom[i].x[1]=Ny+ ((double)rand())/((double)RAND_MAX )* (My-Ny);
         atom[i].x[2]=Nz+ ((double)rand())/((double)RAND_MAX )* (Mz-Nz);

         rejected=0;

         for(j=0;j<i;j++)
         {
            Distance=sqrt(pow(atom[i].x[0]-atom[j].x[0],2)+pow(atom[i].x[1]-atom[j].x[1],2)+pow(atom[i].x[2]-atom[j].x[2],2));
            criterio=atom[i].R+atom[j].R;
            if(Distance<criterio)
            {
               rejected++;
            }
         }

         if(rejected>0)
         {
            accepted=0;
         }
         else
         {
            accepted=1;
            if(strcmp(atom[i].Symbol.c_str(),Symbol_1.c_str()) == 0 )
            {
               cont_S1++;
            }
            else
            {
               cont_S2++;
            }

         }
         if(i<Nat)
         {
            continue;
         }
         else
         {
            break;
         }

      }

   }

   for(i=0;i<Nat;i++)
   {
      atom[i].M=assign_mass(Masas,atom[i].Symbol);
      atom[i].Z=assign_charge(Cargas,atom[i].Symbol);
   }
}


/*

void Cluster::srand_generator_alternativo(string Symbol_1, int N_Symbol_1, string Symbol_2="AAA", int N_Symbol_2=0, float epsilon=1.0)
{

  map<string, double> Radios;
   Radios=radii_dictionary();
   Nat=N_Symbol_1+N_Symbol_2;
   int random;
   int randomS;
   double criterio;
   int accepted=0;
   int rejected=0;
   float Mx,Nx,My,Ny,Mz,Nz;
   float x,y,z;
   double Distance;
   srand(time(NULL));
   atom=new Atom[Nat+1];
   int cont_S1=1;
   int cont_S2=0;

   if(strcmp(Symbol_2.c_str(),"AAA") != 0 )
   {
      type="bimetallic";
   }
   else
   {
      type="monometallic";
   }
///////////////////////////////Coloca el primer
///////////////////////////////átomo en el origen
   atom[0].Symbol=Symbol_1;  //o podría ser que
   atom[0].x[0]=0;           //se genere en un
   atom[0].x[1]=0;           //punto aleatorio
   atom[0].x[2]=0;           //o que admita como
///////////////////////////////argumento este punto
///////////////////////////////por default el origen


   atom[0].R=assign_radii(Radios, atom[0].Symbol);
   for(i=1;i<Nat;i++)
   {
      accepted=0;
      while(accepted==0)
      {
         random=(i-1)+(double)rand()/((double)RAND_MAX/(0-(i-1)+1)+1);

         if(N_Symbol_2!=0)
         {
            randomS=rand()%2;

            if(cont_S1<N_Symbol_1 && randomS==0)
            {
               atom[i].Symbol=Symbol_1;
            }
            else
            {
               if(cont_S2<N_Symbol_2 && randomS==1)
               {
                  atom[i].Symbol=Symbol_2;
               }
               else
               {
                  if(cont_S1>=N_Symbol_1)
                  {
                     atom[i].Symbol=Symbol_2;
                  }
                  else
                  {
                     if(cont_S2>=N_Symbol_2)
                     {
                        atom[i].Symbol=Symbol_1;
                     }
                  }

               }
            }
         }
         else
         {
            atom[i].Symbol=Symbol_1;
         }

         atom[i].R=assign_radii(Radios, atom[i].Symbol);

         double Mrandom =  -1.0;
         double Nrandom =  1.0;

         atom[i].x[0]= atom[random].x[0] + ((atom[i].R) * (Mrandom + (double)rand()/((double)RAND_MAX/(Nrandom-Mrandom+1)+1) ));
         atom[i].x[1]= atom[random].x[1] + ((atom[i].R) * (Mrandom + (double)rand()/((double)RAND_MAX/(Nrandom-Mrandom+1)+1) ));
         atom[i].x[2]= atom[random].x[2] + ((atom[i].R) * (Mrandom + (double)rand()/((double)RAND_MAX/(Nrandom-Mrandom+1)+1) ));

         rejected=0;

         for(j=0;j<i;j++)
         {
            Distance=sqrt(pow(atom[i].x[0]-atom[j].x[0],2)+pow(atom[i].x[1]-atom[j].x[1],2)+pow(atom[i].x[2]-atom[j].x[2],2));
            criterio=atom[i].R+atom[j].R;

            double criterio2 = (2.2 * (1.5928235)*(pow(Nat,1.0/3.0))); // CHECK "as a function of size" e.g ---> (Nat / 2.0)
            if(Distance<criterio || Distance>criterio2)
            {
             cout<<"Criteria not reached"<<endl;
               rejected++;
            }
         }

         if(rejected>0)
         {
            accepted=0;
         }
         else
         {
            accepted=1;
            if(strcmp(atom[i].Symbol.c_str(),Symbol_1.c_str()) == 0 )
            {
               cont_S1++;
            }
            else
            {
               cont_S2++;
            }

         }
         if(i<Nat)
         {
            continue;
         }
         else
         {
            break;
         }

      }

   }

}


*/


/************************** print__xyz ******************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.print_xyz(argv[1]); *****************/
/********************************************************************/

void Atomic_Structure::print_xyz(string outputfile, string tag=" ")
{
   ofstream output_file(outputfile);
   output_file<<Nat<<endl;
   output_file<<tag<<endl;
   for(i=0;i<Nat;i++)
   {
      output_file<<atom[i].Symbol<<" "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<endl;
   }
   output_file.close();
}

/*****************************  move ********************************/
/************** molecule_name.move(DeltaX,DeltaY,DeltaZ) ************/
/********************** Cysteine.move(1,5,3) ************************/
/********************************************************************/

void Atomic_Structure::move(double Dx, double Dy, double Dz)
{
   for(i=0;i<Nat;i++)
   {
       atom[i].x[0]=atom[i].x[0]+Dx;
       atom[i].x[1]=atom[i].x[1]+Dy;
       atom[i].x[2]=atom[i].x[2]+Dz;
   }
}


/*************************** rotate_Rad *****************************/
/********** molecule_name.rotate_Rad(theta(rad), phi(rad)) **********/
/****************** Cysteine.rotate_Rad(pi,pi/2) ********************/
/********************************************************************/

void Molecule::rotate_Rad(float theta, float phi)
{
   float s;
   int l,m;
   float Rot[3][3];

   Rot[0][0]=cos(phi);                   //*************************//
   Rot[0][1]=(-cos(theta)*sin(phi));     //                         //
   Rot[0][2]=sin(phi)*sin(theta);        //                         //
   Rot[1][0]=sin(phi);                   //                         //
   Rot[1][1]=cos(theta)*cos(phi);        //     Rotation Matriz     //
   Rot[1][2]=(-sin(theta)*cos(phi));     //                         //
   Rot[2][0]=0;                          //                         //
   Rot[2][1]=sin(theta);                 //                         //
   Rot[2][2]=cos(theta);                 //*************************//

   for(i=0;i<Nat;i++)
   {
      for(l=0;l<3;l++)
      {
         s=0;
         for(m=0;m<3;m++)
         {
            s=s+(Rot[l][m]*(atom[i].x[m]));
         }

         atom[i].x[l]=s;
      }

   }

}

/*************************** rotate_Rad *****************************/
/********** molecule_name.rotate_Rad(theta(rad), phi(rad)) **********/
/****************** Cysteine.rotate_Rad(pi,pi/2) ********************/
/********************************************************************/

void Cluster::rotate_Rad(float theta, float phi)
{
   float s;
   int l,m;
   float Rot[3][3];

   Rot[0][0]=cos(phi);                   //*************************//
   Rot[0][1]=(-cos(theta)*sin(phi));     //                         //
   Rot[0][2]=sin(phi)*sin(theta);        //                         //
   Rot[1][0]=sin(phi);                   //                         //
   Rot[1][1]=cos(theta)*cos(phi);        //     Rotation Matriz     //
   Rot[1][2]=(-sin(theta)*cos(phi));     //                         //
   Rot[2][0]=0;                          //                         //
   Rot[2][1]=sin(theta);                 //                         //
   Rot[2][2]=cos(theta);                 //*************************//

   for(i=0;i<Nat;i++)
   {
      for(l=0;l<3;l++)
      {
         s=0;
         for(m=0;m<3;m++)
         {
            s=s+(Rot[l][m]*(atom[i].x[m]));
         }

         atom[i].x[l]=s;
      }

   }

}

/*************************** rotate_Deg *****************************/
/********** molecule_name.rotate_Deg(theta(deg), phi(deg)) **********/
/******************* Cysteine.rotate_Deg(45,90) *********************/
/********************************************************************/

void Molecule::rotate_Deg(float thetadeg, float phideg)
{
   float s;
   int l,m;
   float Rot[3][3];
   float theta,phi;

   theta=(thetadeg*(3.1415926535/180));
   phi=(phideg*(3.1415926535/180));

   Rot[0][0]=cos(phi);                   //*************************//
   Rot[0][1]=(-cos(theta)*sin(phi));     //                         //
   Rot[0][2]=sin(phi)*sin(theta);        //                         //
   Rot[1][0]=sin(phi);                   //                         //
   Rot[1][1]=cos(theta)*cos(phi);        //     Rotation Matriz     //
   Rot[1][2]=(-sin(theta)*cos(phi));     //                         //
   Rot[2][0]=0;                          //                         //
   Rot[2][1]=sin(theta);                 //                         //
   Rot[2][2]=cos(theta);                 //*************************//

   for(i=0;i<Nat;i++)
   {
      for(l=0;l<3;l++)
      {
         s=0;
         for(m=0;m<3;m++)
         {
            s=s+(Rot[l][m]*(atom[i].x[m]));
         }
         atom[i].x[l]=s;
      }
   }

}


/*************************** rotate_Deg *****************************/
/********** molecule_name.rotate_Deg(theta(deg), phi(deg)) **********/
/******************* Cysteine.rotate_Deg(45,90) *********************/
/********************************************************************/

void Cluster::rotate_Deg(float thetadeg, float phideg)
{
   float s;
   int l,m;
   float Rot[3][3];
   float theta,phi;

   theta=(thetadeg*(3.1415926535/180));
   phi=(phideg*(3.1415926535/180));

   Rot[0][0]=cos(phi);                   //*************************//
   Rot[0][1]=(-cos(theta)*sin(phi));     //                         //
   Rot[0][2]=sin(phi)*sin(theta);        //                         //
   Rot[1][0]=sin(phi);                   //                         //
   Rot[1][1]=cos(theta)*cos(phi);        //     Rotation Matriz     //
   Rot[1][2]=(-sin(theta)*cos(phi));     //                         //
   Rot[2][0]=0;                          //                         //
   Rot[2][1]=sin(theta);                 //                         //
   Rot[2][2]=cos(theta);                 //*************************//

   for(i=0;i<Nat;i++)
   {
      for(l=0;l<3;l++)
      {
         s=0;
         for(m=0;m<3;m++)
         {
            s=s+(Rot[l][m]*(atom[i].x[m]));
         }
         atom[i].x[l]=s;
      }
   }

}
/******************************  Kick *******************************/
/***************** molecule_name.kick(step_width) *******************/
/*********************** Cysteine.kick(0.8) *************************/
/********************************************************************/

void Cluster::kick(float step_width)
{
   srand(time(NULL));
   double Mrandom =  -1.0;
   double Nrandom =  1.0;

   for(i=0;i<Nat;i++)
   {
      for(j=0;j<3;j++)
      {
         atom[i].x[j]= atom[i].x[j] + ( (Mrandom + (double)rand()/((double)RAND_MAX/(Nrandom-Mrandom+1)+1)) * step_width/2.0 );
      }
   }
}

/************************  Kick_Lennard *****************************/
/************ cluster_name.kick_lennard(step_width) *****************/
/******************* Cysteine.kick_lennard(0.8) *********************/
/********************************************************************/

void Cluster::kick_lennard(float kick_size=1.0)
{

   int time_criterio=0;
   float time=0;
   float time_step=0.0001;
   float G_=-1.0;
   float kick_s=1.5;
   //float kick_size=1.50;  //1.5 funciona bien para 40
   float ax,ay,az;
   for(i=0;i<Nat;i++)
   {
      atom[i].x[0]=atom[i].x[0]*kick_s+random_number(-1,1)*kick_s/2;
      atom[i].x[1]=atom[i].x[1]*kick_s+random_number(-1,1)*kick_s/2;
      atom[i].x[2]=atom[i].x[2]*kick_s+random_number(-1,1)*kick_s/2;
   }
   float Dist;
   float sep;
   double swap;
   while(time_criterio==0)
   {
      for(i=0;i<Nat;i++) //CALCULA LAS FUERZAS
      {
         ax=0;ay=0;az=0;
         for(j=0;j<Nat;j++)
         {
            if(i!=j)
            {
               //Calcula
               Dist=Atomic_Distance(atom[i],atom[j]);
               ax=ax+G_*(4*((12*pow((atom[i].R+atom[j].R)*0.95,12)/pow(Dist,13))-(6*pow((atom[i].R+atom[j].R)*0.95,6)/pow(Dist,7)))*(atom[i].x[0]-atom[j].x[0]));
               ay=ay+G_*(4*((12*pow((atom[i].R+atom[j].R)*0.95,12)/pow(Dist,13))-(6*pow((atom[i].R+atom[j].R)*0.95,6)/pow(Dist,7)))*(atom[i].x[1]-atom[j].x[1]));
               az=az+G_*(4*((12*pow((atom[i].R+atom[j].R)*0.95,12)/pow(Dist,13))-(6*pow((atom[i].R+atom[j].R)*0.95,6)/pow(Dist,7)))*(atom[i].x[2]-atom[j].x[2]));
            }
         }
         atom[i].a[0]=ax*G_;
         atom[i].a[1]=ay*G_;
         atom[i].a[2]=az*G_;
      }
      for(i=0;i<Nat;i++)
      {
         //// Velocidades
         atom[i].v[0]=atom[i].v[0]+(time_step*atom[i].a[0]);//time;
         atom[i].v[1]=atom[i].v[1]+(time_step*atom[i].a[1]);//time;
         atom[i].v[2]=atom[i].v[2]+(time_step*atom[i].a[2]);//time;
         //// Posiciones
         atom[i].x[0]=atom[i].x[0]+(time_step*atom[i].v[0]);
         atom[i].x[1]=atom[i].x[1]+(time_step*atom[i].v[1]);
         atom[i].x[2]=atom[i].x[2]+(time_step*atom[i].v[2]);
      }
      time++;
      if(time>12500*kick_size)
      {
         time_criterio=1;
      }
   }
}

/******************************* swap *******************************/
/************************ swap(TimesToSwap) *************************/
/************************** Cluster.swap(3) *************************/

void Cluster::swap(int Nswap)
{
   int a,b,i;

   string aux, aux2;
   int M=0;
   int N=Nat-1;
   srand(time(NULL));


   i=0;

   if(strcmp(type.c_str(),"bimetallic") == 0 )
   {

     for(i=0;i<Nswap;i++)
     {
        a=M+(double)rand() / ((double)RAND_MAX / (N-M+1)+1);
        b=M+(double)rand() / ((double)RAND_MAX / (N-M+1)+1);

        if(strcmp(atom[a].Symbol.c_str(),atom[b].Symbol.c_str()) != 0 )
        {
           aux=atom[a].Symbol;
           aux2=atom[b].Symbol;

           atom[b].Symbol=aux;
           atom[a].Symbol=aux2;

        }
        else
        {
        i--;
        }
     }



   }

}


/******************************   +   *******************************/
/**************** molecule3 = molecule1 + molecule2 *****************/
/********************************************************************/

Atomic_Structure operator +(Atomic_Structure Mol1, Atomic_Structure Mol2)
{
   int j;
   Atomic_Structure Joined;
   Joined.Nat=Mol1.Nat+Mol2.Nat;

   Joined.atom=new Atom[Joined.Nat+1];

   for(i=0;i<Mol1.Nat;i++)
   {
      Joined.atom[i].Symbol=Mol1.atom[i].Symbol;
      Joined.atom[i].R=Mol1.atom[i].R;
      Joined.atom[i].Z=Mol1.atom[i].Z;
      Joined.atom[i].M=Mol1.atom[i].M;

      for(j=0;j<3;j++)
      {
         Joined.atom[i].x[j]=Mol1.atom[i].x[j];
      }
   }

   for(i=0;i<Mol2.Nat;i++)
   {
   Joined.atom[i+Mol1.Nat].Symbol=Mol2.atom[i].Symbol;
   Joined.atom[i+Mol1.Nat].R=Mol1.atom[i].R;
   Joined.atom[i+Mol1.Nat].Z=Mol1.atom[i].Z;
   Joined.atom[i+Mol1.Nat].M=Mol1.atom[i].M;
      for(j=0;j<3;j++)
      {
         Joined.atom[i+Mol1.Nat].x[j]=Mol2.atom[i].x[j];
      }
   }
   return Joined;
}


/******************************   +   *******************************/
/**************** molecule3 = molecule1 + molecule2 *****************/
/********************************************************************/

Cluster operator +(Cluster Mol1, Cluster Mol2)
{
   int j;
   Cluster Joined;
   Joined.Nat=Mol1.Nat+Mol2.Nat;

   Joined.atom=new Atom[Joined.Nat+1];

   for(i=0;i<Mol1.Nat;i++)
   {
      Joined.atom[i].Symbol=Mol1.atom[i].Symbol;
      Joined.atom[i].R=Mol1.atom[i].R;
      Joined.atom[i].Z=Mol1.atom[i].Z;
      Joined.atom[i].M=Mol1.atom[i].M;
      for(j=0;j<3;j++)
      {
         Joined.atom[i].x[j]=Mol1.atom[i].x[j];
      }
   }

   for(i=0;i<Mol2.Nat;i++)
   {
   Joined.atom[i+Mol1.Nat].Symbol=Mol2.atom[i].Symbol;
   Joined.atom[i+Mol1.Nat].R=Mol2.atom[i].R;
   Joined.atom[i+Mol1.Nat].Z=Mol2.atom[i].Z;
   Joined.atom[i+Mol1.Nat].M=Mol2.atom[i].M;
      for(j=0;j<3;j++)
      {
         Joined.atom[i+Mol1.Nat].x[j]=Mol2.atom[i].x[j];
      }
   }
   return Joined;
}




/********************* Function VASP_to_xyz  ************************/
/************** vasp_to_xyz(inputfile,outputfilet) ******************/
/***************** vasp_to_xyz(argv[1],argv[2]) *********************/

void VASP_to_xyz(string inputfile, string outputfile)
{

  struct Atom
  {
     string Symbol;
     float x[3];
  };

   struct Atom xyz[500];
   struct Atom poscar[500];

   int i, j, l, m;
   string Symbol[500];
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
   string command;
          command="cat "+inputfile+"  >> aux";
   system(command.c_str());
   system("grep \"elective\" aux | wc -l >>sel ; rm aux ");

   ifstream po("sel");
            po>>Sel;
            po.close();
   system("rm sel");
   command.clear();

   if(Sel==1)//Si hay selective dynamics)
   {
      command="tail -$(($(grep -A 500 \"elective\" "+inputfile+" | wc -l )-2)) ";
      command+= inputfile+" | grep .  >> selective  ";

      ofstream tg("commandi");
               tg<<command<<endl;
               tg.close();
      system("chmod +x commandi ");
      command.clear();

      system(" ./commandi ");
      system(" rm commandi");

      system("Nat=$(cat selective | wc -l  ) " );

      system(" awk '{print $4 \" \" $5 \" \" $6}' selective | tr 'T' '1' |tr   'F' '0' >>selectivedynamicsaux ");
      system( " echo \" \" >> selectivedynamics ; echo \" \" >>selectivedynamics ; cat selectivedynamicsaux >> selectivedynamics");
      system("rm selective selectivedynamicsaux");
   }
   else{}  //Creo que esto no importa
/************************************************/
   command = "grep \"irect\" "+inputfile+" | wc -l >> Direct";
   system(command.c_str());

   ifstream x("Direct");
            x>>Direct;
            x.close();
   system("rm Direct");
   command.clear();

   command= "grep \"artesian\" "+inputfile+" | wc -l >> Cartesian";
   system(command.c_str());
   command.clear();
   ifstream xi("Cartesian");
            xi>>Cartesian;
            xi.close();
   system("rm Cartesian");

   if(Direct==1)
   {
      command="head -5 "+ inputfile+ " | tail -3 >> Matriz";
      system(command.c_str());
      command.clear();
      command="head -2 "+inputfile+ " | tail -1 >> Factor";
      system(command.c_str());
      command.clear();
      command= "tail -$(($(grep -A 500 \"irect\" "+inputfile+" | wc -l )-1)) ";
      command+= inputfile+ " | awk '{ print $1 \"  \" $2 \"  \" $3  }' | grep . >> Posiciones ";
      system(command.c_str());
      command.clear();
      command= "head -6 "+inputfile+"  | tail -1 >> aux";
      system(command.c_str());
      command.clear();
      command= "head -7 "+inputfile+"  | tail -1 >> aux2";
      system(command.c_str());
      command.clear();
      system("tr -s '[:blank:]' '\n' < aux >> Simbolos1");
      system("tr -s '[:blank:]' '\n' < aux2 >> Numeros");
      system("grep  . Simbolos1 >> Simbolos");
      system("rm Simbolos1");
      system("cat Simbolos | wc -l >>  Ntyp");

      ifstream f("Matriz");
      ifstream g("Factor");
      ifstream h("Posiciones");
      ifstream q("Simbolos");
      ifstream r("Numeros");
      ifstream p("Ntyp");
      ofstream o(outputfile);

      p>>Ntyp;
      suma=0;

      for(i=0;i<Ntyp;i++)  //Aca imagino podriamos hacer el feof c++
      {
         q>>Symbol[i];
         r>>N_Symbol[i];
         suma=suma+N_Symbol[i];
      }
      Nat=suma;
      l=0;
      for(i=0;i<Ntyp;i++) //Reconoce que tipo de atomo es cada vector posicion
      {
         for(j=0;j<N_Symbol[i];j++)
         {
            xyz[l].Symbol= Symbol[i];
            l=l+1;
         }
      }

      for(i=0;i<3;i++) //Lee los elementos de matriz
      {
         f>>Mi[i][0]>> Mi[i][1]>>Mi[i][2];
      }

      g>>Factor;  //Lee el factor de escala

      for(i=0;i<3;i++)//Multiplica el factor de escala a la matriz
      {
         for(j=0;j<3;j++)
         {
            M[i][j]=Factor*M[j][i];
         }
      }
      for(i=0;i<Nat;i++) //Extrae las coordenadas de Positions
      {
         h>>poscar[i].x[0]>>poscar[i].x[1]>>poscar[i].x[2];
      }

      for(i=0;i<Nat;i++)//Aplica la matriz  de cambio a cada atomo
      {
         for(l=0;l<3;l++)
         {
            for(m=0;m<3;m++)
               {
                  M[l][m]=poscar[i].x[l]*Mi[l][m];
               }
         }
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
   o<<Nat<<endl<<endl;

   for(i=0;i<Nat;i++) //Imprime las lineas de atomos
   {
      o<<xyz[i].Symbol<<" "<< xyz[i].x[0]<<" "<< xyz[i].x[1]<<" "<< xyz[i].x[2]<<endl;
   }
   o.close();

   }//Cierra if/********************************************************/
   else
   {
      if(Cartesian==1)
      {
         command.clear();
         command="tail -$(($(grep -A 500 \"artesian\" "+inputfile+" | wc -l )-1)) "+inputfile+" | awk '{ print $1 \"  \" $2 \"  \" $3  }' | grep . >> Posiciones ";
         system(command.c_str());
         command.clear();
         command= "head -6 "+ inputfile+ "  | tail -1 >> aux";
         system(command.c_str());
         command.clear();
         command= "head -7 "+ inputfile+"  | tail -1 >> aux2";
         system(command.c_str());
         command.clear();
         command="tr -s '[:blank:]' '\n' < aux >> Simbolos1";
         system(command.c_str());
         command.clear();
         command="tr -s '[:blank:]' '\n' < aux2 >> Numeros";
         system(command.c_str());
         command.clear();
         system("grep  . Simbolos1 >> Simbolos");
         system("rm Simbolos1");
         system("cat Simbolos | wc -l >>  Ntyp");//Guarda el numero de especies como archivo


         ifstream hi("Posiciones");
         ifstream qi("Simbolos");
         ifstream ri("Numeros");
         ifstream pi("Ntyp");
         ofstream oi(outputfile);

         pi>>Ntyp;

         for(i=0;i<Ntyp;i++)
         {
            qi>>Symbol[i];
            ri>>N_Symbol[i];
            suma=suma+N_Symbol[i];
         }
         Nat=suma;
         l=0;
         for(i=0;i<Ntyp;i++) //Reconoce que tipo de atomo es cada vector posicion
         {
            for(j=0;j<N_Symbol[i];j++)
            {
               xyz[l].Symbol=Symbol[i];
               l=l+1;
            }
         }
         for(i=0;i<Nat;i++) //Extrae las coordenadas de Positions
         {
            hi>>poscar[i].x[0]>>poscar[i].x[1]>>poscar[i].x[2];
         }
         system("rm  Posiciones Ntyp Simbolos Numeros aux aux2"); //Borra archivos residuales
         oi<<Nat<<endl<<endl;
         for(i=0;i<Nat;i++) //Imprime las lineas de atomos
         {
            oi<<xyz[i].Symbol<<" "<< poscar[i].x[0]<<" "<< poscar[i].x[1]<<" "<< poscar[i].x[2]<<endl;
         }
         oi.clear();
      }
   }

   if (Sel==1)
   {
      command="cat "+outputfile+" >> aux ;  rm "+outputfile+" ; paste aux selectivedynamics  >> "+outputfile;
      system(command.c_str());
      system("rm  aux selectivedynamics");
   }
}//Cierra funcion


/********************* Function xyz_to_VASP  ************************/
/** xyz_to_vasp(InFile,OutFile,Title,ScaleFactor,Matrix(3x3Array) ***/
/******float M[3][3]; xyz_to_vasp(mol.xyz,POSCAR,"Title",1,M);*******/

void xyz_to_VASP(string inputfile,string outputfile,string Titulo, float Factor, double M[3][3])  //Checa la matriz
{
   int i,j;
   int Selective;
   string command;
   ofstream hj(outputfile);
   hj<<Titulo<<endl;
   hj<<Factor<<endl;

   for(i=0;i<3;i++)           //Imprime Matriz de vectores
   {
      for(j=0;j<3;j++)
      {
         hj<<M[i][j]<<"  ";
      }
      hj<<endl;
   }
   hj.close();

   command="cat "+inputfile+" >> selective.txt";
   system(command.c_str());
//   system("head -$(($(head -1 selective.txt )+2)) selective.txt | tail -$(head -1 selective.txt) | awk '{print $5}'");
   command.clear();
   command="awk '{print $5 }'  selective.txt  | grep . | wc -l >> sel";
   system(command.c_str());

   ifstream iy("sel");
   iy>>Selective;
   iy.close();
   system("rm sel");
   if(Selective>0)
   {
      system("head -$(($(head -1 selective.txt )+2)) selective.txt | tail -$(head -1 selective.txt) | awk '{print $5  \" \" $6 \" \" $7 }'  | tr '0' 'T'  >> aux");
      system("rm selective.txt ; cat aux >> selective.txt ; rm aux ");
   }
/////////////////////////////////||||||||||||

   ofstream pl("command");
   pl<< "#!/bin/bash"<<endl;
   pl<<"N=$(head -1 $1)"<<endl;
   pl<<"head -$(($N+2)) $1 | tail -$N |  awk '{print $1 }' | sort -u >> aux"<<endl;
   pl<<"M=$(wc -l aux | awk '{ print $1 }')"<<endl;
   pl<<"for ((i=1;i<$(($M+1));i++))"<<endl;
   pl<<"do"<<endl;
   pl<<"echo -n \"$(head -$i aux | tail -1  )  \" >>$2"<<endl;
   pl<<"done"<<endl;
   pl<<"echo \" \" >>$2"<<endl;

   pl<<"for ((i=1;i<$(($M+1));i++))"<<endl;
   pl<<"do"<<endl;
   pl<<"echo -n \"$( grep -w \"$(head -$i aux | tail -1 )\" $1 | wc -l )  \">> $2"<<endl;
   pl<<"done"<<endl;
   pl<<"echo \" \" >>$2"<<endl;

   if(Selective>0) //Lo que significa que hay al menos una linea seleccionada
   {
      pl<<"echo \"Selective dynamics\" >>$2 "<<endl;
   }
   pl<<"echo \"Cartesian\">>$2 "<<endl;


   pl<<"rm aux"<<endl;
   pl<<"head -$(($N+2)) $1 | tail -$N |  sort -u | awk '{print $2 \"  \"$3 \"  \" $4  }'  >> aux1"<<endl;
   pl<<"head -$(($N+2)) $1 | tail -$N |  sort -u | awk '{print $5 \"  \"$6 \"  \" $7  }'  >> aux2"<<endl;

   if(Selective>0) //Lo que significa que hay al menos una linea seleccionada
   {
      pl<<"Nl=$(cat aux2 | wc -l ) \n for((i=1;i<$(($N+1)); i++)) \n do \n cont=$(head -$i aux2 | tail -1 | tr '0' 'F' | tr '1' 'T' ) "<<endl;
      pl<<" ki=$(echo $cont | wc -c) \n if [ $ki -gt 1 ] \n then echo \"$cont\" >>aux3 \n else \n echo \" T   T   T \" >> aux3  "<<endl;
      pl<<"fi \n done \n paste aux1 aux3 >> $2 \n rm aux1 aux2   aux3 "<<endl;
   }
   else
   {
      pl<<"paste aux1 aux2 >>$2 ; rm aux2 aux1  "<<endl;
   }
   pl.close();
   system("chmod +x command");

   command.clear();
   command="./command  "+inputfile+ " "+outputfile;
   system(command.c_str());
   system("rm selective.txt command");
}


/**************************** read_VASP *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.read_VASP(argv[1]); *****************/
/********************************************************************/

void Atomic_Structure::read_VASP(string file)
{
   VASP_to_xyz(file, "salida");
   read_xyz("salida");
   system("rm salida");
}


/***************************** centroid *****************************/
/********************* molecule_name.centroid() *********************/
/************************ Cysteine.centroid *************************/
/********************************************************************/

void Molecule::centroid()
{
   float Dx, Dy, Dz, sumax=0,sumay=0,sumaz=0;
   for(i=0;i<Nat;i++)
   {
      sumax=atom[i].x[0]+sumax;
      sumay=atom[i].x[1]+sumay;
      sumaz=atom[i].x[2]+sumaz;
   }
   Dx=sumax/Nat;
   Dy=sumay/Nat;
   Dz=sumaz/Nat;

   for(i=0;i<Nat;i++)
   {
       atom[i].x[0]=atom[i].x[0]-Dx;
       atom[i].x[1]=atom[i].x[1]-Dy;
       atom[i].x[2]=atom[i].x[2]-Dz;
   }
}
/***************************** centroid *****************************/
/********************* molecule_name.centroid() *********************/
/************************ Cysteine.centroid *************************/
/********************************************************************/

void Cluster::centroid()
{
   float Dx, Dy, Dz, sumax=0,sumay=0,sumaz=0;
   for(i=0;i<Nat;i++)
   {
      sumax=atom[i].x[0]+sumax;
      sumay=atom[i].x[1]+sumay;
      sumaz=atom[i].x[2]+sumaz;
   }
   Dx=sumax/Nat;
   Dy=sumay/Nat;
   Dz=sumaz/Nat;

   for(i=0;i<Nat;i++)
   {
       atom[i].x[0]=atom[i].x[0]-Dx;
       atom[i].x[1]=atom[i].x[1]-Dy;
       atom[i].x[2]=atom[i].x[2]-Dz;
   }
}

/**************************** print_VASP ****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.print_VASP(argv[1]); ****************/
/********************************************************************/

void Atomic_Structure::print_VASP(string outputfile,string Titulo, float Factor, double M[3][3])
{
   ofstream output_file("output");
   output_file<<Nat<<endl;
   output_file<<" "<<endl;
   for(i=0;i<Nat;i++)
   {
      output_file<<atom[i].Symbol<<" "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<endl;
   }
   output_file.close();
   xyz_to_VASP("output", outputfile, Titulo, Factor, M);
   system("rm output");
}


/**************************** print_VASP ****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.print_VASP(argv[1]); ****************/
/********************************************************************/

void Crystal::print_VASP(string outputfile,string Titulo, float Factor)
{
   ofstream output_file("output");
   output_file<<Nat<<endl;
   output_file<<" "<<endl;
   for(i=0;i<Nat;i++)
   {
      output_file<<atom[i].Symbol<<" "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<endl;
   }
   output_file.close();
   xyz_to_VASP("output", outputfile, Titulo, Factor, lattice);
   system("rm output");
}

/***************************** print_fhi ****************************/
/********************* Molecule  molecule_name; *********************/
/***************** molecule_name.print_fhi(argv[1]); ****************/
/********************************************************************/

void Atomic_Structure::print_fhi(string outputfile)
{
   ofstream output_file(outputfile);
   for(i=0;i<Nat;i++)
   {
      output_file<<"atom "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<" "<<atom[i].Symbol<<endl;
   }
   output_file.close();
}


/***************************** print_fhi ****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.print_fhi(argv[1]); ****************/
/********************************************************************/

void Crystal::print_fhi(string outputfile)
{
   ofstream output_file(outputfile);

   for(i=0;i<3;i++)
   {
      output_file<<"lattice_vector "<<lattice[i][0]<<" "<<lattice[i][1]<<" "<<lattice[i][2]<<endl;
   }
   for(i=0;i<Nat;i++)
   {
      output_file<<"atom "<<atom[i].x[0]<<" "<<atom[i].x[1]<<" "<<atom[i].x[2]<<" "<<atom[i].Symbol<<endl;
   }
   output_file.close();
}


/**************************** read_VASP *****************************/
/********************* Molecule  molecule_name; *********************/
/**************** molecule_name.read_VASP(argv[1]); *****************/
/********************************************************************/

void Crystal::read_VASP(string file)
{
   VASP_to_xyz(file, "salida");
   read_xyz("salida");
   system("rm salida");
   string command="head -5 ";
          command+=file;
          command+=" | tail -3 > matriz ";
   system(command.c_str());
   ifstream fm("matriz");
   for(i=0;i<3;i++) //Lee los elementos de matriz
   {
      fm>>lattice[i][0]>> lattice[i][1]>>lattice[i][2];
   }
   fm.close();
   system("rm matriz");

}

Cluster extract(string estruct, string symbol)
{
   Cluster aux;
   string command="grep \"";
          command+=symbol;
          command+="\" ";
          command+=estruct;
          command+=" >> clus.fhi";
   system(command.c_str());

   aux.read_fhi("clus.fhi");
   aux.print_fhi("impr.xyz");
   system("rm impr.xyz");
   system("rm clus.fhi");
   return aux;
}


/******************************* show *******************************/
/**************** molecule_name.show(visualizer=iqmol); ****************/
/********************************************************************/

void Atomic_Structure::show(string visualizer="iqmol")
{
   print_xyz("visualizer_tmp.xyz");
   if(visualizer=="iqmol")
   {
      system("iqmol visualizer_tmp.xyz 2> /dev/null > /dev/null");
   }
   else
   {
      if(visualizer=="avogadro")
      {
         system("avogadro visualizer_tmp.xyz 2> /dev/null > /dev/null");
      }
      else
      {
         if(visualizer=="vesta" || visualizer=="VESTA")
         {
            system("VESTA visualizer_tmp.xyz 2> /dev/null > /dev/null");
         }
      }
   }
   system("rm visualizer_tmp.xyz");
}


/***************************** make_movie ***************************/
/****************** make_movie(estruct,archivo.xyz); ****************/
/********************************************************************/
/*
void make_movie(Atomic_Structure estructura, string file)
{
    estructura[0].print_xyz("file_movie");
    for(i=1;i<estructura.Nat;i++)
    {
       estructura[i].print_xyz("movie_tmp.xyz");
       system("cat movie_tmp.xyz >> file_movie")
    }
    string command="mv file_movie ";
           command+=file;
    system(command.c_str());
}
*/
/*************************** show_movie  ****************************/
/************* show_movie(visualizer=iqmol, movie.xyz); *************/
/********************************************************************/
/*
void show_movie(string visualizer="iqmol", string file)
{
   string command="cp ";
         command+=file;
         command+="  visualizer_tmp.xyz";
   system(command.c_str());

   if(visualizer=="iqmol")
   {
      system("iqmol visualizer_tmp.xyz 2> /dev/null > /dev/null");
   }
   else
   {
      if(visualizer=="avogadro")
      {
         system("avogadro visualizer_tmp.xyz 2> /dev/null > /dev/null");
      }
      else
      {
         if(visualizer=="vesta")
         {
            system("VESTA visualizer_tmp.xyz 2> /dev/null > /dev/null");
         }
      }
   }
   system("rm visualizer_tmp.xyz");
}
*/
string read_pipe(string command) {
   char buffer[128];
   string result = "";

   FILE* pipe = popen(command.c_str(), "r");
   while (!feof(pipe)) {
      if (fgets(buffer, 128, pipe) != NULL)
         result += buffer;
   }

   pclose(pipe);
   return result;
}

string read_bash(string command)
{
   string variable;
   command+=" > vArIaBlE_to_rEaD";
   system(command.c_str());
   ifstream bar("vArIaBlE_to_rEaD");
   bar>>variable;
   bar.close();
   return variable;
}

string string_pipe(string cmd,string defecto="AAA")
{
//https://www.jeremymorgan.com/tutorials/c-programming/how-to-capture-the-output-of-a-linux-command-in-c/
   string data;
   FILE * stream;
   const int max_buffer = 256;
   char buffer[max_buffer];
   cmd.append(" 2>&1");
   stream = popen(cmd.c_str(), "r");
   if (stream)
   {
      while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
      pclose(stream);
   }
   data.pop_back();
   if(data.length()>0)
   {
      return data;
   }
   else
   {
      return defecto;
   }

}
int int_pipe(string cmd,int defecto=0)
{
   string data;
   FILE * stream;
   const int max_buffer = 256;
   char buffer[max_buffer];
   cmd.append(" 2>&1");
   stream = popen(cmd.c_str(), "r");
   if (stream)
   {
      while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
      pclose(stream);
   }
   if(data.length()>1)
   {
      return stoi(data);
   }
   else
   {
      return defecto;
   }
}

float float_pipe(string cmd,float defecto=0.0)
{
   string data;
   FILE * stream;
   const int max_buffer = 256;
   char buffer[max_buffer];
   cmd.append(" 2>&1");
   stream = popen(cmd.c_str(), "r");
   if (stream)
   {
      while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
      pclose(stream);
   }
   if(data.length()>1)
   {
      return stof(data);
   }
   else
   {
      return defecto;
   }
}


double double_pipe(string cmd,float defecto=0.0)
{
   string data;
   FILE * stream;
   const int max_buffer = 256;
   char buffer[max_buffer];
   cmd.append(" 2>&1");
   stream = popen(cmd.c_str(), "r");
   if (stream)
   {
      while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
      pclose(stream);
   }
   if(data.length()>1)
   {
      return stod(data);
   }
   else
   {
      return defecto;
   }
}



/************************** coulomb_matrix **************************/
/********************************************************************/
/**************** Molecule  molecule_1, molecule_2; *****************/
/*************** c_matrix=coulomb_matrix(molecule_1); ***************/
/********************************************************************/

void coulomb_matrix(Atomic_Structure Molecule1, string energy, string file,bool diagonal)
{
   double  coulomb[Molecule1.Nat][Molecule1.Nat];
   // Calculates Coulomb Matrix
   for(i=0;i<Molecule1.Nat;i++)
   {
      for(j=0;j<Molecule1.Nat;j++)
      {
         if(i==j)
         {
            coulomb[i][j]=0.5*pow(Molecule1.atom[i].Z,2.4);
         }
         else
         {
           coulomb[i][j]=(Molecule1.atom[i].Z)*(Molecule1.atom[j].Z)/Atomic_Distance(Molecule1.atom[i],Molecule1.atom[j]);
         }
      }
   }
   // Flat Coulomb Matrix into a single vector

   ofstream dat(file,std::ios_base::app);
   for(i=0;i<Molecule1.Nat;i++)
   {
      for(j=i;j<Molecule1.Nat;j++)
      {
           if(diagonal && (i==j))
           {
              dat<<coulomb[i][j]<<",";
           }
           else if(j!=i)
           {
              dat<<coulomb[i][j]<<",";
           }
      }
   }
   dat<<energy<<endl;
   dat.close();

}

/********************************************************************/
/************************** coulomb_matrix **************************/
/**************** Molecule  molecule_1, molecule_2; *****************/
/*************** c_matrix=coulomb_matrix(molecule_1); ***************/
/********************************************************************/

void header_csv(Atomic_Structure Molecule1, string file,bool diagonal)
{
   int conter=0;
   ofstream dat(file,std::ios_base::app);
   for(i=0;i<Molecule1.Nat;i++)
   {
      for(j=i;j<Molecule1.Nat;j++)
      {
           if(diagonal && (i==j))
           {
              dat<<conter<<",";
           }
           else if(j!=i)
           {
              dat<<conter<<",";
           }
           conter++;
      }
   }
   dat<<"Eat"<<endl;
   dat.close();

}

void output_qe_to_movie(string file, string file_movie)
{
   string comm=" cp  ";
          comm+=file+" output.temp ; n_steps=$(grep \"ATOMIC_POSITIONS\" output.temp  | wc -l) ; ";
          comm+="nats=$(grep \"number of atoms\" output.temp | awk '{print $5}') ; rep=\"ATOMIC_POSITIONS (angstrom)\" ; ";
          comm+="rep2=\"$nats\\n\" ; ";
          comm+="grep -A $nats \"ATOMIC_POSITIONS\" output.temp   | sed \"s/$rep/$rep2/g\" | sed '/--/d' > "+file_movie;
   system(comm.c_str());
}

void output_qe_split(string file, string dir)
{
   string comm=" #!/bin/bash  \n   cp  ";
          comm+=file+" output.temp ; n_steps=$(grep \"ATOMIC_POSITIONS\" output.temp  | wc -l) ; ";
          comm+="nats=$(grep \"number of atoms\" output.temp | awk '{print $5}') ; rep=\"ATOMIC_POSITIONS (angstrom)\" ; ";
          comm+="rep2=\"$nats\\n\" ; ";
          comm+="grep -A $nats \"ATOMIC_POSITIONS\" output.temp   | sed \"s/$rep/$rep2/g\" | sed '/--/d' > pEli.xyz ; ";
          comm+="to=$(echo $nats+2 | bc ) ; t=$(echo $n_steps+1 | bc )  ;  ";
          comm+=" for ((i=1;i<$t;i++)); do head -$( echo $i*$to | bc ) pEli.xyz | tail -$to > aUxy$i.xyz ;  done ;";
          comm+=" grep \"!\" output.temp > energias ; mkdir "+dir;
          comm+=" ;  for ((i=1;i<$t;i++)); do energia=$(head -$i energias | tail -1 | awk '{print $5}' )  ; ";
          comm+=" sed \"2s/.*/$energia/\" aUxy$i.xyz > "+dir+"/coords$i.xyz ; done ; rm aUxy*.xyz energias  output.temp  pEli.xyz";
          //cout<<comm<<endl;
          ofstream dat("aUxi.sh");
          dat<<comm<<endl;
   system("chmod +x aUxi.sh ; bash aUxi.sh ; rm aUxi.sh ");
   dat.close();
}


double Energy_LJ(Atom a1, Atom a2)
{
   double sigma=(a1.R+a2.R)*0.8163; // Armstrongs
   double epsilon=1; // eV
   double Dist=Atomic_Distance(a1,a2);
   double energia = 4*epsilon*((  pow( (sigma/Dist) , 12 ) - pow( (sigma/Dist) , 6 ) ));
   return energia;
}

double Force_LJ(Atom a1, Atom a2,string calculation)
{
   double sigma=(a1.R+a2.R)*0.8163; // Armstrongs
   double epsilon;
   double f; // auxiliar
   double Dist=Atomic_Distance(a1,a2);
   if(calculation=="go")
   {
      epsilon=300; // eV
   }
   else if(calculation=="md")
   {
      epsilon=(a1.M+a2.M)/2; // eV
   }
   f=48*(epsilon/pow(Dist,2))*((pow( (sigma/Dist) , 12 ) - pow( (sigma/Dist) , 6 )/2 ));
   return f;
}


void Atomic_Structure::molecular_dynamic(Simulation_Cell box=default_box,  int nsteps=1000, string movie_file="false", float temperature=300, double cutoff_radii=10,  string force_field="LJ")
{
//http://phys.ubbcluj.ro/~tbeu/MD/C2_for.pdf
//According with:  Rapaport, D.C., \textit{The Art of Molecular Dynamics Simulation}, Cambridge University Press, Cambridge, 1995.
  int time=0;
  float deltaT=0.01;
  float Dist;
  float sep;
  double swap, distancia, dx, dy, dz;
  double fuerza;
  double xrange=box.M[0][0];
  double yrange=box.M[1][1];
  double zrange=box.M[2][2];
  double F[Nat][Nat];
  double a_aux[Nat][3];
  string command;
  // Para t=0 podemos generar un for que calcule las velocidades iniciales
  // if(Init_velocities==true){} else{ velocidades en 0}
  while(time < nsteps)
  {
     for(i=0;i<Nat;i++)
     {
        for(k=0;k<3;k++)
        {
           // Actualiza las posiciones
           atom[i].x[k] = atom[i].x[k] + atom[i].v[k]*deltaT + atom[i].a[k]*pow(deltaT,2)/2.0;
           if(box.periodicity==true) // AUN FALLAN LAS CONDICIONES PERIODICAS; NECESARIO IMPLEMENTAR COMO EN PACMAN
           {
               if(atom[i].x[k]>box.M[k][k])
               {
                  atom[i].x[k]=atom[i].x[k]-box.M[k][k];
               }
               else if(atom[i].x[k]<0)
               {
                 atom[i].x[k]=atom[i].x[k]+box.M[k][k];
               }
           }
        }
     }
     for(i=0;i<Nat;i++)
     {
        //CALCULA LAS FUERZAS
        for(j=i+1;j<Nat;j++) // Acá se puede hacer j=1 ; pues F[i][j]=-F[j][i]
        {
           if(box.periodicity==false)
           {
              distancia=Atomic_Distance(atom[i],atom[j]);
           }
           else // Periodic Boundary Conditions
           {
              dx=atom[i].x[0]-atom[j].x[0];
              dy=atom[i].x[1]-atom[j].x[1];
              dz=atom[i].x[2]-atom[j].x[2];
              if (abs(dx) > 0.5*xrange)
              {
                 if(dx>0){ dx = xrange + atom[j].x[0]-atom[i].x[0]; }
                 else { dx = xrange + atom[i].x[0]-atom[j].x[0]; }
              }
              if (abs(dy) > 0.5*yrange)
              {
                 if(dy>0){ dy = yrange + atom[j].x[1]-atom[i].x[1]; }
                 else { dy = yrange + atom[i].x[1]-atom[j].x[1]; }
              }
              if (abs(dz) > 0.5*xrange)
              {
                 if(dz>0){ dz = zrange + atom[j].x[2]-atom[i].x[2]; }
                 else { dz = zrange + atom[i].x[2]-atom[j].x[2]; }
              }
              distancia=sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2) ) ;
           }
           if(distancia<cutoff_radii)
           {
              if(force_field=="LJ")
              {
                 F[i][j]=Force_LJ(atom[i],atom[j],"md");
                 F[j][i]=F[i][j];
              }
              else if(force_field=="ML")
              {
                 // IDEM
                 cout<<"ML Force Field Not implemented yet ... coming soon!"<<endl;
                 //F[i][j]=Force_ML(atom[i],atom[j]);
                 //F[j][i]=-F[i][j];
              }
           }
           else // Quiere decir que quedan fuera del cutoff_radii
           {
              F[i][j]=0;
              F[j][i]=0;
           }
        }
     }
     // Esto calcula las aceleraciones del siguiente paso
     for(i=0;i<Nat;i++)
     {
        a_aux[i][0]=0.0; a_aux[i][1]=0.0; a_aux[i][2]=0.0;
        for(j=0;j<Nat;j++)
        {
           if(!(i==j))
           {
              if(force_field=="LJ")
              {
                 for(k=0;k<3;k++)
                 {
                    /////////////////////////////////////////////////
                    a_aux[i][k]=a_aux[i][k]+((F[i][j]/atom[i].M)*(atom[i].x[k]-atom[j].x[k]));
                    /////////////////////////////////////////////////
                 }
              }
              else if(force_field=="ML")
              {
                   // IDEM
                   cout<<"ML Force Field Not implemented yet ... coming soon!"<<endl;
                   //DEBE DEVOLVER a_aux[i][k]; i.e. las aceleraciones
              }
           }
        }
        // En esta parte se podría poner aux[i][k]=aux[i][k]+fuerza_externa[k]
        // if(Fuerza_Externa!=False):
        //      Aux[i,:]=Aux[i,:]+(Fuerza_Externa(e_a.atomo[i]))/e_a.atomo[i].M

     }
     // Calcula las nuevas velocidades
     for(i=0;i<Nat;i++)
     {
        for(k=0;k<3;k++)
        {
           atom[i].v[k] = atom[i].v[k] + (a_aux[i][k]+atom[i].a[k])*deltaT/2.0;
        }
     }
     // Intercambiamos las aceleraciones anteriores por las nuevas
     for(i=0;i<Nat;i++)
     {
        for(k=0;k<3;k++)
        {
           atom[i].a[k] = a_aux[i][k];
        }
     }

     //IMPRIME LAS COORDENADAS PARA HACER PELICULA
     if(!(movie_file=="false"))
     {
        print_xyz("temp.xyz");
        command="cat temp.xyz >> "+movie_file;
        system(command.c_str());
        system("rm temp.xyz");
     }
     //IMPRIME INFORMACION DEL SISTEMA
     cout<<"Step : "<<time<<"/"<<nsteps<<endl;
     if((time+1)>nsteps)
     {
        cout<<" Reached maximum number of steps ... Stopping simulation "<<endl;
     }
     time++;
  }
}


/*

void Atomic_Structure::geometry_optimization(string movie_file="false",int max_steps=1000, Simulation_Cell box=default_box, double tolerance=0.001, double cutoff_radii=10, int optimization_algorithm=0, string force_field="LJ")
{ // optimization_algorithm = 0 : gradient_descent ; 1 : conjugate_gradient ; 2 : quasi_newton
  // gradient_descent is recomended only for Lennard-Jones Potential, conjugate_gradient and quasi_newton might be more suitable for other potentials
   int current_step=0;
   double dE=100.0; // Diferencia energética ; valor inicial irrelevante;
   string command;
   double fuerza, energia=-100000, last_energy;
   double eta=0.01;
   cout<<" Entering geometry optimization loop ... "<<endl;

   while(current_step<max_steps && abs(dE)>tolerance)
   {
      last_energy=energia;
      energia=0;
      switch (optimization_algorithm)
      {
         case 0: ///// / / / /  /  /  /   /    /    /   /  / / / / / ///////  "gradient_descent"
            for(i=0;i<Nat;i++)
            {
               atom[i].a[0]=0.0; atom[i].a[1]=0.0; atom[i].a[2]=0.0;
               for(j=0;j<Nat;j++)
               {
                  if(Atomic_Distance(atom[i],atom[j])<cutoff_radii)
                  {
                     if(i!=j)
                     {
                        for(k=0;k<3;k++)
                        {
                           // CALCULA LAS FUERZAS
                           if(force_field=="LJ")
                           {
                              fuerza=Force_LJ(atom[i],atom[j],"go");
                           }
                         //else if(force_field=="ML")
                           //{
                            //  fuerza=Force_ML(atom[i],atom[j]);
                           //}
                           atom[i].a[k]=atom[i].a[k]+fuerza*(atom[i].x[k]-atom[j].x[k]);
                        }
                        energia=Energy_LJ(atom[i],atom[j])+energia;
                     }
                  }
               }
            }
            //CALCULA Diferencias energéticas
            dE=last_energy-energia;
            //ACTUALIZA LAS  POSICIONES
            for(i=0;i<Nat;i++)
            {
               for(k=0;k<3;k++)
               {
                  if(abs(atom[i].a[k])>1.0)
                  {
                     atom[i].x[k]=atom[i].x[k]+0.005*((atom[i].a[k])/(abs(atom[i].a[k])));
                  }
                  else
                  {
                     atom[i].x[k]=atom[i].x[k]+(eta*atom[i].a[k]);
                     //atom[i].x[k]=atom[i].x[k]+((eta*atom[i].a[k])/(1*abs(eta*atom[i].a[k])));
                  }
               }
            }
            break;
         //case "conjugated_gradient":
          //  ;
          //  break;
      }//end switch
      //IMPRIME LAS COORDENADAS PARA HACER PELICULA
      if(!(movie_file=="false"))
      {
         print_xyz("temp.xyz");
         command="cat temp.xyz >> "+movie_file;
         system(command.c_str());
         system("rm temp.xyz");
      }
      //IMPRIME INFORMACION DEL SISTEMA
    //  if(current_step>0)
      //{
         cout<<"Step : "<<current_step<<"/"<<max_steps<<"  Current energy : "<<energia<<" dE : "<<dE<<endl;
    //  }
      if((current_step+1)>max_steps)
      {
         cout<<" Reached maximum number of steps ... Stopping minimization "<<endl;
      }
      if(abs(dE)<tolerance)
      {
         cout<<" Reached energy accuracy ... Stopping minimization "<<endl;
      }
      current_step++;
   } //end while
} //end function
*/

void Atomic_Structure::geometry_optimization(string movie_file="false",int max_steps=1000, Simulation_Cell box=default_box, double tolerance=0.001, double cutoff_radii=10, int optimization_algorithm=0, string force_field="LJ")
{ // optimization_algorithm = 0 : gradient_descent ; 1 : conjugate_gradient ; 2 : quasi_newton
   int current_step=0;
   double dE=100.0; // Diferencia energética ; valor inicial irrelevante;
   double F[Nat][Nat], fuerza, energia=-100000, last_energy;
   double eta=0.001;
   float Dist;
   float sep;
   double swap, distancia, dx, dy, dz;
   double xrange=box.M[0][0];
   double yrange=box.M[1][1];
   double zrange=box.M[2][2];
   double a_aux[Nat][3];
   string command;
   cout<<" Entering geometry optimization loop ... "<<endl;
   while(current_step<max_steps && abs(dE)>tolerance)
   {
      last_energy=energia;
      switch (optimization_algorithm)
      {
         case 0: ///// / / / /  /  /  /   /    /    /   /  / / / / / ///////  "gradient_descent"
         for(i=0;i<Nat;i++)
         {
         //CALCULA LAS FUERZAS
            for(j=i+1;j<Nat;j++) // Acá se puede hacer j=1 ; pues F[i][j]=-F[j][i]
            {
               if(box.periodicity==false)
               {
                  distancia=Atomic_Distance(atom[i],atom[j]);
               }
               else // Periodic Boundary Conditions
               {
                  dx=atom[i].x[0]-atom[j].x[0];
                  dy=atom[i].x[1]-atom[j].x[1];
                  dz=atom[i].x[2]-atom[j].x[2];
                  if (abs(dx) > 0.5*xrange)
                  {
                     if(dx>0){ dx = xrange + atom[j].x[0]-atom[i].x[0]; }
                     else { dx = xrange + atom[i].x[0]-atom[j].x[0]; }
                  }
                  if (abs(dy) > 0.5*yrange)
                  {
                     if(dy>0){ dy = yrange + atom[j].x[1]-atom[i].x[1]; }
                     else { dy = yrange + atom[i].x[1]-atom[j].x[1]; }
                  }
                  if (abs(dz) > 0.5*xrange)
                  {
                     if(dz>0){ dz = zrange + atom[j].x[2]-atom[i].x[2]; }
                     else { dz = zrange + atom[i].x[2]-atom[j].x[2]; }
                  }
                  distancia=sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2) ) ;
               }
               if(distancia<cutoff_radii)
               {
                  if(force_field=="LJ")
                  {
                     F[i][j]=Force_LJ(atom[i],atom[j],"go");
                     F[j][i]=F[i][j];
                  }
                  else if(force_field=="ML")
                  {
                     // IDEM
                     cout<<"ML Force Field Not implemented yet ... coming soon!"<<endl;
                     //F[i][j]=Force_ML(atom[i],atom[j]);
                     //F[j][i]=-F[i][j];
                  }
               }
               else // Quiere decir que quedan fuera del cutoff_radii
               {
                  F[i][j]=0;
                  F[j][i]=0;
               }
            }
          }
          // Esto calcula las aceleraciones del siguiente paso
          energia=0;
          for(i=0;i<Nat;i++)
          {
               atom[i].a[0]=0.0; atom[i].a[1]=0.0; atom[i].a[2]=0.0;
               for(j=0;j<Nat;j++)
               {
                  if(!(i==j))
                  {
                     if(force_field=="LJ")
                     {
                        for(k=0;k<3;k++)
                        {
                           /////////////////////////////////////////////////
                           atom[i].a[k]=atom[i].a[k]+((F[i][j])*(atom[i].x[k]-atom[j].x[k]));
                           /////////////////////////////////////////////////
                        }
                        energia=Energy_LJ(atom[i],atom[j])+energia;
                     }
                     else if(force_field=="ML")
                     {
                        // IDEM
                        cout<<"ML Force Field Not implemented yet ... coming soon!"<<endl;
                        //DEBE DEVOLVER a_aux[i][k]; i.e. las aceleraciones
                     }
                  }
               }
            }
            //CALCULA Diferencias energéticas
            dE=last_energy-energia;
            //ACTUALIZA LAS  POSICIONES
            for(i=0;i<Nat;i++)
            {
               for(k=0;k<3;k++)
               {
                  if(abs(atom[i].a[k])>1.0)
                  {
                     atom[i].x[k]=atom[i].x[k]+0.005*((atom[i].a[k])/(abs(atom[i].a[k])));
                  }
                  else
                  {
                     atom[i].x[k]=atom[i].x[k]+(eta*atom[i].a[k]);
                     //atom[i].x[k]=atom[i].x[k]+((eta*atom[i].a[k])/(1*abs(eta*atom[i].a[k])));
                  }
               }
            }
            break;
//         case "conjugated_gradient":
  //          ;
    //        break;
      }//end switch
      //IMPRIME LAS COORDENADAS PARA HACER PELICULA
      if(!(movie_file=="false"))
      {
         print_xyz("temp.xyz");
         command="cat temp.xyz >> "+movie_file;
         system(command.c_str());
         system("rm temp.xyz");
      }
      //IMPRIME INFORMACION DEL SISTEMA
    //  if(current_step>0)
      //{
         cout<<"Step : "<<current_step<<"/"<<max_steps<<"  Current energy : "<<energia<<" dE : "<<dE<<endl;
    //  }
      if((current_step+1)>max_steps)
      {
         cout<<" Reached maximum number of steps ... Stopping minimization "<<endl;
      }
      if(abs(dE)<tolerance)
      {
         cout<<" Reached energy accuracy ... Stopping minimization "<<endl;
      }
      current_step++;
   } //end while
} //end function


void Atomic_Structure::initialize_electronic_density(int initialization=0, int charge=0)
{
  // initialization == 0 : random ; initialization == 1 : from file ; et al
   total_charge=0;
   int atomic_charge, _conter=0;
   //Counts total electric charge, sum each atomic charge
   for(i=0;i<Nat;i++)
   {
      atomic_charge=atom[i].Z+atom[i].charge;
      total_charge=total_charge+atomic_charge;
   }
   total_charge=total_charge+charge;
   cout<<" carga total "<<total_charge<<endl;
   //Genera los arreglos
   switch (initialization) {
     case 0:
        for(i=0;i<Nat;i++)
        {
           atomic_charge=atom[i].Z+atom[i].charge;
           atom[i].electron=new Electron[atomic_charge+1];
           for(j=0;j<atomic_charge;j++)
           {
              for(k=0;k<3;k++)
              {
                 atom[i].electron[j].x[k]=atom[i].x[k]+random_number(0,atom[i].R);
                 cout<<atom[i].electron[j].x[k]<<" ";
              }
              cout<<endl;
              _conter++;
              if(total_charge<_conter)
              {
                 break;
              }
           }
        }
        break;
    /*  case 1:
      ifstream coordinates_file(file);
      coordinates_file>>atom[i].electron[j].x[0]>>atom[i].electron[j].x[1]>>atom[i].electron[j].x[2];*/
   }
}



void Cluster::simulated_annealing(float kick_s=1.5,string movie_file="false",int max_steps=1000, double tolerance=0.00001, double cutoff_radii=10, int optimization_algorithm=0, string force_field="LJ")
{
  print_xyz("temporal.xyz");
  for(i=0;i<Nat;i++)
  {
     atom[i].x[0]=atom[i].x[0]*kick_s+random_number(-1,1)*kick_s/2;
     atom[i].x[1]=atom[i].x[1]*kick_s+random_number(-1,1)*kick_s/2;
     atom[i].x[2]=atom[i].x[2]*kick_s+random_number(-1,1)*kick_s/2;
  }
  geometry_optimization("temporal2.xyz");
  if(!(movie_file=="false"))
  {
     string command="cat temporal.xyz > "+movie_file;
            command+=" ; cat temporal2.xyz >> "+movie_file;
     system(command.c_str());
  }
  system("rm temporal.xyz temporal2.xyz");
}
