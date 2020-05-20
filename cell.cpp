a,b,c, longitudes de ejes
cosAB=coseno entre a y b (gamma)
cosAC=coseno entre a y c (beta)
cosBC=coseno entre b y c (alpha)
//Definir estas cosas en términos de celldm
/********************************************/

2  cubic F (fcc)
//v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)

v1[0]=-a/2;  v1[1]=0.0;  v1[2]=a/2;
v2[0]=0.00;  v2[1]=a/2;  v2[2]=a/2;
v3[0]=-a/2;  v3[1]=a/2;  v3[2]=0.0;

3 cubic I (bcc)
  //    v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)
v1[0]= a/2;  v1[1]= a/2;  v1[2]=a/2;
v2[0]=-a/2;  v2[1]= a/2;  v2[2]=a/2;
v3[0]=-a/2;  v3[1]=-a/2;  v3[2]=a/2;
-3  cubic I (bcc), more symmetric axis:
//      v1 = (a/2)(-1,1,1), v2 = (a/2)(1,-1,1),  v3 = (a/2)(1,1,-1)
v1[0]=-a/2;  v1[1]= a/2;  v1[2]= a/2;
v2[0]= a/2;  v2[1]=-a/2;  v2[2]= a/2;
v3[0]= a/2;  v3[1]= a/2;  v3[2]=-a/2;
4  Hexagonal and Trigonal P        celldm(3)=c/a
//      v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)
v1[0]= a;  v1[1]= 0.00;  v1[2]= 0.00;
v2[0]=-a/2;  v2[1]=a*sqrt(3.0)/2.0;  v2[2]=0.0;
v3[0]= 0;  v3[1]= 0;  v3[2]=celldm3;
5       Trigonal R, 3fold axis c        celldm(4)=cos(gamma)
/*      The crystallographic vectors form a three-fold star around
      the z-axis, the primitive cell is a simple rhombohedron:
      v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
      where c=cos(gamma) is the cosine of the angle gamma between
      any pair of crystallographic vectors, tx, ty, tz are:*/
tx=sqrt((1-c)/2.0); ty=sqrt((1-c)/6); tz=sqrt((1+2c)/3);
v1[0]= a*tx;  v1[1]= -a*ty;  v1[2]= a*tz;
v2[0]=0;  v2[1]=2*a*ty;  v2[2]=a*tz;
v3[0]= -a*tx;  v3[1]= -a*ty;  v3[2]=a*tz;
-5         Trigonal R, 3fold axis <111>    celldm(4)=cos(gamma)
/*      The crystallographic vectors form a three-fold star around
      <111>. Defining a' = a/sqrt(3) :
      v1 = a' (u,v,v),   v2 = a' (v,u,v),   v3 = a' (v,v,u)
      where u and v are defined as
        u = tz - 2*sqrt(2)*ty,  v = tz + sqrt(2)*ty
      and tx, ty, tz as for case ibrav=5
      Note: if you prefer x,y,z as axis in the cubic limit,
            set  u = tz + 2*sqrt(2)*ty,  v = tz - sqrt(2)*ty
            See also the note in Modules/latgen.f90 */
tx=sqrt((1-c)/2.0); ty=sqrt((1-c)/6); tz=sqrt((1+2c)/3);
ap = a/sqrt(3.0);
u=tz-2*sqrt(2.0)*ty; v=tz+sqrt(2.0)*ty;
v1[0]= ap*u;  v1[1]= ap*v;  v1[2]= ap*v;
v2[0]=ap*v;  v2[1]=ap*u;  v2[2]=ap*v;
v3[0]= ap*v;  v3[1]= ap*v;  v3[2]=ap*u;

6          Tetragonal P (st)               celldm(3)=c/a
//      v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)
v1[0]= a;  v1[1]= 0.0;  v1[2]= 0.0;
v2[0]=0.0;  v2[1]=a;  v2[2]=0.0;
v3[0]=0.0;  v3[1]= 0.0;  v3[2]=a*celldm3;
7          Tetragonal I (bct)              celldm(3)=c/a
//      v1=(a/2)(1,-1,c/a),  v2=(a/2)(1,1,c/a),  v3=(a/2)(-1,-1,c/a)
v1[0]= a/2;  v1[1]=-a/2;  v1[2]=celldm3*(a/2.0);
v2[0]= a/2;  v2[1]= a/2;  v2[2]=celldm3*(a/2.0);
v3[0]=-a/2;  v3[1]=-a/2;  v3[2]=celldm3*(a/2.0);

  8          Orthorhombic P                  celldm(2)=b/a celldm(3)=c/a
//      v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)
v1[0]= a;  v1[1]= 0.0;  v1[2]= 0.0;
v2[0]=0.0;  v2[1]=b;  v2[2]=0.0;
v3[0]=0.0;  v3[1]= 0.0;  v3[2]=c;

  9          Orthorhombic base-centered(bco) celldm(2)=b/a  celldm(3)=c/a
//      v1 = (a/2, b/2,0),  v2 = (-a/2,b/2,0),  v3 = (0,0,c)
v1[0]= a/2;  v1[1]= b/2;  v1[2]= 0.0;
v2[0]=-a/2;  v2[1]=b/2;  v2[2]=0.0;
v3[0]=0.0;  v3[1]= 0.0;  v3[2]=c;

 -9          as 9, alternate description
//      v1 = (a/2,-b/2,0),  v2 = (a/2, b/2,0),  v3 = (0,0,c)
v1[0]= a/2;  v1[1]=-b/2;  v1[2]= 0.0;
v2[0]=a/2;  v2[1]=b/2;  v2[2]=0.0;
v3[0]=0.0;  v3[1]= 0.0;  v3[2]=c;

 91          Orthorhombic one-face base-centered A-type      celldm(2)=b/a   celldm(3)=c/a
//      v1 = (a, 0, 0),  v2 = (0,b/2,-c/2),  v3 = (0,b/2,c/2)
v1[0]= a;  v1[1]= 0.0;  v1[2]= 0.0;
v2[0]=0.0;  v2[1]=b/2;  v2[2]=-c/2;
v3[0]=0.0;  v3[1]=b/2;  v3[2]=c/2;

 10          Orthorhombic face-centered      celldm(2)=b/a     celldm(3)=c/a
//      v1 = (a/2,0,c/2),  v2 = (a/2,b/2,0),  v3 = (0,b/2,c/2)
v1[0]= a/2;  v1[1]=0v1[0]= a;  v1[1]=0.0;  v1[2]=0.0;
v2[0]=b*cos(gamma); v2[1]=b*sin(gamma);  v2[2]=0.0;
v3[0]=0.0;  v3[1]=0.0;  v3[2]=c;.0;  v1[2]= c/2;
v2[0]=a/2;  v2[1]=b/2;  v2[2]=0.0;
v3[0]=0.0;  v3[1]= b/2;  v3[2]=c/2;

 11          Orthorhombic body-centered      celldm(2)=b/a    celldm(3)=c/a
//      v1=(a/2,b/2,c/2),  v2=(-a/2,b/2,c/2),  v3=(-a/2,-b/2,c/2)
v1[0]= a/2;  v1[1]= b/2;  v1[2]= c/2;
v2[0]=-a/2;  v2[1]=b/2;  v2[2]=c/2;
v3[0]=-a/2;  v3[1]=-b/2;  v3[2]=c/2;

 12          Monoclinic P, unique axis c     celldm(2)=b/a   celldm(3)  celldm(4)=cos(ab)
//      v1=(a,0,0), v2=(b*cos(gamma),b*sin(gamma),0),  v3 = (0,0,c) gamma is angle between a and b.
gamma=acos(celldm4);
v1[0]= a;  v1[1]=0.0;  v1[2]=0.0;
v2[0]=b*cos(gamma); v2[1]=b*sin(gamma);  v2[2]=0.0;
v3[0]=0.0;  v3[1]=0.0;  v3[2]=c;

-12          Monoclinic P, unique axis b     celldm(2)=b/a    celldm(3)=c/a,     celldm(5)=cos(ac)
//      v1 = (a,0,0), v2 = (0,b,0), v3 = (c*cos(beta),0,c*sin(beta))  where beta is the angle between axis a and c
beta=acos(celldm5);
v1[0]= a;  v1[1]=0.0;  v1[2]=0.0;
v2[0]=0.0; v2[1]=b;  v2[2]=0.0;
v3[0]=c*cos(beta);  v3[1]=0.0;  v3[2]=c*sin(beta);

 13          Monoclinic base-centered   celldm(2)=b/a  (unique axis c)  celldm(3)=c/a,    celldm(4)=cos(gamma)
/*      v1 = (  a/2,         0,          -c/2),
        v2 = (b*cos(gamma), b*sin(gamma), 0  ),
        v3 = (  a/2,         0,           c/2), */ gamma=angle between a and b projected on xy plane

-13          Monoclinic base-centered
  /*                                         celldm(2)=b/a
             (unique axis b)                 celldm(3)=c/a,
                                             celldm(5)=cos(beta)
      v1 = (  a/2,      -b/2,             0),
      v2 = (  a/2,       b/2,             0),
      v3 = (c*cos(beta),   0,   c*sin(beta)),
      where beta=angle between axis a and c projected on xz plane*/

 14          Triclinic
                    /*                       celldm(2)= b/a,
                                             celldm(3)= c/a,
                                             celldm(4)= cos(bc),
                                             celldm(5)= cos(ac),
                                             celldm(6)= cos(ab)
      v1 = (a, 0, 0),
      v2 = (b*cos(gamma), b*sin(gamma), 0)
      v3 = (c*cos(beta),  c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma),
           c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
                     - cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) )
      where alpha is the angle between axis b and c
             beta is the angle between axis a and c
            gamma is the angle between axis a and b*/

gamma=acos(celldm4);
beta=acos(celldm5);
alpha=acos(celldm6):
v1[0]= a;  v1[1]=0.0;  v1[2]=0.0;
v2[0]=b*cos(gamma); v2[1]=b*sin(gamma);  v2[2]=0.0;
v3[0]=c*cos(beta);  v3[1]= c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
v3[2]=c*sqrt(1+2*cos(alpha)*cos(beta)*cos(gamma)-cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) );
