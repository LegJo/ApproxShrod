#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include <string.h>





//Equation de Schrod : f''(x)+(2m/hbar²)( E-V(x) )f(x)=0

/*

    dans  I[-a;a]:



    Potentiels V(x) Nul:

       0 : ∀x

    Potentiels V(x) en "rectangle" 

        V0 : -b>x>=b

        0 : sinon

    Potentiels V(x) en "marche":

        0 : -a>x>b

        V0 : b>x>a

    Potentiels V(x) en parabole:

        +inf : x<-b

        0 : -b <= x < -a

        -V0x²+V1 : -a<= x <= a

        0 : a< x <=b

        -inf : x>b

*/



long double E ;
long double V0 ; 
long double V1 ; 
int TypePot;
float a = 0;
float b = 0;

//double m=0.511; //on choisie ici de prendre la masse de l'electron en MeV pour simplifiez les calcules. la masse de l'électron de 9,1 10^-31kg correspond à une énergie de 0,511MeV.

//double hbar=197.326; //on prend hbar en MeV-fm


long double ldabs(long double a) {

    if(a < 0) 
      a=-a;
    return a;
}

long double ldCarre(long double y)
{
  return ((y)*(y));
}

void WriteTabforGnuplot(long double** tab, int n)
{
    FILE* tabfile = fopen("Tableau.tmp","w+");
    if (!tabfile) {
        exit(EXIT_FAILURE);
    }
    for (int i=0; i<=n; i++)
    {
        fprintf(tabfile, "%0.14Lf %0.14Lf\n", tab[i][0], tab[i][1]);
    }
    fclose(tabfile);
} 


void WriteFileTab(long double** tab, int n)
{

    FILE* tabfile = fopen("TableauDeValeurs.txt","w+");

    if (!tabfile) {
        exit(EXIT_FAILURE);
    }
    fprintf(tabfile, "        x           yn          zn          ∫yn² \n");
    for (int i=0; i<=n; i++)
    {
        fprintf(tabfile, "%0.9Le %0.9Le %0.9Le %0.9Le\n", tab[i][0], tab[i][1], tab[i][2], tab[i][3]);
    }
    fclose(tabfile);
} 


void TracerTab(int borneGraphInf, int borneGraphSup, long double z0, long double E, long double prec)

{

    char commande[256]={0};

    snprintf(commande, 256, "./TracerV2.bash %d %d %f %f %d %Lf %Lf %Lf %Lf %Lf", borneGraphInf, borneGraphSup, a, b, TypePot, V0, V1, z0, E, prec);

    system(commande);

}

long double f(int i, long double yn, long double zn, long double tn)
{
    //float k= 2*m / (hbar*hbar);
    /*
            On a f''(x) = (2m/hbar²)[V(x)-E]f(x) => f''(x) = [V(x)-E]f(x)
            et f(x)=y
            On pose y' = z
            z'= (2m/hbar²)[V(x)-E]y => z'= [V(x)-E]y 
            y'= z
    */
    if(i == 1)
    {
         return zn;
    }
    else
    {
        switch(TypePot)
        {
            //V(x)=0 Nul
            case 1 : return ((-E)*yn); break;
            //V(x)=V0 rectangulaire
            case 2 : if(tn > -b && tn < b){return ((V0-E)*yn);}else {return ((-E)*yn );}break;
            //marche  V(x)=V0 x>b et V(x)=0 x<b
            case 3 : if(tn < b){return ((-E)*yn);}else if(tn >= b){return ((V0-E)*yn );}break;
            //V(x)=(-V0*x²)+V1 Parabole
            case 4 : if(tn < -b || tn > b){ return ((-E)*yn);}else if(tn >= -b && tn <= b){return ((-V0*tn*tn+V1-E)*yn );}break;
            default : return ((-E)*yn ); 
        }
   }
}

long double An(int i, float h, long double yn, long double zn, float tn)
{
  return h * f(i, yn, zn, tn);
}

long double Bn(int i, float h, long double yn, long double zn, long double tn, long double A1, long double A2)
{
  return h * f(i, yn + A1/2, zn+ A2/2, tn + h/2);

}

long double Cn(int i, float h, long double yn, long double zn, long double tn, long double B1, long double B2)
{
  return h * f(i, yn + B1/2, zn + B2/2, tn + h/2);

}

long double Dn(int i, float h, long double yn, long double zn, long double tn, long double C1, long double C2)
{

  return h* f(i, yn + C1, zn + C2, tn + h);
}

long double** rk2(int N, float h, long double yn, long double zn, long double tn, long double** tabRK2)

{

  long double A1 = 0;
  long double A2 = 0;
  long double B1 = 0;
  long double B2 = 0;

  for(int i = 1; i <= N; i++)
  {

    A1 = An(1, h, yn, zn, tn);
    A2 = An(2, h, yn, zn, tn);
    B1 = Bn(1, h, yn, zn, tn, A1, A2);
    B2 = Bn(2, h, yn, zn, tn, A1, A2);
    long double yn1 = yn + B1;
    long double zn1 = zn + B2;
    yn = yn1;
    zn =zn1;
    tn += h;
    // printf("Yn : %f, Zn : %f, x: %f \n",  yn1, zn1, tn);

    tabRK2[i][0] = tn;
    tabRK2[i][1] = yn1;
    tabRK2[i][2] = zn1;
    tabRK2[i][3] += h*ldCarre(yn1) ;
  }
  return tabRK2;
}

long double** rk4(int N, float h, long double yn, long double zn, long double tn, long double** tabRK4)
{

  long double A1 = 0;
  long double A2 = 0;
  long double B1 = 0;
  long double B2 = 0;
  long double C1 = 0;
  long double C2 = 0;
  long double D1 = 0;
  long double D2 = 0;
  for(int i = 1; i <= N; i++)
  {

    A1 = An(1, h, yn, zn, tn);
    A2 = An(2, h, yn, zn, tn);
    B1 = Bn(1, h, yn, zn, tn, A1, A2);
    B2 = Bn(2, h, yn, zn, tn, A1, A2);
    C1 = Bn(1, h, yn, zn, tn, B1, B2);
    C2 = Bn(2, h, yn, zn, tn, B1, B2);
    D1 = Bn(1, h, yn, zn, tn, C1, C2);
    D2 = Bn(2, h, yn, zn, tn, C1, C2);

    long double yn1 = yn + (A1 + 2*B1 + 2*C1 + D1)/6;
    long double zn1 = zn + (A2 + 2*B2 + 2*C2 + D2)/6;


    yn = yn1;
    zn = zn1;
    tn += h;
    // printf("Yn : %f, Zn : %f, x: %f \n",  yn1, zn1, tn);


    tabRK4[i][0] = tn;
    tabRK4[i][1] = yn1;
    tabRK4[i][2] = zn1;
    tabRK4[i][3] += h*ldCarre(yn1) ;
  }
  return tabRK4;
} 

int main()
{
    int n;
    printf("Entrez nombre de points n à calculer :");
    scanf("%d",&n);
    int Rk;
    printf("\nEntrez la methode à utiliser \n 1 pour rk4 \n 2 pour rk2 \n=>");
    scanf("%d",&Rk);
    printf("\non considerera  pour simplifier les calculs que 2m/hbar²=1 \n");
    printf("Choississez le type de potentiel V(x): \n➔ 1 pour Potentiel Nul [V(x)=0 ∀ x∈ I] \n➔ 2 pour Potentiel en rectangle [V(x)=V0 -b<x<b et V(x)=0 sinon]"); 
    printf("\n➔ 3 pour Potentiel en marches [V(x)=0 -a<x<=b et V(x)=V0 b<x<a] \n➔ 4 pour Potentiel en parabole \n");
    printf("[+inf : x<-a, 0 : -a <= x < -b, -V0x²+V1 : -b<= x <= b, 0 : b< x <=a,-inf : x>a]\n=>");
    scanf("%d",&TypePot);
    printf("\nEntrez la valeur 'a' definissant l'interval I=[-a;a], Intervale où la Fonction sera approximé :\n=>");
    scanf("%f",&a);
    if(TypePot != 1)
    {  
       printf("\nEntrez la valeur de 'b' compris dans [-a;a]:\n=>");
       scanf("%f",&b);
    }
    
    printf("\nEntrez la valeur de 'V0' :\n=>");
    scanf("%Lf",&V0);
    if(TypePot ==4)
    {
       printf("\nEntrez la valeur de 'V1':\n=>");
       scanf("%Lf",&V1);
    }
    long double y0 = 0;
    float h = ((float)a+(float)a)/((float)n-1) ; //calcul du pas, h 
    long double** Tab  = (long double**) malloc((n+1)*sizeof(long double*));
    for (int i=0;i<=n;i++)// Réservation de l’espace mémoire de chaque ligne
    {
        Tab[i]=(long double*)malloc(4*sizeof(long double));
    }


    
    float bornZ;
    printf("\nEntrez la valeur x (x>0) correspondant à l'écart de depart de recherche pour E \navec E qui variera entre E=V0-x et E=V0 au max:\n=>");
    scanf("%f",&bornZ);
    long double z0;
    printf("\nEntrez la valeur de depart de de z0 :\n=>");
    scanf("%Lf",&z0);
    long double searchStep =0.05;
    long double precision = 0.05;
    
    int bool = 0;
    while(!bool)
    {
        for(E=V0-bornZ;E<V0;E+=searchStep)  //E=-bornE;E<bornE;E+=searchStep     long double i = -a; i < a; i += h
        {
          if(E != 0)
          {
            printf("z0=%0.10Lf E=%0.10Lf\n",z0,E);
            Tab[0][0] = -a; //x
            Tab[0][1] = y0; //f(x)
            Tab[0][2] = z0; //f'(x)
            Tab[0][3] = 0; // ∫f(x)² entre x et x+h
            if(Rk == 1)
                Tab = rk4(n, h, y0 , z0 , -a, Tab );
            else 
                Tab = rk2(n, h, y0 , z0 , -a, Tab );
            if (Tab[n][1]< precision && Tab[n][1]> -precision )
            {
                bool = 1;
                printf("--------------------\nz0=%0.10Lf E=%0.10Lf\n",z0,E);
                printf("f(%f)=%Lf\n",a ,Tab[n][1]);
                printf("∫f(x)² entre %f et %f :%Lf\n",-a , a,Tab[n][3]);
                break;

            }
          }
        }
        z0+= searchStep;  
        if(z0 < 0.3 && z0>-0.3)
          z0 +=0.3-z0;
        if(z0 > 100)
            bool=1;
    }
    WriteFileTab(Tab,n);
    WriteTabforGnuplot(Tab,n);
    //system("cat TableauDeValeurs.txt")
    TracerTab((-a-1), (a+1), z0, E, bornZ);
    return 0 ;
  }