/* TVD scheme */
#include <stdio.h>
#include <math.h>

int main() {
    char name[50];
    double ap,an,as,aw,ae,aww,ass,su;
    double min1,min2,min3,min4,min5,min6,min7,min8,max1,max2,max3,max4;
    int g = 1;
    int s = 3 ;
    double w = 0.1;    
    int iteration = 70000;
    double ROE = 1e-13; 
    /* x=0 2nd upwind TVD (2r) ,x=1 QUICK TVD(2r) ,x=2 2nd upwind (r) TVD ,x=3 QUICK TVD (r) ,x=4 CDS , x=5 upwind , x=6 2nd upwind , x=7 QUICK*/
    int x = 5;
    int steps = 0;  
    double c_diff;    
    int N_[2] ={42,82};
    double h_[2] ={0.025,0.0125};
    double T_exact[N_[g]][N_[g]] = {0};
    double T[N_[g]][N_[g]] = {0};
    double T_new[N_[g]][N_[g]] = {0};
    int Pe[4] = {1,10,100,100000000} ;
    FILE*data2 = fopen( "steps_error.dat" ,"w");
    fprintf(data2,"VARIABLES = \"steps\"\"error\"\n");
    fprintf(data2,"ZONE T=\"UPWIND \" F=POINT\n I = 2000\n");
    for (int i = 0; i < N_[g]; i++) 
    {
        for (int j = 0; j < 2; j++)
        {
            T[i][j] = 1;
            T[i][N_[g] - 1] = 0;    
            T[j][i] = 0;
            T[N_[g] - 1][i] = 1;   
        }   
    }
    
   do 
    {
        c_diff = 0 ;
        for (int i = 2 ; i < N_[g]-1 ; i++)
        {
            for (int j = 2 ; j < N_[g]-1; j++)
            {
                int a1,a2,a3,a4;
                double CDS_N=(T[i+1][j]-T[i][j])/(T[i][j]-T[i-1][j]+1E-10);
                double CDS_S=(T[i][j]-T[i-1][j])/(T[i-1][j]-T[i-2][j]+1E-10);
                double CDS_W=(T[i][j]-T[i][j-1])/(T[i][j-1]-T[i][j-2]+1E-10);
                double CDS_E=(T[i][j+1]-T[i][j])/(T[i][j]-T[i][j-1]+1E-10);
                
                double sinr_N[5] ={CDS_N,0,1,0.75*CDS_N+0.25,2*CDS_N};
                double sinr_S[5] ={CDS_S,0,1,0.75*CDS_S+0.25,2*CDS_S};
                double sinr_W[5] ={CDS_W,0,1,0.75*CDS_W+0.25,2*CDS_W};
                double sinr_E[5] ={CDS_E,0,1,0.75*CDS_E+0.25,2*CDS_E};
                if (x == 0)
                {
                    min1 = fmin(2*CDS_N,1);
                    max1 = fmax(0,min1);
                    min2 = fmin(2*CDS_S,1);
                    max2 = fmax(0,min2);
                    min3 = fmin(2*CDS_W,1);
                    max3 = fmax(0,min3);
                    min4 = fmin(2*CDS_E,1);
                    max4 = fmax(0,min4);
                }
                else if (x == 1)
                {
                    min1 = fmin(2*CDS_N,1);
                    min2 = fmin(0.75*CDS_N+0.25,min1);
                    max1 = fmax(0,min2);
                    min3 = fmin(2*CDS_S,1);
                    min4 = fmin(0.75*CDS_S+0.25,min3);
                    max2 = fmax(0,min4);
                    min5 = fmin(2*CDS_W,1);
                    min6 = fmin(0.75*CDS_W+0.25,min5);
                    max3 = fmax(0,min6);
                    min7 = fmin(2*CDS_E,1);
                    min8 = fmin(0.75*CDS_E+0.25,min7);
                    max4 = fmax(0,min8);   
                }
                else if (x == 2)
                {
                    min1 = fmin(CDS_N,1);
                    max1 = fmax(0,min1);
                    min2 = fmin(CDS_S,1);
                    max2 = fmax(0,min2);
                    min3 = fmin(CDS_W,1);
                    max3 = fmax(0,min3);
                    min4 = fmin(CDS_E,1);
                    max4 = fmax(0,min4);
                }
                else if (x == 3 )
                {
                    min1 = fmin(CDS_N,1);
                    min2 = fmin(0.75*CDS_N+0.25,min1);
                    max1 = fmax(0,min2);
                    min3 = fmin(CDS_S,1);
                    min4 = fmin(0.75*CDS_S+0.25,min3);
                    max2 = fmax(0,min4);
                    min5 = fmin(CDS_W,1);
                    min6 = fmin(0.75*CDS_W+0.25,min5);
                    max3 = fmax(0,min6);
                    min7 = fmin(CDS_E,1);
                    min8 = fmin(0.75*CDS_E+0.25,min7);
                    max4 = fmax(0,min8);                    
                }
                else if (x = 4)
                {
                    a1 = a2 = a3 = a4 =0 ;
                }
                else if (x = 5)
                {
                    a1 = a2 = a3 = a4 =1 ;
                }                
                else if (x = 6)
                {
                    a1 = a2 = a3 = a4 =2 ;
                }
                else if (x = 7)
                {
                    a1 = a2 = a3 = a4 =3 ;
                }  
                     
                for (int m = 0; m < 6; m++)
                {
                    if (max1 == sinr_N[m])
                    {
                        a1 = m ;
                    }
                    if (max2 == sinr_S[m])
                    {
                        a2 = m ;
                    }
                    if (max3 == sinr_W[m])
                    {
                        a3 = m ;
                    }
                    if (max4 == sinr_E[m])
                    {
                        a4 = m ;
                    }                        
                } 
              
                ap = Pe[s]+ Pe[s]+ 4;
                ae = 1;
                an = 1;
                aw = Pe[s]+1;
                as = Pe[s]+1;
                su = Pe[s]*0.5*(sinr_E[a4]*(T[i][j]-T[i][j-1])-sinr_W[a3]*(T[i][j-1]-T[i][j-2])+sinr_N[a1]*(T[i][j]-T[i-1][j])-sinr_S[a2]*(T[i-1][j]-T[i-2][j]));
                T_new[i][j] = (ae*T[i][j+1]+an*T[i+1][j]+aw*T[i][j-1]+as*T[i-1][j]-su)/ap ;
                c_diff = c_diff + fabs(T_new[i][j]-T[i][j]);
                T[i][j] = T[i][j] + w*(T_new[i][j]-T[i][j]);
                
                  
            }
        }
        for (int i = 2; i < N_[g]-1; i++)
        {
            for (int j = 2; j < N_[g]-1; j++)
            {
                T[i][j] = T[i][j] + w*(T_new[i][j]-T[i][j]);           
            }    
        }
        

        if ( steps % 10 == 0 )
        {
            fprintf(data2,"%d\t%e\n",steps,c_diff);
        }
        steps++;
        c_diff = c_diff/((N_[g]-2)*(N_[g]-2));
    } while ( c_diff > ROE && steps < iteration);
    printf("iterations: %d\n", steps);
    printf ("%e",c_diff) ;
    fclose(data2);

        


    FILE *data_2d = fopen( "data_graph.dat", "w");
    fprintf(data_2d, "VARIABLES=\"X\",\"Y\",\"T\"\n");
    fprintf( data_2d, "ZONE T=\"STEP=%d\" F=POINT\n",steps);
    fprintf(data_2d, "I=%d, J=%d\n", N_[g]-1, N_[g]-1);
    sprintf(name , "data_T-Y.dat");
 

    int i, k;
    for (k =1; k < N_[g]; k++) {
        for (i =1; i < N_[g]; i++) {
            fprintf(data_2d, "%d\t%d\t%e\n", i, k, T[k][i]);
        }
    }

    fclose(data_2d); 

    FILE*data = fopen( name ,"w");
    fprintf(data,"VARIABLES = \"Y\"\"T\"\n");
    fprintf(data,"ZONE T=\"UPWIND \" F=POINT\n I = %d\n",(N_[g]-1));
    for (int i = 1; i < N_[g] ;  i++)
    {
        
        fprintf(data,"%lf\t%e\n",((i-1)*h_[g]),T[i][N_[g]/2]);
    }
    fclose(data);
    }




