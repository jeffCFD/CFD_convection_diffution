/* linear equation */
#include <stdio.h>
#include <math.h>

int main() {
    char name[50];
    int g = 0 ;
    int s = 3 ;
    double w = 1;    
    double a1 =0,a2=0,a3=0,a4=0 ;
    int iteration = 10000000;
    double ROE = 1e-13;  
    int steps = 0;  
    double c_diff; 
    double le ;   
    int N_[2] ={42,82};
    double h_[2] ={0.025,0.0125};
    double T_exact[N_[g]][N_[g]] = {0};
    double T[N_[g]][N_[g]] = {0};
    double T_new[N_[g]][N_[g]] = {0};
    int Pe[4] = {1,10,100,100000} ;
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

                int m = 1;
                switch (m)
                {
                case 0/* CDS */:
                    a1 = (1+2*Pe[s])*(T[i-1][j]+T[i][j-1]);
                    a2 = T[i+1][j]+T[i][j+1];
                    a3 = (Pe[s]/2.)*(T[i-2][j]+T[i][j-2]);
                    a4 = (4+3*Pe[s]);    
                    le = (a1+a2-a3)/a4 ;
                    break;

                case 1/* upwind*/:
                    if (Pe[s] = true)
                    {
                        le = (T[i][j-1]+T[i-1][j])/2 ;
                    }
                    else{
                        le = ((1+Pe[s])*T[i][j-1]+T[i][j+1]+T[i+1][j]+(1+Pe[s])*T[i-1][j])/(2*Pe[s]+4);
                    }
                        break;
                case 2/* 2nd order */:
                    if (Pe[s] = true)
                    {
                        a1 = (2/3.)*(T[i][j-1]+T[i-1][j]);
                        a2 = (1/6.)*(T[i-2][j]+T[i][j-2]);
                        le = a1-a2;
                    }
                    else
                    {
                        le = ((1-Pe[s]/2)*T[i][j+1]+(1+Pe[s]/2)*T[i][j-1]+(1-Pe[s]/2)*T[i+1][j]+(1+Pe[s]/2)*T[i-1][j])/4;
                    }
                        break;
                case 3/* QUICK */:
                    if (Pe[s] = true)
                    {
                        a1 = 7*(T[i][j-1]+T[i-1][j])/6.;
                        a2 = (T[i][j+1]+T[i+1][j])/2.;
                        a3 = (T[i][j-2]+T[i-2][j])/6.;                    
                        le = a1-a2-a3;
                    }
                    else
                    {
                        a1 = (1-3/8.*Pe[s])*(T[i][j+1]+T[i+1][j]);
                        a2 = (T[i-2][j]+T[i][j-2])*Pe[s]/8. ;
                        a3 = (1+7/8.*Pe[s])*(T[i][j-1]+T[i-1][j]) ; 
                        a4 = 4+(3/4.)*Pe[s] ;
                        le = (a1 - a2 + a3)/ a4;
                    }

                        break;  
                case 4/* downwind(2r) */:

                    le = ((1-Pe[s])*T[i][j+1]+T[i][j-1]+T[i-1][j]+(1-Pe[s])*T[i+1][j])/(4-2*Pe[s]) ;
                    break;                  
                }
                T_new[i][j] = le ; 
                c_diff = c_diff + fabs(T_new[i][j]-T[i][j]);
                T[i][j] = T[i][j] + w*(T_new[i][j]-T[i][j]);    
            }
        }
        steps++;
        c_diff = c_diff/((N_[g]-2)*(N_[g]-2));
    } while ( c_diff > ROE && steps < iteration);
    printf("iterations: %d\n", steps);
    printf ("%e",c_diff) ; 

        


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
    fprintf(data,"ZONE T=\"UPWIND \" F=POINT\n I = 41\n");
    for (int i = 1; i < N_[g] ;  i++)
    {
        
        fprintf(data,"%lf\t%e\n",((i-1)*h_[g]),T[i][N_[g]/2]);
    }
    fclose(data);
    }

