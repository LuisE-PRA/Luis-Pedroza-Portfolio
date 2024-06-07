#include <stdio.h>
#include <stdlib.h>

double genera_vecinos(int I,int J,int K,int M,int N,double HS[I][K],double HD[J][K])
{
    int i,k,j,k1,count = 0;
    int gli, glj;
    int A[I][K][J][K];
//    double q;
    for (i = 0; i < I; i++){
        for (k = 0; k < K; k++){
            for (j = 0; j < J; j++){
                for (k1 = 0; k1 < K; k1++){
                    if (HD[j][k1] > 0 && HS[i][k] > 0 && k1 >= k){
                        A[i][k][j][k1] = 1; ++count;
                    }else{A[i][k][j][k1] = 0;}
                    printf("%i %i %i %i %i\n",i,k,j,k1,A[i][k][j][k1]);
                }
            }
        }
    }

    gli = count - M;
    count = 0;
    for (j = 0; j < J; j++){
        for (k1 = 0; k1 < K; k1++){
            for (i = 0; i < I; i++){
                for (k = 0; k < K; k++){
                    if (HS[i][k] > 0 && HD[j][k1] > 0 && k1 >= k){
                        A[i][k][j][k1] = 1; ++count;
                    }else{A[i][k][j][k1] = 0;}
                }
            }
        }
    }
    glj = count - N;
 //printf("%i %i", gli,glj);
    return(0);
}


int main()
{
    int I, J, K, M, N;
    double qrec, qhu, qcu;
    FILE *f;
    //LLama a GAMS
    system("gams 06MT lo=2");
    //
    //Abre el archivo Parametros//
    f = fopen("Parametros.c","rt");

        fscanf(f,"%i",&I); /*Conjunto de las corrientes calientes */
        fscanf(f,"%i",&J); /*Conjunto de las corrientes frias */
        fscanf(f,"%i",&K); /*Conjunto de intervalos de temperatura */
        fscanf(f,"%lg",&qrec); /*Calor recuperado*/
        fscanf(f,"%lg",&qhu); /*Calor requerido para servicios de calentamiento*/
        fscanf(f,"%lg",&qcu); /*Calor requerido para servicios de enfriamiento*/
        fscanf(f,"%i",&M); /*Intervalos de entalpia calientes presentes*/
        fscanf(f,"%i",&N); /*Intervalos de entalpia frios presentes*/
    fclose(f);


    double HS[I][K];
    double HD[J][K];
    int i[M],k[M],j[N],k1[N],r,s;
    double hs[M],hd[N];
    FILE *g;

    g = fopen("HS.c","rt");
//HS
    for (r = 0; r < M; r++ ) {
        fscanf(g,"%i",&i[r]);
    }
    for (r = 0; r < M; r++ ) {
        fscanf(g,"%i",&k[r]);
    }
    for (r = 0; r < M; r++ ) {
        fscanf(g,"%lg",&hs[r]);
    }
    for (r = 0; r < I; r++ ) {
        for (s = 0; s < K; s++ ){
            HS[r][s] = 0;
        }
    }
    for (r = 0; r < M; r++ ) {
            HS[i[r]-1][k[r]-1] = hs[r];
    }
    fclose(g);
//HD
    g = fopen("HD.c","rt");
    for (r = 0; r < N; r++ ) {
        fscanf(g,"%i",&j[r]);
    }
    for (r = 0; r < N; r++ ) {
        fscanf(g,"%i",&k1[r]);
    }
    for (r = 0; r < N; r++ ) {
        fscanf(g,"%lg",&hd[r]);
    }
    for (r = 0; r < J; r++ ) {
        for (s = 0; s < K; s++ ){
            HD[r][s] = 0;
        }
    }
    for (r = 0; r < N; r++ ) {
            HD[j[r]-1][k1[r]-1] = hd[r];
    }
    fclose(g);
    //printf("%g",*HS[0][3]);

    int x;
    x = genera_variables(I,J,K,M,N,HS,HD);




    return(0);
}

