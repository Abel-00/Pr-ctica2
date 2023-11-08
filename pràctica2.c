#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define N 512


float V1[N], V2[N], V3[N], V4[N]; // Declarem els vectors
float Mat[N][N], MatDD[N][N]; // Declarem les matrius

void InitData(){
	int i,j;
	srand(4422543);
	//srand(8824553);
	for( i = 0; i < N; i++ )
 		for( j = 0; j < N; j++ ){
 			Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
 			if ( (abs(i - j) <= 3) && (i != j))
 				MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
 			else if ( i == j )
 				MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
 			else MatDD[i][j] = 0.0;
 }
for( i = 0; i < N; i++ ){
 	V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 	V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 	V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
}
}

void PrintVect(float vect[N], int from, int numel) {
    printf("Vector que comença a %d i acaba a %d:\n",from,numel);
    int lon = from + numel; // Calculem la quants nombres del vector volem mostrar
   
    for (int i = from; i < lon; i++) {
        printf("%f ", vect[i]);

    }
    printf("\n");
}
void PrintRow(float matriu[N][N], int row, int from, int numel){
	
	printf("Linea %d de la matriu que comença a %d i acaba a %d: \n",row,from,numel);
	int lon = from + numel;
	for (int i = from; i < lon; i++) {
		printf("%f ",matriu[row][i]);
	}
        printf("\n");	
			}
void MultEscalar(float vect[N], float vectres[N], float alfa){

	for (int i = 0; i < N; i++) {
                vectres[i] = alfa * vect[i];  
	}
}

float Scalar(float vect1[N], float vect2[N]){

	float res = 0;
        for (int i = 0; i < N; i++) {
                res = res + (vect1[i] * vect2[i]);
        }
	printf("El resultat és %f",res);
	return res;
        printf("\n");
}
float Magnitude(float vect1[N]) {
    float mag = 0;

    for (int i = 0; i < N; i++) {
        mag = mag + (vect1[i] * vect1[i]);
    }

    mag = sqrt(mag); // Calculate the square root
    printf("La magnitud del vector és %f",mag);
    return mag; // Return the magnitude
    printf("\n");
}
int Ortogonal( float vect1[N], float vect2[N] ){
	if (Scalar(vect1,vect2)==0){
		return 1;}
	else{
		return 0;}}
void Projection( float vect1[N], float vect2[N], float vectres[N] ){
	float alfa = (Scalar(vect1,vect2))/(Magnitude(vect2));
	MultEscalar(vect2,vectres,alfa);
}

float Infininorm( float M[N][N]){
	float maxim = 0;
	float resultat = 0;
	for (int f=0; f<N; f++){
		for (int c=0; c<N; c++){
			resultat = resultat + fabs(M[f][c]);	
			}
 		if (resultat > maxim){ maxim = resultat;}
                resultat = 0;

	}

	return maxim;

}
float Onenorm( float M[N][N] ){
        float maxim = 0;
        float resultat = 0;
        for (int c=0; c<N; c++){
                for (int f=0; f<N; f++){
                        resultat = resultat + fabs(M[f][c]);}
                if (resultat > maxim){ maxim = resultat;}
                resultat = 0;
                
        }
	return maxim;

}
float NormFrobenius( float M[N][N] ){
	float suma = 0;
	for (int c=0; c<N; c++){
                for (int f=0; f<N; f++){
		suma += (M[f][c]* M[f][c]);
		}
	}
	return sqrt(suma);	
}

int DiagonalDom( float M[N][N] ){
	float suma = 0;
    float diagonal = 0;
    for (int f=0; f<N; f++){
             suma =0;
             diagonal = M[f][f];
             for (int c=0; c<N; c++){
                  suma = suma + fabs(M[f][c]);}
     suma = suma - fabs(diagonal);
     if (suma > fabs(diagonal)) {
            return 0;
        }
    }
    return 1;

}

int Jacobi(float M[N][N] , float vect[N], float vectres[N], float vecterror[N], unsigned iter){


if (DiagonalDom(M) == 0){return 0;}

else{
    float temp[N];
    for (int i =0; i < N; i ++){
	vectres[i] = 0.0;
    	temp[i] = 0.0;
    	vecterror[i] = 0.0;
    }

    for (int i = 0; i < N; i++) {
        float divisor = 0;
        divisor = M[i][i];  // Aseguramos el valor del divisor

        for (int c = 0; c < N; c++) {
            M[i][c] = M[i][c] /( (-1)*divisor);  // Dividir cada elemento de la fila por el valor de la diagonal
        }
        vect[i] = vect[i] / divisor; // Dividir cada elemento del vector por el valor de la diagonal
      }

    for (int i= 0; i<N; i++){ M[i][i] = 0; } //Valors de la diagonal són 0

    for (int i = 0; i < iter; i++) {
	    for(int a =0; a<N; a++){ temp[a] = 0.0; }// Vector temporal inicializado con ceros

        for (int f = 0; f < N; f++) {
            for (int c = 0; c < N; c++) {
                temp[f] += (M[f][c] * vectres[c]);
            }
            temp[f] +=vect[f];
        }

        // Actualizar el vector vectres con los valores temporales
        	for (int f = 0; f < N; f++) {
            	vecterror [f] = fabs(temp[f] - vectres[f]);
            	vectres[f] = temp[f];
        }
		
    }

 return 1;

}

}
int main() {
    InitData();	
    //Apartat A (els resultats són els mateixos si utilizem la llavor srand(4422543)) 
    PrintVect(V1, 0, 10);
    PrintVect(V1, 256, 10);
    PrintVect(V2, 0, 10);
    PrintVect(V2, 256, 10);
    PrintVect(V3, 0, 10);
    PrintVect(V3, 256, 10);
    //Apartat B  (els resultats són els mateixos si utilizem la llavor srand(4422543))
    PrintRow(Mat,0,0,10);
    PrintRow(Mat,100,0,10);
    //Apartat C
    PrintRow(MatDD,0,0,10);
    PrintRow(MatDD,100,95,10); // Perque surti igual que la pràctica hem de visualitzar els elements del 95 al 104 i no del 90 al 99 
    //Apartat D 
    //Matriu Mat
    printf("Mat \n");
    printf("Infininorma: %f \n",Infininorm(Mat));
    printf("Norma u: %f \n",Onenorm(Mat));
    printf("Frobenius: %f \n",NormFrobenius(Mat));//No ens dona del tot exacte
    if(DiagonalDom(Mat)==0){
    	printf("La matriu no és diagonal dominant \n");}
    else{
	printf("La matriu és diagonal dominant \n");}
    //Matriu MatDD
    printf("MatDD \n");
    printf("Infininorma: %f \n",Infininorm(MatDD));
    printf("Norma u: %f \n",Onenorm(MatDD));
    printf("Frobenius: %f \n",NormFrobenius(MatDD));
    if(DiagonalDom(MatDD)==0){
    	printf("La matriu no és diagonal dominant \n");}
    else{
    	printf("La matriu és diagonal dominant \n");}
    //Apartat E
    Scalar(V1,V2);
    printf(" \n");
    Scalar(V1,V3);
    printf("\n");
    Scalar(V2,V3);
    printf("\n");
    //Apartat F
    Magnitude(V1);
    printf("\n");
    Magnitude(V2);
    printf("\n");
    Magnitude(V3);
    printf("\n");
    //Apartat G
    if (Ortogonal(V1,V2)==1){
    	printf("El vector 1 i el vector 2 són ortogonals");
    	printf("\n");}
    if (Ortogonal(V1,V3)==1){
    	printf("El vector 1 i el vector 3 són ortogonals");
    	printf("\n");}
    if (Ortogonal(V2,V3)==1){
    	printf("El vector 2 i el vector 3 són ortogonals");}    
    	printf("\n");
    //Apartat H
    MultEscalar(V3,V2,2);
    PrintVect(V2,0,10);
    PrintVect(V2,256,10);
    //Apartat I
    Projection(V2,V3,V4);
    PrintVect(V4,0,10);
    Projection(V1,V2,V4);
    PrintVect(V4,0,10);
    //Apartat J
    Jacobi(MatDD,V3,V4,V1,1000);
    printf("Resultat:\n"); 
    PrintVect(V4, 0, 10);
    printf("Error:\n");
    PrintVect(V1, 0, 10);
    
    return 0;
}
