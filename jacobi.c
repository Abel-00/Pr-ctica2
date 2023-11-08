#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define M 5


int main() {
    // Declarar y definir la matriz M y el vector vect
    float Ma[M][M] = {
        {10.0, 2.0, -3.0, 1.0, 1.0},
        {2.0, 13.0, -8.0, 1.0, 1.0},
        {3.0, 8.0, 13.0, 1.0, 1.0},
	{1.0, 2.0, 7.0, 20.0, 2.0},
        {2.0, 1.0, 2.0 , 3.0, 12.0}	
    };

    float vect[M] = {1.0, 4.0, 7.0, 12.0,10.0};
    long double vectres[M] = {0, 0, 0,0,0};
    long double temp[M] = {0,0,0,0,0};
    long double vecterror[M] = {0,0,0,0,0};
    // Mostrar la matriz M
   // printf("Matriz M:\n");

//for (int i = 0; i < M; i++) {
     // for


for (int i = 0; i < M; i++) {
        float divisor = 0; 
	divisor = Ma[i][i];  // Aseguramos el valor del divisor

        for (int c = 0; c < M; c++) {
            Ma[i][c] = Ma[i][c] /( (-1)*divisor);  // Dividir cada elemento de la fila por el valor de la diagonal
        }
        vect[i] = vect[i] / divisor; // Dividir cada elemento del vector por el valor de la diagonal
    }

for (int i= 0; i<M; i++){
	Ma[i][i] = 0;}


for (int i = 0; i < 100000; i++) {
        float temp[M] = {0.0};  // Vector temporal inicializado con ceros

        for (int f = 0; f < M; f++) {
            for (int c = 0; c < M; c++) {
                temp[f] += (Ma[f][c] * vectres[c]); 
            }
	    temp[f] +=vect[f];
        }

        // Actualizar el vector vectres con los valores temporales
        for (int f = 0; f < M; f++) {
	    vecterror [f] = fabs(temp[f] - vectres[f]); 
            vectres[f] = temp[f];
        }
    }

// Mostrar el resultado en vectres
    printf("Resultado en vectres:\n");
    for (int f = 0; f < M; f++) {
        printf("Resultat: %Lf\n", vectres[f]);
	printf("Error: %Lf\n",vecterror [f]);
    }
}
