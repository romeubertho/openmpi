// 5248962 Edylson Torikai
// 4471070 Fernanda Tostes
// 8531289 Guilherme Amorim Menegali
// 9266664 Tainá Andrello Piai
// 9293007 Fernando Vinícius Gianini Silva
// 7151905 Romeu Bertho Junior

#include<stdio.h>
#include<stdlib.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>

int mensagem[5];
int npes;
int myrank;
MPI_Status status;

void calc_diferenca(int **matriz, int linha_maior, int coluna_maior, int *linha_menor, int
*coluna_menor, int *valor_maior, int *valor_menor) {
    int dif = INT_MIN, aux_dif = INT_MIN;
    int i, j;
    int mensagem[8];
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    for (i = linha_maior - 1; i < linha_maior + 2; i++) {
        for (j = coluna_maior - 1; j < coluna_maior + 2; j++) {

            if (!(i < 0 || i > linha_maior || j < 0 || j > coluna_maior)) {
                aux_dif = abs(matriz[linha_maior][coluna_maior] - matriz[i][j]);
                if (aux_dif > dif) {
                    dif = aux_dif;
                    *linha_menor = i;
                    *coluna_menor = j;
                    *valor_maior = matriz[linha_maior][coluna_maior];
                    *valor_menor = matriz[i][j];
                }
            }
        }
    }

//	printf("#%d > maior diferenca do elemento %d = %d\n",myrank, matriz[linha_maior][coluna_maior], dif);
    mensagem[0] = dif;
    mensagem[1] = *linha_menor;
    mensagem[2] = *coluna_menor;
    mensagem[3] = *valor_maior;
    mensagem[4] = *valor_menor;
    mensagem[5] = myrank;

//	printf("#%d > Sending...\n",myrank);
    MPI_Send(mensagem, 6, MPI_INT, 0, 1, MPI_COMM_WORLD);
//	MPI_Recv(mensagem, 6, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//	printf("#%d > Received...",mensagem[5]);
}

int main(int argc, char *argv[]) {

    FILE *arquivo_entrada;
    arquivo_entrada = fopen("numeros.txt", "r");
    int linha_menor_aux, coluna_menor_aux, linha_maior_aux, coluna_maior_aux, valor_maior_aux, valor_menor_aux, tam, i = 0, j = 0;
    int linha_menor, coluna_menor, linha_maior, coluna_maior, valor_maior, valor_menor = 0;

    fscanf(arquivo_entrada, "%d\n", &tam);

    int **matriz = ((int **) malloc(tam * sizeof(int *)));
    for (i = 0; i < tam; i++) {
        matriz[i] = ((int *) malloc(tam * sizeof(int)));
    }
    for (i = 0; i < tam; i++) {
        for (j = 0; j < tam; j++) {
            fscanf(arquivo_entrada, "%d\n", &(matriz[i][j]));
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int aux_comp = INT_MIN;
    int diferenca = INT_MIN;
    int pos_i, pos_j;
    int count = 0;

    for (i = 0; i < tam; i++) {
        for (j = 0; j < tam; j++) {
            if (myrank != 0) {
                if (myrank == (((i * tam + j) % (npes - 1)) + 1)) {
//				printf("#%d > (%d,%d)\n",myrank,i,j);
                    calc_diferenca(matriz, i, j, &linha_menor_aux, &coluna_menor_aux, &valor_maior_aux,
                                   &valor_menor_aux);
                }
            } else if (myrank == 0) {
                printf("#%d > Waiting... %d\n", myrank, count);
                count++;
                MPI_Recv(mensagem, 6, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                printf("#%d > Mensagem recebida: mensagem[0]=%d from rank = %d\n", myrank, mensagem[0], mensagem[5]);
//			MPI_Send(mensagem, 6, MPI_INT, mensagem[5], 1, MPI_COMM_WORLD);

                if (diferenca < mensagem[0]) {
                    diferenca = mensagem[0];
                    linha_maior = i;
                    coluna_maior = j;
                    linha_menor = mensagem[1];
                    coluna_menor = mensagem[2];
                    valor_maior = mensagem[3];
                    valor_menor = mensagem[4];
                    printf("#%d > Maior: (%d,%d) = %d | Menor: (%d,%d) = %d\n", myrank, linha_maior, coluna_maior,
                           mensagem[3], linha_menor, coluna_menor, mensagem[4]);
                }
            }
        }
    }

    if (myrank == 0)
        printf("M[%d,%d]=%d M[%d,%d]=%d\n", linha_maior, coluna_maior, valor_maior, linha_menor, coluna_menor,
               valor_menor);

    free(matriz);
    fclose(arquivo_entrada);

    MPI_Finalize();
}