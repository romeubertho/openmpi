// para compilar: gcc medianacol.c -o medianacol -Wall
// para executar: medianacol

/*
Dada uma matriz M de tamanho LxC, composta por números inteiros, construa um algoritmo em C que insira a mediana
dos valores de cada coluna da matriz M em um vetor V de tamanho C e imprima esse vetor de medianas como saída
dessa aplicação.
As dimensoes de M e os números de M serão carregados de um arquivo texto no início da execução.
A primeira linha do arquivo de entrada deve indicar as dimensoes L e C, respectivamente.
As L demais linhas do arquivo indicam os C elementos de cada linha de M.

Considere este exemplo com uma matriz M de 6x4.
Os dados lidos do arquivo são:
6	4
9	3	7	5
4	12	2	40
8	8	4	32
6	14	32	21
33	44	20	1
10	18	17	10

As colunas ordenadas da matriz são (dados temporários):
4	3	2	1
6	8	4	5
8	12	7	10
9	14	17	21
10	18	20	32
33	44	32	40

Saída da aplicação salva em saida.txt (medianas das colunas):
8.5	13.0	12.0	15.5

*/


// Based on Author: Paulo S L Souza
// quick sort code adapted from https://www.geeksforgeeks.org/quick-sort/


#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>

void ordena_colunas(int *, int, int);
void calcula_mediana(int *, float *, int, int);
void quicksort(int *, int, int, int);
int partition (int *, int, int, int);



/* low  --> Starting index,  high  --> Ending index */
//https://www.geeksforgeeks.org/quick-sort/
void quicksort(int *arr, int low, int high, int C)
{
    int pi;

    if (low < high)  {
        /* pi is partitioning index, arr[pi] is now
           at right place */
        pi = partition(arr, low, high, C);

        quicksort(arr, low, pi - 1, C);  // Before pi
        quicksort(arr, pi + 1, high, C); // After pi
    }

} // fim quicksort

/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot
   https://www.geeksforgeeks.org/quick-sort/
*/
int partition (int *arr, int low, int high, int C)
{
    int i, j, swap, pivot;

    // pivot (Element to be placed at right position)
    pivot = arr[high*C];

    i = (low - 1);  // Index of smaller element

    for (j = low; j <= high-1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[j*C] <= pivot)
        {
            i++;    // increment index of smaller element

            // swap arr[i] and arr[j]
            swap = arr[i*C];
            arr[i*C] = arr[j*C];
            arr[j*C] = swap;
        }
    }

    //swap arr[i + 1] and arr[high]
    swap = arr[(i + 1)*C];
    arr[(i + 1)*C] = arr[high*C];
    arr[high*C] = swap;

    return (i + 1);

} // fim partition



void ordena_colunas(int *M, int L, int C)
{
    int j;

    for (j = 0; j < C; j++) {
        //manda o endereco do primeiro elemento da coluna, limites inf e sup e a largura da matriz
        quicksort(&M[j], 0, L-1, C);
    }
} // fim ordena_colunas


void calcula_mediana(int *M, float *vet, int L, int C)
{
    int j;

//  printf("calcula_mediana: L=%d e C=%d \n", L, C);

    for (j = 0; j < C; j++) {
        vet[j] = M[((L/2)*C)+j];
//    printf("vet[%d] = %f, pos = %d, M[]=%d, ", j, vet[j], (L/2)*C+j, M[((L/2)*C)+j]);
        // se o nr de elementos eh par, tem que fazer a media
        // das duas medianas. Tem que pegar a mediana anteior.
        if(!(L%2))  {
            vet[j]+=M[((((L-1)/2)*C)+j)];
//    printf("vet[%d] = %f, pos = %d, M[]=%d, \n", j, vet[j], (((L-1)/2)*C)+j, M[((((L-1)/2)*C)+j)]);
            vet[j]*=0.5;
        } // fim do if
    } // fim do for

    return;
} // fim calcula_mediana



int main(int argc, char *argv[])
{
    int L, C;
    int *M;
    float *vet;

    int i, j;

    FILE *arquivo_entrada,*arquivo_saida;

    if(!(arquivo_entrada=fopen("entrada.txt","r")))
    {
        printf("Erro ao abrir entrada.txt como leitura! Saindo! \n");
        return(-1);
    }

    if(!(arquivo_saida=fopen("saida.txt","w+")))
    {
        printf("Erro ao abrir/criar saida.txt como escrita. Saindo! \n");
        return(-1);
    }

    // Leitura da primeira linha de entrada.txt contendo as dimensoes de M
    fscanf(arquivo_entrada, "%d %d", &L, &C);
    //printf("Ordem da Matriz M: L=%d C=%d\n", L, C);

    // criando M[LxC]
    M = (int *) malloc ( L * C * sizeof(int));

    // criando vet[C] do tipo float. Este vetor terá as medianas das colunas.
    vet = (float *) malloc (C * sizeof (float) );

    // carregando M do arquivo
    for(i=0; i<L; i++)
        for(j=0; j<C; j++)
            fscanf(arquivo_entrada, "%d", &M[i*C+j]);
/*
    // impressao para verificacao apenas
    for(i=0; i<L; i++)
    {
      for(j=0; j<C; j++)
	  printf("%d	", M[i*C+j]);
      printf("\n");
    }
*/
    // ordena as colunas de M usando o quicksort
//    ordena_colunas(M, L, C);
    int npes;
    int myrank;
    float mediana[2];
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    for (j = 0; j < C; j++) {
        //manda o endereco do primeiro elemento da coluna, limites inf e sup e a largura da matriz
        if (myrank != 0) {
            if (myrank == (((j) % (npes - 1)) + 1)) {
                printf("#%d > col=%d\n", myrank, j);
                quicksort(&M[j], 0, L - 1, C);
                calcula_mediana(M, vet, L, C);
                printf("#%d > mediana da col %d = %.1f,\n", myrank, j, vet[j]);
                mediana[0]=vet[j];
                mediana[1]=myrank;
                MPI_Send(mediana, 2, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
            }
        }
            else if(myrank==0)
            {
                MPI_Recv(mediana, 2, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                vet[j]=mediana[0];
                printf("#%d > mediana da col %d = %.1f from #%.0f,\n",myrank,j,vet[j],mediana[1]);
            }
    }



    //printf("Matriz com colunas ordenadas:\n");
    // impressao para verificacao apenas
//    for(i=0; i<L; i++)
//    {
//      for(j=0; j<C; j++)
//	  printf("%d	", M[i*C+j]);
//      printf("\n");
//    }

//    calcula_mediana(M, vet, L, C);
    if(myrank==0) {
        for (j = 0; j < C; j++)
            fprintf(arquivo_saida, "%.1f, ", vet[j]);
        fprintf(arquivo_saida, "\n");
    }
    fclose(arquivo_entrada);
    fclose(arquivo_saida);

    free(vet);
    free(M);

    MPI_Finalize();

    return(0);

}