#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

// Dimensoes da placa e tempo total

#define xa 0

#define xb 1

#define yc 0

#define yd 1

#define TEMPO 1

// Particoes do espaco e tempo

#define X 1440//1500

#define Y 1440//1500

#define T 1E6

// Constante alpha, ela que representara o material da placa

#define alpha 1E-1



// Fonte
void fonte(double w[X][Y], double temp1 ,double temp2)
{

	int i,j;

	
	for(i = 1 ; i < X-1 ; i++)
	{
		for(j = 1; j < Y-1 ; j++)
		{
				
			if((i-X*0.5)*(i-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) == 22500 )
				w[i][j] = temp1; 
						
			if( (i-X*0.5)*(i-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) == 202500)
				w[i][j] = temp2; 					
				
		}
	}
	return;
}
//Condicao inicial
double inicial(double x, double y)
{
	double resultado;
		
	resultado = 0 ;
	
	return resultado;
}

//Funcao auxiliar
void troca(double w[][Y] , double aux[][Y])
{
	int i,j;
	
	#pragma omp parallel shared(w,aux) private(i,j)
	{
		#pragma omp for schedule(dynamic)
		for(i=0; i< X ; i++)
		{
			for(j=0; j < Y ; j++)
			{
				aux[i][j] = w[i][j];
			}
		}
	}

	return;
}

// Função para salvar os dados em arquivos de texto

void salva_dados(double w[][Y],int l)
{
	
	int i,j,v;
	char arquivo[10] = ".dat" , nome[10];	

	FILE *fp;
		
	v = l*0.0001;
	sprintf(nome,"%d",v);
	strcat(nome,arquivo);
	fp=fopen(nome,"w");
			
	for(i=0 ; i < X ; i++)
	{
		for(j=0 ; j < Y ; j++)
		{
			fprintf(fp,"%lf ",w[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	return;


}

void discretizar(double w[][Y], double aux[][Y] , double lambda, int linha, int coluna,int rank,int size)
{

	int i,j,map;
	
	map = rank*linha;
	
	//printf("Estou discretizando RANK: %d\n",rank);
	
        if(size != 1)
        {
		if(rank == 0)
		{
		        #pragma omp parallel shared(w,aux,lambda) private(i,j)
		        {
		        	#pragma omp for schedule(dynamic)
				for(i = 1; i < linha ; i++)
				{
				        for(j = 1; j < coluna-1 ; j++)
				        {
				    		if((i-X*0.5)*(i-X*0.5)+(j-Y*0.5)*(j-Y*0.5) != 22500 && (i-X*0.5)*(i-X*0.5)+(j-Y*0.5)*(j-Y*0.5) != 202500)
							w[i][j] = lambda*(aux[i+1][j]+aux[i-1][j]+aux[i][j+1]+aux[i][j-1]-4*aux[i][j])+aux[i][j];
				        }

				}
			}
		}

		else if(rank == size-1)
		{
			#pragma omp parallel shared(w,aux,lambda) private(i,j)
			{
				#pragma omp for schedule(dynamic)
				for(i = 0; i < linha-1 ; i++)
				{
				        for(j = 1; j < coluna-1 ; j++)
				         {
				    		if((i+map-X*0.5)*(i+map-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) != 22500 && (i+map-X*0.5)*(i+map-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) != 202500)
							w[i][j] = lambda*(aux[i+map+1][j]+aux[i+map-1][j]+aux[i+map][j+1]+aux[i+map][j-1]-4*aux[i+map][j])+aux[i+map][j];
				        }
				}
			}
		}
		       else
			{
		        	#pragma omp parallel shared(w,aux,lambda) private(i,j)
		        	{
		        		#pragma omp for schedule(dynamic)
					for(i = 0; i < linha ; i++)
					{
						for(j = 1; j < coluna-1 ; j++)
						 {
					    		if((i+map-X*0.5)*(i+map-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) != 22500 && (i+map-X*0.5)*(i+map-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) != 202500)
								w[i][j] = lambda*(aux[i+map+1][j]+aux[i+map-1][j]+aux[i+map][j+1]+aux[i+map][j-1]-4*aux[i+map][j])+aux[i+map][j];
						}
					}
				}

			}
        }
        else
        {
               #pragma omp parallel shared(w,aux,lambda) private(i,j)
               {
               	#pragma omp for schedule(dynamic)
		        for(i = 1; i < linha-1 ; i++)
		        {
		                for(j = 1; j < coluna-1 ; j++)
		                 {
				    		if((i-X*0.5)*(i-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) != 22500 && (i-X*0.5)*(i-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) != 202500)
							w[i][j] = lambda*(aux[i+1][j]+aux[i-1][j]+aux[i][j+1]+aux[i][j-1]-4*aux[i][j])+aux[i][j];
				  }
		        }
		}
        }
        
	return;
}


void main(int argc , char **argv)
{
	MPI_Init(&argc,&argv);
	
	int i, j, l, temp1, temp2, rank, size, q_elementos, q_elementos_aux;
	double dx, dy, dt, lambda;
	
	int tiras_row, tiras_col;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	
	
	//Aqui vou estabelecer o tamanho das tiras horizontais que seráo enviadas para cada processo
	// Todas as tiras conterão 1442 colunas, as linhas que irão mudar
	
	tiras_row =(int) X / size; // dividindo o total de linhas pelo total de processos
	
	tiras_col = Y;
	
	q_elementos = tiras_row*tiras_col;
	q_elementos_aux = X*Y;
	
	// cada tira terá um total de X/size linhas e Y = 1442 colunas
	
	/* A quantidade de linhas, aqui são 1442, deverá ser um multiplo da quantidade de processos.
	
	daria para fazer sem ser multiplo, entretanto, o último processo ficaria com mais trabalho que os demais
	
	uma vez que sua tira deverá conter mais linhas para que se complete 1442 */
	
	double placa[X][Y], placa_aux[X][Y], placa_recv[tiras_row][tiras_col];
	
	
	// Passo 1: Determinar dx, dy , dk e lambda. Preencher a matriz com zeros
	
	dx = (xb-xa)/(double) X;
	
	dy = (yd-yc)/(double) Y;
	
	dt = TEMPO/(double) T;
	
	 
	lambda = (alpha*dt)/pow(dx,2); // Lambda deve ser menor ou igual a 1/2  para haver convergencia
	
	// A minha matriz W é a matriz solucao para o tempo K+1, a matriz AUX carregara os valores
	// do tempo K e a matriz W do tempo k+1. Assim consigo fazer uma evolucao no tempo
	
	// Preenchendo a matriz w com zeros
	
	// O processo 0 criará a placa impondo suas condições de contorno e condições iniciais
	
	l=0;
	temp1 = 1000;
	temp2 = 500;
	
	if (rank == 0)
	{
		
		for(i=0 ; i < X ; i++)
		{
			for(j=0 ; j < Y ; j++)
			{
				placa[i][j] = 0.0;
					
			}
		}
		
		//Passo 2: Impor a condicoes de contorno
		
		//Parte de baixo da barra
		for(j=0; j < Y ; j++)
		{	
			placa[0][j] =  293.0;
		}
		
		//Parte de cima da barra	
		for(j=0; j < Y ; j++)
		{
			placa[X-1][j] = 293.0;
		}
		
		//Lateral esquerda da placa
		for(i=0; i < X ; i++)
		{
			placa[i][0] = 293.0;
		}
		
		//Lateral direita da placa
		for(i=0; i < X ; i++)
		{
			placa[i][Y-1] = 293.0;
		}
		
		
		
		
		// Passo 3: Impor condicoes iniciais nos pontos internos da malha
		// 1<= i < X e 1<= j < Y
		
		for(i=1 ; i< X-1; i++)
		{
			for(j=1; j < Y-1 ;j++)
			{
				placa[i][j] = inicial(i,j);
			}
		}
		
		
		// Gerando as fontes na placa
		fonte(placa,temp1,temp2);
		
		// preenchendo a aux
		//troca(placa,placa_aux);
		
		//salvando no tempo 0
		//salva_dados(placa,0);
		
		
	}
	
	
	//Passo 4: Aplicar a discretizacao e resolver para o tempo  t+1
	
	//aquece
	
	for(l=0;l<=500000;l++)
	{
		// o processo zero atualizará a matriz aux e salvara os dados caso esteja na hora
		if(rank == 0)
		{
			//printf("Iteracao: %d\n",l);
			
			troca(placa,placa_aux);
			if(l%10000 == 0)
				salva_dados(placa,l);	
		}		
	
		// reparte para todos os processos a placa 
		MPI_Scatter(placa,q_elementos,MPI_DOUBLE, placa_recv,q_elementos,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		// Envia para todos os processos a placa auxiliar do tempo anterior
		MPI_Bcast(placa_aux, q_elementos_aux,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		// cada processo fará a conta na sua respectiva região
		discretizar(placa_recv,placa_aux,lambda,tiras_row,tiras_col,rank,size);
		
		// cada processo enviara a sua região calculada para o processo root = 0
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(placa_recv,q_elementos,MPI_DOUBLE, placa,q_elementos,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		
	}
	
	// resfria
	if(rank == 0)
	{
		temp1 = temp2 = 400;
		fonte(placa,temp1,temp2);	
	}
	for(l=500001;l<=1000000;l++)
	{
		// o processo zero atualizará a matriz aux e salvara os dados caso esteja na hora
		if(rank == 0)
		{
			//printf("Iteracao: %d\n",l);
			
			troca(placa,placa_aux);
			if(l%10000 == 0)
				salva_dados(placa,l);	
		}		
	
		// reparte para todos os processos a placa 
		MPI_Scatter(placa,q_elementos,MPI_DOUBLE, placa_recv,q_elementos,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		// reparte para todos os processos a placa auxiliar to tepo anterior
		MPI_Bcast(placa_aux, q_elementos_aux,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		// cada processo fará a conta na sua respectiva região
		discretizar(placa_recv,placa_aux,lambda,tiras_row,tiras_col,rank,size);
		
		// cada processo enviara a sua região calculada para o processo root = 0
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(placa_recv,q_elementos,MPI_DOUBLE, placa,q_elementos,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	}
	
	MPI_Finalize();

}
