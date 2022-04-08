#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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
void fonte(double **w, double temp1,double temp2)
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
void troca(double **w , double **aux)
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

void salva_dados(double **w,int l)
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

void discretizar(double **w, double **aux , double lambda)
{

	int i,j;
	
	#pragma omp parallel shared(w,aux,lambda) private(i,j)
	{
		#pragma omp for schedule(dynamic)
		for(i = 1 ; i < X-1 ; i++)
		{		
			for(j = 1; j < Y-1 ; j++)
			{
				
				// Aqui verifica se o o conjunto (i,j) não se encontra sobre uma fonte, caso não seja uma fonte a conta é feita	
				if((i-X*0.5)*(i-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) != 22500 && (i-X*0.5)*(i-X*0.5)+ (j-Y*0.5)*(j-Y*0.5) != 202500)
					w[i][j] = lambda*(aux[i+1][j]+aux[i-1][j]+aux[i][j+1]+aux[i][j-1]-4*aux[i][j])+aux[i][j];
								
			}
		}
	}

	return;
}


void main()
{
	int i,j,l,temp1,temp2;
	double dx,dy,dt,lambda,**w, **aux;
	
	
	//ALOCANDO AS MATRIZES W e AUX
	w =(double**)malloc(X*sizeof(double));
	for (i =0; i < X ; i++)
	{
		w[i]=(double*)malloc(Y*sizeof(double));
	}
	
	aux =(double**)malloc(X*sizeof(double));
	for (i =0; i < X ; i++)
	{
		aux[i]=(double*)malloc(Y*sizeof(double));
	}
	
	
	// Passo 1: Determinar dx, dy , dk e lambda. Preencher a matriz com zeros
	
	dx = (xb-xa)/(double) X;
	
	dy = (yd-yc)/(double) Y;
	
	dt = TEMPO/(double) T;
	
	 
	
	
	lambda = (alpha*dt)/pow(dx,2); // Lambda deve ser menor ou igual a 1/2  para haver convergencia
	
	// A minha matriz W é a matriz solucao para o tempo K+1, a matriz AUX carregara os valores
	// do tempo K e a matriz W do tempo k+1. Assim consigo fazer uma evolucao no tempo
	
	// Preenchendo a matriz w com zeros
	
	for(i=0 ; i < X ; i++)
	{
		for(j=0 ; j < Y ; j++)
		{
			w[i][j] = 0.0;
			
		}
	}
	
	//Passo 2: Impor a condicoes de contorno
	
	//Parte de baixo da barra
	for(j=0; j < Y ; j++)
	{
		w[0][j] =  293.0;
	}
	
	//Parte de cima da barra	
	for(j=0; j < Y ; j++)
	{
		w[X-1][j] = 293.0;
	}
	
	//Lateral esquerda da placa
	for(i=0; i < X ; i++)
	{
		w[i][0] = 293.0;
	}
	
	//Lateral direita da placa
	for(i=0; i < X ; i++)
	{
		w[i][Y-1] = 293.0;
	}
	
	
	
	
	// Passo 3: Impor condicoes iniciais nos pontos internos da malha
	// 1<= i < X e 1<= j < Y
	
	for(i=1 ; i< X-1; i++)
	{
		for(j=1; j < Y-1 ;j++)
			w[i][j]=inicial(i,j);
	}
	
	
	
	
	//Passo 4: Aplicar a discretizacao e resolver para o tempo  t+1
	
	//t=0;
	l=0;
	
	temp1 = 1000;
	temp2 = 500;
	fonte(w,temp1,temp2);
	for(l=0;l<=500000;l++)
	{
		//Gravando a solucao do tempo t de 10 mil em 10 mil iterações
		if(l%10000 == 0)
		{	
			salva_dados(w,l);
		}
		// Inserindo os dados da matriz w na matriz aux
		troca(w,aux);
		
		//Realizando a discretizacao no interior da malha
		// 1<= i < X e 1<= j < Y
		
		discretizar(w,aux,lambda);
		//printf("\nAquecendo a placa \nTempo total: %f Tempo= %f l: %d \n",TEMPO/2.0,t,l);
		//t=t+dt;		
	
	}
	
	temp1 = 400;
	temp2 = 400;
	fonte(w,temp1,temp2);
	for(l=500001;l<=1000000;l++)
	{
		//Gravando a solucao do tempo t de 10 mil em 10 mil iterações
		if(l%10000 == 0)
		{	
			salva_dados(w,l);
		}
		// Inserindo os dados da matriz w na matriz aux
		troca(w,aux);
		
		//Realizando a discretizacao no interior da malha
		// 1<= i < X e 1<= j < Y
		
		discretizar(w,aux,lambda);
		//printf("\nEsfriando a placa \nTempo total: %d Tempo= %f l: %d \n",TEMPO,t,l);
		//t=t+dt;		
	
	}
	

	
}
