import numpy as np

q_arquivos = 100 # quantidade de arquivos que tem

i=0

desvio = 0
erro = []
while(i < q_arquivos):
	arq_serial = 'serial/'+str(i) +'.dat'
	arq_paralelo = 'mpi/'+str(i)+'.dat'
	

	
	serial = np.loadtxt(arq_serial)
	paralelo = np.loadtxt(arq_paralelo)
	
	k=j=0

	for k in range(1440):
		for j in range(1440):
			desvio_quadratico = desvio + (paralelo[k,j] - serial[k,j])**2 
			
	print('\nSerial: ',arq_serial,'\nParalelo: ',arq_paralelo)
	erro.append(desvio_quadratico/2073600)
	i=i+1
	
arquivo_saida = open('Erros.dat','w')

j = 0
for i in erro:
	arquivo_saida.write(str(j)+" "+str(i)+"\n")
	j=j+1

arquivo_saida.close()	 
