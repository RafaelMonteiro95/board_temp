import numpy as np

N = int(input())

temps = {'top'	 : 20,
		 'bottom': 30,
		 'left'	 : 25,
		 'right' : 20}

def generateSystem(N,temps):

	#n eh a "ordem" da matriz com os pontos (nao da matriz real)
	n = N-1

	#a matriz de coeficientes é de tamanho (n-1)**2
	mat_size = n**2

	#criando uma matriz com zeros de dimensões [N*N (altura), N*N(Largura)] para a matriz A do sistema
	matrix = np.zeros([mat_size, mat_size])

	#criando um array com zeros de tamanho N*N para o vetor B do sistema
	b_vector = np.zeros(mat_size)

	#preenchendo a matriz com os valores de cada linha
	for row in range(n):
		for col in range(n):

			#calculando o índice do ponto atual da 
			current = (row*n)+col

			#calculando indices dos pontos vizinhos
			down 	= ((row+1)*n) + col #ponto abaixo
			up 		= ((row-1)*n) + col #ponto acima
			left 	= (row*n) + (col-1) #ponto a direita
			right 	= (row*n) + (col+1) #ponto a esquerda

			#colocando valor do ponto atual
			matrix[current][current] = 4

			#caso eu esteja olhando para um ponto na borda de cima
			if((row-1) < 0):
				b_vector[current] += temps['top']
			else:
				matrix[current][up] = -1

			#caso eu esteja olhando um ponto na borda da esquerda
			if(col-1 < 0):
				b_vector[current] += temps['left']
			else:
				matrix[current][left] = -1

			#caso eu esteja olhando um ponto na borda de baixo
			if((row+1) >= n):
				b_vector[current] += temps['bottom']
			else:
				matrix[current][down] = -1

			#caso eu esteja olhando um ponto na borda da direita
			if(col+1 >= n):
				b_vector[current] += temps['right']
			else:
				matrix[current][right] = -1

	return {'A':matrix,'B':b_vector}