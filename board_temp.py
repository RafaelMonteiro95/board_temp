import numpy as np
import matplotlib.pyplot as plt

def generateSystem(N,temps):
	"""Generates a linear equation system for the board temperature problem
	Given a metal board with 4 different temperatures on its sides, splits this board in NxN points
	and calculate the temperature in each point using the mean temperature of its 4 neighboors
	the generated system is of format
	Ax = b
	
	Args:
		N (int): number of slices in each dimension
		temps (dict): a dictionary containing int values for the keys 'top', 'bottom', 'left', 'right'. Each value represents the temperature of a given edge of the metal board
	
	Returns:
		dict: return a dictionary where the keys 'A' and 'b' represents the A matrix and b vector of the system
	"""

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

	return {'A':matrix,'b':b_vector}

def gaussSeidel(A, b, x = None, iterations = 100):
	"""Solves a linear equation system using gauss-seidel method
	The system should be in format Ax = b
	
	Args:
		A (numpy float matrix): numpy matrix containing the A part of system 
		b (numpy float array): numpy array containing the b part of the system  
		x (numpy float array, optional): numpy array containing the initial x guesses. If none provided, uses a array full of zeroes
		iterations (int, optional): number of iterations. If none provided, uses 100 iterations
	
	Returns:
		numpy float array: an array with the calculated x's values
	"""

	#first, we split A in L and U
	# L is the lower triangle matrix of A
	L = np.tril(A)
	# U is the upper triangle matrix of A without the diagonal
	U = A - L

	#x0 is our initial guess array
	if x == None:
		x0 = np.zeros(len(A))
	else:
		x0 = np.array(x)

	#we need to calculate C and G for our iterations
	# L^-1
	inverse_L = np.linalg.inv(L)
	# -L^-1
	negative_inverse_L = -inverse_L
	# C = (-L^-1 dot U)
	C = np.dot(negative_inverse_L, U)
	# g = (L^-1 dot b)
	g = np.dot(inverse_L, b)

	#then, we iterate k iterations
	for k in range(iterations):
		#for each iteration, x^(k+1) = Cx + g
		x0 = np.dot(C, x0) + g

	#returns our answer
	return x0

def resultMatrix(temps, N, ans):
	""""turns temps and system answer into a matrix with temperatures"
	
	Args:
	    temps (dictionary): temperatures for each edge of the board
	    N (int): number of slices on the board
	    ans (array): answers of the linear system
	
	Returns:
	    TYPE: Description
	"""
	top = temps['top']
	bottom = temps['bottom']
	left = temps['left']
	right = temps['right']

	visu_matrix = np.zeros([N+1,N+1])

	#sides
	for i in range(1,N):
		visu_matrix[0][i] = top
		visu_matrix[i][0] = left
		visu_matrix[i][N] = right
		visu_matrix[N][i] = bottom
	#upper-left
	visu_matrix[0][0] = (top+left)/2.0
	#upper-right
	visu_matrix[0][N] = (top+right)/2.0
	#bottom-left
	visu_matrix[N][0] = (left+bottom)/2.0
	#bottom-right
	visu_matrix[N][N] = (right+bottom)/2.0

	row = 0
	col = 0
	for value in ans:
		visu_matrix[row+1][col+1] = ans[(row*(N-1) + col)]
		col += 1
		if(col == N-1):
			row += 1
			col = 0;

	return visu_matrix

if __name__ == "__main__":

	temps = {'top'	 : 20,
			 'bottom': 30,
			 'left'	 : 25,
			 'right' : 20}

	print('Enter number of slices:')
	#number of slices
	# N = int(input())
	N = 3

	#generating a linear equation system for the given temps and number of slices
	linsys = generateSystem(N,temps)

	#solving system using gaussSeidel method
	ans = gaussSeidel(linsys['A'],linsys['b'])
	print('x:', ans)

	visu = resultMatrix(temps, N, ans)
	print('resulting matrix:')
	print(visu)

	#plotting heatmap
	plt.imshow(visu, cmap='plasma', interpolation='gaussian',origin='lower')
	plt.show()