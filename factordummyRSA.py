#!/usr/bin/env python3
#Example implementation of Dixon's factorization, Pollard's p-1 and Wiener's attack
import math
import random
from operator import add

def dixon(modulus, smoothness):
	#Generating set of smooth primes
	notFound = True
	trivial = True
	if smoothness > 35:
		return "Maximum smoothness is 35."
	primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149]
	B = primes[:smoothness]
	#For our set of primes we simply select the first n primes, this is the modern approach
	lower_bound = int(math.sqrt(modulus))
	#We are searching for smooth squares mod the modulus
	#We are going to repeat this process until we have a correct solution
	while trivial:
		#We are going to repeat this process until we have a square on both sides
		while notFound:
			candidates = []
			counter = 0
			#We are going to look for 1 more candidates than the cardinality of our prime set
			while (counter < len(B)+1):
				candidate = random.randint(lower_bound, modulus)
				terms = [0 for x in range(smoothness)]
				P = candidate ** 2 % modulus
				P_temp = P
				for i in range(len(B)):
					while (P_temp % B[i] == 0):
						P_temp = P_temp/B[i]
						terms[i] = terms[i]+1
				if (P_temp == 1):
					candidates.append([candidate, terms])
					counter = counter + 1
			#At this point we should have a large enough pool of B-smooth squares --> use linear algebra to find nontrivial factors
			L = 1
			R = [0 for x in range(smoothness)]
			for i in candidates:
				L = L * i[0]
				R = list(map(add,R,i[1]))
			if (all((x % 2 == 0) for x in R)):
				#Complete squares on both sides
				notFound = False
		Y = 1
		for i in range(smoothness):
			Y = int( Y * (B[i] ** (R[i]/2)))
		if (math.gcd(Y-L,modulus) != 1 and math.gcd(Y-L,modulus) != modulus):
			trivial = False
			solution = [math.gcd(Y-L,modulus), int(modulus/math.gcd(Y-L,modulus))]
		else:
			notFound = True
	print(solution)

dixon(84923,7)