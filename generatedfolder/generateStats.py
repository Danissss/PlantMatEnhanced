import sys
import os
import math

from rdkit import Chem


def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


def calculateExpectedNumber(num_oh, num_sugar):

	total_num = 0;
	total_num += (num_sugar ** num_oh)
	for i in range(1,num_oh):
		tmp1 = num_sugar ** (num_oh - i)
		tmp2 = nCr(num_oh, num_oh - i)
		total_num += (tmp1 * tmp2)

	return total_num

def main():

	# num_oh = int(sys.argv[1])
	# num_sugar = int(sys.argv[2])
	# 																			# num_oh   , num_sugar
	# print("Expected number                : {0}".format(calculateExpectedNumber(num_oh,num_sugar)))


	suppl = Chem.SDMolSupplier(sys.argv[1])
	mol_array = []
	inds = 0
	for mol in suppl:

		if mol != None:
			mol_array.append(mol)
		else:
			print(inds)
			sys.exit(0)
		inds += 1

			
	# num_oh = int(sys.argv[2])
	# num_sugar = int(sys.argv[3])
	# 																			# num_oh   , num_sugar
	# print("Expected number                : {0}".format(calculateExpectedNumber(num_oh,num_sugar)))
	print("Number of generated metabolites: {0}".format(len(mol_array)))


if __name__ == '__main__':
	main()
