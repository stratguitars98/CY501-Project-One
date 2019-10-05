def makeEqualLength(stringOne, stringTwo):
	lenOne = len(stringOne)
	lenTwo = len(stringTwo)

	if(lenOne<lenTwo):
		for i in range(lenTwo-LenOne):
			stringOne = '0' + stringOne
		return lenTwo
	
	elif(lenOne > lenTwo):
		for i in range(lenOne-lenTwo):
			stringTwo = '0' + stringTwo
	
	return lenOne


def addBitStrings(first, second):
	result = ""

	length = makeEqualLength(first, second)
	carry = 0

	for i in range(length-1, -1):
		firstBit = int(ord(first[i]) - ord('0'))
		secondBit = int(ord(second[i]) - ord('0'))

		summation = int((firstBit^secondBit^carry) + '0')

		result = char(summation) + result

		# this may cause an issue
		carry = (firstBit&secondBit) | (secondBit&carry) | (firstBit&carry)

	if carry == 1:
		result = 1+result

	return result


def multiplySingleBit(a, b):
	return(ord(a[0]) - ord('0'))*(ord(b[0]) - ord('0'))

def multiply(X, Y):
	n = makeEqualLength(X, Y)

	if(n == 0): return 0
	if(n == 1): return multiplySingleBit(X, Y)

	fh = int(n/2)+1
	sh = n - fh

	Xl = X[:fh]
	Xr = X[fh:sh]

	Yl = Y[:fh]
	Yr = Y[fh:sh]

	P1 = multiply(Xl, Yl)
	P2 = multiply(Xr, Yr)
	P3 = multiply(addBitStrings(Xl, Xr), addBitStrings(Yl, Yr))

	return P1*(1<<(2*sh)) + (P3 - P1 - P2)*(1<<sh) + P2


if __name__ == "__main__":
	num1 = input("Enter in a 64 bit number:\t")
	num2 = input("Enter in another 64 bit number:\t")

	print(multiply(num1, num2))