def makeEqualLength(stringOne, stringTwo):
	lenOne = len(stringOne)
	lenTwo = len(stringTwo)

	if(lenOne<lenTwo):
		for i in range(lenTwo-lenOne):
			stringOne = '0' + stringOne
		return (stringOne, stringTwo, lenTwo)
	
	elif(lenOne > lenTwo):
		for i in range(lenOne-lenTwo):
			stringTwo = '0' + stringTwo
	
	return (stringOne, stringTwo, lenOne)


def addBitStrings(first, second):
	result = ""

	first, second, length = makeEqualLength(first, second)
	carry = 0

	for i in range((length-1), -1, -1):
		firstBit = int(int(first[i]) - int('0'))
		secondBit = int(int(second[i]) - int('0'))
		#print(firstBit, secondBit)
		summation = int((firstBit^secondBit^carry) + int('0'))

		if(summation == 1):
			result = '1' + result
		else:
			result = '0'+result

		# this may cause an issue
		carry = (firstBit&secondBit) | (secondBit&carry) | (firstBit&carry)

	if carry == 1:
		result = '1'+result

	return result


def multiplySingleBit(a, b):
	a0 = int(a[0])
	b0 = int(b[0])
	r0 = int('0')
	temp = (a0-r0)*(b0-r0)

	#print('a[0]: {}\tb[0]: {}\t(a[0]-0)*(b[0]-0): {} '.format(a[0], b[0], temp))
	return (temp)

def multiply(X, Y):
	X, Y, n = makeEqualLength(X, Y)
	#print("N:",n)
	#print(X, Y, n)
	if(n == 0): return 0
	if(n == 1): return multiplySingleBit(X, Y)

	fh = int(n/2)
	#print(fh)
	sh = (n - fh)
	#print(sh)

	Xl = X[:fh]
	#print("Xl", Xl)
	Xr = X[-sh:]
	#print("Xr", Xr)
	Yl = Y[:fh]
	#print("Yl", Yl)
	Yr = Y[-sh:]
	#print("Yr", Yr)

	P1 = multiply(Xl, Yl)
	P2 = multiply(Xr, Yr)
	P3 = multiply(addBitStrings(Xl, Xr), addBitStrings(Yl, Yr))

	return P1*(1<<(2*sh)) + (P3 - P1 - P2)*(1<<sh) + P2


if __name__ == "__main__":
	num1 = input("Enter in a 64 bit number:\t")
	num2 = input("Enter in another 64 bit number:\t")

	print(multiply(num1, num2))