import numpy

x = numpy.float64(10**9)

calc = numpy.float64(1.0) / (x - numpy.sqrt(x**2.0-1.0)) # Divide by zero encountered
print(calc) # inf

# Tak bliskie X, że prezyzja Float64 zaokrągla wynik do X 
# print("Sqrt", numpy.sqrt(x**2.0-1.0))

# Wolfram alpha
# 1_999_999_999.9999999994999999