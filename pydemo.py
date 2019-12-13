from ctypes import *

mylib=CDLL("E:/Users/Andrew/eclipse-workspace/clearn/libtest.so")
print(mylib.add_int(5,5))

addFloat=mylib.add_float
addFloat.restype=c_float

print(addFloat(c_float(1.2),c_float(2.3)))

# transfer list to c lib
sum=mylib.sum
sum.restype=c_float
floatArray5=(c_float * 5)    # define a class for array with specific type and length
testArray=floatArray5(1.1,2.2,3.3,4.4,5.5) # to object
print(sum(testArray,5))