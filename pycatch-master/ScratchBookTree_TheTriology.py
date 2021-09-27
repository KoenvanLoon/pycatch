import rpy2
print(type(rpy2.__version__))
print(rpy2.__version__)

print(type('2.9.4'))
print('2.9.4')

if rpy2.__version__ != '2.9.4':
    print("Tested for rpy2 version 2.9.4, current version is", rpy2.__version__)
    print("Please make sure you use the correct version.")
