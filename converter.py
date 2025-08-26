


c01 = 97
c02 = 97
c12 = 56
c_eff = 1/(1/c01 + 1/c02) + c12
print( c_eff )
c = c_eff
Ec = (1.6e-19)**2/(2*c*1e-15)
h = 6.626e-34
Ec_GHz = Ec/h/1e9
print(Ec, Ec_GHz)

Ej = 0.4967*100

print(Ej)