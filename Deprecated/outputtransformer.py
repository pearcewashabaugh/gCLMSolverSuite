#This code takes the ouput of main.py, uout.txt and transforms it into 
#something readily readable by mathematica, uout2.txt

with open('Output/uout.txt', 'r') as f:
	ustr = f.read()

ustr = ustr.split('\n')

ustr = ustr[:-1]

for iindex, i in enumerate(ustr):
	ustr[iindex] = i.split(' ')

for iindex, i in enumerate(ustr):
	for jindex, j in enumerate(i):
		ustr[iindex][jindex] = float(j)

ustr = str(ustr)
ustr = ustr.replace("[","{")
ustr = ustr.replace("]","}")

g = open('Output/uout2.txt', 'w')
g.write(ustr)

f.close()
g.close()

###############################################################################

with open('Output/eta.txt', 'r') as f:
	etastr = f.read()

etastr = etastr.split('\n')

etastr = etastr[:-1]

for iindex, i in enumerate(etastr):
	etastr[iindex] = i.split(' ')

for iindex, i in enumerate(etastr):
	for jindex, j in enumerate(i):
		etastr[iindex][jindex] = float(j)

etastr = str(etastr)
etastr = etastr.replace("[","{")
etastr = etastr.replace("]","}")

g = open('Output/eta2.txt', 'w')
g.write(etastr)

f.close()
g.close()

##############################################################################

with open('Output/etadx.txt', 'r') as f:
	etadxstr = f.read()

etadxstr = etadxstr.split('\n')

etadxstr = etadxstr[:-1]

for iindex, i in enumerate(etadxstr):
	etadxstr[iindex] = i.split(' ')

for iindex, i in enumerate(etadxstr):
	for jindex, j in enumerate(i):
		etadxstr[iindex][jindex] = float(j)

etadxstr = str(etadxstr)
etadxstr = etadxstr.replace("[","{")
etadxstr = etadxstr.replace("]","}")

g = open('Output/etadx2.txt', 'w')
g.write(etadxstr)

f.close()
g.close()
###############################################################################

with open('Output/wout.txt', 'r') as f:
	wstr = f.read()

wstr = wstr.split('\n')

wstr = wstr[:-1]

for iindex, i in enumerate(wstr):
	wstr[iindex] = i.split(' ')

for iindex, i in enumerate(wstr):
	for jindex, j in enumerate(i):
		wstr[iindex][jindex] = float(j)

wstr = str(wstr)
wstr = wstr.replace("[","{")
wstr = wstr.replace("]","}")

g = open('Output/wout2.txt', 'w')
g.write(wstr)

f.close()
g.close()

###############################################################################

with open('Output/etaxint.txt', 'r') as f:
	etaxintstr = f.read()

etaxintstr = etaxintstr.split('\n')

etaxintstr = etaxintstr[:-1]

for iindex, i in enumerate(etaxintstr):
	etaxintstr[iindex] = i.split(' ')

for iindex, i in enumerate(etaxintstr):
	for jindex, j in enumerate(i):
		etaxintstr[iindex][jindex] = float(j)

etaxintstr = str(etaxintstr)
etaxintstr = etaxintstr.replace("[","{")
etaxintstr = etaxintstr.replace("]","}")

g = open('Output/etaxint2.txt', 'w')
g.write(etaxintstr)

f.close()
g.close()