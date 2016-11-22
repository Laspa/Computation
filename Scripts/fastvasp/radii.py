
def getRadii():
	radii = {}
	f=open('radii.txt', 'r')
	while True:
		nextLine = f.readline()
		if not nextLine:
			break
		else:
			nextLine = nextLine.split()
			element = nextLine[0]
			radius = nextLine[1]
			radii[element] = float(radius)
	return radii
