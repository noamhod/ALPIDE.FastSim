import math
me = 0.5109989461/1000.
particles = []

fIn = open("../bppp_000039.out",'r')
for line in fIn:
   if(line.startswith("#")): continue
   words = line.split()
   words = [word for word in words if(len(word.split())!=0)]
   if( abs(int(words[7]))!=11 ): continue
   particles.append( {"E":float(words[0]), "vx":float(words[1]), "vy":float(words[2]), "vz":float(words[3]), "betax":float(words[4]), "betay":float(words[5]), "betaz":float(words[6]), "pdgId":int(words[7]), "wgt":float(words[8])} )
fIn.close()
print("nparticles=%d" % len(particles))

def getE(bx,by,bz,E0):
   b2 = bx*bx+by*by+bz*bz
   gamma = 1./math.sqrt(1.-b2)
   E = gamma*me
   print("E=%8f, E0=%8f" % (E,E0))

for p in particles:
   getE( p["betax"], p["betay"], p["betaz"], p["E"] )
