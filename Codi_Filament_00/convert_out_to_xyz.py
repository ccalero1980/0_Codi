import sys
numin = int(sys.argv[1])
print numin
nparts = 30


nameout = 'trajectory.xyz'
fout = open(nameout,'w')
fout.write(str(nparts) + '\n\n')
fin = open('straight30', 'r')
for k in range(nparts):
  line = fin.readline()
  x=float(line.split()[0])
  y=float(line.split()[1])
  z=float(line.split()[2])
  fout.write('A  %14.8f  %14.8f  %14.8f\n'%(x,y,z))

for n in range(1,numin):
  fout.write(str(nparts) + '\n\n')
  namein = 'out_%i.00'%n
  fin = open(namein, 'r')
  for k in range(nparts):
    line = fin.readline()
    x=float(line.split()[0])
    y=float(line.split()[1])
    z=float(line.split()[2])
  
    fout.write('A  %14.8f  %14.8f  %14.8f\n'%(x,y,z))
 
