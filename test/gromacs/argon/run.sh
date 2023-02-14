#!/bin/bash -v
# From: http://www3.mpibpc.mpg.de/groups/de_groot/compbio/p1/index.html

#wget http://www3.mpibpc.mpg.de/groups/de_groot/compbio1/p1/argon_start.pdb
#wget http://www3.mpibpc.mpg.de/groups/de_groot/compbio1/p1/argon.top
#wget http://www3.mpibpc.mpg.de/groups/de_groot/compbio1/p1/md.mdp

\rm -f \#* 
\rm -f *.gro topol.top
\rm -f argon_1ns*
\rm -f traj*
\rm -f ener*
\rm -f mdout*
\rm -f *log
\rm -f *tpr
\rm -f *.log
\rm -f *.xtc

# minimal
gmx grompp -f md.mdp -c minimal.pdb -p minimal.top
gmx mdrun -s topol.tpr -v -c minimal.gro -nice 0 -x minimal.xtc -g _minimal.log
grep -C 1 -m 1 "LJ (SR)" _minimal.log > minimal.log

# cubic
gmx grompp -f md.mdp -c cubic.pdb -p argon.top
gmx mdrun -s topol.tpr -v -c cubic.gro -nice 0 -x cubic.xtc -g _cubic.log
grep -C 1 -m 1 "LJ (SR)" _cubic.log > cubic.log

# dodecahedron
gmx editconf -f ./cubic.pdb -o ./dodecahedron.pdb -bt dodecahedron -d -0.23
gmx grompp -f md.mdp -c dodecahedron.pdb -p argon.top
gmx mdrun -s topol.tpr -v -c dodecahedron.gro -nice 0 -x dodecahedron.xtc -g _dodecahedron.log
grep -C 1 -m 1 "LJ (SR)" _dodecahedron.log > dodecahedron.log

# octahedron
gmx editconf -f ./cubic.pdb -o ./octahedron.pdb -bt octahedron -d -0.3
gmx grompp -f md.mdp -c octahedron.pdb -p argon.top
gmx mdrun -s topol.tpr -v -c octahedron.gro -nice 0 -x octahedron.xtc -g _octahedron.log
grep -C 1 -m 1 "LJ (SR)" _octahedron.log > octahedron.log

\rm -f  \#* 
\rm -f  *.gro topol.top
\rm -f  argon_1ns*
\rm -f  traj*
\rm -f  ener*
\rm -f  mdout*
\rm -f  *tpr
\rm -f _*

