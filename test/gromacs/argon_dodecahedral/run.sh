# From: http://www3.mpibpc.mpg.de/groups/de_groot/compbio/p1/index.html

#wget http://www3.mpibpc.mpg.de/groups/de_groot/compbio1/p1/argon_start.pdb
#wget http://www3.mpibpc.mpg.de/groups/de_groot/compbio1/p1/argon.top
#wget http://www3.mpibpc.mpg.de/groups/de_groot/compbio1/p1/md.mdp

\rm \#* 
\rm *.gro topol.top
\rm argon_1ns*
\rm traj*
\rm ener*
\rm mdout*
\rm *log
\rm *tpr

gmx grompp -f md.mdp -c argon_start.pdb -p argon.top

gmx mdrun -s topol.tpr -v -c argon.gro -nice 0
