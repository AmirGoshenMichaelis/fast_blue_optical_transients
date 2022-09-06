#! /usr/bin/bash


if [ $BASH_ARGC -ge 1 ]; then
	kappa=$1
else
	kappa=0.3
fi

if [ $BASH_ARGC -ge 2 ]; then
        fs=$1
else
        fs=8
fi

for kappa in 0.1; do 
	cd kappa_${kappa}/
for ang in 10 30 50 80 ; do
	ffmpeg -y -framerate $fs -pattern_type glob -i "${ang}/kappa_${kappa}_dens_photo_${ang}_PN1_*.png" k_${kappa}_a_${ang}.mp4
done
	cd ..
done
# for k in 0.1 0.2 0.3 0.06 1.5; do for p in 10 30 50 80; do echo mkdir -p kappa_$k/$p/ ; done; done;
# for k in 0.1 0.2 0.3 0.06 1.5; do for p in 10 30 50 80; do echo kappa_${k}_dens_photo_${p}_PN1_*.png kappa_$k/$p/ ; done; done;

# ffmpeg -framerate $fs -i kappa_${kappa}_dens_photo_${ang}_PN1_%4d.png k_${kappa}_a_${ang}.mp4
