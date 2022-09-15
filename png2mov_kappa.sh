#! /usr/bin/bash


if [ $BASH_ARGC -ge 1 ]; then
	exp_name=$1
else
	exp_name="exp1"
fi

if [ $BASH_ARGC -ge 2 ]; then
        fs=$2
else
        fs=8
fi

if [ $BASH_ARGC -ge 3 ]; then
	kappa_list=$3
else
	kappa_list="0.1 0.2 0.3 0.06 1.5"
fi

if [ $BASH_ARGC -ge 4 ]; then
	angle_list=$4
else
	angle_list="10 30 50 80"
fi

for kappa in $kappa_list; do 
	cd kappa_${kappa}/
for ang in $angle_list ; do
	ffmpeg -y -framerate $fs -pattern_type glob -i "${ang}/kappa_${kappa}_los_${ang}_dens_*.png" -pix_fmt yuv420p -c:v libx264  ${exp_name}_k_${kappa}_a_${ang}_dens.mp4
	ffmpeg -y -framerate $fs -pattern_type glob -i "${ang}/kappa_${kappa}_los_${ang}_temp_*.png" -pix_fmt yuv420p -c:v libx264  ${exp_name}_k_${kappa}_a_${ang}_temp.mp4
done
	cd ..
done
# for k in 0.1 0.2 0.3 0.06 1.5; do for p in 10 30 50 80; do echo mkdir -p kappa_$k/$p/ ; done; done;
# ffmpeg -framerate $fs -i kappa_${kappa}_los_${ang}_dens_%4d.png  -pix_fmt yuv420p -c:v libx264  ${exp_name}_k_${kappa}_a_${ang}.mp4
