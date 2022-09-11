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

ffmpeg -y -framerate $fs -pattern_type glob -i "pp_dens/${exp_name}/dens_xz_*.png" -pix_fmt yuv420p -c:v libx264  "pp_dens/${exp_name}_dens.mp4"
ffmpeg -y -framerate $fs -pattern_type glob -i "pp_temp/${exp_name}/temp_xz_*.png" -pix_fmt yuv420p -c:v libx264  "pp_temp/${exp_name}_temp.mp4"

