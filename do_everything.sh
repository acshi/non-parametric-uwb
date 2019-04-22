#!/bin/bash
if ! command -v meson ; then
    sudo apt-get install python3 python3-pip ninja-build
    pip3 install --user meson
fi

rm -r builddir

meson builddir
cd builddir
ninja install
cd ../

python3 datasets/dataset_to_binary.py datasets/dataset*.txt
bash collect_evaluations.sh
python3 make_figures.py
