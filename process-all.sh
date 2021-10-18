#!/bin/bash
set -e

files=(
    18dB/LK_cello8.wav
    18dB/LK_flute8.wav
    18dB/LK_horn8.wav
    18dB/PK_bass8.wav
    18dB/PK_bass16.wav
    18dB/UK_clarinet16.wav
    18dB/UK_flute4.wav
    18dB/UK_flute8.wav
    18dB/UK_flute16.wav
    18dB/UK_oboe8.wav
    18dB/UK_string8.wav
    18dB/UK_trombone16.wav
)

mkdir -p processed
for file in ${files[*]}; do
    ./process.d -i "$file" -o processed
done
