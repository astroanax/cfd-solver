#!/usr/bin/env bash

rm -rf build
mkdir build

g++ -std=c++17 -I. -c -O2 euler.cpp -o build/euler.o
g++ build/euler.o -o build/euler

cd build
./euler
cd ..
