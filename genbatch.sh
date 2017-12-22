#!/bin/bash
cp dijet* olddata/
rm dijet*
./generate_dijet low o 100;
./generate_dijet low o 200