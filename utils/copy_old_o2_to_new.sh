#!/bin/bash
h5copy -i $1 -o NewO2Standard$1 -s po2 -d po2
h5copy -i $1 -o NewO2Standard$1 -s recomputed_flow/out0006 -d recomputed_flow/vessels -p
h5copy -i $1 -o NewO2Standard$1 -s po2/out0006/parameters -d parameters/o2 -p 
h5copy -i $1 -o NewO2Standard$1 -s po2/out0006/parameters/D_tissue -d parameters/o2/D_plasma
h5copy -i $1 -o NewO2Standard$1 -s po2/out0006/parameters/S_n -d parameters/o2/sat_curve_exponent
h5copy -i $1 -o NewO2Standard$1 -s po2/out0006/parameters/S_p50 -d parameters/o2/sat_curve_p50
