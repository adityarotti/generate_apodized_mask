#!/bin/bash
./run-gen-apo-mask < submit/gal_ps > out/gal_ps.out &
./run-gen-apo-mask < submit/gal_ps_cl > out/gal_ps_cl.out &
./run-gen-apo-mask < submit/gal_ps_ccl > out/gal_ps_ccl.out &
