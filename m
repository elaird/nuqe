#!/bin/bash

g++ -o foo src/foo.C src/TT_*.C `root-config --cflags --libs` $*
