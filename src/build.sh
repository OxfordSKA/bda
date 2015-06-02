#!/bin/bash

g++-5 \
    -Wall \
    -pedantic \
    -O3 \
    bda.cpp \
    -I/usr/local/include \
    -lcasa_tables \
    -lcasa_casa \
    -lcasa_ms \
    -lcasa_measures \
    -o bda
