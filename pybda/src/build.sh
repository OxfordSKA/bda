#!/bin/bash

g++-5 \
    -Wall \
    -Wextra \
    -pedantic \
    -O3 \
    bda.cpp \
    -I/usr/local/include \
    -lcasa_tables \
    -lcasa_casa \
    -lcasa_ms \
    -lcasa_measures \
    -o bda

g++-5 \
    -Wall \
    -Wextra \
    -pedantic \
    -O3 \
    bda_2.cpp \
    -I/usr/local/include \
    -lcasa_tables \
    -lcasa_casa \
    -o bda_2
