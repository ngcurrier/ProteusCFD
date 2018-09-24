#!/bin/bash

# This fixes weird Debian installation behavior that tries to avoid shadowing
# parallel HDF5 with serial but really just causes problems
ln -s /usr/lib/x86_64-linux-gnu/libhdf5_serial.a /usr/lib/x86_64-linux-gnu/libhdf5.a
ln -s /usr/lib/x86_64-linux-gnu/libhdf5_serial.so /usr/lib/x86_64-linux-gnu/libhdf5.so
