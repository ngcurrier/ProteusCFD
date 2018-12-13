#!/bin/bash

# This format roughly follows the tutorial on using snappyHexMesh from
# http://sourceflux.de/blog/preparing-geometry/

admesh prop.stl

surfaceCheck prop.stl

surfaceOrient prop.stl "(1e10 1e10 1e10)" propFixed.stl

surfaceCheck propFixed.stl | grep -E '^Bounding Box' > BoundingBox.txt
