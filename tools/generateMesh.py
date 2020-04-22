#!/usr/bin/env python3

if __name__ == "__main__":
    import sys, subprocess
    h = 1.0

    if len(sys.argv) != 3:
        print('USAGE: %s <STL filename> <Farfield distance multiplier>' % sys.argv[0])
        print('NOTE: that this script makes use of gmsh and admesh for backend tasks')
        print('both tools must be in your path and accessible via command line call')
        exit()

    filenameSTL = sys.argv[1]
    FFdist = float(sys.argv[2])
    print('STL file: %s' % filenameSTL)
    print('Farfield multiplier: %f' % FFdist)
    
    # make sure our STL is in ascii
    cmd = 'admesh ' + filenameSTL + ' -a ' + filenameSTL
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    
    output = str(p.communicate()[0])
    output = output.split('\\n')
    xmin = 0
    xmax = 0
    ymin = 0
    ymax = 0
    zmin = 0
    zmax = 0
    
    for line in output:
        # find x extents line
        if line.find('Min X') != -1:
            tok = line.split(',')
            xmin = float(tok[0].split('=')[1])
            xmax = float(tok[1].split('=')[1])
        # find y extents line
        elif line.find('Min Y') != -1:
            tok = line.split(',')
            ymin = float(tok[0].split('=')[1])
            ymax = float(tok[1].split('=')[1])
        # find z extents line
        elif line.find('Min Z') != -1:
            tok = line.split(',')
            zmin = float(tok[0].split('=')[1])
            zmax = float(tok[1].split('=')[1])


    xrange = xmax - xmin
    yrange = ymax - ymin
    zrange = zmax - zmin

    maxrange = max(max(xrange, yrange), zrange)
    print('Xrange: %f \tto\t %f' % (xmin, xmax))
    print('Yrange: %f \tto\t %f' % (ymin, ymax))
    print('Zrange: %f \tto\t %f' % (zmin, zmax))
    print('Maximum range detected: %f' % maxrange)
    offsetDist = maxrange*FFdist
    print('Mesh extents from body: %f' % (offsetDist))


    h = maxrange
    
    # generate mesh extent arrays
    x = []
    x.append(xmin - offsetDist)
    x.append(xmax + offsetDist)
    y = []
    y.append(ymin - offsetDist)
    y.append(ymax + offsetDist)
    z = []
    z.append(zmin - offsetDist)
    z.append(zmax + offsetDist)
    
    print(x,y,z)

    fname = filenameSTL.replace('.stl', '')
    fout = open(fname+'.geo', 'w')

    fout.write('Merge \"%s\";\n' % filenameSTL)
    # write corner points
    # minus x plane
    fout.write('Point(1) = {%f, %f, %f, %f};\n' % (x[0], y[0], z[0], h))
    fout.write('Point(2) = {%f, %f, %f, %f};\n' % (x[0], y[0], z[1], h))
    fout.write('Point(3) = {%f, %f, %f, %f};\n' % (x[0], y[1], z[1], h))
    fout.write('Point(4) = {%f, %f, %f, %f};\n' % (x[0], y[1], z[0], h))
    # pos x plane
    fout.write('Point(5) = {%f, %f, %f, %f};\n' % (x[1], y[0], z[0], h))
    fout.write('Point(6) = {%f, %f, %f, %f};\n' % (x[1], y[0], z[1], h))
    fout.write('Point(7) = {%f, %f, %f, %f};\n' % (x[1], y[1], z[1], h))
    fout.write('Point(8) = {%f, %f, %f, %f};\n' % (x[1], y[1], z[0], h))

    # write lines
    # minus x plane
    fout.write('Line(1) = {1,2};\n')
    fout.write('Line(2) = {2,3};\n')
    fout.write('Line(3) = {3,4};\n')
    fout.write('Line(4) = {4,1};\n')
    # pos x plane
    fout.write('Line(5) = {5,6};\n')
    fout.write('Line(6) = {6,7};\n')
    fout.write('Line(7) = {7,8};\n')
    fout.write('Line(8) = {8,5};\n')
    # connection planes
    fout.write('Line(9) = {1,5};\n')
    fout.write('Line(10) = {2,6};\n')
    fout.write('Line(11) = {3,7};\n')
    fout.write('Line(12) = {4,8};\n')
    
    # write surfaces
    fout.write('Curve Loop(1) = {1,2,3,4};\n')
    fout.write('Plane Surface(2) = {1};\n')
    fout.write('Curve Loop(2) = {5,6,7,8};\n')
    fout.write('Plane Surface(3) = {2};\n')
    fout.write('Curve Loop(3) = {2,11,-6,-10};\n')
    fout.write('Plane Surface(4) = {3};\n')
    fout.write('Curve Loop(4) = {3,12,-7,-11};\n')
    fout.write('Plane Surface(5) = {4};\n')
    fout.write('Curve Loop(5) = {9,-8,-12,4};\n')
    fout.write('Plane Surface(6) = {5};\n')
    fout.write('Curve Loop(6) = {5,-10,-1,9};\n')
    fout.write('Plane Surface(7) = {6};\n')

    # write volume
    fout.write('Surface Loop(1) = {1};\n') # this is body
    fout.write('Surface Loop(2) = {2,7,3,4,5,6};\n') # these are farfield
    fout.write('Volume(1) = {1,2};\n')

    # write physical tags
    fout.write('Physical Surface("Body") = {1};\n')
    fout.write('Physical Surface("FarfieldXmin") = {2};\n')
    fout.write('Physical Surface("FarfieldXmax") = {3};\n')
    fout.write('Physical Surface("FarfieldZmax") = {4};\n')
    fout.write('Physical Surface("FarfieldYmax") = {5};\n')
    fout.write('Physical Surface("FarfieldZmin") = {6};\n')
    fout.write('Physical Surface("FarfieldYmin") = {7};\n')

    fout.write('Physical Volume("volume") = {1};\n')
    
    fout.close()

    # generate the mesh with gmsh
    print('Generating mesh with GMSH')
    print('-------------------------------------')
    cmd = 'gmsh ' + fname+'.geo' + ' -3'
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output = str(p.communicate()[0])
    output = output.split('\\n')
    print('Mesh completed: %s' % fname+'.msh')

    print('Running decomp')
    print('-------------------------------------')
    cmd = 'udecomp.x ' + fname+'.msh' + ' 4'
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output = str(p.communicate()[0])
    output = output.split('\\n')
    print('Completed decomp')

    print('Running recomp')
    print('-------------------------------------')
    cmd = 'urecomp.x ' + fname
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output = str(p.communicate()[0])
    output = output.split('\\n')
    print('Completed recomp')
    

    print('Writing bc file')
    print('-------------------------------------')
    fout = open(fname+'.bc','w')

    fout.write('#BC file... comments start with a # symbol\n')
    fout.write('surface #1 = noSlip "Body" twall = [1.0] bleedSteps = [30]\n')
    fout.write('surface #2 = farField "Xmin"\n')
    fout.write('surface #3 = farField "Xmax"\n')
    fout.write('surface #4 = farField "Zmax"\n')
    fout.write('surface #5 = farField "Ymax"\n')
    fout.write('surface #6 = farField "Zmin"\n')
    fout.write('surface #7 = farField "Ymin"\n')
    fout.write('\nbody #1 = [1]\n')
    
