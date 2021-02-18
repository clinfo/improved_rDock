#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import sys

if __name__ == '__main__':

    args = sys.argv
    if len(args) != 2 :
        print 'Usage: python %s gridfile' % args[0]
        quit()

    f = open(args[1], 'r')
    line = f.readline()
    line = f.readline()
    line = f.readline()
    tokens = line.split()
    a = float(tokens[0])
    b = float(tokens[1])
    c = float(tokens[2])

    line = f.readline()
    tokens = line.split()
    nx = int(tokens[0])
    ny = int(tokens[1])
    nz = int(tokens[2])

    line = f.readline()
    tokens = line.split()
    ixs = int(tokens[1])
    ixe = int(tokens[2])
    iys = int(tokens[3])
    iye = int(tokens[4])
    izs = int(tokens[5])
    ize = int(tokens[6])
    
    dx = a / nx
    dy = b / ny
    dz = c / nz

    for iz in range(izs, ize + 1):
        for iy in range(iys, iye + 1):
            for ix in range(ixs, ixe + 1):
                x = ix * dx
                y = iy * dy
                z = iz * dz
                line = f.readline()
                v = float(line)
                if 0.0 < v:
                    print x, y, z, v
                
    f.close()
