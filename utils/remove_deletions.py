#!/usr/bin/env python

import re
import os
import argparse
import sys
import pysam

def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description="This script removes reads with 4+ deletions")
    parser.add_argument("inputBam",help="Input alignment file")
    parser.add_argument("outputBam",help="Output alignment file")
    args = parser.parse_args()
    
    
    #open file handles
    infile = pysam.Samfile(args.inputBam)
    outfile = pysam.Samfile(args.outputBam,"wb",template=infile)
    
    total = 0
    ok = 0
    removed = 0
    
    for read in infile:
        total += 1
        if total % 1000000 == 0:
            print total, ok, removed
        
        read_pass = True
        cigar = read.cigartuples
        for o,l in cigar:
            if o == 2:
                if l >= 4:
                    read_pass = False
                    removed += 1
        if read_pass:
            ok += 1
            outfile.write(read)
    infile.close()
    outfile.close()

    print total, ok, removed

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print "user interrupted, exiting"
        sys.exit(1)
