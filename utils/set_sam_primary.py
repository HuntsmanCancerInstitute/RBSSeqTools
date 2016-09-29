#!/usr/bin/env python

import pysam
import os
import sys
import argparse
import re

def main():
    parser = argparse.ArgumentParser(description="This script sets all mapped reads to primary alignments")
    parser.add_argument("inputBam",help="Alignment file in BAM format")
    parser.add_argument("outputBam",help="Alignment file with primary flags changed")
    
    args = parser.parse_args()
    
    inputBam = pysam.AlignmentFile(args.inputBam,"rb")
    outputBam = pysam.AlignmentFile(args.outputBam,"wb",template=inputBam)
    
    for read in inputBam:
        if not read.is_unmapped and read.is_secondary:
            read.is_secondary = False
        outputBam.write(read)
    
    inputBam.close()
    outputBam.close()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print "user interrupted, exiting"
        sys.exit(1)
