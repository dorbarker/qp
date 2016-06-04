#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import json
import subprocess

def arguments():

    parser = argparse.ArgumentParser()

    return parser.parse_args()

def main():

    args = arguments()

if __name__ == '__main__':
    main()
