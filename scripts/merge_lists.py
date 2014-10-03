#! /usr/bin/env python

import argparse

""" Import PSR names from local data path.
"""
parser = argparse.ArgumentParser()
parser.add_argument("list1")
parser.add_argument("list2")
parser.add_argument("-s", "--savepath", default='config/psrlist.txt')

args = parser.parse_args()

with open(args.list1, 'r') as f:
    psrs1 = f.readlines()

with open(args.list2, 'r') as f:
    psrs2 = f.readlines()

psrlist = set(psrs1 + psrs2)

with open(args.savepath, 'w') as f:
    for psr in psrlist:
        f.write(psr)

print "Merged list saved to %r" % args.savepath