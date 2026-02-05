#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser("one2multi_filter", description="")
parser.add_argument("-m", "--mapping", required=True, help="one2multi mapping")
parser.add_argument("-f", "--file", required=True, help="input file")
parser.add_argument("-1", "--c1", required=True, type=int, help="which column (one)")
parser.add_argument("-2", "--c2", required=True, type=int, help="which column (multi)")
args = parser.parse_args()


mapping = dict()
for line in open(args.mapping, "r"):
	line = line.strip()
	if not line:
		continue
	if line.startswith("#"):
		continue

	key, value = line.split("\t")
	mapping[key] = value.split(",")

for line in sys.stdin if args.file in ("stdin", "-") else open(args.file, "r"):
	line = line.strip()
	if not line:
		continue
	if line.startswith("#"):
		continue

	content = line.split("\t")
	try:
		if content[5] in mapping.keys() and content[0] in mapping[content[5]]:
			print(line)   
	except Exception as e:
		content = line.split(" ") 
		if content[5] in mapping.keys() and content[0] in mapping[content[5]]:
			print(line)   
