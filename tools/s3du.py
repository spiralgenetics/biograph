#!/usr/bin/python3
'''
Display S3 usage as a treemap

First, gather bucket usage data with:

mkdir buckets
for x in `aws s3 ls|cut -f 3 -d ' '`; do aws s3 ls s3://${x} --recursive > buckets/$x; done

Then run this script to turn it into JSON:

s3du.py buckets/* > usage.json

Finally, run Chrome with --allow-file-access-from-files to disable CORS warning:

/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome s3du.html --args --allow-file-access-from-files

D3 treemap code:

http://bl.ocks.org/ganeshv/6a8e9ada3ab7f2d88022
'''
from __future__ import print_function
import sys
import json
from os import path

totals = {}

if not len(sys.argv) > 1:
	raise RuntimeError('Usage: s3du.py bucketfile ...')

import fileinput
for line in fileinput.input():
	bucket = path.basename(fileinput.filename())

	# Ignore invalid entries
	try:
		(_, _, size, fullname) = line.split()
	except ValueError:
		continue

	parts = fullname.split('/', 1)
	if len(parts) > 1:
		subregion = parts[0]
		key = parts[1].split('/', 1)[0]
	else:
		subregion = bucket
		key = parts[0]

	itemkey = '{}/{}'.format(subregion,key)
	if itemkey not in totals:
		totals[itemkey] = {'region': bucket, 'key': key, 'subregion': subregion, 'value': int(size)}
	else:
		prev = totals[itemkey]
		prev['value'] += int(size)
		totals[itemkey] = prev

items = []
summaries = {}
for item in totals:
	entry = totals[item]
	bucket = entry['region']
	# Aggregate all objects < 1GB
	if entry['value'] < 1e9:
		if bucket not in summaries:
			entry['key'] = 'Others'
			summaries[bucket] = entry
		else:
			summaries[bucket]['value'] += entry['value']
	else:
		items.append(totals[item])

for summary in summaries:
	items.append(summaries[summary])

print(json.dumps(items))
