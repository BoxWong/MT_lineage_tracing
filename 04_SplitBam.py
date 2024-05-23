#!/bin/python
import pysam
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Split bam into individual bams')
parser.add_argument('-i', type=str, help='The name of Sam files for parsing')
parser.add_argument('-b', type=str, help='Barcode file for parsing')
args=parser.parse_args()
# file to split
Bam_file = args.i
# where to place output files
out_dir = "./splitted_bam/"
os.system('mkdir -p %s' % (out_dir))
# variable to keep record of barcode index
oldCB = ''
Count = 0

#Read the valid barcode file.
barcode_file = open(args.b, "r")
barcode = []
for i in barcode_file.readlines():
    barcode.append(i.strip())
print(barcode)

# read in upsplit file and loop reads by line
sam_file = pysam.AlignmentFile( Bam_file, "rb")
sam_dict = {}
sam_list = []
no_cb = 0

for read in sam_file.fetch( until_eof=True):#since there is no index for it, until_rof=True is required
    Count += 1
    if "CB:Z" not in read.tostring():
        no_cb += 1
    else:
        CB = read.get_tag('CB')
        if CB == oldCB:
            sam_list.append(read)
        else:
            if oldCB == '':
                sam_list = [read]
                oldCB = CB
            else:
                sam_dict[oldCB] = sam_list
                print("adding %s to dict" %(oldCB))
                sam_list = [read]
                oldCB = CB
#add the last reord to the dict
sam_dict[CB] = sam_list

valid_barcode = [i for i in barcode if i in sam_dict.keys()]

for i in valid_barcode:
    CB_sam_file = pysam.AlignmentFile( out_dir + "{}.bam".format(i), "wb", template=sam_file)
    for j in sam_dict[i]:
        CB_sam_file.write(j)
    CB_sam_file.close()
print("%d alignments with no curated cell barcodes" %(no_cb))
print('%d aligments have been processed and %d have been written into %d separate bam files under %s/splitted_bam' %(Count, Count-no_cb, len(valid_barcode), os.path.abspath('.')))

#CB_sam_file.close()
sam_file.close()
