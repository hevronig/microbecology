#!usr/bin/python

import sys
import os
import uuid
from Bio import SeqIO

# path to fasta file
ffile = sys.argv[1]

# get the directory name. in my case, stands for the name of the sample
direc = os.path.dirname(ffile)
# get the basename of the fasta file. 
base = os.path.basename(ffile)
# create a out file path with the appropriate name. in my case, replace the auto generated "Sample" with "rename"
ofile = str(direc).replace('Sample', '../rename')
print(ofile)

def add_sample_and_uniqid(ffile):
    with open(ffile, 'r') as ff, open(ofile+'.fna', 'w') as out:
        records = SeqIO.parse(ff, 'fasta')
        for record in records:
            # generate a unique id based on user name and current time
            uniqid = str((uuid.uuid1()))
            uniqid = uniqid.split('-')[0]
            # the record.id can be changed based on the relevant dir structure ("direc" or "base")
            record.id = uniqid+" "+str(record.id)+'@'+direc
            SeqIO.write(record, out, 'fasta')

def main():
    add_sample_and_uniqid(ffile)

if __name__ == '__main__':
    main()
