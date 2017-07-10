#!/usr/bin/env python
#from Bio import SeqIO
#count = SeqIO.convert('benchmark/input/uniref2017-01/uniref100.xml', 'uniprot-xml', 'benchmark/uniref2017-01.faa', 'fasta')
#print('Converted {} records'.format(count))

import sys
import xml.etree.ElementTree as ET

def processcontent(content):
    root = ET.fromstring(''.join(content))
    seqid = root.get('id')
    n = Tax = TaxID = ''
    for prop in root:
        if prop.tag == 'property' and prop.get('type') == 'member count':
            n = prop.get('value')
        if prop.tag == 'property' and prop.get('type') == 'common taxon':
            Tax = prop.get('value')
        if prop.tag == 'property' and prop.get('type') == 'common taxon ID':
            TaxID = prop.get('value')
    name = root.find('name').text.replace('Cluster: ', '')


    repr = root.find('representativeMember')
    dbref = repr.find('dbReference')
    accession = proteinname = ncbitax = sourceorg = ''
    for prop in dbref:
        if   prop.tag == 'property' and prop.get('type') == 'UniProtKB accession':
            accession = prop.get('value')
        elif prop.tag == 'property' and prop.get('type') == 'protein name':
            proteinname = prop.get('value')
        elif prop.tag == 'property' and prop.get('type') == 'NCBI taxonomy':
            ncbitax = prop.get('value')
        elif prop.tag == 'property' and prop.get('type') == 'source organism':
            sourceorg = prop.get('value')
    
    print('>{} {} n={} Tax={} TaxID= {} RepID={}'.format(seqid, name, n, Tax, TaxID, dbref.get(id)).strip())
    print(repr.find('sequence').text.strip())

content = []
for line in sys.stdin:
    if line.startswith('<entry'):
        content = [line]
    elif line.startswith('</entry>'):
        processcontent(content + [line])
        content = []
    elif content:
        content += [line]

