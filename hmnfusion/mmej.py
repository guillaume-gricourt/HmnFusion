import hashlib
import logging
import os
import pysam
import sys

import pandas as pd

def extract(finputs, event):
    """Extract event of interest from vcf input files"""
    df = pd.DataFrame()
    for finput in finputs:
        vcf_in = pysam.VariantFile(finput)
        for record in vcf_in.fetch():
            # Check if variant is a deletion.
            if event in ['all', 'deletion'] and len(record.ref) > len(record.alts[0]):
                region_name = '%s %s %s %s'%(record.contig, record.pos, record.ref,
                record.alts[0])
                region_id = hashlib.md5(region_name.encode('utf8')).hexdigest()
                # Check if variant is already seen.
                if region_id in df.index:
                    continue

                df.at[region_id, 'event'] = 'deletion'
                df.at[region_id, 'contig'] = record.contig
                df.at[region_id, 'start'] = record.pos
                df.at[region_id, 'deletion'] = record.ref[1:]

                for sample in record.samples.keys():
                    df.at[region_id, sample] = True

    df['start'] = df['start'].astype(int)
    return df

def signatures(freference, df):
    """Identify MH motif from event"""
    def _signatures(x, freference):
        len_deletion = len(x['deletion'])
        start = int(x['start'])
        region = '%s:%s-%s'%(x['contig'], start, start+(2*len_deletion))
        faidx = pysam.faidx(freference, region, split_lines=True)
        seq = faidx[1]

        left = rec.seq[:len_deletion]
        right = rec.seq[len_deletion:]

        mh_len, mh_seq = 0, ''
        for i in range(len_deletion):
            motif = left[:i]
            if right.startswith(motif):
                mh_len = i
                mh_seq = motif
            if i > mh_len:
                break

        return mh_seq

    df['mmej_sequence'] = df.apply(_signatures, axis=1, args=(freference,))
    return df

def conclude(df):
    """Conclude about the presens of MMEJ signature"""
    def _conclude(x):
        res = ''
        len_deletion = len(x['alt'])
        len_mh = len(x['mmej_sequence'])
        if len_deletion <= 1 or len_mh <= 1:
            res = 'MH deletion'
        elif len_deletion > 1 and len_mh < 5:
            res = 'no clear signature'
        elif len_deletion > 1 and len_mh >= 5:
            res = 'mmej signature'
        return res

    df['mmej_conclusion'] = df.apply(_conclude, axis=1)
    return df

# Write.
def write(filename, finputs, fusions):
    """Write list of fusion to a json file"""
    data = {}
    data['inputs'] = finputs
    data['fusions'] = {}
    for ix, fusion in enumerate(fusions['genefuse']['raw'] + fusions['lumpy']['raw'] + fusions['consensus']):
        data['fusions'][str(ix)] = fusion.to_dict()
    
    with open(filename, 'w') as fod:
        json.dump(data, fod)
