import json
import logging
import pysam

import numpy as np 
import pandas as pd

from .fusion import Fusion
from .region import Region
from .utils import read_json, update_list

# Import
def read_bed(filename):
    """Read a bed file. Return a dataframe"""
    isHeader = 0
    with open(filename) as fid:
        if 'track' in fid.readline():
            isHeader=1
    bed = pd.read_csv(filename, sep='\t', usecols=[0,1,2], names=["chrom", "start", "end"], dtype={"chrom":str, "start":int, "end":int}, skiprows=isHeader)
    return bed

def parse_hmnfusion_json(filename):
    """Read json file construct with exctractfusion command. Return a list of Fusion"""
    data = read_json(filename)
    fusions = []
    for key in sorted(data["fusions"].keys()):
        fusions.append(Fusion.from_dict(data["fusions"][key]))
    return fusions

# Helper functions
def _cigar2position(cigars, start):
    """Construct from a cigar and a position, a position/operation"""
    data = {}
    if cigars[0][0] in [4, 5]:
        for i in range(cigars[0][1]):
            i += 1
            data[start-i] = cigars[0][0]
        del cigars[0]
    for cigar in cigars:
        op, nb = cigar[0], cigar[1]
        if op in [4, 5]:
            for i in range(nb):
                data[start+i] = op
                i += 1
        elif op == 0:
            for i in range(nb):
                data[start] = op
                start+=1
        elif op == 1:
            pass
        else:
            start += 1
    return data


def _select_bed(x, region):
    """Select from a bed file, regions cross over an other region"""
    if x['chrom'] == region.chrom:
        if x['start'] <= region.position and x['end'] >= region.position:
            return True
    return False
        

def run(params, bed, fusions):
    """Main function to quantify fusion"""
    alignment = pysam.AlignmentFile(params['falignment']['path'], params['falignment']['mode'])

    to_delete = []
    for ix, fusion in enumerate(fusions):
        isSkip = False
        if not fusion.isConsensus:
            continue
        # Check fusion against bed.
        sub_first, sub_second = pd.DataFrame(columns=bed.columns), pd.DataFrame(columns=bed.columns)
        if fusion.first.is_init():
            sel = bed.apply(_select_bed, axis=1, args=(fusion.first,))
            sub_first = bed[sel]
        if fusion.second.is_init():
            sel = bed.apply(_select_bed, axis=1, args=(fusion.second,))
            sub_second = bed[sel]
        if len(sub_first) > 1 or len(sub_second) > 1:
            logging.warning('Fusion %s is found multiple times in bed -> skipping'%(fusion,))
            isSkip = True
        if len(sub_first) + len(sub_second) == 2 :
            logging.warning('Fusion %s is found on left and right of breakpoint in the bed -> skipping'%(fusion,))
            isSkip = True
        if len(sub_first) + len(sub_second) == 0 :
            logging.warning("Fusion %s isn't found on left or right of breakpoint in the bed -> skipping"%(fusion,))
            isSkip = True

        if isSkip:
            to_delete.append(ix)
            continue
        # Init.
        bed_sel = pd.DataFrame(columns=bed.columns)
        region = Region()
        if len(sub_first) == 1:
            bed_sel = sub_first
            region = fusion.first
        elif len(sub_second) == 1:
            bed_sel = sub_second
            region = fusion.second
            # Swap.
            fusion.swap_region()
        else:
            logging.warning("Fusion %s, something bad happened -> skipping"%(fusion,))
        count = dict(coverage=0, split=0, mate=0, clipped=0)
        # Run. 
        for aligned_segment in alignment.fetch(bed_sel.iloc[0, 0], bed_sel.iloc[0, 1], bed_sel.iloc[0, 2]):
            # Filtering.
            if aligned_segment.is_unmapped or aligned_segment.is_duplicate or aligned_segment.is_supplementary:
                continue

            cigar2pos = _cigar2position(aligned_segment.cigartuples, aligned_segment.reference_start)                
            if not region.position in cigar2pos.keys(): 
                continue

            count['coverage'] += 1
            # Count split reads.
            if aligned_segment.has_tag("SA"):
                count['split'] += 1
                continue
            
            # Count other Chrom.
            if aligned_segment.is_paired:    
                if not aligned_segment.mate_is_unmapped and not aligned_segment.is_unmapped:
                    if aligned_segment.next_reference_id != aligned_segment.reference_id:
                        count['mate'] += 1
                        continue
        
            # Count reads clipped.
            count_clipped = np.zeros((2, params['clipped']['interval']))
            for i in range(params['clipped']['interval']):
                if cigar2pos.get(region.position-i-1,0) in [4, 5]:
                    count_clipped[0][i] = 1
                if cigar2pos.get(region.position+i+1,0) in [4, 5]:
                    count_clipped[1][i] = 1
            
            if np.max(np.sum(count_clipped, axis=1)) >= params['clipped']['count']:
                count['clipped'] += 1
        fusion.depth = count['coverage']
        fusion.evidence = count['split'] + count['mate'] + count['clipped']
        fusion.evidence_details = count
    fusions = update_list(fusions, to_delete)
    return fusions

# Write.
def _get_header():
    """Provide header to build vcf file"""
    header = []
    header.append('##fileformat=VCFv4.2')
    header.append('##source=HmnFusion')
    header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">')
    header.append('##INFO=<ID=SOFT,Number=1,Type=String,Description=\"Indicated from which software is derived\">')
    header.append('##INFO=<ID=FROM,Number=.,Type=String,Description=\"Indicated from which reference is derived\">')
    header.append('##INFO=<ID=CONS,Number=.,Type=String,Description=\"Is a consensus\">')
    header.append('##INFO=<ID=VAF,Number=.,Type=Float,Description=\"Allelic frequence observed\">')
    header.append('##INFO=<ID=DP,Number=.,Type=Integer,Description=\"Approximate read depth across all samples\">')
    header.append('##INFO=<ID=SU,Number=.,Type=Integer,Description=\"Number of pieces of evidence supporting the variant across all samples\">')
    header.append('##INFO=<ID=PE,Number=.,Type=Integer,Description=\"Number of paired-end reads supporting the variant across all samples\">')
    header.append('##INFO=<ID=SR,Number=.,Type=Integer,Description=\"Number of split reads supporting the variant across all samples\">')
    header.append('##INFO=<ID=SC,Number=.,Type=Integer,Description=\"Number of soft clipped reads supporting the variant across all samples\">')

    header.append('##ALT=<ID=FUS,Description=\"Fusion\">')

    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
    header.append('##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Allelic frequence observed\">')
    header.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">')
    header.append('##FORMAT=<ID=SU,Number=1,Type=Integer,Description=\"Number of pieces of evidence supporting the variant\">')
    header.append('##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Number of paired-end reads supporting the variant\">')
    header.append('##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of split reads supporting the variant\">')
    header.append('##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"Number of soft clipped reads supporting the variant\">')

    return '\n'.join(header)

def write(filename, name, fusions):
    """Write a vcf file from a list of Fusion"""    
    data = {}

    # Header.
    with open(filename, 'w') as fod:
        fod.write(_get_header() + '\n')

    # Fusions.
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [name]
    df = pd.DataFrame(columns=columns)

    for ix, fusion in enumerate(fusions):
        ix += 1

        logging.debug(repr(fusion))
        # First.
        ident = fusion.ident
        ident_1 = ident
        if fusion.second.is_init():
            ident_1 += '_1'
        infos = ['SVTYPE=FUS']
        infos += ['SOFT=%s'%('-'.join(fusion.software),)]
        infos += ['FROM=%s'%('-'.join(fusion.buildFrom),)]
        infos += ['CONS=%s'%(fusion.isConsensus,)]
        infos += ['VAF=%s'%(fusion.get_vaf(),)]
        infos += ['DP=%s'%(fusion.depth,)]
        infos += ['SU=%s'%(fusion.evidence,)]
        ed = fusion.evidence_details
        infos += ['SR=%s'%(ed.get('split',0),)]
        infos += ['PE=%s'%(ed.get('mate',0),)]
        infos += ['SC=%s'%(ed.get('clipped',0),)]
    
        infos = ':'.join(infos)
        values = [fusion.first.chrom, fusion.first.position, ident_1, 'N', '<FUS>', '.', '.', infos, 'GT:VAF:DP:SU:SR:PE:SC', './.:%s:%s:%s:%s:%s'%(fusion.depth, fusion.get_vaf(), fusion.evidence, ed.get('split',0), ed.get('mate',0), ed.get('clipped',0))]
        df = df.append(pd.Series(values, index=columns), ignore_index=True)

        ident_2 = ident
        if fusion.second.is_init():
            ident_2 += '_2'
            # Second.
            infos = ';'.join(['SVTYPE=FUS', 'DP=.', 'SU=.'])
            values = [fusion.second.chrom, fusion.second.position, ident_2, 'N', '<FUS>', '.', '.', infos, 'GT:VAF:DP:SU:SR:PE:SC', './.:.:.']
            df = df.append(pd.Series(values, index=columns), ignore_index=True)
    df.to_csv(filename, mode='a', sep='\t', index=False)
