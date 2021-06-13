import json
import logging
import pysam

import numpy as np 
import pandas as pd

from . import (evidence, fusion, graph, region)

# Import
def read_bed(filename):
    """Read a bed file. Return a dataframe"""
    isHeader = 0
    with open(filename) as fid:
        if 'track' in fid.readline():
            isHeader=1
    bed = pd.read_csv(filename, sep='\t', usecols=[0,1,2], names=["chrom", "start", "end"], dtype={"chrom":str, "start":int, "end":int}, skiprows=isHeader)
    return bed

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

def run(params, bed, graph):
    """Main function to quantify fusion"""
    alignment = pysam.AlignmentFile(params['falignment']['path'], params['falignment']['mode'])

    nodes = graph.nodes
    for n in nodes:
        graph.nodes[n]['is_skip'] = False

        if not graph.nodes[n]['is_consensus']:
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
            logging.warning('Fusion %s is found multiple times in bed -> skipping'%(graph.nodes[n]['fusion'],))
            graph.nodes[n]['is_skip'] = True
        if len(sub_first) + len(sub_second) == 2 :
            logging.warning('Fusion %s is found on left and right of breakpoint in the bed -> skipping'%(graph.nodes[n]['fusion'],))
            graph.nodes[n]['is_skip'] = True
        if len(sub_first) + len(sub_second) == 0 :
            logging.warning("Fusion %s isn't found on left or right of breakpoint in the bed -> skipping"%(graph.nodes[n]['fusion'],))
            graph.nodes[n]['is_skip'] = True

        if graph.nodes[n]['is_skip']:
            continue
        # Init.
        bed_sel = pd.DataFrame(columns=bed.columns)
        region = region.Region()
        if len(sub_first) == 1:
            bed_sel = sub_first
            region = graph.nodes[n]['fusion'].first
        elif len(sub_second) == 1:
            bed_sel = sub_second
            region = graph.nodes[n]['fusion'].second
            # Swap.
            graph.nodes[n]['fusion'].swap_region()
        else:
            logging.warning("Fusion %s, something bad happened -> skipping"%(graph.nodes[n]['fusion'],))
        count = dict(coverage=0, split=0, mate=0, clipped=0)
        # Run. 
        for aligned_segment in alignment.fetch(bed_sel.iloc[0, 0], bed_sel.iloc[0, 1], bed_sel.iloc[0, 2]):
            # Filtering.
            if aligned_segment.is_unmapped or aligned_segment.is_duplicate or aligned_segment.is_supplementary:
                continue

            cigar2pos = _cigar2position(aligned_segment.cigartuples, aligned_segment.reference_start)                
            if not region.position in cigar2pos.keys(): 
                continue

            graph.nodes[n]['fusion'].evidence.depth += 1
            # Count split reads.
            if aligned_segment.has_tag("SA"):
                graph.nodes[n]['fusion'].evidence.split += 1
                continue
            
            # Count other Chrom.
            if aligned_segment.is_paired:    
                if not aligned_segment.mate_is_unmapped and not aligned_segment.is_unmapped:
                    if aligned_segment.next_reference_id != aligned_segment.reference_id:
                        graph.nodes[n]['fusion'].evidence.mate += 1
                        continue
        
            # Count reads clipped.
            count_clipped = np.zeros((2, params['clipped']['interval']))
            for i in range(params['clipped']['interval']):
                if cigar2pos.get(region.position-i-1,0) in [4, 5]:
                    count_clipped[0][i] = 1
                if cigar2pos.get(region.position+i+1,0) in [4, 5]:
                    count_clipped[1][i] = 1
            
            if np.max(np.sum(count_clipped, axis=1)) >= params['clipped']['count']:
                graph.nodes[n]['fusion'].evidence.clipped += 1

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

def write(filename, name, graph):
    """Write a vcf file from a list of Fusion"""    
    data = {}

    # Header.
    with open(filename, 'w') as fod:
        fod.write(_get_header() + '\n')

    # Fusions.
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [name]
    df = pd.DataFrame(columns=columns)

    nodes = graph.nodes
    for n in nodes:
        logging.debug(graph.nodes[n]['fusion'])
        # First.
        ident = graph.nodes[n]['fusion'].get_name()
        ident_1 = ident
        if graph.nodes[n]['fusion'].second.is_init():
            ident_1 += '-1'
        infos = ['SVTYPE=FUS']
        infos += ['SOFT=%s'%(graph.nodes[n]['fusion'].software,)]
        infos += ['FROM=%s'%('-'.join(graph.label_build_from(n)),)]
        infos += ['CONS=%s'%(graph.nodes[n]['is_consensus'],)]
        infos += ['VAF=%s'%(graph.nodes[n]['fusion'].evidence.get_vaf(),)]
        infos += ['DP=%s'%(graph.nodes[n]['fusion'].evidence.depth,)]
        infos += ['SU=%s'%(graph.nodes[n]['fusion'].evidence.get_sum(),)]
        infos += ['SR=%s'%(graph.nodes[n]['fusion'].evidence.split,)]
        infos += ['PE=%s'%(graph.nodes[n]['fusion'].evidence.mate,)]
        infos += ['SC=%s'%(graph.nodes[n]['fusion'].evidence.clipped,)]
    
        infos = ':'.join(infos)
        values = [graph.nodes[n]['fusion'].first.chrom, graph.nodes[n]['fusion'].first.position]
        values += [ident_1, 'N', '<FUS>', '.', '.', infos]
        values += ['GT:VAF:DP:SU:SR:PE:SC', './.:%s:%s:%s:%s:%s'%(graph.nodes[n]['fusion'].evidence.depth, graph.nodes[n]['fusion'].evidence.get_vaf(), graph.nodes[n]['fusion'].evidence.get_sum(), graph.nodes[n]['fusion'].evidence.split, graph.nodes[n]['fusion'].evidence.mate, graph.nodes[n]['fusion'].evidence.clipped)]
        df = df.append(pd.Series(values, index=columns), ignore_index=True)

        ident_2 = ident
        if graph.nodes[n]['fusion'].second.is_init():
            ident_2 += '-2'
            # Second.
            infos = ';'.join(['SVTYPE=FUS', 'DP=.', 'SU=.'])
            values = [graph.nodes[n]['fusion'].second.chrom, graph.nodes[n]['fusion'].second.position, ident_2, 'N', '<FUS>', '.', '.', infos, 'GT:VAF:DP:SU:SR:PE:SC', './.:.:.']
            df = df.append(pd.Series(values, index=columns), ignore_index=True)
    df.to_csv(filename, mode='a', sep='\t', index=False)
