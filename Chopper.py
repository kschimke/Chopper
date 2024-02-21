#!/usr/bin/env python3
# Kayla Schimke
# Christopher Vollmers

import sys
import os
import argparse
import numpy as np
import mappy as mp
VERSION = '1.0.0'

def commandline():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('mode',nargs='?',default='detect,assemble', help='Comma separated list of modes to run. detect: Bins consensus reads and detects plasmid length based on peak bins; assemble: Assembles reads around plasmid length(s); plot: runs detect and reports the bins and peaks within a range. Default: detect,assemble.')
    parser.add_argument('-c','--consensus',type=str, required=True, help='FASTA file containing basecalled consensus reads')
    parser.add_argument('-s','--subreads',type=str, help='FASTQ file containing basecalled R2C2 subreads')
    parser.add_argument('-o','--output',type=str,default=os.getcwd(), help='Directory to create and store output files. Defaults to current directory.')
    parser.add_argument('-b','--binsize',type=int,default=100, help='Bin size when sorting reads by length for detect, also used to define a range around --plasmidsize. Defaults to 100.')
    parser.add_argument('-p','--plasmidsize',type=str, help='Individual or comma separated list of sizes corresponding with plasmids in data. If running detect, set in addition to automatically detected peaks.')
    parser.add_argument('-r','--range',type=str,default='0,20000', help='The inclusive range of bins to print to the console. Defaults to 0,20000.')
    parser.add_argument('-t','--threshold',type=int,default=5, help='Minimum percentage of total consensus reads required to call a peak. Defaults to 5%%.') # a single % causes issue with parsing help text
    parser.add_argument('-su','--subsample',type=int, help='Max number of subreads used to make an assembly.')
    parser.add_argument('-V','--verbose',action='store_true', help='Print log text to console.')
    parser.add_argument('-v','--version',action='version', version='Chopper.py '+VERSION, help='Prints the version and exits.')
    args = parser.parse_args()
    return args

def parseFastX(fastx,qual=False):
    reads = {}
    for name,seq,q in mp.fastx_read(fastx):
        if qual:
            reads[name] = [seq,q]
        else:
            reads[name] =  seq

    return reads

def writeFastX(file_name,reads,sample):
    file = open(file_name,'w')
    ext = file_name.split('.')[-1]
    if ext == 'fasta':
        for name in sample:
            seq = reads[name]
            file.write(f'>{name}\n{seq}\n')
    elif ext == 'fastq':
        for name in sample:
            seq,q = reads[name]
            file.write(f'@{name}\n{seq}\n+\n{q}\n')
    file.close()

def rounddown(x,size):
    return int(np.floor(x / size)) * size

def binReads(reads,binsize):
    bins = {}
    for name in reads:
        bin = rounddown(len(reads[name]),binsize)
        if bin not in bins:
            bins[bin]={}
        bins[bin][name] = reads[name]

    return bins

def findPeaks(bins,binsize,totalReads,threshold):
    prev = 0
    peaks = []
    for i in range(min(bins)+binsize,max(bins),binsize):
        if i not in bins:
            bins[i] = []
        current = len(bins[i-binsize])
        next = len(bins[i])
        if prev < current > next and (current/totalReads)*100 >= threshold:
            peaks.append((i-binsize,i,'automatically_detected',bins[i-binsize]))
        prev = current

    return peaks

def printPeaks(binsize,printrange,bins,peaks):
    minimum,maximum=[int(i) for i in printrange.split(',')]
    for step in range(minimum,maximum+binsize,binsize):
        if step in bins:
            string=str(step)+'\t'+str(len(bins[step]))
        else:
            string=str(step)+'\t0'
        if step in peaks:
            string += ' <'+'-'*(20-len(string))+'peak'
        print(string)

def assembleReads(output_path,reads,subreads,subsample):
    medaka_path = output_path+'medaka/'
    if os.path.exists(medaka_path):
        os.system(f'rm -r {medaka_path}')
    os.mkdir(medaka_path)

    # create reference fasta
    max_coverage = 0
    for name in reads:
        try:
            coverage=float(name.split('_')[3])
        except:
            coverage=float(name.split('_')[2])

        if coverage<=5: # why are we limiting the coverage to 5?
            if coverage>=max_coverage:
                 best=(name,reads[name])
                 max_coverage=coverage

    duplicate_read=best[1]+best[1]
    # print(f'Reference Read Length: {len(duplicate_read)}')

    ref_fasta= medaka_path+'refread.fasta'
    reference = open(ref_fasta, 'w')
    reference.write(f'>Duplicated_High_Coverage_Consensus_Read\n{duplicate_read}\n')
    reference.close()

    if reads == subreads:
        peak_fastx = medaka_path+'peak_reads.fasta'
    else:
        peak_fastx = medaka_path+'peak_reads.fastq'

    if subsample:
        sample = np.random.choice(list(subreads.keys()),size=min(len(subreads),subsample),replace=False)
        writeFastX(peak_fastx,subreads,sample)
    else:
        writeFastX(peak_fastx,subreads,list(subreads.keys()))

    medaka_messages = medaka_path+'medaka_messages.txt'
    medaka_errors = medaka_path+'medaka_errors.txt'
    minimapSAM = medaka_path+'alignment.sam'
    minimapBAM = medaka_path+'alignment.bam'
    minimapSortedBAM = medaka_path+'alignment_sorted.bam'
    medakaHDF = medaka_path+'medaka.hdf'
    medaka_fasta = medaka_path+'medaka.fasta'
    os.system(f'minimap2 -ax map-ont --secondary=no -t 1 {ref_fasta} {peak_fastx} >{minimapSAM} 2> {medaka_path}minimap2.messages')
    os.system(f'samtools view -b {minimapSAM} >{minimapBAM}')
    os.system(f'samtools sort {minimapBAM} >{minimapSortedBAM}')
    os.system(f'samtools index {minimapSortedBAM}')
#    os.system(f'medaka consensus {minimapSortedBAM} {medakaHDF} --model r941_min_sup_g507 > {medaka_messages} 2> {medaka_errors}')
    os.system(f'/home/k2so/.local/bin/medaka consensus {minimapSortedBAM} {medakaHDF} --model r1041_e82_400bps_sup_v4.0.0 > {medaka_messages} 2> {medaka_errors}')
    os.system(f'/home/k2so/.local/bin/medaka stitch {medakaHDF} {ref_fasta} {medaka_fasta} >> {medaka_messages} 2>> {medaka_errors}')

    corrected_consensus = ''
    for name,seq,qual,comment in mp.fastx_read(medaka_fasta,read_comment=True):
        for i in range(int(len(seq)/2)-50,len(seq),1):
            if seq[:20] == seq[i:i+20]:
                corrected_consensus = seq[:i]
                break

    return corrected_consensus

def writeLog(output_path,args):
    log = open(output_path+'/Chopper.log','w')
    log.write(f'Chopper {VERSION}\n\n')
    log.write(f'Mode(s): {args.mode}\n')
    log.write(f'Output Directory: {output_path}\n')
    log.write(f'Consensus Reads: {args.consensus}\n')
    if args.subreads:
        log.write(f'Subreads: {args.subreads}\n')
    else:
        log.write(f'Subreads: not provided\n')
    log.write(f'Bin Size: {args.binsize}\n')
    if 'plot' in args.mode:
        log.write(f'Plot Range: {args.range}\n')
    if 'detect' in args.mode:
        log.write(f'Peak Threshold: {args.threshold}\n')
    if 'assemble' in args.mode:
        log.write(f'Plasmid Size: {args.plasmidsize}\n')
        log.write(f'Subsample: {args.subsample}\n')
    log.close()

def main():
    args = commandline()

    output_path = args.output
    if not output_path.endswith('/'):
        output_path += '/'
    if not os.path.exists(output_path):
        if args.verbose: sys.stderr.write(f'Making output directory: {output_path}\n')
        os.mkdir(output_path)

    writeLog(output_path, args)

    consensus_fasta = args.consensus
    if not os.path.exists(consensus_fasta):
        sys.stderr.write('Consensus FASTA does not exist.\n')
        sys.exit(1)
    if args.verbose: sys.stderr.write(f'Reading consensus FASTA: {consensus_fasta}\n')
    consensus_reads = parseFastX(consensus_fasta)

    detected=[]
    if 'detect' in args.mode or 'plot' in args.mode:
        if args.verbose: sys.stderr.write(f'Detecting plasmid lengths from consensus reads.\n')
        bins = binReads(consensus_reads,args.binsize)
        detected = findPeaks(bins,args.binsize,len(consensus_reads),args.threshold)

        if 'plot' in args.mode:
            printPeaks(args.binsize,args.range,bins,[i[0] for i in detected])

    if 'assemble' in args.mode:
        peaks=detected
        if not peaks and not args.plasmidsize:
            sys.stderr.write(f'Must specify --plasmidsize when running assemble without detect.\n')
            sys.exit(1)

        if args.plasmidsize:
            plasmidsizes=[int(i) for i in args.plasmidsize.split(',')]
            for p in plasmidsizes:
                minimum = p - args.binsize/2
                maximum = p + args.binsize/2
                manual_bin = {}
                for name in consensus_reads:
                    length=len(consensus_reads[name])
                    if minimum<=length<=maximum:
                        manual_bin[name] = consensus_reads[name]
                peaks.append((minimum,maximum,'manually_entered',manual_bin))

        asm_fasta = output_path+'R2C2_full_length_plasmids.fasta'
        af = open(asm_fasta,'w')

        if not args.subreads:
            if args.verbose: sys.stderr.write('No subreads provided. Using only consensus reads to assemble.\n')
            for minimum,maximum,type1,inRange in peaks:
                if args.verbose: sys.stderr.write(f'Assembling reads from {type1} range {minimum}-{maximum}\n')
                assembly = assembleReads(output_path,inRange,inRange,args.subsample) #inRange doesn't have qual
                if args.verbose: sys.stderr.write(f'Assembly Length: {len(assembly)}\n')
                af.write(f'>{type1}_{minimum}-{maximum}_assembly\n{assembly}\n')

        elif not os.path.exists(args.subreads):
            sys.stderr.write('Subreads file does not exist.\n')
            sys.exit(1)

        else:
            if args.verbose: sys.stderr.write(f'Reading subreads FASTQ: {args.subreads}\n')
            subreads = parseFastX(args.subreads, True)

            for minimum,maximum,type1,inRange in peaks:
                if args.verbose: sys.stderr.write(f'Assembling reads from {type1} range {minimum}-{maximum}\n')
                names = set([i.split('_')[0] for i in inRange])
                inRangeSub={}
                for name in subreads:
                    if name.split('_')[0] in names:
                        length=len(subreads[name][0])
                        if minimum<=length<=maximum:
                            inRangeSub[name] = subreads[name]

                assembly = assembleReads(output_path,inRange,inRangeSub,args.subsample)
                if args.verbose: sys.stderr.write(f'Assembly Length: {len(assembly)}\n')
                af.write(f'>{type1}_{minimum}-{maximum}_assembly\n{assembly}\n')
        af.close()

main()
