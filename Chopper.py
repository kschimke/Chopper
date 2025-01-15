#!/usr/bin/env python3
# Kayla Schimke
# Christopher Vollmers

import sys
import os
import argparse
import numpy as np
import mappy as mp
racon='racon'
abpoa='abpoa'
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
    parser.add_argument('-i','--iteration',type=int,default=1, help='How many polished plasmid sequences to generate')
    parser.add_argument('-V','--verbose',action='store_true', help='Print log text to console.')
    parser.add_argument('-v','--version',action='version', version='Chopper.py '+VERSION, help='Prints the version and exits.')
    args = parser.parse_args()
    return args

def q_score_to_value(q_score_char):
    # Convert Q score character to its numerical value
    return ord(q_score_char) - 33

def average_quality(q_scores):
    numerical_values=[]
    # Convert each Q score character to its numerical value
    for i in range(0,len(q_scores),50):
        char=q_scores[i]
        value=ord(char) - 33
        numerical_values.append(value)
    # Calculate the average of these numerical values
    return sum(numerical_values) / len(numerical_values)


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

def assembleReads(output_path,reads,subreads,subsample,iteration):
    medaka_path = output_path+'medaka/'
    if os.path.exists(medaka_path):
        os.system(f'rm -r {medaka_path}')
    os.mkdir(medaka_path)

    if reads == subreads:
        peak_fastx1 = medaka_path+'peak_reads1.fasta'
        peak_fastx2 = medaka_path+'peak_reads2.fasta'
    else:
        peak_fastx1 = medaka_path+'peak_reads1.fastq'
        peak_fastx2 = medaka_path+'peak_reads2.fastq'

    print('reads',len(reads),'subreads',len(subreads))

    if subsample:
        qSubs=[]
        if peak_fastx1 == medaka_path+'peak_reads1.fastq':
            for name1 in subreads:
                seq,q1=subreads[name1]
                aveQual=average_quality(q1)
                qSubs.append((aveQual,name1))
            sortedqSubs=sorted(qSubs,reverse=True)
        allSamples=[]
        for q1,name1 in sortedqSubs:
            allSamples.append(name1)
        sub_sample1=allSamples[:subsample]
        sub_sample2=allSamples[:subsample]
#        sub_sample1 = np.random.choice(list(subreads.keys()),size=min(len(subreads),subsample),replace=False)
#        sub_sample2 = np.random.choice(list(subreads.keys()),size=min(len(subreads),subsample),replace=False)
        writeFastX(peak_fastx1,subreads,sub_sample1)
        writeFastX(peak_fastx2,subreads,sub_sample2)
    else:
        writeFastX(peak_fastx1,subreads,list(subreads.keys()))
        writeFastX(peak_fastx2,subreads,list(subreads.keys()))


    # create reference fasta
    max_qscore = 0
    max_repeats = 0
    r2c2_reads=[]
    for name,sequence1 in reads.items():
        repeats=int(name.split('_')[3])
        qscore=float(name.split('_')[1])
        if repeats<=5:
            r2c2_reads.append((repeats,qscore,sequence1))
#        if repeats>max_repeats:
#             best=(name,sequence1)
#             max_repeats=repeats
#             max_qscore=qscore
#        elif repeats==max_repeats:
#             if qscore>max_qscore:
#                 best=(name,sequence1)
#                 max_qscore=qscore
    sorted_r2c2_reads=sorted(r2c2_reads,reverse=True)
#    print(sorted_r2c2_reads[iteration][:2])
    sample=np.random.choice(list(reads.keys()),size=min(len(reads),1),replace=False)
    best=(iteration,sorted_r2c2_reads[iteration][2])
#    best=(sample[0],reads[sample[0]])
    plasmidLength=len(best[1])

#    poa_reads = output_path+'/poa_reads.fasta'
#    poa_output = output_path+'/poa_output.fasta'
    poa_cons = output_path+'/poa_consensus.fasta'
#    oriented_reads={}
#    mm_align = mp.Aligner(seq=best[1], preset='map-ont')
#    for name,sequence in reads.items():
#        repeats=int(name.split('_')[3])
#        qscore=float(name.split('_')[1])
#        for hit in mm_align.map(sequence):
#            if hit.strand==1:
#                oriented_reads[name]=sequence
#            elif hit.strand==-1:
#                oriented_reads[name]=mp.revcomp(sequence)

#    if not oriented_reads:
#        print('no sequences','reads')
#        consensus_sequence = best[1]+best[1]
#    else:
#        poa_fh=open(poa_reads,'w')
#        poaCounter=0
#     sample = np.random.choice(list(reads.keys()),size=min(len(reads),50),replace=False)
#     writeFastX(poa_reads,reads,sample)
#        sortedSequences=sorted(sequences,reverse=True)
#        for qscore,sequence in sortedSequences[:50]:
#            poaCounter+=1
#            poa_fh.write(f'>{poaCounter}\n{sequence}\n')
#        poa_fh.close()

#        os.system(f'{abpoa} -r 2 -S {poa_reads} > {poa_output} 2> apboa.messages')
#        cons_seq=''
#        for consName,consSeq,consQ in mp.fastx_read(f'{poa_output}'):
#            cons_seq=consSeq
#        if cons_seq:
#            consensus_sequence = cons_seq

#    if len(consensus_sequence)>plasmidLength*1.5:
#        center = int(len(consensus_sequence)/2)
#        start = center-int(plasmidLength*0.75)
#        end=start+int(plasmidLength*1.5)
#        consensus_sequence=consensus_sequence[start:end]
    consensus_sequence=best[1]+best[1][:int(len(best[1])/2)]
    windowLength=len(consensus_sequence)
#    print(len(consensus_sequence),plasmidLength)

    out_cons_file = open(poa_cons, 'w')
    out_cons_file.write(f'>Consensus\n{consensus_sequence}\n')
    out_cons_file.close()



    duplicate_read=consensus_sequence
    # print(f'Reference Read Length: {len(duplicate_read)}')
    overlap=output_path+'overlap.paf'
    overlap_fh=open(overlap,'w')

    output_cons=output_path+'racon.fasta'

    mm_align = mp.Aligner(seq=duplicate_read, preset='map-ont')
    for name1 in sub_sample1:
        sequence1,q1 = subreads[name1]
        for hit in mm_align.map(sequence1):
            overlap_fh.write(f'{name1}\t{str(len(sequence1))}\t{hit.q_st}\t{hit.q_en}\t{hit.strand}\tConsensus\t{hit.ctg_len}\t{hit.r_st}\t{hit.r_en}\t{hit.mlen}\t{hit.blen}\t{hit.mapq}\n')

    overlap_fh.close()
    os.system(f'{racon} --no-trimming -q 20 -t 1 --window-length {windowLength} {peak_fastx1} {overlap} {poa_cons} > {output_cons} 2>./racon_messages.txt')

    raconreads=parseFastX(output_cons)

#    duplicate_read=raconreads['Consensus']

    ref_fasta= medaka_path+'refread.fasta'
    reference = open(ref_fasta, 'w')

    reference.write(f'>Duplicated_High_Coverage_Consensus_Read\n{duplicate_read}\n')
    reference.close()

    medaka_messages = medaka_path+'medaka_messages.txt'
    medaka_errors = medaka_path+'medaka_errors.txt'
    minimapSAM = medaka_path+'alignment.sam'
    minimapBAM = medaka_path+'alignment.bam'
    minimapSortedBAM = medaka_path+'alignment_sorted.bam'
    medakaHDF = medaka_path+'medaka.hdf'
    medaka_fasta = medaka_path+'medaka.fasta'
    os.system(f'minimap2 -ax map-ont --secondary=no -t 1 {ref_fasta} {peak_fastx2} >{minimapSAM} 2> {medaka_path}minimap2.messages')
    os.system(f'samtools view -b {minimapSAM} >{minimapBAM}')
    os.system(f'samtools sort {minimapBAM} >{minimapSortedBAM}')
    os.system(f'samtools index {minimapSortedBAM}')
#    os.system(f'medaka consensus {minimapSortedBAM} {medakaHDF} --model r941_min_sup_g507 > {medaka_messages} 2> {medaka_errors}')
    os.system(f'CUDA_VISIBLE_DEVICES="" medaka inference {minimapSortedBAM} {medakaHDF} --model r1041_e82_400bps_sup_v5.0.0 > {medaka_messages} 2> {medaka_errors}')
    os.system(f'CUDA_VISIBLE_DEVICES="" medaka sequence {medakaHDF} {ref_fasta} {medaka_fasta} >> {medaka_messages} 2>> {medaka_errors}')
#    os.system(f'medaka_consensus -i {peak_fastx} -d {ref_fasta} -o {medaka_path} >> {medaka_messages} 2>> {medaka_errors}')
    corrected_consensus = ''
    for name,seq,qual,comment in mp.fastx_read(medaka_fasta,read_comment=True):
#        center = int(len(seq)/2)
#        start = center-int(plasmidLength/2)
        start=10
        screenrange = (start+plasmidLength-100,start+plasmidLength+100)
#        print(len(seq),start,screenrange)
        for i in range(screenrange[0],screenrange[1],1):
#            print(seq[start:start+20],seq[i:i+20])
            if seq[start:start+20] == seq[i:i+20]:
                corrected_consensus = seq[start:i]
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
                for iteration in range(args.iteration):
                    assembly = assembleReads(output_path,inRange,inRangeSub,args.subsample,iteration)
                    if args.verbose: sys.stderr.write(f'Assembly Length: {len(assembly)}\n')
                    af.write(f'>{type1}_{minimum}-{maximum}_assembly_{iteration}\n{assembly}\n')
        af.close()

main()
