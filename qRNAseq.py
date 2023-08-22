#!/usr/bin/env python3

####author: CY G
####nongxinshengxin

import re
import os
import argparse
from textwrap import dedent
import sys


class Hats():
    
    def __init__(self,directory,genome,annotation,output,read,strand,intro,thread):
        self.directory=directory
        self.genome=genome
        self.annotation=annotation
        self.output=output
        self.read=read
        self.strand=strand
        self.intro=intro
        self.thread=thread

    
    def find_fastq(self,x):
        fqlist=os.listdir(x)
        fqs={}
        for fqfile in fqlist:
            a=re.search(r'(.+)1.(clean.fq.gz|clean.fastq.gz|clean.fq|clean.fastq|fq.gz|fastq.gz|fq|fastq)',fqfile)
            if a!=None:
                part1=str(a.group(1))
                part2=str(a.group(2))
                fastq1=part1+'1.'+part2
                fastq2=part1+'2.'+part2
                fqs[part1]=[fastq1,fastq2]
                sys.stderr.write('+%s\n' %part1)
                sys.stderr.write('- fastq file 1:%s\n'%fastq1)
                sys.stderr.write('- fastq file 1:%s\n'%fastq2)
        return fqs

                


    def rnaseq(self):
        sys.stderr.write('Parsing fastq file(s) in %s \n' %self.directory)
        fqs=self.find_fastq(self.directory)
        num=len(fqs.keys())
        if num < 1:
            sys.stderr.write('ERROR:Fastq files could not be found!\n')
            sys.exit()
        sys.stderr.write('Creating new output-dir:%s\n'%self.output)
        outdir=self.make_dir(self.output)

        jug1=os.path.exists(self.genome+'.1.ht2')
        if jug1 is not True:
            sys.stderr.write("No index in this directory, build a new index.")
            self.build_db('hisat2-build',self.genome,outdir)
        sys.stderr.write('Make RNAseq flow\n')
        rnaseq_strand_op=''
        if self.strand=='1':
            rnaseq_strand_op='--rna-strandness RF'
        if self.read=="p":
            for key,value in fqs.items():
                sys.stderr.write('Mapping %s'%key)
                hisat2_option='%s --max-intronlen %s --no-unal --no-temp-splicesite --novel-splicesite-outfile %s/bam/%s.splicesite --novel-splicesite-infile %s/bam/%s.splicesite'%(rnaseq_strand_op,self.intro,outdir,key,outdir,key)
                mapping_cmd='hisat2 -p %s %s -x %s -1 %s/%s -2 %s/%s'%(self.thread,hisat2_option,self.genome,self.directory,value[0],self.directory,value[1])
                sam2bam_cmd='samtools view -bS -'
                sortbam_cmd='samtools sort -o %s/bam/%s.bam -'%(outdir,key)
                os.system('%s|%s|%s'%(mapping_cmd,sam2bam_cmd,sortbam_cmd))
        elif self.read=="s":
            for fq in os.listdir(self.directory):

                sys.stderr.write('Mapping %s'%key)
                hisat2_option='%s --max-intronlen %s --no-unal --no-temp-splicesite --novel-splicesite-outfile %s/bam/%s.splicesite --novel-splicesite-infile %s/bam/%s.splicesite'%(rnaseq_strand_op,self.intro,outdir,key,outdir,key)
                mapping_cmd='hisat2 -p %s %s -x %s -U %s/%s '%(self.thread,hisat2_option,self.genome,self.directory,fq)
                sam2bam_cmd='samtools view -bS -'
                sortbam_cmd='samtools sort -o %s/bam/%s.bam -'%(outdir,key)
                os.system('%s|%s|%s'%(mapping_cmd,sam2bam_cmd,sortbam_cmd))
        else:
            sys.stderr.write('ERROR:Parameter error!\n')
        sys.stderr.write('Start feature counting...\n')
        featureCounts_cmd='featureCounts -p -s 2 -t exon -g gene_id -a %s -o %s/featureCounts/rnaseq.raw.tsv %s/bam/*.bam'%(self.annotation,outdir,outdir)
        os.system(featureCounts_cmd)
        sys.stderr.write('RNAseq has been done!\n')
        self.format_fctable(outdir)
        sys.stderr.write('DONE!\n')
            

    def build_db(self,app,file,dir):
        os.system('cp %s %s/refseq.fa'%(file,dir))
        os.system('%s %s %s/refseq.fa'%(app,file,dir))
        self.genome='%s/refseq.fa'%(dir)

    
    def make_dir(self,outdir):
        os.system('mkdir -p %s/bam'%outdir)
        os.system('mkdir -p %s/featureCounts'%outdir)
        return outdir
        
    def format_fctable(self,dir):
        raw=open('%s/featureCounts/rnaseq.raw.tsv'%dir,'r')
        tsv=open('%s/featureCounts/rnaseq.tsv'%dir,'w')
        for line in raw.readlines():
            if len(line.strip().split('\t'))==1:
                tsv.write('')
            else:
                lin=line.strip().split('\t',6)
                print(lin[0]+'\t'+lin[6],file=tsv)
        raw.close()
        tsv.close()


        
    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-d', '--directory',type=str,help='directory that containing FASTQ file(s), .fq/.fastq(.gz) fq文件所在文件夹')
    parser.add_argument('-g', '--genome',type=str,help='reference genome 参考基因组文件或参考基因组索引文件的前缀')
    parser.add_argument('-a', '--annotation',type=str,help='annotation file in GTF format 注释GTF路径')
    parser.add_argument('-o', '--output',type=str,help='output file path 结果输出文件存放路径')
    parser.add_argument('-r', '--read',type=str,default='p',help=dedent('\
            Single-read or Paired-End 单端测序或双端测序 \
                -r p Paired-End [default] 双端测序(默认)\
                -r s Single-read 单端测序'))
    parser.add_argument('-s', '--strand',type=str,default='1',help=dedent('\
            strand-specific information 链特异性信息 \
                -s 1 stranded [default] 链特异性(默认)\
                -s 0 unstranded 非链特异性'))
    parser.add_argument('-i','--intro',type=str,default='4000',help='max intron length, default 4000 bp 最大内含子长度，默认4000bp')
    parser.add_argument('-t','--thread',type=str,default='1',help='Running threads, default 1 运行线程，默认1')
    args=parser.parse_args()
    hat=Hats(args.directory,args.genome,args.annotation,args.output,args.read,args.strand,args.intro,args.thread)
    hat.rnaseq()







if __name__=='__main__':
    main()
