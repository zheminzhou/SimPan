# SimPan
A simulator of pan-genome for bacterial population

```$ python SimPan.py -h
usage: SimPan.py [-h] [-p PREFIX] [--genomeNum GENOMENUM] [--geneLen GENELEN]
                 [--igrLen IGRLEN] [--backboneBlock BACKBONEBLOCK]
                 [--mobileBlock MOBILEBLOCK] [--operonBlock OPERONBLOCK]
                 [--aveSize AVESIZE] [--nBackbone NBACKBONE] [--nCore NCORE]
                 [--nMobile NMOBILE] [--pBackbone PBACKBONE]
                 [--pMobile PMOBILE] [--tipAccelerate TIPACCELERATE]
                 [--rec REC] [--recLen RECLEN] [--seqRec SEQREC]
                 [--insRec INSREC] [--delRec DELREC] [--noSeq]
                 [--idenOrtholog IDENORTHOLOG] [--idenParalog IDENPARALOG]
                 [--idenDuplication IDENDUPLICATION] [--indelRate INDELRATE]
                 [--indelMax INDELMAX] [--freqStart FREQSTART]
                 [--freqStop FREQSTOP]

SimPan is a simulator for bacterial pan-genome.
Global phylogeny and tree distortions are derived from SimBac and the gene and intergenic sequences are simulated using indel-seq-gen.

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        prefix for all intermediate files and outputs
  --genomeNum GENOMENUM
                        No of genome in population [DEFAULT: 20]
  --geneLen GENELEN     [negative bionomial with r=2] mean,min,max sizes of genes [DEFAULT: 900,150,6000]
  --igrLen IGRLEN       [negative bionomial] mean,min,max sizes of intergenic regions [DEFAULT: 50,0,300]
  --backboneBlock BACKBONEBLOCK
                        [geometric] mean,min,max number of backbone genes per block [DEFAULT: 3,0,30]
  --mobileBlock MOBILEBLOCK
                        [geometric] mean,min,max number of mobile genes per block [DEFAULT: 10,0,100]
  --operonBlock OPERONBLOCK
                        [geometric] mean,min,max number of continuous genes that share the same coding strand [DEFAULT: 3,0,15]
  --aveSize AVESIZE     average gene number per genome (greater than nBackbone). [DEFAULT: 4500]
  --nBackbone NBACKBONE
                        number of backbone genes (present in common ancestor) per genome. [DEFAULT: 4000]
  --nCore NCORE         sizea of core gene (smaller than the size of backbone genes). [DEFAULT: 3500]
  --nMobile NMOBILE     size of mobile gene pool for accessory genome. [DEFAULT: 20000]
  --pBackbone PBACKBONE
                        propotion of paralogs in backbone (core) genes. [DEFAULT: 0.05]
  --pMobile PMOBILE     propotion of paralogs in mobile (accessory) genes. [DEFAULT: 0.4]
  --tipAccelerate TIPACCELERATE
                        grandient increasing of gene indels in recent times. [DEFAULT: 100]
  --rec REC             expected coverage of homoplastic events in pairwise comparisons. [DEFAULT: 0.05]
  --recLen RECLEN       expected size of homoplastic events. [DEFAULT: 1000]
  --seqRec SEQREC       Use homoplastic events to infer sequences. Use 0 to disable [DEFAULT: 1]
  --insRec INSREC       Use homoplastic events to infer gene insertions. Use 0 to disable [DEFAULT: 1]
  --delRec DELREC       Use homoplastic events to infer gene deletions. Use 0 to disable [DEFAULT: 1]
  --noSeq               Do not infer sequence but only the gene presence/absence. [DEFAULT: False]
  --idenOrtholog IDENORTHOLOG
                        average nucleotide identities for orthologous genes. [DEFAULT: 0.98]
  --idenParalog IDENPARALOG
                        average nucleotide identities for paralogous genes. [DEFAULT: 0.6]
  --idenDuplication IDENDUPLICATION
                        average nucleotide identities for recent gene duplications. [DEFAULT: 0.995]
  --indelRate INDELRATE
                        summarised average frequency of indel events. [DEFAULT: 0.01]
  --indelMax INDELMAX   maximum size of short indel events within each gene (<=300). [DEFAULT: 30]
  --freqStart FREQSTART
                        frequencies of start codons of ATG,GTG,TTG. DEFAULT: 0.83,0.14,0.03
  --freqStop FREQSTOP   frequencies of stop codons of TAA,TAG,TGA. DEFAULT: 0.63,0.08,0.29
```
