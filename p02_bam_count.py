'''
count reads at each position of the genome including strand information. 
The format of the dataframe is: [count of reads, chromosome, start, end,strand,length of mapped protion of the read]
'''
import os,argparse
from natsort import natsorted
from f02_parse_bam import bam_parse
import pysam,glob
from multiprocessing import Pool


parser = argparse.ArgumentParser(description='Prepare read count at each position, each chromosome, in both strands, each mapping length')
parser.add_argument('-b','--bam',action='store',dest='bam_path',help='path to bam files')
parser.add_argument('-t','--thread',action='store',dest='thread',type=int,help='number of thread',default=1)
parser.add_argument('-m','--maxl',action='store',dest='max_len',type=int,help='maximum mapping length to consider',default=37)
parser.add_argument('-r','--readl',action='store',dest='read_len',type=int,help='read length in sequencing',default=50)

args = parser.parse_args()

bam_path = args.bam_path
if bam_path.endswith('/'): bam_path = bam_path[:-1]
thread   = args.thread
max_len  = args.max_len
read_len = args.read_len

def fwd_rev_cov(bamFile,outpath,max_len=max_len,seq_len=75):
    """
    This function calculates mapping 5' and 3' end position of each read and count the number for those with the same positions
    * bamFile: str. Bam file
    output a file with cover dataframe. The columns are: ['num','chr','end5','end3','strand','len'].
    num is number of reads. len is length of mapping in reads
    """
    if not os.path.exists(cov_path):
        os.mkdir(cov_path)
    
    bamHandle = pysam.AlignmentFile(bamFile,'rb')
    Handle = bam_parse(bamHandle)
    count_dic = Handle.bam_fwd_rev_count(max_len,seq_len)
    covFile = outpath+'/'+ os.path.basename(bamFile)[:3] + '_cov.txt'
    
    def write_dic(dic,File):
        """
        This function write the dictionary to a file with 'key tab value'
        * dic: dict. Dictionary that need to be write.
        * File: filename. File that stores the output.
        """
        Handle = open(File,'w')
        for key in natsorted(dic):
            Handle.write(str(dic[key])+'\t'+key+'\n')
        Handle.close()
    
    write_dic(count_dic,covFile)




thread = 6
# run fwd_rev_cov in parallel
bamFiles = natsorted(glob.glob(bam_path+'/*.sort.bam'))
cov_path = os.path.dirname(bam_path) + '/02_cov'
if not os.path.exists(cov_path): os.mkdir(cov_path)

p = Pool(processes=thread)
for bam in bamFiles:
    p.apply_async(fwd_rev_cov,args=(bam,cov_path,max_len,read_len,))
    p.close()
    p.join()



