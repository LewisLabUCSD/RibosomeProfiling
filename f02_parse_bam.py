import pysam


class bam_parse(object):
    """
    This is for single end data analysis.
    Input is a bamHandle read by pysam.AlignmentFile(bamFile,'rb')
    """
    def __init__(self,bamHandle):
        self.bamHandle = bamHandle
        
    def bam_fwd_rev_count(self,max_len=37,seq_len=50):  # get fwd and reverse count data
        """
        Return two dictioanry: {name\tpos:count} for fwd and reverse alignment
        """
        count_dic = {}
        for align in self.bamHandle.fetch():
            if align.alen > seq_len:
                length = align.qlen   # query length: length of the nt in the read that maps
            else:
                length = align.alen   # alignment length
            if length > max_len: continue
            tid = align.reference_id
            ref_name = self.bamHandle.getrname(tid)
            align5 = align.pos + 1
            align3 = align.aend
            if align.is_reverse:
                key = '\t'.join([ref_name,str(align5),str(align3),'-',str(length)])
                if key in count_dic:
                    count_dic[key] = count_dic[key] + 1
                else:
                    count_dic[key] = 1
            else:  # forward strand
                key = '\t'.join([ref_name,str(align5),str(align3),'+',str(length)])
                if key in count_dic:
                    count_dic[key] = count_dic[key] + 1
                else:
                    count_dic[key] = 1
        return count_dic
    
    def align_len_distribution(self,seq_len=50):
        """
        This function returns a dictioanry of alignment length distribution. eg: {29:30%}
        """
        len_dist = {}
        for align in self.bamHandle.fetch():
            if align.alen > seq_len:
                length = align.qlen
            else:
                length = align.alen
            if length in len_dist:
                len_dist[length] = len_dist[length]+1
            else:
                len_dist[length] = 1
        total = self.bamHandle.mapped
        for key in len_dist:
            len_dist[key] = len_dist[key]/float(total)
        return len_dist
        
    