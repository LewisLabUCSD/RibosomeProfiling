class trpr(object):
    """
    Input is a pandas dataframe from 01_pr_cds.txt. Chr    cds_start    cds_end    GeneID    PrAccess    Strand
    """
    def __init__(self,df):
        self.df = df
        self.df.columns = ['chr','start','end','geneid','access','strand']
        self.df['geneid'] = self.df['geneid'].astype(str)
    

    def get_longest_trprs(self):
        '''this function extracts longest rna or protein'''
        df = self.df
        df['len'] = df['end'] - df['start']

        longest = df[['geneid','access','len']].groupby(['geneid','access']).sum()
        longest['gid'] = longest.index.get_level_values(0)
        longest['feid']   = longest.index.get_level_values(1)
        longest = longest.groupby(level='geneid').max()
        longest['gid'] = longest.index
        return longest

    def get_trpr_pos(self,fid):
        df = self.df
        pos_df = df[df['access'].values==fid]
        pos = []
        strd = set(pos_df['strand'])
        if len(strd) > 1:
            assert False,'map to different strand'
        for s,e in zip(pos_df['start'],pos_df['end']):
            pos.extend(range(s,e))
        if strd == {'-'}:
            pos = pos[::-1]
        return pos
    # id part
#     def get_chrom(self,g,id_type='gene'):
#         """
#         This function get chromosome from gene id
#         """
#         if id_type == 'gene':
#             q_id = 'geneid'
#         else:
#             q_id = 'access'
#         g_df = self.df[self.df[q_id].values==g]
#         chrom = list(set(g_df['chr'].tolist()))
#         if len(chrom)==1:
#             return chrom[0]
#         else:
#             assert False,chrom + 'is wrong,more than one mapping'
    
#     def all_gene_ids(self):
#         """get all the gene ids in the dataframe"""
#         genes = list(set(self.df['geneid'].tolist()))
#         return genes
    
#     def fwd_gene_ids(self):
#         """get all gene ids in the forward strand"""
#         genes = list(set(self.df[self.df['strand'].values=='+'].tolist()))
#         return genes
    
#     def rev_gene_ids(self):
#         """get all gene ids in the reverse strand"""
#         genes = list(set(self.df[self.df['strand'].values=='-'].tolist()))
#         return genes
    
#     def genes2multi_chr(self):
#         """This function finds genes that map to multiple chromosomes"""
#         gene_chr_df = self.df[['chr','geneid']].drop_duplicates()
#         gene_chr_dic = gene_chr_df.groupby('geneid').groups
#         genes = []
#         for key in gene_chr_dic:
#             if len(gene_chr_dic[key])>1:
#                 genes.append(key)
#         return genes
        
#     def trpr_gene_dic(self):
#         """This function generates a function of {praccess:geneid}"""
#         pr_g_dic = self.df.set_index('access')['geneid'].to_dict()
#         return pr_g_dic


#     ##### position part
#     def get_longest_trpr(self,gene):
#         """get the longest protein access or transcript access"""
#         gene_df = self.df[self.df['geneid'].values==gene].copy()
#         gene_df.loc[:,'len'] = gene_df['end'] - gene_df['start']
#         trprs = list(set(gene_df['access'].tolist()))
#         if len(trprs)!=1:
#             lens = {}
#             for trpr in trprs:
#                 trpr_len = gene_df[gene_df['access'].values==trpr]['len'].sum()
#                 lens[trpr_len] = trpr
#             return lens[max(lens.keys())]
#         else:
#             return trprs[0]    
    
#     def get_gene_trpr_len(self,gene,multi_chrom='N'):
#         """
#         This function gets the lenth of all exon/CDS in a gene.
#         * gene: str. gene id
#         * multi_chrom: N or Y. Whether to include the length of multiple chromosome
#         """
#         trpr_df = self.df[self.df['geneid'].values==gene]
#         chrom = list(set(trpr_df['chr'].tolist()))
#         if multi_chrom == 'N':
#             if len(chrom) > 1:
#                 assert False,'gene map to multiple chromosome'
#         # get total length
#         pos = []
#         for start,end in zip(trpr_df['start'],trpr_df['end']):
#             inter = range(int(start)+1,int(end)+1)
#             pos = inter+pos
#         return len(list(set(pos)))
            
    
#     def get_trpr_pos(self,trpr,level='trpr'):
#         """
#         This function gets all position of transcript or protein or gene
#         """
#         # decide whether to proceed
#         trpr = str(trpr)
#         if level == 'trpr':
#             trpr_df = self.df[self.df['access'].values==trpr]
#         elif level == 'gene':
#             trpr_df = self.df[self.df['geneid'].values==trpr]
#         chrom = list(set(trpr_df['chr'].tolist()))
#         if len(chrom) !=1:
#             assert False,'Protein {trpr} map to multiple chromosome'.format(trpr=trpr)
#         strand = list(set(trpr_df['strand'].tolist()))
#         if len(strand)!=1: 
#             assert False,'Protein {trpr} map to both strand'.format(trpr=trpr)
#         # get position
#         pos = [] # store the cds positions
#         for start,end,stra in zip(trpr_df['start'],trpr_df['end'],trpr_df['strand']):
#             inter = range(int(start)+1,int(end)+1)
#             if stra == '-':
#                 inter.reverse()
#                 pos = inter + pos
#             else:
#                 pos.extend(inter)
#         if pos[0]>pos[1]: # - strand
#             new_pos = natsorted(list(set(pos)))[::-1]
#         else:
#             new_pos = natsorted(list(set(pos)))
#         return new_pos
    
#     def get_ribo_pr_pos(self,trpr,rna_pos,offset,level='trpr'):
#         """get position of A site for ribosome profiling.
#         includes the offset of the distance between 5' end and the A site.
#         * trpr: pr accession
#         * rna_pos: a list of rna positions
#         * offset: offset between 5' end and start of A site in ribosome 
#         """
#         pos = self.get_trpr_pos(trpr,level=level)
#         off = offset
#         index = rna_pos.index(pos[0]) # index of pr start pos in rna positions
#         if index >= off:
#                 new_pos = rna_pos[index-off:index] + pos[:-off]
#         else:  # index < offset
#             if pos[0]>pos[1]: # - strand
#                 step = -1
#             elif pos[0]<pos[1]:
#                 step = 1
#             new_pos = range(rna_pos[0]-step*(off-index),rna_pos[0],step) + rna_pos[:index]+pos[:-off]
#         return new_pos
    
#     def htseq_genome_array_set(self,stranded=False):
#         """
#         This function generates a genomic array of set for finding the position that are ovelapped by genes
#         * stranded: logical. 
#         """
#         gene_array_set = ht.GenomicArrayOfSets('auto',stranded=stranded)
#         for row in self.df.itertuples(index=False):
#             try:
#                 gene_array_set[ht.GenomicInterval(row[0],row[1],row[2],row[5])] += row[3]
#             except:
#                 print row,'has problem'
#                 continue
#         return gene_array_set
       
#     def genes_cov_target_pos(self,chrom,pos):
#         """
#         This function returns the genes that cover the position.
#         * chrom: str. chromosome.
#         * pos: int. position. 
#         return list of genes containing that position
#         """
#         chr_df = self.df[self.df['chr'].values==chrom]
#         chr_df = chr_df[(chr_df['start']< pos) & (chr_df['end']>=pos)]
#         genes = list(set(chr_df['geneid'].tolist()))
#         return genes
    
#     def all_pr_start_end_pos(self):
#         """This function gets all the start and end position of all proteins.
#         For end position it is the start of the stop codon.
#         """
#         proteins = list(set(self.df['access'].tolist()))
#         chr_pr_start_end_sites = pd.DataFrame()
#         Chr = [];GeneID=[];Start=[];Stop=[];Pr=[];Strand=[]
#         for pr in proteins:
#             pr_df = self.df[self.df['access'].values==pr]
#             pr_obj = trpr(pr_df)
#             try:
#                 pos = pr_obj.get_trpr_pos(pr)
#             except:
#                 print pr,'map to multiple chromosome'
#                 continue
#             strand = pr_df['strand'].tolist()
#             # assign
#             Chr.append(pr_df['chr'].tolist()[0])
#             GeneID.append(pr_df['geneid'].tolist()[0])
#             Start.append(pos[0])
#             Stop.append(pos[-3])
#             Pr.append(pr)
#             Strand.append(strand[0])
#         chr_pr_start_end_sites['Chr']=pd.Series(Chr)
#         chr_pr_start_end_sites['GeneID']=pd.Series(GeneID)
#         chr_pr_start_end_sites['Start'] = pd.Series(Start)
#         chr_pr_start_end_sites['Stop'] = pd.Series(Stop)
#         chr_pr_start_end_sites['Pr_Access'] = pd.Series(Pr)
#         chr_pr_start_end_sites['Strand'] = pd.Series(Strand)
        
#         return chr_pr_start_end_sites
    
#     def get_exn_cds_start_end_pos(self,gene):
#         """
#         This function gets the start and end position of exons or cds of a gene provided
#         returns 2 list of start position and end position
#         """
#         gene_df = self.df[self.df['geneid'].values==gene]
#         strand = list(set(gene_df['strand'].tolist()))[0]
#         if strand == '+':
#             exn_start = gene_df['start'].tolist()
#             exn_end = gene_df['end'].tolist()
# #             ex_start = [abs(n-exn_start[0]) for n in exn_start]
# #             ex_end = [abs(n-exn_start[0]) for n in exn_end]
#         else:
#             exn_start = gene_df['end'].tolist()
#             exn_end = gene_df['start'].tolist()
#             exn_start.reverse()
#             exn_end.reverse()
# #             ex_start = [abs(n-exn_end[-1]) for n in exn_start]
# #             ex_end = [abs(n-exn_end[-1]) for n in exn_end]
#         return exn_start,exn_end