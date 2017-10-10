'''
this file extracts cds and exn position into bed file format

'''

import argparse,os
from f01_parse_gff import ncbi_gff
from f03_parse_trpr_df import trpr
import pandas as pd

parser = argparse.ArgumentParser(description='Prepare related annotaion file for downstream analysis')
parser.add_argument('-g','--gff',action='store',dest='gff_fn',help='path to gff file')
args = parser.parse_args()


gff_fn = args.gff_fn
if not gff_fn.endswith('gff'):
    assert False,'input file is not in gff format'
db_path = os.path.dirname(gff_fn)


#################### 1. get CDS and exn position bed file ###############
def generate_feature_pos_df(gff_obj,outFile,feature):
    '''outFile is the file stores the file'''
    if not os.path.exists(outFile):
        df = gff_obj.get_feature_pos_df(feature)
        if feature == 'CDS':
            df = df.sort_values(['Chr','cds_s'])
        else:
            df = df.sort_values(['Chr','exn_s'])
        df.to_csv(outFile,sep='\t',index=False)

# extract cds and exon positions separately. These two dataframe is used for 
# extracting protein or rna sequence positions for downstream analysis.
gff_df = pd.read_csv(gff_fn,sep='\t',header=None,comment='#')
gff_obj = ncbi_gff(gff_df)
if not os.path.exists(db_path): os.mkdir(db_path)
exnFile = db_path + '/01_pr_rna.bed'
if not os.path.exists(exnFile):
    generate_feature_pos_df(gff_obj,exnFile,'exon')
cdsFile = db_path + '/01_pr_cds.bed'
if not os.path.exists(cdsFile):
    generate_feature_pos_df(gff_obj,cdsFile,'CDS')
print('exn and cds position file done')


################## 2. get all gene id files ##################

all_id_fn = db_path + '/02_gff_allid.txt'
if not os.path.exists(all_id_fn):
    all_id_df = gff_obj.get_all_id()
    all_id_df.to_csv(all_id_fn,sep='\t',index=False)

print('id file preparation done')
################# 3. remove proteins that map to multiple scaffolds ##############
# Get proteins that map to multiple scaffolds, we remove them in the downstream analysis
multi_chr_prs = gff_obj.multi_chr_protein()
cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
cds_df = cds_df[~cds_df['PrAccess'].isin(multi_chr_prs)]
cds_df.to_csv(cdsFile,sep='\t',index=False)
print('remove proteins with multiple chromosomes in cds position file done')

multi_chr_prs_fn = db_path + '/01_multi_chr_proteins.txt'
if not os.path.exists(multi_chr_prs_fn):
    with open(multi_chr_prs_fn,'w') as f:
        f.write('\n'.join(multi_chr_prs))
################# 4. get length of 5UTR and 3UTR ##############
'''
Get length of 5'UTR and 3'UTR for all proteins. This information 
is used to get the offset of P sites in riboseq. Sometimes there 
are splice sites between TSS of rna and TSS of protein sequence, 
if we want to get upstream sequence of CDS and upstream base pairs
pass the splice sites, we need to skip the intron sequence upstream 
of CDS to get the correct upstream rna base pairs of CDS. 
The format of this information is in the following format: 
[geneid,traccess,praccess,t_start,p_start,dist_start,t_end,p_end,dist_end]'''


def get_pos_dic(bed):
    '''this function builds dictioanry for protein or rna {id:[(s1,e1),(s2,e2)]}'''
    dic = {}
    with open(bed) as f:
        for line in f:
            if line.startswith('Chr'): continue
            item = line.strip().split('\t')
            if item[4] in dic:
                dic[item[4]].append((item[1],item[2]))
            else:
                dic[item[4]] = [item[-1],(item[1],item[2])]
    return dic


def get_pos(dic,access):
    strd = dic[access][0]
    pos = []
    for p in dic[access][1:]:
        pos.extend(range(int(p[0])+1,int(p[1])+1))
    if strd == '-':
        pos = pos[::-1]
    return pos

def get_utr_len(row,exn_dic,cds_dic):
    tr = row['TrAccess']
    pr = row['PrAccess']
    if pr not in cds_dic:
        return pd.Series([0,0])
    if tr not in exn_dic or tr == '-':
        return pd.Series([0,0])
    tr_pos = get_pos(exn_dic,tr)
    pr_pos = get_pos(cds_dic,pr)
    utr5 = tr_pos.index(pr_pos[0])
    utr3 = len(tr_pos) - tr_pos.index(pr_pos[-1]) - 1
    if tr_pos[0] < tr_pos[1]: 
        strd = '+'
    else:
        strd = '-'
    if pr_pos[0] < pr_pos[1]:
        p_strd = '+'
    else:
        p_strd = '-'
    if strd != p_strd: 
        assert False,tr+' and ' + pr + ' strand inconsistant'
    return pd.Series([utr5,utr3,strd])




utr_len_fn = db_path +'/03_utr_len.txt'
if not os.path.exists(utr_len_fn):    
    cdsFile = db_path + '/01_pr_cds.bed'
    exnFile = db_path + '/01_pr_rna.bed'

    cds_pos_dic = get_pos_dic(cdsFile)
    exn_pos_dic = get_pos_dic(exnFile)
    all_id_fn =  db_path + '/02_gff_allid.txt'
    all_id_df = pd.read_csv(all_id_fn,sep='\t',header=0)
    utr_df = all_id_df[all_id_df['PrAccess'].values != '-']
    utr_df = utr_df.reset_index(drop=True)
    utr_df[['utr5_len','utr3_len','strand']] = utr_df.apply(lambda row: get_utr_len(row,exn_pos_dic,cds_pos_dic),axis=1)
    utr_df = utr_df[~utr_df['PrAccess'].isin(multi_chr_prs)]
    utr_df.to_csv(utr_len_fn,sep='\t',index=False)
print('utr5 and utr3 length done')

############### 5. get 5UTR and 3UTR for longest proteins ##############
long_utr_fn = db_path + '/04_long_utr_len.txt'
utr_len_fn = db_path +'/03_utr_len.txt'
if not  os.path.exists(long_utr_fn):
    # get longest proteins
    cds_df = pd.read_csv(cdsFile,sep='\t',header=0,low_memory=False)
    cds_obj = trpr(cds_df)
    long_pr = cds_obj.get_longest_trprs()
    long_prs = list(set(long_pr['feid']))

    # extract longest proteins utr
    utr_df = pd.read_csv(utr_len_fn,sep='\t',header=0,low_memory=False)
    utr_df = utr_df[utr_df['PrAccess'].isin(long_prs)]
    utr_df.to_csv(long_utr_fn,sep='\t',index=False)

