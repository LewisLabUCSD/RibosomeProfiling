import argparse
import os,glob
import pandas as pd
from natsort import natsorted
from multiprocessing import Pool



parser = argparse.ArgumentParser(description='Get coverage at each position of transcript including 50 bps up and downstream bps')
parser.add_argument('-c','--cov_path',action='store',dest='cov_path',help='coverage path that stores all the coverage files')
parser.add_argument('-g','--gff_file',action='store',dest='gff',help='full path to gff file')
parser.add_argument('-o','--offset',action='store',dest='offset',help='full path to offset file')
parser.add_argument('-t','--thread',action='store',dest='thread',type=int,help='thread to use',default=1)
parser.add_argument('-e','--end',action='store',dest='end',type=int,help='which end of read to use',default=5)
args = parser.parse_args()

cov_path = args.cov_path
if cov_path.endswith('/'): cov_path = cov_path[:-1]
gff_fn   = args.gff
offset_fn= args.offset
thread   = args.thread
end      = args.end




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
    if access not in dic:
        pos = []
    else:
        strd = dic[access][0]
        pos = []
        for p in dic[access][1:]:
            pos.extend(range(int(p[0])+1,int(p[1])+1))
        if strd == '-':
            pos = pos[::-1]
    return pos


def cov5_3_dic(covFile,m_lens):
    '''prepare two dictionaries for mapping of 5end and 3end of the reads.
    format {chr:{pos+/-:count}}'''
    cov_5dic = {}
    cov_3dic = {}
    with open(covFile) as cov:
        for line in cov:
            item = line.strip().split('\t')
            if int(item[-1]) not in m_lens: continue
            count = int(item[0])
            chrom = item[1]
            end5  = item[2]
            end3  = item[3]
            strd  = item[4]
            if chrom in cov_5dic:
                if end5+strd in cov_5dic[chrom]:
                    cov_5dic[chrom][end5+strd] += count
                else:
                    cov_5dic[chrom][end5+strd] = count
                if end3+strd in cov_3dic[chrom]:
                    cov_3dic[chrom][end3+strd] += count
                else:
                    cov_3dic[chrom][end3+strd] = count
            else:
                cov_5dic[chrom] = {}
                cov_3dic[chrom] = {}
    return cov_5dic,cov_3dic


def get_pos_cov(dic5,dic3,chrom,pos,end):
    end = str(end)
    if pos[0]<pos[1]:
        strd = '+'
    else:
        strd = '-'
    pos_cov = []
    for p in pos:
        try:
            if end == '5' and strd == '+':
                pos_cov.append(dic5[chrom][str(p)+strd])
            elif end == '5' and strd == '-':
                pos_cov.append(dic3[chrom][str(p)+strd])
            elif end == '3' and strd == '+':
                pos_cov.append(dic3[chrom][str(p)+strd])
            elif end == '3' and strd == '-':
                pos_cov.append(dic5[chrom][str(p)+strd])
        except:
            pos_cov.append(0)
    return pos_cov



def get_full_tr_pos_cov(covFile,outpath,exnFile,cdsFile,long_utr_fn,offset_fn,end):
    # prepare rna {id:position} dictionary
    cds_pos_dic = get_pos_dic(cdsFile)
    exn_pos_dic = get_pos_dic(exnFile)
    # define outpath
    cov_fn_folder = outpath + '/' + os.path.basename(covFile)[:-4]
    if not os.path.exists(cov_fn_folder): os.mkdir(cov_fn_folder)
    # prepare offset dictionary
    off_df = pd.read_csv(offset_fn,sep='\t',header=0)
    offset_dic = {k:list(v) for k,v in off_df.groupby('offset')['len']}
    # prepare pr tr ids
    utr_df = pd.read_csv(long_utr_fn,sep='\t',header=0)
    utr_df.index = utr_df['PrAccess']
    pr_tr_dic = utr_df.set_index('PrAccess')['TrAccess'].to_dict()
    pr_chr_dic= utr_df.set_index('PrAccess')['Chrom'].to_dict()
    prs = pr_tr_dic.keys()
    
    for key in offset_dic:
        out_fn = cov_fn_folder + '/rna_pos_cov_off_' + str(key) + '.txt'
        with open(out_fn,'w') as f:
            m_lens = offset_dic[key]
            # prepare coverage dictionary
            dic5,dic3 = cov5_3_dic(covFile,m_lens)
            for pr in prs:
                tr = pr_tr_dic[pr]
                chrom = pr_chr_dic[pr]
                if tr in exn_pos_dic:
                    tr_pos = get_pos(exn_pos_dic,tr)
                else:
                    tr_pos = get_pos(cds_pos_dic,pr)
                    tr = pr
                if tr_pos[0] < tr_pos[1]:
                    tr_pos = list(range(tr_pos[0]-50,tr_pos[0]))+tr_pos+list(range(tr_pos[-1]+1,tr_pos[-1]+51))
                else:
                    tr_pos = list(range(tr_pos[0]+50,tr_pos[0],-1))+tr_pos+list(range(tr_pos[-1]-1,tr_pos[-1]-51,-1))
                if len(set(tr_pos)) != len(tr_pos):
                    assert False,'duplicate position in tr position'
                # get coverage
                pos_cov = get_pos_cov(dic5,dic3,chrom,tr_pos,str(end))
                pos_cov = [str(p) for p in pos_cov]
                f.write('\t'.join([tr,pr]+pos_cov)+'\n')



import time
start = time.time()

db_path  = os.path.dirname(gff_fn)
path = os.path.dirname(cov_path)

long_utr_fn = db_path + '/04_long_utr_len.txt'
exnFile = db_path + '/01_pr_rna.bed'
cdsFile = db_path + '/01_pr_cds.bed'
tr_pos_cov = path + '/04_tr_pos_cov'
if not os.path.exists(tr_pos_cov): os.mkdir(tr_pos_cov)



covFiles = natsorted(glob.glob(cov_path+'/*.txt'))

p = Pool(processes=thread)
for covFile in covFiles:
    p.apply_async(get_full_tr_pos_cov,args=(covFile,tr_pos_cov,exnFile,cdsFile,long_utr_fn,offset_fn,end,))
p.close()
p.join()
# get_full_tr_pos_cov(covFiles[0],tr_pos_cov,exnFile,cdsFile,long_utr_fn,offset_fn,end)
print('get position coverage for transcript finished')
print('total run time is: '+str((time.time()-start)/60) +' minutes')

