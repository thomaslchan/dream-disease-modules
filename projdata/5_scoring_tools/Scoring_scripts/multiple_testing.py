'''
Created on Oct 24, 2016

@author: sarvenaz choobdar

Functions for multiple testing correction of modules scored  by PASCAL tool.

To use this script  provide the following input at the command line:
    - the path to the PASCAL output containing the module scores for all the modules.
    - the name of the geneset file  scored by PASCAL


'''
import os,csv,sys
import numpy as np
from statsmodels.stats.multitest import fdrcorrection

#### text file containing the list of gwas for multiple testing correction
_gwas_collection_filelist ='gwas_final.txt'

def multiple_testing_correction(ps, alpha=0.05,
        method='benjamini-hochberg', **kwargs):
    """ correct pvalues for multiple testing and add corrected `q` value
    
    :param ps: list of pvalues
    :param alpha: significance level default : 0.05
    :param method: multiple testing correction method [bonferroni|benjamini-hochberg]
    :returns rej: rejected nodes
    """
    _p = np.array(ps)
    q = _p.copy()
    rej = _p.copy()
    mask = ~np.isnan(_p)
    p = _p[mask]
    if method == 'bonferroni':
        q[mask] = p / len(p)
        rej[mask] = q[mask] < alpha
    elif method == 'benjamini-hochberg':
        _rej, _q = fdrcorrection(p, alpha)
        rej[mask] = _rej
        q[mask] = _q
    else:
        raise ValueError(method)
    return rej



def multiple_testing(path_to_pascal_output,gene_set_fileName,selected_gwas,method='fdr',alpha = 0.05):
    '''
    input: 
        - path_to_pascal_output: path to the pascal output
        - gene_set_fileName: name of file contating the geneset
        - selected_gwas: list of gwas for which the  gene_set_fileName is scores
        - method: correction method: fdr (benjamini-hochberg') , bn (bonferroni)
    output: 
        - N: Number of modules in geneSet file
        - NS: Number of modules that showed enrichment at least for one of the gwas in the selected_gwas list
    '''
    gene_set_fileName = os.path.splitext(os.path.basename(gene_set_fileName))[0]
    N=0
    onlyfiles=[]
    for g in selected_gwas:
        onlyfiles += [ f for f in os.listdir(path_to_pascal_output) if g+'.PathwaySet--' +gene_set_fileName in f ]
    
    if len(onlyfiles)==0:
        print '\nNo PASCAL results found'
        return None
    sig_mids = []
    for filename in onlyfiles: # phenotype.PathwasySet.teamI(/userID).submissionID.networkFileName.sum.txt: 3342188.6971288.6_homology_anonym_v2.txt
        path2 =os.path.join(path_to_pascal_output,filename)  
        L = csv.reader( open(path2, 'rU'),delimiter='\t')
        L.next()
        pval =[]
        mids = []
        _NG = 0
        for a in L:
            _NG+=1
            if a[1]!='NA':
                pval.append(float(a[1]))
                mids.append(int(a[0]))
        N =_NG
        if method=='bn':
            NS_corected = list((np.array(pval)* N< alpha)*1)
        if method=='fdr':
            NS_corected = list(multiple_testing_correction(np.array(pval)))
        sig_mids +=[mids[i] for i in xrange(len(mids)) if  NS_corected[i]==1]
        NS =len(set(sig_mids))
    return NS

if __name__ == '__main__':
    selected_gwas=[]
    L = csv.reader( open(_gwas_collection_filelist ,'rU'),delimiter='\t')
    for a in L:
        f=a[0]
        selected_gwas.append('.'.join(f.split('.')[:-2])) #EUR.IBDGenetics.UC.txt.gz
    ns = multiple_testing(sys.argv[1],sys.argv[2],selected_gwas)
    print '\nMultiple testing correction of PASCAL output for genesetfile:',sys.argv[2]
    print 'Number of Significant modules:',ns
