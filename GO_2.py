## Script for analysing the gene enrichment of a set of differentially expressed genes.
## Input: file with differentially expressed genes, file with identified genes in experiment
##		  file with process category-terms for all genes in species. 

#Last edited: 09/01/2018

#input paths
path_category = "../data/openMS/results/goterm.csv"
path_sign = "../data/openMS/statistics/sign_genes.csv"
path_iden = "../data/openMS/statistics/diffacto_q_values.csv"

#output path
path_sign_out = "../data/openMS/enrichment/GO/sign_process.csv"
path_iden_out = "../data/openMS/enrichment/GO/identified_process.csv"
path_process = "../data/openMS/enrichment/GO/enriched_GO.csv"



import pandas as pd
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import interpolate


#Use q-value estimation by nfusi:
def estimate(pv, m=None, verbose=False, lowmem=False, pi0=None):
    """
    Estimates q-values from p-values

    Args
    =====

    m: number of tests. If not specified m = pv.size
    verbose: print verbose messages? (default False)
    lowmem: use memory-efficient in-place algorithm
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
         For most GWAS this is not necessary, since pi0 is extremely likely to be
         1

    """
    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = sp.arange(0, 0.90, 0.01)
        counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = sp.array(pi0)

        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        if verbose:
            print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)

        if pi0 > 1:
            if verbose:
                print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
            pi0 = 1.0

    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0

    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = sp.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
        pv[last_pv] = -sp.inf
        prev_qv = last_pv
        for i in xrange(int(len(pv))-2, -1, -1):
            cur_max = pv.argmax()
            qv_i = (pi0*m*pv[cur_max]/float(i+1))
            pv[cur_max] = -sp.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]

    else:
        p_ordered = sp.argsort(pv)
        pv = pv[p_ordered]
        qv = pi0 * m/len(pv) * pv
        qv[-1] = min(qv[-1], 1.0)

        for i in range(len(pv)-2, -1, -1):
            qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])

        # reorder qvalues
        qv_temp = qv.copy()
        qv = sp.zeros_like(qv)
        qv[p_ordered] = qv_temp

    # reshape qvalues
    qv = qv.reshape(original_shape)

    return qv
	
#import file with significantly expressed genes:
sign_df = pd.read_csv(path_sign, header=0, index_col=0, usecols=["protein","f_value","p_value","q_value"])
cg_df = pd.read_csv(path_category, header=0, index_col=0, usecols=["gene_id","go_description","go_id"])
iden_df = pd.read_csv(path_iden, header=0, index_col=0, usecols=["protein","f_value","p_value","q_value"])


###calc qvalue by me###
def get_q_value(pvals):
	q_sign = 0.05
	num = len(pvals)
	q_vals = []
	sign = []
	##starting from the lowest ranked p-value:
	prev_q = pvals[-1]
	q_vals.append(prev_q)
	if prev_q < q_sign:
		sign.append("Yes")
	else:
		sign.append("No")
	##minimum for those with higher rank
	for rank in range (num-1,0,-1):
		q_val = min(pvals[rank-1]*num/rank , prev_q)
		q_vals.append(q_val)
		prev_q = q_val
		##Find significant results
		if q_val < q_sign:
			sign.append("Yes")
		else:
			sign.append("No")
	return q_vals[::-1], sign[::-1]


#new part for GO
##############################fold change = (sign/tot_sign)/(iden/tot_sign)
s_df = pd.DataFrame(columns=["gene_id","go_description","go_id"])
s_df.set_index("gene_id", inplace=True)
no_go = list()
print(len(sign_df.index.tolist()))
for gene in sign_df.index.tolist():
	s_ = cg_df.loc[cg_df.index == gene]
	s_df=s_df.append(s_)
	if gene not in cg_df.index:
		no_go.append(gene)
print(s_df)

s_df.rename(columns={"go_id":"go_significant"},inplace=True)
#print(no_go ,len(no_go), sep="\n") #191 not in go_term list ! 533 present.
GO_counts_sign = s_df["go_significant"].value_counts()	#counts for all GOs for sign set
N = len(s_df.index.value_counts())

##IDENTIFIED##
i_df = pd.DataFrame(columns=["gene_id","go_description","go_id"])
i_df.set_index("gene_id", inplace=True)
no_go_i = list()
for gene in iden_df.index.tolist():
	i_ = cg_df.loc[cg_df.index == gene]
	i_df=i_df.append(i_)
	if gene not in cg_df.index:
		no_go_i.append(gene)

i_df.rename(columns={"go_id":"go_identified"},inplace=True)
#print(no_go_i ,len(no_go_i),sep="\n") #437 not in go_term list ! 958 present.
GO_counts_iden = i_df["go_identified"].value_counts()	#counts for all GOs for iden set
M = len(i_df.index.value_counts())
GOsign_df= pd.DataFrame(GO_counts_sign)
GOiden_df = pd.DataFrame(GO_counts_iden)
df = GOsign_df.merge(GOiden_df, left_on = ["go_significant"],right_on=["go_identified"],left_index=True,right_index=True)
fc = [(sign/N)/(iden/M) for sign,iden in zip(df["go_significant"],df["go_identified"])]

df = df.assign(fold_change=fc, missing_GO_significant=[len(no_go)]*len(df), \
missing_GO_identified=[len(no_go_i)]*len(df), total_significant_N=[N]*len(df), \
total_identified_M=[M]*len(df))



###Statistics: Hypergeometric test###
x=df["go_significant"]	#GO counts in significant
n=df["go_identified"]	#GO counts in identified
pvalsH = stats.hypergeom.sf(x-1, M, n, N, loc=0)	#perform hypergeometric test
pvalsL = []
#cdf wont accept my input as array/list/df. avoiding the struggle:
for xL,nL in zip(x,n):
	pvalL = stats.hypergeom.cdf(xL, M, nL, N, loc=0)
	pvalsL.append(pvalL)

df = df.assign(p_value_enrichment=pvalsH,p_value_depletion=pvalsL)


#Sort according to p-value_depletion and add q-values
df = df.sort_values(by=["p_value_depletion"],ascending=True)
qvalsL, signL = get_q_value(df["p_value_depletion"])
qvL = estimate(df["p_value_depletion"])
df = df.assign(depletion_significant=signL,q_value_depletion=qvalsL,q_value_pi0_depletion=qvL)


#sort according to p-value enrichment and add q-values
df = df.sort_values(by=["p_value_enrichment"],ascending=True)
qvalsH, signH = get_q_value(df["p_value_enrichment"])
qvH = estimate(df["p_value_enrichment"])
df = df.assign(enrichment_significant=signH,q_value_enrichment=qvalsH,q_value_pi0_enrichment=qvH)

#output process counts dataframe
df=df[["go_significant","go_identified","total_identified_M","total_significant_N",\
"fold_change","missing_GO_identified","missing_GO_significant", \
"p_value_depletion","q_value_depletion","depletion_significant", "q_value_pi0_depletion",\
"p_value_enrichment","q_value_enrichment","enrichment_significant","q_value_pi0_enrichment"]]	
df.to_csv(path_process)




"""
###Plots########################################################################################################
counts_df_sign = counts_df.sort_values(by=["counts_significant"],ascending=True) #Sort according to #significant
num= len(counts_df)
fig, ax = plt.subplots(figsize=(13,7))
ind = np.arange(num)  # the x locations for the groups
width = 0.3       # the width of the bars
rects1 = ax.barh(ind, counts_df_sign["counts_identified"], width, color='green')
rects2 = ax.barh(ind + width, counts_df_sign["counts_significant"], width, color='orange')

# add some text for labels, title and axes ticks
ax.set_yticklabels(counts_df_sign.index)
ax.set_xlabel('Number of genes annotated')
ax.set_title('Enrichment of Categories')
ax.set_yticks(ind + width/2)

#save plot
ax.legend((rects1[0], rects2[0]), ('Identified', 'ANOVA q<0.05'))	
plt.subplots_adjust(left=0.4)
plt.savefig("../data/openMS/enrichment/new/enrichment_process_plot_new.png")
plt.close(fig)

###normalized plot###
norm=[x/y for x, y in zip( counts_df_pval["counts_significant"], counts_df_pval["counts_identified"])]
counts_df_pval = counts_df_pval.assign(significant_normalized=norm)
counts_df_norm = counts_df_pval.sort_values(by=["significant_normalized"],ascending=True) #sort according to norm-value
fig, ax = plt.subplots(figsize=(10,6))
rects = ax.barh(ind, counts_df_norm["significant_normalized"], width, color='#335B8E')
for i, bar in enumerate(rects):
	if counts_df_norm['enrichment_significant'][i] == 'Yes':
			bar.set_color("#B7DBDB")
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
#         fancybox=True, shadow=True, ncol=2)

# add some text for labels, title and axes ticks
ax.set_yticklabels(counts_df_norm.index)
ax.set_xlabel('Proportion of genes annotated')
ax.set_title('Normalised Enrichment of Categories')
ax.set_yticks(ind + width/2)

#save plot
colour_blue = matplotlib.patches.Patch(color='#335B8E', label='Nonsignificant')
colour_yel = matplotlib.patches.Patch(color="#B7DBDB", label='Significant enrichment')
#plt.legend(handles=[colour_blue,colour_yel], loc=8)
plt.legend(handles=[colour_blue,colour_yel], loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
plt.subplots_adjust(left=0.5, bottom=0.2)
plt.savefig("../data/openMS/enrichment/new/enrichment_process_plot_new_norm.png")
plt.close(fig)
"""
"""
###normalized plot###V2
norm=[x/M for x in (counts_df_pval["counts_identified"])]
counts_df_pval = counts_df_pval.assign(significant_normalized=norm)
counts_df_norm = counts_df_pval.sort_values(by=["significant_normalized"],ascending=True) #sort according to norm-value
fig, ax = plt.subplots(figsize=(13,6))
rects = ax.barh(ind, counts_df_norm["significant_normalized"], width, color='blue')
for i, bar in enumerate(rects):
	if counts_df_norm['enrichment_significant'][i] == 'Yes':
			bar.set_color("#edd012")

# add some text for labels, title and axes ticks
ax.set_yticklabels(counts_df_norm.index)
ax.set_xlabel('Proportion of genes annotated')
ax.set_title('Normalised Enrichment of Categories')
ax.set_yticks(ind + width/2)

#save plot
colour_blue = matplotlib.patches.Patch(color='blue', label='Nonsignificant')
colour_yel = matplotlib.patches.Patch(color="#edd012", label='Significant enrichment')
plt.legend(handles=[colour_blue,colour_yel])
plt.subplots_adjust(left=0.4)
plt.savefig("../data/openMS/enrichment/new/enrichment_process_plot_new_norm_V2.png")
plt.close(fig)
"""