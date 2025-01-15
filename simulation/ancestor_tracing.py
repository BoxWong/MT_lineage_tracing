import pandas as pd
import re
import numpy as np
import pickle
from Bio import Phylo
import sys
from collections import Counter
from io import StringIO

def read_tree(file):
    with open(file, 'r') as f:
        tree_nwk = f.readline()
    tree_nwk = StringIO(tree_nwk)
    tree = Phylo.read(tree_nwk, format='newick')
    return tree

def ancestor_tracing_tg(cells, ancestors, tree, count=True):
    '''
    from Bio import Phylo
    tree = Phylo.read('xxx.nwk', format='newick')
    cells = ['<19_1015>', '<25_1878>', '<25_586>']
    ancestors = ['<3_0>', '<3_1>', '<3_4>', '<3_3>', '<3_5>', '<3_6>', '<3_7>', '<3_8>']
    ancestor_tracing(cells, ancestors, tree)
    '''
    ancestors = [tree.find_any(i) for i in ancestors]
    res = np.array([None]*len(cells))
    for anc in ancestors:
        res[np.isin(cells, [i.name for i in anc.get_terminals()])] = anc.name
    if count:
        return Counter(res)
    else:
        return res


def find_ancestor(tree, clone, cutoff):
    res = tree.common_ancestor(clone)
    level = int(res.name.split('_')[0][1:])
    coverage = []
    best_coverage = 1
    ress = []
    for i in Phylo.BaseTree._level_traverse(res, lambda elem: elem.clades):
        curr_level = int(i.name.split('_')[0][1:])
        if curr_level > level:
            level = curr_level
            if np.sum(np.array(coverage)>cutoff)==0:
                break
            else:
                res = ress[np.where(np.array(coverage)==max(coverage))[0][0]]
            coverage = []
            ress = []
        termi = [i.name for i in i.get_terminals()]
        coverage.append(np.sum(np.isin(clone, termi))/len(clone))
        ress.append(i)
        
    return res.name


sim_id = sys.argv[1]
gen = sys.argv[2]
clone_df = pd.read_csv(sys.argv[3])
#tree = Phylo.read(sys.argv[4], format= "newick")
tree = read_tree(sys.argv[4])
output_name = sys.argv[5]
ancestor_output_name = sys.argv[6]
anc_count_ls = []
most_anc_ls = []
most_anc_gen_ls = []
#define a pattern for parsing the cell generation
pattern = r"<(\d+)_(\d+)>"




df1 = pd.DataFrame(columns=["counts", "clone", "ancestor"])
for i in clone_df["Clone"].unique():
    #cell_id = pd.Series([k.split(">")[0]+">" for k in clone_df["Cell"][clone_df["Clone"]==i]]).drop_duplicates()
    cell_id = pd.Series([k.split(">")[0]+">" for k in clone_df["Cell"][clone_df["Clone"]==i]])
    
    #print(ancestor_tracing(cells=cell_id, ancestors = ['<3_0>', '<3_1>', '<3_4>', '<3_3>', '<3_5>', '<3_6>', '<3_7>', '<3_8>'], tree=tree))
    anc_count_ls.append(len(ancestor_tracing_tg(cells=cell_id, ancestors = ['<3_0>', '<3_1>', '<3_4>', '<3_3>', '<3_5>', '<3_6>', '<3_7>', '<3_8>'], tree=tree)))
    tmp_ser = pd.Series(ancestor_tracing_tg(cells=cell_id, ancestors = ['<3_0>', '<3_1>', '<3_4>', '<3_3>', '<3_5>', '<3_6>', '<3_7>', '<3_8>'], tree=tree), name="counts")
    tdf = tmp_ser.to_frame()
    tdf["clone"] = i
    tdf['ancestor'] = tdf.index
    df1 = pd.concat([df1, tdf], axis=0)
    
    ancestor = find_ancestor(tree, cell_id, 0.8)
    most_anc_ls.append(ancestor)
    most_anc_gen_ls.append(re.findall(pattern, ancestor)[0][0])
df1.to_csv(output_name)

df2 = pd.DataFrame({"Ancestor":most_anc_ls, "Ancestor_gen":most_anc_gen_ls}, index = clone_df["Clone"].unique())
df2.to_csv(ancestor_output_name)
