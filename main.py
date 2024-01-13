# A small tool for obtaining chemical space based on functional groups
# Author: Chemromi
# Description: From an initial molecule and a given functional group, 
#              all molecules after gradually increasing functional groups 
#              are obtained and stored as a bidirectional graph.
# Copyright: (C) 2022 Chemromi
# License: Apache License 2.0

import csv
import random
import copy
import numpy as np 
import pandas as pd 

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw, MolFromSmiles, MolToSmiles
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, MolFromSmiles, MolToSmiles, AllChem, MACCSkeys, Fragments 

import dgl
from dgl import DGLGraph

from chemutils import (attach_mols_nx, copy_edit_mol, decode_stereo,
                        enum_assemble_nx, set_atommap)

import warnings

warnings.filterwarnings("ignore")
rdkit.RDLogger.DisableLog('rdApp.*')

def worker_init_fn(id_):
    lg = rdkit.RDLogger.logger()
    lg.setLevel(rdkit.RDLogger.CRITICAL)



class CombineMolTree(DGLGraph):
    def __init__(self, smiles1,smiles2):
        DGLGraph.__init__(self)
        self.nodes_dict = {}
        self.add_nodes(2)
        mol1 = MolFromSmiles(smiles1)
        Chem.Kekulize(mol1)
        smiles1 = Chem.MolToSmiles(mol1, kekuleSmiles=True)
        mol2 = MolFromSmiles(smiles2)
        Chem.Kekulize(mol2)
        smiles2 = Chem.MolToSmiles(mol2, kekuleSmiles=True)
        self.nodes_dict[0] = dict(
                smiles=smiles1,
                mol=mol1,
                clique=[], 
                nid = 0,
            )
        self.nodes_dict[1] = dict(
                smiles=smiles2,
                mol=mol2,
                clique=[], 
                nid = 0,
            )
        self.add_edges(np.array([0,1]), np.array([1,0]))

def dfs_assemble(mol_tree_msg, mol_vec, cur_mol,
                 global_amap, fa_amap, cur_node_id, fa_node_id):
    nodes_dict = mol_tree_msg.nodes_dict
    fa_node = nodes_dict[fa_node_id] if fa_node_id is not None else None
    cur_node = nodes_dict[cur_node_id]

    fa_nid = fa_node['nid'] if fa_node is not None else -1
    prev_nodes = [fa_node] if fa_node is not None else []

    children_node_id = [v for v in mol_tree_msg.successors(cur_node_id).tolist()
                        if nodes_dict[v]['nid'] != fa_nid]
    children = [nodes_dict[v] for v in children_node_id]
    neighbors = [nei for nei in children if nei['mol'].GetNumAtoms() > 1]
    neighbors = sorted(
        neighbors, key=lambda x: x['mol'].GetNumAtoms(), reverse=True)
    singletons = [nei for nei in children if nei['mol'].GetNumAtoms() == 1]
    neighbors = singletons + neighbors

    cur_amap = [(fa_nid, a2, a1)
                for nid, a1, a2 in fa_amap if nid == cur_node['nid']]
    cands = enum_assemble_nx(cur_node, neighbors, prev_nodes, cur_amap)
    if len(cands) == 0:
        return None
    cand_smiles, cand_mols, cand_amap = list(zip(*cands))
    backup_mol = Chem.RWMol(cur_mol)
    ret_mols = []
    for i in range(len(cand_smiles)):
        cur_mol = Chem.RWMol(backup_mol)
        pred_amap = cand_amap[i]
        new_global_amap = copy.deepcopy(global_amap)

        for nei_id, ctr_atom, nei_atom in pred_amap:
            if nei_id == fa_nid:
                continue
            new_global_amap[nei_id][nei_atom] = new_global_amap[cur_node['nid']][ctr_atom]

        cur_mol = attach_mols_nx(cur_mol, children, [], new_global_amap)
        new_mol = cur_mol.GetMol()
        new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(new_mol))

        if new_mol is None:
            continue
        result = True
        
        end = True
        for nei_node in children:
            if not nei_node['is_leaf']:
                end = False
        if end:
            ret_mols = ret_mols + [cur_mol]
        
        for nei_node_id, nei_node in zip(children_node_id, children):
            if nei_node['is_leaf']:
                continue
                
            cur_mols = dfs_assemble(
                mol_tree_msg, mol_vec, cur_mol, new_global_amap, pred_amap,
                nei_node_id, cur_node_id)
            ret_mols = ret_mols + cur_mols
            if len(cur_mols)==0:
                return []
    if result:
        return ret_mols
    else:
        return []

def tree_restruct(mol_tree, nodes_dict, effective_nodes):
    effective_nodes_list = effective_nodes
    nodes_dict = [nodes_dict[v] for v in effective_nodes_list]
    for i, (node_id, node) in enumerate(zip(effective_nodes_list, nodes_dict)):
            node['idx'] = i
            node['nid'] = i + 1
            node['is_leaf'] = True
            if mol_tree.in_degrees(node_id) > 1:
                node['is_leaf'] = False
                set_atommap(node['mol'], node['nid'])

    mol_tree_sg = mol_tree.subgraph(effective_nodes)
    #mol_tree_sg.copy_from_parents()
    mol_tree_msg = mol_tree_sg
    #mol_tree_msg = unbatch(mol_tree_msg)[0]
    mol_tree_msg.nodes_dict = nodes_dict

    cur_mol = copy_edit_mol(nodes_dict[0]['mol'])
    global_amap = [{}] + [{} for node in nodes_dict]
    global_amap[1] = {atom.GetIdx(): atom.GetIdx()
                      for atom in cur_mol.GetAtoms()}
    mol_vec = 0
    cur_mol = dfs_assemble(
        mol_tree_msg, mol_vec, cur_mol, global_amap, [], 0, None)
    #print(cur_mol)
    if cur_mol is None:
        return []
    return cur_mol 

def mol_dispersion(smiles):
    mol_tree = DGLMolTree(smiles)
    #mol_tree.recover()
    #mol_tree.assemble()
    nodes_dict = mol_tree.nodes_dict
    effective_nodes = list(mol_tree.nodes_dict)
    mols = tree_restruct(mol_tree, nodes_dict, effective_nodes)
    slist = []
    for m in mols:
        for atom in m.GetAtoms():
            atom.SetAtomMapNum(0)
        smi = MolToSmiles(m)
        slist.append(smi)
    
    return slist


def mol_dispersion2(smiles1,smiles2):
    mol_tree = CombineMolTree(smiles1,smiles2)
    #mol_tree.recover()
    #mol_tree.assemble()
    nodes_dict = mol_tree.nodes_dict
    effective_nodes = list(mol_tree.nodes_dict)
    mols = tree_restruct(mol_tree, nodes_dict, effective_nodes)
    slist = []
    for m in mols:
        for atom in m.GetAtoms():
            atom.SetAtomMapNum(0)
        smi = MolToSmiles(m)
        slist.append(smi)
    
    return slist


def get_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        #print(smiles)
        return None
    #Chem.Kekulize(mol)
    return mol

def get_mols_by_gen(df,gen):
    # 根据代际获得分子
    df1 = df[(df['gen']==gen )]
    df1 = df1[['index','gen','smiles']]
    back_dict = df1.to_dict('list')
    return back_dict
def get_csv_index(df):
    index_num = df.shape[0]
    return index_num
def write_line_csv(df,content):
    idx_now = get_csv_index(df)
    content = [idx_now] + content
    df.loc[idx_now] = content
    #print(idx_now,end=" ")
    return idx_now

def mol_fps(smi):
    m = Chem.MolFromSmiles(smi)
    fps1 = AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=1024) 
    return fps1

def check_same_fps(fps1,fps2):
    return DataStructs.FingerprintSimilarity(fps1,fps2,metric=DataStructs.DiceSimilarity)

def check_same(fomu_k, fps_k):
    global fp
    global fomu
    index_list = [index for index, value in enumerate(fomu) if value == fomu_k]
    matching_indices = [index for index in index_list if check_same_fps(fp[index], fps_k) == 1]

    if matching_indices:
        return matching_indices
    else:
        return []
    
def mol_sim(smi1,smi2):
    """molecular Similarity"""
    m1 = Chem.MolFromSmiles(smi1)
    m2 = Chem.MolFromSmiles(smi2)
    fps1 = AllChem.GetMorganFingerprintAsBitVect(m1,2,nBits=1024) 
    fps2 = AllChem.GetMorganFingerprintAsBitVect(m2,2,nBits=1024)
    r=DataStructs.FingerprintSimilarity(fps1,fps2,metric=DataStructs.DiceSimilarity)
    return r

def mol_same(smi1,smi2):
    """"Returns whether two molecules are the same
    
    Except:
        Sometimes somethings go wrong
    """
    try:
        s = mol_sim(smi1,smi2)
        if s == 1: 
            return True
        return False
    except:
        print("warning : error in mol_same func")
        return False
def get_mol_fomula(smi):
    """Function as name"""
    mol = MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except:
        pass
    #atoms = mol.GetAtoms()
    return (Chem.rdMolDescriptors.CalcMolFormula(mol))

def initial_data(mol_csv,vec_csv,smi0):
    csvheader = ["index","gen",'fomula',"smiles",'value','rank']
    csvfile  = open(mol_csv, "w",newline='')
    writer = csv.writer(csvfile)
    writer.writerow(csvheader)
    global fomu
    global smilist
    global fp
    fomu = []
    smilist = []
    fp = []
    for i in range(len(smi0)):
        idx = i
        f = get_mol_fomula(smi0[i])
        content = [i,1,f,smi0[i]]
        writer.writerow(content)
        fomu.append(f)
        smilist.append(smi0[i])
        fp.append(mol_fps(smi0[i]))
    csvfile.close()
    csvheader = ["index","start","end",'smiles_start','smiles_end',"smiles_delta",'value']
    csvfile  = open(vec_csv, "w",newline='')
    writer = csv.writer(csvfile)
    writer.writerow(csvheader)
    csvfile.close()
    print(get_mols_by_gen(pd.read_csv(mol_csv),1))
    
def gen_map(mol_csv,vec_csv,gen,lines):
    global fomu
    global smilist
    global fp 
    dfm = pd.read_csv(mol_csv)
    dfv = pd.read_csv(vec_csv)
    dictm = dfm.to_dict(orient='records')
    dictv = dfv.to_dict(orient='records')
    smidict = get_mols_by_gen(dfm,gen)
    for i in range(len(smidict['index'])):
        print('gen: ',gen," ",i," in ",len(smidict['smiles']))
        for j in range(len(lines)):
            s = lines[j]
            try:
                slist = mol_dispersion2(smidict['smiles'][i],s)
            except:
                continue
            slist = list(set(slist))
            
            for k in range(len(slist)):
                fps_k = mol_fps(slist[k])
                fomu_k = get_mol_fomula(slist[k])
                same = check_same(fomu_k, fps_k)
                if same:
                    #print("same",same)
                    content = {'index': len(dictv), 'start': smidict['index'][i], 'end': same[0], 'smiles_start': smidict['smiles'][i], 
                                   'smiles_end': slist[k], 'smiles_delta': "+"+s, 'value': 0}
                    dictv.append(content)
                    #write_line_csv(dfv,content)
                    content = {'index': len(dictv), 'start': same[0], 'end': smidict['index'][i], 'smiles_start': slist[k], 
                                   'smiles_end': smidict['smiles'][i], 'smiles_delta': "-"+s, 'value': 0}
                    dictv.append(content)
                    #content = [same[0],smidict['index'][i],slist[k],smidict['smiles'][i],"-"+s,0]
                    #write_line_csv(dfv,content)
                else: 
                    #print(get_mol(slist[k]))
                    if not get_mol(slist[k]):
                        print("###########None############")
                        continue
                    content = {'index': len(dictm), 'gen': gen+1, 'fomula': fomu_k, 'smiles': slist[k], 'value': 0.0, 'rank': 0.0}
                    #content = [gen+1,get_mol_fomula(slist[k]),slist[k],0,0]
                    dictm.append(content)
                    new_idx = len(dictm)
                    content = {'index': len(dictv), 'start': smidict['index'][i], 'end': new_idx, 'smiles_start':smidict['smiles'][i], 
                                   'smiles_end': slist[k], 'smiles_delta': "+"+s, 'value': 0}
                    #content = [smidict['index'][i],new_idx,smidict['smiles'][i],slist[k],"+"+s,0]
                    dictv.append(content)
                    content = {'index': len(dictv), 'start': new_idx, 'end': smidict['index'][i], 'smiles_start': slist[k], 
                                   'smiles_end': smidict['smiles'][i], 'smiles_delta': "-"+s, 'value': 0}
                    #content = [new_idx,smidict['index'][i],slist[k],smidict['smiles'][i],"-"+s,0]
                    dictv.append(content)
                    fomu.append(fomu_k)
                    smilist.append(slist[k])
                    fp.append(fps_k)
    
    dfm = pd.DataFrame(dictm)
    dfv = pd.DataFrame(dictv)
    dfm.to_csv(mol_csv,index=False)
    dfv.to_csv(vec_csv,index=False)
    
if __name__=='__main__':
    import argparse
    # Initialize the parser
    parser = argparse.ArgumentParser(description='Process some paths and lists.')

    # Add arguments
    parser.add_argument(
        '--mpath',
        type=str,
        required=True,
        help='path to the mol file'
    )

    parser.add_argument(
        '--vpath',
        type=str,
        required=True,
        help='path to the vec file'
    )

    parser.add_argument(
        '--input',
        type=str,
        nargs='+',
        required=True,
        help='list of functional groups'
    )

    parser.add_argument(
        '--lines',
        type=str,
        nargs='+',
        required=True,
        help='list of functional groups'
    )

    parser.add_argument(
        '--num',
        type=int,
        required=True,
        help='number of generation iters'
    )

    # Parse the arguments
    args = parser.parse_args()

    # Access the arguments
    mpath = args.mpath
    vpath = args.vpath
    input_list = args.input
    lines_list = args.lines
    num = args.num

    # Below this line, you can add the logic to use the parsed arguments
    # For example, just printing them:
    print(f"Mol CSV Path: {mpath}")
    print(f"Vec CSV Path: {vpath}")
    print(f"Input List: {input_list}")
    print(f"Functional Groups List: {lines_list}")
    print(f"Number of generation iters: {num}")

    # If you want to run this script, save it to a file like `parse_args.py`
    # and run it using the command line, for example:
    # python parse_args.py --mpath "model.pth" --vpath "video.mp4" --input "input1" "input2" --lines "line1" "line2" --num 42
    
    initial_data(mpath,vpath,input_list)
    for i in range(num-1):
        gen_map(mpath,vpath,i+1,lines_list)