from shutil import copyfile
import numpy as np
import networkx as nx
import os
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

plt.rcParams.update({'font.size': 20})
rc('text', usetex=True)

def remove_hydrogens(pdb_file, save_copy=False):
    """Remove hydrogens from a pdb file
    
    Parameters
    ----------
    pdb_file : str
        pdb file path
    save_copy : boolean, optional
        if True, save a copy of the original file (default is False)
    """

    pdb_id = pdb_file.rsplit(".", 1)[0]
    copyfile(pdb_file, pdb_file.rsplit(".", 1)[0] + "_ORIG.pdb")

    with open(pdb_file, "r") as f:
        pdb = [line.rsplit("\n")[0] for line in f]
    pdb_new = []
    for line in pdb:
        if line[0:4] == "ATOM" and (line[13] == "H" or line[12] == "H"):
            pass
        else:
            pdb_new.append(line)
    with open(pdb_file, "w") as f:
        f.writelines([line + "\n" for line in pdb_new])

def local_partitioning_difference(variant1, variant2):
    lp1 = variant1.get_local_partitioning()
    lp2 = variant2.get_local_partitioning()

    if lp1.keys() != lp2.keys():
        raise Exception("The two protein variants must have the same sequence positions.")

    return {node: np.array(lp2[node]) - np.array(lp1[node]) for node in lp1.keys()}

def plot_local_partitioning_difference(variant1, variant2, folder_path, name=None, figsize=(20, 5), save_csv=False):
    chains = variant1.chains

    oligomer=False
    if len(chains) > 1:
        oligomer = True
    diff = local_partitioning_difference(variant1, variant2)

    for c in chains:
        net1 = variant1.net
        net2 = variant2.net
        nodes_chain = [n for n in net1.nodes if n[0] == c]
        positions = np.array(sorted([int(n[1::]) for n in nodes_chain]))
        residues1 = nx.get_node_attributes(net1, "residue_type")
        residues2 = nx.get_node_attributes(net2, "residue_type")
        seq1 = [residues1[c + str(n)] for n in positions]
        seq2 = [residues2[c + str(n)] for n in positions]

        xticks = [aa1 + str(pos) + aa2 if aa2 != aa1 else aa1 + str(pos) for aa1, pos, aa2 in zip(seq1, positions, seq2)]

        plt.figure(figsize=figsize)
        diff1D = [diff[c + str(n)][0] for n in positions]
        diff2D = [diff[c + str(n)][1] for n in positions]
        diff3D = [diff[c + str(n)][2] for n in positions]
        diff4D = [diff[c + str(n)][3] for n in positions]

        plt.plot(positions, [0] * len(positions), c="gray")
        w = .2
        plt.bar(positions - 2*w, diff1D, width=w, label=r"$\Delta f_{w_{1D}}$")
        plt.bar(positions - w, diff2D, width=w, label=r"$\Delta f_{w_{2D}}$")
        plt.bar(positions + w, diff3D, width=w, label=r"$\Delta f_{w_{3D}}$")
        if oligomer:
            plt.bar(positions + 2*w, diff4D, width=w, label=r"$\Delta f_{w_{4D}}$")
        plt.xticks(positions, xticks,
                rotation=90, fontsize=8)
        plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
        plt.ylim(-.5, .5)
        plt.tight_layout()
        plt.grid()
        if name:
            plt.savefig(os.path.join(folder_path, "%s_plot_differences_chain%s.png" %(name, c)))
            plt.savefig(os.path.join(folder_path, "%s_plot_differences_chain%s.pdf" %(name, c)))
        else:
            plt.savefig(os.path.join(folder_path, "plot_differences_chain%s.png" %c))
            plt.savefig(os.path.join(folder_path, "plot_differences_chain%s.pdf" %c))
        plt.show()

        if save_csv:
            df = pd.DataFrame()
            df["position"] = positions
            df["label"]= xticks
            df["diff1D"] = diff1D
            df["diff2D"] = diff2D
            df["diff3D"] = diff3D
            if oligomer:
                df["diff4D"] = diff4D
            df.to_csv(os.path.join(folder_path, "differences_chain%s.csv" %c), index=False)







