from biographs import Pmolecule
from biopandas.pdb import PandasPdb
import networkx as nx
import numpy as np
import amino_acids_conversion as aaconv

class Protein(Pmolecule):
    
    def __init__(self, pdb, pos_start=None, pos_stop=None, selected_chains=None):
        super().__init__(pdb)
        self.net = self.network()
        self.assign_residue_type()

        if selected_chains:
            self.remove_selected_chains
            self.chains = selected_chains
        else:
            self.chains = set([n[0] for n in self.net.nodes])

        self.pos_start = pos_start
        self.pos_stop = pos_stop
        self.remove_unwanted_positions()

        self.assign_secondary_structure()
        self.assign_relation_links()
        self.assign_local_partitioning()
    
    def remove_unwanted_chains(self, selected_chains):
        net = self.net
        nodes = self.nodes
        to_remove = [n for n in nodes if n[0] not in selected_chains]
        net.remove_nodes_from(to_remove)

    def remove_unwanted_positions(self):
        pos_start = self.pos_start
        pos_stop = self.pos_stop
        if not pos_start:
            pos_start = int(sorted(self.net.nodes)[0][1::])
        if not pos_stop:
            pos_stop = int(sorted(self.net.nodes)[-1][1::])
        net = self.net
        nodes = net.nodes
        to_remove = [n for n in nodes if int(n[1::]) < pos_start or int(n[1::]) > pos_stop]
        net.remove_nodes_from(to_remove)
    
    def assign_residue_type(self):
        net = self.net
        residues = self.model.get_residues()
        residues_dict = {}
        for residue in residues:
            res_type = residue.resname.strip()
            if len(res_type) < 3:
                res_type = aaconv.one2three(res_type)
            res_pos = residue.parent.id + str(residue.id[1])
            residues_dict[res_pos] = res_type
        nx.set_node_attributes(net, residues_dict, "residue_type")
    
    
    def assign_secondary_structure(self):

        pdb = self.path_to_file
        ppdb = PandasPdb().read_pdb(pdb)
        secondary_structure = {}
        helices_from_pdb = ppdb.df["OTHERS"][ppdb.df["OTHERS"]["record_name"] == "HELIX"][
        "entry"
        ]
        for helix in helices_from_pdb:
            identifier_h = helix[5:8].strip()
            initial_chain_h = helix[13].strip()
            initial_pos_h = helix[16:19].strip()
            final_pos_h = helix[28:31].strip()
            for i in range(int(initial_pos_h), int(final_pos_h) + 1):
                secondary_structure[initial_chain_h + str(i)] = (
                    "helix" + identifier_h + "-" + initial_chain_h
                )

        sheets_from_pdb = ppdb.df["OTHERS"][ppdb.df["OTHERS"]["record_name"] == "SHEET"][
            "entry"
        ]
        for sheet in sheets_from_pdb:
            identifier_s = sheet[6:8].strip()
            initial_chain_s = sheet[15].strip()
            initial_pos_s = sheet[17:20].strip()
            final_pos_s = sheet[28:31].strip()
            for i in range(int(initial_pos_s), int(final_pos_s) + 1):
                secondary_structure[initial_chain_s + str(i)] = (
                    "sheet" + identifier_s + "-" + initial_chain_s
                )
        
            net = self.net

            residues = list(net.nodes)  # assume they are ordered
            last_structure = None
            last_chain = None
            i = 0
            for residue in residues:
                chain = residue[0]
                try:
                    structure = secondary_structure[residue]
                    if structure != last_structure:
                        i += 1
                except KeyError:
                    if chain != last_chain:
                        i += 1
                    structure = "loop" + str(i)
                    secondary_structure[residue] = structure
                last_structure = structure
                last_chain = chain

            nx.set_node_attributes(net, secondary_structure, "secondary_structure")

    def assign_relation_links(self):
        net = self.net
        edges = net.edges
        secondary_structure = nx.get_node_attributes(self.net, "secondary_structure")
        relations = {}
        for u, v in edges:
            chain_u = u[0]
            chain_v = v[0]
            pos_u = u[1::]
            pos_v = v[1::]
            struct_u = secondary_structure[u]
            struct_v = secondary_structure[v]

            if chain_u == chain_v:
                dist = int(pos_v) - int(pos_u)
                if np.abs(dist) == 1:
                    relation = "1D"
                elif struct_u == struct_v:
                    if np.abs(dist) < 5:
                        relation = "2D"
                    else:
                        relation = "3D"
                else:
                    relation = "3D"
            else:
                relation = "4D"
            
            relations[(u, v)] = relation

        nx.set_edge_attributes(net, relations, "relation")

    def assign_local_partitioning(self):
        net = self.net
        nodes = net.nodes
        edges = net.edges
        relations = nx.get_edge_attributes(self.net, "relation")
        local_partitioning = {}

        for u in nodes:
            w = nx.degree(net, u, weight="weight")
            neighbors = nx.neighbors(net, u)
            separated_weights = {c: 0 for c in ["1D", "2D", "3D", "4D"]}
            for v in neighbors:
                try:
                    relation = relations[(u, v)]
                    link_weight = edges[(u, v)]["weight"]
                except KeyError:
                    relation = relations[(v, u)]
                    link_weight = edges[(v, u)]["weight"]
                separated_weights[relation] += link_weight
            local_partitioning[u] = [float(separated_weights[c]) / w for c in ["1D", "2D", "3D", "4D"]]
        
        nx.set_node_attributes(net, local_partitioning, "local_partitioning")

    def get_local_partitioning(self):
        return nx.get_node_attributes(self.net, "local_partitioning")
