import os
import shutil
import pandas as pd
import pickle

class Feature:
    def __init__(self, name, type, site, strand):
        self.type = type
        self.site = site
        self.name = name
        self.strand = strand
    def __str__(self):
        return str(["Feature object:", self.name, self.type, self.site, self.strand])
    
class Genetic_Object:
    def __init__(self, sequence, annotation, features):
        self.sequence = sequence
        self.annotation = annotation
        self.features = features
    def __str__(self):
        return str(self.annotation) + " " + str(self.sequence)
    def add_feature(self, new_feature):
        self.features.append(new_feature)

vector_structure = []
for i in os.listdir("D:\Dev\cluster_seeker\motifs"):
    vector_structure.append(i[:-5])

vecs = os.listdir("D:\Dev\cluster_seeker\storage\\vectors")
chum = [[0,0,0,0,0,0,1,1,0], [0,0,0,0,0,0,0,1,0], [0,0,0,0,0,0,1,0,0]]
for i in vecs:
    with open(f"D:\Dev\cluster_seeker\storage\\vectors\{i}", "rb") as f:
        vectors = pickle.load(f)
        filtered = []
    for j in vectors:
        if j in chum:
            continue
        filtered.append(j)
    with open(f"D:\Dev\cluster_seeker\storage\\filtered_vecs\{i}", "wb") as f:
        pickle.dump(filtered, f)
    print(len(vectors))
    print(len(filtered))