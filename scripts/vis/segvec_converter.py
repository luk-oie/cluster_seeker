from sklearn.manifold import TSNE
#from keras.datasets import mnist
from sklearn.datasets import load_iris
from numpy import reshape
import seaborn
import pandas
import os
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


def segvec_converter(segment, vector_structure):
    vector = []
    for i in vector_structure:
        counter = 0
        for j in segment.features:
            if j.name == i:
                counter += 1
        vector.append(counter)
    return vector

segs = os.listdir("D:\Dev\cluster_seeker\storage\segments")
vector_structure = []
for i in os.listdir("D:\Dev\cluster_seeker\motifs"):
    vector_structure.append(i[:-5])

for i in segs:
    with open(f"D:\Dev\cluster_seeker\storage\segments\{i}", "rb") as f:
        segments_of_i = pickle.load(f)
        vectors = []
    for j in segments_of_i:
        vector = segvec_converter(j, vector_structure)
        vectors.append(vector)
    with open(f"D:\Dev\cluster_seeker\storage\\vectors\{i}", "wb") as f:
        pickle.dump(vectors, f)


