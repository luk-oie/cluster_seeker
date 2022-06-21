from Bio import SeqIO
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

with open("D:\Dev\cluster_seeker\ATG5_go", "rb") as f:
    d = pickle.load(f)

for i in d.features:
    print(i)

memes = os.listdir("D:\Dev\cluster_seeker\motifs")
print(memes)


GGAACCCGCTGAATGGCTGGGAAACTGTTTTTAATAAACAAGAGTTCAAACGCTGAGTCTCCATAGCCAACAGTCACGCCTTAAACTCTAAACTTTACATGAAAACTCTCCCAGTCCCTTTGAACTCTCCTTCCCTCTCCGACGGAGCATTCCAGTGTTTATGCATTTTTCGAGATTGCTGGCAGGATTGCGAGGCGCTTTGAATACTTTCCCTCTCTTTAGCAATCTCCTGCCTTCAGAACCCATTAATGCCCCACTCAGCAACATCAAGGTAATGCTTTGAAGTCCCTTCTGGCTACAGGCCTCTCACACTCTTTGGTTACTGTGCAATCAAATAATTAAATGGATTCTCAGGGAAAAAAAAAATCATTCTTTTCCTTGGTAATTAAATACTCCAGGTTGGACAGGAGCGCTGTGTTG
