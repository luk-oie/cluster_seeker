from sklearn import cluster
from sklearn.manifold import TSNE
#from keras.datasets import mnist
from sklearn.datasets import load_iris
from numpy import reshape
import seaborn as sns
import pandas as pd
import os
import pickle
import matplotlib.pyplot as plt
import mplcursors as mpc


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

def filter(chum, data, labels):
    new_data = []
    new_labels = []
    for i in range(len(data)):
        if data[i] in chum:
            continue
        new_data.append(data[i])
        new_labels.append(labels[i])
    
    return(new_data, new_labels)

vecs = os.listdir("D:\Dev\cluster_seeker\storage\\vectors")
segs = os.listdir("D:\Dev\cluster_seeker\storage\\segments")
dataset = []
labels = []
motifs = []
for i in vecs:
    with open(f"D:\Dev\cluster_seeker\storage\\vectors\{i}", "rb") as f:
        vectors = pickle.load(f)
    counter = 0
    for j in vectors:
        dataset.append(j)
        labels.append(i+" segment #" + str(counter))
        counter += 1

for i in segs:
    with open(f"D:\Dev\cluster_seeker\storage\\segments\{i}", "rb") as f:
        segments = pickle.load(f)
    counter = 0
    for j in segments:
        t = i+" segment #" + str(counter) + "   "
        for k in j.features:
            t += str(k.name + " ")
        motifs.append(t)
        counter += 1

print(dataset)
print(labels)
print(len(dataset), len(labels))

chum = [[0,0,0,0,0,0,1,1,0], [0,0,0,0,0,0,0,1,0], [0,0,0,0,0,0,1,0,0]]
fdata, flabels = filter(chum, dataset, motifs)

tsne = TSNE(n_components=2, verbose=1, random_state=123)
coordinates = tsne.fit_transform(fdata)

df = pd.DataFrame()

df["comp-1"] = coordinates[:,0]
df["comp-2"] = coordinates[:,1]
x = coordinates[:,0]
y = coordinates[:,1]


fig, ax = plt.subplots(1, figsize=(8,6))
sc = ax.scatter(x, y)
ax.scatter(x[272], y[272])
ax.scatter(x[160:184], y[160:184])
plt.xlabel("component-1")
plt.ylabel("component-2")

cursor = mpc.cursor(sc, hover = "True")

cursor.connect("add", lambda sel: sel.annotation.set_text(flabels[sel.index]))


plt.show()

hlyst = []
c = 0
for i in flabels:
    if "AXIN2" in i:
        hlyst.append(c)
    c+=1
print(hlyst)
