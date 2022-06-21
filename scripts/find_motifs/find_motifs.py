from Bio import SeqIO
import os
import shutil
import pandas as pd
import pickle


def motifs_in_temp_fasta(meme):
    shutil.copy(meme, "D:\Dev\cluster_seeker\scripts\\find_motifs\\temp.meme")
    os.system('cmd /c "docker run -v D:\Dev\cluster_seeker\scripts\\find_motifs:/home/meme memesuite/memesuite fimo temp.meme temp.fasta"')

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

def child_of_genobj(parent, start, stop):
    child = Genetic_Object(parent.sequence[start:stop], parent.annotation + str([start, stop]), [])
    for i in parent.features:
        if start < i.site[0] and i.site[1] < stop:
            new_start = i.site[0] - start
            new_stop = i.site[1] - start
            child.features.append(Feature(i.name, i.type, [new_start, new_stop], i.strand))
    return child

def genobj_from_fasta(fasta):
    handle = open(fasta, "r")
    lines = handle.readlines()
    annotation = lines[0]
    sequence = "".join(lines[1:]).replace("\n", "")
    handle.close()
    return Genetic_Object(sequence, annotation, [])

def features_from_meme_output(df, type):
    features = []
    for i, j in df.iterrows():
        features.append(Feature(j["motif_alt_id"], type, [int(j["start"]), int(j["stop"])], j["strand"]))
    return features

def segments_from_genobj(gen_obj, core, radius):
    segments = []
    for i in gen_obj.features:
        segment = child_of_genobj(gen_obj, i.site[0], i.site[1])

    return segments

"""

SeqIO.convert(gb, "genbank", "temp.fasta", "fasta")
goi = genobj_from_fasta("temp.fasta")

motifs_in_temp_fasta(meme)
print(goi)
df = pd.read_csv("D:\Dev\cluster_seeker\scripts\\find_motifs\\fimo_out\\fimo.tsv", sep='\\t', header=0)
df.drop(df.tail(3).index,inplace=True)
#df = df.sort_values(by = "start")
print(df)
print(len(df))

#x = segments_from_df(df, goi.sequence, "ATG5", 500)
#print(x[0].features[0])
features = features_from_meme_output(df, "motif")

for i in features:
    goi.add_feature(i)

with open("D:\Dev\cluster_seeker\ATG5_go", "wb") as f:
    pickle.dump(goi, f)
"""

gene = "LEF1"
gb = f"D:\Dev\cluster_seeker\input_seqi\{gene}.gb"
meme = "D:\Dev\cluster_seeker\motifs\TCF7L2.meme"
memes = os.listdir("D:\Dev\cluster_seeker\motifs")

SeqIO.convert(gb, "genbank", "temp.fasta", "fasta")
goi = genobj_from_fasta("temp.fasta")

for i in memes:
    meme = f"D:\Dev\cluster_seeker\motifs\{i}"
    motifs_in_temp_fasta(meme)
    df = pd.read_csv("D:\Dev\cluster_seeker\scripts\\find_motifs\\fimo_out\\fimo.tsv", sep='\\t', header=0)
    df.drop(df.tail(3).index,inplace=True)
    if df.empty:
        continue
    df = df.sort_values(by = "start")
    features = features_from_meme_output(df, "motif")

    for j in features:
        goi.add_feature(j)


print(len(goi.features))
for i in goi.features:
    print(i)
with open(f"D:\Dev\cluster_seeker\storage\seqi_with_motifs_GO\{gene}", "wb") as f:
    pickle.dump(goi, f)

