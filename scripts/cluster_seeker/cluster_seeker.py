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

def extend_gen_obj(extendee, extension, direction):
    if direction == "+":
        extendee.sequence = extendee.sequence + extension.sequence
        extendee.annotation = "Merged objects: " + extendee.annotation + extension.annotation
        extendee.features = extendee.features + extension.features
    if direction == "-":
        extendee.sequence = extension.sequence + extendee.sequence
        extendee.annotation = "Merged ojects: " + extension.annotation + extendee.annotation
        extendee.features = extension.features + extendee.features
    return extendee

def segments_from_genobj(gen_obj, core, radius):
    segments = []
    prev_core_pos = -2*radius
    for i in gen_obj.features:
        if i.name != core:
            continue
        if i.site[0] < prev_core_pos + 2*radius:
            print(i.name, i)
            extension = child_of_genobj(gen_obj, prev_core_pos + radius, i.site[0] + radius)
            extend_gen_obj(segments[-1], extension, "+")
            prev_core_pos = i.site[1]
            continue

        seg_start = i.site[0] - radius
        #we gonna have a problem here when we encounter a case with i.site[0] < radius
        seg_end = i.site[1] + radius 
        segment = child_of_genobj(gen_obj, seg_start, seg_end)
        #add annotation - segment X-Y (X - core, Y - number)
        #may encounter issues with cut off features at the borders of child - check behaviour of child_of_genobj()
        segments.append(segment)
        prev_core_pos = i.site[1]
    return segments



def segments_from_genobj_new(gen_obj, core, radius):
    segments = []
    prev_core_pos = -2*radius
    for i in gen_obj.features:
        if i.name not in core:
            continue
        if i.site[0] < prev_core_pos + 2*radius:
            print(i.name, i)
            extension = child_of_genobj(gen_obj, prev_core_pos + radius, i.site[0] + radius)
            extend_gen_obj(segments[-1], extension, "+")
            prev_core_pos = i.site[1]
            continue

        seg_start = i.site[0] - radius
        #we gonna have a problem here when we encounter a case with i.site[0] < radius
        seg_end = i.site[1] + radius 
        segment = child_of_genobj(gen_obj, seg_start, seg_end)
        #add annotation - segment X-Y (X - core, Y - number)
        #may encounter issues with cut off features at the borders of child - check behaviour of child_of_genobj()
        segments.append(segment)
        prev_core_pos = i.site[1]
    return segments




core = ["TCF7L2","TCF7L1"]
gen_objs = os.listdir("D:\Dev\cluster_seeker\storage\seqi_with_motifs_GO")
for i in gen_objs:
    gen_obj = i
    with open(f"D:\Dev\cluster_seeker\storage\seqi_with_motifs_GO\{gen_obj}", "rb") as f:
        go = pickle.load(f)
    segs = segments_from_genobj_new(go, core, 250)
    with open(f"D:\Dev\cluster_seeker\storage\segments\{gen_obj}", "wb") as f:
        pickle.dump(segs, f)

"""
gen_obj = "AXIN2"
core = "TCF7L2"
with open(f"D:\Dev\cluster_seeker\storage\seqi_with_motifs_GO\{gen_obj}", "rb") as f:
    go = pickle.load(f)

segs = segments_from_genobj(go, core, 250)
print(segs)

"""