# -*- coding: utf-8 -*-
import os
from collections import defaultdict
import argparse
import random
import string
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description='This script is to perform UniCluster for those cds without UniRef90 annotation.')
    #parser.add_argument("-o","--ouput_UniRefUnknown",
    #                    required=True,
    #                    type=str,
    #                    help="This is the filename that temporary store the UniRef90_unknown seqs, which will be used for MMSeq clustering.")
    parser.add_argument("-f","--folder_UniRefannotation",
                        required=True,
                        type=str,
                        help="The folder name containning all UniRef90 ananotations results.")
    args = parser.parse_args()
    return args

current_wd = os.getcwd()#get current work directory.

def generate_unique_random_strings(n, length=11):
    #"""generate n组不重复的length位随机字符串。"""
    characters = string.ascii_letters + string.digits
    unique_strings = set()
    while len(unique_strings) < n:
        new_string = ''.join(random.choice(characters) for _ in range(length))
        unique_strings.add(new_string)
    return list(unique_strings) #n 个不重复的长度为11的随机字符串（包括数字和大小写字母）。

def extract_unknown_Uniref(extracted_seqs,annotation_folder):   # store Uniref_unknown seqs to new file
    output = open(current_wd+"/"+extracted_seqs,"w")
    for root,dirs,files in os.walk(current_wd+"/"+annotation_folder):
        for file in files:
            if str(file).endswith("annotated"):
                f = open(os.path.join(root,file),"r")
                genelist = []
                for record in SeqIO.parse(f,"fasta"):
                    if "UniRef90_unknown" in record.description:
                        genelist.append(str(record.description).split("|")[0])
                faa = open(os.path.join(root,str(file).replace("annotated","translated")),"r")
                for record in SeqIO.parse(faa,"fasta"):
                    if str(record.description) in genelist:
                        output.write(">" + str(record.description).strip() + "\n")
                        output.write(str(record.seq) + "\n")
    output.close()


def parse_mmseq(MMseqresult):
    tmpin = open(current_wd+"/"+MMseqresult,"r")
    represent = [] # all clusters' representative gene id
    NewUnirefID = defaultdict(list) # representative gene id vs members
    for line in tmpin:
        cont = str(line).strip().split("\t")
        represent.append(cont[0])
        NewUnirefID[cont[0]].append(cont[1])
    Unirep = list(set(represent))
    GeneClusLib = {} # a dictionary: gene ids (used for clustering) vs UniCluster_****
    n = len(Unirep) # number of unique clusters
    ClusterIDs = generate_unique_random_strings(n)
    #print(ClusterIDs)
    i = 0
    for key,value in NewUnirefID.items():
        ID = ClusterIDs[i]
        i += 1
        for item in value:
            GeneClusLib[item] = "UniClust90_"+ID
    return GeneClusLib

def rename_ffn(MMseqresult,annotation_folder):
    GeneLib = parse_mmseq(MMseqresult)
    #print(GeneLib)
    for root,dirs,files in os.walk(current_wd+"/"+annotation_folder):
        for file in files:
            if str(file).endswith("annotated"):
                f = open(os.path.join(root,file),"r")
                fout = open(os.path.join(root,str(file)+"-UniCluster-updated.ffn"),"w")
                for record in SeqIO.parse(f,"fasta"):
                    if "UniRef90_unknown" not in str(record.description):
                        fout.write(">"+str(record.description).strip()+"\n"+str(record.seq).strip()+"\n")
                    else:
                        cont = str(record.description).strip().split("|")
                        geneid = cont[0].split(" ")[0]
                        try:
                            fout.write(">"+cont[0]+"|"+GeneLib[geneid]+"|"+cont[-1]+"\n"+str(record.seq).strip()+"\n")
                        except:
                            print(geneid," Not found in the Cluster dictionary.")
                fout.close()

def main():
    args = get_args()
    ouput_UniRefUnknown = "UniRef_unknown.fasta"
    extract_unknown_Uniref(ouput_UniRefUnknown,args.folder_UniRefannotation)
    os.system("mmseqs easy-cluster --cov-mode 0 -c 0.8 --min-seq-id 0.9 " + ouput_UniRefUnknown + " result tmp")
    rename_ffn("result_cluster.tsv", args.folder_UniRefannotation)

if __name__ == '__main__':
    main()