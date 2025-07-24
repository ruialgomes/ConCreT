import os
from src.utils import open_fasta
import subprocess
from keras.utils import Sequence
import numpy as np


def get_prediction_classes(path):
    prediction_classes_list1 = []
    classes_list1 = []
    with open(path, "r") as classfile:
        for line in classfile:
            line2 = line.split(" , ")
            if line2[0] not in classes_list1:
                classes_list1.append(line2[0])
                class_list = (((line2[1].replace("\n", "")).replace("[", "")).replace("]", "")).split(", ")
                prediction_classes_list1.append((line2[0], class_list))
    return prediction_classes_list1


def get_singletons_test_identifiers_list(fasta_test_ident_list,blast_results_dic):
    singletons_test_identifiers_list = []
    for identifier in fasta_test_ident_list:
        if len(blast_results_dic[identifier]) == 0:
            singletons_test_identifiers_list.append(identifier)
    return singletons_test_identifiers_list


def get_tot_taxa_lvls(path):
    if not os.path.isdir(path):
        return 0
    subfolders = os.listdir(path)
    levels = []
    for subfol in subfolders:
        subfol_path = os.path.join(path, subfol)
        if os.path.isdir(subfol_path):
            levels.append(get_tot_taxa_lvls(subfol_path))
    if not levels:
        return 0
    return 1 + max(levels)


class BatchGenerator(Sequence):
    def __init__(self, data, batch_size=128):
        self.data = data
        self.batch_size = batch_size

    def __len__(self):
        return int(np.ceil(len(self.data) / self.batch_size))

    def __getitem__(self, idx):
        batch = self.data[idx * self.batch_size:(idx + 1) * self.batch_size]
        return batch


def parse_singletons(dir_path, singletons_taxa_dic, singletons_dic):
    singletondic = open_fasta(os.path.join(dir_path, "Singletons_taxa.fasta"))
    taxa_path = os.path.relpath(dir_path)
    with open(os.path.join(dir_path, "Singletons_taxa.txt"), "r") as singleton_taxa:
        for line in singleton_taxa:
            line = line.replace("\n", "")
            line_list = (line[1:]).split(" : ")
            new_ident = line_list[0]
            singletons_taxa_dic[new_ident] = taxa_path + "/" + line_list[1]
            singletons_dic[new_ident] = singletondic[line_list[0]]
    return singletons_taxa_dic, singletons_dic


def get_model_files(dir_path, model_info_dic, singletons_taxa_dic, singletons_dic):
    modelfile = ""
    classfile = ""
    dir_list = []
    fasta_list = []
    last_dir = (dir_path.split(os.path.sep))[-1]
    files_model_list = os.listdir(dir_path)
    for file in files_model_list:
        if os.path.isdir(os.path.join(dir_path, file)):
            dir_list.append(file)
        if file.split(".")[-1] == "h5":
            modelfile = file
        if file.split(".")[-1] == "csv":
            classfile = file
        if file.split(".")[-1] == "fasta":
            if file == "Singletons_taxa.fasta":
                singletons_taxa_dic, singletons_dic = parse_singletons(dir_path, singletons_taxa_dic, singletons_dic)
            else:
                fasta_list.append(file)

    model_info_dic[last_dir] = [modelfile, classfile, fasta_list, dir_list, dir_path]

    for curdir in dir_list:
        subfol_path = os.path.join(dir_path, curdir)
        model_info_dic, singletons_taxa_dic, singletons_dic = get_model_files(subfol_path, model_info_dic,
                                                                              singletons_taxa_dic, singletons_dic)

    return model_info_dic, singletons_taxa_dic, singletons_dic


def blast_sequences_list(sequence_list_path, genus_path, seq_result_file_path):
    blast_cl = ["blastn",
                "-query", sequence_list_path,
                "-subject", genus_path,
                "-outfmt", "6 qseqid sseqid pident evalue qcovhsp",
                "-out", seq_result_file_path
                ]

    process = subprocess.run(blast_cl, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

    results_list = []
    with open(seq_result_file_path, "r") as blast_file:
        for line in blast_file:
            line = line.replace("\n", "")
            line = line.split("\t")
            results_list.append(line)

    return results_list, process.stderr.decode()
