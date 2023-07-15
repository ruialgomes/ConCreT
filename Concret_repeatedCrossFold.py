# -*- coding: utf-8 -*-
"""Autotaxa_modelling_noextra.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1z4GdoCiDBF8DNhginPjLQrWb8TuumlmB

#Modelling

##Libraries
"""

import os
import time

import numpy as np
import math
import random

import multiprocessing


from sklearn.metrics import f1_score
from sklearn import metrics
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import MaxPooling2D
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Flatten
from tensorflow.keras.layers import Dropout

from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping

import keras

from sklearn.model_selection import RepeatedStratifiedKFold

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# ######################################################################################################################
# Functions
########################################################################################################################
########################################################################################################################
# write to an output file
def write_output(path, output):
    with open(path, "a+") as outfile:
        outfile.write(output)


########################################################################################################################
# create sequence in fasta format
def create_fasta_seq(seq_in):
    seq_out = ""
    for index_let, letter in enumerate(seq_in):
        seq_out += letter
        if index_let % 73 == 0 and index_let != 0:
            seq_out += "\n"
    return seq_out


########################################################################################################################
########################################################################################################################
def open_fasta(path):
    fasta_dic = {}
    fasta_seq = ""
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if len(fasta_dic) == 0:
                    seq_id = line[1:].replace("\n", "")
                    fasta_dic[seq_id] = ""
                    continue
                else:
                    fasta_dic[seq_id] = fasta_seq
                    fasta_seq = ""
                seq_id = line[1:].replace("\n", "")
                fasta_dic[seq_id] = fasta_seq
            else:
                fasta_seq += line.replace("\n", "")
        else:
            fasta_dic[seq_id] = fasta_seq
    return fasta_dic


#######################################################################################################################
#######################################################################################################################
#One-hot encode
def onehot_sequences(xDataset, yDataset):
    enconded_x = []
    nYdataset = []
    for indx1, sequence in enumerate(xDataset):

        if len(sequence) <= INPUT_LENGTH:
            nseq = repeat_to_length(sequence, INPUT_LENGTH)

        else:
            nseq = sequence[:INPUT_LENGTH]
            print(len(nseq))
            input("Maior 10k (-1?)")

        mapping = dict(zip("ATCG-KNSMYRWDVHB", range(16)))
        seq2 = [mapping[i] for i in nseq]
        seq_list = []
        for val in seq2:
            if val in (4, 5, 6, 7, 8, 9, 10, 11,12,13,14,15):
                seq_list.append(np.asmatrix(np.zeros(4)))
            else:
                seq_list.append(np.asmatrix(np.eye(4)[val]))
        seq_matrix = np.concatenate(tuple(seq_list))
        enconded_x.append(seq_matrix)
        nYdataset.append(yDataset[indx1])
        print('\rEncoding... %d%%' % (indx1 * (100 / len(xDataset))), end="")
    print("\n")
    return np.array(enconded_x), np.array(nYdataset)


########################################################################################################################
########################################################################################################################
def repeat_to_length(string_to_expand, length):
    return (string_to_expand * (int(length / len(string_to_expand)) + 1))[:length]


########################################################################################################################
########################################################################################################################
# enconde/decode one-hot classes
def onehot_classes(yDataset, mode):
    enconded_Y = []
    if mode == "e":
        for class_num in yDataset:
            file_class = []
            for index in range(TOT_FILES):
                if index == class_num:
                    file_class.append(1)
                else:
                    file_class.append(0)
            enconded_Y.append(np.array(file_class))

    elif mode == "d":
        for class_num in yDataset:
            enconded_Y.append(np.argmax(class_num))

    return np.array(enconded_Y)

# ######################################################################################################################
# ######################################################################################################################
def model_parse():


    class_file_dic = {}
    data_list = []
    data_ident_dic = {}

    global TOT_FILES
    TOT_FILES = len(files_dic)

    class_num = 0
    biggest_class_tot = 0
    for cur_class, cur_class_dic in files_dic.items():
        file_class = []
        print("-------------------")
        print(len(cur_class_dic))
        for index in range(TOT_FILES):
            if index == class_num:
                file_class.append(1)
            else:
                file_class.append(0)
        print(file_class)
        class_file_dic[cur_class] = file_class

        for ident1, seq1 in cur_class_dic.items():

            data_list.append([seq1, class_num])
            data_ident_dic[seq1] = ident1
        class_num += 1
        if biggest_class_tot < len(cur_class_dic):
            biggest_class_tot = len(cur_class_dic)


    category_names = list(files_dic.keys())
    categories_len = []
    for i1 in category_names:
      categories_len.append(len(files_dic[i1]))

    c=0
    for fileclass, listclass in class_file_dic.items():
        print(fileclass," - ",c," - ",str(listclass))
        write_output(dir_path + dir_path_name + "_classfile.csv", fileclass + " , " + str(listclass) + "\n")
        c+=1

    traintime = time.time()
    scores1 = []
    Fscores1 = []
    histories1 = []
    scores2 = []
    Fscores2 = []
    histories2 = []

    X_train1 = np.array([x for x, y in data_list])
    y_train1 = np.array([y for x, y in data_list])

    rskfold = RepeatedStratifiedKFold(n_splits=2,n_repeats=5, random_state=1)
    print("Total splits: ", rskfold.get_n_splits(X_train1, y_train1))

    # ##################################################################################################################
    print("Model 1 using 16 kernels in the first layer")
    for i, (train_ix, test_ix) in enumerate(rskfold.split(X_train1, y_train1)):
        fold_start_time = time.time()

        enconded_Y = onehot_classes(y_train1, "e")
        X_train, y_train, X_val, y_val = X_train1[train_ix], enconded_Y[train_ix], \
            X_train1[test_ix], enconded_Y[test_ix]

        aug_x_train = X_train
        aug_y_train = y_train

        print("Fold ", i)

        enconded_aug_x_train, aug_y_train = onehot_sequences(aug_x_train, aug_y_train)

        enconded_x_test, y_test2 = onehot_sequences(X_val, y_val)
        model1 = define_model1()

        custom_early_stopping = EarlyStopping(monitor="val_loss", patience=150, verbose=1, restore_best_weights=True)

        history = model1.fit(enconded_aug_x_train, aug_y_train, epochs=10000, batch_size=2,
                            validation_data=(enconded_x_test, y_test2),
                            verbose=1, callbacks=[custom_early_stopping])

        _, auc = model1.evaluate(enconded_x_test, y_test2, verbose=0)

        preds = model1.predict(enconded_x_test)

        decoded_y_pred = onehot_classes(preds, "d")
        decoded_y_test = onehot_classes(y_test2, "d")
        print("\n>Auc-PR:", auc)
        print(">F1 score: ", f1_score(decoded_y_test, decoded_y_pred, average='weighted'))
        print("Time took: " + str((time.time() - fold_start_time) / 60) + " min\n")

        scores1.append(auc)
        Fscores1.append(f1_score(decoded_y_test, decoded_y_pred, average='weighted'))
        histories1.append(history)

    # ##################################################################################################################
    print("Model 2, using 64 kernels in the first layer")
    for i, (train_ix, test_ix) in enumerate(rskfold.split(X_train1, y_train1)):
        fold_start_time = time.time()
        print("Fold ", i)
        enconded_Y = onehot_classes(y_train1, "e")
        X_train, y_train, X_val, y_val = X_train1[train_ix], enconded_Y[train_ix], \
            X_train1[test_ix], enconded_Y[test_ix]

        aug_x_train = X_train
        aug_y_train = y_train

        enconded_aug_x_train, aug_y_train = onehot_sequences(aug_x_train, aug_y_train)

        enconded_x_test, y_test2 = onehot_sequences(X_val, y_val)
        model2 = define_model2()

        custom_early_stopping = EarlyStopping(monitor="val_loss", patience=150, verbose=1, restore_best_weights=True)

        history = model2.fit(enconded_aug_x_train, aug_y_train, epochs=10000, batch_size=2,
                            validation_data=(enconded_x_test, y_test2),
                            verbose=1, callbacks=[custom_early_stopping])

        _, auc = model2.evaluate(enconded_x_test, y_test2, verbose=0)

        preds = model2.predict(enconded_x_test)

        decoded_y_pred = onehot_classes(preds, "d")
        decoded_y_test = onehot_classes(y_test2, "d")
        print("\n>Auc-PR:", auc)
        print(">F1 score: ", f1_score(decoded_y_test, decoded_y_pred, average='weighted'))
        print("Time took: " + str((time.time() - fold_start_time) / 60) + " min\n")

        scores2.append(auc)
        Fscores2.append(f1_score(decoded_y_test, decoded_y_pred, average='weighted'))
        histories2.append(history)

    # ##################################################################################################################
    print('F1 score: média=%.3f desvio=%.3f' % (np.mean(Fscores1) * 100, np.std(Fscores1) * 100))
    print(Fscores1)
    print('Auc: média=%.3f desvio=%.3f' % (np.mean(scores1) * 100, np.std(scores1) * 100))
    print(scores1)

    print('F1 score: média=%.3f desvio=%.3f' % (np.mean(Fscores2) * 100, np.std(Fscores2) * 100))
    print(Fscores2)
    print('Auc: média=%.3f desvio=%.3f' % (np.mean(scores2) * 100, np.std(scores2) * 100))
    print(scores2)

    print("Time took: " + str((time.time() - traintime) / 60) + " min\n")

    t_statistic, p_value = metrics.ttest_rel(Fscores1, Fscores2)
    print(f"The p value: {p_value}, using F1 score." )
    if p_value < 0.05:
        print("There are statistically significant differences between models.")
    else:
        print("There are  NO statistically significant differences between models.")

    t_statistic, p_value = metrics.ttest_rel(scores1, scores2)
    print(f"The p value: {p_value}, using auc-pr score." )
    if p_value < 0.05:
        print("There are statistically significant differences between models.")
    else:
        print("There are  NO statistically significant differences between models.")


# ######################################################################################################################
# ######################################################################################################################
def define_model1():
    model = Sequential()
    model.add(Conv2D(16, (150, 4), activation='relu', padding='same', input_shape=(INPUT_LENGTH, INPUT_HEIGHT, 1)))
    model.add(MaxPooling2D((50, 1), padding='same'))
    model.add(Dropout(0.5))
    model.add(Conv2D(32, (50, 4), activation='relu', padding='same'))
    model.add(MaxPooling2D((25, 1), padding='same'))
    model.add(Dropout(0.5))

    model.add(Flatten())
    model.add(Dense(64, activation='relu'))
    model.add(Dropout(0.4))
    model.add(Dense(TOT_FILES, activation='softmax'))

    # compilando modelo
    opt = Adam(learning_rate=0.0001)

    if TOT_FILES == 2:
        METRICS = [keras.metrics.AUC(name='auc', curve="PR")]
        model.compile(optimizer=opt, loss='binary_crossentropy', metrics=METRICS)
    else:
        METRICS = [keras.metrics.AUC(name='auc', curve="PR",multi_label = True)]
        model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=METRICS)
    print(model.summary())
    return model

# ######################################################################################################################
# ######################################################################################################################
def define_model2():
    model = Sequential()
    model.add(Conv2D(64, (150, 4), activation='relu', padding='same', input_shape=(INPUT_LENGTH, INPUT_HEIGHT, 1)))
    model.add(MaxPooling2D((50, 1), padding='same'))
    model.add(Dropout(0.5))
    model.add(Conv2D(32, (50, 4), activation='relu', padding='same'))
    model.add(MaxPooling2D((25, 1), padding='same'))
    model.add(Dropout(0.5))


    model.add(Flatten())
    model.add(Dense(64, activation='relu'))
    model.add(Dropout(0.4))
    model.add(Dense(TOT_FILES, activation='softmax'))

    # compilando modelo
    opt = Adam(learning_rate=0.0001)

    if TOT_FILES == 2:
        METRICS = [keras.metrics.AUC(name='auc', curve="PR")]
        model.compile(optimizer=opt, loss='binary_crossentropy', metrics=METRICS)
    else:
        METRICS = [keras.metrics.AUC(name='auc', curve="PR",multi_label = True)]
        model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=METRICS)
    print(model.summary())
    return model

########################################################################################################################
########################################################################################################################
# vars

BATCH_SIZE = multiprocessing.cpu_count() - 1

INPUT_LENGTH = 15000
INPUT_HEIGHT = 4

stime = time.time()
########################################################################################################################
dir_path = "DataVMR/Data/model_data/Cressdnaviricota/"  #model

dir_path_name = dir_path.split("/")[-2]

files_dic = {}
max_len = 0
maxlendic = {}

for file in os.listdir(dir_path):
    if os.path.isdir(dir_path+file) or (dir_path+file).split(".")[-1] != "fasta" or file == "Singletons_taxa.fasta":
        continue
    seqs_dic = open_fasta(dir_path+file)

    for ident,seq in seqs_dic.items():
        if len(seq) > max_len:
            maxlendic[ident] = seq
            max_len = len(seq)
    files_dic[file.split(".")[0]] = seqs_dic

resp = model_parse()