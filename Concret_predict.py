import os
import statistics

import time
import subprocess
import tensorflow as tf
import numpy as np
import argparse

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# ######################################################################################################################
# Functions
########################################################################################################################
########################################################################################################################
# Auxiliar functions
########################################################################################################################
########################################################################################################################
# write to an output file
def write_output(path, output):
    with open(path, "a+") as outfile:
        outfile.write(output)


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


########################################################################################################################
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
def encode_sequence(seq):
    nseq_list = []
    if len(seq) <= INPUT_LENGTH:
        nseq = repeat_to_length(seq, INPUT_LENGTH)
        nseq_list.append(nseq)
    else:
        nseq = seq[:INPUT_LENGTH]
        nseq_list.append(nseq)
        print("Sequence bigger than the maximum length, is it Cressdnaviricota? Removing excess...")
        write_output(logpath, "Sequence bigger than the maximum length, is it Cressdnaviricota? Removing excess...\n")

    # for seq3 in nseq_list:
    mapping = dict(zip("ATCG-KNSMYRWD", range(13)))
    seq2 = [mapping[i] for i in nseq]
    seq_list = []
    for val in seq2:
        if val in (4, 5, 6, 7, 8, 9, 10, 11,12):
            seq_list.append(np.asmatrix(np.zeros(4)))
        else:
            seq_list.append(np.asmatrix(np.eye(4)[val]))
    sequence_matrix = np.concatenate(tuple(seq_list))

    return sequence_matrix
########################################################################################################################
########################################################################################################################
def repeat_to_length(string_to_expand, length):
    return (string_to_expand * (int(length / len(string_to_expand)) + 1))[:length]

########################################################################################################################
########################################################################################################################
def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

########################################################################################################################
########################################################################################################################
def has_lower(inputString):
    return any(char.islower() for char in inputString)

########################################################################################################################
########################################################################################################################
def has_letter(inputString):
    return any(char.isalpha() for char in inputString)


########################################################################################################################
########################################################################################################################
def get_gpu_memory_usage():
    result = subprocess.run(['nvidia-smi', '--query-gpu=memory.used', '--format=csv,nounits,noheader'],
                            capture_output=True, text=True)
    output = result.stdout.strip().split('\n')
    gpu_memory_usage = [int(line) for line in output]
    return gpu_memory_usage


########################################################################################################################
########################################################################################################################
## Main functions

#######################################################################################################################
########################################################################################################################
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


########################################################################################################################
########################################################################################################################
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


########################################################################################################################
########################################################################################################################
def get_model_files(dir_path,last_dir):
    model_file, class_file, dir_list,fasta_list = get_files_names(dir_path)
    #for genuses with not directories for species, saves the fasta files for BLAST
    if dir_list == []:
        MODEL_INFO_DIC[last_dir] = [model_file, class_file, fasta_list, dir_path]
    else:
        MODEL_INFO_DIC[last_dir] = [model_file, class_file, dir_list, dir_path]

    for curdir in dir_list:
        subfol_path = os.path.join(dir_path, curdir)
        get_model_files(subfol_path,curdir)


########################################################################################################################
########################################################################################################################
def get_files_names(cur_path):
    modelfile = ""
    classfile = ""
    dir_list = []
    fasta_list = []
    files_model_list = os.listdir(cur_path)
    for file in files_model_list:
        if os.path.isdir(os.path.join(cur_path,file)):
            dir_list.append(file)
        if file.split(".")[-1] == "h5":
            modelfile = file
        if file.split(".")[-1] == "csv":
            classfile = file
        if file.split(".")[-1] == "fasta":
            #separe singletons
            if file == "Singletons_taxa.fasta": # save acc taxa in ident
                singletondic = open_fasta(os.path.join(cur_path,"Singletons_taxa.fasta"))
                taxa_path = cur_path.split(firstTaxa_model_name)[1]
                with open(os.path.join(cur_path,"Singletons_taxa.txt"),"r") as singleton_taxa:
                    for line in singleton_taxa:
                        line_list = (line[1:]).split(" : ")  #by standard: '>MN100000.1_Opuntia_virus_1 : Opunvirus'
                        new_ident = line_list[0]
                        SINGLETONS_TAXA_DIC[new_ident]=taxa_path+"/"+line_list[1]
                        SINGLETONS_DIC[new_ident] = singletondic[line_list[0]]
                continue
            fasta_list.append(file)
    return modelfile,classfile,dir_list,fasta_list


########################################################################################################################
########################################################################################################################
def get_blast_results(cur_predictions_list1,cur_predictions_values_dic1,cur_type_list1):

    ############################################################################################################
    # BLAST for singletons and each genus in the predicted family

    singletons_blast_resp_list = blast_Singletons_predict(cur_sequence)

    #if there is BLAST results for the singletons, uses then as the current best result
    if singletons_blast_resp_list != []:
        best_allgenus_ident_info = (singletons_blast_resp_list[0][0],
                 [float(singletons_blast_resp_list[0][1]),float(singletons_blast_resp_list[0][2])],singletons_blast_resp_list)
    else:
        best_allgenus_ident_info = ("", [0, 0], [])

    for genus_fasta_path in fasta_list:
        ident_list = blast_genus(cur_sequence, genus_fasta_path, dir_path)
        # save higher identity results
        if ident_list != []:

            if statistics.mean([float(ident_list[0][1]),float(ident_list[0][2])]) > statistics.mean(best_allgenus_ident_info[1]):
                best_allgenus_ident_info = ((genus_fasta_path.split("."))[0], [float(ident_list[0][1]),float(ident_list[0][2])], ident_list)

    ############################################################################################################
    # parse blast results
    final_blast_results = []
    # if there is a blast result
    if best_allgenus_ident_info[0] != "" :
        # check if the genus with best blast result match the predicted genus
        if best_allgenus_ident_info[0] != cur_predictions_list1[-1]:
            #when don't match, change the predictions, checking if is a singleton or other genus in the family
            #first for when exists singletons results
            if singletons_blast_resp_list and best_allgenus_ident_info[0] == singletons_blast_resp_list[0][0]:

                singleton_taxa_info_list = SINGLETONS_TAXA_DIC[singletons_blast_resp_list[0][0]].replace("\n","").split("/")
                singleton_taxa_info_list = singleton_taxa_info_list[1:]
                new_pred_type = ["Singleton" for _ in singleton_taxa_info_list]

                new_predValue_dic = {v1: ([[1.00]],[]) for i1,v1 in enumerate(singleton_taxa_info_list)
                                     if singleton_taxa_info_list[i1] != singleton_taxa_info_list[-1]}

                new_predValue_dic[singleton_taxa_info_list[-1]] = ([[float(best_allgenus_ident_info[1][0]) / 100]],[])

                cur_predictions_list1 = singleton_taxa_info_list.copy()
                cur_type_list1 = new_pred_type.copy()
                cur_predictions_values_dic1 = new_predValue_dic.copy()
                final_blast_results.append("Accession with sequence similar to singleton!")

            else:
                genus_pred = cur_predictions_list1.pop()
                cur_predictions_list1.append(best_allgenus_ident_info[0])

                pred_categories = list(cur_predictions_values_dic1[genus_pred][1])
                del cur_predictions_values_dic1[genus_pred]
                cur_predictions_values_dic1[best_allgenus_ident_info[0]] = ([[float(best_allgenus_ident_info[1][0]) / 100]],pred_categories)

                cur_type_list1.pop()
                cur_type_list1.append("Blast")
                final_blast_results.append("Accession with sequence similar to other genus!")
        else:
            final_blast_results.append("Accession with sequence similar to the predicted genus!")

        # get blast results info to the final output
        final_blast_results.extend(best_allgenus_ident_info[2])

    else:
        final_blast_results.append("Accession with NO sequence similarity using blast to the predicted family genuses members!")

    return final_blast_results,cur_predictions_list1, cur_predictions_values_dic1,cur_type_list1


########################################################################################################################
########################################################################################################################
def blast_Singletons_predict(sequence):
    seq_fasta_file_path = os.path.join(files_path,"seq_file.fasta")
    acc_id_forPath = cur_seq_identifier.replace("/","_")
    seq_result_file_path = os.path.join(files_path,"blast_results_Singletons_" + acc_id_forPath + ".txt")
    with open(seq_fasta_file_path, "w") as f:
        f.write(">" + cur_seq_identifier + "\n" + sequence + "\n")

    ################################################################################################################
    process = subprocess.Popen(["blastn", "-query", seq_fasta_file_path, "-subject", SINGLETONS_FASTA_PATH,
                                "-outfmt", "6 qseqid sseqid pident evalue qcovhsp", "-evalue", "1000",
                                "-out", seq_result_file_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    if stderrdata.decode() != "":
        print("Run Blast Error: \n")
        print(stderrdata.decode())

    ################################################################################################################
    resp_file_blast = open(seq_result_file_path, "r")
    resp_list = []
    for line in resp_file_blast:
        line = line.replace("\n", "")
        line = line.split("\t")
        resp_list.append((line[1],line[2],line[4]))
    resp_file_blast.close()

    os.remove(seq_fasta_file_path)
    if resp_list == []:
        os.remove(seq_result_file_path)
    return resp_list


########################################################################################################################
########################################################################################################################
def blast_genus(sequence,genus_path,dir_path1):
    seq_fasta_file_path = os.path.join(files_path, "seq_file.fasta")
    acc_id_forPath = cur_seq_identifier.replace("/", "_")
    genusname = genus_path.split(".")[0]
    seq_result_file_path = os.path.join(files_path, "blast_results_" + cur_seq_identifier + "_gen_"
                                        + genusname + ".txt")

    with open(seq_fasta_file_path, "w") as f:
        f.write(">" + cur_seq_identifier + "\n" + sequence + "\n")

    ################################################################################################################
    process = subprocess.Popen(["blastn", "-query", seq_fasta_file_path, "-subject",
                                os.path.join(dir_path1 , genus_path), "-outfmt",
                                "6 qseqid sseqid pident evalue qcovhsp",
                                "-out", seq_result_file_path],    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdoutdata, stderrdata = process.communicate()
    if stderrdata.decode() != "":
        print("Blast Error: \n")
        print(stderrdata.decode())
    ################################################################################################################
    resp_file_blast = open(seq_result_file_path, "r")
    resp_list = []
    for line in resp_file_blast:
        line = line.replace("\n","")
        line = line.split("\t")
        resp_list.append((line[1],line[2],line[4]))
    resp_file_blast.close()

    os.remove(seq_fasta_file_path)
    if resp_list == []:
        os.remove(seq_result_file_path)
    return resp_list

########################################################################################################################
########################################################################################################################

inputfile = "DataVMR/Data/test_data/classes_fasta/Genus/all_genus.fasta"

model_dir_path = "DataVMR/Data/model_data/Cressdnaviricota_borderAugn64_phylnoaugn"

output_path = "DataVMR/Data/test_data/output_predictions"
# output_path = os.getcwd()

parser = argparse.ArgumentParser()

parser.add_argument('--in', dest='inputfile', default=inputfile, help='The input file or path to file. Default is '
                                                                      'the example file.')
parser.add_argument('--model', dest='modelpath', default=model_dir_path, help='A path to the model directory. '
                                                                           'Default is the Cressdnaviricota model.')
parser.add_argument('--out', dest='outputpath', default=output_path, help='A path for the output directory, where a '
                        'output folder will be created. Default is the directory where the script is being executed.')
parser.add_argument('--output-name', dest='outputname', default="Results-name-default", help='A name for the result '
                                                                             'file. Default is the input file name.')
parser.add_argument('--max-blast-hits', dest='max_blast_hits', type=int, default=3, help='Maximum of BLAST hits that '
                                                                    'will be printed in the output. Default is 3.')

args = parser.parse_args()

init_test_file_path = args.inputfile
model_dir_path = args.modelpath
output_path = args.outputpath
output_name = args.outputname
BLAST_GENUS_MAX_HITS = args.max_blast_hits

if output_name == "Results-name-default":
    directory, file_name_aux = os.path.split(init_test_file_path)
    filename = (file_name_aux.split("."))[0]
else:
    filename = output_name

########################################################################################################################
#output
#dir and files for output information


filename_path = os.path.join(output_path ,filename+"_Results")
os.makedirs(filename_path, exist_ok=True)
files_path = os.path.join(filename_path,"files")
os.makedirs(files_path, exist_ok=True)

#Clean if files already exists, no restart possibility for now.
#cleaning main csv file
predpath = os.path.join(filename_path , filename + "_predictions.csv")
if os.path.exists(predpath):
    os.remove(predpath)

#cleaning secondary csv with full information
fullinfopath = os.path.join(filename_path , filename + "_FullInfo.csv")
if os.path.exists(fullinfopath):
    os.remove(fullinfopath)

#cleaning file for occurrences during prediction
logpath = os.path.join(filename_path , filename + ".log")
if os.path.exists(logpath):
    os.remove(logpath)

#creating files
write_output(logpath,"Starting...\n")
print("Starting...")
write_output(predpath,"Identifier,Length,Predidction,Pred. conf.(%),Pred. result,Blast Ids,Blast(ident/cov)\n")
write_output(fullinfopath,"Identifier,Class (pred), pred(type/perc),"
                                                 "Order(pred), pred(type/perc),Family(pred), pred(type/perc),"
                                                 "Genus(pred), pred(type/perc),Blast Ids,""Blast(ident/cov)\n")

########################################################################################################################
#Load model information

firstTaxa_model_name = (model_dir_path.split(os.path.sep))[-1]

taxa_lvl_dic = {0:"Genus",1:"Family",2:"Order",3:"Class",4:"Phylum",5:"Kingdom",6:"Realm"}
total_taxa_levels = get_tot_taxa_lvls(model_dir_path)

print("Classifing from the '"+taxa_lvl_dic[total_taxa_levels]+ "' Taxonomic level.")
write_output(logpath,"Classifing from the "+taxa_lvl_dic[total_taxa_levels]+ " Taxonomic level.\n")

MODEL_INFO_DIC = {}
SINGLETONS_DIC = {}
SINGLETONS_TAXA_DIC = {}
get_model_files(model_dir_path, firstTaxa_model_name)
########################################################################################################################
#create singletons file for BLAST
SINGLETONS_FASTA_PATH = os.path.join(files_path, "singletons_len" + str(len(SINGLETONS_DIC)) + ".fasta")
if not os.path.exists(SINGLETONS_FASTA_PATH):
    with open(SINGLETONS_FASTA_PATH, "w") as f2:
        for identifier, sequence2 in SINGLETONS_DIC.items():
            idn2 = identifier.replace("\n", "")
            seq2 = create_fasta_seq(sequence2.replace("\n", ""))
            f2.write(">" + idn2 + "\n" + create_fasta_seq(sequence2.replace("\n", "")) + "\n")

########################################################################################################################
#loading the models files with TensorFlow
write_output(logpath,"Loading models...\n")
print("Loading models files...")
load_time = time.time()
MODELS_MAIN_DIC = {}
#Allow gpu growth while opening multiple models in tensorflow
tf.config.experimental.set_memory_growth(tf.config.list_physical_devices('GPU')[0], True)

last_gpu_usage = get_gpu_memory_usage()

disk_space_list = []
for taxa,infolist in MODEL_INFO_DIC.items():
    if infolist[0]!="":
        write_output(logpath,"Loading model: "+ str(infolist[0])+"\n")
        print("Loading model: ", infolist[0],"\n")
        cur_model_file = tf.keras.models.load_model(os.path.join(infolist[3], infolist[0]))
        MODELS_MAIN_DIC[infolist[0]] = cur_model_file
        disk_space_usage = round((os.path.getsize(os.path.join(infolist[3], infolist[0])))/ (1000 * 1000),1)
        disk_space_list.append(disk_space_usage)

cur_gpu_memory_usage = get_gpu_memory_usage()
write_output(logpath,"Models loaded in "+ str(round(time.time() - load_time)) +
             " seconds. Models Disk Size: "+str(sum(disk_space_list))+" MBs.\nApproximated Used GPU memory: "+
             str(cur_gpu_memory_usage[0]-last_gpu_usage[0])+" MBs.\n")
print("Models loaded in "+ str(round(time.time() - load_time)) + " seconds. Models Disk Size: "+str(sum(disk_space_list))+
      " MBs.")
print("Approximated Used GPU memory: "+str(cur_gpu_memory_usage[0]-last_gpu_usage[0])+" MBs.")
#using last loaded model to get input sizes
layers = cur_model_file.layers
first_layer = layers[0]
INPUT_LENGTH = (first_layer.input_shape)[1]

########################################################################################################################
#start parsing

fasta_test_dic = open_fasta(init_test_file_path)

count_tot_seqs = 0
file_time = time.time()
med_time_list = []
write_output(logpath, "Parsing...\n")
print("Parsing...\n###")
for cur_seq_identifier,cur_sequence in fasta_test_dic.items():
    #Remove chars that cant be used in files names
    cur_seq_identifier = ((cur_seq_identifier.replace("\n","")).replace(" ","_")).replace("\t","")
    idtime = time.time()
    count_tot_seqs+=1

    write_output(logpath, "###\n###\nACCESSION: " + cur_seq_identifier + "\n")
    print("ACCESSION: " + cur_seq_identifier)

    ################################################################################################################
    ################################################################################################################
    # Model prediction
    cur_predictions_list = []
    cur_type_list = []
    cur_predictions_values_dic = {}

    seq_matrix = encode_sequence(cur_sequence)

    taxa_ops_count = total_taxa_levels
    last_taxa = firstTaxa_model_name
    while taxa_ops_count >= 0:
        taxa_ops_count -=1
        model_file = MODEL_INFO_DIC[last_taxa][0]
        categories_file = MODEL_INFO_DIC[last_taxa][1]
        fasta_list = MODEL_INFO_DIC[last_taxa][2]
        dir_path = MODEL_INFO_DIC[last_taxa][3]

        if model_file == "":

            single_name = ((MODEL_INFO_DIC[cur_predictions_list[-1]][2][0]).split(".f"))[0]

            cur_predictions_list.append(single_name)
            cur_predictions_values_dic[single_name] = ([[1.00]], [(single_name,['1'])])
            cur_type_list.append("Single")
            last_taxa = single_name

        else:
            tf.keras.backend.clear_session()
            prediction_classes_list = get_prediction_classes(os.path.join(dir_path, categories_file))
            model = MODELS_MAIN_DIC[model_file]

            predict1 = model.predict(np.array([seq_matrix]), verbose=0)
            sorted_indices = (np.argsort(predict1))[0][::-1]

            pre_order = []
            for si in sorted_indices:
                for category_name, category_list in prediction_classes_list:
                    if category_list[si] == "1":
                        pre_order.append(category_name)

            last_taxa = pre_order[0] #get the highest prediction probability category
            cur_predictions_list.append(pre_order[0])
            cur_predictions_values_dic[pre_order[0]] = (predict1, prediction_classes_list)
            cur_type_list.append("Model")

    ################################################################################################################
    ################################################################################################################
    #blast against all singletons and all genus in predicted family

    genus_blast_results,cur_predictions_list,cur_predictions_values_dic,cur_type_list = \
        get_blast_results(cur_predictions_list,cur_predictions_values_dic,cur_type_list)

    ####################################################################################################################
    ####################################################################################################################
    # Format and write prediction, create csv lines, with ',' for separation of columns
    main_csv_line = cur_seq_identifier + "," + str(len(cur_sequence)) + ","
    fullInfo_csv_line = cur_seq_identifier + ","

    pred_output_list = []
    last_taxa2 = ""
    for result_idx,taxa in enumerate(cur_predictions_list):
        pred_perc_idx = int(np.argmax(np.array(cur_predictions_values_dic[taxa][0])))
        pred_perc = round((cur_predictions_values_dic[taxa][0][0][pred_perc_idx]) * 100, 2)
        last_taxa2 = taxa
        pred_output_list.append(taxa)

        if cur_type_list[result_idx] == "Singleton":
            pred_output_list.append(cur_type_list[result_idx] + "/" + str(genus_blast_results[1][1]) + "%")
        else:
            pred_output_list.append(cur_type_list[result_idx] + "/" + str(pred_perc) + "%")

    write_output(logpath, "Pred: " + str(last_taxa2) + "\n")
    print("Pred: " + str(last_taxa2))

    for k,v in cur_predictions_values_dic.items():
        write_output(logpath, "Taxa: " + str(k) +"; Values: "+str(v)+ "\n")

    for joined_taxa in pred_output_list:
        fullInfo_csv_line += joined_taxa + ","

    main_csv_line += pred_output_list[-2]+","
    main_csv_line += pred_output_list[-1]+","
    main_csv_line += genus_blast_results[0]+","


    if len(genus_blast_results) > 1:
        main_csv_line += genus_blast_results[1][0] + "," + genus_blast_results[1][1] + "/" + genus_blast_results[1][2]
        for idx, blast_results_set in enumerate(genus_blast_results):
            if idx == 0:
                continue
            fullInfo_csv_line += blast_results_set[0] + "," + blast_results_set[1] + "/" + blast_results_set[2] + ","

    write_output(predpath, main_csv_line+"\n")
    write_output(fullinfopath, fullInfo_csv_line + "\n")

    print("Total id time: " + str((time.time() - idtime)) + " s;")
    write_output(logpath,"Total id time " + str((time.time() - idtime)) + " s;\n")
    print("Total time: " + str((time.time() - file_time) / 60) + " min;")
    write_output(logpath,"Total time: " + str((time.time() - file_time) / 60) + " min;\n")
    print("current: ", count_tot_seqs, ". Total: ", len(fasta_test_dic), "\n----------------------")

print("tot accs: ",len(fasta_test_dic))
write_output(logpath,"tot accs: "+str(len(fasta_test_dic))+ "\n")
print("Total time: " + str((time.time() - file_time) / 60) + " min;")










