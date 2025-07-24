import time
import os
import numpy as np
import tensorflow as tf
from src.utils import encode_sequence, write_output, open_fasta, create_fasta_seq
from src.parsing_utils import get_prediction_classes, get_tot_taxa_lvls, get_model_files, blast_sequences_list, \
    BatchGenerator


def parse_fasta_file(
        init_test_file_path,
        model_path,
        output_path,
        output_name,
        save_files,
        use_gpu,
        blast_genus_max_hits,
        batch_size
):
    # create paths and folders
    if output_name == "Results-name-default":
        directory, file_name_aux = os.path.split(init_test_file_path)
        filename = (file_name_aux.split("."))[0]
    else:
        filename = output_name

    filename_dir_path = str(os.path.join(output_path, filename + "_Results"))
    os.makedirs(filename_dir_path, exist_ok=True)

    files_dir_path = os.path.join(filename_dir_path, "files")
    os.makedirs(files_dir_path, exist_ok=True)

    # Clean if files already exists
    predpath = os.path.join(filename_dir_path, filename + "_predictions.csv")
    if os.path.exists(predpath):
        os.remove(predpath)

    fullinfopath = os.path.join(filename_dir_path, filename + "_FullInfo.csv")
    if os.path.exists(fullinfopath):
        os.remove(fullinfopath)

    logpath = os.path.join(filename_dir_path, filename + ".log")
    if os.path.exists(logpath):
        os.remove(logpath)

    tottime = time.time()
    write_output(logpath, "Starting...\n")
    write_output(predpath,
                 "Identifier,Length,Predidction,Pred. conf.(%),Pred. result,Blast Ids,Blast(ident/cov)\n")
    write_output(fullinfopath,
                 "Identifier,Class (pred), pred(type/perc),"
                 "Order(pred), pred(type/perc),Family(pred), pred(type/perc),"
                 "Genus(pred), pred(type/perc),Blast Ids,""Blast(ident/cov)\n")

    model_info_dic = {}
    singletons_taxa_dic = {}
    singletons_dic = {}

    # # load model saved information
    model_info_dic, singletons_taxa_dic, singletons_dic, = get_model_files(model_path, model_info_dic,
                                                                           singletons_taxa_dic, singletons_dic)
    # open the input fasta file in a dic
    fasta_test_dic = open_fasta(init_test_file_path)
    fasta_test_ident_list = list(fasta_test_dic.keys())

    #get prediction with the model
    write_output(logpath, "Getting model predictions...\n")
    pred_time = time.time()
    model_predict_results_dic, taxa_prediction_idents_dic = model_predict_fasta_dic(
        fasta_test_dic,
        model_path,
        model_info_dic,
        logpath,
        batch_size,
        use_gpu
    )
    write_output(logpath, f"Time took for predictions: {round((time.time() - pred_time) / 60, 2)} minutes.\n")

    #blast against each predicted genus
    write_output(logpath, f"Running blast against each predicted genus.\n")
    blastgenus_time = time.time()
    blast_results_dic = blast_genus_search(
        fasta_test_dic,
        model_predict_results_dic,
        model_info_dic,
        files_dir_path,
        save_files,
        blast_genus_max_hits
    )
    write_output(logpath, f"Time took for genus blast: {round((time.time() - blastgenus_time) / 60, 2)} minutes.\n")

    #blast against singletons
    write_output(logpath,
                 f"Running blast against singletons for tested sequences with no predicted genus identity.\n")
    singletons_time = time.time()
    singletons_test_identifiers_list = []
    for identifier in fasta_test_ident_list:
        if len(blast_results_dic[identifier]) == 0:
            singletons_test_identifiers_list.append(identifier)
    sub_fasta_test_dic = {ident: fasta_test_dic[ident] for ident in singletons_test_identifiers_list}


    singletons_blast_results = blast_singletons_search(
        sub_fasta_test_dic,
        singletons_dic,
        files_dir_path,
        batch_size,
        logpath,
        save_files,
        blast_genus_max_hits
    )
    write_output(logpath, f"Time took for singletons blast: {round(time.time() - singletons_time, 2)} seconds.\n")

    #compare and reorganize blast results
    final_blast_results_dic = {}
    for identifier in fasta_test_ident_list:

        if blast_results_dic[identifier]:

            if identifier not in singletons_blast_results.keys():
                final_blast_results_dic[identifier] = blast_results_dic[identifier]
            elif singletons_blast_results[identifier]:
                for results in blast_results_dic[identifier]:

                    if float(results[0][2]) > singletons_blast_results[identifier][0][2]:
                        final_blast_results_dic[identifier] = results

                    else:
                        ident = singletons_blast_results[identifier][0][1]
                        taxa_list1 = singletons_taxa_dic[ident].split(os.sep)[1:]
                        taxa_list = [(taxa_list1[i], None) for i in range(1, len(taxa_list1))]
                        taxa_prediction_idents_dic[identifier] = taxa_list
                        final_blast_results_dic[identifier] = singletons_blast_results[identifier]

            else:
                final_blast_results_dic[identifier] = blast_results_dic[identifier]

        elif singletons_blast_results[identifier]:
            ident = singletons_blast_results[identifier][0][1]
            taxa_list1 = singletons_taxa_dic[ident].split(os.sep)[1:]
            taxa_list = [(taxa_list1[i], None) for i in range(1, len(taxa_list1))]
            taxa_prediction_idents_dic[identifier] = taxa_list
            final_blast_results_dic[identifier] = singletons_blast_results[identifier]

        else:
            final_blast_results_dic[identifier] = []

    # format the final results
    write_output(logpath, f"Saving results...\n")
    csv_pred_file = ["Identifier,Genus prediction,Blast best hit\n"]
    csv_full_info_file = [
        "Predicted Identifier,"
        " Class Pred.,"
        " Order Pred.,"
        " Family Pred.,"
        " Genus Pred.,"
        " Blast best hit,"
        " Ident.,"
        " eval.,"
        " cov.\n"
    ]
    for ident, rest_list in final_blast_results_dic.items():
        if rest_list:
            results_list_sorted = sorted(rest_list, key=lambda x: float(x[2]), reverse=True)
            csv_pred_line = f"{ident},{taxa_prediction_idents_dic[ident][-1][0]},{results_list_sorted[0][1]}\n"
            csv_fullinfo_line = f"{ident},"

            for pred in taxa_prediction_idents_dic[ident]:
                csv_fullinfo_line += f"{pred[0]}/{round(pred[1][np.argmax(pred[1])] * 100, 3) if pred[1] else 'NA'},"

            for blast_result_list in results_list_sorted:

                for idx, blast_result_info in enumerate(blast_result_list):
                    if idx == 0: continue
                    if idx == 2: blast_result_info = str(float(blast_result_info))
                    csv_fullinfo_line += f"{blast_result_info},"

            csv_fullinfo_line += "\n"

            csv_pred_file.append(csv_pred_line)
            csv_full_info_file.append(csv_fullinfo_line)
        else:
            csv_pred_line = f">{ident},{taxa_prediction_idents_dic[ident][-1][0]}\n"
            csv_fullinfo_line = f">{ident},"

            for pred in taxa_prediction_idents_dic[ident]:
                csv_fullinfo_line += f"{pred[0]}/{round(pred[1][np.argmax(pred[1])] * 100, 3) if pred[1] else 'NA'},"
            csv_fullinfo_line += "\n"

            csv_pred_file.append(csv_pred_line)
            csv_full_info_file.append(csv_fullinfo_line)

    #save results in the files
    with open(predpath, "w") as outfile:
        for line in csv_pred_file:
            outfile.write(line)

    with open(fullinfopath, "w") as outfile:
        for line in csv_full_info_file:
            outfile.write(line)
    write_output(logpath, f"Total time taken: {round((time.time() - tottime) / 60, 2)} minutes.\n")


def model_predict_fasta_dic(
        fasta_dic,
        model_path,
        model_info_dic,
        logpath,
        batch_size,
        use_gpu
):
    # loading gpu
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        try:
            write_output(logpath, "GPU is available for usage.\n")
            tf.config.experimental.set_memory_growth(gpus[0], True)
        except RuntimeError as e:
            write_output(logpath, "error: " + str(e) + "\n")

    if use_gpu and gpus:
        write_output(logpath, "Using GPU.\n")
        tf.config.set_visible_devices(gpus[0], 'GPU')
    elif use_gpu:
        write_output(logpath, "GPU available, but set to 'False', using CPU.\n")
        tf.config.set_visible_devices([], 'CPU')
    else:
        write_output(logpath, "GPU not available, using CPU.\n")
        tf.config.set_visible_devices([], 'CPU')

    # loading models
    write_output(logpath, "Loading models...\n")
    load_time = time.time()
    loaded_models_dic = {}
    for taxa, infolist in model_info_dic.items():
        if infolist[0] != "":
            write_output(logpath, "Loading model for " + str(taxa) + "\n")
            cur_model_file = tf.keras.models.load_model(os.path.join(infolist[4], infolist[0]))
            loaded_models_dic[infolist[0]] = cur_model_file
        else:
            write_output(logpath, "No model found for " + str(taxa) + ", will be defined as a single class rank.\n")
    write_output(logpath, f"Total loading time: {round(time.time() - load_time, 2)} seconds\n")

    #get the number of taxonomic levels that exists
    total_taxa_levels = get_tot_taxa_lvls(model_path)

    #get input shape from the first model
    first_tax_class_file = model_path.split(os.sep)[-1]
    first_taxa = model_info_dic[first_tax_class_file][0]
    input_shape = loaded_models_dic[first_taxa].layers[0].input_shape[1]

    # encode the test sequences using the input shape
    write_output(logpath, "Encoding...\n")
    encd_time = time.time()
    enconded_fasta_dic = {}
    for ident, seq in fasta_dic.items():
        seq_matrix = encode_sequence(seq, input_shape)
        enconded_fasta_dic[ident] = seq_matrix
    write_output(logpath, f"Time took encoding: {round(time.time() - encd_time, 2)} seconds\n")

    #first predictions round using all sequences against the first model taxa, i.e. Cressdnavicota in the ex
    #dictionary for each taxonomic class and related identifiers as a list with prediction values when available
    test_identifiers_dic = {first_tax_class_file: [(identifier, None) for identifier in fasta_dic.keys()]}
    prediction_idents_dic = {}
    ranks_done_list = []

    # start predictions against each level of the taxonomic ranks
    write_output(logpath, f"Starting predictions against each level of the taxonomic ranks.\n")
    taxa_lvl_dic = {0: "Genus", 1: "Family", 2: "Order", 3: "Class", 4: "Phylum", 5: "Kingdom", 6: "Realm"}
    while total_taxa_levels >= 0:
        write_output(logpath, f"### Predicting for the {taxa_lvl_dic[total_taxa_levels + 1]} level. ###\n")
        aux_test_ident_dic = {}

        for taxa_rank, taxa_ident_list in test_identifiers_dic.items():
            if taxa_rank in ranks_done_list: continue
            write_output(logpath, f"Current taxa: {taxa_rank}; total sequences on it:{len(taxa_ident_list)};\n")

            #get taxa information
            rank_time = time.time()
            model_file = model_info_dic[taxa_rank][0]
            categories_file = model_info_dic[taxa_rank][1]
            fasta_taxa_subranks_list = model_info_dic[taxa_rank][2]
            dir_path = model_info_dic[taxa_rank][4]

            if model_file != "":
                #get all encoded sequences for the taxa rank identifiers list
                encoded_sequences = np.array([enconded_fasta_dic[k[0]] for k in taxa_ident_list])

                #start the prediction batches
                tf.keras.backend.clear_session()
                model = loaded_models_dic[model_file]
                gen = BatchGenerator(encoded_sequences, batch_size=batch_size)
                predict1 = model.predict(gen, verbose=0)

                #parse predictions saving in the 2 results dics
                prediction_classes_list = get_prediction_classes(os.path.join(dir_path, categories_file))
                pred_dic = {}
                for idx, pred in enumerate(predict1):
                    for category_name, category_encoded_list in prediction_classes_list:
                        if np.argmax(np.array(category_encoded_list)) == np.argmax(pred):
                            pred_dic[taxa_ident_list[idx][0]] = (category_name, list(pred))

                for identifier, prediction in pred_dic.items():

                    if identifier not in prediction_idents_dic.keys():
                        prediction_idents_dic[identifier] = [prediction]
                    else:
                        prediction_idents_dic[identifier].append(prediction)

                    if prediction[0] not in aux_test_ident_dic.keys():
                        aux_test_ident_dic[prediction[0]] = [(identifier, prediction[1])]
                    else:
                        aux_test_ident_dic[prediction[0]].append((identifier, prediction[1]))

            else:
                aux_test_ident_dic[(fasta_taxa_subranks_list[0].split("."))[0]] = []
                for ident in taxa_ident_list:
                    aux_test_ident_dic[(fasta_taxa_subranks_list[0].split("."))[0]].append((ident[0], None))
                    if ident[0] not in prediction_idents_dic.keys():
                        prediction_idents_dic[ident[0]] = [((fasta_taxa_subranks_list[0].split("."))[0], None)]
                    else:
                        prediction_idents_dic[ident[0]].append(((fasta_taxa_subranks_list[0].split("."))[0], None))

            #save as done
            ranks_done_list.append(taxa_rank)
            write_output(logpath,
                         f"Time took on the taxa:{round(time.time() - rank_time, 2)} seconds\n")

        #save in the pred dic
        test_identifiers_dic.update(aux_test_ident_dic)
        total_taxa_levels -= 1

    return test_identifiers_dic, prediction_idents_dic


def blast_genus_search(
        fasta_dic,
        model_predict_results_dic,
        model_info_dic,
        files_path,
        save_files,
        blast_genus_max_hits
):
    blast_results_dic = {}

    for tax_class, data_info in model_info_dic.items():
        #when there is no other directory the class is a family
        if not data_info[3]:
            for genus in model_info_dic[tax_class][2]:
                if genus.split(".")[0] in model_predict_results_dic.keys():
                    identifiers_results_dic = {}
                    #get tested sequences for the predicted genus

                    idents_list = [ident for ident, pred_val in model_predict_results_dic[genus.split(".")[0]]]
                    fasta_list = [(identifier, fasta_dic[identifier]) for identifier in idents_list]
                    fasta_list_path = os.path.join(files_path, genus.split(".")[0] + "_blast_test.fasta")
                    with open(fasta_list_path, "w+") as file:
                        for ident1, seq1 in fasta_list:
                            file.write(">" + ident1 + "\n" + seq1 + "\n")

                    #get the genus fasta sequences file
                    genus_fasta_path = os.path.join(model_info_dic[tax_class][4], genus)

                    #create blast output file
                    seq_result_file_path = os.path.join(files_path, genus.split(".")[0] + "_blast_results.txt")

                    #run blast
                    blast_test_results = blast_sequences_list(fasta_list_path, genus_fasta_path, seq_result_file_path)

                    #if no result is found
                    if not blast_test_results[0]:
                        for ident2 in idents_list:
                            identifiers_results_dic[ident2] = []
                        if os.path.exists(seq_result_file_path):
                            os.remove(seq_result_file_path)

                        blast_results_dic.update(identifiers_results_dic)

                        continue

                    if not save_files and os.path.exists(seq_result_file_path):
                        os.remove(fasta_list_path)
                        os.remove(seq_result_file_path)

                    for result in blast_test_results[0]:
                        if result[0] in identifiers_results_dic.keys():
                            if len(identifiers_results_dic[result[0]]) < blast_genus_max_hits:
                                identifiers_results_dic[result[0]].append(result)
                        else:
                            identifiers_results_dic[result[0]] = [result]

                    for identifier in idents_list:
                        if identifier not in identifiers_results_dic.keys():
                            identifiers_results_dic[identifier] = []

                    #save results for each genus
                    blast_results_dic.update(identifiers_results_dic)
    return blast_results_dic


def blast_singletons_search(
        fasta_dic,
        singletons_dic,
        files_path,
        batch_size,
        log_path,
        save_files,
        blast_genus_max_hits
):
    #create singletons file, if not exists
    singletons_fasta_path = os.path.join(files_path, "singletons_len" + str(len(singletons_dic)) + ".fasta")
    if os.path.exists(singletons_fasta_path):
        os.remove(singletons_fasta_path)

    with open(singletons_fasta_path, "w") as f2:
        for identifier, sequence2 in singletons_dic.items():
            f2.write(
                ">" + identifier.replace("\n", "") + "\n" + create_fasta_seq(sequence2.replace("\n", "")) + "\n")

    idents_list = list(fasta_dic.keys())

    identifiers_results_dic = {}
    seq_result_file_path = os.path.join(files_path, "singletons_blast_results.txt")
    for idx in range(0, len(idents_list), batch_size):
        cur_idents_batch = idents_list[idx:idx + batch_size]
        fasta_list = [(identifier, fasta_dic[identifier]) for identifier in cur_idents_batch]
        fasta_list_path = os.path.join(files_path, "batchs_blast_test.fasta")
        with open(fasta_list_path, "w+") as file:
            for ident1, seq1 in fasta_list:
                file.write(">" + ident1 + "\n" + seq1 + "\n")

        blast_test_results = blast_sequences_list(fasta_list_path, singletons_fasta_path, seq_result_file_path, )
        if not blast_test_results[0]:
            write_output(log_path,
                         "Error while testing for: Singletons blast, err:" + blast_test_results[1] + "\n")
            write_output(log_path, "Error identifiers:\n")
            for ident2 in cur_idents_batch:
                write_output(log_path, ">" + ident2 + "\n")
            if os.path.exists(fasta_list_path):
                os.remove(fasta_list_path)
            continue

        if not save_files:
            os.remove(fasta_list_path)
            os.remove(singletons_fasta_path)
            if os.path.exists(seq_result_file_path):
                os.remove(seq_result_file_path)

        # parse results
        aux_identifiers_results_dic = {}
        for result in blast_test_results[0]:
            if result[0] in aux_identifiers_results_dic.keys():
                if len(aux_identifiers_results_dic[result[0]]) < blast_genus_max_hits:
                    aux_identifiers_results_dic[result[0]].append(result)
            else:
                aux_identifiers_results_dic[result[0]] = [result]

        for identifier in idents_list:
            if identifier not in aux_identifiers_results_dic.keys():
                aux_identifiers_results_dic[identifier] = []

        identifiers_results_dic.update(aux_identifiers_results_dic)

    return identifiers_results_dic
