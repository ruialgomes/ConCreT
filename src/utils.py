import numpy as np

def write_output(path, output):
    with open(path, "a+") as outfile:
        outfile.write(output)

def open_fasta(path):
    fasta_dic = {}
    fasta_seq = ""
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if len(fasta_dic) == 0:
                    seq_id = line[1:].replace("\n", "").replace(" ","_")
                    fasta_dic[seq_id] = ""
                    continue
                else:
                    fasta_dic[seq_id] = fasta_seq
                    fasta_seq = ""
                seq_id = line[1:].replace("\n", "").replace(" ","_")
                fasta_dic[seq_id] = fasta_seq
            else:
                fasta_seq += line.replace("\n", "")
        else:
            fasta_dic[seq_id] = fasta_seq
    return fasta_dic

def create_fasta_seq(seq_in):
    seq_out = ""
    for index_let, letter in enumerate(seq_in):
        seq_out += letter
        if index_let % 73 == 0 and index_let != 0:
            seq_out += "\n"
    return seq_out

def repeat_to_length(string_to_expand, length):
    return (string_to_expand * (int(length / len(string_to_expand)) + 1))[:length]

def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def has_lower(inputString):
    return any(char.islower() for char in inputString)

def has_letter(inputString):
    return any(char.isalpha() for char in inputString)

def encode_sequence(seq, input_length):
    nseq_list = []
    if len(seq) <= input_length:
        nseq = repeat_to_length(seq, input_length)
        nseq_list.append(nseq)
    else:
        nseq = seq[:input_length]
        nseq_list.append(nseq)

    mapping = dict(zip("ATCG-KNSMYRWD", range(13)))
    seq2 = [mapping[i] for i in nseq]
    seq_list = []
    for val in seq2:
        if val in (4, 5, 6, 7, 8, 9, 10, 11, 12):
            seq_list.append(np.asmatrix(np.zeros(4)))
        else:
            seq_list.append(np.asmatrix(np.eye(4)[val]))
    sequence_matrix = np.concatenate(tuple(seq_list))

    return sequence_matrix