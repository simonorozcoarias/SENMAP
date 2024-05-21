import os
import sys
import time
import argparse

import numpy as np
from numpy import argmax
from Bio import SeqIO
import tensorflow as tf
from tensorflow.keras import backend as K
from focal_loss import BinaryFocalLoss


def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall


def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision


def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2 * ((precision * recall) / (precision + recall + K.epsilon()))


def write_sequences_file(sequences, filename):
    try:
        SeqIO.write(sequences, filename, "fasta")
    except FileNotFoundError:
        print("FATAL ERROR: I couldn't find the file, please check: '" + filename + "'. Path not found")
        sys.exit(0)
    except PermissionError:
        print("FATAL ERROR: I couldn't access the files, please check: '" + filename + "'. I don't have permissions.")
        sys.exit(0)
    except Exception as exp:
        print("FATAL ERROR: There is a unknown problem writing sequences in : '" + filename + "'.")
        print(exp)
        sys.exit(0)


def create_output_folders(folder_path):
    try:
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
    except FileNotFoundError:
        print("FATAL ERROR: I couldn't create the folder " + folder_path + ". Path not found")
        sys.exit(0)
    except PermissionError:
        print("FATAL ERROR: I couldn't create the folder " + folder_path + ". I don't have permissions.")
        sys.exit(0)


def delete_files(file_path):
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
        except PermissionError:
            print("WARNING: The file " + file_path + " couldn't be removed because I don't have permissions.")


"""
This converts a fasta sequences (in nucleotides) to one-hot representation
"""
def fasta2onehot(library, max_len):
    sequences = [x for x in SeqIO.parse(library, "fasta")]
    IDs = [x.id for x in sequences]
    langu = ['A', 'C', 'G', 'T', 'N']
    rep2d = np.zeros((len(sequences), len(langu), max_len), dtype=bool)

    idx_seq = 0
    for sequence in sequences:
        posNucl = 0
        if len(sequence.seq) < max_len:
            # to complete TEs with NNs centering the sequence
            times = int((max_len-len(sequence.seq))/2)
            seq = str('N'*times+str(sequence.seq)+'N'*(times+1))[0:max_len]
        else:
            seq = sequence.seq[0:max_len]

        for nucl in seq:
            posLang = langu.index(nucl.upper())

            rep2d[idx_seq][posLang][posNucl] = 1
            posNucl += 1

        idx_seq += 1
    
    return rep2d, IDs


"""
This converts a one-hot representation of a sequence to its fasta form (in nucleotides)
"""
def onehot2fasta(dataset):
    langu = ['A', 'C', 'G', 'T', 'N']
    fasta_seqs = ""
    for j in range(dataset.shape[1]):
        if sum(dataset[:, j]) > 0:
            pos = argmax(dataset[:, j])
            fasta_seqs += langu[pos]
    return fasta_seqs



def model_predict(model_path, dataX):
    model = tf.keras.models.load_model(model_path,
                                       custom_objects={'LeakyReLU': tf.keras.layers.LeakyReLU(0.01),
                                                       'f1_m': f1_m})
    predictions = model.predict(dataX, verbose=0)
    return predictions


def curate_libs(predictions, IDs, kept_thr, outputdir, library):
    kept_seqs = []
    remove_seqs = []
    predictions_out = open(outputdir + "/predictions_log.txt", "w")
    label_predicted = [argmax(x) for x in predictions]
    for i in range(len(label_predicted)):
        removed = ""
        if label_predicted[i] < kept_thr:
            remove_seqs.append(IDs[i])
            removed = "Removed"
        else:
            kept_seqs.append(IDs[i])
            removed = "Kept"
        predictions_out.write(IDs[i] + "\t" + str(predictions[i, 0]) + "\t" + str(predictions[i, 1]) + "\t" + removed + "\n")

    predictions_out.close()
    write_sequences_file([x for x in SeqIO.parse(library, "fasta") if x.id in kept_seqs], outputdir+"/kept_sequences.fa")
    write_sequences_file([x for x in SeqIO.parse(library, "fasta") if x.id in remove_seqs], outputdir + "/removed_sequences.fa")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', required=True, dest='library',
                        help='Path to the library in fasta format. Required*')
    parser.add_argument('-v', required=False, dest='verbose', default=True,
                        help='Verbose: activate (True) or desactivate (False)')
    parser.add_argument('-o', '--output', required=True, dest='outputdir',
                        help='Path to the output directory. Required*')
    options = parser.parse_args()
    library = options.library
    outputdir = options.outputdir
    verbose = options.verbose

    if verbose is None:
        verbose = True

    installation_path = os.path.dirname(os.path.realpath(__file__))
    model_path = installation_path + "/trained_model.hdf5"
    max_len = 23200
    kept_thr = 0.7

    # Step 1: Convert fasta libs into one_hot 
    start_time = time.time()
    dataX, IDs = fasta2onehot(library, max_len)
    end_time = time.time()
    if verbose:
        print("Sequences converted into one-hot coding scheme successfully [" + str(
            end_time - start_time) + " seconds]")

    # Step 2: SENMAP predictions
    start_time = time.time()
    predictions = model_predict(model_path, dataX)
    end_time = time.time()
    if verbose:
        print("Predictions done by the Neural Network [" + str(
            end_time - start_time) + " seconds]")

    # Step 3: Generate final files
    start_time = time.time()
    curate_libs(predictions, IDs, kept_thr, outputdir, library)
    end_time = time.time()
    if verbose:
        print("Final files created successfully [" + str(
            end_time - start_time) + " seconds]")
    
