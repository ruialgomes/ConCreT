import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import argparse
import tensorflow as tf
from src.predictor import parse_fasta_file

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        tf.config.experimental.set_memory_growth(gpus[0], True)
    except RuntimeError as e:
        print("Error:", e)

parser = argparse.ArgumentParser()
parser.add_argument(
    '--in',
    dest='inputfile',
    default="examples/example_test_allgenus.fasta",
    help='The input path to query file. Default is the example file.'
)
parser.add_argument(
    '--model',
    dest='modelpath',
    default="models/Cressdnaviricota",
    help='A path to the model directory. Default is the Cressdnaviricota model.'
)
parser.add_argument(
    '--out',
    dest='outputpath',
    default="output",
    help='A path for the output directory, where a output folder will be created. Default is the output directory.'
)
parser.add_argument(
    '--outputname',
    default="Results-name-default",
    help='A name for the result file. Default is the input file name.'
)
parser.add_argument(
    '--save',
    default=True,
    help='To select if you want to save the blast results files. Default is no.'
)
parser.add_argument(
    '--max-blast-hits',
    dest='max_blast_hits',
    default=2,
    help='Maximum of BLAST hits that will be printed in the output. Default is 3.'
)
parser.add_argument(
    '--gpu',
    help='Use GPU for the analysis when available. Default is using CPU.',
    default=True
)
parser.add_argument(
    '--batch_size',
    help='The number of sequences predicted at a time, reducing memory usage. Default is 128',
    default=128
)


def main():
    args = parser.parse_args()

    init_test_file_path = args.inputfile
    model_path = args.modelpath
    output_path = args.outputpath
    output_name = args.outputname
    save_files = args.save
    use_gpu = args.gpu
    blast_genus_max_hits = args.max_blast_hits
    batch_size = args.batch_size

    print("Starting...")
    #parsing the test input
    parse_fasta_file(
        init_test_file_path,
        model_path,
        output_path,
        output_name,
        save_files,
        use_gpu,
        blast_genus_max_hits,
        batch_size
    )
    print("Ended.")


if __name__ == "__main__":
    main()
