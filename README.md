# ConCreT

Welcome, this is a script for a predictor based on convolutional neural networks called Convolutional Neural Networks for Cressdnaviricota Taxonomy. Using this script you can analyze fasta/multi-fasta files of Cressdnaviricota viruses, and quickly get a taxonomy prediction with blast confirmation. Currently, the dataset used for predictions is based on the ICTV VMR release (https://ictv.global/vmr) from 08/31/2022. 

The faster way to use this tool is to run the from the Colab directory, which already has Python3 and TensorFlow for the main script ‘Concret_predictor.ipynb’. If you wish to run it locally is necessary to install Python and TensorFlow in one environment, and then use Python3 to run the ‘concret_prediction.py’, as better explained below. Also, as BLAST is used for each sample, for running locally you must install it as explained for different OS in the link: (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). For Colab this is already done in the script.

For those who prefer graphical interfaces to select the files and parameters, there is also a GUI version. For GPU predictions, with a really faster execution time, it is necessary to be available on the computer, which is possible for some graphic cards, directly for NVIDIA, but you might need to install some extra libraries in TensorFlow for Intel or AMD.  This tool is in version 1.0, and still under development, in case of errors, please contact us.

If you use our tool, please remember to cite: https://pubmed.ncbi.nlm.nih.gov/37536450/ 

Any questions or suggestions, feel free to contact me at the email address: ruitherarthur@gmail.com 

**Usage**<br/>

Download the folder “ConCreT”, and extract it from the zip file at the desired directory. To run the script open the terminal, and navigate to the extraction directory. With Python3 already installed, inside the folder, you will create a new environment for TensorFlow installation. You can do that by using the command:

_python -m venv ’environment_name’_  	for Windows.
or
_python3 -m venv ’environment_name’_	for Unix/Mac.

Where you will change  “environment_name” to a name of your preference, like “concret_venv”. Then you need to activate the environment using:

_’environment_name’/Scripts/activate_	for Windows
or
_source ’environment_name’/bin/activate_	for Unix/Mac

Now you can install the libraries required to run (TensorFlow  and Numpy) using:

_pip install -r requirements.txt_

Then you can use the command line options, with more options below:

_python concret.py --in “yourFasta.fasta”_	for Windows
or
_python3 concret.py --in “yourFasta.fasta”_ 	for  Unix/Mac

When you finish, deactivate the venv with:

_deactivate_

**Command line**<br/>

--in &nbsp;&nbsp;&nbsp;&nbsp;The input path to the query file. The default is the example file.

--model &nbsp;&nbsp;&nbsp;&nbsp; A path to the model directory. The default is the Cressdnaviricota model.

--out &nbsp;&nbsp;&nbsp;&nbsp; A path for the output directory, where an output folder will be created. The default is the output directory.

'--output-name &nbsp;&nbsp;&nbsp;&nbsp;A name for the result file. The default is the input file name.

--save &nbsp;&nbsp;&nbsp;&nbsp; To select if you want to save the blast results files. The default is no.

--max-blast-hits &nbsp;&nbsp;&nbsp;&nbsp; Maximum of BLAST hits that 'will be printed in the output. The default is 3.

--gpu &nbsp;&nbsp;&nbsp;&nbsp;Use GPU (if available). The default is using CPU.


**Command line example**

python3 Concret.py  --in ’path/to/your/file.fasta’




