# ConCreT


*** THIS TOOL IS UNDER DEVELOPMENT***

Welcome, this is a notebook script for a predictor based on convolutional neural networks called Convolutional neural networks for Cressdnaviricota Taxonomy. Using this script you can analyze fasta/multifasta files of Cressdnaviricota viruses, and quickly get a taxonomy prediction with blast confirmation. Currently, the dataset used for predictions is based on the ICTV VMR release (https://ictv.global/vmr) from 08/31/2022. 

The faster way to use this tool is to run the from the Colab directory, which already has Python3 and TensorFlow for the main script ‘Concret_predictor.ipynb’. If you wish to run it locally is necessary to install Python and TensorFlow in one environment, and then use Python3 to run the ‘concret_prediction.py’, as better explained below. Also, as BLAST is used for each sample, for running locally you must install it as explained for different OS in the link: (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html). For Colab this is done already in the script.

This tool is in version 1.0, and still under development, in case o errors, please contact us.

Also available on Colab: [https://github.com/ruialgomes/ConCre](https://drive.google.com/drive/u/0/folders/1xL2dN3VTwxi9XZCs7oPM4rbdVOjaHTBl)T 

Any questions or suggestions, feel free to contact me at the email address: ruitherarthur@gmail.com 

**Usage**

**Google Colab**
Download the entire directory, extract it to a folder, and upload the entire folder to your drive. From there, open the .ipynb file with Google Colaboratory can change the paths to match your drive’s path and run the script.

**Script .py**
Download the folder “ConCreTv1.0”, and extract it from the zip file at the desired directory. To run the script open the terminal, and navigate to the extraction directory. With Python3 already installed, inside the folder, you will create a new environment for TensorFlow installation. You can do that by using the command:

python -m venv ’environment_name’  	for Windows.
or
python3 -m venv ’environment_name’ 	for Unix/Mac.

Where you will change  “environment_name” to a name of your preference, like “concret_venv”. Then you need activate the environment using:

’environment_name’/Scripts/activate		for Windows
or
source ’environment_name’/bin/activate 		for Unix/Mac

Now you can install the libraries required to run (TensorFlow  and Numpy) using:

pip install -r requirements.txt

Then you can use the command line options, with more options below:

python concret_prediction.py --in “yourFasta.fasta” 		for Windows
or
python3 concret_prediction.py --in “yourFasta.fasta” 		for  Unix/Mac

Or use the GUI option with:

python concret_prediction_GUI.py 		for Windows
or
python3 concret_prediction_GUI.py 		for  Unix/Mac

When you finish, deactivate the venv with:

deactivate

**Execultables**
**
Linux64X**  ***Under tests***
Download the folder “ConCreTv1.0_linux64”, and extract it from the zip file at the desired directory. To run the script open the terminal, navigate to the extraction directory, and run the line command, explained below, directly from it. The output will be saved in the output folder by default.


**Command line options**

--in			The input path to the query file. The default is the example file.
--mode	l 		A path to the model directory. The default is the Cressdnaviricota model.
--out			A path for the output directory, where an output folder will be created. The default is the output directory.
'--output-name	 	A name for the result file. The default is the input file name.
--save			To select if you want to save the blast results files. The default is no.
--max-blast-hits	Maximum of BLAST hits that 'will be printed in the output. The default is 3.
--gpu			Use GPU (if available). The default is using CPU.
**Command line examples**

“concret_prediction  --in ’path/to/your/file.fasta’”




