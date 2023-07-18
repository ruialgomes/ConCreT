import tkinter as tk
from tkinter import filedialog

class GUIApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Concret")

        # Variables to store the user inputs
        self.inputfile = tk.StringVar()
        self.modelpath = tk.StringVar()
        self.outputpath = tk.StringVar()
        self.outputname = tk.StringVar()
        self.save = tk.BooleanVar()
        self.max_blast_hits = tk.IntVar()
        self.gpu = tk.BooleanVar()

        # Set default values for max_blast_hits and save options
        self.max_blast_hits.set(3)
        self.save.set(False)

        # Labels on the left side
        label_inputfile = tk.Label(self.root, text="Path to input file:")
        label_inputfile.grid(row=0, column=0, sticky="w")
        label_modelpath = tk.Label(self.root, text="Path to model directory:")
        label_modelpath.grid(row=1, column=0, sticky="w")
        label_outputpath = tk.Label(self.root, text="Path to output directory:")
        label_outputpath.grid(row=2, column=0, sticky="w")
        label_outputname = tk.Label(self.root, text="Results-name-default:")
        label_outputname.grid(row=3, column=0, sticky="w")
        label_max_blast_hits = tk.Label(self.root, text="Max blast hits:")
        label_max_blast_hits.grid(row=4, column=0, sticky="w")
        label_save = tk.Label(self.root, text="Save result files:")
        label_save.grid(row=5, column=0, sticky="w")
        label_gpu = tk.Label(self.root, text="Use GPU (if available):")
        label_gpu.grid(row=6, column=0, sticky="w")

        # Entry and Button widgets
        tk.Entry(self.root, textvariable=self.inputfile).grid(row=0, column=1, padx=5, pady=5)
        tk.Button(self.root, text="Browse", command=self.browse_inputfile).grid(row=0, column=2)

        tk.Entry(self.root, textvariable=self.modelpath).grid(row=1, column=1, padx=5, pady=5)
        tk.Button(self.root, text="Browse", command=self.browse_modelpath).grid(row=1, column=2)

        tk.Entry(self.root, textvariable=self.outputpath).grid(row=2, column=1, padx=5, pady=5)
        tk.Button(self.root, text="Browse", command=self.browse_outputpath).grid(row=2, column=2)

        tk.Entry(self.root, textvariable=self.outputname).grid(row=3, column=1, padx=5, pady=5)

        tk.Spinbox(self.root, from_=1, to=15, textvariable=self.max_blast_hits).grid(row=4, column=1, padx=5, pady=5)

        tk.Checkbutton(self.root, variable=self.save).grid(row=5, column=1, padx=5, pady=5)

        tk.Checkbutton(self.root, variable=self.gpu).grid(row=6, column=1, padx=5, pady=5)

        # Submit button
        submit_button = tk.Button(self.root, text="Submit", command=self.on_submit)
        submit_button.grid(row=7, column=1, padx=5, pady=10)

    def browse_inputfile(self):
        file_path = filedialog.askopenfilename()
        self.inputfile.set(file_path)

    def browse_modelpath(self):
        dir_path = filedialog.askdirectory()
        self.modelpath.set(dir_path)

    def browse_outputpath(self):
        dir_path = filedialog.askdirectory()
        self.outputpath.set(dir_path)

    def on_submit(self):
        # Close the GUI window
        self.root.destroy()

    def run(self):
        self.root.geometry("400x300")  # Set the window size here
        self.root.mainloop()

if __name__ == "__main__":
    app = GUIApp()
    app.run()
