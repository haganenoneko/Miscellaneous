"""
Extract audio from a video file using ffmpeg

Instructions:
1. Open a command prompt window in the folder where this script is saved 
2. Run `python _ex.py [input] [outut]`, where `input` and `output` are 
	optional values for input and output filenames, respectively. 

Copyright 2021 Delbert Yip (Github: haganenoneko)
"""

import os 
import argparse 
import tkinter as tk
from tkinter import filedialog

def open_fd() -> str:
	"""
	Open tkinter file dialog to select input 
	video file to extract audio from
	"""
	root = tk.Tk()
	root.withdraw()
	file_path = filedialog.askopenfilename(initialdir=os.getcwd())
	
	if os.path.isfile(file_path):
		return file_path 
	else:
		raise FileExistsError(f"{file_path} is not a valid file")

parser = argparse.ArgumentParser()

# optional arguments
parser.add_argument("--filename", help="filename", type=str)
parser.add_argument("--dest", help="name of output file", type=str)
args = parser.parse_args()

fname = open_fd() if not args.filename else args.filename 
dest = fname if not args.dest else args.dest 

print(f"Extracting from {fname}. Output will be saved to {dest}.m4a")

cmd = f"ffmpeg -i {fname} -vn -acodec copy {dest}.m4a"
print(f"Running...\n{cmd}")

os.system(cmd)

"""
References

Tkinter file dialog: 
https://stackoverflow.com/questions/9319317/quick-and-easy-file-dialog-in-python

Initial tkinter file dialog path: 
https://stackoverflow.com/questions/31122704/specify-file-path-in-tkinter-file-dialog

Extract audio with ffmpeg: 
https://stackoverflow.com/questions/9913032/how-can-i-extract-audio-from-video-with-ffmpeg
"""

