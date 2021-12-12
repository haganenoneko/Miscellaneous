"""Move all PDF files with suffix `prob.pdf`, `sol.pdf`, `sum.pdf`, or `ex*.pdf` into one folder above, then delete all remaining files"""
import os, glob 
from tkinter import filedialog, Tk 
from typing import List 

def yes_or_no(question: str) -> bool:
	"""https://gist.github.com/garrettdreyfus/8153571"""
    reply = str(raw_input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return True
    if reply[0] == 'n':
        return False
    else:
        return yes_or_no("Uhhhh... please enter ")

class mover:
	def __init__(self, folder: str=None):
		self._suffixes = ['prob.pdf', 'sol.pdf', 'sum.pdf', 'ex*.pdf']
		
		self._getFolderPath(folder)
		paths = self._getFilePaths()
		self._moveFiles(paths)
		self._deleteRemainder()
		
	def _getFolderPath(self, folder: str):
		"""start a tkinter instance and retrieve the filename with a tkinter gui"""
		
		if folder is not None:
			self._folderPath = folder 
			return 		
		
		try:
			root = Tk() 
			root.filename =  filedialog.askdirectory(
				initialdir = "/", title = "Select directory", 
			)
		except KeyboardInterrupt:
			print("Dialog closed without setting a folder path. Exiting...")
			exit()
		
		print(f"Selected folder: {root.filename}")
		self._folderPath = root.filename 
		
	def _moveFiles(self, files: List[str]):
		"""Move files from their current path to `self._folderPath`"""
		
		new_file_name = self._folderPath + "/{0}"
		
		for f in files:
			fname = os.path.basename(f)
			try: 
				os.rename(f, new_file_name.format(fname))
			except Exception as e:
				print(f"Failed to move file\n {f}\n due to exception:\n{e}.")
		
		print(f"{len(files)} were successfully moved to {new_file_name[:-3]}")
		
	def _deleteRemainder(self):
		"""
		Delete all files in `self._folderPath` that do not have a .pdf extension
		https://stackoverflow.com/questions/185936/how-to-delete-the-contents-of-a-folder
		"""
		
		delete_files = yes_or_no("Do you want to delete remaining files?")
		if not delete_files:
			return 			
		
		for filename in os.scandir(self._folderPath):
			if ".pdf" == filename[-4:]: continue 
			print(f"Removing... {filename}")
			os.remove(filename)
	
	def _getFilePaths(self) -> List[str]:
		"""Search through all subdirectories of `self._folderPath` for files matching `self._suffixes`"""
		pattern = self._folderPath + "/**/*{0}"
		
		all_paths = [] 
		
		for suffix in self._suffixes:
			files = glob.glob( pattern.format(suffix), recursive=True ) 
			print(suffix)
			if not files: continue 
			
			all_paths.extend(files)
		
		if not all_paths: 
			raise ValueError(f"No files were found with the pattern:\n {pattern}")
		
		return all_paths 
			
m = mover()

			
	
	