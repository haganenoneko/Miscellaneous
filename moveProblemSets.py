"""Move all PDF files with suffix `prob.pdf`, `sol.pdf`, `sum.pdf`, or `ex*.pdf` into one folder above, then delete all remaining files"""
import os, glob, shutil
from tkinter import filedialog, Tk 
from typing import List 

def yes_or_no(question: str) -> bool:
	"""https://gist.github.com/garrettdreyfus/8153571"""
	reply = str(input(question+' (y/n): ')).lower().strip()
	if reply[0] == 'y':
		return True
	if reply[0] == 'n':
		return False
	else:
		return yes_or_no("Uhhhh... please enter ")

class mover:
	def __init__(self, folder: str=None, suffixes: List[str]=None):
		
		self._getSuffixes(suffixes)		
		self._getFolderPath(folder)		
		paths = self._getFilePaths()		
		self._moveFiles(paths)		
		self._deleteRemainder()
	
	def _getSuffixes(self, suffixes: List[str]):
		if suffixes is None:
			self._suffixes = ['prob.pdf', 'sol.pdf', 'sum.pdf', 'ex*.pdf']
		else:
			for i, s in enumerate(suffixes):
				if '.pdf' == s[:-4]: continue 
				else:
					suffixes[i] += '.pdf'
			
			self._suffixes = suffixes 
		
	def _getFolderPath(self, folder: str):
		"""start a tkinter instance and retrieve the filename with a tkinter gui"""
		
		if folder is not None:
			self._folderPath = folder 
			print(f"Selected folder: {folder}")
			return 		

		root = Tk() 
		root.filename =  filedialog.askdirectory(
			initialdir = "/", title = "Select directory", 
		)
			
		if root.filename is None:
			raise ValueError("No directory selected.")
			
		print(f"Selected folder: {root.filename}")
		self._folderPath = root.filename 
		
	def _getFilePaths(self) -> List[str]:
		"""Search through all subdirectories of `self._folderPath` for files matching `self._suffixes`"""
		pattern = self._folderPath + "/**/*{0}"
		
		all_paths = [] 
		
		for suffix in self._suffixes:
			files = glob.glob( pattern.format(suffix), recursive=True ) 
			if not files: continue 
			
			all_paths.extend(files)
		
		if not all_paths: 
			print(f"No files were found with the pattern:\n {pattern}")
			return [] 
				
		return all_paths 
	
	def _moveFiles(self, files: List[str]):
		"""Move files from their current path to `self._folderPath`"""
		
		if not files:
			print("No files to move")
			return
		
		new_file_name = self._folderPath + "/{0}"
		
		for f in files:
			fname = os.path.basename(f)
			try: 
				os.rename(f, new_file_name.format(fname))
			except WindowsError as e:
				if e.winerror == 183:
					print("File already exists in destination. Removing original.")
					os.remove(f)
				else:
					print(f"Failed to move file\n {f}\n due to exception:\n{e}.")
		
		print(f"{len(files)} were successfully moved to {new_file_name[:-3]}")
		return 
				
	def _deleteRemainder(self):
		"""
		Delete all files in `self._folderPath` that do not have a .pdf extension
		https://stackoverflow.com/questions/185936/how-to-delete-the-contents-of-a-folder
		"""
		
		delete_files = yes_or_no("Do you want to delete remaining files?")
		if not delete_files:
			return 			
		
		with os.scandir(self._folderPath) as files:
				
			for f in files:
				if ".pdf" == f.name[-4:]: continue 
				print(f"Removing... {f.name}")

				if os.path.isdir(f.path):
					shutil.rmtree(f.path)
				else:
					os.remove(f.path)
	
# defaults are for MIT linalg 
# p6.041: suffixes = ['sol', 'final', 'f09', 's09', '_assn*', '_L*', '_rec*', '_tut*', '_qu*']
	
m = mover(
	folder = r"C:/Users/delbe/Documents/Code_Repositories/MIT 6.041/6-041sc-fall-2013/final-exam",
	suffixes = ['sol', 'final', 'f09', 's09', '_assn*', '_L*', '_rec*', '_tut*', '_qu*']
)

			
	
	