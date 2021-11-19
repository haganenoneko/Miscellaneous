# pymol publication quality
# https://bionerdnotes.wordpress.com/2018/11/12/getting-high-quality-pictures-in-pymol/

from pymol import cmd
import os 

"""
### Ideally, we could set the below, but there are no 'cmd' for these functions, so do them manually

set cartoon_highlight_color, grey75
set antialias, 4

set ray_trace_mode, 1
set ray_trace_gain, 0.005
set cartoon_discrete_colors, 1
set ray_shadows, 0   
"""

_u = [
    "set cartoon_highlight_color, grey75",
    "set antialias, 4"n
    "set ray_trace_mode, 1",
    "set ray_trace_gain, 0.005",
    "set cartoon_discrete_colors, 1",
    "set ray_shadows, 0"
]

print("Please set the following: \n", "\n".join(_u))

def _reset():
    cmd.reinitialize("everything")
    cmd.reinitialize("original_settings")

class load_structures:
    def __init__(self, path="", names=[]):
        
        if not os.path.isdir(path):
            path = r"C:/Users/delbe/Downloads/wut/wut/Post_grad/UBC/Research/lab/structure files/pdb_files/HCN2_MOUSE_O88703_Potassiumsodium_hyperpolarizati/models/"
            
        if len(names) < 1:
            names = ["7nmn/7nmn", "6uqf/6uqf", "7np3_mHCN2"]
        
        self._paths = [path + n + ".pdb" for n in names]
        
        self._views = [
            # dummy 
            (0.),
            
            # 5 close ups for multi-model 
            (0.927794456,    0.365347803,   -0.075607292,
            0.357938945,   -0.928815305,   -0.095773198,
            -0.105213568,    0.061797936,   -0.992521822,
            -0.001014138,    0.001799855,  -40.068046570,
            144.100341797,  134.904006958,  115.075866699,
            14.193284035,   66.120056152,  -20.000000000 ),
            (0.918314517,    0.308308542,    0.248276174,
            0.325078785,   -0.945250809,   -0.028589115,
            0.225866720,    0.106965333,   -0.968261242,
            -0.000293555,    0.001118235,  -50.398750305,
            146.105010986,  136.049407959,  114.180130005,
            -12.006711960,  113.528327942,  -20.000000000),
            (0.624206722,   -0.776324987,    0.087679751,
            -0.753706276,   -0.627924204,   -0.193991154,
            0.205658153,    0.055005763,   -0.977073669,
            0.000912417,    0.000822168,  -38.829242706,
            144.630264282,  132.983993530,  115.019546509,
            -3.563867569,   81.349311829,  -20.000000000),
            (0.968158007,   -0.250144631,    0.009050513,
            -0.250012606,   -0.968132794,   -0.013895181,
            0.012239689,    0.011192130,   -0.999857962,
            -0.001026906,    0.003676833,  -35.612518311,
            143.959976196,  133.969512939,  114.222045898,
            -20.115106583,   91.698898315,  -20.000000000 ),
            (0.630304396,    0.764664412,    0.134201080,
            0.775923669,   -0.614777684,   -0.141375616,
            -0.025603920,    0.193240792,   -0.980812788,
            -0.000249907,   -0.000595212,  -40.734756470,
            142.911941528,  134.563812256,  113.872619629,
            -9.208276749,   90.970481873,  -20.000000000 ),
            
            # side (single or multi-model)
            (0.853055477,   -0.190655425,    0.485733509,
            -0.521587849,   -0.338814467,    0.783033252,
            0.015285396,   -0.921327472,   -0.388475388,
            -0.000381979,    0.002122164, -176.534942627,
            134.342285156,  132.092651367,  110.693161011,
            128.609085083,  224.218704224,  -20.000001907 ),
            
            # zoom (single model)
            (0.835139990,    0.537639856,   -0.116132051,
            0.545744836,   -0.836256802,    0.053179424,
            -0.068523534,   -0.107785024,   -0.991804361,
            -0.001693513,    0.002039829,  -69.380439758,
            142.934219360,  132.703445435,  115.728683472,
            35.292522430,  104.049064636,  -20.000000000 )
        ]
        
    def _load(self):
        """Load structures in `self._paths`"""
        if len(self._paths) < 1:
            raise ValueError("No file paths defined.")
        
        model_names = [] 
        for p in self._paths:
            n = os.path.basename(p)[:-4]
            model_names.append(n)
                
            cmd.load(p)    

        self._modelNames = model_names
        
    def _align(self, reference="7nmn"):
        
        if len(self._modelNames) < 1:
            raise ValueError("Model names are not defined.")
        
        if reference not in self._modelNames:
            raise ValueError("Reference < %s > is not a valid model name" % reference)
        
        print("Aligning to reference: < %s >" % reference)
        cmd.select("reference", "m. %s" % reference)
        
        # 6uqf                 RMSD =    3.842 (13168 atoms)
        # 7np3_mHCN2           RMSD =    2.223 (12992 atoms)
        for m in self._modelNames:        
            if m == reference:
                continue 
            
            r = cmd.align("model %s" % m, "reference")[0]
            print("RMSD for %s:     %.4f" % (m, r))
    
    def _multiColour(
        self, clrs=["gray75", "orange", "slate"], 
        models=["7nmn", "7np3_mHCN2", "6uqf"], 
    ):
        
        cmd.bg_color("white")
        
        if len(models) < 1:
            models = self._modelNames
        
        for i, m in enumerate(models):
            if i == len(clrs) - 1:
                print("Number of colours provided < number of models. Colouring stopped.")
                break 
            
            cmd.color(clrs[i], "m. %s" % m)
            
    
    def _select(self, idx, radius, chain="A"):
        """Select central residue with index `idx` and residues within `radius` of it

        Args:
            center_idx (int): index of central residue
            radius (float or int): distance in Angstroms
            chain (str): chain to select
        """

        # main residue 
        _root = "c. %s and i. %d" % (chain, idx)
        cmd.select("root", _root)

        # select nearby residues
        _nearby1 = "(not c. %s)" % chain
        _nearby2 = "(c. %s and i. 1-%d)" % (chain, idx-3)
        _nearby = "%s or %s" % (_nearby1, _nearby2)
        
        # select residues within `radius` of `root`
        cmd.select("nearby", "br. (%s) within %.1f of root" % (_nearby, radius))
        # select backbone carbons
        cmd.select("nearby", "nearby and not name CA")
        
        cmd.show("sticks", "root")
        cmd.show("sticks", "nearby")
        
        
    # saved view 
    def set_view(self, num):
        cmd.set_view(self._views[int(num)])
                
def hide(how=0)
    
    # multi close up 
    if how == 0:
        to_hide = ["c. C or c. D", "i. 365-411"]
    
    else:
        to_hide = ["c. B or c. D"]
    
    to_hide.extend([
        # remove all but s5 and s6 of c. B
        "c. B and i. 1-412",
        "c. B and i. 434-999",
        
        "i. 1-291",     # S1-S3
        "i. 460-999"    # CTD 
    ])
    
    for t in to_hide:
        cmd.remove(t)
    
    if how > 1:
        cmd.color("lightblue", "reference")
        cmd.color("ruby", "nearby")

def save(fname, width="8in", height="9in", path=""):
    
    if not os.path.isdir(path):
        path = r"C:/Users/delbe/Downloads/wut/wut/Post_grad/UBC/Research/lab/structure files/figures/"
    
    out = path + fname + ".png"
        
    cmd.png(out, width, height, 300, 1)
    print("Figure saved to < %s >" % out)
