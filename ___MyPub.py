# pymol publication quality
# https://bionerdnotes.wordpress.com/2018/11/12/getting-high-quality-pictures-in-pymol/

from pymol import cmd
import numpy as np
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
    "set dash_radius, 0.0500",
    "set dash_width, 3",
    
    "bg_color white", 
    
    "set sphere_transparency, 0.45",
    "set cartoon_highlight_color, grey75",
    "set antialias, 4",
    "set ray_trace_mode, 1",
    "set ray_trace_gain, 0.005",
    "set cartoon_discrete_colors, 1",
    "set ray_shadows, 0",
    
    "set fog, on",
    "set ray_trace_fog, 1"
]


def _reset():
    cmd.reinitialize("everything")
    cmd.reinitialize("original_settings")


class load_structures:
    def __init__(self, path="", names=[]):

        if not os.path.isdir(path):
            path = r"C:/Users/delbe/Downloads/wut/wut/Post_grad/UBC/Research/lab/structure files/pdb_files/HCN2_MOUSE_O88703_Potassiumsodium_hyperpolarizati/models/"

        if len(names) < 1:
            names = ["7nmn/7nmn", "5u6o/5u6o", "6uqf/6uqf",
                     "6uqg/6uqg", "7np3_mHCN2"]

        self.out_path = path
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
             14.193284035,   66.120056152,  -20.000000000),
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
             -20.115106583,   91.698898315,  -20.000000000),
            (0.630304396,    0.764664412,    0.134201080,
             0.775923669,   -0.614777684,   -0.141375616,
             -0.025603920,    0.193240792,   -0.980812788,
             -0.000249907,   -0.000595212,  -40.734756470,
             142.911941528,  134.563812256,  113.872619629,
             -9.208276749,   90.970481873,  -20.000000000),

            # side (single or multi-model)
            (0.853055477,   -0.190655425,    0.485733509,
             -0.521587849,   -0.338814467,    0.783033252,
             0.015285396,   -0.921327472,   -0.388475388,
             -0.000381979,    0.002122164, -176.534942627,
             134.342285156,  132.092651367,  110.693161011,
             128.609085083,  224.218704224,  -20.000001907),

            # zoom (single model)
            (0.835139990,    0.537639856,   -0.116132051,
             0.545744836,   -0.836256802,    0.053179424,
             -0.068523534,   -0.107785024,   -0.991804361,
             -0.001693513,    0.002039829,  -69.380439758,
             142.934219360,  132.703445435,  115.728683472,
             35.292522430,  104.049064636,  -20.000000000)
        ]

        self._load()
        self._align()

    def _load(self):
        """Load structures in `self._paths`"""
        if len(self._paths) < 1:
            raise ValueError("No file paths defined.")

        model_names = []
        for p in self._paths:
            n = os.path.basename(p)[:-4]
            model_names.append(n)

            print("Loading... ", n)
            cmd.load(p)

        self._modelNames = model_names

    def _align(self, reference="7nmn"):

        if len(self._modelNames) < 1:
            raise ValueError("Model names are not defined.")

        if reference not in self._modelNames:
            raise ValueError(
                "Reference < %s > is not a valid model name" % reference)

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
                print(
                    "Number of colours provided < number of models. Colouring stopped.")
                break

            cmd.color(clrs[i], "m. %s" % m)

    def unpack_atom_info(self, ids):
        """
        Read chain, resn, etc. for a list of atom indices
        https://pymolwiki.org/index.php/Iterate
        """

        _names = []

        def _add(chain, resn, resi, name):
            _names.append("%s %s %s %s" % (chain, resn, resi, name))

        _space = {"_add": _add}

        if isinstance(ids, str):
            cmd.iterate(ids, "_add(chain, resn, resi, name)", space=_space)
        else:
            for a in ids:
                cmd.iterate("%s`%d" %
                            a, "_add(chain, resn, resi, name)", space=_space)

        return _names

    def write_dists(self, model, root, sel):
        """
        Get distance between `root` and `sel`, and write the distances to a text file 
        https://pymolwiki.org/index.php/Distance

        model: str = name of the current model (primarily for bookkeeping)
        root: str = string corresponding to a PyMol selection object 
        sel: str = string corresponding to another PyMol selection object

        Returns 
            mean: float = mean minimum distance 
            sd: float = standard deviation of distance 
        """

        # iterator over atoms
        # e.g.  cmd.index: (7np3_mHCN2`14248)
        atoms1 = []
        dists1 = []
        for at1 in cmd.index(sel):
            dists2 = []
            atoms2 = []

            for at2 in cmd.index(root):
                _dx = cmd.get_distance(atom1="%s`%d" %
                                       at1, atom2="%s`%d" % at2)
                dists2.append(_dx)
                atoms2.append(at2)

            i = np.argmin(dists2)
            atoms1.append(atoms2[i])
            dists1.append(dists2[i])

        # compute summary statistics on distances
        mu = np.mean(dists1)
        sd = np.std(dists1)
        header2 = "\n \n Average distance: \nmin = %.3f \nmu = %.3f \nsd = %.3f" % (
            min(dists1), mu, sd)
        print(header2)

        # convert distances to strings to add to output txt file
        dists1 = ["%.2f" % d for d in dists1]

        atoms1 = self.unpack_atom_info(atoms1)
        atoms2 = self.unpack_atom_info(sel)

        f_out = self.out_path + model + "_dists.txt"
        with open(f_out, "w") as io:
            header = "Distances for model < %s > and root < %s > and selection < %s >\n" % (
                model, root, sel)
            io.write(header)
            io.write(header2 + "\n")

            io.write("Columns: Root atom    Selection atom  Distance (A)\n")

            for i, atom1 in enumerate(atoms1):

                # despite explicitly telling the selection to exclude all residues 400-460,
                # 431 somehow slips through in the contacts, so we need to
                # filter out intra-residue contacts here
                if "431" in atoms2[i]:
                    continue

                line = "  ".join([atom1, atoms2[i], dists1[i]])
                io.write(line + "\n")
        print("Successfully wrote distances to < %s > " % f_out)

        return mu, sd

    def _select(self, model, idx=431, radius=4.5, chain="A"):
        """Select central residue with index `idx` and residues within `radius` of it

        Args:
            center_idx (int): index of central residue
            radius (float or int): distance in Angstroms
            chain (str): chain to select
        """

        pre = "(m. %s) and " % model

        # main residue (side chain only)
        _root = pre + "c. %s and i. %d and sc." % (chain, idx)

        # select residues to exclude
        # roughly most of the S6 on the same chain
        ex = pre + ("c. %s and i. 400-460" % chain)

        sels = [model + s for s in ['_ex', '_root', '_nearby', '_nearbyFull']]
        cmd.select(sels[0], ex)
        cmd.select(sels[1], _root)
        cmd.select(sels[2], pre + "(not %s) within %.1f of (%s)" %
                   (sels[0], radius, sels[1]))

        # refine 'nearby' again to exclude 'root'
        cmd.select(sels[2], pre + "(%s) and (not %s)" % (sels[2], sels[1]))

        # expand 'nearby' selction to full residues
        cmd.select(model + "_nearbyFull", "br. (%s)" % sels[2])

        # save distances for all 'nearby' atoms and the atoms in 'root'
        # calculates the minimum distance between each atom in `nearby` and all atoms in `root`
        # i.e. finds the closest `root` atom for each `nearby` atom
        # mu = mean of minimum distances, sd = standard deviation
        self.write_dists(model, sels[1], sels[2])

        # the above doesn't generate distance objects, so we do so here
        cmd.distance("dx_%s" % model, sels[1], sels[2], radius)
        # remove all labels
        cmd.hide("labels", "all")

        # show full residues
        cmd.show("sticks", sels[3])
        cmd.show("sticks", sels[1])

    def selectAll(self):

        for m in self._modelNames:
            self._select(m)

    # saved view
    def set_view(self, num):
        cmd.set_view(self._views[int(num)])


def hide(how=0):

    # zoom ins, S5-S6 c.A and S6 of c. B
    if how == 0:
        to_hide = [
            "c. C or c. D", 
            "i. 365-411",
            
            # remove all but s5 and s6 of c. B
            "c. B and i. 1-410",
            "c. B and i. 440-999",

            "i. 1-291",     # S1-S3
            "i. 460-999"    # CTD
        ]

    # full shots of c. A and c. B
    elif how == 1:
        to_hide = ["c. C or c. D"]
    
    # truncation of c. A and c. B (S4-S6)
    elif how == 2:
        to_hide = []

    for t in to_hide:
        # cmd.remove(t)
        cmd.hide("cartoon", t)

    if how > 1:
        cmd.color("lightblue", "reference")
        cmd.color("ruby", "nearby")


def showBonds(model, dashes=True):

    pre = "(m. %s) and " % model
    
    cmd.show("sticks", pre + "(br. %s_root)" % model)
    cmd.show("sticks", pre + "%s_nearbyFull" % model)
    cmd.show("sphere", pre + "%s_nearby" % model)

    if dashes:
        cmd.show("dashes", "dx_%s" % model)
        
    
def presetViews(model, set_view=True):

    if model not in ['7nmn', '7np3_mHCN2', '5u6o', '6uqf', '6uqg']:
        raise ValueError("%s is not supported in presetViews" % model)

    print("Reminder: you may want to: \nset sphere_transparency=0.6")
    
    cmd.hide("all")
    
    pre = "(m. %s) and " % model
    if int(model[0]) < 7:
        # cmd.hide("cartoon", pre + "c. B")
        cmd.show("cartoon", pre + "(c. A or c. B)")
    else:
        cmd.show("cartoon", pre + "(c. A)")
    
    cmd.hide("everything", "not m. %s" % model)
    
    hide(0)        
    showBonds(model, dashes=True)
    cmd.color("black", "dx_%s" % model)

    if not set_view: return None 

    if model in ['7nmn', '7np3_mHCN2']:
        cmd.set_view(
            (0.735802114,    0.128497943,   -0.664881766,
            0.658835709,   -0.362820596,    0.658981502,
            -0.156554669,   -0.922942996,   -0.351636201,
            -0.003335670,    0.003044491,  -47.385192871,
            145.109497070,  130.834396362,  117.095565796,
            -61.299564362,  155.086624146,  -20.000000000 ))

    elif model == '5u6o':
        cmd.set_view((
            0.743370116,   -0.040746111,    0.667627156,
            -0.337418705,   -0.884657085,    0.321714938,
            0.577528596,   -0.464431494,   -0.671387017,
            0.000611171,    0.003140345,  -52.546115875,
            140.624206543,  133.207962036,  115.828216553,
            -229.189453125,  334.605041504,  -20.000000000 ))

    elif model == '6uqf':
        cmd.set_view(
            (0.883197486,    0.456724256,   -0.106543846,
             0.384107351,   -0.834776402,   -0.394442677,
             -0.269092858,    0.307454586,   -0.912718534,
             -0.001593671,    0.001918254,  -53.864448547,
             142.160278320,  135.679916382,  115.980529785,
             13.909035683,   94.075057983,  -20.000000000))

    elif model == '6uqg':
        cmd.set_view(
            (0.732726216,   -0.608215094,    0.305242926,
             -0.556425214,   -0.793694794,   -0.245789051,
             0.391773760,    0.010250090,   -0.920000672,
             0.001189526,    0.003080618,  -36.255870819,
             140.349334717,  136.380432129,  117.012245178,
             12.242701530,   60.759925842,  -20.000001907))

# subunit A-B in wide view
set_view (\
    -0.997289836,    0.059505343,   -0.043070592,\
     0.032249175,   -0.172066659,   -0.984544337,\
    -0.066005029,   -0.983280301,    0.169689313,\
    -0.000528872,    0.000868529, -380.889251709,\
   134.466247559,  138.314025879,  128.699386597,\
   319.289703369,  443.352264404,  -20.000000000 )

def wideView(model='5u6o', clrs=["0x1a3673", "0xd7e1f7"], view_ind=0):

    if model is not None:
        pre = "(m. %s) and " % model

        cmd.hide("cartoon", "all")
        cmd.show("cartoon", "m. %s and (c. A or c. B)" % model)
        cmd.hide("everything", "not m. %s" % model)
        
        # molecular details 
        showBonds(model, dashes=False)
        
        cmd.color(clrs[0], pre + "c. A")
        cmd.color(clrs[1], pre + "c. B")

    else:
        
        cmd.show("cartoon", "all")
        hide(1)
        
    if view_ind == 0:
        cmd.set_view((
            -0.940882266,    0.098684147,    0.324018717,\
            -0.333495766,   -0.102597132,   -0.937139690,\
            -0.059243105,   -0.989811838,    0.129451454,\
            0.001483327,    0.000707418, -372.652679443,\
            138.741653442,  139.680511475,  130.467956543,\
            332.032318115,  413.853546143,  -20.000000000 ))
    elif view_ind == 1:
        cmd.set_view((
            -0.986978233,    0.156574160,    0.036624614,
            -0.016419291,    0.128467515,   -0.991565764,
            -0.159967989,   -0.979271472,   -0.124220252,
            -0.000528872,    0.000868529, -350.829803467,
            134.466247559,  138.314025879,  128.699386597,
            289.230255127,  413.292816162,  -20.000000000 
        ))
    elif view_ind == 2:
        cmd.set_view((
            -0.997289836,    0.059505343,   -0.043070592,
            0.032249175,   -0.172066659,   -0.984544337,
            -0.066005029,   -0.983280301,    0.169689313,
            -0.000528872,    0.000868529, -380.889251709,
            134.466247559,  138.314025879,  128.699386597,
            319.289703369,  443.352264404,  -20.000000000 )
        )
    else:
        print("view_ind = 0-2")
    
def colorSegments(model='5u6o', hide_rest=True):
    
    clrs = {
        "S6" : [411, 442, '0x7de274'],
        "S5" : [333, 364, '0xc25b5b'],
        "S4" : [292, 327, '0x7f8acc']
    }
    
    if hide_rest:
        cmd.select("tmp", "m. %s and i. 292-411 and c. A" % model)
        cmd.hide("everything", "m. %s and not tmp" % model)
    
    pre = "(m. %s) and " % model 
    for k, v in clrs.items():
        cmd.select(k, pre + "i. %d-%d" % (v[0], v[1]))
        cmd.color(v[2], k)

def colorModels(show_dashes=True):
    clrs = ["0x6a51f6", "raspberry", "0x21db00", "0x9ad4d6", "wheat"]
    modelNames = ["5u6o", "6uqf", "6uqg", "7nmn", "7np3_mHCN2"]
    # clrs = {k : v for k, v in zip(modelNames, clrs)}
    
    cmd.show("cartoon", "all")
    hide(1)
        
    for model, clr in zip(modelNames, clrs):
        
        cmd.color(clr, "m. %s" % model)
        showBonds(model, dashes=show_dashes)
        cmd.color(clr, "dx_%s" % model)
        

def save(fname, width="8in", height="9in", path=""):

    if not os.path.isdir(path):
        path = r"C:/Users/delbe/Downloads/wut/wut/Post_grad/UBC/Research/lab/structure files/figures/"

    out = path + fname + ".png"

    cmd.png(out, width, height, 300, 1)
    print("Figure saved to < %s >" % out)


##############################

_reset()

u = load_structures()

u.selectAll()

hide(1)

print("Please set the following (you can just copy and paste the whole block into the terminal): ")
for x in _u:
    print(x)

print("Use `presetViews(...)` to engage preset views for a given model.")
print("Use `hide(0)` for zoom-ins and `hide(1)` for multi-subunit shots")
presetViews('7nmn')


"""
# view for zoom ins (4in x 4in) of 6uqf and 6uqg
set_view (\
     0.732726216,   -0.608215094,    0.305242926,\
    -0.556425214,   -0.793694794,   -0.245789051,\
     0.391773760,    0.010250090,   -0.920000672,\
     0.001189526,    0.003080618,  -52.275936127,\
   140.349334717,  136.380432129,  117.012245178,\
     2.082305908,  102.960418701,  -20.000001907 )
"""