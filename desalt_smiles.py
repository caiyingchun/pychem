"""
Removes common salts and solvents from SMILES strings.

This script will remove common salts and solvents from SMILES strings using
an internal library. SMILES, SDF and Maestro input and output formats can be
used. However, the script will internally convert all structures to canonical
SMILES prior to removing the salt.

If a salt(s) is detected and the salt(s) is unknown, the largest molecule in
the SMILES string will be returned as the desalted compound. If a single
entry contains duplicate molecules, the first molecule will be retained. If
writing out to either a Maestro or SDF file, the structure will most likely
be distorted and so should be prepared with LigPrep.

Copyright Schrodinger, LLC. All rights reserved.
"""

import argparse
import os
import sys

from schrodinger import structure
from schrodinger.infra import canvas
from schrodinger.structutils import build
from schrodinger.structutils import smiles
from schrodinger.utils import cmdline
from schrodinger.utils import fileutils

###############################################################################
# Globals and Constants

FORMATS = [fileutils.SMILES, fileutils.SD, fileutils.MAESTRO]


###############################################################################
def salt_set():
    """
    Create a list of all salt SMILES strings.

    :return:  Set of all salt SMILES strings.
    :rtype:  set
    """

    anion_salts = {
        "O=C(C)O",  # Acetic acid
        "[O-]C(=O)C",  # Acetic acid
        "O=C(Oc1ccccc1C(=O)O)C",  # 2-Acetoxybenzoic acid*
        "O=C(O)C12CC3CC(C1)CC(C2)C3",  # Adamantanecarboxylic acid
        "O=C(O)CCCCC(=O)O",  # Adipic acid
        "O=C(O)c1ccc(cc1O)N",  # 4-Aminosalicylic acid*
        "OCC(O)C(OC1=O)C(O)=C1O",  # Ascorbic acid
        "O=C(O)CC(N)C(=O)O",  # Aspartic acid*
        "O=C(O)CCCCCCCC(=O)O",  # Azelaic acid*
        "O=S(=O)(O)c1ccccc1",  # Benzenesulfonic acid
        "[O-]S(=O)(=O)c1ccccc1",  # Benzenesulfonic acid
        "O=C(O)c1ccccc1",  # Benzoic acid
        "[O-]C(=O)c1ccccc1",  # Benzoic acid
        "[O-]C(=O)C(C)C",  # Butyric acid
        "O=C([O-])O",  # Carbonate
        "OC(=O)O",  # Carbonate
        "[O-]C([O-])=O",  # Carbonate
        r"O=C(O)\C=C\c1ccccc1",  # Cinnamic acid*
        "O=C(O)CC(O)(C(=O)O)CC(=O)O",  # Citric acid
        "O=C(O)C1CCCCC1",  # Cyclohexancarboxylic acid*
        "O=C(O)CCCCCCCCC",  # Decanoic acid*
        "O=C(O)CCCCCCCCCCC",  # Dodecanoic acid*
        "O=S(=O)(OCCCCCCCCCCCC)O",  # Dodecylsulfuric acid*
        "O=S(=O)(O)CCS(O)(=O)=O",  # Ethane-1,2-disulfonic acid*
        "CCOS([O-])(=O)=O",  # Ethylsulfuric acid
        "CCOS(=O)(=O)O",  # Ethylsulfuric acid
        "CCO[S+]([O-])([O-])=O",  # Ethylsulfuric acid
        "O=C[O-]",  # Formic acid
        "OC=O",  # Formic acid
        "O=C(O)/C=C/C(=O)O",  # Fumaric acid
        "[O-]C(=O)/C=C/C(=O)O",  # Fumaric acid
        "[O-]C(=O)/C=C/C([O-])=O",  # Fumaric acid
        "[O-]C(=O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)CO",  # Gluconic acid
        "[O-]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",  # Gluconic acid
        "O=C(O)CCC(N)C(=O)O",  # Glutamic acid*
        "O=C(O)CO",  # Glycolic acid
        "F",  # Halogens
        "[F-]",  # Halogens
        "Cl",  # Halogens
        "[Cl-]"  # Halogens
        "Br",  # Halogens
        "[Br-]",  # Halogens
        "I",  # Halogens
        "[I-]",  # Halogens
        "O=S(=O)(O)CCO",  # 2-Hydroxyethanesulfonic acid
        r"OC(=O)\C=C(\O)C(O)=O",  # Hydroxymaleic acid*
        "[O-]C(=O)c1ccncc1",  # Isonicotinic acid
        "O=C(O)C(C)O",  # Lactic acid
        "[O-]C(=O)[C@H](C)O",  # Lactic acid
        "[O-]C(=O)[C@@H](C)O",  # Lactic acid
        r"O=C(O)/C=C\C(=O)O",  # Maleic acid
        "[O-]C(=O)CC(O)C(=O)O",  # Malic acid
        "[O-]C(=O)C(O)CC(=O)O",  # Malic acid
        "O=C(O)[C@@H](O)CC(=O)O",  # Malic acid
        "[O-]C(=O)CC(=O)O",  # Malonic acid
        "O=C(O)CC(=O)O",  # Malonic acid
        "O=C(O)C(O)c1ccccc1",  # Mandelic acid*
        "Cc1c(S(=O)(=O)O)c(C)cc(C)c1",  # Mesitylenesulfonic acid
        "CS(=O)(=O)O",  # Methansulfonic acid
        "C[S-](=O)(=O)=O",  # Methansulfonic acid
        "CS([O-])(=O)=O",  # Methansulfonic acid
        "O=S(=O)(O)c1ccccc1C",  # 2-Methylbenzenesulfonic
        "O=S(=O)(O)c1cccc(c1)C",  # 3-Methylbenzenesulfonic acid*
        r"O=C(O)\C=C(/C(=O)O)C",  # Methylmaleic acid*
        "COS([O-])(=O)=O",  # Methylsulfuric acid
        "COS(=O)(=O)O",  # Methylsulfuric acid
        "O=S(=O)(O)c1cccc2c1cccc2S(=O)(=O)O",  # 1,5-Naphthalene-disulfonic acid*
        "c1cccc(c12)ccc(c2)S([O-])(=O)=O",  # 2-Naphthalenesulfonic acid
        "O=S(=O)(O)NC1CCCCC1",  # N-cyclohexylsulfamic acid*
        "O=S(=O)(O)NCC",  # N-ethylsulfamic acid*
        "[O-]C(=O)c1cccnc1",  # Nicotinic acid*
        "O=C(O)c1cccnc1",  # Nicotinic acid*
        "O=[N+]([O-])O",  # Nitrate
        "O=N(=O)O",  # Nitrate
        "[O-][N+]([O-])=O",  # Nitrate
        "O=N([O-])=O",  # Nitrate
        "O=S(=O)(O)NC",  # N-methylsulfamic acid*
        "O=C(O)c2ccccc2Oc1ccccc1",  # 2-Phenoxybenzoic acid*
        "CCCNS(O)(=O)=O",  # N-propylsulfamic acid*
        "O=C(O)CCCCCCC",  # Octanoic acid
        "[O-]C(=O)c1cc(=O)[nH]c(=O)[nH]1",  # Orotic acid
        "O=C(O)c1cc(=O)[nH]c(=O)[nH]1",  # Orotic acid
        "O=C1CC(C(=O)O)NC(=O)N1",  # Orotic acid
        "O=C(O)C(=O)O",  # Oxalic acid
        "[O-]C(=O)C(=O)O",  # Oxalic acid
        "[O-]C(=O)C([O-])=O",  # Oxalic acid
        "O=C(O)C(=O)CCC(=O)O",  # Oxoglutaric acid
        "[O-][Cl](=O)(=O)=O",  # Perchlorate
        "O=[Cl](=O)(=O)O",  # Perchlorate
        "O=[Cl-](=O)(=O)=O",  # Perchlorate
        "O=C(O)Cc1ccccc1",  # Phenylacetic acid
        "O=P(O)(O)O",  # Phosphate
        "O=[P-](O)(O)O",  # Phosphate
        "[O-]P(=O)(O)O",  # Phosphate
        "O=[P](=O)(=O)O",  # Phosphate
        "[O-]P([O-])(=O)O",  # Phosphate
        "[O-][N+](=O)c(c1O)cc([N+]([O-])=O)cc1[N+]([O-])=O",  # Picric acid
        "[O-]c1c([N+]([O-])=O)cc([N+]([O-])=O)cc1[N+]([O-])=O",  # Picric acid
        "O=C(O)CCCCCC(=O)O",  # Pimelic acid*
        "O=C(O)CC",  # Propionic acid*
        "O=C(O)c1ccccc1C(=O)O",  # Phthalic acid*
        "O=C(O)C(N1)CCC1=O",  # Pyroglutamic acid
        "O=C(O)[C@@H](N1)CCC1=O",  # Pyroglutamic acid
        "O=C(O)c1c(O)cccc1",  # Salicylic acid
        "[O-]C(=O)c1c(O)cccc1",  # Salicylic acid
        "O=C(O)CCCCCCC(=O)O",  # Suberic acid*
        "O=C(O)CCC(=O)O",  # Succinic acid
        "[O-]C(=O)CCC(=O)O",  # Succinic acid
        "O=S(=O)(O)N",  # Sulfamic acid*
        "Nc1ccc(cc1)S(=O)(=O)O",  # Sulfanilic acid
        "O=S(=O)(O)O",  # Sulfuric acid
        "[O-]S(=O)(=O)O",  # Sulfuric acid
        "[O-]S([O-])(=O)=O",  # Sulfuric acid
        "O=[S-2](=O)(=O)=O",  # Sulfuric acid
        "[O-][S](=O)(=O)=O",  # Sulfuric acid
        "[O-][S+](=O)(O)O",  # Sulfuric acid
        "[O-][S+]([O-])([O-])=O",  # Sulfuric acid
        "[O-][S+]([O-])(=O)O",  # Sulfuric acid
        "[O-]C(=O)C(O)C(O)C(=O)O",  # Tartaric acid
        "O=C(O)C(O)C(O)C(=O)O",  # Tartaric acid
        "O=C(O)[C@H](O)[C@@H](O)C(=O)O",  # Tartaric acid
        "O=C(O)[C@@H](O)[C@H](O)C(=O)O",  # Tartaric acid
        "O=C(O)[C@@H](O)[C@@H](O)C(=O)O",  # Tartaric acid
        "[O-]C(=O)[C@H](O)[C@@H](O)C(=O)O",  # Tartaric acid
        "F[B-](F)(F)F",  # Tetrafluroborate
        "Cc1ccc(cc1)S(=O)(=O)O",  # p-toluenesulfonic acid
        "Cc1ccc(S([O-])(=O)=O)cc1",  # p-toluenesulfonic acid
        "Cc1ccc([S+]([O-])([O-])=O)cc1",  # p-toluenesulfonic acid
        "O=C(C(F)(F)F)O",  # TFA
        "[O-]C(=O)C(F)(F)F",  # TFA
        "SC#N",  # Thiocyanate
        "S=C=[N-]",  # Thiocyanate
        "[S-]C#N"
    }  # Thiocyanate

    cation_salts = {
        "[NH4+]",  # Ammonia
        "N",  # Ammonia
        "O=C(O)C(N)CCCNC(=N)N",  # Arginine
        "NC1CC1",  # Cyclopropylamine
        "NCCCN",  # Diaminopropane
        "N1(C)CCN(C)CC1",  # N.N'-dimethylpiperazine
        "C1CCCCC1NC2CCCCC2",  # Dicyclohexylamine
        "NCCN",  # Ethylenediamine
        "[Li+]",
        "[LiH]",
        "[Na+]",
        "[NaH]",
        "[Mg+2]",
        "[MgH+]",
        "[MgH2]",
        "[Al+3]",
        "[Al]",
        "[Al-]",
        "[SiH6+2]",
        "[K+]",
        "[KH]",
        "[KH2-]",
        "[Ca+2]",
        "[CaH+]",
        "[CaH2]",
        "[Sc+3]",
        "[Sc]",
        "[Ti+2]",
        "[Ti]",
        "[V]",
        "[Cr]",
        "[Cr+3]",
        "[Mn+2]",
        "[Mn]",
        "[Fe]",
        "[Fe+2]",
        "[Fe+3]",
        "[Fe+]",
        "[Co+2]",
        "[Co+3]",
        "[Co]",
        "[Ni+2]",
        "[Ni]",
        "[Cu+2]",
        "[Cu+]",
        "[Cu]",
        "[Zn+2]",
        "[ZnH2]",
        "[Ga]",
        "[Ga+3]",
        "[AsH6+3]",
        "[Rb+]",
        "[Rb]",
        "[Sr+2]",
        "[Y]",
        "[Y+3]",
        "[Zr+2]",
        "[Zr]",
        "[Mo]",
        "[RuH+]",
        "[RuH2+2]",
        "[Ru+2]",
        "[Ru]",
        "[Rh+]",
        "[Rh]",
        "[RhH]",
        "[PdH4]",
        "[PdH2+2]",
        "[Ag+]",
        "[Ag-]",
        "[Ag]",
        "[Cd+2]",
        "[CdH2]",
        "[In]",
        "[In+3]",
        "[SnH4]",
        "[SnH2+2]",
        "[Sb+3]",
        "[Cs+]",
        "[Cs]",
        "[Ba+2]",
        "[Ba]",
        "[La]",
        "[La+3]",
        "[Ce+3]",
        "[Ce]",
        "[Pr+3]",
        "[Pr]",
        "[Nd]",
        "[Nd+3]",
        "[Sm+3]",
        "[Sm]",
        "[Eu+3]",
        "[Eu]",
        "[Gd]",
        "[Gd+3]",
        "[Tb+3]",
        "[Tb]",
        "[Dy+3]",
        "[Dy]",
        "[Ho]",
        "[Ho+3]",
        "[Er+3]",
        "[Er]",
        "[Tm+3]",
        "[Tm]",
        "[Yb+3]",
        "[Lu]",
        "[Lu+3]",
        "[Ta]",
        "[W]",
        "[Ir+]",
        "[PtH4]",
        "[PtH2+2]",
        "[Au+]",
        "[Hg+2]",
        "[HgH2]",
        "[HgH+]",
        "[Tl+]",
        "[Tl]",
        "[PbH4]",
        "[PbH2+2]",
        "[Bi+3]",
        "[BiH3]",  # Elements
        "N1(CC)CCCCC1",  # N-ethylpiperidine*
        "O=C(O)C(N)CCCCN",  # Lysine
        "CC[N+](CC)(CC)CC"
    }  # TEA

    solvents = {
        "CC(=O)C",  # Acetone
        "N#CC",  # Acetonitrile*
        "[NH3+]c1ccccc1",  # Aniline
        "c1ccccc1",  # Benzene
        "OCCCC",  # 1-butanol*
        "OC(C)CC",  # 2-butanol*
        "O=C(C)CC",  # 2-butanone*
        "OC(C)(C)C",  # t-butyl alcohol*
        "N#Cc1ccccc1",  # Benzonitrile
        "OCCCCO",  # 1,4-Butanediol
        "ClC(Cl)(Cl)Cl",  # Carbon tetrachloride*
        "Clc1ccccc1",  # Chlorobenzene*
        "ClC(Cl)Cl",  # Chloroform*
        "C1CCCCC1",  # Cyclohexane
        "ClCCCl",  # 1,2-dichloroethane*
        "NC1CCCCC1",
        "[NH3+]C1CCCCC1",  # Cyclohexylamine
        "CC(=O)N(C)C",  # Dimethylacetamide
        "O=CN(C)C",  # Dimethylformamide
        "CCOCC",  # Diethylether
        "COCCOCCOC",  # Diglyme*
        "COCCOC",  # DME*
        "O(C)C",  # Dimethylether*
        "O=CN(C)C",  # DMF*
        "C1COCCO1",  # 1,4-Dioxane
        "CS(=O)C",  # DMSO
        "OCC",  # Ethanol
        "OCCO",  # Ethylene glycol*
        "OCC(O)CO",  # Glycerin*
        "CCCCCCC",  # Heptane
        "O=P(N(C)C)(N(C)C)N(C)C",  # HMPA*
        "N(P(N(C)C)N(C)C)(C)C",  # HMPT*
        "CCCCCC",  # Hexane*
        "NCCO",
        "[NH3+]CCO",  # MEA
        "CO",  # Methanol
        "ClCCl",  # Methylene Chloride*
        "CN1CCOCC1",
        "C[NH+]1CCOCC1",  # N-Methylmorpholine
        "C1COCCN1",
        "C1COCC[NH2+]1",  # Morpholine
        "O(C(C)(C)C)C",  # MTBE*
        "[O-][N+](=O)C",  # Nitromethane*
        "O=C1N(C)CCC1",  # NMP*
        "CCCCC",  # Pentane*
        "Cc1ccncc1",  # 4-Picoline
        "C1CNCCN1",  # Piperazine
        "C1CC[NH2+]CC1",  # Piperidine
        "C1CCNCC1",  # Piperidine
        "OCCCO",  # 1,3-Propanediol
        "CCCO",  # 1-Propanol
        "CC(C)O",  # 2-Propanol
        "c1cc[nH+]cc1",  # Pyridine
        "c1ccncc1",  # Pyridine
        "C1CCOC1",  # THF
        "Cc1ccccc1",  # Toluene
        "CCN(CC)CC",  # Triethylamine
        "CC[NH+](CC)CC",  # Triethylamine
        "O",  # Water
        "[OH-]",  # Water
        "Cc1c(C)cccc1",  # Xylene (ortho)
        "Cc1cc(C)ccc1",  # Xylene (meta)
        "Cc1ccc(C)cc1"
    }  # Xylene (para)

    misc = {"[H][H]", "[H+]"}

    salts = anion_salts | cation_salts | solvents | misc

    return salts


###############################################################################
def parse_args():
    """
    Parse the command line options.

    :return:  All script arguments and options, the input and output file
        formats and the salt output file name.
    :rtype:  class:`argparse.Namespace`, str, str, str
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("infile",
                        help="Input file. Valid formats are SMILES, Maestro and"
                        " SDF formats.")

    parser.add_argument("outfile",
                        help="Desalted structure output file. Valid formats "
                        "are SMILES, Maestro or SDF format.")

    parser.add_argument("-name",
                        default="",
                        help="Field in a SD or Maestro file that can be used "
                        "to overwrite the name of molecules.")

    parser.add_argument("-neutralize",
                        action='store_true',
                        help="Neutralize the output structures.")

    args = parser.parse_args()

    if not os.path.isfile(args.infile):
        print("Input file %s cannot be found." % args.infile)
        sys.exit(1)

    in_format = fileutils.get_structure_file_format(args.infile)
    out_format = fileutils.get_structure_file_format(args.outfile)
    salt_outfile = fileutils.get_basename(args.outfile) + "_salts.usmi"

    if in_format not in FORMATS:
        print("Unsupported input file format.")
        sys.exit(1)

    if out_format not in FORMATS:
        print("Unsupported output file format.")
        sys.exit(1)

    if in_format == fileutils.SMILES and args.name:
        print("Cannot specify an alternate name field in a SMILES file.")
        sys.exit(1)

    if os.path.isfile(salt_outfile):
        os.remove(salt_outfile)

    return args, in_format, out_format, salt_outfile


###############################################################################
def canonical_smiles(split_smiles):
    """
    For each SMILES string, convert them to canonical SMILES and add them to
    a list. Then ensure all entries in the list are unique.

    :param split_smiles:  List of SMILES strings.
    :type split_smiles:   list

    :return:  list of canonical SMILES strings
    :rtype:  list
    """

    sg = canvas.ChmMmctSmilesGenerator()

    unique_smiles = [sg.canonicalize(smi) for smi in split_smiles]
    # Remove duplicates while preserving sort order:
    unique_smiles_list = list(dict.fromkeys(unique_smiles))

    return unique_smiles_list


###############################################################################
def remove_known_salts(unique_smiles_list, salts):
    """
    For each SMILES string in the list, test to see if it is in the list of
    known salts.
    Returns either a single unsalted SMILES or a list of unsalted SMILES along
    with the list of known salts -- either the final SMILES will be None or
    the list of unsalted SMILES will be empty.

    :param unique_smiles_list:  List of unique SMILES strings.
    :type unique_smiles_list:   list

    :param salts:  Set of known salts in SMILES format.
    :type salts:   set

    :return:  Final unsalted SMILES, SMILES containing unknown
        salts and a list of identified known salts.
    :rtype:  str, list, list
    """

    known_salts = []
    tmp_smiles = []
    final_smiles = None

    unsalted_smiles = None

    for smi in unique_smiles_list:
        if smi in salts:
            known_salts.append(smi)
        else:
            tmp_smiles.append(smi)

    if len(tmp_smiles) == 1:
        final_smiles = tmp_smiles[0]
    elif len(tmp_smiles) >= 1:
        unsalted_smiles = tmp_smiles

    return final_smiles, unsalted_smiles, known_salts


#########################################################################
def process_identified_salts(known_salts, title, salt_file):
    """
    For each identified salt SMILES, write the salts to a SMILES file.

    :param known_salts:  List of identified salt SMILES strings.
    :type known_salts:   list

    :param title:  Title of the input structure.
    :type title:   str

    :param salt_file:  File name for the file containing the identified salts.
    :type salt_file:   str
    """

    with open(salt_file, 'a') as id_salts:
        for salt in sorted(known_salts):
            id_salts.write(str(salt) + " " + str(title) + "\n")


###############################################################################
def retain_largest_molecule(unsalted_smiles):
    """
    For those structures where a salt cannot be identified, retain the
    largest molecule as the main compound and return the rest as salts.

    :param unsalted_smiles:  List of unsalted SMILES.
    :type unsalted_smiles:   list

    :return: Largest unsalted SMILES, Removed salts
    :rtype: str, list(str)
    """

    unsalted_smi = None
    removed_salts = []

    largest = None

    min_weight = 0

    for i, smi in enumerate(unsalted_smiles):
        mol = canvas.ChmMol_fromSMILES(smi)
        weight = round(mol.getMW(), 2)

        if weight > min_weight:
            largest = i
            min_weight = weight

    if largest is not None:
        unsalted_smi = unsalted_smiles.pop(largest)
        removed_salts = unsalted_smiles

    return unsalted_smi, removed_salts


###############################################################################
def desalt_smiles(salted_smiles, salts):
    """
    Get the desalted SMILES and all salts.

    :param salted_smiles:  salted SMILES.
    :type salted_smiles:   string

    :param salts:  known salts to search for
    :type salted_smiles:   list(str)

    :return:  Final unsalted SMILES, a list of known and unknown salts.
    :rtype:   str, list
    """

    split_smiles = salted_smiles.split(".")

    if len(split_smiles) == 1:
        # Only one molecule present
        return salted_smiles, []

    try:
        unique_smiles_list = canonical_smiles(split_smiles)
    except RuntimeError:
        # If the structures cannot be converted to unique SMILES, skip them
        return None, None

    final_smiles, unsalted_smiles, removed_salts = remove_known_salts(
        unique_smiles_list, salts)

    if unsalted_smiles is not None:
        final_smiles, removed_unknown_salts = retain_largest_molecule(
            unsalted_smiles)
        removed_salts += removed_unknown_salts

    return final_smiles, removed_salts


def main():
    """
    Main body of the script.
    """

    salts = salt_set()

    final_smiles = None

    sg = smiles.SmilesGenerator(smiles.STEREO_FROM_ANNOTATION_AND_GEOM,
                                unique=True)

    cmd_args, input_format, output_format, salt_file = parse_args()

    reader = structure.SmilesReader(cmd_args.infile) \
        if input_format == fileutils.SMILES else \
        structure.StructureReader(cmd_args.infile)

    writer = structure.SmilesWriter(cmd_args.outfile) \
        if output_format == fileutils.SMILES else \
        structure.StructureWriter(cmd_args.outfile)

    for i, st in enumerate(reader):
        salted_smiles = st.smiles if input_format == fileutils.SMILES else \
            sg.getSmiles(st)

        if cmd_args.name:
            title = st.property.get(cmd_args.name,
                                    "Untitled Molecule_" + str(i))
        else:
            title = st.title
            if not title:
                title = "Untitled Molecule_" + str(i)

        final_smiles, removed_salts = desalt_smiles(salted_smiles, salts)

        if final_smiles is None:
            # Error; skip structure
            continue

        if removed_salts:
            process_identified_salts(removed_salts, title, salt_file)

        smi_struct = structure.SmilesStructure(final_smiles)
        smi_struct.title = title

        if output_format is not fileutils.SMILES or cmd_args.neutralize:
            smi_struct = smi_struct.get2dStructure()

            if cmd_args.neutralize:
                smi_struct = build.neutralize_structure(smi_struct)

        writer.append(smi_struct)
    writer.close()


if __name__ == '__main__':
    cmdline.main_wrapper(main)
