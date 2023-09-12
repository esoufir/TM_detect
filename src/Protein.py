import argparse
from Bio.PDB import *
import numpy as np
import pymol
import copy

import Vector
import sys


class AminoAcid:
    """
    A class to represent a amino acid.

    ...

    Attributes
    ----------
    id : int
        Position in sequence
    code : str
        3 letters code
    is_hydrophobic : bool
        
    point : Point
        Coordinates of the c-alpha

    
    Methods
    -------
    find_best_axis():
        Prints the person's name and age.
    """

    hydrophobics_amino_acids = ['PHE', 'GLY', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR']
    
    def __init__(self, code, id_aa, x, y, z):
        self.id = id_aa
        self.code = code
        self.is_hydrophobic = False
        if self.code in AminoAcid.hydrophobics_amino_acids:
            self.is_hydrophobic = True
        self.point = Vector.Point(x, y, z)

    def __str__(self):
        """Redifining print() comportement."""
        return f"AA n°{self.id} is {self.code} at ({self.point.get_x():.3f}," \
               f"{self.point.get_y():.3f},{self.point.get_z():.3f}) \n"


class Protein:
    """
    A class to represent a Protein with one chain.

    ...

    Attributes
    ----------
    name : str
        Name of the protein
    mass_center : Point
        Coordinates of the mass center
    amino_acid_sequence : list <AminoAcid>
        List of the amino acids accessible to solvant
    full_sequence : list <AminoAcid>
        Full sequence
    best_positions : list <Axis>
        List of every axis with its own best position 
    
    Methods
    -------
    find_best_axis():
        Prints the person's name and age.
    """

    def __init__(self, name):
        self.name = name
        self.mass_center = Vector.Point(0, 0, 0)
        self.amino_acid_sequence = []
        self.full_sequence = []
        self.best_positions = []

    def find_best_axis(self):
        """
        Finds the best axis among the list self.best_positions. 
        The best axis is the one with the best value of hydrophobicity. 

        Returns
        -------
        best_axis : Axis
            Axis with the highest value of hydrophobicity
        """
        best_axis_val = 0
        best_axis_found = None
        for axis in self.best_positions:
            if axis.best_hydrophobicity > best_axis_val:
                best_axis_val = axis.best_hydrophobicity
                best_axis_found = copy.deepcopy(axis)
        return best_axis_found

    def compute_mass_center(self, chain):
        """
        Compute the coordinates of the mass center.
        
        Parameters
        ----------
        chain : Chain (Bio.PDB)
            Chain to select. 
        
        Returns
        -------
        best_axis : Point
            Coordinates of the mass center
        """
        result = chain.center_of_mass()
        x, y, z = result[0], result[1], result[2]
        self.mass_center = Vector.Point(x, y, z)


def check_input_file(input_file):
    """
    Does some verifications on the input file

        Parameter:
        ----------
            input_file : str
                Path to the input PDB file. 
        Returns:
        -------
            None
    """
    if not input_file.lower().endswith(".pdb"):
        print("Wrong extension. A PDB file is required.")
        sys.exit(0)
    try:
        # Attempt to parse the PDB file
        parser = PDBParser()
        structure = parser.get_structure("temp", input_file)
        # Check if the structure contains at least one chain labeled as "A"
        for model in structure:
            for chain in model:
                if chain.id == 'A':
                    return True
    except Exception:
        print("Issue when attempting to read the PDB file. Please check the content of your input file.")
        sys.exit(0)


def caculate_solvant_accessibility(structure, input_file):
    """
    Calculate the solvent accessibility of residues using DSSP (https://swift.cmbi.umcn.nl/gv/dssp/).

        Parameter:
        ----------
            input_file : str
                Path to the input PDB file. 
        Returns:
        -------
            None
    """
    print("Computing solvant accessibility...")
    model = structure[0]
    dssp = DSSP(model, input_file, dssp='mkdssp')
    return dssp


def parse_pdb(input_file,chain):
    """
    Parse the PDB file to construct a Protein object with the
    useful informations.

        Parameter:
        ----------
            input_file : str
                Path to the input PDB file. 
            chain : str
                Chain to consider. 
        Returns:
        -------
            protein : Protein
                A Protein initialised with the infos parsed in the PDB.
    """
    print("Parsing the PDB file ...")
    p = PDBParser()

    structure = p.get_structure("structure", input_file)

    id_amino_acid = 1  # Counting the accessible residues
    id_full_amino_acid = 1  # Counting all residues
    protein = Protein(name=input_file[8:-4])
    # Compute the solvant accessibility and set it for each amino acid of the protein
    dssp_res = caculate_solvant_accessibility(structure, input_file=input_file)
    print("Trimming to get only exposed residues...")
    exposed_residues = []
    # Getting only the residue id and its solvant accessibility (0.3)
    for res_id, aa, bb, asa, _, _, _, _, _, _, _, _, _, _ in dssp_res:
        if asa > 0.3:
            exposed_residues.append(res_id)
    print(f"Found {len(exposed_residues)} exposed residues.")
    # Getting only one chain
    model = structure[0]
    chain_selected = model["A"]

    # Calculating mass_center of the chain :
    protein.compute_mass_center(chain_selected)
    print("Mass center is", protein.mass_center)

    # Creating the residues of the protein:
    for residue in chain_selected:
        if residue.has_id("CA"):
            atom = residue['CA']
            x, y, z = atom.get_coord()
            # Creating and setting the AminoAcid object
            new_amino_acid = AminoAcid(code=residue.get_resname(), id_aa=residue.get_id()[1], x=x, y=y, z=z)
            # Linking the amino acid to its protein object
            protein.full_sequence.append(new_amino_acid)
            id_full_amino_acid += 1
            # If the residue is exposed, adding it to a dedicated list :
            if residue.get_id()[1] in exposed_residues:
                if residue.has_id("CA"):
                    atom = residue['CA']
                    x, y, z = atom.get_coord()
                    # Creating and setting the AminoAcid object
                    new_amino_acid = AminoAcid(code=residue.get_resname(), id_aa=residue.get_id()[1], x=x, y=y, z=z)
                    # Linking the amino acid to its protein object
                    protein.amino_acid_sequence.append(new_amino_acid)
                    id_amino_acid += 1
    return protein


def show_in_pymol(plane1, plane2, pdb_file, mass_center):
    """
    Use PyMol GUI to display the protein with the predicted membrane.

        Parameter:
        ----------
            pdb_file : str
                Path to the input PDB file. 
            plane1 : Plane
                Upper plane of the membrane.
            plane2 : Plane
                Lower plane of the membrane.
            mass_center : Point
                Coordinates of the mass center of the protein.

        Returns:
        -------
            None
    """
    # Run in quiet mode to avoid problem with argument parsing
    pymol.finish_launching(['pymol', '-q'])
    # Load the PDB file
    pymol.cmd.load(pdb_file, "protein")  # TODO
    pymol.cmd.remove("solvent")

    # Determine the dimensions of the molecule in the X-axis
    x_min, x_max = float("inf"), float("-inf")
    for atom in pymol.cmd.get_model("protein").atom:
        x = atom.coord[0]
        if x < x_min:
            x_min = x
        if x > x_max:
            x_max = x

    # Determine the dimensions of the molecule in the Y-axis
    y_min, y_max = float("inf"), float("-inf")
    for atom in pymol.cmd.get_model("protein").atom:
        y = atom.coord[1]
        if y < y_min:
            y_min = y
        if y > y_max:
            y_max = y

    # If point_x is not provided, use the midpoint between x_min and x_max
    center_x = (x_min + x_max) / 2.0

    # Define the step size for sampling points
    step = 3
    
    # Generate points on plane 1
    points_on_plane1 = []
    for x in np.arange(center_x - (x_max - x_min) / 2, center_x + (x_max - x_min) / 2 + step, step):
        for y in np.arange(y_min, y_max + step, step):
            # Calculate z coordinate using the plane equation for plane 1
            z1 = (-plane1.a * x - plane1.b * y - plane1.d) / plane1.c
            points_on_plane1.append((x, y, z1))

    # Generate points on plane 2
    points_on_plane2 = []
    for x in np.arange(x_min, x_max + step, step):
        for y in np.arange(y_min, y_max + step, step):
            # Calculate z coordinate using the plane equation for plane 2
            z2 = (-plane2.a * x - plane2.b * y - plane2.d) / plane2.c
            points_on_plane2.append((x, y, z2))

    # Create pseudoatoms for points on plane 1
    for idx, point in enumerate(points_on_plane1):
        x, y, z = point
        atom_name = f"plane1_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=[x, y, z], color="yellow")
        pymol.cmd.show("spheres", f"plane1_{idx}")

    # Create pseudoatoms for points on plane 2
    for idx, point in enumerate(points_on_plane2):
        x, y, z = point
        atom_name = f"plane2_{idx}"
        pymol.cmd.pseudoatom(atom_name, pos=[x, y, z], color="yellow")
        pymol.cmd.show("spheres", f"plane2_{idx}")

    # Mass center
    if mass_center is not None:
        atom_name = "mass_center"
        pymol.cmd.pseudoatom(atom_name, pos=[mass_center.get_x(), mass_center.get_y(), mass_center.get_z()],
                             color="magenta")
        pymol.cmd.show("spheres", atom_name)

    # Show the protein structure
    pymol.cmd.show("cartoon", "protein")


def optimizing_width_membrane(gap_membrane, axis_init, amino_acid_sequence,
                              plane_to_consider):  # plane to consider is plane1 or plane 2
    """
    Optimise the mebrane width by exploring a little above/below a plane.

        Parameter:
        ----------
            gap_membrane : int
                Size of sliding to adjust the membrane.
            axis_init : Axis
                Axis to adjust.
            amino_acid_sequence : List <AminoAcid>
                Residues exposed.
            plane_to_consider : int
                Either to consider plane 1 or plane 2.

        Returns:
        -------
            best_axis : Axis
                Best axis with the planes maximizing hydrophobicity.
    """
    # Keeping in mind what was the best axis with the best planes yet :
    best_axis = copy.deepcopy(axis_init)
    # First slide to be able to compare two diffent axis : 
    if plane_to_consider == 1:
        # Plane 1 is slided
        axis_init.plane1.slide_plane(gap_membrane)
    else:
        # Plane 2 is slided
        axis_init.plane2.slide_plane(gap_membrane)

    while axis_init.explore_axe_bis(amino_acid_sequence, ref=best_axis):
        if plane_to_consider == 1:
            # Plane 1 is slided
            axis_init.plane1.slide_plane(gap_membrane)
        else:
            # Plane 2 is slided
            axis_init.plane2.slide_plane(gap_membrane)
    return best_axis
