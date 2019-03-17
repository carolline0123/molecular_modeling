"""Simulator of atoms movement and interactions

This script allows the user to generate simulation of atoms moving in
box with dimensions (30 A x 30 A x 30 A), when periodic boundary
conditions applied. Moreover, the formation of bonds between atoms is
possible with probability of 50% when distance between atoms does not
exceed 2 A. Default number of simulation frames was set to 1000.

As an input, script requires (.xyz) file in working directory. The atom
symbols are determined as 'C' (carbon element) by default.

This script generates file 'simulation.xyz' which contains coordinates
for every atom for every of 1000 frames. This file can be further opened
in visualisation program like VMD to enable observation of simulation.

This file can be also imported as a module and contains the following
functions:

    * read_coordinates - reads data from input file and returns
                         coordinates in the form of namedtuple
    * save_coordinates - saves data for every frame of simulation to the
                         output file
    * calc_distance - calculates the distance between atoms and checks
                      whether bond formation between them is possible
    * generate_positions - introduces shifts of atoms by generation of
                           random numbers
    * check_pbc - checks whether particular atom or set of atoms is not
                  outside periodic box
    * main - the main function of the script
"""

from collections import namedtuple
from math import sqrt
from random import uniform, choice

def read_coordinates(filename):
    """Reads data from .xyz file.

    Parameters
    ----------
    filename : str
        name of input file

    Returns
    -------
    namedtuple
        namedtuple of coordinates of atoms
    """

    with open(filename) as f:

        atom_number = int(f.readline())
        f.readline()
        coordinates = []

        for i in range(atom_number):

            line = f.readline().split()

            symbol = line[0]
            x = float(line[1])
            y = float(line[2])
            z = float(line[2])

            Atom = namedtuple(symbol, ['x', 'y', 'z'])
            atom = Atom(x, y, z)
            coordinates.append(atom)

        return coordinates

def save_coordinates(filename, coordinates, frames):
    """Saves coordinates of every atom in all frames to output file.

    Parameters
    ----------
    filename : str
        name of the output file
    coordinates : namedtuple
        namedtuple of coordinates of atoms
    frames : int
        a number of simulation frames
    """

    with open(filename, 'w') as f:

        coordinate_list = []
        atom_number = len(coordinates)
        f.write(str(atom_number))
        f.write('\n')

        for i in range(atom_number):
        	f.write('\n{} {:.3f} {:.3f} {:.3f}'.format(type \
                   (coordinates[i]).__name__, coordinates[i].x,
                    coordinates[i].y, coordinates[i].z))

        #conversion of namedtuple to list of coordinates
        for i in range(len(coordinates)):
            coordinate_list.append([coordinates[i].x,
                                    coordinates[i].y, coordinates[i].z])

        #creating empty dictionary for information which atoms are bonded
        #bonded = {key : value} = {atom_number : [list of atoms bonded to it]}
        bonded = {}

        while(frames):

            f.write('\n%s' %(str(atom_number)))
            f.write('\n')
            bonded = calc_distance(coordinate_list, atom_number, bonded)
            coordinate_list = generate_positions(coordinate_list, atom_number,
                                                 bonded)

            for i in range(atom_number):
                f.write('\nC {:.3f} {:.3f} {:.3f}'.format(coordinate_list[i][0],
                                                          coordinate_list[i][1],
                                                          coordinate_list[i][2])
                                                          )
            frames -= 1

def calc_distance(coordinate_list, atom_number, bonded):
    """Calculates distance between atoms and checkes whether bond
    formation is possible (when distance smaller than 2 A).

    Parameters
    ----------
    coordinate_list : list
        a list of lists of coordinates for every atom
    atom_number : int
        a number of atoms
    bonded : dict
        a dictionary of bonded atoms

    Returns
    -------
    bonded : dict
        a dictionary of bonded atoms
    """

    #Distance matrix formation. Value of 1 will be assigned to pairs of
    #atoms which are within the distance of 2 and value of 0 will be
    #assigned to ones too far from each other to form a bond.
    distmatrix = []

    for l in range(atom_number):
        distmatrixrow = []

        for m in range(atom_number):
            distance = sqrt((coordinate_list[m][0] - coordinate_list[l][0])**2
                          + (coordinate_list[m][1] - coordinate_list[l][1])**2
                          + (coordinate_list[m][2] - coordinate_list[l][2])**2)
            probability = choice([0, 1])

            if distance != 0:
                if (distance <= 2 and probability == 1):
                    distmatrixrow.append(1)
                else:
                    distmatrixrow.append(0)
            else:
                distmatrixrow.append(0)

        distmatrix.append(distmatrixrow)

    for m in range(atom_number):
        for n in range(atom_number):
            distmatrix[m][n] = distmatrix[n][m]

    #Creating dictionary is based on values for particular pairs of
    #atoms in distance matrix
    for m in range(atom_number):
        for n in range(atom_number):
            if distmatrix[m][n] == 1:
                if m + 1 in bonded:
                    bonded[m + 1].append(n + 1)
                else:
                    bonded[m + 1] = [n + 1]

    return bonded

def generate_positions(coordinate_list, atom_number, bonded):
    """Generates positions for atoms in particular frame.

    Parameters
    ----------
    coordinate_list : list
        a list of lists of coordinates for every atom
    atom_number : int
        a number of atoms
    bonded : dict
        a dictionary of bonded atoms

    Returns
    -------
    coordinate_list : list
        a refreshed list of lists of coordinates for every atom
    """

    bonded_dict = bonded

    for j in range(atom_number):
        if j + 1 in bonded_dict:
            for k in range(3):
                random_number = uniform(-1, 1)
                coordinate_list[j][k] += random_number
                for ii in set(bonded_dict[j + 1]):
                    coordinate_list[ii - 1][k] += random_number
                    coordinate_list[ii - 1][k] = \
                    check_pbc(coordinate_list[ii - 1][k])
                coordinate_list[j][k] = check_pbc(coordinate_list[j][k])

        else:
            for k in range(3):
                coordinate_list[j][k] += uniform(-1, 1)
                coordinate_list[j][k] = check_pbc(coordinate_list[j][k])

    return coordinate_list

def check_pbc(coordinate):
    """Checks whether particular atom or set of atoms is not outside
    periodic box.

    Parameters
    ----------
    coordinate : float
        a single x or y or z coordinate of atom

    Returns
    -------
    coordinate : float
        a single x or y or z coordinate of atom after PBC correction
    """

    if coordinate > 30:
        coordinate -= 30

    elif coordinate < 0:
        coordinate += 30

    return coordinate

def main():
    frames = 1000
    coordinates = read_coordinates('example.xyz')
    save_coordinates('simulation.xyz', coordinates, frames)

if __name__ == "__main__":
    main()
