from __future__ import print_function
import sys
import os
sys.executable

import numpy
numpy.set_printoptions(threshold=sys.maxsize)
import pandas as pd
import csv
from math import factorial

import MDAnalysis as mda
import argparse


def calculate_distances(trj_prot, aas_range):
    all_distances = []

    for ts in trj_prot.trajectory:
        com_coordinates = []
        dist_list = []

        for i in aas_range:
            atoms_i = trj_prot.select_atoms(f'resid {i} and name N DN CA DCA C DC O DO LPOA LPOB')
            com_coordinates.append(atoms_i.center_of_mass())

        for i in range(len(aas_range)):
            for j in range(i + 1, len(aas_range)):
                dist = numpy.linalg.norm(com_coordinates[j] - com_coordinates[i])
                dist_list.append(dist)

        all_distances.append(dist_list)

    return all_distances

def calculate_angles(trj_prot, aas_range):
    all_angles = []

    for ts in trj_prot.trajectory:
        com_coordinates = []
        angle_list = []

        for i in aas_range:
            atoms_i = trj_prot.select_atoms(f'resid {i} and name N DN CA DCA C DC O DO LPOA LPOB')
            com_coordinates.append(atoms_i.center_of_mass())

        for i in range(len(aas_range)):
            for j in range(i + 1, len(aas_range)):
                for k in range(j + 1, len(aas_range)):
                    vec1 = com_coordinates[i] - com_coordinates[j]
                    vec2 = com_coordinates[k] - com_coordinates[j]

                    angle = numpy.arccos(numpy.dot(vec1, vec2) / (numpy.linalg.norm(vec1) * numpy.linalg.norm(vec2)))
                    angle_list.append(numpy.degrees(angle))
        all_angles.append(angle_list)

    return all_angles

def calculate_dist_hp(trj_prot, s1, s2):
    all_distances = []

    for ts in trj_prot.trajectory:
        com_coordinates_i = []
        com_coordinates_j = []
        dist_list = []

        for i in s1:
            atoms_i = trj_prot.select_atoms(f"resid {i} and name P DP O1P DO1P O2P DO2P O5' DO5' C5' DC5' C4' DC4' O4' DO4' LPRA LPRB LPX C1' DC1' C2' DC2' C3' DC3' O3' DO3'")
            com_coordinates_i.append(atoms_i.center_of_mass())
        for j in s2:
            atoms_j = trj_prot.select_atoms(f"resid {j} and name P DP O1P DO1P O2P DO2P O5' DO5' C5' DC5' C4' DC4' O4' DO4' LPRA LPRB LPX C1' DC1' C2' DC2' C3' DC3' O3' DO3'")
            com_coordinates_j.append(atoms_j.center_of_mass())

        for i in range(len(s1)):
            for j in range(len(s2)):
                dist = numpy.linalg.norm(com_coordinates_i[i] - com_coordinates_j[j])
                dist_list.append(dist)
        all_distances.append(dist_list) 
            
    return all_distances


def calculate_dist_hb(trj_prot):
    all_distances_hb = []

    for ts in trj_prot.trajectory:
        dist_list = []

        # pair resid 2 H3 to resid 6 O2
        atoms_i = trj_prot.select_atoms("resid 2 and name H3")
        atoms_j = trj_prot.select_atoms("resid 6 and name O2")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # Pair 2: resid 2 C2 to resid 6 H3
        atoms_i = trj_prot.select_atoms("resid 2 and name C2")
        atoms_j = trj_prot.select_atoms("resid 6 and name H3")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # Pair 3: resid 2 H3 to resid 6 O4
        atoms_i = trj_prot.select_atoms("resid 2 and name H3")
        atoms_j = trj_prot.select_atoms("resid 6 and name O4")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # Pair 4: resid 2 O2 to resid 6 H3
        atoms_i = trj_prot.select_atoms("resid 2 and name O2")
        atoms_j = trj_prot.select_atoms("resid 6 and name H3")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # Pair 5: resid 5 H3 to resid 8 N1
        atoms_i = trj_prot.select_atoms("resid 5 and name H3")
        atoms_j = trj_prot.select_atoms("resid 8 and name N1")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # Append the distances for the current timestep
        all_distances_hb.append(dist_list)

    return all_distances_hb

def calculate_dist_stack(trj_prot):
    all_distances_stack = []
    base = "(name N9 DN9 C4 DC4 N3 DN3 LP3 C2 DC2 H2 N1 DN1 LP1 C6 DC6 N6 DN6 H61 H62 C5 DC5 N7 DN7 LP7 C8 DC8 H8 H1 H6 N2 DN2 H21 H22 C6 DC6 O6 DO6 LP6A LP6B O2 DO2 LP2A LP2B N3 DN3 H3 LP3 C4 DC4 O4 DO4 LP4A LP4B N4 DN4 H41 H42 C5M DC5M H51 H52 H53 H5)"
    for ts in trj_prot.trajectory:
        dist_list = []

        # pair resid 1 to resid 2
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 2 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 1 to resid 3
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 3 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 1 to resid 4
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 4 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 1 to resid 5
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 5 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 1 to resid 6
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 6 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)


        # pair resid 1 to resid 7
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 7 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 1 to resid 8
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 8 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)


        # pair resid 1 to resid 9
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 9 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 1 to resid 10
        atoms_i = trj_prot.select_atoms(f"resid 1 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 10 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 2 to resid 3
        atoms_i = trj_prot.select_atoms(f"resid 2 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 3 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 2 to resid 5
        atoms_i = trj_prot.select_atoms(f"resid 2 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 5 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 2 to resid 6
        atoms_i = trj_prot.select_atoms(f"resid 2 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 6 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 2 to resid 7
        atoms_i = trj_prot.select_atoms(f"resid 2 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 7 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 2 to resid 8
        atoms_i = trj_prot.select_atoms(f"resid 2 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 8 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 2 to resid 10
        atoms_i = trj_prot.select_atoms(f"resid 2 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 10 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 3 to resid 5
        atoms_i = trj_prot.select_atoms(f"resid 3 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 5 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 3 to resid 6
        atoms_i = trj_prot.select_atoms(f"resid 3 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 6 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 3 to resid 7
        atoms_i = trj_prot.select_atoms(f"resid 3 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 7 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 3 to resid 8
        atoms_i = trj_prot.select_atoms(f"resid 3 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 8 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 3 to resid 10
        atoms_i = trj_prot.select_atoms(f"resid 3 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 10 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 4 to resid 8
        atoms_i = trj_prot.select_atoms(f"resid 4 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 8 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 5 to resid 6
        atoms_i = trj_prot.select_atoms(f"resid 5 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 6 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 5 to resid 7
        atoms_i = trj_prot.select_atoms(f"resid 5 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 7 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 5 to resid 8
        atoms_i = trj_prot.select_atoms(f"resid 5 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 8 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 5 to resid 10
        atoms_i = trj_prot.select_atoms(f"resid 5 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 10 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 6 to resid 7
        atoms_i = trj_prot.select_atoms(f"resid 6 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 7 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 6 to resid 8
        atoms_i = trj_prot.select_atoms(f"resid 6 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 8 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 6 to resid 10
        atoms_i = trj_prot.select_atoms(f"resid 6 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 10 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 7 to resid 8
        atoms_i = trj_prot.select_atoms(f"resid 7 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 8 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 7 to resid 10
        atoms_i = trj_prot.select_atoms(f"resid 7 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 10 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)

        # pair resid 8 to resid 10
        atoms_i = trj_prot.select_atoms(f"resid 8 and {base}")
        atoms_j = trj_prot.select_atoms(f"resid 10 and {base}")
        if len(atoms_i) > 0 and len(atoms_j) > 0:
            dist = numpy.linalg.norm(atoms_i.center_of_mass() - atoms_j.center_of_mass())
            dist_list.append(dist)


        # Append the distances for the current timestep
        all_distances_stack.append(dist_list)

    return all_distances_stack
