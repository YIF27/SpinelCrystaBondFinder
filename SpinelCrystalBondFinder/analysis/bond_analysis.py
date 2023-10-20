from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core import Structure
import numpy as np

class BondAnalyzer:
  def __init__(self, input_data):
    if isinstance(input_data, Structure):
      self.structure = input_data
    elif isinstance(input_data, str):  # Assuming a file path is given as a string
        self.structure = Poscar.from_file(input_data).structure
    else:
        raise ValueError("Input must be a pymatgen Structure or a path to a POSCAR file.")

  def compute_distances_and_displacements(self, ox, met):
    self.ox = ox
    self.met = met
    x_offsets = [-1, 0, 1]
    y_offsets = [-1, 0, 1]
    z_offsets = [-1, 0, 1]

    r = []
    d = []

    for x_offset in x_offsets:
        for y_offset in y_offsets:
            for z_offset in z_offsets:
                dx = self.ox[0] - self.met[0] + x_offset
                dy = self.ox[1] - self.met[1] + y_offset
                dz = self.ox[2] - self.met[2] + z_offset

                r.append(np.sqrt(dx ** 2 + dy ** 2 + dz ** 2))
                d.append([dx, dy, dz])

    return r, d

  def count_bond_number(self, list_in, min_degree, max_degree):
    self.list_in = list_in
    self.min_degree = min_degree
    self.max_degree = max_degree
    count = 0
    for i in range(len(self.list_in)):
        if self.min_degree < self.list_in[i][-1]:
            if self.list_in[i][-1] <= self.max_degree:
                count = count + 1
    return count



  def generate_bond_info(self):
    composition_values = list(self.structure.composition.as_dict().values())
    metal1_number = int(composition_values[0])
    metal2_number = int(composition_values[1])
    oxygen_number = int(composition_values[2])
    metal1_coords = self.structure.frac_coords[0: metal1_number]
    metal2_coords = self.structure.frac_coords[metal1_number: metal2_number + metal1_number]
    oxygen_coords = self.structure.frac_coords[metal2_number + metal1_number:]

    R_metal1_oxygen = np.zeros((oxygen_number, metal1_number))
    PBD_O_metal1_type = np.zeros((oxygen_number, metal1_number))
    c_O_metal1 = np.zeros((oxygen_number * metal1_number, 3))
    t = 0
    for i in range(oxygen_number):
        for j in range(metal1_number):
            r, d = self.compute_distances_and_displacements(oxygen_coords[i], metal1_coords[j])
            minimum_r = min(r)
            arg_minimum_r = int(np.argmin(r))
            R_metal1_oxygen[i, j] = minimum_r
            PBD_O_metal1_type[i, j] = arg_minimum_r
            c_O_metal1[t] = d[arg_minimum_r]
            t = t + 1

    R_metal2_oxygen = np.zeros((oxygen_number, metal2_number))
    PBD_O_metal2_type = np.zeros((oxygen_number, metal2_number))
    c_O_metal2 = np.zeros((oxygen_number * metal2_number, 3))
    t = 0
    for i in range(oxygen_number):
        for j in range(metal2_number):
            r, d = self.compute_distances_and_displacements(oxygen_coords[i], metal2_coords[j])
            minimum_r = min(r)
            arg_minimum_r = int(np.argmin(r))
            R_metal2_oxygen[i, j] = minimum_r
            PBD_O_metal2_type[i, j] = arg_minimum_r
            c_O_metal2[t] = d[arg_minimum_r]
            t = t + 1

    metal1_oxygen_metal2_1 = []
    for i in range(oxygen_number):
        for j in range(metal1_number):
            for k in range(metal2_number):
                if R_metal1_oxygen[i][j] < 0.255 and R_metal2_oxygen[i][k] < 0.255:
                    metal1_oxygen_metal2_1_tmp = [i, j, k]
                    metal1_oxygen_metal2_1_tmp.append(R_metal1_oxygen[i][j])
                    d1 = R_metal1_oxygen[i][j]
                    metal1_oxygen_metal2_1_tmp.append(R_metal2_oxygen[i][k])
                    d2 = R_metal2_oxygen[i][k]
                    c1 = c_O_metal1[i * metal1_number + j][:]
                    c2 = c_O_metal2[i * metal2_number + k][:]
                    p = np.dot(c1, c2) / (d1 * d2)
                    metal1_oxygen_metal2_1_tmp.append(np.degrees(np.arccos(p)))
                    metal1_oxygen_metal2_1.append(metal1_oxygen_metal2_1_tmp)

    metal1_oxygen_metal2_2 = []
    for i in range(oxygen_number):
        for j in range(metal1_number):
            for k in range(metal2_number):
                if R_metal1_oxygen[i][j] > 0.402:
                    if R_metal1_oxygen[i][j] < 0.46:
                        if R_metal2_oxygen[i][k] < 0.255:
                            metal1_oxygen_metal2_2_tmp = [i, j, k]
                            metal1_oxygen_metal2_2_tmp.append(R_metal1_oxygen[i][j])
                            d1 = R_metal1_oxygen[i][j]
                            metal1_oxygen_metal2_2_tmp.append(R_metal2_oxygen[i][k])
                            d2 = R_metal2_oxygen[i][k]
                            c1 = c_O_metal1[i * metal1_number + j][:]
                            c2 = c_O_metal2[i * metal2_number + k][:]
                            p = np.dot(c1, c2) / (d1 * d2)
                            if p < 0.2588:
                              metal1_oxygen_metal2_2_tmp.append(np.degrees(np.arccos(p)))
                              metal1_oxygen_metal2_2.append(metal1_oxygen_metal2_2_tmp)

    metal1_oxygen_metal2_3 = []
    for i in range(oxygen_number):
        for j in range(metal1_number):
            for k in range(metal2_number):
                if R_metal1_oxygen[i][j] < 0.255:
                    if R_metal2_oxygen[i][k] > 0.40:
                        if R_metal2_oxygen[i][k] < 0.46:
                            metal1_oxygen_metal2_3_tmp = []
                            metal1_oxygen_metal2_3_tmp.append(i)
                            metal1_oxygen_metal2_3_tmp.append(j)
                            metal1_oxygen_metal2_3_tmp.append(k)
                            metal1_oxygen_metal2_3_tmp.append(R_metal1_oxygen[i][j])
                            d1 = R_metal1_oxygen[i][j]
                            metal1_oxygen_metal2_3_tmp.append(R_metal2_oxygen[i][k])
                            d2 = R_metal2_oxygen[i][k]
                            c1 = c_O_metal1[i * metal1_number + j][:]
                            c2 = c_O_metal2[i * metal2_number + k][:]
                            p = np.dot(c1, c2) / (d1 * d2)
                            if p < 0.2588:
                              metal1_oxygen_metal2_3_tmp.append(np.degrees(np.arccos(p)))
                              metal1_oxygen_metal2_3.append(metal1_oxygen_metal2_3_tmp)

    metal2_oxygen_metal2_1 = []
    for i in range(oxygen_number):
        for j in range(metal2_number - 1):
            for k in range(j + 1, metal2_number):
                if R_metal2_oxygen[i][j] < 0.255:
                    if R_metal2_oxygen[i][k] < 0.255:
                        metal2_oxygen_metal2_1_tmp = []
                        metal2_oxygen_metal2_1_tmp.append(i)
                        metal2_oxygen_metal2_1_tmp.append(j)
                        metal2_oxygen_metal2_1_tmp.append(k)
                        metal2_oxygen_metal2_1_tmp.append(R_metal2_oxygen[i][j])
                        d1 = R_metal2_oxygen[i][j]
                        metal2_oxygen_metal2_1_tmp.append(R_metal2_oxygen[i][k])
                        d2 = R_metal2_oxygen[i][k]
                        c1 = c_O_metal2[i * metal2_number + j][:]
                        c2 = c_O_metal2[i * metal2_number + k][:]
                        p = np.dot(c1, c2) / (d1 * d2)
                        metal2_oxygen_metal2_1_tmp.append(np.degrees(np.arccos(p)))
                        metal2_oxygen_metal2_1.append(metal2_oxygen_metal2_1_tmp)

    metal2_oxygen_metal2_2 = []
    R_metal2_oxygen_1 = []
    R_metal2_oxygen_2 = []
    for i in range(oxygen_number):
        for j in range(metal2_number - 1):
            for k in range(j + 1, metal2_number):
                if (R_metal2_oxygen[i][j] < 0.255 and 0.40 < R_metal2_oxygen[i][k] < 0.46) or (
                        R_metal2_oxygen[i][k] < 0.255 and 0.40 < R_metal2_oxygen[i][j] < 0.46):
                    metal2_oxygen_metal2_2_tmp = [i, j, k]
                    metal2_oxygen_metal2_2_tmp.append(R_metal2_oxygen[i][j])
                    d1 = R_metal2_oxygen[i][j]
                    metal2_oxygen_metal2_2_tmp.append(R_metal2_oxygen[i][k])
                    d2 = R_metal2_oxygen[i][k]
                    c1 = c_O_metal2[i * metal2_number + j][:]
                    c2 = c_O_metal2[i * metal2_number + k][:]
                    p = np.dot(c1, c2) / (d1 * d2)
                    if p < 0.2588:
                        metal2_oxygen_metal2_2_tmp.append(np.degrees(np.arccos(p)))
                        metal2_oxygen_metal2_2.append(metal2_oxygen_metal2_2_tmp)
                        if d1 > d2:
                            R_metal2_oxygen_1.append(d2)
                            R_metal2_oxygen_2.append(d1)
                        if d2 > d1:
                            R_metal2_oxygen_1.append(d1)
                            R_metal2_oxygen_2.append(d2)

    metal1_oxygen_metal1_1 = []
    R_metal1_oxygen_1 = []
    R_metal1_oxygen_2 = []
    for i in range(oxygen_number):
        for j in range(metal1_number - 1):
            for k in range(j + 1, metal1_number):
                if (R_metal1_oxygen[i][j] < 0.255 and 0.40 < R_metal1_oxygen[i][k] < 0.46) or (
                        R_metal1_oxygen[i][k] < 0.255 and 0.40 < R_metal1_oxygen[i][j] < 0.46):
                    metal1_oxygen_metal1_1_tmp = [i,j,k]
                    metal1_oxygen_metal1_1_tmp.append(R_metal1_oxygen[i][j])
                    d1 = R_metal1_oxygen[i][j]
                    metal1_oxygen_metal1_1_tmp.append(R_metal1_oxygen[i][k])
                    d2 = R_metal1_oxygen[i][k]
                    c1 = c_O_metal1[i * metal1_number + j][:]
                    c2 = c_O_metal1[i * metal1_number + k][:]
                    p = np.dot(c1, c2) / (d1 * d2)
                    if p < 0.2588:
                        metal1_oxygen_metal1_1_tmp.append(np.degrees(np.arccos(p)))
                        metal1_oxygen_metal1_1.append(metal1_oxygen_metal1_1_tmp)
                        if d1 > d2:
                            R_metal1_oxygen_1.append(d2)
                            R_metal1_oxygen_2.append(d1)
                        if d2 > d1:
                            R_metal1_oxygen_1.append(d1)
                            R_metal1_oxygen_2.append(d2)

    metal1_oxygen_metal1_2 = []
    for i in range(oxygen_number):
        for j in range(metal1_number - 1):
            for k in range(j + 1, metal1_number):
                if R_metal1_oxygen[i][j] < 0.255:
                    if R_metal1_oxygen[i][k] < 0.255:
                        d1 = R_metal1_oxygen[i][j]
                        d2 = R_metal1_oxygen[i][k]
                        c1 = c_O_metal1[i * metal1_number + j][:]
                        c2 = c_O_metal1[i * metal1_number + k][:]
                        p = np.dot(c1, c2) / (d1 * d2)
                        angle = np.arccos(p) * 180 / np.pi
                        metal1_oxygen_metal1_2_tmp = [i, j, k, d1, d2, angle]
                        metal1_oxygen_metal1_2.append(metal1_oxygen_metal1_2_tmp)

    metal1_oxygen_metal2_1 = np.array(metal1_oxygen_metal2_1)
    metal1_oxygen_metal2_1_new = metal1_oxygen_metal2_1[metal1_oxygen_metal2_1[:, 5].argsort()]

    metal1_oxygen_metal2_2 = np.array(metal1_oxygen_metal2_2)
    metal1_oxygen_metal2_2_new = metal1_oxygen_metal2_2[metal1_oxygen_metal2_2[:, 5].argsort()]

    metal1_oxygen_metal2_3 = np.array(metal1_oxygen_metal2_3)
    metal1_oxygen_metal2_3_new = metal1_oxygen_metal2_3[metal1_oxygen_metal2_3[:, 5].argsort()]

    metal1_oxygen_metal1_1 = np.array(metal1_oxygen_metal1_1)
    metal1_oxygen_metal1_1_new = metal1_oxygen_metal1_1[metal1_oxygen_metal1_1[:, 5].argsort()]

    metal1_oxygen_metal1_2 = np.array(metal1_oxygen_metal1_2)
    metal1_oxygen_metal1_2_new = metal1_oxygen_metal1_2[metal1_oxygen_metal1_2[:, 5].argsort()]

    metal2_oxygen_metal2_1 = np.array(metal2_oxygen_metal2_1)
    metal2_oxygen_metal2_1_new = metal2_oxygen_metal2_1[metal2_oxygen_metal2_1[:, 5].argsort()]

    metal2_oxygen_metal2_2 = np.array(metal2_oxygen_metal2_2)
    metal2_oxygen_metal2_2_new = metal2_oxygen_metal2_2[metal2_oxygen_metal2_2[:, 5].argsort()]

    metal2_oxygen_metal2_79 = self.count_bond_number(metal2_oxygen_metal2_2_new, 76, 80)
    metal2_oxygen_metal2_93 = self.count_bond_number(metal2_oxygen_metal2_1_new, 0, 94)
    metal2_oxygen_metal2_123 = self.count_bond_number(metal2_oxygen_metal2_1_new, 120, 126)
    metal2_oxygen_metal2_126 = self.count_bond_number(metal2_oxygen_metal2_2_new, 124, 127)
    metal2_oxygen_metal2_158 = self.count_bond_number(metal2_oxygen_metal2_2_new, 157, 160)
    metal2_oxygen_metal2_179 = self.count_bond_number(metal2_oxygen_metal2_2_new, 177, 180)
    metal1_oxygen_metal2_79 = self.count_bond_number(metal1_oxygen_metal2_3_new, 0, 80)
    metal1_oxygen_metal2_93 = self.count_bond_number(metal1_oxygen_metal2_1_new, 92, 94)
    metal1_oxygen_metal2_123 = self.count_bond_number(metal1_oxygen_metal2_1_new, 120, 126)
    metal1_oxygen_metal2_126 = self.count_bond_number(metal1_oxygen_metal2_3_new, 124, 127)
    metal1_oxygen_metal2_158 = self.count_bond_number(metal1_oxygen_metal2_3_new, 157, 160)
    metal1_oxygen_metal2_179 = self.count_bond_number(metal1_oxygen_metal2_3_new, 177, 180)
    metal2_oxygen_metal1_79 = self.count_bond_number(metal1_oxygen_metal2_2_new, 77, 81)
    metal2_oxygen_metal1_126 = self.count_bond_number(metal1_oxygen_metal2_2_new, 124, 127)
    metal2_oxygen_metal1_158 = self.count_bond_number(metal1_oxygen_metal2_2_new, 157, 160)
    metal2_oxygen_metal1_179 = self.count_bond_number(metal1_oxygen_metal2_2_new, 177, 180)
    metal1_oxygen_metal1_79 = self.count_bond_number(metal1_oxygen_metal1_1_new, 77, 81)
    metal1_oxygen_metal1_93 = self.count_bond_number(metal1_oxygen_metal1_2_new, 90, 95)
    metal1_oxygen_metal1_123 = self.count_bond_number(metal1_oxygen_metal1_2_new, 120, 125)
    metal1_oxygen_metal1_126 = self.count_bond_number(metal1_oxygen_metal1_1_new, 124, 127)
    metal1_oxygen_metal1_158 = self.count_bond_number(metal1_oxygen_metal1_1_new, 157, 160)
    metal1_oxygen_metal1_179 = self.count_bond_number(metal1_oxygen_metal1_1_new, 177, 180)
    bond = [metal1_oxygen_metal2_123, metal2_oxygen_metal2_123, metal1_oxygen_metal1_123, metal2_oxygen_metal1_158,
            metal2_oxygen_metal2_158, metal1_oxygen_metal2_158, metal1_oxygen_metal1_158, metal1_oxygen_metal2_179,
            metal2_oxygen_metal2_179, metal2_oxygen_metal1_179, metal1_oxygen_metal1_179,
            metal1_oxygen_metal2_93, metal2_oxygen_metal2_93, metal1_oxygen_metal1_93, metal2_oxygen_metal1_126,
            metal1_oxygen_metal2_126, metal2_oxygen_metal2_126, metal1_oxygen_metal1_126, metal2_oxygen_metal2_79,
            metal2_oxygen_metal1_79, metal1_oxygen_metal2_79, metal1_oxygen_metal1_79]

    return bond


