from math import acos, pi, ceil
from scipy.spatial.distance import cdist
import numpy as np
import matplotlib.pyplot as plt
import re
from math import sin,cos,sqrt
from scipy.ndimage.filters import gaussian_filter1d

def angle(a,b):
    """ calculate the angle between vector a and b """
    return acos(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))

def smear(data, sigma):
    """
    Apply Gaussian smearing to spectrum y value.

    Args:
        sigma: Std dev for Gaussian smear function
    """
    diff = [data[i + 1,0] - data[i,0] for i in range(len(data) - 1)]
    avg_x_per_step = np.sum(diff) / len(diff)
    data[:,1] = gaussian_filter1d(data[:,1], sigma / avg_x_per_step)
    return data

class Element:
    def __init__(self, input_value):
        self.input = input_value

        # list with atomic number z, short name, full name, valence, 
                # valence electrons, covalent radius, good bonds, Maximum CN:
        self.elements_list = [
            (1, 'H', 'Hydrogen',    1.0, 1, 0.31),
            (2, 'He', 'Helium',     0.5, 2, 0.28),
            (3, 'Li', 'Lithium',    1.0, 1, 1.28),
            (4, 'Be', 'Beryllium',  2.0, 2, 0.96),
            (5, 'B', 'Boron',       3.0, 3, 0.84),
            (6, 'C', 'Carbon',      4.0, 4, 0.70),
            (7, 'N', 'Nitrogen',    3.0, 5, 0.71),
            (8, 'O', 'Oxygen',      2.0, 6, 0.66),
            (9, 'F', 'Fluorine',    1.0, 7, 0.57),
            (10, 'Ne', 'Neon',      0.5, 8, 0.58),
            (11, 'Na', 'Sodium',    1.0, 1, 1.66),
            (12, 'Mg', 'Magnesium', 2.0, 2, 1.41),
            (13, 'Al', 'Aluminium', 3.0, 3, 1.21),
            (14, 'Si', 'Silicon',   4.0, 4, 1.11),
            (15, 'P', 'Phosphorus', 3.0, 5, 1.07),
            (16, 'S', 'Sulfur',     2.0, 6, 1.05),
            (17, 'Cl', 'Chlorine',  1.0, 7, 1.02),
            (18, 'Ar', 'Argon',     0.5, 8, 1.06),
            (19, 'K', 'Potassium',  1.0, 1, 2.03),
            (20, 'Ca', 'Calcium',   2.0, 2, 1.76),
            (21, 'Sc', 'Scandium',  3.0, 3, 1.70),
            (22, 'Ti', 'Titanium',  4.0, 4, 1.60),
            (23, 'V', 'Vanadium',   4.0, 5, 1.53),
            (24, 'Cr', 'Chromium',  3.0, 6, 1.39),
            (25, 'Mn', 'Manganese', 4.0, 5, 1.39),
            (26, 'Fe', 'Iron',      3.0, 3, 1.32),
            (27, 'Co', 'Cobalt',    3.0, 3, 1.26),
            (28, 'Ni', 'Nickel',    2.0, 3, 1.24),
            (29, 'Cu', 'Copper',    2.0, 2, 1.32),
            (30, 'Zn', 'Zinc',      2.0, 2, 1.22),
            (31, 'Ga', 'Gallium',   3.0, 3, 1.22),
            (32, 'Ge', 'Germanium', 4.0, 4, 1.20),
            (33, 'As', 'Arsenic',   3.0, 5, 1.19),
            (34, 'Se', 'Selenium',  2.0, 6, 1.20),
            (35, 'Br', 'Bromine',   1.0, 7, 1.20),
            (36, 'Kr', 'Krypton',   0.5, 8, 1.16),
            (37, 'Rb', 'Rubidium',  1.0, 1, 2.20),
            (38, 'Sr', 'Strontium', 2.0, 2, 1.95),
            (39, 'Y', 'Yttrium',    3.0, 3, 1.90),
            (40, 'Zr', 'Zirconium', 4.0, 4, 1.75),
            (41, 'Nb', 'Niobium',   5.0, 5, 1.64),
            (42, 'Mo', 'Molybdenum',4.0, 6, 1.54),
            (43, 'Tc', 'Technetium',4.0, 5, 1.47),
            (44, 'Ru', 'Ruthenium', 4.0, 3, 1.46),
            (45, 'Rh', 'Rhodium',   4.0, 3, 1.42),
            (46, 'Pd', 'Palladium', 4.0, 3, 1.39),
            (47, 'Ag', 'Silver',    1.0, 2, 1.45),
            (48, 'Cd', 'Cadmium',   2.0, 2, 1.44),
            (49, 'In', 'Indium',    3.0, 3, 1.42),
            (50, 'Sn', 'Tin',       4.0, 4, 1.39),
            (51, 'Sb', 'Antimony',  3.0, 5, 1.39),
            (52, 'Te', 'Tellurium', 2.0, 6, 1.38),
            (53, 'I', 'Iodine',     1.0, 7, 1.39),
            (54, 'Xe', 'Xenon',     0.5, 8, 1.40),
            (55, 'Cs', 'Caesium',   1.0, 1, 2.44),
            (56, 'Ba', 'Barium',    2.0, 2, 2.15),
            (57, 'La', 'Lanthanum', 3.0, 3, 2.07),
            (58, 'Ce', 'Cerium',    4.0, 3, 2.04),
            (59,'Pr','Praseodymium',3.0, 3, 2.03),
            (60, 'Nd', 'Neodymium', 3.0, 3, 2.01),
            (61, 'Pm', 'Promethium',3.0, 3, 1.99),
            (62, 'Sm', 'Samarium',  3.0, 3, 1.98),
            (63, 'Eu', 'Europium',  3.0, 3, 1.98),
            (64, 'Gd', 'Gadolinium',3.0, 3, 1.96),
            (65, 'Tb', 'Terbium',   3.0, 3, 1.94),
            (66, 'Dy', 'Dysprosium',3.0, 3, 1.92),
            (67, 'Ho', 'Holmium',   3.0, 3, 1.92),
            (68, 'Er', 'Erbium',    3.0, 3, 1.89),
            (69, 'Tm', 'Thulium',   3.0, 3, 1.90),
            (70, 'Yb', 'Ytterbium', 3.0, 3, 1.87),
            (71, 'Lu', 'Lutetium',  3.0, 3, 1.87),
            (72, 'Hf', 'Hafnium',   4.0, 3, 1.75),
            (73, 'Ta', 'Tantalum',  5.0, 3, 1.70),
            (74, 'W', 'Tungsten',   4.0, 3, 1.62),
            (75, 'Re', 'Rhenium',   4.0, 3, 1.51),
            (76, 'Os', 'Osmium',    4.0, 3, 1.44),
            (77, 'Ir', 'Iridium',   4.0, 3, 1.41),
            (78, 'Pt', 'Platinum',  4.0, 3, 1.36),
            (79, 'Au', 'Gold',      1.0, 3, 1.36),
            (80, 'Hg', 'Mercury',   2.0, 3, 1.32),
            (81, 'Tl', 'Thallium',  3.0, 3, 1.45),
            (82, 'Pb', 'Lead',      4.0, 4, 1.46),
            (83, 'Bi', 'Bismuth',   3.0, 5, 1.48),
            (84, 'Po', 'Polonium',  2.0, 6, 1.40),
            (85, 'At', 'Astatine',  1.0, 7, 1.50),
            (86, 'Rn', 'Radon',     0.5, 8, 1.50),
            (87, 'Fr', 'Francium',  1.0, 1, 2.60),
            (88, 'Ra', 'Radium',    2.0, 2, 2.21),
            (89, 'Ac', 'Actinium',  3.0, 3, 2.15),
            (90, 'Th', 'Thorium',   4.0, 3, 2.06),
            (91,'Pa','Protactinium',4.0, 3, 2.00),
            (92, 'U', 'Uranium',    4.0, 3, 1.96),
            (93, 'Np', 'Neptunium', 4.0, 3, 1.90),
            (94, 'Pu', 'Plutonium', 4.0, 3, 1.87),
            (95, 'Am', 'Americium', 4.0, 3, 1.80),
            (96, 'Cm', 'Curium',    4.0, 3, 1.69),
            (97, 'Bk', 'Berkelium', 4.0, 3, None),
            (98,'Cf','Californium', 4.0, 3, None),
            (99,'Es','Einsteinium', 4.0, 3, None),
            (100, 'Fm', 'Fermium',  4.0, 3, None),
            (101,'Md','Mendelevium',4.0, 3, None),
            (102, 'No', 'Nobelium', 4.0, 3, None),
            (103, 'Lr','Lawrencium',4.0, 3, None),
            (104,'Rf','Rutherfordium',4.0,3,None),
            (105, 'Db', 'Dubnium',  2.0, 3, None),
        ]
        
        #scatter factor
        self.sf=[
            [  0.493,  0.323,  0.140,  0.041, 10.511, 26.126,  3.142, 57.800,  0.003],
            [  0.873,  0.631,  0.311,  0.178,  9.104,  3.357, 22.928,  0.982,  0.006],
            [  1.128,  0.751,  0.618,  0.465,  3.955,  1.052, 85.391,168.261,  0.038],
            [  1.592,  1.128,  0.539,  0.703, 43.643,  1.862,103.483,  0.542,  0.038],
            [  2.055,  1.333,  1.098,  0.707, 23.219,  1.021, 60.350,  0.140, -0.193],
            [  2.310,  1.020,  1.589,  0.865, 20.844, 10.208,  0.569, 51.651,  0.216],
            [ 12.213,  3.132,  2.013,  1.166,  0.006,  9.893, 28.997,  0.583,-11.529],
            [  3.049,  2.287,  1.546,  0.867, 13.277,  5.701,  0.324, 32.909,  0.251],
            [  3.539,  2.641,  1.517,  1.024, 10.283,  4.294,  0.262, 26.148,  0.278],
            [  3.955,  3.112,  1.455,  1.125,  8.404,  3.426,  0.231, 21.718,  0.352],
            [  4.763,  3.174,  1.267,  1.113,  3.285,  8.842,  0.314,129.424,  0.676],
            [  5.420,  2.174,  1.227,  2.307,  2.828, 79.261,  0.381,  7.194,  0.858],
            [  6.420,  1.900,  1.594,  1.965,  3.039,  0.743, 31.547, 85.089,  1.115],
            [  6.292,  3.035,  1.989,  1.541,  2.439, 32.334,  0.678, 81.694,  1.141],
            [  6.435,  4.179,  1.780,  1.491,  1.907, 27.157,  0.526, 68.164,  1.115],
            [  6.905,  5.203,  1.438,  1.586,  1.468, 22.215,  0.254, 56.172,  0.867],
            [ 11.460,  7.196,  6.256,  1.645,  0.010,  1.166, 18.519, 47.778, -9.557],
            [  7.484,  6.772,  0.654,  1.644,  0.907, 14.841, 43.898, 33.393,  1.444],
            [  8.219,  7.440,  1.052,  0.866, 12.795,  0.775,213.187, 41.684,  1.423],
            [  8.627,  7.387,  1.590,  1.021, 10.442,  0.660, 85.748,178.437,  1.375],
            [  9.189,  7.368,  1.641,  1.468,  9.021,  0.573,136.108, 51.353,  1.333],
            [  9.759,  7.356,  1.699,  1.902,  7.851,  0.500, 35.634,116.105,  1.281],
            [ 10.297,  7.351,  2.070,  2.057,  6.866,  0.438, 26.894,102.478,  1.220],
            [ 10.641,  7.354,  3.324,  1.492,  6.104,  0.392, 20.263, 98.740,  1.183],
            [ 11.282,  7.357,  3.019,  2.244,  5.341,  0.343, 17.867, 83.754,  1.090],
            [ 11.769,  7.357,  3.522,  2.305,  4.761,  0.307, 15.354, 76.881,  1.037],
            [ 12.284,  7.341,  4.003,  2.349,  4.279,  0.278, 13.536, 71.169,  1.012],
            [ 12.838,  7.292,  4.444,  2.380,  3.878,  0.257, 12.176, 66.342,  1.034],
            [ 13.338,  7.168,  5.616,  1.673,  3.583,  0.247, 11.397, 64.831,  1.191],
            [ 14.074,  7.032,  5.165,  2.410,  3.266,  0.233, 10.316, 58.710,  1.304],
            [ 15.235,  6.701,  4.359,  2.962,  3.067,  0.241, 10.781, 61.414,  1.719],
            [ 16.082,  6.375,  3.707,  3.683,  2.851,  0.252, 11.447, 54.763,  2.131],
            [ 16.672,  6.070,  3.431,  4.278,  2.635,  0.265, 12.948, 47.797,  2.531],
            [ 17.001,  5.820,  3.973,  4.354,  2.410,  0.273, 15.237, 43.816,  2.841],
            [ 17.179,  5.236,  5.638,  3.985,  2.172, 16.580,  0.261, 41.433,  2.956],
            [ 17.355,  6.729,  5.549,  3.537,  1.938, 16.562,  0.226, 39.397,  2.825],
            [ 17.178,  9.644,  5.140,  1.529,  1.789, 17.315,  0.275,164.934,  3.487],
            [ 17.566,  9.818,  5.422,  2.669,  1.556, 14.099,  0.166,132.376,  2.506],
            [ 17.776, 10.295,  5.726,  3.266,  1.403, 12.801,  0.261,104.354,  1.912],
            [ 17.876, 10.948,  5.417,  3.657,  1.276, 11.916,  0.118, 87.663,  2.069],
            [ 17.614, 12.014,  4.042,  3.533,  1.189, 11.766,  0.205, 69.796,  3.756],
            [  3.703, 17.236, 12.888,  3.743,  0.277,  1.096, 11.004, 61.658,  4.387],
            [ 19.130, 11.095,  4.649,  2.713,  0.864,  8.145, 21.571, 86.847,  5.404],
            [ 19.267, 12.918,  4.863,  1.568,  0.809,  8.435, 24.800, 94.293,  5.379],
            [ 19.296, 14.350,  4.734,  1.289,  0.752,  8.218, 25.875, 98.606,  5.328],
            [ 19.332, 15.502,  5.295,  0.606,  0.699,  7.989, 25.205, 76.899,  5.266],
            [ 19.281, 16.688,  4.805,  1.046,  0.645,  7.473, 24.660, 99.816,  5.179],
            [ 19.221, 17.644,  4.461,  1.603,  0.595,  6.909, 24.701, 87.482,  5.069],
            [ 19.162, 18.560,  4.295,  2.040,  0.548,  6.378, 25.850, 92.803,  4.939],
            [ 19.189, 19.101,  4.458,  2.466,  5.830,  0.503, 26.891, 83.957,  4.782],
            [ 19.642, 19.045,  5.037,  2.683,  5.303,  0.461, 27.907, 75.283,  4.591],
            [ 19.964, 19.014,  6.145,  2.524,  4.817,  0.421, 28.528, 70.840,  4.352],
            [ 20.147, 18.995,  7.514,  2.273,  4.347,  0.381, 27.766, 66.878,  4.071],
            [ 20.293, 19.030,  8.977,  1.990,  3.928,  0.344, 26.466, 64.266,  3.712],
            [ 20.389, 19.106, 10.662,  1.495,  3.569,  0.311, 24.388,213.904,  3.335],
            [ 20.336, 19.297, 10.888,  2.696,  3.216,  0.276, 20.207,167.202,  2.773],
            [ 20.578, 19.599, 11.373,  3.287,  2.948,  0.244, 18.773,133.124,  2.147],
            [ 21.167, 19.770, 11.851,  3.330,  2.812,  0.227, 17.608,127.113,  1.863],
            [ 22.044, 19.670, 12.386,  2.824,  2.774,  0.222, 16.767,143.644,  2.058],
            [ 22.684, 19.685, 12.774,  2.851,  2.662,  0.211, 15.885,137.903,  1.985],
            [ 23.340, 19.610, 13.123,  2.875,  2.563,  0.202, 15.101,132.721,  2.029],
            [ 24.004, 19.426, 13.440,  2.896,  2.473,  0.196, 14.400,128.007,  2.210],
            [ 24.627, 19.089, 13.760,  2.293,  2.388,  0.194, 13.755,123.174,  2.575],
            [ 25.071, 19.080, 13.852,  3.545,  2.253,  0.182, 12.933,101.398,  2.420],
            [ 25.898, 18.219, 14.317,  2.954,  2.243,  0.196, 12.665,115.362,  3.583],
            [ 26.507, 17.638, 14.560,  2.966,  2.180,  0.202, 12.190,111.874,  4.297],
            [ 26.905, 17.294, 14.558,  3.638,  2.071,  0.198, 11.441, 92.657,  4.568],
            [ 27.656, 16.428, 14.978,  2.982,  2.074,  0.224, 11.360,105.703,  5.920],
            [ 28.182, 15.885, 15.154,  2.987,  2.029,  0.239, 10.998,102.961,  6.756],
            [ 28.664, 15.434, 15.309,  2.990,  1.989,  0.257, 10.665,100.417,  7.567],
            [ 28.948, 15.221, 15.100,  3.716,  1.902,  9.985,  0.261, 84.330,  7.976],
            [ 29.144, 15.173, 14.759,  4.300,  1.833,  9.600,  0.275, 72.029,  8.582],
            [ 29.202, 15.229, 14.514,  4.765,  1.773,  9.370,  0.296, 63.364,  9.244],
            [  0.000,  0.000,  0.000,  0.000,  1.722,  9.231,  0.323, 57.725,  9.858],
            [ 28.762, 15.719, 14.556,  5.442,  1.672,  9.092,  0.350, 52.086, 10.472],
            [ 28.189, 16.155, 14.931,  5.676,  1.629,  8.979,  0.383, 48.165, 11.000],
            [ 27.305, 16.730, 15.611,  5.834,  1.593,  8.866,  0.418, 45.001, 11.472],
            [ 27.006, 17.764, 15.713,  5.784,  1.513,  8.812,  0.425, 38.610, 11.688],
            [ 16.882, 18.591, 25.558,  5.860,  0.461,  8.622,  1.483, 36.396, 12.066],
            [ 20.681, 19.042, 21.657,  5.968,  0.545,  8.448,  1.573, 38.325, 12.609],
            [ 27.545, 19.158, 15.538,  5.526,  0.655,  8.708,  1.963, 45.815, 13.175],
            [ 31.062, 13.064, 18.442,  5.970,  0.690,  2.358,  8.618, 47.258, 13.412],
            [ 33.369, 12.951, 16.588,  6.469,  0.704,  2.924,  8.794, 48.009, 13.578],
            [  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
            [  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
            [  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
            [  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
            [  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000],
            [ 35.660, 23.103, 12.598,  4.087,  0.589,  3.652, 18.599,117.020, 13.527],
            [ 35.564, 23.422, 12.747,  4.807,  0.563,  3.462, 17.831, 99.172, 13.431],
            [ 35.885, 23.295, 14.189,  4.173,  0.548,  3.415, 16.924,105.251, 13.429],
            [  0.000,  0.000,  0.000,  0.000,  0.530,  3.335, 16.143,101.371, 13.393],
            [ 36.187, 23.596, 15.640,  4.186,  0.512,  3.254, 15.362, 97.491, 13.357],
            [ 36.526, 23.808, 16.771,  3.479,  0.499,  3.264, 14.946,105.980, 13.381]]

        self.z = None
        self.short_name = None
        self.long_name = None
        self.valence = None
        self.valence_electrons = None
        self.covalent_radius = None

        pos = None

        try:
            int(self.input)
            self.z = self.input

            for i, el in enumerate(self.elements_list):
                if el[0] == self.z:
                    pos = i
                    self.short_name = el[1]
                    self.long_name = el[2]
                    break
        except ValueError:
            self.short_name = self.input
            for i, el in enumerate(self.elements_list):
                if el[1] == self.short_name:
                    pos = i
                    self.z = el[0]
                    self.long_name = el[2]
                    break

            if not self.z:
                self.short_name = None
                self.long_name = self.input
                for i, el in enumerate(self.elements_list):
                    if el[2] == self.long_name:
                        pos = i
                        self.z = el[0]
                        self.short_name = el[1]
                        break
                if not self.z:
                    self.long_name = None

        if pos is not None:
            self.valence = self.elements_list[pos][3]
            self.valence_electrons = self.elements_list[pos][4]
            self.covalent_radius = self.elements_list[pos][5]
            self.scatter = self.sf[pos]

    def get_all(self, pos):
        els = []
        for el in self.elements_list:
            els.append(el[pos])
        return els

    def get_sf(self, pos):
        els = []
        for el in self.sf:
            els.append(el[pos])
        return els

    def all_z(self):
        return self.get_all(0)
    def all_short_names(self):
        return self.get_all(1)
    def all_long_names(self):
        return self.get_all(2)
    def all_valences(self):
        return self.get_all(3)
    def all_valence_electrons(self):
        return self.get_all(4)
    def all_covalent_radii(self):
        return self.get_all(5)
    def get_sf(self):
        return self.get_sf()

class crystal(object):
    """a class of crystal structure. 
    Attributes:
        cell_para: a,b,c, alpha, beta, gamma
        cell_matrix: 3*3 matrix
        rec_matrix: reciprocal of cell matrix
        atom_type:  elemental type (e.g. Na Cl)
        composition: chemical composition (e.g., [1,1])
        coordinate: atomic positions (e.g., [[0,0,0],[0.5,0.5,0.5]])
    """

    def __init__(self, fileformat='POSCAR', filename=None, \
                 lattice=None, atom_type=None, composition=None, coordinate=None):
        """Return a structure object with the proper structures info"""
        if fileformat == 'POSCAR':
           self.from_POSCAR(filename)
        elif fileformat == 'cif':
           self.from_cif(filename)
        else:
           self.from_dict(lattice, atom_type, composition, coordinate)
    
    def from_cif(self, filename):
        cif_struc = cif(filename)
        lattice = self.para2matrix(cif_struc.cell_para)
        composition = cif_struc.composition
        coordinate = cif_struc.coordinate
        atom_type = cif_struc.atom_type
        self.from_dict(lattice, atom_type, composition, coordinate)

    def from_POSCAR(self, filename):

        f = open(filename)

        tag = f.readline()
        lattice_constant = float(f.readline().split()[0])

        # Now the lattice vectors
        a = []
        for ii in range(3):
            s = f.readline().split()
            floatvect = float(s[0]), float(s[1]), float(s[2])
            a.append(floatvect)
        lattice = np.array(a) * lattice_constant

        # Number of atoms. 
        atom_type = f.readline().split()
        comp = f.readline().split()
        composition = []
        if len(atom_type)==len(comp):
           for num in comp:
               composition.append(int(num))
        else:
           print('Value Error POSCAR symbol and composition is inconsistent')
        ac_type = f.readline().split()
        # Check if atom coordinates are cartesian or direct
        cartesian = ac_type[0].lower() == "c" or ac_type[0].lower() == "k"
        tot_natoms = sum(composition)
        coordinate = np.empty((tot_natoms, 3))
        for atom in range(tot_natoms):
            ac = f.readline().split()
            coordinate[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
        # Done with all reading
        f.close()
        if cartesian:
            coordinate *= lattice_constant
        cell_para = []
        self.from_dict(lattice, atom_type, composition, coordinate)

    def from_dict(self, lattice, atom_type, composition, coordinate):
        self.cell_matrix = np.array(lattice) 
        self.atom_type = atom_type
        self.composition = np.array(composition)
        self.coordinate = np.array(coordinate)
        self.cell_para = self.matrix2para(self.cell_matrix)
        self.rec_matrix = self.rec_lat(self.cell_matrix)
        self.volume = np.linalg.det(self.cell_matrix)
        self.name = ''
        for ele, num in zip(self.atom_type, self.composition):
            self.name += ele
            if num > 1:
               self.name += str(num)
   
    #def show(self, L=2):
    #    """show crystal structure"""
    #    
    #    for i in range(-L, L+1):
    #        for j in range(-L, L+1):
    #            for k in range(-L, L+1):
    #                for m in self.coordinate:
    #                    
    #                    sphere(pos=vector(m[0], m[1], m[2]), radius=R)
    
    @staticmethod
    def rec_lat(matrix):
        """ calculate the reciprocal lattice """
        rec_lat = np.zeros([3,3])
        V = np.linalg.det(matrix)
        rec_lat[0] = np.cross(matrix[1], matrix[2])/V
        rec_lat[1] = np.cross(matrix[2], matrix[0])/V
        rec_lat[2] = np.cross(matrix[0], matrix[1])/V
        return  rec_lat #* 2 * pi

    @staticmethod
    def matrix2para(matrix):
        """ 3x3 representation -> 1x6 (a, b, c, alpha, beta, gamma)"""
        cell_para = np.zeros(6)
        cell_para[0] = np.linalg.norm(matrix[0])
        cell_para[1] = np.linalg.norm(matrix[1])
        cell_para[2] = np.linalg.norm(matrix[2])
    
        cell_para[5] = angle(matrix[0], matrix[1])
        cell_para[4] = angle(matrix[0], matrix[2])
        cell_para[3] = angle(matrix[1], matrix[2])

        return cell_para

    @staticmethod
    def para2matrix(cell_para):
        """ 1x6 (a, b, c, alpha, beta, gamma) -> 3x3 representation -> """
        matrix = np.zeros([3,3])
        matrix[0][0] = cell_para[0]
        matrix[1][0] = cell_para[1]*cos(cell_para[5])
        matrix[1][1] = cell_para[1]*sin(cell_para[5])
        matrix[2][0] = cell_para[2]*cos(cell_para[4])
        matrix[2][1] = cell_para[2]*cos(cell_para[3])*sin(cell_para[4])
        matrix[2][2] = sqrt(cell_para[2]**2 - matrix[2][0]**2 - matrix[2][1]**2)
        
        return matrix
    
class cif(object):
    """a class of cif reader
    Attributes:
        wavelength: default: 1.54181a, namely Cu-Ka
        max2theta: the range of 2theta angle
        intensity: intensities for all hkl planes
        pxrd: powder diffraction data
    """

    def __init__(self, filename):
        """Return a XRD object with the proper info"""
        self.from_file(filename)
        self.parse_cell()
        self.parse_atom()
        self.apply_symops()

    def from_file(self, filename):
        cif = np.genfromtxt(filename, dtype=str, delimiter='\n')
        
        # 3 modes in each flag:  
        # 0: not started; 
        # 1: reading; 
        # 2: done
        flags = {'cell':0, 'symops':0, 'atom':0}

        atom = {}
        cell = {}
        symops = {'string':[], 'matrix':[]}

        for lines in cif:

            if 'loop_' in lines:  
                #if a _loop lines starts, the current reading flag switch to 0
                for item in flags.keys():
                    if flags[item] == 1:
                        flags[item] = 2

            elif '_cell_length_' in lines or '_cell_angle_' in lines:
                #_cell_length_a          4.77985

                flags['cell'] = 1
                cell_str = lines.split()
                item = cell_str[0].replace(' ','')
                value = float(cell_str[1].split("(")[0])
                cell[item] = value

            elif '_symmetry_equiv_pos_as_xyz' in lines:
                #_symmetry_equiv_pos_as_xyz
                flags['symops'] = 1
      
            elif '_space_group_symop_operation_xyz' in lines:
                #_space_group_symop_operation_xyz
                flags['symops'] = 1
                
            elif flags['symops'] == 1:
                #1, 'x, y, z'
                #    x, -y, z
                raw_line = lines.strip().strip("'").split(' ', 1)
                if raw_line[0].isdigit():     
                    sym_str = raw_line[1].strip("'")
                else:
                    sym_str = lines.strip().strip("'").replace(' ', '')
                sym_str = sym_str.replace("'","")
                symops['string'].append(sym_str)
                symops['matrix'].append(self.xyz2sym_ops(sym_str))

            elif '_atom_site' in lines: 
                flags['atom'] = 1
                atom_str = lines.replace(' ','')
                item = atom_str
                atom[item] = []

            elif flags['atom'] == 1:
                raw_line = lines.split()
                for i, item in enumerate(atom.keys()):
                    raw_text = raw_line[i]
                    
                    if item.find('fract')>0:
                       value = float(raw_text.split("(")[0])
                    elif item.find('symbol')>0:
                       m_symbol = re.compile("([A-Z]+[a-z]*)")
                       value = str(m_symbol.findall(raw_text)).strip("[]").strip("''")
                       #print(raw_text, value)
                    else:
                       value = raw_text
                       
                    atom[item].append(value)

            elif flags['cell'] + flags['symops'] + flags['atom'] == 6:
                break

        self.cell = cell
        self.atom = atom
        self.symops = symops
   
    def parse_cell(self):
        cell_para = np.zeros(6)
        cell = self.cell
        for item in cell.keys():
            if item.find('_length_a') > 0:
                cell_para[0] = cell[item]
            elif item.find('_length_b') > 0:
                cell_para[1] = cell[item]
            elif item.find('_length_c') > 0:
                cell_para[2] = cell[item]
            elif item.find('_angle_alpha') > 0:
                cell_para[3] = np.radians(cell[item])
            elif item.find('_angle_beta') > 0:
                cell_para[4] = np.radians(cell[item])
            elif item.find('_angle_gamma') > 0:
                cell_para[5] = np.radians(cell[item])
        self.cell_para = cell_para

    def parse_atom(self):
        atom = self.atom
        N_atom = len(atom['_atom_site_fract_x'])
        cif_xyz = np.zeros([N_atom, 3])

        for item in atom.keys():
            if item.find('_fract_x') > 0:
                cif_xyz[:,0] = np.array(atom[item])
            elif item.find('_fract_y') > 0:
                cif_xyz[:,1] = np.array(atom[item])
            elif item.find('_fract_z') > 0:
                cif_xyz[:,2] = np.array(atom[item])

        self.cif_xyz = cif_xyz

    #generates all coordinates from rotation matrices and translation vectors
    def apply_symops(self):
        fract_xyz = self.cif_xyz
        symops_matrix = self.symops['matrix']
        atom_type = self.atom['_atom_site_type_symbol']
        sym_coordinates = {}
        
        for item in atom_type:
            sym_coordinates[item] = []


        for ii,item in enumerate(atom_type):
            for mat_vec in symops_matrix:
                sym_temp = np.dot(mat_vec[0], fract_xyz[ii].transpose()) + mat_vec[1]
                sym_coordinates[item].append(sym_temp)
        self.coordinate, self.composition, self.atom_type = \
                      self.remove_duplicate(sym_coordinates)

    #remove equivalent points and keep the unique ones
    #get the numbers of atoms per species
    @staticmethod
    def remove_duplicate(sym_coordinates):
        coordinate = []
        composition = []
        atom_type = []
        for item in sym_coordinates.keys():
            atom_type.append(item)
            raw_equiv = np.array(sym_coordinates[item])
            raw_equiv = raw_equiv - np.floor(raw_equiv)
            raw_equiv = np.around(raw_equiv, 4)
            raw_equiv = np.unique(raw_equiv, axis=0)
            composition.append(len(raw_equiv))
            if coordinate == []:
                coordinate = raw_equiv
            else:
                coordinate = np.concatenate((coordinate,raw_equiv),axis=0)

        return coordinate, composition, atom_type


    #function generates rotation matrices and translation vectors from equivalent points
    @staticmethod
    def xyz2sym_ops(string):
        #rotational matrix dictionary
        rot_dic = {}
        rot_dic['x'] = np.array([1.0,0,0])
        rot_dic['y'] = np.array([0,1.0,0])
        rot_dic['z'] = np.array([0,0,1.0])
        parts = string.strip().replace(' ','').lower().split(',')
        rot_mat = []
        rot_temp = np.array([0.,0.,0.])
        trans_vec = np.array([0.,0.,0.])
        #use re module to read xyz strings
        m_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
        m_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")
        for jj,item in enumerate(parts):
            #rotation matrix
            for ii,m in enumerate(m_rot.finditer(item)):
                coef = -1 if m.group(1) == '-' else 1
                if m.group(2) != '':
                    if m.group(3) != '':
                        coef *= float(m.group(2))/float(m.group(3))
                    else:
                        coef *= float(m.group(2))
                if ii == 0:                  
                    rot_temp = rot_dic[m.group(4)]*coef
                else:
                    rot_temp += rot_dic[m.group(4)]*coef
            rot_mat.append(rot_temp)
            #translation vector
            for m in m_trans.finditer(item):
                coef = -1 if m.group(1) == '-' else 1
                if m.group(3) != '':
                    coef = float(m.group(2))/float(m.group(3))
                else:
                    coef = float(m.group(2))
                trans_vec[jj] = 1.0*coef
        return (rot_mat, trans_vec)
         

class RDF(object):
    """a class of crystal structure. 
    Attributes:
        R_max: maximum cutoff distance
        R_bin: length of each bin when computing the RDF 
        width:
    """

    def __init__(self, crystal, R_max=12, R_bin=0.2, sigma=0.2):
        """Return a RDF object with the proper info"""
        self.R_max = R_max
        self.R_bin = R_bin
        self.sigma = sigma
        self.compute_RDF(crystal)
        self.plot_RDF()

    def compute_RDF(self, crystal):
        R_max = self.R_max
        R_bin = self.R_bin

        R = np.zeros(round(R_max/R_bin))
        R1 = self.RDF_in_box(crystal, R_max, R_bin)
        rij_dot = self.find_supercell(crystal, R_max)

        for atom in crystal.coordinate:
            origin = [np.dot(atom, crystal.cell_matrix)]
            #print(np.shape(origin))
            #print(np.shape(rij_dot))
            #print(rij_dot)
            rij_dist = cdist(rij_dot, origin)
            #print(rij_dist[rij_dist<10])
            for i in range(len(R)):
                d_min, d_max = R_bin*(i+0.5), R_bin*(i+1.5)
                R[i] += len([x for x in rij_dist if d_min<=x<d_max])

        self.RDF = np.ones([len(R), 2])
        rho = len(crystal.coordinate)/crystal.volume
        for i in range(len(R)):
            dr = R_bin
            r = (i+1)*dr
            self.RDF[i,0] = r
            self.RDF[i,1] = R[i]/(4*pi*dr*rho*r**2)
        #print(self.RDF)
        #print(R)

    def plot_RDF(self, filename=None):
        """ plot PXRD """
        datax = smear(self.RDF, self.sigma)
        #plt.plot(datax[:,0]/3.822, datax[:,1])
        plt.plot(datax[:,0], datax[:,1])
        ax=plt.gca()
        plt.grid()
        plt.xlabel(r"$r (\AA)$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)
        if filename is None:
           plt.show()
        else:
           plt.savefig(filename)
           plt.close()

    def RDF_in_box(self, crystal, R_max, R_bin):
        """ PDF for the box """
        R = np.zeros(round(R_max/R_bin))

        rij_dot = np.dot(crystal.coordinate, crystal.cell_matrix)
        for atom in crystal.coordinate:
            origin = [np.dot(atom, crystal.cell_matrix)]
            rij_dist = cdist(rij_dot, origin)
            for i in range(len(R)):
                d_min, d_max = R_bin*(i+0.5), R_bin*(i+1.5)
                R[i] += len([x for x in rij_dist if d_min<=x<d_max])

        return R/2

    def find_supercell(self, crystal, R_max, atom=[0.5, 0.5, 0.5]):
        # make supercell to cover the big sphere centered at atom with a radius of R_max
        # for all basic directions (1,0,0), (0,1,0), (1,1,0), (1,-1,0) and so on, 26 in total 
        hkl_max = np.array([1,1,1])
        for h1 in [-1, 0, 1]:
            for k1 in [-1, 0, 1]:
                for l1 in [-1, 0, 1]:
                    hkl_index = np.array([[h1,k1,l1]])
                    R = float(np.linalg.norm( np.dot(hkl_index-atom, crystal.cell_matrix), axis=1))
                    #print(hkl_index, R)
                    multiple = round(R_max/R)
                    hkl_index *= multiple
                    for i in range(len(hkl_max)):
                        if hkl_max[i] < hkl_index[0,i]:
                           hkl_max[i] = hkl_index[0,i]
        h, k, l = hkl_max
        coors = []
        for i in range(-h, h+1):
            for j in range(-k, k+1):
                for k in range(-l, l+1):
                    for coor in crystal.coordinate:
                        coors.append(coor+[i,j,k])
        supercell = np.asarray(coors)#.reshape(len(coors), 3)
        rij_dot = np.dot(supercell, crystal.cell_matrix)
        return rij_dot

from optparse import OptionParser
import pandas as pd
from tabulate import tabulate

if __name__ == "__main__":
    #-------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-c", "--crystal", dest="structure",default='',
                      help="crystal from file, cif or poscar, REQUIRED", metavar="crystal")
    parser.add_option("-r", "--Rmax", dest="Rmax", default=12, type='float',
                      help="Rmax, default: 12 A", metavar="Rmax")
    parser.add_option("-s", "--sigma", dest="sigma",default=0.20, type='float',
            help="sigma width for Gaussian smearing, default: 0.20", metavar="sigma")
    parser.add_option("-d", "--delta", dest="delta",default=0.08, type='float',
            help="step length, default: 0.08", metavar="R_bin")
    parser.add_option("-o", "--output", dest="mstyle",default='bmh',
            help="matplotlib style, fivethirtyeight, bmh, grayscale, dark_background, ggplot", metavar="mstyle")

    (options, args) = parser.parse_args()    
    if options.structure.find('cif') > 0:
       fileformat = 'cif'
    else:
       fileformat = 'POSCAR'

    plt.style.use(options.mstyle)
    test = crystal(fileformat, filename=options.structure)
    RDF(test, R_max=options.Rmax, R_bin=options.delta, sigma=options.sigma)   
