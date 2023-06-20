import tkinter as tk
import tkinter.filedialog
import tkinter.messagebox as msgbox
from tkinter import *
import numpy as np
import pandas as pd
import fnmatch
import os
import ntpath
from PIL import ImageTk,Image
import matplotlib.pyplot as plt

# diamagnetic corrections - dictionaries
# 1 - sample
# 2 - holder
# 3 - medium

# magnetic susceptibility for cations (-1x10^(-6) emu / mol)
cations = {'Ag1': 28, 'Ag2': 24, 'Al3': 2, 'As3': 9, 'As5': 6, 'Au1': 40, 'Au3': 32, 'B3': 0.2, 'Ba2': 26.5, 'Be2': 0.4,
           'Bi3': 25, 'Bi5': 23, 'Br5': 6, 'C4': 0.1, 'Ca2': 10.4, 'Cd2': 24, 'Ce3': 20, 'Ce4': 17, 'Cl5': 2, 'Co2': 12,
           'Co3': 10, 'Cr2': 15, 'Cr3': 11, 'Cr4': 8, 'Cr5': 5, 'Cr6': 3, 'Cs1': 35, 'Cu1': 12, 'Cu2': 11, 'Dy3': 19, 'Er3': 18,
          'Eu2': 22, 'Eu3': 20, 'Fe2': 13, 'Fe3': 10, 'Ga3': 8, 'Ge4': 7, 'Gd3': 20, 'H1': 0, 'Hf4': 16, 'Hg2': 40, 'Ho3': 19,
          'I5': 12, 'I7': 10, 'In3': 19, 'Ir1': 50, 'Ir2': 42, 'Ir3': 35, 'Ir4': 29, 'Ir5': 20, 'K1': 14.9, 'La3': 20, 'Li1': 1,
           'Lu3': 17, 'Mg2': 5, 'Mn2': 14, 'Mn3': 10, 'Mn4': 8, 'Mn6': 4, 'Mn7': 3, 'Mo2': 31, 'Mo3': 23, 'Mo4': 17, 'Mo5': 12,
           'Mo6': 7, 'N5': 0.1, '(NH4)1': 13.3, '(N(CH3)4)1': 52, '(N(C2H5)4)1': 101, 'Na1': 6.8, 'Nb5': 9, 'Nd3': 20, 'Ni2': 12,
           'Os2': 44, 'Os3': 36, 'Os4': 29, 'Os6': 18, 'Os8': 11, 'P3': 4, 'P5': 1, 'Pb2': 32, 'Pb4': 26, 'Pd2': 25,
           'Pd4': 18, 'Pm3': 27, 'Pr3': 20, 'Pr4': 18, 'Pt2': 40, 'Pt3': 3, 'Pt4': 28, 'Rb1': 22.5, 'Re3': 36, 'Re4': 28,
           'Re6': 16, 'Re7': 12, 'Rh3': 22, 'Rh4': 18, 'Ru3': 23, 'Ru4': 18, 'S4': 3, 'S6': 1, 'Sb3': 17, 'Sb5': 14,
           'Sc3': 6, 'Se4': 8, 'Se6': 5, 'Si4': 1, 'Sm2': 23, 'Sm3': 20, 'Sn2': 20, 'Sn4': 16, 'Sr2': 19, 'Ta5': 14, 'Tb3': 19,
           'Tb4': 17, 'Te4': 14, 'Te6': 12, 'Th4': 23, 'Ti3': 9, 'Ti4': 5, 'Tl1': 35.7, 'Tl3': 31, 'Tm3': 18,
          'U3': 46, 'U4': 35, 'U5': 26, 'U6': 19, 'V2': 15, 'V3': 10, 'V4': 7, 'V5': 4, 'VO2': 12.5, 'W2': 41, 'W3': 36,
           'W4': 23, 'W5': 19, 'W6': 13, 'Y3': 12, 'Yb2': 20, 'Yb3': 18, 'Zn2': 15, 'Zr4': 10}

# Atoms in covalent species (-1x10^(-6) emu / mol)
atomscov = {'Ag': 31, 'Al': 13, 'AsIII': 20.9, 'AsV': 43, 'B': 7, 'Bi': 192, 'Br': 30.6, 'C': 6, 'Cring': 6.24, 'Ca': 15.9,
            'Cl': 20.1, 'F': 6.3, 'H': 2.93, 'HgII': 33, 'I': 44.6, 'K': 18.5, 'Li': 4.2, 'Mg': 10, 'Nring': 4.61,
            'Nchain': 5.57, 'Na': 9.2, 'O': 4.6, 'P': 26.3, 'PbII': 46, 'S': 15, 'SbIII': 74, 'Se': 23, 'Si': 13,
            'SnIV': 30, 'Te': 37.3, 'TlI': 40, 'Zn': 13.5}

# specific bonds types
# d - doble, s - single, t - triple
# Ar - aromatic ring (1x10^(-6) emu / mol)
bonds = {'CdC': 5.5, 'CtC': 0, 'CdCsCdC': 10.6, 'ArsCtCsAr': 3.85, 'CH2dCHsCH2s(allyl)': 4.5, 'CdO': 6.3, 'COOH': -5,
         'COOR': -5, 'C(dO)NH2': -3.5, 'NdN': 1.85, 'CdNs': 8.15, 'sCtN': 0.8, 'sNtC': 0, 'NdO': 1.7, 'sNO2': -2, 'CsCl': 3.1,
        'ClsCR2CR2sCl': 4.3, 'R2CCl2': 1.44, 'RCHCL2': 6.43, 'CsBr': 4.1, 'BrsCR2CR2sBr': 6.24, 'CsI': 4.1, 'ArsOH': -1,
        'ArsNR2': 1, 'ArsC(dO)R': -1.5, 'ArsCOOR': -1.5, 'ArsCdC': -1, 'ArsCtC': -1.5, 'ArsOR': -1,'ArsCHO': -1.5,
         'ArsAr': -0.5, 'ArsNO2': -0.5, 'ArsBr': -3.5, 'ArsCl': -2.5, 'ArsI': -3.5, 'ArsCOOH': -1.5, 'ArsC(dO)NH2': -1.5,
        'R2CdNsNdCR2': 10.2, 'RCdCsC(dO)R': 0.8, 'Benzene': -1.4, 'Cyclobutane': 7.2, 'Cyclohexadiene': 10.56,
         'Cyclohexane': 3, 'Cyclohexene': 6.9, 'Cyclopentane': 0, 'Cyclopropane': 7.2, 'Dioxane': 5.5, 'Furan': -2.5,
        'Imidazole': 8, 'Isoxazole': 1, 'Morpholine': 5.5, 'Piperazine': 7, 'Piperidine': 3, 'Pyrazine': 9, 'Pyridine': 0.5,
        'Pyrimidine': 6.5, 'alpha-Pyrone': -1.4, 'gamma-pyrone': -1.4, 'Pyrrole': -3.5, 'Pyrrolidine': 0, 'Tetrahydrofuran': 0,
        'Thiazole': -3, 'Thiophene': -7, 'Triazine': -1.4}

#Anions (-1x10^(-6) emu / mol)
anions = {'(AsO3)3': 51, '(AsO4)3': 60, '(BF4)1': 37, '(BO3)3': 35, 'Br1': 34.6, '(BrO3)1': 40, 'Cl1': 23.4, '(Clo3)1': 30.2,
          '(Clo4)1': 32, 'CN1': 13, '(C5H5)1': 65, 'C6H5COO1': 71, '(CO3)2': 28, '(C2O4)2': 34, 'F1': 9.1, 'HCOO1': 17,
         'I1': 50.6, '(IO3)1': 51, '(IO4)1': 51.9, '(NO2)1': 10, '(NO3)1': 18.9, 'NCO1': 23, 'NCS1': 31, 'O2': 12, 'OAc1': 31.5,
         'OH1': 12, '(PO3)3': 42, '(PtCl6)2': 148, 'S2': 30, '(SO3)2': 38, '(SO4)2': 40.1, '(S2O3)2': 46, '(S2O8)2': 78,
          '(HSO4)1': 35, 'Se2': 48, '(SeO3)2': 44, '(SeO4)2': 51, '(SiO3)2': 36, 'Te2': 70, '(TeO3)2': 63, '(TeO4)2': 55}

# common ligands (-1x10^(-6) emu / mol)
ligands = {'Acac1': 52, 'Bipy': 105, 'CO': 10, '(C5H5)1': 65, 'En': 46.5, 'Ethylene': 15, 'Glycinate': 37, 'H20': 13,
           'Hydrazine': 20, 'Malonate': 45, 'NH3': 18, 'Phen': 128, 'o-PBMA': 194, 'Phthalocyanine': 442, 'PPh3': 167,
          'Pyrazine': 50, 'Pyridine': 49, 'Salen2': 182, 'Urea': 34}

# common solvents (-1x10^(-6) emu / mol)
solvents = {'CCl4': 66.8, 'CHCl3': 58.9, 'CH2Cl2': 46.6, 'CH3Cl': 32, 'CH3NO2': 21, 'CH3OH': 21.4, 'CCl3COOH': 73,
            'CF3COOH': 43.3, 'CH3CN': 27.8, '1,2-C2H4Cl2': 59.6, 'CH3COOH': 31.8, 'CH3CH2OH': 33.7, 'HOCH2CH2OH': 38.9,
           'CH3CH2SH': 44.9, 'CH3C(dO)CH3': 33.8, 'CH3C(dO)OC(dO)CH3': 52.8, 'CH3CH2CH2CN': 50.4, 'CH3C(dO)OCH2CH3': 54.1,
           'CH3CH2CH2CH2OH': 56.4, 'CH3CH2OCH2CH3': 55.5, 'Pentane': 61.5, 'o-Dichlorobenzene': 84.4, 'Benzene': 54.8,
           'Cyclohexane': 68, 'Hexane': 74.1, 'Triethyloamine': 83.3, 'Benzonitrile': 65.2, 'Toluene': 65.6,
            'Isooctane': 99.1, 'Naphthalene': 91.6}

#holders (emu / g)
holders = {'capsule': -0.4074*0.000000001, 'delrin': -16.46*0.000000001/30, 'PET': -11.86*0.000000001/14}

#nujol, wosk (emu / g)
nujol = {'nujol': -11.86*0.000000001/14, 'wosk': -243.06*0.000000001/282.55}

# function returning a path to a file
def path(self):
    start, end = ntpath.split(self)
    return end or ntpath.basename(start)

# definition of lists for calculating diamgnetic correcctions for dc data
DC = []
all_corr = []

class Window(tk.Tk):
    """Class for magnetic data processing with graphic interface"""

    def __init__(self):
        super().__init__()
        self.title('magnetic data processing')
        self.label_text = tk.StringVar()
        self.label_text2 = tk.StringVar()
        self.label_text3 = tk.StringVar()
        self.label_text.set('Magnetic data processing')

        self.label = tk.Label(self, textvar=self.label_text)
        self.label.pack(fill=tk.BOTH, expand=1, padx=100, pady=10)

        # click the button to choose the file for processing
        b1 = tk.Button(self, text='Choose file', command=self.choose_file)
        b1.pack(fill='x')

    def choose_file(self):
        global f
        global b4
        # choosing a file for processing
        f = tkinter.filedialog.askopenfilename(
            parent=window,
            initialdir='/',
            title='Choose file for processing.',
            filetypes=[('dat files', '.dat')])

        ntpath.basename(f)

        filename = path(f)
        if fnmatch.fnmatch(path(filename), '*.ac.dat'):
            # ac data
            # show the name of the currently processed file
            self.label_text.set(filename)
            # selction of needed columns and rows
            ac = pd.read_table(f, skiprows=28, header=None, engine='python',
                               sep='\s\s+|,', usecols=[2, 3, 4, 5, 6, 7, 8, 10, 13, 14])
            # just data
            ac2 = pd.DataFrame(ac)
            # information from comments (mass, molar mmass)
            data3 = pd.read_table(f, skiprows=[i for i in range(0, 8)] + [j for j in range(11, 31)],
                                  header=None, engine='python', sep='\s\s+|,', usecols=[2])
            # number of moles
            n = data3.iat[0, 0] * 0.001 / data3.iat[1, 0]
            # information from comments (holders etc.)
            data4 = pd.read_table(f, skiprows=[i for i in range(0, 12)] +
                                              [j for j in range(14, 31)], header=None, engine='python', sep='\s\s+|,',
                                  usecols=[2, 3])
            # the form of the output file: 1. row - headers, 2. row - units, 3. row - comments
            ac1 = pd.DataFrame(
                [['', '', 'H', 'T', 'm1', 'dm1', 'm2', 'dm2', 'amp', '', 'phase', '', '', 'Hamp', 'f', 'chi1', 'chi2',
                  'Tround'],
                 ['', '', 'Oe', 'K', 'emu', 'emu', 'emu', 'emu', 'emu', '', 'deg', '', '', 'Oe', 'Hz', 'cm3/mol',
                  'cm3/mol',
                  'K'],
                 ['', '', data3.iat[0, 0], n, data4.iat[0, 0], data4.iat[0, 1], '', '', '', '', '', '', '',
                  '', '', '']])
            # claculation of chi', chi'' (susceptibility) and round temperature
            ac2[15] = ac2[4] / (ac2[13] * n)
            ac2[16] = ac2[6] / (ac2[13] * n)
            ac2[17] = ac2[3].round(2)
            # data for file construction
            ac3 = ac1.append(ac2)

            ac4 = ac3.drop([0, 1, 9, 11, 12], axis=1)
            # saving a file
            ac4.to_csv(f + '.txt', sep=',', index=False, header=0)
            # click the button to save the pictures with graps
            # link to the function for drawing graphs and saving them as .png files
            b4 = tk.Button(self, text='Draw graphs', command=self.draw_graphs)
            b4.pack(fill='x')

        elif fnmatch.fnmatch(path(filename), '*.dc.dat') or fnmatch.fnmatch(path(filename), '*.rso.dat'):
            # dc data
            global n1
            global dc1
            global dc2
            global dc3
            global dc4
            # show the name of the currently processed file
            self.label_text.set(filename)
            # selction of needed columns and rows
            dc = pd.read_table(f, skiprows=31, header=None, engine='python', sep='\s\s+|,', usecols=[2, 3, 4, 5])
            # just data
            dc2 = pd.DataFrame(dc)
            # information from comments (mass, molar mmass)
            data3 = pd.read_table(f, skiprows=[i for i in range(0, 8)] +
                                              [j for j in range(10, 31)], header=None, engine='python', sep='\s\s+|,',
                                  usecols=[2])
            # number of moles
            n1 = data3.iat[0, 0] * 0.001 / data3.iat[1, 0]
            # information from comments (holders etc.)
            data4 = pd.read_table(f, skiprows=[i for i in range(0, 12)] +
                                              [j for j in range(14, 31)], header=None, engine='python', sep='\s\s+|,',
                                  usecols=[2, 3])
            # the form of the output file: 1. row - headers, 2. row - units, 3. row - comments
            dc1 = pd.DataFrame([['', '', 'H', 'T', 'm', 'dm', 'chi', 'chiT', 'M'],
                                ['', '', 'Oe', 'K', 'emu', 'emu', 'cm3/mol', 'cm3K/mol', 'uB/mol'],
                                ['', '', data3.iat[0, 0], n1, data4.iat[0, 0], data4.iat[0, 0], '', '', '']])
            # link to the function asking for diamagnetic corrections
            self.diamagnetic_corrections()

        else:
            # in case of other files, chosen by accident, the information of incorrect file displayed
            self.label_text.set('Incorrect file.')

    def diamagnetic_corrections(self):
        # box asking if the user want to include diamagnetic corrections
        if msgbox.askyesno('Diamagnetic corrections', 'Do you want to include diamagnetic corrections?'):
            # if the user want to include diamagnetic corrections, the image with symbols is displayed
            # (keys to the dictionaries with values of diamagentic corrections)
            global pic
            global name_entry
            global b3
            global b6

            load = Image.open("keys.png")
            render = ImageTk.PhotoImage(load)
            pic = Label(self, image=render)
            pic.image = render
            pic.pack(side="top", fill="both", expand="yes")
            # information how to enter the components of the diamagentic correction
            self.label_text2.set('Enter the sample components: cations, atoms in covalent species,\n \
            specific bonds, anions ligands or solvents, holder materials or another medium (if was added) \
            \n Nomenclature: (sample item) (number of sample items) (holder material) (mass of holder in mg) \n \
            (medium) (mass of medium in mg). \n \
            Example: Co2 1 NCS1 2 Pyridine 2 C 2 H 6')

            self.name_text = tk.StringVar()

            self.label = tk.Label(self, textvar=self.label_text2)
            self.label.pack(fill=tk.BOTH, expand=1, padx=100, pady=10)
            # field to enter the components
            self.name_entry = tk.Entry(self, textvar=self.name_text)
            self.name_entry.pack(fill=tk.BOTH, expand=1, padx=20, pady=10)
            # click the button to calculate the diamagnetic correction,
            # link to the function calculating diamagnetic corrections on the basis of the entered components
            b3 = tk.Button(self, text='Calculate', command=self.diamag_corr)
            b3.pack(fill='x')

        else:
            # if the user does not want to include the diamagnetic correction,
            # ask if the user want to include ferromagnetic correction
            msgbox.showinfo('Ok', 'The data will not be corrected for diamagnetism')

            if msgbox.askyesno('Ferromagnetic corrections', 'Do you want to include ferromagnetic corrections?'):
                # if yes, ask for file to calculate the ferromagnetic correction
                # the measurement of magnetic moment vs. field in high tmeperature is needed
                # (only for samples displaying the magnetic properties in low temperature)
                self.label_text2.set('Choose another file to calculate ferromagnetic corrections.')
                self.label = tk.Label(self, textvar=self.label_text2)
                self.label.pack(fill=tk.BOTH, expand=1, padx=100, pady=10)
                # a path to the function to calculate the ferromagnetic correction and include it in the data
                self.ferro_path()

            else:
                msgbox.showinfo('Ok', 'The data will not be corrected for ferromagnetic impurities.')
                # if no, no diamagnetic or ferromagnetic correction will be included.
                # calulation of chi (susceptibility), chiT, and M (magnetisation)
                dc2[6] = dc2[4] / (dc2[2] * n1)
                dc2[7] = dc2[6] * dc2[3]
                dc2[8] = dc2[4] / (n1 * 5580.54)
                # data for file construction
                dc3 = dc1.append(dc2)

                dc4 = dc3.drop([0, 1], axis=1)
                # saving a file
                dc4.to_csv(f + '.txt', sep=',', index=False, header=0)
                # click the button to save the pictures with graps
                # link to the function for drawing graphs and saving them as .png files
                b6 = tk.Button(self, text='Draw graphs', command=self.draw_graphs)
                b6.pack(fill='x')

    def diamag_corr(self):
        const = 0
        global b5
        # function calculating diamagnetic corrections as diamagnetic susceptibility
        # (to obtain it in emu multiple by filed in Oe)
        while True:
            diamagwhich = self.name_entry.get()
            diamagwhich = diamagwhich.split()
            La = []
            Lb = []
            # division the list of entries into two lists: first one with components, second one with the quantities
            for i in diamagwhich:
                if diamagwhich.index(i) % 2 == 0:
                    La.append(i)
                if diamagwhich.index(i) % 2 != 0:
                    Lb.append(i)
            try:
                L = []
                for i in La:
                    if i in cations.keys():
                        x = cations[i] * -1 * 0.000001 * float(Lb[La.index(i)]) * n1
                        L.append(x)
                        DC.append(sum(L))
                    elif i in atomscov.keys():
                        x = atomscov[i] * -1 * 0.000001 * float(Lb[La.index(i)]) * n1
                        L.append(x)
                        DC.append(sum(L))
                    elif i in bonds.keys():
                        x = bonds[i] * 0.000001 * float(Lb[La.index(i)]) * n1
                        L.append(x)
                        DC.append(sum(L))
                    elif i in anions.keys():
                        x = anions[i] * -1 * 0.000001 * float(Lb[La.index(i)]) * n1
                        L.append(x)
                        DC.append(sum(L))
                    elif i in ligands.keys():
                        x = ligands[i] * -1 * 0.000001 * float(Lb[La.index(i)]) * n1
                        L.append(x)
                        DC.append(sum(L))
                    elif i in solvents.keys():
                        x = solvents[i] * -1 * 0.000001 * float(Lb[La.index(i)]) * n1
                        L.append(x)
                        DC.append(sum(L))
                    elif i in holders.keys():
                        x = holders[i] * float(Lb[La.index(i)])
                        L.append(x)
                        DC.append(sum(L))
                    elif i in nujol.keys():
                        x = nujol[i] * float(Lb[La.index(i)])
                        L.append(x)
                        DC.append(sum(L))
                    else:
                        const = 1
                # the diamagnetic corrections are summed up, the sum is appended to the all_coor list
                # there is always the last value on the list taken for corrections -
                # in case the user entered incorrectly the components, it is always possibility to enter them
                # again and recalculate the correction.
                all_corr.append(DC[-1])
                if const == 0:
                    self.label_text2.set(all_corr[-1])
                    # to accept the correction click accept
                    # link to the function asking for ferromagnetic corrections
                    b5 = tk.Button(self, text='Accept', command=self.ferro)
                    b5.pack(fill='x')
                    break
                else:
                    self.label_text2.set('Something went wrong')
                    if b5.winfo_exists() == 1:
                        b5.destroy()
                    break
            # errors handling - just in case
            except KeyError:
                self.label_text2.set('Key Error')
                break
            except ValueError:
                self.label_text2.set('Value Error')
                break
            except IndexError:
                self.label_text2.set('Index Error')
                break

    def ferro(self):
        global b6
        # function asking if the user wants to include ferromagnetic correction
        # the unnecessary elements in the interface disappear
        pic.destroy()
        self.name_entry.destroy()
        b3.destroy()
        b5.destroy()
        if msgbox.askyesno('Ferromagnetic corrections', 'Do you want to include ferromagnetic corrections?'):
            # if the users wants to include ferromagnetic corrections, the choice of the appropriate file is necessary
            self.label_text3.set('Choose another file to calculate ferromagnetic corrections.')
            self.label = tk.Label(self, textvar=self.label_text3)
            self.label.pack(fill=tk.BOTH, expand=1, padx=100, pady=10)
            # function calculating ferromagnetic correction and simoultaneuosly including
            # both dia- and ferromagnetic corrections to data
            self.ferro_path2()
        else:
            # if the user does not want to include ferromagnetic correction,
            # the data will be correctoed only for diamagnetism
            msgbox.showinfo('Ok', 'The data will not be corrected for ferromagnetic impurities.')
            dc2[6] = (dc2[4] - all_corr[-1] * dc2[2]) / (dc2[2] * n1)
            dc2[7] = dc2[6] * dc2[3]
            dc2[8] = (dc2[4] - all_corr[-1] * dc2[2]) / (n1 * 5580.54)

            dc3 = dc1.append(dc2)

            dc4 = dc3.drop([0, 1], axis=1)

            dc4.to_csv(f + '.txt', sep=',', index=False, header=0)
            # click the button to draw graphs
            b6 = tk.Button(self, text='Draw graphs', command=self.draw_graphs)
            b6.pack(fill='x')

    def ferro_path(self):
        # ferromagnetic correction calculated in emu
        global f2
        global b6
        f2 = tkinter.filedialog.askopenfilename(
            parent=window,
            initialdir='C:/Users/Magdalena Ceglarska/Desktop',
            title='Choose file',
            filetypes=[('dat files', '.dat')])
        # the file chosen for calculating ferromagnetic correction is converted the same the dc data
        # in case the user want to use the file for other purposes. It is also saved to .txt file
        ntpath.basename(f2)

        filename = path(f2)

        fr = pd.read_table(f2, skiprows=31, header=None, engine='python', sep='\s\s+|,', usecols=[2, 3, 4, 5])

        fr2 = pd.DataFrame(fr)

        fr3 = pd.read_table(f2, skiprows=[i for i in range(0, 8)] +
                                         [j for j in range(10, 31)], header=None, engine='python', sep='\s\s+|,',
                            usecols=[2])

        n2 = fr3.iat[0, 0] * 0.001 / fr3.iat[1, 0]

        fr4 = pd.read_table(f2, skiprows=[i for i in range(0, 12)] +
                                         [j for j in range(14, 31)], header=None, engine='python', sep='\s\s+|,',
                            usecols=[2, 3])

        fr1 = pd.DataFrame([['', '', 'H', 'T', 'm', 'dm', 'chi', 'chiT', 'M']])

        fr2[6] = (fr2[4]) / (fr2[2] * n2)
        fr2[7] = fr2[6] * fr2[3]
        fr[8] = (fr2[4]) / (n2 * 5580.54)

        fr3 = fr1.append(fr2)

        fr4 = fr3.drop([0, 1], axis=1)

        fr4.to_csv(f2 + '.txt', sep=',', index=False, header=0)

        graphdata = pd.read_table(f2 + '.txt', engine='python',
                                  sep='\s\s+|,')
        # To calculate the ferromagnetic correction the M(H) data
        # (few point at high, positive fields, few points at high, negative fields) are fitted.
        x = graphdata['H']
        y = graphdata['m']

        X = np.array(x)
        Y = np.array(y)
        x1 = X[0:4]
        y1 = Y[0:4]
        x2 = X[len(X) - 4:len(X)]
        y2 = Y[len(Y) - 4:len(Y)]
        a, b = np.polyfit(x1, y1, 1)
        c, d = np.polyfit(x2, y2, 1)

        global ferro
        # From the intercepts (b and d) the correction is calculated.
        ferro = np.absolute(b - d) / 2
        self.label_text2.set(ferro)
        # only ferromagnetic correction included to the data
        dc2[6] = (dc2[4] - ferro) / (dc2[2] * n1)
        dc2[7] = dc2[6] * dc2[3]
        dc2[8] = (dc2[4] - ferro) / (n1 * 5580.54)

        dc3 = dc1.append(dc2)

        dc4 = dc3.drop([0, 1], axis=1)

        dc4.to_csv(f + '.txt', sep=',', index=False, header=0)
        # click the button to draw graphs
        b6 = tk.Button(self, text='Draw graphs', command=self.draw_graphs)
        b6.pack(fill='x')

    def ferro_path2(self):
        # the same as ferro_path, but in the data both diamagnetic and ferromagnetic corrections are included
        global f2
        global b6
        f2 = tkinter.filedialog.askopenfilename(
            parent=window,
            initialdir='/',
            title='Choose file',
            filetypes=[('dat files', '.dat')])

        ntpath.basename(f2)

        filename = path(f2)

        fr = pd.read_table(f2, skiprows=31, header=None, engine='python', sep='\s\s+|,', usecols=[2, 3, 4, 5])

        fr2 = pd.DataFrame(fr)

        fr3 = pd.read_table(f2, skiprows=[i for i in range(0, 8)] +
                                         [j for j in range(10, 31)], header=None, engine='python', sep='\s\s+|,',
                            usecols=[2])

        n2 = fr3.iat[0, 0] * 0.001 / fr3.iat[1, 0]

        fr4 = pd.read_table(f2, skiprows=[i for i in range(0, 12)] +
                                         [j for j in range(14, 31)], header=None, engine='python', sep='\s\s+|,',
                            usecols=[2, 3])

        fr1 = pd.DataFrame([['', '', 'H', 'T', 'm', 'dm', 'chi', 'chiT', 'M']])

        fr2[6] = (fr2[4]) / (fr2[2] * n2)
        fr2[7] = fr2[6] * fr2[3]
        fr[8] = (fr2[4]) / (n2 * 5580.54)

        fr3 = fr1.append(fr2)

        fr4 = fr3.drop([0, 1], axis=1)

        fr4.to_csv(f2 + '.txt', sep=',', index=False, header=0)

        graphdata = pd.read_table(f2 + '.txt', engine='python',
                                  sep='\s\s+|,')

        x = graphdata['H']
        y = graphdata['m']

        X = np.array(x)
        Y = np.array(y)
        x1 = X[0:4]
        y1 = Y[0:4]
        x2 = X[len(X) - 4:len(X)]
        y2 = Y[len(Y) - 4:len(Y)]
        a, b = np.polyfit(x1, y1, 1)
        c, d = np.polyfit(x2, y2, 1)

        global ferro1
        ferro1 = np.absolute(b - d) / 2
        self.label_text3.set(ferro1)
        dc2[6] = (dc2[4] - all_corr[-1] * dc2[2] - ferro1) / (dc2[2] * n1)
        dc2[7] = dc2[6] * dc2[3]
        dc2[8] = (dc2[4] - all_corr[-1] * dc2[2] - ferro1) / (n1 * 5580.54)

        dc3 = dc1.append(dc2)

        dc4 = dc3.drop([0, 1], axis=1)

        dc4.to_csv(f + '.txt', sep=',', index=False, header=0)

        b6 = tk.Button(self, text='Draw graphs', command=self.draw_graphs)
        b6.pack(fill='x')

    def draw_graphs(self):
        # function for drawing the appropriate graphs, according to the identificator of data in the file name
        if (fnmatch.fnmatch(path(f + '.txt'), '*h1000*') or fnmatch.fnmatch(path(f + '.txt'), '*h5000*') or \
            fnmatch.fnmatch(path(f + '.txt'), '*h100*') or fnmatch.fnmatch(path(f + '.txt'), '*h2000*')) and \
                (fnmatch.fnmatch(path(f + '.txt'), '*dc*') or fnmatch.fnmatch(path(f + '.txt'), '*rso*')):
            graphdata = pd.read_table(f + '.txt', engine='python', skiprows=[i for i in range(1, 3)],
                                      sep='\s\s+|,')

            x = np.array(graphdata['T'])
            y = np.array(graphdata['chi'])
            y2 = np.array(graphdata['chiT'])

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
            ax1.plot(x, y, 'o')
            ax2.plot(x, y2, 'o')

            ax1.set(xlabel='T (K)', ylabel='$\chi$ (cm$^{3}$/mol)')
            ax2.set(xlabel='T (K)', ylabel='$\chi$T (cm$^{3}$K/mol)')

            plt.savefig(f + '.png')
            b6.destroy()
            self.label_text2.set('')
            self.label_text3.set('')

        if (fnmatch.fnmatch(path(f + '.txt'), '*t1k8*') or fnmatch.fnmatch(path(f + '.txt'), '*t2k*') or \
            fnmatch.fnmatch(path(f + '.txt'), '*t2k5*') or fnmatch.fnmatch(path(f + '.txt'), '*t3k*') or \
            fnmatch.fnmatch(path(f + '.txt'), '*mht*')) and \
                (fnmatch.fnmatch(path(f + '.txt'), '*dc*') or fnmatch.fnmatch(path(f + '.txt'), '*rso*')):
            graphdata = pd.read_table(f + '.txt', engine='python', skiprows=[i for i in range(1, 3)],
                                      sep='\s\s+|,')

            x = np.array(graphdata['H'])
            y = np.array(graphdata['M'])

            fig, ax1 = plt.subplots(1, figsize=(15, 5))
            ax1.plot(x, y, 'o')

            ax1.set(xlabel='H (Oe)', ylabel='M ($\mu$$_{B}$)')

            plt.savefig(f + '.png')
            b6.destroy()
            self.label_text2.set('')
            self.label_text3.set('')

        if fnmatch.fnmatch(path(f + '.txt'), '*zfcfc*') and \
                (fnmatch.fnmatch(path(f + '.txt'), '*dc*') or fnmatch.fnmatch(path(f + '.txt'), '*rso*')):
            graphdata = pd.read_table(f + '.txt', engine='python', skiprows=[i for i in range(1, 3)],
                                      sep='\s\s+|,')

            x = np.array(graphdata['T'])
            y = np.array(graphdata['chi'])

            fig, ax1 = plt.subplots(1, figsize=(15, 5))
            ax1.plot(x, y, 'o')

            ax1.set(xlabel='T (K)', ylabel='$\chi$ (cm$^{3}$/mol)')

            plt.savefig(f + '.png')
            b6.destroy()
            self.label_text2.set('')
            self.label_text3.set('')

        if fnmatch.fnmatch(path(f + '.txt'), '*act*') and \
                fnmatch.fnmatch(path(f + '.txt'), '*ac*'):
            graphdata = pd.read_table(f + '.txt', engine='python', skiprows=[i for i in range(1, 3)],
                                      sep='\s\s+|,')

            x = np.array(graphdata['T'])
            y = np.array(graphdata['chi1'])
            y2 = np.array(graphdata['chi2'])

            fig, ax1 = plt.subplots(1, figsize=(15, 5))
            ax1.plot(x, y, 'o')
            ax1.plot(x, y2, 'o')

            ax1.set(xlabel='T (K)', ylabel='$\chi$$_{AC}$ (cm$^{3}$/mol)')

            plt.savefig(f + '.png')
            b6.destroy()
            self.label_text2.set('')
            self.label_text3.set('')

        if fnmatch.fnmatch(path(f + '.txt'), '*acf-h*') and \
                fnmatch.fnmatch(path(f + '.txt'), '*ac*'):
            graphdata = pd.read_table(f + '.txt', engine='python', skiprows=[i for i in range(1, 3)],
                                      sep='\s\s+|,')
            x = np.array(graphdata['T'])
            y = np.array(graphdata['chi1'])
            y2 = np.array(graphdata['chi2'])

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
            ax1.plot(x, y, 'o')
            ax1.plot(x, y2, 'o')
            ax2.plot(y, y2, 'o')

            ax1.set(xlabel='T (K)', ylabel='$\chi$$_{AC}$ (cm$^{3}$/mol)')
            ax2.set(xlabel='$\chi$\' (cm$^{3}$/mol)', ylabel='$\chi$\'\' (cm$^{3}$/mol)')

            plt.savefig(f + '.png')
            b4.destroy()

        if fnmatch.fnmatch(path(f + '.txt'), '*acf-t*') and \
                fnmatch.fnmatch(path(f + '.txt'), '*ac*'):
            graphdata = pd.read_table(f + '.txt', engine='python', skiprows=[i for i in range(1, 3)],
                                      sep='\s\s+|,')
            x = np.array(graphdata['H'])
            y = np.array(graphdata['chi1'])
            y2 = np.array(graphdata['chi2'])

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
            ax1.plot(x, y, 'o')
            ax1.plot(x, y2, 'o')
            ax2.plot(y, y2, 'o')

            ax1.set(xlabel='H (Oe)', ylabel='$\chi$$_{AC}$ (cm$^{3}$/mol)')
            ax2.set(xlabel='$\chi$\' (cm$^{3}$/mol)', ylabel='$\chi$\'\' (cm$^{3}$/mol)')

            plt.savefig(f + '.png')
            b4.destroy()

# After the graphs are drawn, the program is ready to use it for another file (choose another file for processing)
if __name__ == "__main__":
    window = Window()
    close = Button(window, text='Exit', command=window.quit)
    close.pack()
    window.mainloop()
