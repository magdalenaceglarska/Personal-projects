import os
import pandas as pd
import matplotlib.pyplot as plt
import unittest
import os.path
from os import path
# path to the file
filename = os.path.abspath(
    "C:/Users/magda/Documents/GitHub/Personal-projects/heat capacity project/HC_AJ1396_Co-morph.dat")

# selction of needed columns and rows
# hc = heat capacity
hc = pd.read_table(filename, skiprows=15, header=None, engine='python',
                   sep='\s\s+|,', usecols=[5, 7, 9, 11])
# averaging
for i in range(0, len(hc)):
    if i % 2 == 0:
        hc.loc[i / 2, 5] = (hc[5][i] + hc[5][i + 1]) / 2
        hc.loc[i / 2, 7] = (hc[7][i] + hc[7][i + 1]) / 2
        hc.loc[i / 2, 9] = (hc[9][i] + hc[9][i + 1]) / 2
        hc.loc[i / 2, 11] = (hc[11][i] + hc[11][i + 1]) / 2
hc = hc.drop([i for i in range(int(len(hc) / 2), len(hc))], axis=0)

# data for mass, molar mass
hc1 = pd.read_table(filename, skiprows=[i for i in range(6, len(hc) + 15)]
                    + [j for j in range(0, 3)] + [k for k in range(4, 5)],
                    header=None, engine='python', sep='\s\s+|,',
                    usecols=[1])
# number of moles
n = hc1.iat[0, 0] * 0.001 / hc1.iat[1, 0]

# the form of the output file: 1. row - headers, 2. row - units, 3. row - comments

hc2 = pd.DataFrame([['', '', '', '', '', 'H', '', 'T', '', 'HC', '', 'Add', 'C', 'C/T'],
                    ['', '', '', '', '', 'Oe', '', 'K', '', 'uJ/K', '', 'uJ/K', 'J/molK', 'J/molK2'],
                    ['', '', '', '', '', 'm = ' + str(hc1.iat[0, 0]) + ' mg', '',
                     'M = ' + str(hc1.iat[1, 0]) + ' g/mole',
                     '', str(n) + ' mole', '', 'HC avg', '', '']])
# calculation of C, C/T
hc[12] = hc[9] * 0.000001 / n
hc[13] = hc[12] / hc[7]
# construction of the output data for file
hc3 = hc2.append(hc)
hc3 = hc3.drop([0, 1, 2, 3, 4, 6, 8, 10], axis=1)
hc3.to_csv(filename + '.txt', sep=',', index=False, header=0)
# data for graph to remove outliers
graphdata = pd.read_table(filename + '.txt', engine='python', skiprows=[i for i in range(1, 3)],
                          sep='\s\s+|,')

X = graphdata['T']
Y = graphdata['C']
Y1 = graphdata['C/T']

fig, ax1 = plt.subplots(1, figsize=(15, 15))
ax1.plot(X, Y, 'o', picker=True, pickradius=2)
ax1.set(xlabel='T (K)', ylabel='C (J/(mol K))')
plt.title('Pick outliers')
plt.savefig(filename + '.png')

picked = []

# function to pick a point
def out_points(point):
    line = point.artist
    xdata = line.get_xdata()
    ydata = line.get_ydata()
    ind = point.ind
    point.canvas.draw()
    picked.append(tuple(zip(xdata[ind], ydata[ind])))
    points = list(zip(xdata[ind], ydata[ind]))
    print(points)


fig.canvas.mpl_connect('pick_event', out_points)
plt.show()
# rows to remove
to_remove = []
for i in range(0,len(picked)):
    for j in range(0,len(picked[i])):
        to_remove.append(hc[hc[7] == picked[i][j][0]].index[0])
print(to_remove)
# removing
#if i in hc.index:
hc = hc.drop(to_remove, axis=0)
#else:
#    pass
# construction of the output data without outliers for file
hc3 = hc2.append(hc)
hc3 = hc3.drop([0, 1, 2, 3, 4, 6, 8, 10], axis=1)
hc3.to_csv(filename + '-out' + '.txt', sep=',', index=False, header=0)
# data for graphs drawing: C(T) and C/T(T)
graphdata1 = pd.read_table(filename + '-out' + '.txt', engine='python', skiprows=[i for i in range(1, 3)],sep='\s\s+|,')

X = graphdata1['T']
Y = graphdata1['C']
Y1 = graphdata1['C/T']

fig, ax1 = plt.subplots(2, figsize=(15, 15))
ax1[0].plot(X, Y, 'o')
ax1[0].set(xlabel='T (K)', ylabel='C (J/(mol K))')
ax1[1].plot(X, Y1, 'o')
ax1[1].set(xlabel='T (K)', ylabel='C/T (J/(mol K2))')
plt.savefig(filename + '-out' + '.png')

plt.show()
# function to find the length of the file - for tests
def file_length(file):
    with open(file) as f:
        for i, k in enumerate(f):
            pass
    return i + 1


class TestMethods(unittest.TestCase):

    def test_file(self):
        self.assertEqual(path.exists(filename + '.txt'), True)
    def test_fileout(self):
        self.assertEqual(path.exists(filename + '-out' + '.txt'), True)
    def test_png(self):
        self.assertEqual(path.exists(filename + '.png'), True)
    def test_pngout(self):
        self.assertEqual(path.exists(filename + '-out' + '.png'), True)

    def test_file_empty(self):
        fsize = os.path.getsize(filename + '.txt')
        self.assertEqual(fsize != 0, True)

    def test_file_empty2(self):
        fsize2 = os.path.getsize(filename + '-out' + '.txt')
        self.assertEqual(fsize2 != 0, True)

    def test_file_empty_png(self):
        psize = os.path.getsize(filename + '.png')
        self.assertEqual(psize != 0, True)

    def test_file_empty_png2(self):
        psize2 = os.path.getsize(filename + '.png')
        self.assertEqual(psize2 != 0, True)

    def test_null(self):
        self.assertEqual(hc3.isnull().all().all(), False)

    def test_nan(self):
        self.assertEqual(hc3[3::].isna().all().all(), False)

    def test_columns(self):
        data = pd.read_table(filename + '.txt', engine='python', sep='\s\s+|,')
        self.assertEqual(len(data.columns), 6)

    def test_rows(self):
        data = pd.read_table(filename + '.txt', engine='python', sep='\s\s+|,')
        self.assertEqual(len(data), (file_length(filename)-15)/2+3-1)

    def test_columns2(self):
        data = pd.read_table(filename + '-out' + '.txt', engine='python', sep='\s\s+|,')
        self.assertEqual(len(data.columns), 6)

    def test_rows2(self):
        data = pd.read_table(filename + '-out' + '.txt', engine='python', sep='\s\s+|,')
        self.assertEqual(len(data), (file_length(filename)-15)/2+3-1-len(set(to_remove)))

    #def test_values(self):
    #    data = pd.read_table(filename + '-out' + '.txt', engine='python', sep='\s\s+|,')
    #    self.assertEqual(round(float(data.iloc[2,4]), 4), 106.8183)

    #def test_values2(self):
    #    data = pd.read_table(filename + '-out' + '.txt', engine='python', sep='\s\s+|,')
    #    self.assertEqual(round(float(data.iloc[4,5]), 4), 1.7092)


if __name__ == '__main__':
    unittest.main()