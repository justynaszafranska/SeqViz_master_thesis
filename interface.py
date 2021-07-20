import wx
import wx.grid as grid
from Bio import AlignIO
import os
import matplotlib.pyplot as plt 
from dicts import *
import numpy as np


def Exit(evt):
    dialog = wx.MessageDialog(Window1, 'Are you sure?', 'Closing the program',
                              style = wx.OK | wx.CANCEL)
    x = dialog.ShowModal()
    dialog.Destroy()
    if x == wx.ID_OK:
        Window1.Close()

def Read_Fasta(evt):
    global n
    global path_file
    global new
    global id_seq
    new = []
    Names = []
    id_seq = []
    path_file = []
    dialog = wx.FileDialog(Window1, message = 'Select .fasta/.fa file',
                           defaultFile = '', wildcard = '*.FASTA',
                           style = wx.FD_OPEN, pos=(10,10))
    if dialog.ShowModal() == wx.ID_OK:
        selected_file = dialog.GetPaths()
        file = open(selected_file[0],'r')
        align = AlignIO.read(file, "fasta")
        sequences = []
        desc = [] #description
        for x in align:
            sequences.append(x.seq)
            id_seq.append(x.id)

            desc.append(x.description)
            
        new = (list(map(str,sequences)))
        path_file = os.path.basename(file.name)
        Names.append(path_file)
        List.InsertItems(Names,n)
        n=n+1
        List.Show()
        dialog.Destroy()
        
    return new, id_seq, path_file


def Read_Clustal(evt):
    global n
    global path_file
    global file
    global new
    global id_seq
    new = []
    Names = []
    id_seq = []
    path_file = []
    dialog = wx.FileDialog(Window1, message = 'Select .aln file',
                           defaultFile = '', wildcard = '*.aln',
                           style = wx.FD_OPEN, pos=(10,10))
    if dialog.ShowModal() == wx.ID_OK:
        selected_file = dialog.GetPaths()
        file = open(selected_file[0],'r')
        align = AlignIO.read(file, "clustal")
        sequences = []
        desc = [] #description
        for x in align:
            sequences.append(x.seq)
            id_seq.append(x.id)

            desc.append(x.description)
            
        new = (list(map(str,sequences)))
        path_file = os.path.basename(file.name)
        Names.append(path_file)
        List.InsertItems(Names,n)
        n=n+1
        List.Show()
        dialog.Destroy()

    return new, id_seq, path_file

###############################
def Plots(evt):

    global matrix
    #choosing the best path
    def get_crawl_matrix(matrix):
        x, y = len(matrix) - 1, len(matrix[0]) - 1
        result = matrix[x][y]
        while (x, y) != (0, 0):
            if x == 0:
                y -= 1
            elif y == 0:
                x -= 1
            else:
                vecs = [(-1, 0), (0, -1), (-1, -1)]
                poss = [(x + vec_x, y + vec_y) for vec_x, vec_y in vecs]
                poss = [(matrix[i][j], i, j) for i, j in poss]
                _, x, y = max(poss)
            result += matrix[x][y]
        return result

    #creating a scoring matrix for each column
    def get_matrix_score(A, B, penalty = -1):
        matrix = [[0 for e in B] for e in A]
        for i in range(len(A)):
            for j in range(len(B)):
                try:
                    matrix[i][j] = blosum62[(A[i], A[j])]
                except KeyError: #when it is - we give the point -1
                    matrix[i][j] = penalty
        return get_crawl_matrix(matrix)



    global get_data_matrix
    #we get the matrix - here we have to give sequences (depending on what we are loading)
    def get_data_matrix():
        global new
        global id_seq
        id_seq = id_seq
        new = new
        return new, id_seq


    def get_column(matrix, j):
        result = []
        for i in range(len(matrix)):
            result.append(matrix[i][j])
        return result


    def get_results(matrix):
        result = []
        for j in range(len(matrix[0])):
            word = get_column(matrix, j)
            word = ''.join(word)
            
            ga = get_matrix_score(word, word)
            result.append(ga)
        return result

    #create a plot
    def make_plot(data):
        plt.xlabel('column index')
        plt.ylabel('matrix score')
        x_seq = list(range(len(data)))
        y_seq = data
        ax = plt.gca()
        plt.ylim(min(y_seq)-5, max(y_seq)+5)
        ax.plot(x_seq, y_seq)
        
        #get updated x value after zoom
        def on_xlims_change(event_ax):
            global new_x
            new_x = event_ax.get_xlim()
            if new_x[0] < 0:
                new_x = list(new_x)
                new_x[0] = 0
                new_x = tuple(new_x)
                
            return new_x

        ax.callbacks.connect('xlim_changed', on_xlims_change)
        plt.show()

    def solve():
        filename = path_file
        matrix, _ = get_data_matrix()
        result = get_results(matrix)
        make_plot(result)

        
    if __name__ == '__main__':
        solve()

    return get_data_matrix(), new, id_seq

###############################
def Visualization_Of_Sequences(evt):

    protein = {
        '-': 'white',
        "A": "#F08080",
        "R": "#E9967A",
        "N": "#FFC0CB",
        "D": "#FFA07A",
        "C": "#FFEFD5",
        "Q": "#FFE4B5",
        "E": "#EEE8AA",
        "G": "#E6E6FA",
        "H": "#D8BFD8",
        "I": "#DDA0DD",
        "L": "#EE82EE",
        "K": "#AFEEEE",
        "M": "#5F9EA0",
        "F": "#ADD8E6",
        "P": "#A9A9A9",
        "S": "#DCDCDC",
        "T": "#6A5ACD",
        "W": "#FF7F50",
        "Y": "#FFDEAD",
        "V": "#DAA520",
        "B": "6495ED",
        "X": "00BFFF",
        "Z": "008000"
    }

    def get_cell(char):
        return f'<td style=\"background-color:{protein[char]};\">{char}</td>'


    def get_line(name, line):
        return (
            '<tr>\n\t\t' + 
            f'<td style=\"width:200px;\">{name}</td>' + 
            '\n\t\t'.join([get_cell(e) for e in line]) + 
            '\n\t</tr>'
        )


    def gen_table(matrix, id_seq, filename = 'file.html'):
        beg = '''<html>
    <head>
    <style>
    table, th, td {
      border: 1px solid black;
      border-collapse: collapse;
    }
    </style>
    </head>
    <body>'''
        en = '''</body>
    </html>'''

        res = (
            '<table style=\"border=\"3\";\">\n\t' + 
            '\n\t'.join([get_line(idd, line) for idd, line in zip(id_seq, matrix)]) + 
            '\n</table>'
        )
        open(filename, 'w+').write(beg + res + en)


    if __name__ == '__main__':
        gen_table(new, id_seq)

    os.startfile("file.html")

###############################
def Fragment_visualization(evt):
    
    from_x = int(new_x[0])
    to_x = int(new_x[1])
    global fragment
    
    fragment = []
    for i in range(1, len(new)):
        one = new[i]
        two = one[from_x:to_x]
        fragment.append(two)
  
    protein = {
        '-': 'white',
        "A": "#F08080",
        "R": "#E9967A",
        "N": "#FFC0CB",
        "D": "#FFA07A",
        "C": "#FFEFD5",
        "Q": "#FFE4B5",
        "E": "#EEE8AA",
        "G": "#E6E6FA",
        "H": "#D8BFD8",
        "I": "#DDA0DD",
        "L": "#EE82EE",
        "K": "#AFEEEE",
        "M": "#5F9EA0",
        "F": "#ADD8E6",
        "P": "#A9A9A9",
        "S": "#DCDCDC",
        "T": "#6A5ACD",
        "W": "#FF7F50",
        "Y": "#FFDEAD",
        "V": "#DAA520",
        "B": "6495ED",
        "X": "00BFFF",
        "Z": "008000"
    }

    def get_cell(char):
        return f'<td style=\"background-color:{protein[char]};\">{char}</td>'


    def get_line(name, line):
        return (
            '<tr>\n\t\t' + 
            f'<td style=\"width:200px;\">{name}</td>' + 
            '\n\t\t'.join([get_cell(e) for e in line]) + 
            '\n\t</tr>'
        )


    def gen_table(matrix, id_seq, filename = 'fragment_visualization.html'):
        beg = '''<html>
    <head>
    <style>
    table, th, td {
      border: 1px solid black;
      border-collapse: collapse;
    }
    </style>
    </head>
    <body>'''
        en = '''</body>
    </html>'''

        res = (
            '<table style=\"border=\"3\";\">\n\t' + 
            '\n\t'.join([get_line(idd, line) for idd, line in zip(id_seq, matrix)]) + 
            '\n</table>'
        )
        open(filename, 'w+').write(beg + res + en)


    if __name__ == '__main__':
        gen_table(fragment, id_seq)

    os.startfile("fragment_visualization.html")


### WX ###
Prog = wx.App()
Sequences = []
Names = []
n=0
Window1 = wx.Frame(None, title = 'Menu', size = (800,600), pos = (200,50))
MenuBar = wx.MenuBar()
ProgMenu = wx.Menu()
ProgMenuItem1 = ProgMenu.Append(wx.ID_ANY, 'Fasta', 'Read Data')
Window1.Bind(wx.EVT_MENU, Read_Fasta, ProgMenuItem1)
ProgMenuItem2 = ProgMenu.Append(wx.ID_ANY, 'Clustal', 'Read Data')
Window1.Bind(wx.EVT_MENU, Read_Clustal, ProgMenuItem2)
MenuBar.Append(ProgMenu, 'Data')
Window1.SetMenuBar(MenuBar)

ProgMenu = wx.Menu()

ProgMenuItem1 = ProgMenu.Append(wx.ID_ANY, 'Plot')
Window1.Bind(wx.EVT_MENU, Plots, ProgMenuItem1)
ProgMenuItem2 = ProgMenu.Append(wx.ID_ANY, 'Visualization of sequences')
Window1.Bind(wx.EVT_MENU, Visualization_Of_Sequences, ProgMenuItem2)
ProgMenuItem3 = ProgMenu.Append(wx.ID_ANY, 'Fragment visualization')
Window1.Bind(wx.EVT_MENU, Fragment_visualization, ProgMenuItem3)
MenuBar.Append(ProgMenu, 'Options')


ProgMenu = wx.Menu()
ProgMenuItem1 = ProgMenu.Append(wx.ID_EXIT, 'Close', 'Close the program')
MenuBar.Append(ProgMenu,'Exit')
Window1.Bind(wx.EVT_MENU, Exit, ProgMenuItem1)

panel = wx.Panel(parent = Window1)
List = wx.ListBox(parent = panel, pos = (20,20), size = (400,400), style=wx.LB_MULTIPLE,)
List.Hide()

Result = wx.ListBox(parent = panel, pos=(440,20), size=(200,400))
Result.Hide()

Window1.Show()

Prog.MainLoop()

