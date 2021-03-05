#!/usr/bin/env python3
################################################################################
# Name:    vector_tracer.py
# Purpose: This Python script traces SVG and CSV vector files using the
#          discrete Fourier transform.
# Author:  Huidae Cho
# GitHub:  https://github.com/HuidaeCho/vector_tracer.py
# Since:   March 4, 2021
#
# Copyright (C) 2021, Huidae Cho <https://idea.isnew.info>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
################################################################################

import sys
import pathlib
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def fft_1d(f):
    '''
    1-dimensional fast Fourier transform
    f:  Original image
    '''
    M = len(f)
    K = int(M/2)
    feven = f[::2]
    fodd = f[1::2]
    if K == 1:
        Feven = feven[0]/2
        Fodd = fodd[0]/2
        F = np.array((Feven+Fodd, Feven-Fodd))
    else:
        Feven = fft_1d(feven)/2
        Fodd = fft_1d(fodd)/2*np.exp(-2j*np.pi*np.arange(0,K)/(2*K))
        F = np.concatenate((Feven+Fodd, Feven-Fodd))
    return F

def ifft_1d(F):
    '''
    1-dimensional inverse fast Fourier transform
    F:  DFT
    '''
    M = len(F)
    K = M/2
    Feven = F[::2]
    Fodd = F[1::2]
    if K == 1:
        feven = Feven[0]/2
        fodd = Fodd[0]/2
        f = np.array((feven+fodd, feven-fodd))
    else:
        feven = ifft_1d(Feven)/2
        fodd = ifft_1d(Fodd)/2*np.exp(2j*np.pi*np.arange(0,K)/(2*K))
        f = np.concatenate((feven+fodd, feven-fodd))
    return f

def pad_1d(x):
    '''
    Pad zeros to x so len(x) becomes a power of two
    x:  1-dimensional array
    '''
    m = len(x)
    M = 2**math.ceil(math.log2(m))
    if m != M:
        M -= m
        x = np.pad(x, (0, M), 'wrap')
    return x

def fft(f):
    '''
    2-dimensional fast Fourier transform
    f:  Original image
    '''
    M, N = f.shape
    Fx = np.zeros(f.shape, dtype=complex)
    for x in range(0, M):
        Fx[x] = fft_1d(f[x])
    F = np.zeros(f.shape, dtype=complex)
    for v in range(0, N):
        F[...,v] = fft_1d(Fx[...,v])
    return F

def ifft(F):
    '''
    2-dimensional inverse fast Fourier transform
    F:  DFT
    '''
    M, N = F.shape
    fx = np.zeros(F.shape, dtype=complex)
    for x in range(0, M):
        fx[x] = ifft_1d(F[x])
    f = np.zeros(F.shape, dtype=complex)
    for v in range(0, N):
        f[...,v] = ifft_1d(fx[...,v])
    return f

def pad(x):
    '''
    Pad zeros to x so x.shape becomes a power of two
    x:  2-dimensional array
    '''
    m, n = x.shape
    M = 2**math.ceil(math.log2(m))
    N = 2**math.ceil(math.log2(n))
    if m != M or n != N:
        M -= m
        N -= n
        x = np.pad(x, ((int(M/2), M-int(M/2)), (int(N/2), N-int(N/2))),
                   'constant')
    return x

def shift(x):
    '''
    Shift x by (-1)**(row+column)
    x:  2-dimensional array
    '''
    m, n = x.shape
    y = x.copy()
    for r in range(0, m):
        for c in range(0, n):
            y[r,c] *= (-1)**(r+c)

    return y

def spectrum(F):
    '''
    Calculate the Fourier spectrum
    F:  Fourier transform
    '''
    return abs(F)

def phase_angle(F):
    '''
    Calculate the Fourier phase angle
    F:  Fourier transform
    '''
    return np.arctan2(F.imag, F.real)

def power_spectrum(F):
    '''
    Calculate the Fourier power spectrum
    F:  Fourier transform
    '''
    return spectrum(F)**2

def plot(x, y, label):
    plt.plot(x, y, label=label)
    plt.legend()
    plt.show()

def pulse(x, params):
    K = params[0]
    A = params[1]
    return np.where(x < K, A, 0)

def line(x, params):
    slope = params[0]
    intercept = params[1]
    return slope*x+intercept

def trace_fourier(f, color='black', draw_circles=True, gif_path=None,
                  cuts=None, buf=0.2):
    # https://andymac-2.github.io/fourier-polygon/writeup/
    # https://github.com/andymac-2/fourier-polygon
    # https://dsp.stackexchange.com/a/59133/41245
    n = len(f)
    f = pad_1d(f)
    M = len(f)

    F = fft_1d(f)

    f_ui = []
    for u in range(M):
        f_ui.append([])
        for i in range(M):
            f_ui[u].append(F[u]*np.exp(2j*np.pi*u*i/M))

    x_min = min(f.real)
    x_max = max(f.real)
    y_min = min(f.imag)
    y_max = max(f.imag)
    x_range = x_max-x_min
    y_range = y_max-y_min
    xlim_min = x_min-buf*x_range
    xlim_max = x_max+buf*x_range
    ylim_min = y_min-buf*y_range
    ylim_max = y_max+buf*y_range

#    if cuts:
#        for i in range(len(cuts)-1):
#            s, e = cuts[i], cuts[i+1]
#            plt.plot(f[s:e].real, f[s:e].imag, 'g')
#    else:
#        plt.plot(f.real, f.imag, 'g')
#    plt.xlim(xlim_min, xlim_max)
#    plt.ylim(ylim_min, ylim_max)
#    plt.axis('square')
#    plt.show()

#    plt.plot(F.real, F.imag)
#    plt.axis('square')
#    plt.show()

#    plt.plot(X.real, 'g', label='X(u)')
#    plt.plot(Y.real, 'r--', label='Y(u)')
#    plt.title('Fourier Transform')
#    plt.xlabel('u')
#    plt.legend()
#    plt.show()

#    f2 = np.zeros(M, dtype=complex)
#    for i in range(M):
#        for u in np.arange(M):
#            f2[i] += F[u]*np.exp(2j*np.pi*u*i/M)
#    plt.plot(f2.real, f2.imag)
#    plt.show()

    fig = plt.figure()
    ax = plt.axes(xlim=(xlim_min, xlim_max), ylim=(ylim_min, ylim_max))
    ax.set_aspect('equal')
    ax.axis('off')
    f2 = []
    lines = []
    for u in range(2*M+1):
        if u < M:
            # circles
            line, = ax.plot([], [], ':', color='lightgray', linewidth=0.5)
        elif u < 2*M:
            # arms
            line, = ax.plot([], [], color='gray', linewidth=0.5)
        else:
            # figure
            line, = ax.plot([], [], color=color, linewidth=1)
        lines.append(line)

    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def animate(i):
        if i == 0:
            f2.clear()
        f_sum = 0
        f_prev = f_sum
        for u in range(M):
            if draw_circles:
                # circle
                f_u = np.zeros(M, dtype=complex)
                for j in range(M):
                    f_u[j] = f_prev+f_ui[u][j]
                lines[u].set_data(f_u.real, f_u.imag)

            # arm
            f_sum += f_ui[u][i]
            lines[M+u].set_data((f_prev.real, f_sum.real),
                                (f_prev.imag, f_sum.imag))
            f_prev = f_sum

        if cuts and len(f2) in cuts:
            f2.append(np.nan)
        else:
            f2.append(f_sum)
        x = np.array(f2).real
        y = np.array(f2).imag
        lines[2*M].set_data(x, y)
        return lines

    anim = FuncAnimation(fig, animate, init_func=init,
                         frames=n+1, interval=5, blit=True)

    if gif_path:
        anim.save(gif_path)
    else:
        plt.show()

def trace_svg(svg_path, color='black', draw_circles=True, gif_path=None):
    from xml.dom import minidom
    from svg.path import parse_path
    from svg.path.path import Line, Move, Close

    # https://stackoverflow.com/a/56913776
    doc = minidom.parse(svg_path)
    path_strings = []
    for path in doc.getElementsByTagName('path'):
        g = path.parentNode
        # XXX: for grasslogo.svg
        # these paths are part of GRASS GIS text with bezier curves, which are
        # not supported yet
        if g.getAttribute('id').startswith('flowRoot'):
            # in ('flowRoot5356', 'flowRoot6517'):
            continue
        # https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/transform
        transform = g.getAttribute('transform')
        m_matrix = re.match(
                'matrix\(([^,]*),([^,]*),([^,]*),([^,]*),([^,]*),([^,]*)\)',
                transform)
        m_translate = re.match('translate\(([^,]*),([^,]*)\)', transform)
        if m_matrix:
            m = m_matrix
            t = [float(m[1]), float(m[2]), float(m[3]),
                 float(m[4]), float(m[5]), float(m[6])]
        elif m_translate:
            m = m_translate
            t = [1, 0, 0, 1, float(m[1]), float(m[2])]
        else:
            t = [1, 0, 0, 1, 0, 0]
        path_strings.append({'transform':t, 'd':path.getAttribute('d')})
    doc.unlink()

    f = []
    cuts = []

    for i in range(len(path_strings)):
        t = path_strings[i]['transform']
        path_string = path_strings[i]['d']
        path = parse_path(path_string)
        first = True
        for e in path:
            if isinstance(e, Line):
                if first:
                    x = t[0]*e.start.real+t[2]*e.start.imag+t[4]
                    y = t[1]*e.start.real+t[3]*e.start.imag+t[5]
                    f.append(x-y*1j)
                    first = False
                x = t[0]*e.end.real+t[2]*e.end.imag+t[4]
                y = t[1]*e.end.real+t[3]*e.end.imag+t[5]
                f.append(x-y*1j)
            elif isinstance(e, Move):
                first = True
                cuts.append(len(f))
                # for cut in animation
                x = t[0]*e.start.real+t[2]*e.start.imag+t[4]
                y = t[1]*e.start.real+t[3]*e.start.imag+t[5]
                f.append(x-y*1j)
    cuts.append(len(f))

    trace_fourier(f, color, draw_circles, gif_path, cuts)

def trace_csv(csv_path, color='black', draw_circles=True, gif_path=None):
    import csv

    f = []
    with open(csv_path) as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            f.append(float(row[0])+float(row[1])*1j)

    trace_fourier(f, color, draw_circles, gif_path)

def trace_vector(input_path, color='black', draw_circles=True, gif_path=None):
    if input_path.lower().endswith('.svg'):
        trace_svg(input_path, color, draw_circles, gif_path)
    elif input_path.lower().endswith('.csv'):
        trace_csv(input_path, color, draw_circles, gif_path)
    else:
        print(f'{input_path}: Unsupported file format')
        exit(1)

def print_usage():
    print(f'''
Usage: vector_tracer.py [-h] [-c] input_path [color] [output_gif_path]

  -h               display this help and exit
  -c               draw circles

  input_path       SVG or CSV path
  color            outline color (default: black)
  output_gif_path  output GIF path (default: display)

  For color names, please check
  https://matplotlib.org/stable/gallery/color/named_colors.html
''')

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print_usage()
        exit(0)

    draw_circles = False
    input_path = None
    color = None
    gif_path = None

    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-h':
            print_usage()
            exit(0)
        elif sys.argv[i] == '-c':
            draw_circles = True
        elif not input_path:
            input_path = sys.argv[i]
            if not pathlib.Path(input_path).exists():
                print(f'{input_path}: No such file found')
                exit(1)
        elif not color:
            color = sys.argv[i]
        elif not gif_path:
            gif_path = sys.argv[i]
            if pathlib.Path(gif_path).exists():
                print(f'{gif_path}: File already exists')
                exit(1)
    if not color:
        color = 'black'

    trace_vector(input_path, color, draw_circles, gif_path)
