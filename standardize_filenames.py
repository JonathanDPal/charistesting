from glob import glob
import os
import sys
from copy import deepcopy


def valuefinder(filename, param):
    param = str.lower(param)
    paramlengths = {'annuli': 6, 'subsections': 11, 'movement': 8, 'spectrum': 8, 'smooth': 6, 'highpass': 8, 'kl': 2}
    paramlength = paramlengths[param]
    startingindex = None  # will be defined soon
    for j in range(len(filename)):
        if str.lower(filename[j: j + paramlength]) == param:
            if param == 'kl':
                startingindex = j + paramlength
            else:
                startingindex = j - 1

    if startingindex is not None:
        if param != 'kl':
            valuelength = 0
            while startingindex >= 0 and filename[startingindex] != '_':
                startingindex -= 1
                valuelength += 1
            end_index = startingindex + 1 + valuelength
            value = filename[startingindex + 1: end_index]
        else:
            value = filename[startingindex: startingindex + 2]
    else:
        value, startingindex, end_index = None, None, None

    return value, startingindex, end_index


def standardize(filename):
    fix_mov, fix_cs, fix_hp = False, False, False
    movement, mov_si, mov_ei = valuefinder(filename, 'movement')
    smooth, cs_si, cs_ei = valuefinder(filename, 'smooth')
    highpass, hp_si, hp_ei = valuefinder(filename, 'highpass')
    if movement is not None and smooth is not None and highpass is not None:
        movement, smooth, highpass = str.lower(movement), str.lower(smooth), str.lower(highpass)
        if '.' not in movement:
            movement = float(int(movement))
            fix_mov = True
        if '.' not in smooth:
            smooth = float(int(smooth))
            fix_cs = True
        if '.' not in highpass:
            if highpass == 'true' or highpass == 'false':
                fix_hp = False
            else:
                highpass = float(int(highpass))
                fix_hp = True
        needs_fixing = fix_mov or fix_cs or fix_hp
        if needs_fixing:
            new_filename = deepcopy(filename)
            if fix_mov:
                new_filename = f'{new_filename[:mov_si]}{movement}{new_filename[mov_ei:]}'
                offset_mov = len(str(movement)) - (mov_ei - mov_si)
            if fix_cs:
                if fix_mov:
                    cs_si += offset_mov
                    cs_ei += offset_mov
                    new_filename = f'{new_filename[:cs_si]}{smooth}{new_filename[cs_ei:]}'
                else:
                    new_filename = f'{new_filename[:cs_si]}{smooth}{new_filename[cs_ei:]}'
                offset_cs = len(str(smooth)) - (cs_ei - cs_si)
            if fix_hp:
                if fix_mov and fix_cs:
                    hp_si += (offset_mov + offset_cs)
                    hp_ei += (offset_mov + offset_cs)
                    new_filename = f'{new_filename[:hp_si]}{highpass}{new_filename[hp_ei:]}'
                elif fix_mov:
                    hp_si += offset_mov
                    hp_ei += offset_mov
                    new_filename = f'{new_filename[:hp_si]}{highpass}{new_filename[hp_ei:]}'
                elif fix_cs:
                    hp_si += offset_cs
                    hp_ei += offset_cs
                    new_filename = f'{new_filename[:hp_si]}{highpass}{new_filename[hp_ei:]}'
                else:
                    new_filename = f'{new_filename[:hp_si]}{highpass}{new_filename[hp_ei:]}'
            os.rename(filename, new_filename)


direc = sys.argv[1]
stuff = glob(f'{direc}/*')
directories = []
filenames = []
for thing in stuff:
    if os.path.isdir(thing):
        directories.append(thing)
    elif os.path.isfile(thing):
        filenames.append(thing)

more_subdirectories = len(directories) != 0
directory_length = [0, len(directories)]
i = 0
while more_subdirectories:
    for directory in directories[directory_length[i]:]:
        for item in glob(f'{directory}/*'):
            if os.path.isdir(item):
                directories.append(item)
            elif os.path.isfile(item):
                filenames.append(item)
    directory_length.append(len(directories))
    i += 1
    more_subdirectories = not directory_length[i + 1] == directory_length[i]

for file in filenames:
    standardize(file)

# condor record = 3344 (209 KLIP runs simultaneously)
