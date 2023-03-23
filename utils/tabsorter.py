import os, sys

filename = sys.argv[1]
outname = 'test.py'

f = open(filename, 'r')
content = f.readlines()
f.close()

word = ''
for line in content:
    line = line.replace('    ','\t')
    word += line

g = open(filename,'w')
g.write(word)
g.close()
    