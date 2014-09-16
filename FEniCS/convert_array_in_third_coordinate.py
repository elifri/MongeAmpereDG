#!/usr/bin/python

import sys
import numpy

from lxml import etree

print len(sys.argv)

if len(sys.argv) != 3:
  print 'Error, please specify an input and an outputfile!'
  sys.exit(-1)

filename_input = sys.argv[1]
filename_output = sys.argv[2]

print filename_input

tree = etree.parse(filename_input)

#get number of points
piece_data = tree.find('.//Piece')
number_of_points = int(piece_data.attrib['NumberOfPoints'])
print 'number of points', number_of_points

#extract solution values
for elem in tree.findall('.//DataArray'):
    if 'Name' in elem.attrib.keys():
      if elem.attrib['Name'] == 'u':
        sol_data = elem.text.split('  ')

x = numpy.zeros(number_of_points)
y = numpy.zeros(number_of_points)
z = numpy.zeros(number_of_points)

for points in tree.findall('.//Points'):
  for elem in points.findall('.//DataArray'):
    point_data = elem.text.split(' ')

    for i in range(0, number_of_points):
      x[i]=point_data[4*i]
      y[i]=point_data[4*i+1]
      z[i]=sol_data[i]
 
    #combine point data
    combi = []
    for e in zip(x,y,z):
      combi.extend(e)
    
    #write point data
    elem.text = ' '.join([str(number) for number in combi])

#write data to file
with open(filename_output, 'w') as file_handle:
    file_handle.write(etree.tostring(tree, pretty_print=True, encoding='utf8'))