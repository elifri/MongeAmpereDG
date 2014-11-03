#!/usr/bin/python

import sys
import re
import os

# Replacement Rules
# http://docs.python.org/2/library/re.html#re.sub

def locate(pattern = r'\d+[_]', root=os.curdir):
    for files in os.listdir(os.path.abspath(root)):
        for filename in re.findall(pattern, ''.join(files)):
            print(str(filename))
            yield filename

filepath = "."

no = 1
deg =2

for f_str in locate('MA_template.cfg',filepath):
    with open(filepath+"/"+f_str, "r") as f:
        print("opten "+filepath+"/"+f_str)
        inputstr = f.read()
        for p_grad in [0,10,30,50,70,100,200]:
            for alpha in [1, 3, 5, 7]:
                
                tempstr = re.sub("\<p\>", str(p_grad), inputstr)
                tempstr = re.sub("\<alpha\>", str(alpha), tempstr)
                tempstr = re.sub("\<no\>", str(no), tempstr)
                tempstr = re.sub("\<deg\>", str(deg), tempstr)

                # Write result to temp file
                outfile = "cases/MA"+str(no)+"_deg"+str(deg)+"_p"+str(p_grad)+"_alpha0"+str(alpha)+".cfg"
                output = open(outfile, "w+")
                output.write(tempstr)
                output.close()
