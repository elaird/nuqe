#!/usr/bin/python

import os,sys
filename=sys.argv[1]
first_line_to_save=int(sys.argv[2])
last_line_to_save=int(sys.argv[3])
#cmd="cat "+file+" | head -"+str(last_line_to_save)+" | tail -"+str(last_line_to_save-first_line_to_save+1)+" | grep -v S | grep -v N | grep -v set > out.txt"

file=open(filename)
line_index=0
graph_index=0
print "graph0=TGraph()"

keep_going=True
for line in file :
    line_index=line_index+1
    if (keep_going) :
        if (line_index<first_line_to_save) : continue
        if (line_index>last_line_to_save) : continue
        if (line.find("S")>-1) : continue
        if (line.find("N")>-1) : continue
        if (line.find("set")>-1) : continue
        if (line.find("showpage")>-1) : 
            keep_going=False
            continue

        list=line.split()
        print "graph0.SetPoint("+str(graph_index)+","+list[0]+","+list[1]+")"
        graph_index=graph_index+1

print "graph0.Draw(\"ap\")"
