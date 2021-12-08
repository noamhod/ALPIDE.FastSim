import os, sys, time

processList      = ["elaser", "glaser"]
phaseList        = ["phase0", "phase1"]
xiList           = [10.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 0.5]
polarizationList = ["gpc", "ppw"]

for proc in processList:
    for phase in phaseList:
        for xiStr in xiList:
            for pol in polarizationList:
                if(proc=="glaser"):
                    additional=' -g'
                else:
                    additional=''
                command = "time python h5Format2root.py -process "+proc+" -phase "+phase+" -xi "+str(xiStr)+" -polarization "+pol+additional
                print(">>>>>> The command: ", command)
                os.system(command)