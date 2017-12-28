import os
import glob
import time
import sys

def main():
  #WhatToCount=["Total Events","TotalNuInteractions","TotalNuInteractionsinVertex","CCNuInteractions","CC-COH_Only pi- in TPC","CC-COH_Only mu+ in TPC","CC-COH_Only mu- in TPC","CC-COH_Only pi+ in TPC","CC-COH_Both  Particles Out of TPC","CC-COH_Both  Particles in TPC","CC-COH_13","CC-COH_211","CCQE_Only mu- in TPC","CCQE_Only p+ in TPC","CCQE_Both  Particles Out of TPC","CCQE_Both  Particles in TPC","CCQE_13","CCQE_2212","CCQE_ProCount","CCQRes_Only mu- in TPC","CCQRes_Only Pion+ in TPC","CCQRes_Both  Particles Out of TPC","CCQRes_13","CCQRes_211","CCQECount","CCResCount_","CC-COH"]
  #Count=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  Files=[]
  counter=0
  what={}
 

  for file in glob.glob(r'*/*.out'):
    Files.insert(counter,os.path.abspath(file))
    counter+=1

  for i in range(len(Files)):
    readFile=open(Files[i],"r")
    readLine=readFile.readlines()
    status(i,counter)    
    for k in (readLine):
      x=k.split("<>")
      if (x[0].strip() in what and len(x)>1 and x[1].strip().isdigit()) :
        what[x[0].strip()]+=int(x[1].strip())
      elif (x[0].strip() not in what and len(x)>1 and x[1].strip().isdigit()) :
	      what[x[0].strip()]=int(x[1].strip())  
    readFile.close()
    
  file1=raw_input("please name the file->  ");
  WriteFile=open("/uboone/app/users/ilker228/v06_26_00/Grid/" + file1,"w+")
  
  for i in sorted(what):
    WriteFile.write(str(i) + " = " + str(what[i]) + "\n")
  WriteFile.close()
 

def status(r,n):
    sys.stdout.write('\r')
    percent=r/float(n)*100
    sys.stdout.write("Process ->>>>>  "+ str(round(percent)) + "%")
    sys.stdout.flush()



if __name__=="__main__":
  print 'Starting to Collect the information..',
  print '\b'*12
  main()
  print '\n Completed!' 
  print('\a')
  print('\a')
  print('\a')
  print('\a')

 

