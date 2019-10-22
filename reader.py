# Designed to work with python2.
# reader.py
# Create ROOT TTree using the ODE integration data
#
import sys
import numpy
import ROOT
from math import sqrt, tan, atan2, atan, sin, cos, pi
import math
import random
import argparse

parser = argparse.ArgumentParser(description='Pendulum ODE tree maker')
parser.add_argument('-c','--case',type=int,default=2, help='case ')
args = parser.parse_args()
case = args.case
print 'case ',case

ROOT.gROOT.ProcessLine(
"struct MyStruct {\
   Int_t        fi;\
   Double_t     ft;\
   Double_t     ftheta;\
   Double_t     fomega;\
   Double_t     fenergy;\
};" );

mystruct = ROOT.MyStruct()

if case == 0:
   name = 'euler_0p01'
elif case == 1:
   name = 'rk2_0p01'
elif case == 2:
   name = 'rk4_0p01'

infile = name+'.dat'
outfile = name+'.root'

# Initialize histo file
f = ROOT.TFile(outfile, 'update')
tree = ROOT.TTree( 'T', 'Just A Tree' )
tree.Branch( 'i', ROOT.AddressOf( mystruct, 'fi' ), 'i/I' ) 
tree.Branch( 't', ROOT.AddressOf( mystruct, 'ft' ), 't/D' ) 
tree.Branch( 'theta', ROOT.AddressOf( mystruct, 'ftheta' ), 'theta/D' ) 
tree.Branch( 'omega', ROOT.AddressOf( mystruct, 'fomega' ), 'omega/D' )
tree.Branch( 'energy', ROOT.AddressOf( mystruct, 'fenergy' ), 'energy/D' )
htheta = ROOT.TH1D("htheta","theta vs time; time(s); theta(rad)",1001,-0.005,10.005)
homega = ROOT.TH1D("homega","omega vs time; time(s); omega(rad/s)",1001,-0.005,10.005)

#print('using file ',infile)
print 'reading as input file  ',infile
print 'writing as output file ',outfile
with open(infile) as fp:
   j = 0
   for line in fp:
      data = line.split()
      t = float(data[0])
      theta = float(data[1])
      omega = float(data[2])
      energy = float(data[3])
      mystruct.fi = j
      mystruct.ft = t
      mystruct.ftheta = theta
      mystruct.fomega = omega
      mystruct.fenergy = energy
      tree.Fill()
      htheta.Fill(t,theta)
      homega.Fill(t,omega)
      j += 1 

f.Write()
