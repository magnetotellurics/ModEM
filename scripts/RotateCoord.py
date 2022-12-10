import sys
import math
import readdatafile
import readZdatafile


def RotateCoord(y, x, Angle):
    ynew = y*math.cos(Angle*math.pi/180.0)+x*math.sin(Angle*math.pi/180.0)
    xnew = x*math.cos(Angle*math.pi/180.0)-y*math.sin(Angle*math.pi/180.0)
    return ynew, xnew


def RotateCoord_Sites(sites, Angle):
    for item in sites.items:
        yout, xout = RotateCoord(item.y, item.x, Angle)
        item.y = yout
        item.x = xout
        item.coordori = item.coordori+Angle


if len(sys.argv) < 5:
    print('Input Parameters: [datatype] [datafile_Input] [datafile_Output] [rotation angle(deg)]')
    print('datatype: rp, z')
    print('No input parameters or parameters not enough.')
    exit()

datatype = sys.argv[1]
datafile_input = sys.argv[2]
datafile_output = sys.argv[3]
rotangle = float(sys.argv[4])

# Read ModEM datafile and output rotated datafile
if datatype == 'z':
    sites = readZdatafile.Tsites()
    readZdatafile.readfiles(datafile_input, sites)
    RotateCoord_Sites(sites, rotangle)
    readZdatafile.writefiles(datafile_output, sites)
elif datatype == 'rp':
    sites = readdatafile.Tsites()
    readdatafile.readfiles(datafile_input, sites)
    RotateCoord_Sites(sites, rotangle)
    readdatafile.writefiles(datafile_output, sites)
