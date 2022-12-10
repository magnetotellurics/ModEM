# filecoding=utf-8

import re
import copy
import math

D2R = math.pi/180
R2D = 180/math.pi


class Treal(object):
    def __init__(self):
        self.re = 0.0
        self.err = 0.0
        self.s_re = False


class Tcmplx(object):
    def __init__(self):
        self.re = 0.0
        self.im = 0.0
        self.err = 0.0
        self.s_re = False


class Tsite(object):

    def __init__(self):
        self.id = ''
        self.lon = 0.0
        self.lat = 0.0
        self.x = 0.0
        self.y = 0.0
        self.elev = 0.0
        self.coordori = 0.0
        self.freqs = []
        self.hxori = []
        self.hyori = []
        self.exori = []
        self.eyori = []
        self.hxori_ref = []
        self.hyori_ref = []
        self.exori_ref = []
        self.eyori_ref = []
        self.zxx = []
        self.zxy = []
        self.zyx = []
        self.zyy = []
        self.tzx = []
        self.tzy = []
        self.rxx = []
        self.pxx = []
        self.rxy = []
        self.pxy = []
        self.ryx = []
        self.pyx = []
        self.ryy = []
        self.pyy = []

    def Assign(self, source):
        self.id = source.id
        self.lon = source.lon
        self.lat = source.lat
        self.x = source.x
        self.y = source.y
        self.elev = source.elev
        self.coordori = source.coordori
        self.freqs = copy.deepcopy(source.freqs)
        self.hxori = copy.deepcopy(source.hxori)
        self.hyori = copy.deepcopy(source.hyori)
        self.exori = copy.deepcopy(source.exori)
        self.eyori = copy.deepcopy(source.eyori)
        self.hxori_ref = copy.deepcopy(source.hxori_ref)
        self.hyori_ref = copy.deepcopy(source.hyori_ref)
        self.exori_ref = copy.deepcopy(source.exori_ref)
        self.eyori_ref = copy.deepcopy(source.eyori_ref)
        self.zxx = copy.deepcopy(source.zxx)
        self.zxy = copy.deepcopy(source.zxy)
        self.zyx = copy.deepcopy(source.zyx)
        self.zyy = copy.deepcopy(source.zyy)
        self.tzx = copy.deepcopy(source.tzx)
        self.tzy = copy.deepcopy(source.tzy)
        self.rxx = copy.deepcopy(source.rxx)
        self.pxx = copy.deepcopy(source.pxx)
        self.rxy = copy.deepcopy(source.rxy)
        self.pxy = copy.deepcopy(source.pxy)
        self.ryx = copy.deepcopy(source.ryx)
        self.pyx = copy.deepcopy(source.pyx)
        self.ryy = copy.deepcopy(source.ryy)
        self.pyy = copy.deepcopy(source.pyy)

    def AddAngle(self, Angle):
        for i in range(len(self.freqs)):
            self.hxori[i] = self.hxori[i]+Angle
            self.hyori[i] = self.hyori[i]+Angle
            self.exori[i] = self.exori[i]+Angle
            self.eyori[i] = self.eyori[i]+Angle
            self.hxori_ref[i] = self.hxori_ref[i]+Angle
            self.hyori_ref[i] = self.hyori_ref[i]+Angle
            self.exori_ref[i] = self.exori_ref[i]+Angle
            self.eyori_ref[i] = self.eyori_ref[i]+Angle


class Tsites(object):

    def __init__(self):
        self.items = []
        self.rhomax = -1e20
        self.rhomin = 1e20
        self.phasemax = -1e20
        self.phasemin = 1e20
        self.freqmax = -1e20
        self.freqmin = 1e20
        self.latmin = 1e20
        self.latmax = -1e20
        self.lonmin = 1e20
        self.lonmax = -1e20
        self.xmin = 1e20
        self.xmax = -1e20
        self.ymin = 1e20
        self.ymax = -1e20

    def AddSite(self, site):
        self.items.append(site)

    def Clear(self):
        self.items.clear()
        self.rhomax = -1e20
        self.rhomin = 1e20
        self.phasemax = -1e20
        self.phasemin = 1e20
        self.freqmax = -1e20
        self.freqmin = 1e20
        self.latmin = 1e20
        self.latmax = -1e20
        self.lonmin = 1e20
        self.lonmax = -1e20
        self.xmin = 1e20
        self.xmax = -1e20
        self.ymin = 1e20
        self.ymax = -1e20

    def DelSite(self, index):
        self.items.remove(index)

    def IndexofSite(self, site):
        index = -1
        for item in self.items:
            if site == item.id:
                index = self.items.index(item)
                break
        return index

    def WriteFile(self, filename):
        f = open(filename, 'w')
        print('site freq  x  y  rxy pxy ryx pyx', file=f)
        for item in self.items:
            for i in range(len(item.freqs)):
                print(item.id, end=' ', file=f)
                print(item.freqs[i], end=' ', file=f)
                print(str(item.x), end=' ', file=f)
                print(str(item.y), end=' ', file=f)
                print(str(item.rxy[i].re), end=' ', file=f)
                print(str(item.pxy[i].re), end=' ', file=f)
                print(str(item.ryx[i].re), end=' ', file=f)
                print(str(item.pyx[i].re), file=f)

        f.close()
