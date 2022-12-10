# filecoding=utf-8

import re
import math
import sitesclass


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False


Treal = sitesclass.Treal
Tcmplx = sitesclass.Tcmplx
Tsite = sitesclass.Tsite
Tsites = sitesclass.Tsites

headers = []
datatype = ''
timefactor = ''
dataunit = ''
coordorient = ''
geoorigin = ''
nPnS = ''


def readfiles(dir, sites):
    global datatype
    global timefactor
    global dataunit
    global coordorient
    global geoorigin
    global nPnS
    headers.clear()
    freqlist = []
    # 读取观测数据
    freqlist.clear()
    openfile = open(dir, 'r')
    lineindex = 0
    ImpUnitsFactor = 1.
    for line in openfile:
        items = re.split(r'[:|\s]+', line.strip())
        if items[0] != '>':
            if (items[0] == '#') or (not is_number(items[0])):
                headers.append(line.strip())
                continue
        else:
            lineindex = lineindex+1
            if lineindex == 1:
                datatype = items[1]
            elif lineindex == 2:
                timefactor = line.strip()[1:]
            elif lineindex == 3:
                dataunit = items[1]
                if items[1] == '[V/m]/[T]':
                    ImpUnitsFactor = 1./1000.
                elif items[1] == '[mV/km]/[nT]':
                    ImpUnitsFactor = 1.
                elif (items[1] == '[V/m]/[A/m]') or (items[1] == 'Ohm'):
                    ImpUnitsFactor = (4.*math.pi)/10000.
            elif lineindex == 4:
                coordorient = items[1]
            elif lineindex == 5:
                geoorigin = line.strip()[1:]
            elif lineindex == 6:
                nPnS = line.strip()[1:]
            continue

        index = sites.IndexofSite(items[1])
        if index < 0:
            site = Tsite()
            sites.items.append(site)

            site.id = items[1]
            site.lat = float(items[2])
            site.lon = float(items[3])
            site.x = float(items[4])
            site.y = float(items[5])
            site.elev = float(items[6])
            site.coordori = float(coordorient)
        else:
            site = sites.items[index]

        if items[0] not in freqlist:
            freqlist.append(items[0])
            site.freqs.append(items[0])
            site.hxori.append(0.0)
            site.hyori.append(0.0)
            site.exori.append(0.0)
            site.eyori.append(0.0)
            site.hxori_ref.append(0.0)
            site.hyori_ref.append(0.0)
            site.exori_ref.append(0.0)
            site.eyori_ref.append(0.0)
            site.zxx.append(Tcmplx())
            site.zxy.append(Tcmplx())
            site.zyx.append(Tcmplx())
            site.zyy.append(Tcmplx())
            site.tzx.append(Tcmplx())
            site.tzy.append(Tcmplx())
            site.rxx.append(Treal())
            site.pxx.append(Treal())
            site.rxy.append(Treal())
            site.pxy.append(Treal())
            site.ryx.append(Treal())
            site.pyx.append(Treal())
            site.ryy.append(Treal())
            site.pyy.append(Treal())
        if items[0] not in site.freqs:
            site.freqs.append(items[0])
            site.hxori.append(0.0)
            site.hyori.append(0.0)
            site.exori.append(0.0)
            site.eyori.append(0.0)
            site.hxori_ref.append(0.0)
            site.hyori_ref.append(0.0)
            site.exori_ref.append(0.0)
            site.eyori_ref.append(0.0)
            site.zxx.append(Tcmplx())
            site.zxy.append(Tcmplx())
            site.zyx.append(Tcmplx())
            site.zyy.append(Tcmplx())
            site.tzx.append(Tcmplx())
            site.tzy.append(Tcmplx())
            site.rxx.append(Treal())
            site.pxx.append(Treal())
            site.rxy.append(Treal())
            site.pxy.append(Treal())
            site.ryx.append(Treal())
            site.pyx.append(Treal())
            site.ryy.append(Treal())
            site.pyy.append(Treal())
        index = site.freqs.index(items[0])

        if items[7].lower() == 'zxx':
            if len(items) == 11:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 12:
                site.hxori[index] = float(items[11])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori[index] = site.hxori[index]
                site.eyori[index] = site.exori[index]+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 15:
                site.hxori[index] = float(items[11])
                site.hyori[index] = float(items[12])
                site.exori[index] = float(items[13])
                site.eyori[index] = float(items[14])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.zxx[index].re = float(items[8])
            site.zxx[index].im = float(items[9])
            site.zxx[index].err = float(items[10])
            site.zxx[index].s_re = True

            site.rxx[index].re = 0.2*(pow(float(items[8]), 2) +
                                      pow(float(items[9]), 2))*float(items[0])
            site.rxx[index].err = math.sqrt(
                0.4*site.rxx[index].re*pow(float(items[10]), 2)*float(items[0]))
            site.rxx[index].s_re = True

            site.pxx[index].re = 180/math.pi*math.atan2(float(items[9]), float(items[8]))
            site.pxx[index].err = 180/math.pi * \
                math.sqrt(8*pow(math.cos(site.pxx[index].re*math.pi/180), 4)*pow(
                    float(items[10]), 2)*(pow(float(items[8]), 2) + pow(float(items[9]), 2)) /
                    (pow(2*float(items[8]), 4)+pow(10, -8)))
            site.pxx[index].s_re = True

            if site.rxx[index].re > 0:
                sites.rhomax = max(sites.rhomax, site.rxx[index].re)
                sites.rhomin = min(sites.rhomin, site.rxx[index].re)
            sites.phasemax = max(sites.phasemax, site.pxx[index].re)
            sites.phasemin = min(sites.phasemin, site.pxx[index].re)

        elif items[7].lower() == 'zxy':
            if len(items) == 11:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 12:
                site.hxori[index] = float(items[11])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori[index] = site.hxori[index]
                site.eyori[index] = site.exori[index]+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 15:
                site.hxori[index] = float(items[11])
                site.hyori[index] = float(items[12])
                site.exori[index] = float(items[13])
                site.eyori[index] = float(items[14])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.zxy[index].re = float(items[8])
            site.zxy[index].im = float(items[9])
            site.zxy[index].err = float(items[10])
            site.zxy[index].s_re = True

            site.rxy[index].re = 0.2*(pow(float(items[8]), 2) +
                                      pow(float(items[9]), 2))*float(items[0])
            site.rxy[index].err = math.sqrt(
                0.4*site.rxy[index].re*pow(float(items[10]), 2)*float(items[0]))
            site.rxy[index].s_re = True

            site.pxy[index].re = 180/math.pi*math.atan2(float(items[9]), float(items[8]))
            site.pxy[index].err = 180/math.pi * \
                math.sqrt(8*pow(math.cos(site.pxy[index].re*math.pi/180), 4)*pow(
                    float(items[10]), 2)*(pow(float(items[8]), 2) + pow(float(items[9]), 2)) /
                    (pow(2*float(items[8]), 4)+pow(10, -8)))
            site.pxy[index].s_re = True

            if site.rxy[index].re > 0:
                sites.rhomax = max(sites.rhomax, site.rxy[index].re)
                sites.rhomin = min(sites.rhomin, site.rxy[index].re)
            sites.phasemax = max(sites.phasemax, site.pxy[index].re)
            sites.phasemin = min(sites.phasemin, site.pxy[index].re)

        elif items[7].lower() == 'zyx':
            if len(items) == 11:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 12:
                site.hxori[index] = float(items[11])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori[index] = site.hxori[index]
                site.eyori[index] = site.exori[index]+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 15:
                site.hxori[index] = float(items[11])
                site.hyori[index] = float(items[12])
                site.exori[index] = float(items[13])
                site.eyori[index] = float(items[14])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.zyx[index].re = float(items[8])
            site.zyx[index].im = float(items[9])
            site.zyx[index].err = float(items[10])
            site.zyx[index].s_re = True

            site.ryx[index].re = 0.2*(pow(float(items[8]), 2) +
                                      pow(float(items[9]), 2))*float(items[0])
            site.ryx[index].err = math.sqrt(
                0.4*site.ryx[index].re*pow(float(items[10]), 2)*float(items[0]))
            site.ryx[index].s_re = True

            site.pyx[index].re = 180/math.pi*math.atan2(float(items[9]), float(items[8]))
            site.pyx[index].err = 180/math.pi * \
                math.sqrt(8*pow(math.cos(site.pyx[index].re*math.pi/180), 4)*pow(
                    float(items[10]), 2)*(pow(float(items[8]), 2) + pow(float(items[9]), 2)) /
                    (pow(2*float(items[8]), 4)+pow(10, -8)))
            site.pyx[index].s_re = True

            if site.ryx[index].re > 0:
                sites.rhomax = max(sites.rhomax, site.ryx[index].re)
                sites.rhomin = min(sites.rhomin, site.ryx[index].re)
            sites.phasemax = max(sites.phasemax, site.pyx[index].re)
            sites.phasemin = min(sites.phasemin, site.pyx[index].re)

        elif items[7].lower() == 'zyy':
            if len(items) == 11:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 12:
                site.hxori[index] = float(items[11])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori[index] = site.hxori[index]
                site.eyori[index] = site.exori[index]+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 15:
                site.hxori[index] = float(items[11])
                site.hyori[index] = float(items[12])
                site.exori[index] = float(items[13])
                site.eyori[index] = float(items[14])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.zyy[index].re = float(items[8])
            site.zyy[index].im = float(items[9])
            site.zyy[index].err = float(items[10])
            site.zyy[index].s_re = True

            site.ryy[index].re = 0.2*(pow(float(items[8]), 2) +
                                      pow(float(items[9]), 2))*float(items[0])
            site.ryy[index].err = math.sqrt(
                0.4*site.ryy[index].re*pow(float(items[10]), 2)*float(items[0]))
            site.ryy[index].s_re = True

            site.pyy[index].re = 180/math.pi*math.atan2(float(items[9]), float(items[8]))
            site.pyy[index].err = 180/math.pi * \
                math.sqrt(8*pow(math.cos(site.pyy[index].re*math.pi/180), 4)*pow(
                    float(items[10]), 2)*(pow(float(items[8]), 2) + pow(float(items[9]), 2)) /
                    (pow(2*float(items[8]), 4)+pow(10, -8)))
            site.pyy[index].s_re = True

            if site.ryy[index].re > 0:
                sites.rhomax = max(sites.rhomax, site.ryy[index].re)
                sites.rhomin = min(sites.rhomin, site.ryy[index].re)
            sites.phasemax = max(sites.phasemax, site.pyy[index].re)
            sites.phasemin = min(sites.phasemin, site.pyy[index].re)

        elif items[7].lower() == 'tx':
            if len(items) == 11:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 12:
                site.hxori[index] = float(items[11])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori[index] = site.hxori[index]
                site.eyori[index] = site.exori[index]+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 15:
                site.hxori[index] = float(items[11])
                site.hyori[index] = float(items[12])
                site.exori[index] = float(items[13])
                site.eyori[index] = float(items[14])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.tzx[index].re = float(items[8])
            site.tzx[index].im = float(items[9])
            site.tzx[index].err = float(items[10])
            site.tzx[index].s_re = True

        elif items[7].lower() == 'ty':
            if len(items) == 11:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 12:
                site.hxori[index] = float(items[11])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori[index] = site.hxori[index]
                site.eyori[index] = site.exori[index]+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 15:
                site.hxori[index] = float(items[11])
                site.hyori[index] = float(items[12])
                site.exori[index] = float(items[13])
                site.eyori[index] = float(items[14])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.tzy[index].re = float(items[8])
            site.tzy[index].im = float(items[9])
            site.tzy[index].err = float(items[10])
            site.tzy[index].s_re = True

        items.clear()

    openfile.close()
    for site in sites.items:
        site.freqs = [1/float(i) for i in site.freqs]
    freqlist = [1/float(i) for i in freqlist]
    sites.freqmax = max(freqlist)
    sites.freqmin = min(freqlist)


def writefiles(dir, sites):
    openfile = open(dir, 'w')
    for line in headers:
        openfile.write(line + '\n')
    openfile.write('> ' + datatype + '\n')
    openfile.write('> ' + timefactor + '\n')
    openfile.write('> ' + dataunit + '\n')
    openfile.write('> ' + str(sites.items[0].coordori) + '\n')
    openfile.write('> ' + geoorigin + '\n')
    openfile.write('> ' + nPnS + '\n')
    for item in sites.items:
        for i in range(len(item.freqs)):
            if item.zxx[i].s_re:
                tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
                    1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
                tmpstr2 = '{} {:15.6e} {:15.6e} {:15.6e}'.format(
                    'ZXX', item.zxx[i].re, item.zxx[i].im, item.zxx[i].err)
                tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
                    item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
                tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
                openfile.write(tmpstr + '\n')
            if item.zxy[i].s_re:
                tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
                    1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
                tmpstr2 = '{} {:15.6e} {:15.6e} {:15.6e}'.format(
                    'ZXY', item.zxy[i].re, item.zxy[i].im, item.zxy[i].err)
                tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
                    item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
                tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
                openfile.write(tmpstr + '\n')
            if item.zyx[i].s_re:
                tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
                    1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
                tmpstr2 = '{} {:15.6e} {:15.6e} {:15.6e}'.format(
                    'ZYX', item.zyx[i].re, item.zyx[i].im, item.zyx[i].err)
                tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
                    item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
                tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
                openfile.write(tmpstr + '\n')
            if item.zyy[i].s_re:
                tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
                    1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
                tmpstr2 = '{} {:15.6e} {:15.6e} {:15.6e}'.format(
                    'ZYY', item.zyy[i].re, item.zyy[i].im, item.zyy[i].err)
                tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
                    item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
                tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
                openfile.write(tmpstr + '\n')
            # if item.tzx[i].s_re:
            #     tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
            #         1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
            #     tmpstr2 = '{} {:15.6e} {:15.6e} {:15.6e}'.format(
            #         'TX', item.tzx[i].re, item.tzx[i].im, item.tzx[i].err)
            #     tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
            #         item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
            #     tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
            #     openfile.write(tmpstr + '\n')
            # if item.tzy[i].s_re:
            #     tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
            #         1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
            #     tmpstr2 = '{} {:15.6e} {:15.6e} {:15.6e}'.format(
            #         'TY', item.tzy[i].re, item.tzy[i].im, item.tzy[i].err)
            #     tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
            #         item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
            #     tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
            #     openfile.write(tmpstr + '\n')
    openfile.close()


def UpdateErrFlr(ObseSites, RespSites):
    # 更新门槛误差
    for obs in ObseSites.items:
        for res in RespSites.items:
            if obs.id == res.id:
                for obsf in obs.freqs:
                    for resf in res.freqs:
                        if math.isclose(obsf, resf, rel_tol=0.00001):
                            io = obs.freqs.index(obsf)
                            ir = res.freqs.index(resf)
                            res.rxx[ir].err = obs.rxx[io].err
                            res.pxx[ir].err = obs.pxx[io].err
                            res.rxy[ir].err = obs.rxy[io].err
                            res.pxy[ir].err = obs.pxy[io].err
                            res.ryx[ir].err = obs.ryx[io].err
                            res.pyx[ir].err = obs.pyx[io].err
                            res.ryy[ir].err = obs.ryy[io].err
                            res.pyy[ir].err = obs.pyy[io].err
                            break
                break


def UpdateZErrFlr(ObseSites, RespSites):
    # 更新门槛误差
    for obs in ObseSites.items:
        for res in RespSites.items:
            if obs.id == res.id:
                for obsf in obs.freqs:
                    for resf in res.freqs:
                        if math.isclose(obsf, resf, rel_tol=0.00000001):
                            io = obs.freqs.index(obsf)
                            ir = res.freqs.index(resf)
                            res.zxx[ir].err = obs.zxx[io].err
                            res.zxy[ir].err = obs.zxy[io].err
                            res.zyx[ir].err = obs.zyx[io].err
                            res.zyy[ir].err = obs.zyy[io].err
                            break
                break
