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
    headers.clear()
    global datatype
    global timefactor
    global dataunit
    global coordorient
    global geoorigin
    global nPnS
    freqlist = []
    # 读取观测数据
    freqlist.clear()
    openfile = open(dir, 'r')
    lineindex = 0
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

            sites.xmin = min(sites.xmin, site.x)
            sites.xmax = max(sites.xmax, site.x)
            sites.ymin = min(sites.ymin, site.y)
            sites.ymax = max(sites.ymax, site.y)
            sites.latmin = min(sites.latmin, site.lat)
            sites.latmax = max(sites.latmax, site.lat)
            sites.lonmin = min(sites.lonmin, site.lon)
            sites.lonmax = max(sites.lonmax, site.lon)
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
            site.rxx.append(Treal())
            site.pxx.append(Treal())
            site.rxy.append(Treal())
            site.pxy.append(Treal())
            site.ryx.append(Treal())
            site.pyx.append(Treal())
            site.ryy.append(Treal())
            site.pyy.append(Treal())
        index = site.freqs.index(items[0])

        if items[7].lower() == 'rhoxy':
            if len(items) == 10:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 11:
                site.hxori[index] = float(items[10])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori = site.hxori
                site.eyori = site.exori+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 14:
                site.hxori[index] = float(items[10])
                site.hyori[index] = float(items[11])
                site.exori[index] = float(items[12])
                site.eyori[index] = float(items[13])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.rxy[index].re = float(items[8])
            site.rxy[index].err = float(items[9])
            site.rxy[index].s_re = True
            sites.rhomax = max(sites.rhomax, site.rxy[index].re)
            sites.rhomin = min(sites.rhomin, site.rxy[index].re)
        elif items[7].lower() == 'phsxy':
            if len(items) == 10:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 11:
                site.hxori[index] = float(items[10])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori = site.hxori
                site.eyori = site.exori+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 14:
                site.hxori[index] = float(items[10])
                site.hyori[index] = float(items[11])
                site.exori[index] = float(items[12])
                site.eyori[index] = float(items[13])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.pxy[index].re = float(items[8])
            site.pxy[index].err = float(items[9])
            site.pxy[index].s_re = True
            sites.phasemax = max(sites.phasemax, site.pxy[index].re)
            sites.phasemin = min(sites.phasemin, site.pxy[index].re)
        elif items[7].lower() == 'rhoyx':
            if len(items) == 10:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 11:
                site.hxori[index] = float(items[10])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori = site.hxori
                site.eyori = site.exori+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 14:
                site.hxori[index] = float(items[10])
                site.hyori[index] = float(items[11])
                site.exori[index] = float(items[12])
                site.eyori[index] = float(items[13])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.ryx[index].re = float(items[8])
            site.ryx[index].err = float(items[9])
            site.ryx[index].s_re = True
            sites.rhomax = max(sites.rhomax, site.ryx[index].re)
            sites.rhomin = min(sites.rhomin, site.ryx[index].re)
        elif items[7].lower() == 'phsyx':
            if len(items) == 10:
                site.hxori[index] = 0.0
                site.hyori[index] = 90.0
                site.exori[index] = 0.0
                site.eyori[index] = 90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 11:
                site.hxori[index] = float(items[10])
                site.hyori[index] = site.hxori[index]+90.0
                site.exori = site.hxori
                site.eyori = site.exori+90.0
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            elif len(items) == 14:
                site.hxori[index] = float(items[10])
                site.hyori[index] = float(items[11])
                site.exori[index] = float(items[12])
                site.eyori[index] = float(items[13])
                site.hxori_ref[index] = site.hxori[index]
                site.hyori_ref[index] = site.hyori[index]
                site.exori_ref[index] = site.exori[index]
                site.exori_ref[index] = site.exori[index]
            else:
                raise RuntimeError('data block is not enough.')

            site.pyx[index].re = float(items[8])
            site.pyx[index].err = float(items[9])
            site.pyx[index].s_re = True
            sites.phasemax = max(sites.phasemax, site.pyx[index].re)
            sites.phasemin = min(sites.phasemin, site.pyx[index].re)

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
            if item.rxy[i].s_re:
                tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
                    1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
                tmpstr2 = '{} {:15.6e} {:15.6e}'.format('RHOXY', item.rxy[i].re,  item.rxy[i].err)
                tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
                    item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
                tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
                openfile.write(tmpstr + '\n')
            if item.pxy[i].s_re:
                tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
                    1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
                tmpstr2 = '{} {:15.6e} {:15.6e}'.format('PHSXY', item.pxy[i].re, item.pxy[i].err)
                tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
                    item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
                tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
                openfile.write(tmpstr + '\n')
            if item.ryx[i].s_re:
                tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
                    1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
                tmpstr2 = '{} {:15.6e} {:15.6e}'.format('RHOYX', item.ryx[i].re, item.ryx[i].err)
                tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
                    item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
                tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
                openfile.write(tmpstr + '\n')
            if item.pyx[i].s_re:
                tmpstr1 = '{:12.6e} {} {:.3f} {:.3f} {:15.3f} {:15.3f} {:15.3f}'.format(
                    1/item.freqs[i], item.id, item.lat, item.lon, item.x, item.y, item.elev)
                tmpstr2 = '{} {:15.6e} {:15.6e}'.format('PHSYX', item.pyx[i].re, item.pyx[i].err)
                tmpstr3 = '{:9.3f} {:9.3f} {:9.3f} {:9.3f}'.format(
                    item.hxori[i], item.hyori[i], item.exori[i], item.eyori[i])
                tmpstr = tmpstr1 + ' ' + tmpstr2 + ' ' + tmpstr3
                openfile.write(tmpstr + '\n')
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
