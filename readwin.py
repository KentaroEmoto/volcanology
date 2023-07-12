import numpy as np
import warnings
import glob
from obspy import Stream, Trace, UTCDateTime

#
# last update: 11 July 2023 
#


def __s4(value):
    return -(value & 0b1000) | (value & 0b0111)

def __read_win_1(filename):
    output = {}
    srates = {}

    with open(filename, "rb") as fpin:
        fpin.seek(0, 2)
        sz = fpin.tell()
        fpin.seek(0)
        leng = 0
        status0 = 0
        start = 0
        while leng < sz:
            truelen = int.from_bytes(fpin.read(4), 'big')
            leng = 4
            if truelen == 0:
                break
            buff = fpin.read(6)
            leng += 6

            yy = "20%02x" % ord(buff[0:1])
            mm = "%02x" % ord(buff[1:2])
            dd = "%02x" % ord(buff[2:3])
            hh = "%02x" % ord(buff[3:4])
            mi = "%02x" % ord(buff[4:5])
            sec = "%02x" % ord(buff[5:6])

            date = UTCDateTime(int(yy), int(mm), int(dd), int(hh), int(mi), int(sec))

            if start == 0:
                start = date
            if status0 == 0:
                sdata = None
            while leng < truelen:
                buff = fpin.read(4)
                leng += 4
                flag = '%02x' % ord(buff[0:1])
                chanum = '%02x' % ord(buff[1:2])
                chanum = "%02s%02s" % (flag, chanum)
                datawide = int('%x' % (ord(buff[2:3]) >> 4))
                srate = ord(buff[3:4])
                xlen = (srate - 1) * datawide
                if datawide == 0:
                    xlen = srate // 2
                    datawide = 0.5

                leng += 4
                idata22 = int.from_bytes(fpin.read(4), 'big', signed=True)

                if chanum in output:
                    output[chanum].append(idata22)
                else:
                    output[chanum] = [idata22, ]
                    srates[chanum] = srate
                sdata = fpin.read(xlen)
                leng += xlen

                if len(sdata) < xlen:
                    fpin.seek(-(xlen - len(sdata)), 1)
                    sdata += fpin.read(xlen - len(sdata))
                    msg = "This shouldn't happen, it's weird..."
                    warnings.warn(msg)

                if datawide == 0.5:
                    for i in range(xlen):
                        idata2 = output[chanum][-1] + (int.from_bytes(sdata[i:i + 1], 'big', signed=True) >> 4)
                        output[chanum].append(idata2)
                        if i == (xlen -1):
                            break
                        idata2 += __s4(int.from_bytes(sdata[i:i+1],'big') & 0b1111)
                        output[chanum].append(idata2)
                elif datawide == 1:
                    for i in range((xlen // datawide)):
                        idata2 = output[chanum][-1] + int.from_bytes(sdata[i:i + 1], 'big', signed=True)
                        output[chanum].append(idata2)
                elif datawide == 2:
                    for i in range((xlen // datawide)):
                        idata2 = output[chanum][-1] +\
                                int.from_bytes(sdata[i*2:(i+1)*2], 'big', signed=True)
                        output[chanum].append(idata2)
                elif datawide == 3:
                    for i in range((xlen // datawide)):
                        idata2 = output[chanum][-1] +\
                            int.from_bytes(sdata[i*3:(i+1)*3], 'big', signed=True)
                        output[chanum].append(idata2)
                elif datawide == 4:
                    for i in range((xlen // datawide)):
                        idata2 = output[chanum][-1] +\
                                int.from_bytes(sdata[i*4:(i+1)*4], 'big', signed=True)
                        output[chanum].append(idata2)
                else:
                    msg = "DATAWIDE is %s " % datawide + \
                        "but only values of 0.5, 1, 2, 3 or 4 are supported."
                    raise NotImplementedError(msg)
    traces = []
    for i in output.keys():
        t = Trace(data=np.array(output[i]))
        t.stats.channel = str(i)
        t.stats.sampling_rate = float(srates[i])
        t.stats.starttime = start
        traces.append(t)
    return Stream(traces=traces)

def read_win(*winfiles):
    winfile = []
    st = Stream()
    for path in winfiles:
        filenames = sorted(glob.glob(path))
        for filename in filenames:
            winfile.append(filename)
    for filename in winfile:
        st += __read_win_1(filename)
    st.merge(method=1)
    return st

