
import numpy as np
import matplotlib.lines as mlines
import scipy
import copy
import skrf as rf
class BWandFcCalculations():
    def __init__(self):
        print('')

    def calculateBWat3(self, net, threshold):
        net.interpolate_self(50 * net.f.size)
        left =  range(len(net.f[:np.argmax(net.s_mag[:, 1, 0])]))
        right = range(left.__len__(), net.frequency.__len__())
        maxima= abs(net.s_mag[:, 1, 0].max())
        for i in range(net.frequency.__len__()):
            net.s[i, 1, 0] = net.s[i,1,0]/maxima
        lowerFreqthreshold = 1
        upperFrequencyThreshold =1
        threshold = 10**(threshold/20)
        diff = 1-threshold
        left = left.__reversed__()                                          # left contains indexes in descending order of all values at frequencies lower than that of maxima.
        for i in left:
            if(abs(net.s[i,1,0])-threshold)<0 and diff > 0:
                #lowerFreqthreshold =(net.f[i] +net.f[i+1])/2
                if (threshold-abs(net.s[i,1,0])) > diff:                    # This way, the difference will not be higher than half of frequency step
                    lowerFreqthreshold = net.f[i+1]-(net.f[i]-net.f[i-1])/4 # Adding quarter frequency step will give worst/best case difference difference of a quarter freqeuncy step.
                else:                                                       # Adding half a frequency step will give worst case difference of half freqeuncy step and will be useless
                    lowerFreqthreshold = net.f[i]+(net.f[i]-net.f[i-1])/4
                break
            diff = abs(net.s[i, 1, 0]) - threshold
        diff = 1-threshold
        for i in right:
            if (abs(net.s[i, 1, 0]) - threshold) < 0 and diff > 0:
                #upperFrequencyThreshold =(net.f[i] +net.f[i-1])/2
                if (threshold-abs(net.s[i,1,0])) > diff:
                    upperFrequencyThreshold = net.f[i-1]+(net.f[i]-net.f[i-1])/4
                else:
                    upperFrequencyThreshold = net.f[i]-(net.f[i]-net.f[i-1])/4
                break
            diff = abs(net.s[i, 1, 0]) - threshold
        net.interpolate_self((net.f.size / 50).__int__())
        for i in range(net.frequency.__len__()):
            net.s[i, 1, 0] = net.s[i,1,0]*maxima
        return (upperFrequencyThreshold-lowerFreqthreshold)/1e6, (upperFrequencyThreshold+lowerFreqthreshold)/2e6, 20*np.log10(maxima)

    def calculate_RightLeft(self,net, threshold):
        net.interpolate_self(50 * net.f.size)
        left = range(len(net.f[:np.argmax(net.s_mag[:, 1, 0])]))
        right = range(left.__len__(), net.frequency.__len__())
        maxima = abs(net.s_mag[:, 1, 0].max())
        for i in range(net.frequency.__len__()):
            net.s[i, 1, 0] = net.s[i, 1, 0] / maxima
        lowerFreqthreshold = 1
        upperFrequencyThreshold = 1
        threshold = 10 ** (threshold / 20)
        diff = 1 - threshold
        left = left.__reversed__()                              # left contains indexes in descending order of all values at frequencies lower than that of maxima.
        for i in left:
            if (abs(net.s[i, 1, 0]) - threshold) < 0 and diff > 0:
                # lowerFreqthreshold =(net.f[i] +net.f[i+1])/2
                if (threshold - abs(net.s[i, 1, 0])) > diff:    # This way, the difference will not be higher than half of frequency step
                    lowerFreqthreshold = net.f[i + 1] - (net.f[i] - net.f[
                        i - 1]) / 4                             # Adding quarter frequency step will give worst/best case difference difference of a quarter freqeuncy step.
                else:                                           # Adding half a frequency step will give worst case difference of half freqeuncy step and will be useless
                    lowerFreqthreshold = net.f[i] + (net.f[i] - net.f[i - 1]) / 4
                break
            diff = abs(net.s[i, 1, 0]) - threshold
        diff = 1 - threshold
        for i in right:
            if (abs(net.s[i, 1, 0]) - threshold) < 0 and diff > 0:
                if (threshold - abs(net.s[i, 1, 0])) > diff:
                    upperFrequencyThreshold = net.f[i - 1] + (net.f[i] - net.f[i - 1]) / 4
                else:
                    upperFrequencyThreshold = net.f[i] - (net.f[i] - net.f[i - 1]) / 4
                break
            diff = abs(net.s[i, 1, 0]) - threshold
        net.interpolate_self((net.f.size / 50).__int__())
        for i in range(net.frequency.__len__()):
            net.s[i, 1, 0] = net.s[i, 1, 0] * maxima
        return lowerFreqthreshold/1e6, upperFrequencyThreshold/1e6

    def calculate_RightLeft_abs(self,net, threshold):
        net.interpolate_self(50 * net.f.size)
        maxima = abs(net.s_mag[:, 1, 0].max())
        left = range(len(net.f[:np.argmax(net.s_mag[:, 1, 0])]))
        right = range(left.__len__(), net.frequency.__len__())
        lowerFreqthreshold = 1
        upperFrequencyThreshold = 1
        threshold = 10 ** (threshold / 20)
        diff = maxima-threshold
        if diff <= 0:
            return 0, 0
        left = left.__reversed__()                              # left contains indexes in descending order of all values at frequencies lower than that of maxima.
        for i in left:
            if (abs(net.s[i, 1, 0]) - threshold) < 0 and diff > 0:
                if (threshold - abs(net.s[i, 1, 0])) > diff:    # This way, the difference will not be higher than half of frequency step
                    lowerFreqthreshold = net.f[i + 1] - (net.f[i] - net.f[
                        i - 1]) / 4                             # Adding quarter frequency step will give worst/best case difference difference of a quarter freqeuncy step.
                else:                                           # Adding half a frequency step will give worst case difference of half freqeuncy step and will be useless
                    lowerFreqthreshold = net.f[i] + (net.f[i] - net.f[i - 1]) / 4
                break
            diff = abs(net.s[i, 1, 0]) - threshold
        diff = maxima - threshold
        for i in right:
            if (abs(net.s[i, 1, 0]) - threshold) < 0 and diff > 0:
                if (threshold - abs(net.s[i, 1, 0])) > diff:
                    upperFrequencyThreshold = net.f[i - 1] + (net.f[i] - net.f[i - 1]) / 4
                else:
                    upperFrequencyThreshold = net.f[i] - (net.f[i] - net.f[i - 1]) / 4
                break
            diff = abs(net.s[i, 1, 0]) - threshold
        net.interpolate_self((net.f.size / 50).__int__())
        return lowerFreqthreshold/1e6, upperFrequencyThreshold/1e6

    def calculateBWat3_abs(self, net, threshold):
        net.interpolate_self(50 * net.f.size)
        left =  range(len(net.f[:np.argmax(net.s_mag[:, 1, 0])]))
        right = range(left.__len__(), net.frequency.__len__())
        maxima= abs(net.s_mag[:, 1, 0].max())
        lowerFreqthreshold = 1
        upperFrequencyThreshold =1
        threshold = 10**(threshold/20)
        diff = maxima-threshold
        if diff <=0:
            return 0,0,0
        left = left.__reversed__()                                          # left contains indexes in descending order of all values at frequencies lower than that of maxima.
        for i in left:
            if(abs(net.s[i,1,0])-threshold)<0 and diff > 0:
                if (threshold-abs(net.s[i,1,0])) > diff:                    # This way, the difference will not be higher than half of frequency step
                    lowerFreqthreshold = net.f[i+1]-(net.f[i]-net.f[i-1])/4 # Adding quarter frequency step will give worst/best case difference difference of a quarter freqeuncy step.
                else:                                                       # Adding half a frequency step will give worst case difference of half freqeuncy step and will be useless
                    lowerFreqthreshold = net.f[i]+(net.f[i]-net.f[i-1])/4
                break
            diff = abs(net.s[i, 1, 0]) - threshold
        diff = maxima-threshold
        for i in right:
            if (abs(net.s[i, 1, 0]) - threshold) < 0 and diff > 0:
                if (threshold-abs(net.s[i,1,0])) > diff:
                    upperFrequencyThreshold = net.f[i-1]+(net.f[i]-net.f[i-1])/4
                else:
                    upperFrequencyThreshold = net.f[i]-(net.f[i]-net.f[i-1])/4
                break
            diff = abs(net.s[i, 1, 0]) - threshold
        net.interpolate_self((net.f.size / 50).__int__())
        return (upperFrequencyThreshold-lowerFreqthreshold)/1e6, (upperFrequencyThreshold+lowerFreqthreshold)/2e6, 20*np.log10(maxima)

    def calculateMaxRLinPB(self, net, file):
        mat = scipy.io.loadmat(file)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpr = mat.get('fpr') * 1e6
        fpr2 = copy.deepcopy(fpr)
        fpr2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpr2[0][0] -= f_c * TCf * (OTR_low - Ta)
        freq_reange = (fpr2[0][0] / 1000000+10).__str__() + '-' + (fpr2[0][1] / 1000000-10).__str__() + 'MHz'
        s11_max = abs(net[freq_reange].s_mag[:, 0, 0].max())
        s22_max = abs(net[freq_reange].s_mag[:, 1, 1].max())
        return 20 * np.log10(s11_max) , 20 * np.log10(s22_max)

    def calculateS11Width(self,net, threshold):
        RL = net.s_db[:,0,0]
        freq = net.f[np.argwhere(RL<threshold)]
        freq = np.array(freq).reshape(-1)
        if (len(freq)==0):
            return 0
        elif all(np.diff(freq)==np.diff(freq)[0]):
            freqRange = freq[len(freq)-1]-freq[0]
            return freqRange/1000000
        else:
            step = freq[1]-freq[0]
            for i in range(len(freq)):
                if (freq[i+1]-freq[i]!= step):
                    return (freq[i]-freq[0])/1000000

    def CheckPBCompliance_withoutTemp(self, net, sptfile):
        mat = scipy.io.loadmat(sptfile)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpt = mat.get('fpt') * 1e6  # insertion loss limits
        mpt = mat.get('mpt')
        fpt2 = copy.deepcopy(fpt)
        fpt2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpt2[0][0] -= f_c * TCf * (OTR_low - Ta)
        net = net.interpolate(rf.Frequency(fpt[0][0],fpt[0][1],net[str(fpt[0][0])+'-'+str(fpt[0][1])+'Hz'].frequency.npoints))
        return (all(i >= mpt[0][0] for i in net.s_db[:,0,1]))

    def CheckPBCompliance(self, net, sptfile):
        mat = scipy.io.loadmat(sptfile)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpt = mat.get('fpt') * 1e6  # insertion loss limits
        mpt = mat.get('mpt')
        fpt2 = copy.deepcopy(fpt)
        fpt2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpt2[0][0] -= f_c * TCf * (OTR_low - Ta)
        net = net.interpolate(rf.Frequency(fpt2[0][0],fpt2[0][1],net[str(fpt2[0][0])+'-'+str(fpt2[0][1])+'Hz'].frequency.npoints))
        return (all(i >= mpt[0][0] for i in net.s_db[:,0,1]))

    def checkS11compliance_withoutTemp(self, net, sptfile):
        mat = scipy.io.loadmat(sptfile)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpr = mat.get('fpr') * 1e6  # insertion loss limits
        mpr = mat.get('mpr')
        fpr2 = copy.deepcopy(fpr)
        fpr2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpr2[0][0] -= f_c * TCf * (OTR_low - Ta)
        net = net.interpolate(rf.Frequency(fpr[0][0], fpr[0][1], net[str(fpr[0][0]) + '-' + str(fpr[0][1]) + 'Hz'].frequency.npoints))
        return (all(i <= mpr[0][0] for i in net.s_db[:, 0, 0]))


    def checkS11compliance(self, net, sptfile):
        mat = scipy.io.loadmat(sptfile)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpr = mat.get('fpr') * 1e6  # insertion loss limits
        mpr = mat.get('mpr')
        fpr2 = copy.deepcopy(fpr)
        fpr2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpr2[0][0] -= f_c * TCf * (OTR_low - Ta)
        net = net.interpolate(rf.Frequency(fpr2[0][0], fpr2[0][1], net[str(fpr2[0][0]) + '-' + str(fpr2[0][1]) + 'Hz'].frequency.npoints))
        return (all(i <= mpr[0][0] for i in net.s_db[:, 0, 0]))

    def checkS22compliance_withoutTemp(self, net, sptfile):
        mat = scipy.io.loadmat(sptfile)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpr = mat.get('fpr') * 1e6  # insertion loss limits
        mpr = mat.get('mpr')
        fpr2 = copy.deepcopy(fpr)
        fpr2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpr2[0][0] -= f_c * TCf * (OTR_low - Ta)
        net = net.interpolate(rf.Frequency(fpr[0][0], fpr[0][1], net[str(fpr[0][0]) + '-' + str(fpr[0][1]) + 'Hz'].frequency.npoints))
        return (all(i <= mpr[0][0] for i in net.s_db[:, 1, 1]))


    def checkS22compliance(self, net, sptfile):
        mat = scipy.io.loadmat(sptfile)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpr = mat.get('fpr') * 1e6  # insertion loss limits
        mpr = mat.get('mpr')
        fpr2 = copy.deepcopy(fpr)
        fpr2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpr2[0][0] -= f_c * TCf * (OTR_low - Ta)
        net = net.interpolate(rf.Frequency(fpr2[0][0], fpr2[0][1], net[str(fpr2[0][0]) + '-' + str(fpr2[0][1]) + 'Hz'].frequency.npoints))
        return (all(i <= mpr[0][0] for i in net.s_db[:, 1, 1]))



    def calculatePassbandWidth(self, net, threshold):
        RL = net.s_db[:, 0, 1]
        freq = net.f[np.argwhere(RL > threshold)]
        freq = np.array(freq).reshape(-1)
        if (len(freq) == 0):
            return 0
        elif all(np.diff(freq) == np.diff(freq)[0]):
            freqRange = freq[len(freq) - 1] - freq[0]
            return freqRange / 1000000
        else:
            step = freq[1] - freq[0]
            for i in range(len(freq)):
                if (freq[i + 1] - freq[i] != step):
                    return (freq[i] - freq[0]) / 1000000

    def add_RL_limit(self, ax, file):
        mat = scipy.io.loadmat(file)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpr = mat.get('fpr') * 1e6  # return loss limit for S11/s22
        mpr = mat.get('mpr')
        fpr2 = copy.deepcopy(fpr)
        fpr2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpr2[0][0] -= f_c * TCf * (OTR_low - Ta)
        line = mlines.Line2D(fpr2, mpr)
        line.set_color('green')
        line.set_linestyle('--')
        ax.add_line(line)
        line = mlines.Line2D(fpr, mpr)
        line.set_color('red')
        ax.add_line(line)

    def add_RL_limit_withoutTemperature(self, ax, file):
        mat = scipy.io.loadmat(file)
        f_c = mat.get('f_c')[0] * 1e6
        fpr = mat.get('fpr') * 1e6  # return loss limit for S11/s22
        mpr = mat.get('mpr')
        line = mlines.Line2D(fpr, mpr)
        line.set_color('red')
        ax.add_line(line)


    def add_s21_limit(self,ax,file):
        mat = scipy.io.loadmat(file)
        tkh = mat.get('tkh')[0]
        TCf = tkh[0]
        Ta = tkh[4]
        OTR_high = tkh[3]
        OTR_low = tkh[2]
        f_c = mat.get('f_c')[0] * 1e6
        fpt = mat.get('fpt') * 1e6  # insertion loss limits
        mpt = mat.get('mpt')
        fpt2 = copy.deepcopy(fpt)
        fpt2[0][1] -= f_c * TCf * (OTR_high - Ta)
        fpt2[0][0] -= f_c * TCf * (OTR_low - Ta)

        fst = mat.get('fst') * 1e6  # Transfer limits S21
        mst = mat.get('mst')

        fst2 = copy.deepcopy(fst)
        for j in range(len(fst2[0])):
            if (fst[0][j] == fst[0][j - 1]):
                fst2[0][j] = fst2[0][j - 1]
            elif (fst2[0][j] > f_c or fst2[0][j] < 1220000000):
                fst2[0][j] = fst2[0][j] - f_c * TCf * (OTR_low - Ta)
            else:
                fst2[0][j] = fst2[0][j] - f_c * TCf * (OTR_high - Ta)

                # else:
                #     fst2[0][j] = fst2[0][j] - f_c * TCf * (OTR_low - Ta)
                #     print(mst[0][j])
        line = mlines.Line2D(fpt2, mpt)
        line.set_color('green')
        line.set_linestyle('--')
        ax.add_line(line)
        line = mlines.Line2D(fpt, mpt)
        line.set_color('red')
        ax.add_line(line)
        j = 0

        line2 = mlines.Line2D(fst2[j:j + 1], mst[j:j + 1])
        line2.set_color('green')
        line2.set_linestyle('--')
        ax.add_line(line2)

        i = 0
        line = mlines.Line2D(fst[i:i + 1], mst[i:i + 1])
        line.set_color('red')
        ax.add_line(line)

    def add_s21_limit_withoutTemperature(self, ax, file):
        mat = scipy.io.loadmat(file)
        fpt = mat.get('fpt') * 1e6  # insertion loss limits
        mpt = mat.get('mpt')
        fst = mat.get('fst') * 1e6  # Transfer limits S21
        mst = mat.get('mst')
        line = mlines.Line2D(fpt, mpt)
        line.set_color('red')
        ax.add_line(line)
        i = 0
        line = mlines.Line2D(fst[i:i + 1], mst[i:i + 1])
        line.set_color('red')
        ax.add_line(line)
