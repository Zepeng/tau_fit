import ROOT
#binning of fijs 15 bins equally from -1 to 1 in Cos Zenith Angle
# 15 bins equally from -0.1 to 1.1 in tau NN output.
def convert(skx, sample, error_name):
    tfile = ROOT.TFile('./sk' + str(skx) + '.' + sample + '.root','READ')
    fij_1D = tfile.Get(error_name + '_fij')
    fijs = []
    for i in range(15*15):
        fijs.append(fij_1D.GetBinContent(i+1))
    return fijs

def get_errors(skx, sample):
    tfile = ROOT.TFile('./sk' + str(skx) + '.' + sample + '.root','READ')
    keys = tfile.GetListOfKeys()
    ErrorNames = []
    for key in keys:
        if key.GetClassName() == 'TH1D' and key.GetName()[:-4] not in ErrorNames:
            ErrorNames.append(key.GetName()[:-4])
    tfile.Close()
    return ErrorNames

def get_max(skx, sample, error):
    pos_max = max(convert(skx, sample, error))
    neg_max = min(convert(skx, sample, error))
    return max(pos_max, abs(neg_max))

def write_root(skx, sample):
    errors = get_errors(skx, sample)
    print errors
    error_file = ROOT.TFile('./sys_pdf/error.sk' + str(skx) + '.' + sample + '.root','RECREATE')
    sk_type = ['I', 'II','III','IV']
    event_file = ROOT.TFile('./data/SK-' + sk_type[skx - 1] + '.root','READ')
    dict_sample = {'fcmc':'bkg','tau':'tau'}
    events_2D = event_file.Get(dict_sample[sample] + 'HistoZenith2D')
    th1 = ROOT.TH1F('th1', 'th1', 230, 0,230)
    for i in range(225):
        th1.SetBinContent(i+1,0)
    for error in errors:
        fijs = convert(skx, sample, error)
        if get_max(skx, sample, error) < 0.02 or 't2k' in error : #and 'pid_multi_ring' not in error and 'polfit_1_ring' not in error: #and get_max(skx, 'tau', error) < 0.2:
            continue
        for m in range(len(fijs)):
            th1.SetBinContent(m+1, fijs[m] + th1.GetBinContent(m+1))
        th2_pos = ROOT.TH2F(error + '_pos',error + '_pos', 15, -1.0, 1.0, 15, -0.1, 1.1 )
        th2_pos.SetDirectory(error_file)
        th2_neg = ROOT.TH2F(error + '_neg', error + '_neg', 15, -1.0, 1.0, 15, -0.1, 1.1)
        th2_neg.SetDirectory(error_file)
        error_file.cd()
        for i in range(15):
            for j in range(15):
                if fijs[15*i + j] > 0:
                    th2_pos.SetBinContent(i+1, j+1, fijs[15*i + j]*events_2D.GetBinContent(i+1, j+1) )
                else:
                    th2_neg.SetBinContent(i+1, j+1, -fijs[15*i + j]*events_2D.GetBinContent(i+1, j+1) )
        print error, th2_pos.GetMaximum()
        print error, th2_neg.GetMaximum()
        th2_pos.Write()
        th2_neg.Write()

    error_file.Close()

def mc_pdf(skx, sample):
    errors = get_errors(skx, sample)
    print errors
    error_file = ROOT.TFile('./mc_pdf/error.sk' + str(skx) + '.' + sample + '.root','RECREATE')
    sk_type = ['I', 'II','III','IV']
    event_file = ROOT.TFile('./data/SK-' + sk_type[skx - 1] + '.root','READ')
    dict_sample = {'fcmc':'bkg','tau':'tau'}
    events_2D = event_file.Get(dict_sample[sample] + 'HistoZenith2D')
    th1 = ROOT.TH1F('th1', 'th1', 230, 0,230)
    for i in range(225):
        th1.SetBinContent(i+1,0)
    for error in errors:
        #if error not in ['up_down_ratio', 'horizontal_vertical_ratio', 'K_pi_ratio', 'NC_CC_ratio', 'up_down_energy_calibration']:
        #    continue
        fijs = convert(skx, sample, error)
        if get_max(skx, sample, error) < 0.02:
            continue
        for m in range(len(fijs)):
            th1.SetBinContent(m+1, fijs[m] + th1.GetBinContent(m+1))
        th2_mc = ROOT.TH2F(error ,error , 15, -1.0, 1.0, 15, -0.1, 1.1 )
        th2_mc.SetDirectory(error_file)
        error_file.cd()
        for i in range(15):
            for j in range(15):
                th2_mc.SetBinContent(i+1, j+1, fijs[15*i + j]*events_2D.GetBinContent(i+1, j+1) )
        print error, th2_mc.GetMaximum()
        th2_mc.Write()

    error_file.Close()

if __name__ == '__main__':
    for i in range(1,5):
        write_root(i, 'fcmc')
        write_root(i,'tau')
