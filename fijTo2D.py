import ROOT
#binning of fijs 15 bins equally from -1 to 1 in Cos Zenith Angle
# 15 bins equally from -0.1 to 1.1 in tau NN output.
def convert(skx, sample, error_name):
    tfile = ROOT.TFile('./data/fij.sk' + str(skx) + '.' + sample + '.root','READ')
    fij_1D = tfile.Get(error_name + '_fij')
    fijs = []
    for i in range(15*15):
        fijs.append(fij_1D.GetBinContent(i+1))
    return fijs

def get_errors(skx, sample):
    tfile = ROOT.TFile('./data/fij.sk' + str(skx) + '.' + sample + '.root','READ')
    keys = tfile.GetListOfKeys()
    ErrorNames = []
    for key in keys:
        if key.GetClassName() == 'TH1D' and key.GetName()[:-4] not in ErrorNames:
            ErrorNames.append(key.GetName()[:-4])
    tfile.Close()
    return ErrorNames

def write_root(skx, sample):
    errors = get_errors(skx, sample)
    print errors
    error_file = ROOT.TFile('./sys_pdf/error.sk' + str(skx) + '.' + sample + '.root','RECREATE')
    sk_type = ['I', 'II','III','IV']
    event_file = ROOT.TFile('./objects_extra/SK-' + sk_type[skx - 1] + '.root','READ')
    dict_sample = {'fcmc':'bkg','tau':'tau'}
    events_2D = event_file.Get(dict_sample[sample] + 'HistoZenith2D')
    th1 = ROOT.TH1F('th1', 'th1', 230, 0,230)
    for i in range(225):
        th1.SetBinContent(i+1,0)
    for error in errors:
        fijs = convert(skx, sample, error)
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

if __name__ == '__main__':
    for i in range(4,5):
        write_root(i, 'fcmc')
        write_root(i,'tau')
