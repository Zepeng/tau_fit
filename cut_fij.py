import ROOT, os
from fijTo2D import get_errors, convert

def get_max(skx, sample):
    max_ratio = []
    error_names = get_errors(skx, sample)
    for error in error_names:
        pos_max = max(convert(skx, sample, error))
        neg_max = min(convert(skx, sample, error))
        max_ratio.append(max(pos_max, abs(neg_max)))
    return max_ratio

def get_sum(skx, sample):
    sum_ratio = []
    error_names = get_errors(skx, sample)
    for error in error_names:
        sum = 0
        for bin in convert(skx, sample, error):
            sum= sum + (bin)
        sum_ratio.append(sum)
    return sum_ratio

def get_sys_count(skx):
    tfile = ROOT.TFile('./sys_pdf/error.sk' + str(skx) + '.root', 'READ')
    keys = tfile.GetListOfKeys()
    sum = 0
    for key in keys:
        th2f = tfile.Get(key.GetName())
        if 'neg' in key.GetName():
            sum = sum - th2f.Integral()
        else:
            sum = sum + th2f.Integral()
    print sum

def plot_max():
    ROOT.gROOT.Macro( os.path.expanduser( '~/Root/rootlogon.C' ) )
    th1 = ROOT.TH1F('fij max','fij max', 60, 0, 0.6)
    th2 = ROOT.TH1F('fij sum','fij sum', 60, 0, 30)
    th2d = ROOT.TH2F('fij 2D', 'fij 2D', 60, 0, 0.6, 60, 0, 30)

    c1 = ROOT.TCanvas('mycanvas','mycanvas', 800, 600)
    for sk in range(1,5):
        for sample in ['fcmc', 'tau']:
            ratio_list = get_max(sk, sample)
            sum_list = get_sum(sk, sample)
            for ratio , sum in zip(ratio_list, sum_list):
                th1.Fill(ratio)
                th2.Fill(sum)
                th2d.Fill(ratio, sum)
    th1.Draw()
    th1.GetXaxis().SetTitle('max bin(ratio) in 2D fijs')
    print th1.GetEntries()
    c1.Print('fij_max.png','png')
    th2.Draw()
    th2.GetXaxis().SetTitle('sum of 2D fijs')
    c1.Print('fij_sum.png','png')
    th2d.Draw('colz')
    th2d.GetXaxis().SetTitle('max bin of 2D fijs')
    th2d.GetYaxis().SetTitle('sum of 2D fijs')
    c1.Print('fij_max_sum.png','png')

if __name__ == '__main__':
    for skx in range(1,5):
        get_sys_count(skx)
