import ROOT, os

def plot_fit():
    ROOT.gROOT.Macro( os.path.expanduser( '~/Root/rootlogon.C' ) )
    infile = open('apr16.txt','r')
    i = 1
    flag = 0
    fit_sys = []
    for line in infile:
        if 'NAME' in line:
            flag = flag + 1
        if 'CONVERGED' in line:
            flag = flag + 1
        if 'EXTERNAL' in line:
            flag = flag + 1
        if flag == 4 and 'cc_nutau' not in line and 'NAME' not in line:
            if len(line.strip()) < 30:
                continue
            words = line.strip().split()
            i = i+1
            fit_sys.append(words)

    data_file = ROOT.TFile('./data/SK-merged.root','READ')
    bkg_2D = data_file.Get('bkgHistoZenith2D')
    tau_2D = data_file.Get('tauHistoZenith2D')
    data_2D = data_file.Get('dataHistoZenith2D')

    sys_file = ROOT.TFile('./sys_pdf/error.merged.root','READ')
    keys = sys_file.GetListOfKeys()
    for key in keys:
        for i in range(len(fit_sys)):
            if key.GetName()[0:-4] == fit_sys[i][1]:
                _hist = sys_file.Get(key.GetName())
                fit_value = float(fit_sys[i][2])
                if key.GetName()[-3:] == 'pos':
                    bkg_2D.Add(_hist, fit_value)
                else:
                    bkg_2D.Add(_hist,-fit_value)

    tau_taulike = ROOT.TH1F('tau_taulike', 'tau_taulike', 15, -1, 1)
    bkg_taulike = ROOT.TH1F('bkg_taulike', 'bkg_taulike', 15, -1, 1)
    data_taulike = ROOT.TH1F('data_taulike', 'data_taulike', 15, -1, 1)
    for i in range(9,16):
        for j in range(15):
            tau_taulike.SetBinContent(j+1, tau_taulike.GetBinContent(j+1) + 1.52*tau_2D.GetBinContent(j+1, i))
            bkg_taulike.SetBinContent(j+1, bkg_taulike.GetBinContent(j+1) + bkg_2D.GetBinContent(j+1, i))
            data_taulike.SetBinContent(j+1, data_taulike.GetBinContent(j+1) + data_2D.GetBinContent(j+1, i))

    tau_nontau = ROOT.TH1F('tau_nontau', 'tau_nontau', 15, -1, 1)
    bkg_nontau = ROOT.TH1F('bkg_nontau', 'bkg_nontau', 15, -1, 1)
    data_nontau = ROOT.TH1F('data_nontau', 'data_nontau', 15, -1, 1)
    for i in range(1,9):
        for j in range(15):
            tau_nontau.SetBinContent(j+1, tau_nontau.GetBinContent(j+1) + 1.52*tau_2D.GetBinContent(j+1, i))
            bkg_nontau.SetBinContent(j+1, bkg_nontau.GetBinContent(j+1) + bkg_2D.GetBinContent(j+1, i))
            data_nontau.SetBinContent(j+1, data_nontau.GetBinContent(j+1) + data_2D.GetBinContent(j+1, i))

    tau_upward = ROOT.TH1F('tau_upward', 'tau_upward', 15, -0.1, 1.1)
    bkg_upward = ROOT.TH1F('bkg_upward', 'bkg_upward', 15, -0.1, 1.1)
    data_upward = ROOT.TH1F('data_upward', 'data_upward', 15, -0.1, 1.1)
    for i in range(1,9):
        for j in range(15):
            tau_upward.SetBinContent(j+1, tau_upward.GetBinContent(j+1) + 1.52*tau_2D.GetBinContent(i, j+1))
            bkg_upward.SetBinContent(j+1, bkg_upward.GetBinContent(j+1) + bkg_2D.GetBinContent(i, j+1))
            data_upward.SetBinContent(j+1, data_upward.GetBinContent(j+1) + data_2D.GetBinContent(i, j+1))

    tau_downward = ROOT.TH1F('tau_downward', 'tau_downward', 15, -0.1, 1.1)
    bkg_downward = ROOT.TH1F('bkg_downward', 'bkg_downward', 15, -0.1, 1.1)
    data_downward = ROOT.TH1F('data_downward', 'data_downward', 15, -0.1, 1.1)
    for i in range(9,15):
        for j in range(15):
            tau_downward.SetBinContent(j+1, tau_downward.GetBinContent(j+1) + 1.52*tau_2D.GetBinContent(i, j+1))
            bkg_downward.SetBinContent(j+1, bkg_downward.GetBinContent(j+1) + bkg_2D.GetBinContent(i, j+1))
            data_downward.SetBinContent(j+1, data_downward.GetBinContent(j+1) + data_2D.GetBinContent(i, j+1))

    c1 = ROOT.TCanvas('mycanvas', 'mycanvas', 800,600)
    c1.Divide(2,2)
    c1.cd(1)
    leg = ROOT.TLegend(0.2, 0.7,0.45,0.9)
    leg.SetBorderSize(0)
    leg.AddEntry(bkg_taulike, 'Background after fit', 'f')
    leg.AddEntry(tau_taulike, 'Tau after fit', 'f')
    leg.AddEntry(data_taulike, 'Data', 'lep')
    ts1 = ROOT.THStack('taulike', 'taulike')
    ts1.Add(bkg_taulike)
    ts1.Add(tau_taulike)
    ts1.SetMaximum(data_taulike.GetMaximum()*1.2)
    tau_taulike.SetFillColor(ROOT.kRed)
    ts1.Draw('hist')
    ts1.GetXaxis().SetTitle('Cos Zenith Angle')
    ts1.GetYaxis().SetTitle('# of events')
    data_taulike.Draw('E1same')
    leg.Draw()
    data_taulike.SetLineColor(ROOT.kBlack)

    c1.cd(2)
    ts2 = ROOT.THStack('upward', 'upward')
    ts2.Add(bkg_upward)
    ts2.Add(tau_upward)
    tau_upward.SetFillColor(ROOT.kRed)
    ts2.Draw('hist')
    ts2.SetMaximum(data_upward.GetMaximum()*1.2)
    ts2.GetXaxis().SetTitle('NN output')
    ts2.GetYaxis().SetTitle('# of events')
    data_upward.Draw('E1same')
    data_upward.SetLineColor(ROOT.kBlack)

    c1.cd(3)
    ts3 = ROOT.THStack('nontau', 'nontau')
    ts3.Add(bkg_nontau)
    ts3.Add(tau_nontau)
    tau_nontau.SetFillColor(ROOT.kRed)
    ts3.Draw('hist')
    ts3.SetMaximum(data_nontau.GetMaximum()*1.2)
    ts3.GetXaxis().SetTitle('Cos Zenith Angle')
    ts3.GetYaxis().SetTitle('# of events')
    data_nontau.Draw('E1same')
    data_nontau.SetLineColor(ROOT.kBlack)

    c1.cd(4)
    ts4 = ROOT.THStack('downward', 'downward')
    ts4.Add(bkg_downward)
    ts4.Add(tau_downward)
    tau_downward.SetFillColor(ROOT.kRed)
    ts4.Draw('hist')
    ts4.SetMaximum(data_downward.GetMaximum()*1.2)
    ts4.GetXaxis().SetTitle('NN output')
    ts4.GetYaxis().SetTitle('# of events')
    data_downward.Draw('E1same')
    data_downward.SetLineColor(ROOT.kBlack)
    c1.Print('test.png','png')

def plot_sys():
    ROOT.gROOT.Macro( os.path.expanduser( '~/Root/rootlogon.C' ) )
    infile = open('apr16.txt','r')
    th1 = ROOT.TH1F('sys_error', 'Systematic errors', 40, -10, 10)
    th2 = ROOT.TH1F('Systematic errors', 'Systematic errors', 20, -1, 1)
    th3 = ROOT.TH1F('errors', 'errors', 40, 0, 40)
    th3.SetTitle('')
    i = 1
    flag = 0
    for line in infile:
        if 'NAME' in line:
            flag = flag + 1
        if 'CONVERGED' in line:
            flag = flag + 1
        if 'EXTERNAL' in line:
            flag = flag + 1
        if flag == 4 and 'cc_nutau' not in line and 'NAME' not in line:
            if len(line.strip()) < 30:
                continue
            words = line.strip().split()
            th3.GetXaxis().SetBinLabel(i, words[1])
            th3.SetLabelSize(0.02)
            th3.SetBinContent(i, float(words[2]))
            th1.Fill((float(words[2]) )/(1 - float(words[3]) ))
            th2.Fill((float(words[2]) ))
            th3.GetYaxis().SetTitle('unit of #sigma')
            i = i+1

    #print th1.GetMean(), th1.GetRMS()
    #print th2.GetMean(), th2.GetRMS()
    c1 = ROOT.TCanvas('mycanvas', 'mycanvas', 800,600)
    c1.SetBottomMargin(0.35)
    th3.Draw()
    th3.GetXaxis().SetLabelSize(0.035)
    th3.GetXaxis().LabelsOption('v')
    c1.Print('error_terms.png','png')
if __name__ == '__main__':
    plot_fit()
