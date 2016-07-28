import ROOT, os
import array

def plot_fit():
    ROOT.gROOT.Macro( os.path.expanduser( '~/Root/rootlogon.C' ) )
    infile = open('refit.txt','r')
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

    ROOT.gStyle.SetOptTitle(0)
    c1 = ROOT.TCanvas('mycanvas', 'mycanvas', 800,600)
    c1.Divide(2,2)
    c1.cd(1)
    leg = ROOT.TLegend(0.2, 0.7,0.45,0.9)
    leg.SetBorderSize(0)
    leg.AddEntry(bkg_taulike, 'BG after fit', 'f')
    leg.AddEntry(tau_taulike, 'Tau after fit', 'f')
    leg.AddEntry(data_taulike, 'Data', 'lep')
    latex1 = ROOT.TLatex(-0.3, 0.2*bkg_taulike.GetMaximum(), 'Tau-like')
    latex1.SetTextSize(0.1)
    latex1.SetTextAlign(13)
    ts1 = ROOT.THStack('taulike', 'taulike')
    ts1.Add(bkg_taulike)
    ts1.Add(tau_taulike)
    bkg_taulike.SetLineColor(ROOT.kBlack)
    ts1.SetMaximum(data_taulike.GetMaximum()*1.2)
    tau_taulike.SetFillColor(ROOT.kRed)
    tau_taulike.SetLineColor(ROOT.kRed)
    ts1.Draw('hist')
    ts1.GetXaxis().SetTitle('Cos#Theta')
    ts1.GetXaxis().CenterTitle()
    ts1.GetYaxis().SetTitle('Events')
    data_taulike.Draw('E1same')
    data_taulike.SetMarkerStyle(24)
    data_taulike.SetMarkerSize(0.9)
    leg.Draw()
    latex1.Draw()
    data_taulike.SetLineColor(ROOT.kBlack)

    c1.cd(2)
    ts2 = ROOT.THStack('upward', 'upward')
    ts2.Add(bkg_upward)
    ts2.Add(tau_upward)
    leg2 = ROOT.TLegend(0.6, 0.7,0.85,0.9)
    leg2.SetBorderSize(0)
    leg2.AddEntry(bkg_upward, 'BG after fit', 'f')
    leg2.AddEntry(tau_upward, 'Tau after fit', 'f')
    leg2.AddEntry(data_upward, 'Data', 'lep')
    tau_upward.SetFillColor(ROOT.kGray)
    tau_upward.SetLineColor(ROOT.kGray)
    bkg_upward.SetLineColor(ROOT.kBlack)
    latex2 = ROOT.TLatex(0.2,0.2*bkg_upward.GetMaximum(), 'Upward')
    latex2.SetTextSize(0.1)
    latex2.SetTextAlign(13)
    ts2.Draw('hist')
    ts2.SetMaximum(data_upward.GetMaximum()*1.2)
    ts2.GetXaxis().SetTitle('NN output')
    ts2.GetXaxis().CenterTitle()
    ts2.GetYaxis().SetTitle('Events')
    ts2.GetYaxis().SetNoExponent()
    data_upward.Draw('E1same')
    data_upward.SetMarkerSize(0.9)
    data_upward.SetMarkerStyle(24)
    #leg2.Draw()
    latex2.Draw()
    data_upward.SetLineColor(ROOT.kBlack)

    c1.cd(3)
    ts3 = ROOT.THStack('nontau', 'nontau')
    leg3 = ROOT.TLegend(0.6, 0.7,0.85,0.9)
    leg3.SetBorderSize(0)
    leg3.AddEntry(bkg_nontau, 'BG after fit', 'f')
    leg3.AddEntry(tau_nontau, 'Tau after fit', 'f')
    leg3.AddEntry(data_nontau, 'Data', 'lep')
    ts3.Add(bkg_nontau)
    ts3.Add(tau_nontau)
    tau_nontau.SetFillColor(ROOT.kGray)
    tau_nontau.SetLineColor(ROOT.kGray)
    bkg_nontau.SetLineColor(ROOT.kBlack)
    latex3 = ROOT.TLatex(-0.3,0.2*bkg_nontau.GetMaximum(), 'Non tau-like')
    latex3.SetTextSize(0.1)
    latex3.SetTextAlign(13)
    ts3.Draw('hist')
    ts3.SetMaximum(data_nontau.GetMaximum()*1.2)
    ts3.GetXaxis().SetTitle('Cos#Theta')
    ts3.GetYaxis().SetTitle('Events')
    data_nontau.Draw('E1same')
    data_nontau.SetMarkerSize(0.9)
    data_nontau.SetMarkerStyle(24)
    data_nontau.SetLineColor(ROOT.kBlack)
    latex3.Draw()

    c1.cd(4)
    ts4 = ROOT.THStack('downward', 'downward')
    ts4.Add(bkg_downward)
    ts4.Add(tau_downward)
    latex4 = ROOT.TLatex(0.2,0.2*bkg_downward.GetMaximum(), 'Downward')
    latex4.SetTextSize(0.1)
    latex4.SetTextAlign(13)
    tau_downward.SetFillColor(ROOT.kGray)
    tau_downward.SetLineColor(ROOT.kGray)
    bkg_downward.SetLineColor(ROOT.kBlack)
    ts4.Draw('hist')
    ts4.SetMaximum(data_downward.GetMaximum()*1.2)
    ts4.GetXaxis().SetTitle('NN output')
    ts4.GetXaxis().CenterTitle()
    ts4.GetYaxis().SetTitle('Events')
    ts4.GetYaxis().SetNoExponent()
    data_downward.Draw('E1same')
    data_downward.SetMarkerSize(0.9)
    data_downward.SetMarkerStyle(24)
    data_downward.SetLineColor(ROOT.kBlack)
    latex4.Draw()
    c1.Print('fourpanel.png','png')

    c2 = ROOT.TCanvas('mycanvas2', 'mycanvas2', 800,600)
    ts1.Draw('hist')
    ts1.GetXaxis().SetTitle('Cos#Theta')
    ts1.GetYaxis().SetTitle('Events')
    data_taulike.Draw('E1same')
    data_taulike.SetMarkerStyle(24)
    data_taulike.SetMarkerSize(0.9)
    leg.Draw()
    data_taulike.SetLineColor(ROOT.kBlack)
    c2.Print('tau_like_red.png','png')

def store_sys():
    infile = open('refit.txt','r')
    i = 1
    flag = 0
    fit_sys = ROOT.TTree('fit_sys','fit_sys')
    strings = ROOT.vector('string')()
    fit_mean = ROOT.vector('float')()
    fit_sigma = ROOT.vector('float')()
    fit_sys.Branch('sys_term', strings)
    fit_sys.Branch('sys_mean', fit_mean)
    fit_sys.Branch('sys_sigma', fit_sigma)
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
            strings.push_back(words[1])
            fit_mean.push_back(float(words[2]))
            fit_sigma.push_back(float(words[3]))
            i = i+1
    fit_sys.Fill()
    fit_sys.Scan()
    root_sys = ROOT.TFile('fit_sys.root','RECREATE')
    fit_sys.Write()
    root_sys.Close()

def plot_sys():
    ROOT.gROOT.Macro( os.path.expanduser( '~/Root/rootlogon.C' ) )
    ROOT.gStyle.SetOptTitle(0)
    infile = open('refit.txt','r')
    th1 = ROOT.TH1F('sys_error', 'Systematic errors', 20, -3, 3)
    th2 = ROOT.TH1F('Systematic errors', 'Systematic errors', 20, -2, 2)
    th3 = ROOT.TH1F('errors', 'errors', 43, 0, 43)
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

    print th1.GetMean(), th1.GetRMS()
    print th2.GetMean(), th2.GetRMS()
    c1 = ROOT.TCanvas('mycanvas', 'mycanvas', 800,600)
    c1.SetBottomMargin(0.35)
    th2.Draw()
    th2.SetXTitle('Fit Error')
    #th3.GetXaxis().SetLabelSize(0.035)
    #th3.GetXaxis().LabelsOption('v')
    c1.Print('error.png','png')
    th1.Draw()
    #th1.Fit('gaus')
    th1.SetXTitle('Corrected Fit Error')
    c1.Print('error_corrected.png','png')
    th3.Draw()
    th3.GetXaxis().SetLabelSize(0.035)
    th3.GetXaxis().LabelsOption('v')
    c1.Print('error_terms.png','png')
if __name__ == '__main__':
    store_sys()
