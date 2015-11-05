import ROOT, os

def plot_fit():
    ROOT.gROOT.Macro( os.path.expanduser( '~/Root/rootlogon.C' ) )
    infile = open('fit.txt','r')
    th1 = ROOT.TH1F('Systematic errors', 'Systematic errors', 102, 0, 102)
    i = 1
    errors = []
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
            th1.GetXaxis().SetBinLabel(i, words[1])
            th1.SetLabelSize(0.02)
            th1.SetBinContent(i, float(words[2]))
            th1.SetBinError(i, float(words[3]))
            th1.GetYaxis().SetTitle('unit of #sigma')
            i = i+1

    c1 = ROOT.TCanvas('mycanvas', 'mycanvas', 800,600)
    th1.Draw()
    c1.Print('sys_fit.pdf','pdf')

if __name__ == '__main__':
    plot_fit()
