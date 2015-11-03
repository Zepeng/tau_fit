#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TParameter.h>
#include <TKey.h>
#include <TH1.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTree.h>

#include <RooWorkspace.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooKeysPdf.h>
#include <RooRandom.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooNDKeysPdf.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooGaussian.h>
#include <RooMinuit.h>
#include <RooAddition.h>

using namespace RooFit;
using namespace std;

void joint_fit_sys(char* job)
{

  // Define variables used in the fit
  RooRealVar cosTheta      ("cosTheta",    "Zenith of the Cosine Angle", -1.0, 1.0);
  RooRealVar NN_output     ("NN_output",   "NeuralNet Output",           -0.1, 1.1);
  RooRealVar NN_selected   ("NN_selected", "NeuralNet selected",          0.0, 1.0);
  RooRealVar eventWeight   ("eventWeight", "Event Weight",                0.0, 10.0);
  RooArgSet  chainVars     (cosTheta, NN_output, NN_selected, eventWeight);
  RooArgList axisVariables (cosTheta, NN_output);

  TString objDir = "./data/" ;

  // Open SKI ,SKII ,SKIII and SKIV files
  TFile *skIfile    = new TFile((objDir +"SK-I.root"    ).Data(),  "READ");
  TFile *skIIfile   = new TFile((objDir +"SK-II.root"   ).Data(),  "READ");
  TFile *skIIIfile  = new TFile((objDir +"SK-III.root"  ).Data(),  "READ");
  TFile *skIVfile   = new TFile((objDir +"SK-IV.root"   ).Data(),  "READ");

  // Grab the trees and Histograms we need.
  // Data chain could be used in building KeysPdf, but is not used in this
  // fitting now.
  //SK1
  TTree *SKIChain              = (TTree*)skIfile->Get("dataHistoTMVAOutputTree");
  //TTree *SKIBkgChain           = (TTree*)skIfile->Get("bkgHistoTMVAOutputTree");
  //TTree *SKISigChain           = (TTree*)skIIIfile->Get("tauHistoTMVAOutputTree");
  TH2F  *tauHisto2DZenithSKI   = (TH2F*) skIfile->Get("tauHistoZenith2D");
  TH2F  *bkgHisto2DZenithSKI   = (TH2F*) skIfile->Get("bkgHistoZenith2D");

  //SK2
  TTree *SKIIChain              = (TTree*)skIIfile->Get("dataHistoTMVAOutputTree");
  //TTree *SKIIBkgChain           = (TTree*)skIIfile->Get("bkgHistoTMVAOutputTree");
  //TTree *SKIISigChain           = (TTree*)skIIfile->Get("tauHistoTMVAOutputTree");
  TH2F  *tauHisto2DZenithSKII   = (TH2F*) skIIfile->Get("tauHistoZenith2D");
  TH2F  *bkgHisto2DZenithSKII   = (TH2F*) skIIfile->Get("bkgHistoZenith2D");
  //SK3
  TTree *SKIIIChain              = (TTree*)skIIIfile->Get("dataHistoTMVAOutputTree");
  //TTree *SKIIIBkgChain           = (TTree*)skIIIfile->Get("bkgHistoTMVAOutputTree");
  //TTree *SKIIISigChain           = (TTree*)skIIIfile->Get("tauHistoTMVAOutputTree");
  TH2F  *tauHisto2DZenithSKIII   = (TH2F*) skIIIfile->Get("tauHistoZenith2D");
  TH2F  *bkgHisto2DZenithSKIII   = (TH2F*) skIIIfile->Get("bkgHistoZenith2D");
  //SK4
  TTree *SKIVChain              = (TTree*)skIVfile->Get("dataHistoTMVAOutputTree");
  //TTree *SKIVBkgChain           = (TTree*)skIVfile->Get("bkgHistoTMVAOutputTree");
  //TTree *SKIVSigChain           = (TTree*)skIVfile->Get("tauHistoTMVAOutputTree");
  TH2F  *tauHisto2DZenithSKIV   = (TH2F*) skIVfile->Get("tauHistoZenith2D");
  TH2F  *bkgHisto2DZenithSKIV   = (TH2F*) skIVfile->Get("bkgHistoZenith2D");

  //Kernel estimation PDF, not used temporarily.
  /*
  RooDataSet RoofitBkgSetI("bkgTreeI", "Bkg ", chainVars, RooFit::Import(*SKIBkgChain), 
			   RooFit::Cut("NN_selected==1"), RooFit::WeightVar("eventWeight"));
  RooNDKeysPdf bkgKeysPdfI("bkgpdfI",     "Bkg I",   axisVariables, RoofitBkgSetI ,"m");
  RooDataSet RoofitSigSetI("sigTreeI", "Sig ", chainVars, RooFit::Import(*SKISigChain), 
			   RooFit::Cut("NN_selected==1"), RooFit::WeightVar("eventWeight"));
  RooNDKeysPdf sigKeysPdfI("sigpdfI", "Sig I", axisVariables, RoofitSigSetI, "m");
  */

  // Read the data set from TTree, and cut the events by NN_selected.
  // SK1
  RooDataSet  dataSetI  ("dset I", "SKI Data", chainVars, RooFit::Import(*SKIChain), 
			 RooFit::Cut("NN_selected==1"));
  RooDataHist bkgHistoI("bkgI", "Bkg Histo", 
                               axisVariables, bkgHisto2DZenithSKI);
  RooDataHist sigHistoI("sigI", "Sig Histo", 
                               axisVariables, tauHisto2DZenithSKI);
  // SK2
  RooDataSet  dataSetII  ("dset II", "SKII Data", chainVars, RooFit::Import(*SKIIChain), 
			 RooFit::Cut("NN_selected==1"));
  RooDataHist bkgHistoII("bkgII", "Bkg Histo", 
                               axisVariables, bkgHisto2DZenithSKII);
  RooDataHist sigHistoII("sigII", "Sig Histo", 
                               axisVariables, tauHisto2DZenithSKII);
  // SK3
  RooDataSet  dataSetIII  ("dset III", "SKIII Data", chainVars, RooFit::Import(*SKIIIChain), 
			 RooFit::Cut("NN_selected==1"));
  RooDataHist bkgHistoIII("bkgIII", "Bkg Histo", 
                               axisVariables, bkgHisto2DZenithSKIII);
  RooDataHist sigHistoIII("sigIII", "Sig Histo", 
                               axisVariables, tauHisto2DZenithSKIII);
  // SK4
  RooDataSet  dataSetIV  ("dset IV", "SKIV Data", chainVars, RooFit::Import(*SKIVChain), 
			 RooFit::Cut("NN_selected==1"));
  RooDataHist bkgHistoIV("bkgIV", "Bkg Histo", 
                               axisVariables, bkgHisto2DZenithSKIV);
  RooDataHist sigHistoIV("sigIV", "Sig Histo", 
                               axisVariables, tauHisto2DZenithSKIV);

  // Set the expected number of signal/background events (these are constants)
  // These are not normalized by any fraction (signal, background or DIS).
  RooRealVar cc_nutau_xsec("cc_nutau_xsec", "cc_nutau_xsec",  0, 3);
  //SK1
  RooRealVar expbkgI("expbkgI", "expbkgI",  bkgHisto2DZenithSKI->Integral());
  RooRealVar expsigI("expsigI", "expsigI",  tauHisto2DZenithSKI->Integral());
  RooFormulaVar tau_normI("tau_normI", "expsigI*cc_nutau_xsec",  RooArgList(expsigI, cc_nutau_xsec));
  //SK2
  RooRealVar expbkgII("expbkgII", "expbkgII",  bkgHisto2DZenithSKII->Integral());
  RooRealVar expsigII("expsigII", "expsigII",  tauHisto2DZenithSKII->Integral());
  RooFormulaVar tau_normII("tau_normII", "expsigII*cc_nutau_xsec",  RooArgList(expsigII, cc_nutau_xsec));
  //SK3
  RooRealVar expbkgIII("expbkgIII", "expbkgIII",  bkgHisto2DZenithSKIII->Integral());
  RooRealVar expsigIII("expsigIII", "expsigIII",  tauHisto2DZenithSKIII->Integral());
  RooFormulaVar tau_normIII("tau_normIII", "expsigIII*cc_nutau_xsec",  RooArgList(expsigIII, cc_nutau_xsec));
  //SK4
  RooRealVar expbkgIV("expbkgIV", "expbkgIV",  bkgHisto2DZenithSKIV->Integral());
  RooRealVar expsigIV("expsigIV", "expsigIV",  tauHisto2DZenithSKIV->Integral());
  RooFormulaVar tau_normIV("tau_normIV", "expsigIV*cc_nutau_xsec",  RooArgList(expsigIV, cc_nutau_xsec));

  // Make base PDFs with background simulation(normal ATM neutrinos, except for tau), and tau neutrino
  // simulation. Both PDFs have to be normalized to livetime, because they are fixed in the final PDF,
  // only the systematic errors are varied in the fitting.
  // SK1
  RooHistPdf  bkgPdfI("bkgpdfI", "Background PDF ", 
                             axisVariables, bkgHistoI, 0);
  RooHistPdf  sigPdfI("sigpdfI", "Signal PDF ",     
                             axisVariables, sigHistoI, 0);
  // SK2
  RooHistPdf  bkgPdfII("bkgpdfII", "Background PDF ", 
                             axisVariables, bkgHistoII, 0);
  RooHistPdf  sigPdfII("sigpdfII", "Signal PDF ",     
                             axisVariables, sigHistoII, 0);
  // SK3
  RooHistPdf  bkgPdfIII("bkgpdfIII", "Background PDF ", 
                             axisVariables, bkgHistoIII, 0);
  RooHistPdf  sigPdfIII("sigpdfIII", "Signal PDF ",     
                             axisVariables, sigHistoIII, 0);
  // SK4
  RooHistPdf  bkgPdfIV("bkgpdfIV", "Background PDF ", 
                             axisVariables, bkgHistoIV, 0);
  RooHistPdf  sigPdfIV("sigpdfIV", "Signal PDF ",     
                             axisVariables, sigHistoIV, 0);

  // Create a list of PDFs and coefficients that will be used to construct the final PDF.
  // SK1
  RooArgList pdf_listI("pdf_listI");
  pdf_listI.add(bkgPdfI);
  pdf_listI.add(sigPdfI);
  RooArgList coeff_listI("coeff_listI");
  coeff_listI.add(expbkgI);
  coeff_listI.add(tau_normI);
  // SK2
  RooArgList pdf_listII("pdf_listII");
  pdf_listII.add(bkgPdfII);
  pdf_listII.add(sigPdfII);
  RooArgList coeff_listII("coeff_listII");
  coeff_listII.add(expbkgII);
  coeff_listII.add(tau_normII);
  // SK3
  RooArgList pdf_listIII("pdf_listIII");
  pdf_listIII.add(bkgPdfIII);
  pdf_listIII.add(sigPdfIII);
  RooArgList coeff_listIII("coeff_listIII");
  coeff_listIII.add(expbkgIII);
  coeff_listIII.add(tau_normIII);
  // SK4
  RooArgList pdf_listIV("pdf_listIV");
  pdf_listIV.add(bkgPdfIV);
  pdf_listIV.add(sigPdfIV);
  RooArgList coeff_listIV("coeff_listIV");
  coeff_listIV.add(expbkgIV);
  coeff_listIV.add(tau_normIV);

  //SK1-4 MC study - statistical only
  if(kFALSE)
  {
    RooRandom::randomGenerator()->SetSeed(0);
    TFile* mc_study = new TFile("mc_tau.root","recreate");
    TH1F* mc_xsecI   = new TH1F("xsec1","xsec",40,0,4);
    TH1F* mc_xsecII  = new TH1F("xsec2","xsec",40,0,4);
    TH1F* mc_xsecIII = new TH1F("xsec3","xsec",40,0,4);
    TH1F* mc_xsecIV  = new TH1F("xsec4","xsec",40,0,4);
    TH1F* mc_xsec_com  = new TH1F("xsec_com","xsec",40,0,4);
    TRandom3* r3 = new TRandom3();
    RooAddPdf modelI_stat ("modelI",  "signal+bkgd SKI",  pdf_listI, coeff_listI );
    RooAddPdf modelII_stat ("modelII",  "signal+bkgd SKII",  pdf_listII, coeff_listII );
    RooAddPdf modelIII_stat ("modelIII",  "signal+bkgd SKIII",  pdf_listIII, coeff_listIII );
    RooAddPdf modelIV_stat ("modelIV",  "signal+bkgd SKIV",  pdf_listIV, coeff_listIV );
    for(int i = 0; i < 1000; i++)// Create MC sample for study.
    {
        RooDataSet  *mc_bkgI = bkgPdfI.generate(axisVariables, TMath::Nint(r3->PoissonD(2817)));
        RooDataSet  *mc_sigI = sigPdfI.generate(axisVariables, TMath::Nint(r3->PoissonD(56)));
        mc_sigI->append(*mc_bkgI);
        RooDataSet  *mc_bkgII = bkgPdfII.generate(axisVariables, TMath::Nint(r3->PoissonD(1496)));
        RooDataSet  *mc_sigII = sigPdfII.generate(axisVariables, TMath::Nint(r3->PoissonD(32)));
        mc_sigII->append(*mc_bkgII);
        RooDataSet  *mc_bkgIII = bkgPdfIII.generate(axisVariables, TMath::Nint(r3->PoissonD(988)));
        RooDataSet  *mc_sigIII = sigPdfIII.generate(axisVariables, TMath::Nint(r3->PoissonD(19)));
        mc_sigIII->append(*mc_bkgIII);
        RooDataSet  *mc_bkgIV = bkgPdfIV.generate(axisVariables, TMath::Nint(r3->PoissonD(3343)));
        RooDataSet  *mc_sigIV = sigPdfIV.generate(axisVariables, TMath::Nint(r3->PoissonD(68)));
        mc_sigIV->append(*mc_bkgIV);
        //Create likelihood function for each SK period.
        RooAbsReal* nll_sk1 = modelI_stat.createNLL(*mc_sigI,  RooFit::Extended(kTRUE));
        RooAbsReal* nll_sk2 = modelII_stat.createNLL(*mc_sigII, RooFit::Extended(kTRUE));
        RooAbsReal* nll_sk3 = modelIII_stat.createNLL(*mc_sigIII, RooFit::Extended(kTRUE));
        RooAbsReal* nll_sk4 = modelIV_stat.createNLL(*mc_sigIV, RooFit::Extended(kTRUE));
        // Simultaneous fit of SK1-4 with combined likelihood.
        RooAddition nllCombined("nll_combined", "nll_combined", RooArgSet(*nll_sk1, *nll_sk2, *nll_sk3, *nll_sk4));
       
        //Mininize the likelihood for SK1-4 individually and simultaneously.
        RooMinuit minit1(*nll_sk1);
        minit1.migrad();
        mc_xsecI->Fill(cc_nutau_xsec.getVal());
        RooMinuit minit2(*nll_sk2);
        minit2.migrad();
        mc_xsecII->Fill(cc_nutau_xsec.getVal());
        RooMinuit minit3(*nll_sk3);
        minit3.migrad();
        mc_xsecIII->Fill(cc_nutau_xsec.getVal());
        RooMinuit minit4(*nll_sk4);
        minit4.migrad();
        mc_xsecIV->Fill(cc_nutau_xsec.getVal());
        RooMinuit minit_com(nllCombined);
        minit_com.migrad();
        mc_xsec_com->Fill(cc_nutau_xsec.getVal());
    }
    mc_study->Write();
    return;
  }
  // Make PDFs for systematic errors. The PDF is built with the change in each bin after shifting
  // systematic error by one sigma. Here, we assume the change corresponding to the systematic error
  // is linear. Each systematic error consists of two histograms, positive and negative part. This is
  // due to the fact that PDF has to be positive. The PDF of systematic errors is added with a varied
  // parameter to the base PDF. 
  // PDF_final = PDF_bkg + PDF_tau + Sum_i(var_sys_i*PDF_pos - var_sys_i*PDF_neg)
  //
  //systematic errors that are used in the SK Osc3++ analysis. skip t2k related errors. Each error is
  //constrained with a Gaussian function centered at 0 and sigma=1. CC_nutau_xsec is not constrained.
  RooWorkspace* tau_fit= new RooWorkspace("tau_fit");
  
  // Read the file of histograms for SK1 systematic errors.
  TFile* fijs_sk1 = new TFile("./sys_pdf/error.sk1.root", "READ");

  // Read the names of systematic errors from fijs file, and store in a vector of string.
  // create Gaussian constraint for each systematic error.
  std::vector<std::string> sk1_errors;
  TIter next(fijs_sk1->GetListOfKeys());
  TKey* key;
  RooArgSet parts_pdfI("parts_pdfI");
  while ((key = (TKey*)next())) {
      std::string fij_name = key->GetName();
      std::string error_term = fij_name.substr(0,fij_name.length()-4) ;
      std::vector<std::string>::iterator iter;
      iter = find(sk1_errors.begin(), sk1_errors.end(), error_term);
      TH2F* th2 = (TH2F*)fijs_sk1->Get(fij_name.c_str());
      if(iter == sk1_errors.end() && error_term != "CC_nutau_xsec" && th2->Integral() > 0)
      // read fijs for each error, and use the untrivial fijs. DO NOT constrain CC_nutau_xsec.
      {
        sk1_errors.push_back(error_term);
        // Define the variables for systematic errors, 1 means shifting the systematic error by 1 sigma
        // // The PDFs of shifting 1 sigma is built with Osc3++. 
        TString constraint_expr = TString::Format("RooGaussian::%s_sys(%s[-3,3], 0, 1.0)", 
                error_term.c_str(), error_term.c_str());//Gaussian constraint.
        tau_fit->factory(constraint_expr);
        TString constraint_name = TString::Format("%s_sys", error_term.c_str());
        parts_pdfI.add(*tau_fit->pdf(constraint_name));
      }
  }
 
  //Create RooHistPDF for each systematic error. Divide each term to positive and
  //and negative parts.
  //SK1 positive
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_posI;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_posI;
  for(unsigned int i = 0; i < sk1_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk1_errors[i] + "_pos";
    TH2F* hist_temp = (TH2F*)fijs_sk1->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_posI[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_posI[i].reset(new RooFormulaVar(coeff_name.c_str(), "@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk1_errors[i].c_str())) ));
    if(hist_temp->Integral() > 0)// Only keep the non-trivial PDFs.
    {
        pdf_listI.add(*hist_pdfs_posI[i]);
        coeff_listI.add(*coeffs_posI[i]);
    }
  }
  //SK1 negative
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_negI;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_negI;
  for(unsigned int i = 0; i < sk1_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk1_errors[i] + "_neg";
    TH2F* hist_temp = (TH2F*)fijs_sk1->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_negI[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_negI[i].reset(new RooFormulaVar(coeff_name.c_str(), "-1*@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk1_errors[i].c_str())) ));
    if(hist_temp->Integral()>0)
    {
        pdf_listI.add(*hist_pdfs_negI[i]);
        coeff_listI.add(*coeffs_negI[i]);
    }
  }
  //stitch the parts of PDF to the final PDF.
  RooAddPdf modelI ("modelI",  "signal+bkgd SKI",  pdf_listI, coeff_listI );//base PDFs + PDFs of sys errors
  RooProdPdf modelI_sys("modelI with sys errors", "modelI with sys errors", parts_pdfI);//Gaussian constraint of sys errors
  
  // Read the file of histograms for sk2 systematic errors.
  TFile* fijs_sk2 = new TFile("./sys_pdf/error.sk2.root", "READ");

  // Read the names of systematic errors from fijs file, and store in a vector of string.
  // create Gaussian constraint for each systematic error.
  std::vector<std::string> sk2_errors;
  TIter next2(fijs_sk2->GetListOfKeys());
  RooArgSet parts_pdfII("parts_pdfII");
  while ((key = (TKey*)next2())) {
      std::string fij_name = key->GetName();
      std::string error_term = fij_name.substr(0,fij_name.length()-4) ;
      std::vector<std::string>::iterator iter;
      iter = find(sk2_errors.begin(), sk2_errors.end(), error_term);
      std::vector<std::string>::iterator iter1;
      iter1 = find(sk1_errors.begin(), sk1_errors.end(), error_term);
      TH2F* th2 = (TH2F*)fijs_sk2->Get(fij_name.c_str());
      if(iter == sk2_errors.end() && error_term != "CC_nutau_xsec" && th2->Integral() > 0)
      // read fijs for each error, and use the untrivial fijs. DO NOT constrain CC_nutau_xsec.
      {
        sk2_errors.push_back(error_term);
        // Define the variables for systematic errors, 1 means shifting the systematic error by 1 sigma
        // // The PDFs of shifting 1 sigma is built with Osc3++. 
        TString constraint_expr = TString::Format("RooGaussian::%s_sys(%s[-3,3], 0, 1.0)", 
                error_term.c_str(), error_term.c_str());//Gaussian constraint.
        if(iter1 == sk1_errors.end()) tau_fit->factory(constraint_expr);
        TString constraint_name = TString::Format("%s_sys", error_term.c_str());
        parts_pdfII.add(*tau_fit->pdf(constraint_name));
      }
  }
 
  //Create RooHistPDF for each systematic error. Divide each term to positive and
  //and negative parts.
  //SKII positive
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_posII;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_posII;
  for(unsigned int i = 0; i < sk2_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk2_errors[i] + "_pos";
    TH2F* hist_temp = (TH2F*)fijs_sk2->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_posII[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_posII[i].reset(new RooFormulaVar(coeff_name.c_str(), "@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk2_errors[i].c_str())) ));
    if(hist_temp->Integral() > 0)// Only keep the non-trivial PDFs.
    {
        pdf_listII.add(*hist_pdfs_posII[i]);
        coeff_listII.add(*coeffs_posII[i]);
    }
  }
  //SKII negative
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_negII;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_negII;
  for(unsigned int i = 0; i < sk2_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk2_errors[i] + "_neg";
    TH2F* hist_temp = (TH2F*)fijs_sk2->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_negII[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_negII[i].reset(new RooFormulaVar(coeff_name.c_str(), "-1*@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk2_errors[i].c_str())) ));
    if(hist_temp->Integral()>0)
    {
        pdf_listII.add(*hist_pdfs_negII[i]);
        coeff_listII.add(*coeffs_negII[i]);
    }
  }
  //stitch the parts of PDF to the final PDF.
  //Create RooHistPDF for each systematic error. Divide each term to positive and
  RooAddPdf modelII ("modelII",  "signal+bkgd SKII",  pdf_listII, coeff_listII );
  RooProdPdf modelII_sys("modelII with sys errors", "modelII with sys errors", parts_pdfII);
 
  //SK2 MC study with systematic errors. Turned off by default.
  if(kFALSE)
  {
    RooRandom::randomGenerator()->SetSeed(0);
    TFile* mc_study = new TFile("mc_tau_sk2.root","recreate");
    TH1F* mc_xsec = new TH1F("xsec","xsec",40,0,4);
    TRandom3* r3 = new TRandom3();
    for(int i = 0; i < 1; i++)// Create MC sample for study.
    {
        RooDataSet  *mc_bkgII = bkgPdfII.generate(axisVariables, TMath::Nint(r3->PoissonD(1496)));
        RooDataSet  *mc_sigII = sigPdfII.generate(axisVariables, TMath::Nint(r3->PoissonD(32)));
        mc_sigII->append(*mc_bkgII);
        //Create likelihood function for each SK period.
        RooAbsReal* nll_sk2 = modelII.createNLL(*mc_sigII,RooFit::ExternalConstraints(modelII_sys),RooFit::Extended(kTRUE));
  
        //RooAddition nllCombined("nll_combined", "nll_combined", RooArgSet(*nll_sk1, *nll_sk2, *nll_sk3, *nll_sk4));
        RooMinuit minit(*nll_sk2);
        minit.migrad();
        mc_xsec->Fill(cc_nutau_xsec.getVal());
    }
    mc_study->Write();
    return;
  }
  
  // Read the file of histograms for sk3 systematic errors.
  TFile* fijs_sk3 = new TFile("./sys_pdf/error.sk3.root", "READ");

  // Read the names of systematic errors from fijs file, and store in a vector of string.
  // create Gaussian constraint for each systematic error.
  std::vector<std::string> sk3_errors;
  TIter next3(fijs_sk3->GetListOfKeys());
  RooArgSet parts_pdfIII("parts_pdfIII");
  while ((key = (TKey*)next3())) {
      std::string fij_name = key->GetName();
      std::string error_term = fij_name.substr(0,fij_name.length()-4) ;
      std::vector<std::string>::iterator iter;
      iter = find(sk3_errors.begin(), sk3_errors.end(), error_term);
      std::vector<std::string>::iterator iter1;
      iter1 = find(sk1_errors.begin(), sk1_errors.end(), error_term);
      TH2F* th2 = (TH2F*)fijs_sk3->Get(fij_name.c_str());
      if(iter == sk3_errors.end() &&  error_term != "CC_nutau_xsec" && th2->Integral() > 0)
      // read fijs for each error, and use the untrivial fijs. DO NOT constrain CC_nutau_xsec.
      {
        sk3_errors.push_back(error_term);
        // Define the variables for systematic errors, 1 means shifting the systematic error by 1 sigma
        // // The PDFs of shifting 1 sigma is built with Osc3++. 
        TString constraint_expr = TString::Format("RooGaussian::%s_sys(%s[-3,3], 0, 1.0)", 
                error_term.c_str(), error_term.c_str());//Gaussian constraint.
        if(iter1 == sk1_errors.end()) tau_fit->factory(constraint_expr);
        TString constraint_name = TString::Format("%s_sys", error_term.c_str());
        parts_pdfIII.add(*tau_fit->pdf(constraint_name));
      }
  }
 
  //Create RooHistPDF for each systematic error. Divide each term to positive and
  //and negative parts.
  //SKIII positive
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_posIII;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_posIII;
  for(unsigned int i = 0; i < sk3_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk3_errors[i] + "_pos";
    TH2F* hist_temp = (TH2F*)fijs_sk3->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_posIII[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_posIII[i].reset(new RooFormulaVar(coeff_name.c_str(), "@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk3_errors[i].c_str())) ));
    if(hist_temp->Integral() > 0)// Only keep the non-trivial PDFs.
    {
        pdf_listIII.add(*hist_pdfs_posIII[i]);
        coeff_listIII.add(*coeffs_posIII[i]);
    }
  }
  //SKIII negative
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_negIII;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_negIII;
  for(unsigned int i = 0; i < sk3_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk3_errors[i] + "_neg";
    TH2F* hist_temp = (TH2F*)fijs_sk3->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_negIII[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_negIII[i].reset(new RooFormulaVar(coeff_name.c_str(), "-1*@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk3_errors[i].c_str())) ));
    if(hist_temp->Integral()>0)
    {
        pdf_listIII.add(*hist_pdfs_negIII[i]);
        coeff_listIII.add(*coeffs_negIII[i]);
    }
  }
  //stitch the parts of PDF to the final PDF.
  //Create RooHistPDF for each systematic error. Divide each term to positive and
  RooAddPdf modelIII ("modelIII",  "signal+bkgd SKIII",  pdf_listIII, coeff_listIII );
  RooProdPdf modelIII_sys("modelIII with sys errors", "modelIII with sys errors", parts_pdfIII);
 
  //SK3 MC study with systematic errors. Turned off by default.
  if(kFALSE)
  {
    RooRandom::randomGenerator()->SetSeed(0);
    TFile* mc_study = new TFile("mc_tau_sk3.root","recreate");
    TH1F* mc_xsec = new TH1F("xsec","xsec",40,0,4);
    TRandom3* r3 = new TRandom3();
    for(int i = 0; i < 1000; i++)// Create MC sample for study.
    {
        RooDataSet  *mc_bkgIII = bkgPdfIII.generate(axisVariables, TMath::Nint(r3->PoissonD(988)));
        RooDataSet  *mc_sigIII = sigPdfIII.generate(axisVariables, TMath::Nint(r3->PoissonD(19)));
        mc_sigIII->append(*mc_bkgIII);
        //Create likelihood function for each SK period.
        RooAbsReal* nll_sk3 = modelIII.createNLL(*mc_sigIII, RooFit::ExternalConstraints(modelIII_sys),RooFit::Extended(kTRUE));
  
        RooMinuit minit(*nll_sk3);
        minit.migrad();
        mc_xsec->Fill(cc_nutau_xsec.getVal());
    }
    mc_study->Write();
  }

  // Read the file of histograms for sk4 systematic errors.
  TFile* fijs_sk4 = new TFile("./sys_pdf/error.sk4.root", "READ");

  // Read the names of systematic errors from fijs file, and store in a vector of string.
  // create Gaussian constraint for each systematic error.
  std::vector<std::string> sk4_errors;
  TIter next4(fijs_sk4->GetListOfKeys());
  RooArgSet parts_pdfIV("parts_pdfIV");
  while ((key = (TKey*)next4())) {
      std::string fij_name = key->GetName();
      std::string error_term = fij_name.substr(0,fij_name.length()-4) ;
      std::vector<std::string>::iterator iter;
      iter = find(sk4_errors.begin(), sk4_errors.end(), error_term);
      std::vector<std::string>::iterator iter1;
      iter1 = find(sk1_errors.begin(), sk1_errors.end(), error_term);
      TH2F* th2 = (TH2F*)fijs_sk4->Get(fij_name.c_str());
      if(iter == sk4_errors.end() && error_term != "CC_nutau_xsec" && th2->Integral() > 0)
      // read fijs for each error, and use the untrivial fijs. DO NOT constrain CC_nutau_xsec.
      {
        sk4_errors.push_back(error_term);
        // Define the variables for systematic errors, 1 means shifting the systematic error by 1 sigma
        // // The PDFs of shifting 1 sigma is built with Osc3++. 
        TString constraint_expr = TString::Format("RooGaussian::%s_sys(%s[-3,3], 0, 1.0)", 
                error_term.c_str(), error_term.c_str());//Gaussian constraint.
        if(iter1 == sk1_errors.end()) tau_fit->factory(constraint_expr);
        TString constraint_name = TString::Format("%s_sys", error_term.c_str());
        parts_pdfIV.add(*tau_fit->pdf(constraint_name));
      }
  }
 
  //Create RooHistPDF for each systematic error. Divide each term to positive and
  //and negative parts.
  //SKIV positive
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_posIV;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_posIV;
  for(unsigned int i = 0; i < sk4_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk4_errors[i] + "_pos";
    TH2F* hist_temp = (TH2F*)fijs_sk4->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_posIV[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_posIV[i].reset(new RooFormulaVar(coeff_name.c_str(), "@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk4_errors[i].c_str())) ));
    if(hist_temp->Integral() > 0)// Only keep the non-trivial PDFs.
    {
        pdf_listIV.add(*hist_pdfs_posIV[i]);
        coeff_listIV.add(*coeffs_posIV[i]);
    }
  }
  //SKIV negative
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_negIV;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_negIV;
  for(unsigned int i = 0; i < sk4_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk4_errors[i] + "_neg";
    TH2F* hist_temp = (TH2F*)fijs_sk4->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_negIV[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_negIV[i].reset(new RooFormulaVar(coeff_name.c_str(), "-1*@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk4_errors[i].c_str())) ));
    if(hist_temp->Integral()>0)
    {
        pdf_listIV.add(*hist_pdfs_negIV[i]);
        coeff_listIV.add(*coeffs_negIV[i]);
    }
  }
  //stitch the parts of PDF to the final PDF.
  //Create RooHistPDF for each systematic error. Divide each term to positive and
  RooAddPdf modelIV ("modelIV",  "signal+bkgd SKIV",  pdf_listIV, coeff_listIV );
  RooProdPdf modelIV_sys("modelIV with sys errors", "modelIV with sys errors", parts_pdfIV);

  //SK4 MC study with systematic errors. Turned off by default.
  if(kFALSE)
  {
    RooRandom::randomGenerator()->SetSeed(0);
    TFile* mc_study = new TFile("mc_tau_sk4.root","recreate");
    TH1F* mc_xsec = new TH1F("xsec","xsec",40,0,4);
    TRandom3* r3 = new TRandom3();
    for(int i = 0; i < 1000; i++)// Create MC sample for study.
    {
        RooDataSet  *mc_bkgIV = bkgPdfIV.generate(axisVariables, TMath::Nint(r3->PoissonD(3343)));
        RooDataSet  *mc_sigIV = sigPdfIV.generate(axisVariables, TMath::Nint(r3->PoissonD(68)));
        mc_sigIV->append(*mc_bkgIV);
        //Create likelihood function for each SK period.
        RooAbsReal* nll_sk4 = modelIV.createNLL(*mc_sigIV, RooFit::ExternalConstraints(modelIV_sys),RooFit::Extended(kTRUE));
  
        RooMinuit minit(*nll_sk4);
        minit.migrad();
        mc_xsec->Fill(cc_nutau_xsec.getVal());
    }
    mc_study->Write();
  }

  //SK1-4 MC study with systematic errors. Turned off by default. Because a single fit takes 20mins, this part
  //need to be processed in parallel. The code has been written for job submission in condor. A paramter is read
  //from the main function, and attached to the output file name.
  if(kFALSE)
  {
    RooRandom::randomGenerator()->SetSeed(0);
    TString root_file = "./mc_extra/mc_tau_" ;
    root_file += job;
    root_file += ".root";
    TFile* mc_study = new TFile(root_file,"recreate");
    TH1F* mc_xsec = new TH1F("xsec","xsec",40,0,4);
    TRandom3* r3 = new TRandom3();
    r3->SetSeed(0);
    for(int i = 0; i < 5; i++)// Create MC sample for study.
    {
        RooDataSet  *mc_bkgI = bkgPdfI.generate(axisVariables, TMath::Nint(r3->PoissonD(2817)));
        RooDataSet  *mc_sigI = sigPdfI.generate(axisVariables, TMath::Nint(r3->PoissonD(56*1.)));
        mc_sigI->append(*mc_bkgI);
        RooDataSet  *mc_bkgII = bkgPdfII.generate(axisVariables, TMath::Nint(r3->PoissonD(1496)));
        RooDataSet  *mc_sigII = sigPdfII.generate(axisVariables, TMath::Nint(r3->PoissonD(32*1.)));
        mc_sigII->append(*mc_bkgII);
        RooDataSet  *mc_bkgIII = bkgPdfIII.generate(axisVariables, TMath::Nint(r3->PoissonD(988)));
        RooDataSet  *mc_sigIII = sigPdfIII.generate(axisVariables, TMath::Nint(r3->PoissonD(19*1.)));
        mc_sigIII->append(*mc_bkgIII);
        RooDataSet  *mc_bkgIV = bkgPdfIV.generate(axisVariables, TMath::Nint(r3->PoissonD(3343)));
        RooDataSet  *mc_sigIV = sigPdfIV.generate(axisVariables, TMath::Nint(r3->PoissonD(68*1.)));//89 for adding extra data.
        mc_sigIV->append(*mc_bkgIV);
        //Create likelihood function for each SK period.
        RooAbsReal* nll_sk1 = modelI.createNLL(*mc_sigI,  RooFit::ExternalConstraints(modelI_sys),RooFit::Extended(kTRUE));
        RooAbsReal* nll_sk2 = modelII.createNLL(*mc_sigII,  RooFit::ExternalConstraints(modelII_sys),RooFit::Extended(kTRUE));
        RooAbsReal* nll_sk3 = modelIII.createNLL(*mc_sigIII,  RooFit::ExternalConstraints(modelIII_sys),RooFit::Extended(kTRUE));
        RooAbsReal* nll_sk4 = modelIV.createNLL(*mc_sigIV,  RooFit::ExternalConstraints(modelIV_sys),RooFit::Extended(kTRUE));
        RooAddition nllCombined("nll_combined", "nll_combined", RooArgSet(*nll_sk1, *nll_sk2, *nll_sk3, *nll_sk4));
        RooMinuit minit(nllCombined);
  
        minit.migrad();
        mc_xsec->Fill(cc_nutau_xsec.getVal());
    }
    mc_study->Write();
  }

  //Fit the data against PDFs.
  if(kTRUE)
  {
      //Create likelihood function for each SK period.
      RooAbsReal* nll_sk1 = modelI.createNLL(dataSetI,  RooFit::ExternalConstraints(modelI_sys),RooFit::Extended(kTRUE));
      RooAbsReal* nll_sk2 = modelII.createNLL(dataSetII,  RooFit::ExternalConstraints(modelII_sys),RooFit::Extended(kTRUE));
      RooAbsReal* nll_sk3 = modelIII.createNLL(dataSetIII,  RooFit::ExternalConstraints(modelIII_sys),RooFit::Extended(kTRUE));
      RooAbsReal* nll_sk4 = modelIV.createNLL(dataSetIV,  RooFit::ExternalConstraints(modelIV_sys),RooFit::Extended(kTRUE));
      RooAddition nllCombined("nll_combined", "nll_combined", RooArgSet(*nll_sk1, *nll_sk2, *nll_sk3, *nll_sk4));
      RooMinuit minit(nllCombined);
  
      minit.migrad();
      //cc_nutau_xsec.Print();
  }
}

int main(int argc, char** argv) {
    std::cout << argc << argv[1] << std::endl;// Read paramters from command line.
    joint_fit_sys(argv[1]);// The first paramter is passed to the function to specify the file name.
}
