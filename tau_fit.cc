#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>
#include <memory>
#include <chrono>
#include <unistd.h>
#include <fcntl.h>
#include <vector>
#include <string>

#include <tclap/CmdLine.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <TCanvas.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>
#include <TParameter.h>
#include <TKey.h>

#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooKeysPdf.h>
#include <RooRandom.h>
#include <RooBinning.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooNDKeysPdf.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooProdPdf.h>
#include <RooGaussian.h>

void joint_fit_sys()
{

  RooWorkspace* tau_fit= new RooWorkspace("tau_fit");

  // Define variables used in the fit
  RooRealVar cosTheta      ("cosTheta",    "Zenith of the Cosine Angle", -1.0, 1.0);
  RooRealVar NN_output     ("NN_output",   "NeuralNet Output",           -0.1, 1.1);
  RooRealVar NN_selected   ("NN_selected", "NeuralNet selected",          0.0, 1.0);
  RooRealVar eventWeight   ("eventWeight", "Event Weight",                0.0, 10.0);
  RooArgSet  chainVars     (cosTheta, NN_output, NN_selected, eventWeight);
  RooArgList axisVariables (cosTheta, NN_output);

  // Define Categories
  RooCategory dataPeriod("dataPeriod","SK Run Period");
  dataPeriod.defineType("SK1");

  TString objDir = "./data/" ;
  // Define variable used in the output tree

  // Open SKI and SKII and SKIII files
  TFile *skIfile   = new TFile((objDir +"SK-I.root"  ).Data(),  "READ");
  // Grab the trees and Histograms we need.
  //SK1
  TTree *SKIChain              = (TTree*)skIfile->Get("dataHistoTMVAOutputTree");
  //TTree *SKIBkgChain           = (TTree*)skIfile->Get("bkgHistoTMVAOutputTree");
  //TTree *SKISigChain           = (TTree*)skIfile->Get("tauHistoTMVAOutputTree");
  TH2F  *tauHisto2DZenithSKI   = (TH2F*) skIfile->Get("tauHistoZenith2D");
  TH2F  *bkgHisto2DZenithSKI   = (TH2F*) skIfile->Get("bkgHistoZenith2D");
/*
  RooDataSet RoofitBkgSetI("bkgTreeI", "Bkg ", chainVars, RooFit::Import(*SKIBkgChain), 
			   RooFit::Cut("NN_selected==1"), RooFit::WeightVar("eventWeight"));
  RooNDKeysPdf bkgKeysPdfI("bkgpdfI",     "Bkg I",   axisVariables, RoofitBkgSetI ,"m");
  RooDataSet RoofitSigSetI("sigTreeI", "Sig ", chainVars, RooFit::Import(*SKISigChain), 
			   RooFit::Cut("NN_selected==1"), RooFit::WeightVar("eventWeight"));
  RooNDKeysPdf sigKeysPdfI("sigpdfI", "Sig I", axisVariables, RoofitSigSetI, "m");
  */
  RooDataSet  dataSetI  ("dset", "Data", chainVars, RooFit::Import(*SKIChain), 
			 RooFit::Cut("NN_selected==1"));

  RooDataHist bkgHistoI("bkgI", "Bkg Histo", 
                               axisVariables, bkgHisto2DZenithSKI);
  RooDataHist sigHistoI("sigI", "Sig Histo", 
                               axisVariables, tauHisto2DZenithSKI);

  // Set the expected number of signal/background events (these are constants)
  // These are not normalized by any fraction (signal, background or DIS).
  RooRealVar expbkgI("expbkgI", "expbkgI",  bkgHisto2DZenithSKI->Integral());
  RooRealVar expsigI("expsigI", "expsigI",  tauHisto2DZenithSKI->Integral());
  RooRealVar cc_nutau_xsec("cc_nutau_xsec", "cc_nutau_xsec", -3, 3);
  RooFormulaVar tau_normI("tau_normI", "expsigI*cc_nutau_xsec",  RooArgList(expsigI, cc_nutau_xsec));

  // Make base PDFs with background simulation(normal ATM neutrinos, except for tau), and tau neutrino
  // simulation. Both PDFs have to be normalized to livetime, because they are fixed in the final PDF,
  // only the systematic errors are varied in the fitting.
  RooHistPdf  bkgPdfI("bkgpdfI", "Background PDF ", 
                             axisVariables, bkgHistoI, 0);
  RooHistPdf  sigPdfI("sigpdfI", "Signal PDF ",     
                             axisVariables, sigHistoI, 0);

  // Create a list of PDFs and coefficients that will be used to construct the final PDF.
  RooArgList pdf_list("pdf_list");
  pdf_list.add(bkgPdfI);
  pdf_list.add(sigPdfI);
  RooArgList coeff_list("coeff_list");
  coeff_list.add(expbkgI);
  coeff_list.add(tau_normI);

  // Make PDFs for systematic errors. The PDF is built with the change in each bin after shifting
  // systematic error by one sigma. Here, we assume the change corresponding to the systematic error
  // is linear. Each systematic error consists of two histograms, positive and negative part. This is
  // due to the fact that PDF has to be positive. The PDF of systematic errors is added with a varied
  // parameter to the base PDF. 
  // PDF_final = PDF_bkg + PDF_tau + Sum_i(var_sys_i*PDF_pos - var_sys_i*PDF_neg)
  //
  //systematic errors that are used in the SK Osc3++ analysis. skip t2k related errors. Each error is
  //constrained with a Gaussian function centered at 0 and sigma=1. CC_nutau_xsec is not constrained.
  //
  
  // SK1 systematic errors, SK2 3&4 share part of the terms.
  /*
   'abs_norm_E_lt_1GeV', 'nu_nubar_ratio_E_lt_1GeV', 'nu_nubar_ratio_E_gt_10GeV', 'nuebar_nue_E_lt_1GeV', 'nuebar_nue_E_gt_10GeV', 'numubar_numu_E_lt_1GeV', 'numubar_numu_E_gt_10GeV', 'up_down_ratio', 'horizontal_vertical_ratio', 'K_pi_ratio', 'nu_path', 'abs_norm_E_gt_1GeV', 'relative_norm_FC_multi_GeV', 'axial_mass_QE_and_1pi', 'CCQE_xsec_ratio', 'single_meson_xsec', 'DIS_model_difference', 'DIS_xsec', 'coherent_pi_xsec', 'NC_CC_ratio', 'CC_nutau_xsec', 'FC_reduction_sk1', 'FC_PC_separation_sk1', 'hadron_simulation', 'non_nue_bg_elike_sk1', 'non_nu_bg_mulike_sk1', 'ring_separation_sk1', 'pid_1ring_sk1', 'pid_multi_ring_sk1', 'energy_calibration_sk1', 'up_down_energy_calibration_sk1', 'non_nue_bg_mg_1ring_elike_sk1', 'non_nue_bg_mg_multiring_elike_sk1', 'mgmre_nue_nuebar_separation_sk1', 'solar_activity_sk1', 'polfit_1_ring_pi0_sk1', 'decay_e_tagging_sk1', 'CCQE_nu_nubar_ratio', 'CCQE_numu_nue_ratio', 'pi0_qpi_ratio', 'nubar_nu_1pi_ratio', 'nu_nubar_ratio_1_E_10GeV', 'nuebar_nue_1_E_10GeV', 'numubar_numu_1_E_10GeV', 'pi_decay_tagging_error', 'fiducial_volume_sk1', 'dis_q2_high_W', 'dis_q2_low_W', 'mgmre_other_separation_sk1', 'theta13_error', 'solDm2_error', 'theta12_error', 't2k_m23_error', 't2k_s23_error', 'matter_effect_error', 'neut_axial_mass'
   */
  
  // Define the variables for systematic errors, 1 means shifting the systematic error by 1 sigma
  // The PDFs of shifting 1 sigma is built with Osc3++. 
  // these variables should be independent of SK periods.
  
  // Read the file of histograms for systematic errors.
  TFile* fijs_fc = new TFile("./sys_pdf/error.sk1.root", "READ");

  // Read the names of systematic errors from fijs file, and store in a vector of string.
  std::vector<std::string> sk1_errors;
  TIter next(fijs_fc->GetListOfKeys());
  TKey* key;
  RooArgSet parts_final_pdf("parts_final_pdf");
  
  while ((key = (TKey*)next())) {
      std::string fij_name = key->GetName();
      std::string error_term = fij_name.substr(0,fij_name.length()-4) ;
      std::vector<std::string>::iterator iter;
      iter = find(sk1_errors.begin(), sk1_errors.end(), error_term);
      TH2F* th2 = (TH2F*)fijs_fc->Get(fij_name.c_str());
      if(iter == sk1_errors.end() && error_term != "CC_nutau_xsec" && th2->Integral()>0)
      {
        sk1_errors.push_back(error_term);
        std::string expr = error_term + "[-3,3]";
        TString constraint_expr = TString::Format("RooGaussian::%s_sys(%s[-3,3], 0, 1.0)", 
                error_term.c_str(), error_term.c_str());
        tau_fit->factory(constraint_expr);
        TString constraint_name = TString::Format("%s_sys", error_term.c_str());
        parts_final_pdf.add(*tau_fit->pdf(constraint_name));
      }
  }
  
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_pos;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_pos;

  for(unsigned int i = 0; i < sk1_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk1_errors[i] + "_pos";
    TH2F* hist_temp = (TH2F*)fijs_fc->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_pos[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_pos[i].reset(new RooFormulaVar(coeff_name.c_str(), "@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk1_errors[i].c_str())) ));
    pdf_list.add(*hist_pdfs_pos[i]);
    coeff_list.add(*coeffs_pos[i]);
  }
  
  std::map<int, std::shared_ptr<RooHistPdf>> hist_pdfs_neg;
  std::map<int, std::shared_ptr<RooFormulaVar>> coeffs_neg;

  for(unsigned int i = 0; i < sk1_errors.size() ; i++) 
  {
    //build PDFs for each systematic error.
    string name = sk1_errors[i] + "_neg";
    TH2F* hist_temp = (TH2F*)fijs_fc->Get(name.c_str());
    std::string datahist_name = "hist_" + name;
    RooDataHist* datahist_temp = new RooDataHist(datahist_name.c_str(), datahist_name.c_str(), axisVariables, hist_temp);
    std::string pdf_name = "pdf_" + name;
    hist_pdfs_neg[i].reset(new RooHistPdf(pdf_name.c_str(), pdf_name.c_str(), axisVariables, *datahist_temp));
    std::string coeff_name = "coeff_" + name; 
    coeffs_neg[i].reset(new RooFormulaVar(coeff_name.c_str(), "-1*@0*@1", RooArgList(RooFit::RooConst(hist_temp->Integral()), *tau_fit->var(sk1_errors[i].c_str())) ));
    pdf_list.add(*hist_pdfs_neg[i]);
    coeff_list.add(*coeffs_neg[i]);
  }
  
  RooAddPdf model1I ("modelI",  "signal+bkgd SKI",  pdf_list, coeff_list );
  parts_final_pdf.add(model1I);
  RooProdPdf model_sys("model with sys errors", "model with sys errors", parts_final_pdf);
  
  RooDataSet  *mc_bkg = bkgPdfI.generate(axisVariables, 2817);
  RooDataSet  *mc_sig = sigPdfI.generate(axisVariables, 56);
  mc_sig->append(*mc_bkg);
  RooFitResult *r1 = model_sys.fitTo(*mc_sig, RooFit::Save(), RooFit::Strategy(1), RooFit::Minimizer("Minuit2", "Migrad"));
  r1->Print();
}

int main(int argc, char** argv) {
    joint_fit_sys();
}
