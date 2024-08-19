from ROOT import *
import numpy as np
import random

myfile = TFile("pi0_invarmass_phdef.root","RECREATE")
f = TFile("1pi0_6May.root")
ana_tree = f.Get("ana_tree")

#calibration factors
cf_ecal_mean = 125.
cf_ecal_median = 134.7
#cf_hcal_mean = 10.45
cf_hcal_mean = 9.83
cf_hcal_median = 10.59
cf_lw = 1.60
#cf_byslices = [13.1,10.2,9.6,9.6,9.3,9.8,8.7,13.9]
#cf_byslices500keV = [9.1,9.9,9.7,9.7,9.3,9.8,9.2,10.8]
#cf_byslices500keV = [9.25,10.27,9.81,9.83,9.74,9.57,8.99,9.91]
cf_byslices = [22., 11.1, 8.9, 8.4, 8.]

pi0_mass = 135.

ecal_mean = 129.51112899462882
ecal_min = pi0_mass - 30.
ecal_max = pi0_mass + 30.

hcal_mean = 126.76624206466326
hcal_min = pi0_mass - 50.
hcal_max = pi0_mass + 50.

eh_mean = 124.41244081563899
eh_min = pi0_mass - 50.
eh_max = pi0_mass + 50.

# lines
pi0_mass_ecal = TLine(135., 0, 135, 0.141)
pi0_mass_ecal.SetLineColor(2)
ecal_lb = TLine(105., 0, 105, 0.141)
ecal_lb.SetLineColor(12)
ecal_lb.SetLineStyle(2)
ecal_ub = TLine(165., 0, 165, 0.141)
ecal_ub.SetLineColor(12)
ecal_ub.SetLineStyle(2)

pi0_mass_hcal = TLine(135.,0,135,0.141)
pi0_mass_hcal.SetLineColor(2)
hcal_lb = TLine(85.,0,85,0.141)
hcal_lb.SetLineColor(12)
hcal_lb.SetLineStyle(9)
hcal_ub = TLine(185,0,185,0.141)
hcal_ub.SetLineColor(12)
hcal_ub.SetLineStyle(9)

pi0_mass_1e1h = TLine(135,0,135,0.095)
pi0_mass_1e1h.SetLineColor(2)
eh_lb = TLine(85,0,85,0.095)
eh_lb.SetLineColor(12)
eh_lb.SetLineStyle(2)
eh_ub = TLine(185,0,185,0.095)
eh_ub.SetLineColor(12)
eh_ub.SetLineStyle(2)

def cos_theta(px, py, pz, ph1, ph2):
    return (px[ph1[0]]*px[ph2[0]]+py[ph1[0]]*py[ph2[0]]+pz[ph1[0]]*pz[ph2[0]])/(np.sqrt(px[ph1[0]]*px[ph1[0]]+py[ph1[0]]*py[ph1[0]]+pz[ph1[0]]*pz[ph1[0]])*np.sqrt(px[ph2[0]]*px[ph2[0]]+py[ph2[0]]*py[ph2[0]]+pz[ph2[0]]*pz[ph2[0]]))

#def hcal_byslices_energy500keV(e,ph):
#    if(e[ph[0]]>=0 and e[ph[0]]<100):
#        return cf_byslices500keV[0]*e[ph[0]]
#    if(e[ph[0]]>=100 and e[ph[0]]<200):
#        return cf_byslices500keV[1]*e[ph[0]]
#    if(e[ph[0]]>=200 and e[ph[0]]<300):
#        return cf_byslices500keV[2]*e[ph[0]]
#    if(e[ph[0]]>=300 and e[ph[0]]<400):
#        return cf_byslices500keV[3]*e[ph[0]]
#    if(e[ph[0]]>=400 and e[ph[0]]<500):
#        return cf_byslices500keV[4]*e[ph[0]]
#    if(e[ph[0]]>=500 and e[ph[0]]<600):
#        return cf_byslices500keV[5]*e[ph[0]]
#    if(e[ph[0]]>=600 and e[ph[0]]<800):
#        return cf_byslices500keV[6]*e[ph[0]]
#    if(e[ph[0]]>=800 and e[ph[0]]<1400):
#        return cf_byslices500keV[7]*e[ph[0]]

def hcal_byslices_energy500keV(e,ph):
    if(e[ph[0]]<1.5):
        return 25.*e[ph[0]]
    if(e[ph[0]]>=1.5 and e[ph[0]]<10):
        return cf_byslices[0]*e[ph[0]]
    if(e[ph[0]]>=10 and e[ph[0]]<20):
        return cf_byslices[1]*e[ph[0]]
    if(e[ph[0]]>=20 and e[ph[0]]<30):
        return cf_byslices[2]*e[ph[0]]
    if(e[ph[0]]>=30 and e[ph[0]]<40):
        return cf_byslices[3]*e[ph[0]]
    if(e[ph[0]]>=40):
        return cf_byslices[4]*e[ph[0]]

# total reco energy
def total_reco_energy_per_ph(ee, eh, ph):
    return cf_lw*ee[ph[0]]+hcal_byslices_energy500keV(eh,ph)

#hcal    
def invar_mass_hcal_byslices500keV(e,px,py,pz,ph1,ph2):
    return np.sqrt(2.*hcal_byslices_energy500keV(e,ph1)*hcal_byslices_energy500keV(e,ph2)*(1.-cos_theta(px,py,pz,ph1,ph2)))

#def invar_mass_hcal_calibrated_mean(e,px,py,pz,ph1,ph2):
#    return np.sqrt(2.*cf_hcal_mean*cf_hcal_mean*e[ph1[0]]*e[ph2[0]]*(1.-cos_theta(px,py,pz,ph1,ph2)))

#1e1h cases
#def invar_mass_1e1h_calibrated_mean_lw(ee,eh,px,py,pz,ph1,ph2):
#    return np.sqrt(2.*cf_lw*ee[ph1[0]]*cf_hcal_mean*eh[ph2[0]]*(1.-cos_theta(px,py,pz,ph1,ph2)))

def invar_mass_1e1h_byslices500keV(ee,eh,px,py,pz,ph1,ph2):
    return np.sqrt(2.*cf_lw*ee[ph1[0]]*hcal_byslices_energy500keV(eh,ph2)*(1.-cos_theta(px,py,pz,ph1,ph2)))

#ecal
def invar_mass_ecal_calibrated_lw(e,px,py,pz,ph1,ph2):
    return np.sqrt(2.*cf_lw*cf_lw*e[ph1[0]]*e[ph2[0]]*(1.-cos_theta(px,py,pz,ph1,ph2)))

# total reco invariant mass
def invar_mass_total_e(ee,eh,px,py,pz,ph1,ph2):
    return np.sqrt(2.*total_reco_energy_per_ph(ee,eh,ph1)*total_reco_energy_per_ph(ee,eh,ph2)*(1.-cos_theta(px,py,pz,ph1,ph2)))

#historgram definitions
h_invarmass_hcal_comb = TH1F("h_invarmass_hcal_comb","'Reco' #pi^{0} mass with bleed (hcal-hcal)",40,0,400) ##takes "visible/sim" energy
h_invarmass_hcal_tot = TH1F("h_invarmass_hcal_tot","'Reco' #pi^{0} mass with bleed (hcal-hcal)",40,0,400) ##takes "visible/sim" energy

h_invarmass_1e1h_comb_1 = TH1F("h_invarmass_1e1h_comb_1","'Reco' #pi^{0} mass with bleed (ecal-hcal)",40,0,400)
h_invarmass_1e1h_comb_2 = TH1F("h_invarmass_1e1h_comb_2","#pi^{0} invariant mass 1e1h case",40,0,400)
h_invarmass_1e1h_comb_tot = TH1F("h_invarmass_1e1h_comb_tot","'Reco' #pi^{0} mass with bleed (ecal-hcal)",40,0,400)
h_invarmass_1e1h_tot_1 = TH1F("h_invarmass_1e1h_tot_1","#pi^{0} invariant mass 1e1h case",40,0,400)
h_invarmass_1e1h_tot_2 = TH1F("h_invarmass_1e1h_tot_2","#pi^{0} invariant mass 1e1h case",40,0,400)
h_invarmass_1e1h_tot_tot = TH1F("h_invarmass_1e1h_tot_tot","#pi^{0} invariant mass 1e1h case",40,0,400)

h_invarmass_ecal_lw = TH1F("h_invarmass_ecal_lw","#pi^{0} invariant mass ecal case (layer-weight)",40,0,400)
h_invarmass_ecal_lw_tot = TH1F("h_invarmass_ecal_lw_tot","'Reco' #pi^{0} mass with bleed (ecal-ecal)",40,0,400)

h_n_ecal = TH1F("h_n_ecal","",40,0,400)
h_n_1e1h1 = TH1F("h_n_1e1h1","",40,0,400)
h_n_1e1h2 = TH1F("h_n_1e1h2","",40,0,400)
h_n_hcal = TH1F("h_n_hcal","",40,0,400)

for entryNum in range(0,ana_tree.GetEntries()):
    ana_tree.GetEntry(entryNum)
    ph1_idx = getattr(ana_tree,"pi0_photon1_idx")
    ph2_idx = getattr(ana_tree,"pi0_photon2_idx")
    ph1_det = getattr(ana_tree,"pi0_photon1_det")
    ph2_det = getattr(ana_tree,"pi0_photon2_det")
    ph_e = getattr(ana_tree,"sim_p_ecal_e")
    ph_h = getattr(ana_tree,"sim_p_hcal_e")
    ph_px = getattr(ana_tree,"sim_p_px")
    ph_py = getattr(ana_tree,"sim_p_py")
    ph_pz = getattr(ana_tree,"sim_p_pz")
    ph_e_true = getattr(ana_tree,"sim_p_e")
    ph_e_lw = getattr(ana_tree,"sim_p_ecal_e_lw")
    if(len(ph_e)>ph1_idx[0] and len(ph_h)>ph1_idx[0] and len(ph_e)>ph2_idx[0] and len(ph_h)>ph2_idx[0]):
        
        # ecal
        if(ph1_det[0] == 1 and ph2_det[0] == 1):
            h_invarmass_ecal_lw.Fill(invar_mass_ecal_calibrated_lw(ph_e_lw,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))
            h_invarmass_ecal_lw_tot.Fill(invar_mass_total_e(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))
            if(invar_mass_ecal_calibrated_lw(ph_e_lw,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<ecal_max and invar_mass_ecal_calibrated_lw(ph_e_lw,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>ecal_min):
                h_n_ecal.Fill(invar_mass_ecal_calibrated_lw(ph_e_lw,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))

        # hcal
        if(ph1_det[0] == 2 and ph2_det[0] == 2):
            h_invarmass_hcal_comb.Fill(invar_mass_hcal_byslices500keV(ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))
            h_invarmass_hcal_tot.Fill(invar_mass_total_e(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))
            if(invar_mass_hcal_byslices500keV(ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<hcal_max and invar_mass_hcal_byslices500keV(ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>hcal_min):
                h_n_hcal.Fill(invar_mass_hcal_byslices500keV(ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))

        # 1e1h 1
        if(ph1_det[0] == 1 and ph2_det[0] == 2):
            h_invarmass_1e1h_comb_1.Fill(invar_mass_1e1h_byslices500keV(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))
            h_invarmass_1e1h_tot_1.Fill(invar_mass_total_e(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))
            if(invar_mass_1e1h_byslices500keV(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<hcal_max and invar_mass_1e1h_byslices500keV(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>hcal_min):
                h_n_1e1h1.Fill(invar_mass_1e1h_byslices500keV(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))

        # 1e1h 2
        if(ph1_det[0] ==2 and ph2_det[0] == 1):
            h_invarmass_1e1h_comb_2.Fill(invar_mass_1e1h_byslices500keV(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph2_idx,ph1_idx))
            h_invarmass_1e1h_tot_2.Fill(invar_mass_total_e(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph1_idx,ph2_idx))
            if(invar_mass_1e1h_byslices500keV(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph2_idx,ph1_idx)<hcal_max and invar_mass_1e1h_byslices500keV(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph2_idx,ph1_idx)>hcal_min):
                h_n_1e1h2.Fill(invar_mass_1e1h_byslices500keV(ph_e_lw,ph_h,ph_px,ph_py,ph_pz,ph2_idx,ph1_idx))

# sum histos
h_invarmass_1e1h_comb_tot = h_invarmass_1e1h_comb_1.Clone()
h_invarmass_1e1h_comb_tot.Add(h_invarmass_1e1h_comb_2)
h_invarmass_1e1h_tot_tot = h_invarmass_1e1h_tot_1.Clone()
h_invarmass_1e1h_tot_tot.Add(h_invarmass_1e1h_comb_2)

print(h_invarmass_hcal_comb.GetEntries())
print(h_invarmass_1e1h_comb_tot.GetEntries())
print(h_invarmass_ecal_lw.GetEntries())

print("pass selection: ecal, hcal, 1e1h")
print(h_n_ecal.GetEntries())
print(h_n_hcal.GetEntries())
print(h_n_1e1h1.GetEntries())
print(h_n_1e1h2.GetEntries())

h_invarmass_hcal_comb.Scale(1./h_invarmass_hcal_comb.Integral())
h_invarmass_1e1h_comb_tot.Scale(1./h_invarmass_1e1h_tot_tot.Integral())
h_invarmass_ecal_lw.Scale(1./h_invarmass_ecal_lw.Integral())

c1 = TCanvas("c1")
h_invarmass_hcal_comb.GetXaxis().SetTitle("m_{#gamma#gamma} [MeV]")
h_invarmass_hcal_comb.GetYaxis().SetTitle("A.U.")
h_invarmass_hcal_comb.SetStats(0)
h_invarmass_hcal_comb.SetLineColor(1)
#h_invarmass_hcal_tot.SetLineColor(1)
#h_invarmass_hcal_comb.Draw("hist")
h_invarmass_hcal_comb.Draw("hist")
pi0_mass_hcal.Draw("hist same")
hcal_lb.Draw("hist same")
hcal_ub.Draw("hist same")
#l1 = TLegend(0.6,0.6,0.95,0.85)
#l1.AddEntry(h_invarmass_hcal_comb,"only hcal energy")
#l1.AddEntry(h_invarmass_hcal_tot,"total reco energy")
#l1.Draw("hist same")
print(h_invarmass_hcal_comb.GetMean())
#print(h_invarmass_hcal_tot.GetMean())
#c1.SaveAs("invarmass_hcal.png")
myfile.WriteObject(c1,"pi0_hh_invarmass")

c2 = TCanvas("c2")
h_invarmass_1e1h_comb_tot.GetXaxis().SetTitle("m_{#gamma#gamma} [MeV]")
h_invarmass_1e1h_comb_tot.GetYaxis().SetTitle("A.U.")
h_invarmass_1e1h_comb_tot.SetStats(0)
h_invarmass_1e1h_comb_tot.SetLineColor(1)
#h_invarmass_1e1h_tot_tot.SetLineColor(1)
h_invarmass_1e1h_comb_tot.Draw("hist")
pi0_mass_1e1h.Draw("hist same")
eh_lb.Draw("hist same")
eh_ub.Draw("hist same")
#l3 = TLegend(0.6,0.6,0.95,0.85)
#l3.AddEntry(h_invarmass_1e1h_comb_tot,"only primary det energy")
#l3.AddEntry(h_invarmass_1e1h_tot_tot,"total reco energy")
#l3.Draw("hist same")
print(h_invarmass_1e1h_comb_tot.GetMean())
#print(h_invarmass_1e1h_tot_tot.GetMean())
#c2.SaveAs("invarmass_1e1h.png")
myfile.WriteObject(c2,"pi0_eh_invarmass")

c3 = TCanvas("c3")
h_invarmass_ecal_lw.GetXaxis().SetTitle("m_{#gamma#gamma} [MeV]")
h_invarmass_ecal_lw.GetYaxis().SetTitle("A.U.")
h_invarmass_ecal_lw.SetStats(0)
#h_invarmass_ecal_lw.SetLineColor(9)
h_invarmass_ecal_lw.SetLineColor(1)
h_invarmass_ecal_lw.Draw("hist")
pi0_mass_ecal.Draw("hist same")
ecal_lb.Draw("hist same")
ecal_ub.Draw("hist same")
#l2 = TLegend(0.6,0.6,0.95,0.85)
#l2.AddEntry(h_invarmass_ecal_lw,"only ecal energy")
#l2.AddEntry(h_invarmass_ecal_lw_tot,"total reco energy")
#l2.Draw("hist same")
print(h_invarmass_ecal_lw.GetMean())
#print(h_invarmass_ecal_lw_tot.GetMean())
#c3.SaveAs("invarmass_ecal.png")
myfile.WriteObject(c3,"pi0_ee_invarmass")

c4 = TCanvas("c4")
h_invarmass_1e1h_comb_tot.SetLineColor(4)
h_invarmass_hcal_comb.SetLineColor(2)
h_invarmass_ecal_lw.Draw("hist")
h_invarmass_hcal_comb.Draw("hist same")
h_invarmass_1e1h_comb_tot.Draw("hist same")
pi0_mass_ecal.Draw("hist same")
ecal_lb.Draw("hist same")
ecal_ub.Draw("hist same")
hcal_lb.Draw("hist same")
hcal_ub.Draw("hist same")
l4 = TLegend(0.6,0.3,0.9,0.6)
l4.SetBorderSize(0)
l4.AddEntry(h_invarmass_ecal_lw,"ecal-ecal")
l4.AddEntry(h_invarmass_hcal_comb,"hcal-hcal")
l4.AddEntry(h_invarmass_1e1h_comb_tot,"ecal-hcal")
l4.AddEntry(ecal_lb,"ecal-ecal event sel")
l4.AddEntry(hcal_lb,"hcal event sel")
l4.Draw("hist same")
c4.SaveAs("pi0_invarmass_phdef.png")

myfile.WriteObject(c4,"pi0_invarmass_phdef")
