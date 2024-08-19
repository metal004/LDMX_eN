from ROOT import *
import numpy as np
import random

myfile = TFile("hists.root","RECREATE")
f = TFile("output_file_sidehcaledep_1M.root","UPDATE")
ana_tree = f.Get("ana_tree")

#calibration factors
cf_lw = 1.60
cf_byslices = [50,21.85,10.92,8.95,8.25,8.05]

# event selection
hcal_min = 134.97 - 50.
hcal_max = 134.97 + 50.

# lines
pi0_mass_ecal = TLine(135., 0, 135, 0.1385)
pi0_mass_ecal.SetLineColor(2)
ecal_lb = TLine(105., 0, 105, 0.1385)
ecal_lb.SetLineColor(12)
ecal_lb.SetLineStyle(2)
ecal_ub = TLine(165., 0, 165, 0.1385)
ecal_ub.SetLineColor(12)
ecal_ub.SetLineStyle(2)

pi0_mass_hcal = TLine(135.,0,135,0.102)
pi0_mass_hcal.SetLineColor(2)
hcal_lb = TLine(85.,0,85,0.102)
hcal_lb.SetLineColor(12)
hcal_lb.SetLineStyle(2)
hcal_ub = TLine(185,0,185,0.102)
hcal_ub.SetLineColor(12)
hcal_ub.SetLineStyle(2)

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

def ecal_energy(e,ph):
    return cf_lw*e[ph[0]]

def hcal_energy(e,ph):
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

def invarmass(e1,e2,px,py,pz,ph1,ph2):
    return np.sqrt(2.*e1*e2*(1-cos_theta(px,py,pz,ph1,ph2)))

def hcal_theta_max(thz,ph1,ph2):
    if(thz[ph1[0]]>thz[ph2[0]]):
        return thz[ph1[0]]
    else:
        return thz[ph2[0]]

def sidehcal_edep(p):
    total_edep = hcal_top_e[p] + hcal_bottom_e[p] + hcal_left_e[p] + hcal_right_e[p]
    return total_edep

# charged hadron multiplicity
h_prot = TH1D("h_prot","",10,0,10)
h_pi = TH1D("h_pi","",10,0,10)
h_ch = TH1D("h_ch","",10,0,10)
h_ch_dif = TH1D("h_ch_dif","",10,0,10)

# side hcal reduction. [40,80] is full extent; 75% is [40,70]; 80% is [40,80.4]; 50% is [40,60]
h_prot_75 = TH1D("h_prot_75","",10,0,10)
h_pi_75 = TH1D("h_pi_75","",10,0,10)
h_prot_66 = TH1D("h_prot_66","",10,0,10)
h_pi_66 = TH1D("h_pi_66","",10,0,10)
h_prot_50 = TH1D("h_prot_50","",10,0,10)
h_pi_50 = TH1D("h_pi_50","",10,0,10)

h_ch_90 = TH1D("h_ch_90","",10,0,10)
h_ch_80 = TH1D("h_ch_80","",10,0,10)
h_ch_70 = TH1D("h_ch_70","",10,0,10)
h_ch_60 = TH1D("h_ch_60","",10,0,10)
h_ch_50 = TH1D("h_ch_50","",10,0,10)
h_ch_dif_90 = TH1D("h_ch_dif_90","",10,0,10)
h_ch_dif_80 = TH1D("h_ch_dif_80","",10,0,10)
h_ch_dif_70 = TH1D("h_ch_dif_70","",10,0,10)
h_ch_dif_60 = TH1D("h_ch_dif_60","",10,0,10)
h_ch_dif_50 = TH1D("h_ch_dif_50","",10,0,10)

for entryNum in range(0,ana_tree.GetEntries()):
    ana_tree.GetEntry(entryNum)
    n_sim_p = getattr(ana_tree,"n_sim_p")
    p_pdg = getattr(ana_tree,"sim_p_pdg")
    pi0_idx = getattr(ana_tree,"pi0_idx")
    ph1_idx = getattr(ana_tree,"pi0_photon1_idx")
    ph2_idx = getattr(ana_tree,"pi0_photon2_idx")
    ph1_det = getattr(ana_tree,"pi0_photon1_det")
    ph2_det = getattr(ana_tree,"pi0_photon2_det")
    ph_e_true = getattr(ana_tree,"sim_p_e")
    ph_e_lw = getattr(ana_tree,"sim_p_ecal_e_lw")
    ph_h = getattr(ana_tree,"sim_p_hcal_e")
    ph_px = getattr(ana_tree,"sim_p_px")
    ph_py = getattr(ana_tree,"sim_p_py")
    ph_pz = getattr(ana_tree,"sim_p_pz")
    ph_thetaz = getattr(ana_tree,"sim_p_thetaz")
    hcal_top_e = getattr(ana_tree,"hcal_top_e")
    hcal_bottom_e = getattr(ana_tree,"hcal_bottom_e")
    hcal_left_e = getattr(ana_tree,"hcal_left_e")
    hcal_right_e = getattr(ana_tree,"hcal_right_e")
    hcal_back_e = getattr(ana_tree,"hcal_back_e")
    n_prot = getattr(ana_tree,"n_sim_prot")

    n_ch = 0
    n_ch_90 = 0
    n_ch_80 = 0
    n_ch_70 = 0
    n_ch_60 = 0
    n_ch_50 = 0

    # 100% is 80 deg, 50% is 60 deg; 90% is 76 deg, 80% is 72 deg, 70% is 68 deg, 60% is 64 deg
    for particle in range(n_sim_p):
        if((p_pdg[particle]==2212 or p_pdg[particle]==211 or p_pdg[particle]==-211) and sidehcal_edep(particle)>0):
            n_ch += 1
            if(ph_thetaz[particle]<76*3.1415/180):
                n_ch_90 += 1
            if(ph_thetaz[particle]<72*3.1415/180):
                n_ch_80 += 1
            if(ph_thetaz[particle]<68*3.1415/180):
                n_ch_70 += 1
            if(ph_thetaz[particle]<64*3.1415/180):
                n_ch_60 += 1
            if(ph_thetaz[particle]<60*3.1415/180):
                n_ch_50 += 1
        h_ch.Fill(n_ch)
        h_ch_90.Fill(n_ch_90)
        h_ch_80.Fill(n_ch_80)
        h_ch_70.Fill(n_ch_70)
        h_ch_60.Fill(n_ch_60)
        h_ch_50.Fill(n_ch_50)
        h_ch_dif_90.Fill(n_ch-n_ch_90)
        h_ch_dif_80.Fill(n_ch-n_ch_80)
        h_ch_dif_70.Fill(n_ch-n_ch_70)
        h_ch_dif_60.Fill(n_ch-n_ch_60)
        h_ch_dif_50.Fill(n_ch-n_ch_50)

c4 = TCanvas("c4")
#c4.SetLogy()
h_ch_50.SetStats(0)
h_ch_50.GetXaxis().SetTitle("multiplicity")
h_ch_50.GetYaxis().SetTitle("Count")
h_ch_50.Draw("hist")
h_ch.SetLineColor(1)
h_ch_90.SetLineColor(8)
h_ch_80.SetLineColor(87)
h_ch_70.SetLineColor(91)
h_ch_60.SetLineColor(96)
h_ch_50.SetLineColor(2)
h_ch_90.Draw("hist same")
h_ch_80.Draw("hist same")
h_ch_60.Draw("hist same")
h_ch.Draw("hist same")
h_ch_70.Draw("hist same")
l_ch = TLegend(0.6,0.6,0.85,0.85)
l_ch.SetBorderSize(0)
l_ch.AddEntry(h_ch,"Full side hcal")
l_ch.AddEntry(h_ch_90,"90% of side hcal")
l_ch.AddEntry(h_ch_80,"80% of side hcal")
l_ch.AddEntry(h_ch_70,"70% of side hcal")
l_ch.AddEntry(h_ch_60,"60% of side hcal")
l_ch.AddEntry(h_ch_50,"50% of side hcal")
l_ch.Draw("hist same")
c4.SaveAs("ch_mult_reduction.png")

c3 = TCanvas("c3")
c3.SetLogy()
h_ch_70.SetStats(0)
h_ch_70.GetXaxis().SetTitle("multiplicity")
h_ch_70.GetYaxis().SetTitle("Count")
h_ch_70.Draw("hist")
h_ch.SetLineColor(1)
h_ch_90.SetLineColor(8)
h_ch_80.SetLineColor(87)
h_ch_70.SetLineColor(91)
h_ch_60.SetLineColor(96)
h_ch_50.SetLineColor(2)
h_ch_90.Draw("hist same")
h_ch_80.Draw("hist same")
h_ch_60.Draw("hist same")
h_ch_50.Draw("hist same")
h_ch.Draw("hist same")
l_ch = TLegend(0.6,0.6,0.85,0.85)
l_ch.SetBorderSize(0)
l_ch.AddEntry(h_ch,"Full side hcal")
l_ch.AddEntry(h_ch_90,"90% of side hcal")
l_ch.AddEntry(h_ch_80,"80% of side hcal")
l_ch.AddEntry(h_ch_70,"70% of side hcal")
l_ch.AddEntry(h_ch_60,"60% of side hcal")
l_ch.AddEntry(h_ch_50,"50% of side hcal")
l_ch.Draw("hist same")
c3.SaveAs("ch_mult_reduction_log.png")

c5 = TCanvas("c5")
h_ch_dif_90.SetStats(0)
h_ch_dif_90.GetXaxis().SetTitle("multiplicity")
h_ch_dif_90.GetYaxis().SetTitle("Count")
h_ch_dif_90.Draw("hist")
h_ch_dif_90.SetLineColor(8)
h_ch_dif_80.SetLineColor(87)
h_ch_dif_70.SetLineColor(91)
h_ch_dif_60.SetLineColor(96)
h_ch_dif_50.SetLineColor(2)
h_ch_dif_70.Draw("hist same")
h_ch_dif_80.Draw("hist same")
h_ch_dif_60.Draw("hist same")
h_ch_dif_50.Draw("hist same")
l_ch_dif = TLegend(0.6,0.6,0.85,0.85)
l_ch_dif.SetBorderSize(0)
l_ch_dif.AddEntry(h_ch_dif_90,"90% of side hcal")
l_ch_dif.AddEntry(h_ch_dif_80,"80% of side hcal")
l_ch_dif.AddEntry(h_ch_dif_70,"70% of side hcal")
l_ch_dif.AddEntry(h_ch_dif_60,"60% of side hcal")
l_ch_dif.AddEntry(h_ch_dif_50,"50% of side hcal")
l_ch_dif.Draw("hist same")
c5.SaveAs("ch_mult_reduction_dif.png")
