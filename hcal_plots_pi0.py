from ROOT import *
import numpy as np
import random

myfile = TFile("hcal_pi0_plots_int_eff_hists.root","RECREATE")
f_pi0 = TFile("1pi0_30July.root","UPDATE")
ana_tree = f_pi0.Get("ana_tree")

#calibration factors
cf_lw = 1.60
cf_reduced = [50, 28.55, 15.25, 12.75, 11.55, 10.95]
cf_nominal = [50,21.85,10.92,8.95,8.25,8.05]

# event selection
hcal_min = 134.97 - 50.
hcal_max = 134.97 + 50.

ecal_min = 134.97 - 30.
ecal_max = 134.97 + 30.

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

def hcal_energy_nominal(e,ph):
    if(e[ph[0]]<1.5):
        return cf_nominal[0]*e[ph[0]]
    if(e[ph[0]]>=1.5 and e[ph[0]]<10):
        return cf_nominal[1]*e[ph[0]]
    if(e[ph[0]]>=10 and e[ph[0]]<20):
        return cf_nominal[2]*e[ph[0]]
    if(e[ph[0]]>=20 and e[ph[0]]<30):
        return cf_nominal[3]*e[ph[0]]
    if(e[ph[0]]>=30 and e[ph[0]]<40):
        return cf_nominal[4]*e[ph[0]]
    if(e[ph[0]]>=40):
        return cf_nominal[5]*e[ph[0]]

def hcal_energy_reduced(e,ph):
    if(e[ph[0]]<1.5):
        return cf_reduced[0]*e[ph[0]]
    if(e[ph[0]]>=1.5 and e[ph[0]]<10):
        return cf_reduced[1]*e[ph[0]]
    if(e[ph[0]]>=10 and e[ph[0]]<20):
        return cf_reduced[2]*e[ph[0]]
    if(e[ph[0]]>=20 and e[ph[0]]<30):
        return cf_reduced[3]*e[ph[0]]
    if(e[ph[0]]>=30 and e[ph[0]]<40):
        return cf_reduced[4]*e[ph[0]]
    if(e[ph[0]]>=40):
        return cf_reduced[5]*e[ph[0]]

def invarmass(e1,e2,px,py,pz,ph1,ph2):
    return np.sqrt(2.*e1*e2*(1-cos_theta(px,py,pz,ph1,ph2)))

def hcal_theta_max(thetaz,ph1,ph2):
    #print(f'event:{entryNum}')
    #print(len(ph1_idx))
    #print(ph1_idx[0])
    #print(ph2_idx[0])
    #print(len(thetaz))
    #print(thetaz[ph1_idx[0]])
    if(thetaz[ph1[0]]>thetaz[ph2[0]]):
        return thetaz[ph1[0]]
    else:
        return thetaz[ph2[0]]

def sidehcal_edep(ph):
    total_edep = hcal_top_e[ph] + hcal_bottom_e[ph] + hcal_left_e[ph] + hcal_right_e[ph]
    return total_edep

#def ph_sep(ph1_idx,ph2_idx,ph_px,ph_py,ph_pz,ph_thetaz):
#    phi_ph1 = np.arccos(ph_px[ph1_idx[0]]/np.sqrt(ph_px[ph1_idx[0]]*ph_px[ph1_idx[0]]+ph_py[ph1_idx[0]]*ph_py[ph1_idx[0]]))
#    phi_ph2 = np.arccos(ph_px[ph2_idx[0]]/np.sqrt(ph_px[ph2_idx[0]]*ph_px[ph2_idx[0]]+ph_py[ph2_idx[0]]*ph_py[ph2_idx[0]]))
#    r = 240.
#    x_ph1 = r*np.sin(ph_thetaz[ph1_idx[0]])*np.cos(phi_ph1)
#    y_ph1 = r*np.sin(ph_thetaz[ph1_idx[0]])*np.sin(phi_ph1)
#    x_ph2 = r*np.sin(ph_thetaz[ph2_idx[0]])*np.cos(phi_ph2)
#    y_ph2 = r*np.sin(ph_thetaz[ph2_idx[0]])*np.sin(phi_ph2)
#    dist = np.sqrt((x_ph2 - x_ph1)*(x_ph2 - x_ph1)+(y_ph2 - y_ph1)*(y_ph2 - y_ph1))
#    return dist

def ph_sep(ph1_idx,ph2_idx,px,py,pz,theta):
    z = 240.
    x1 = px[ph1_idx[0]] * z / pz[ph1_idx[0]]
    y1 = py[ph1_idx[0]] * z / pz[ph1_idx[0]]
    x2 = px[ph2_idx[0]] * z / pz[ph2_idx[0]]
    y2 = py[ph2_idx[0]] * z / pz[ph2_idx[0]]
    sep = np.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
    #print(x1)
    #print(y1)
    #print(x2)
    #print(y2)
    #print(sep)
    return sep

def ph_sep_x(ph1_idx,ph2_idx,px,py,pz,theta):
    z = 240.
    x1 = px[ph1_idx[0]] * z / pz[ph1_idx[0]]
    x2 = px[ph2_idx[0]] * z / pz[ph2_idx[0]]
    sep = np.abs(x2-x1)
    return sep

def ph_sep_y(ph1_idx,ph2_idx,px,py,pz,theta):
    z = 240.
    y1 = py[ph1_idx[0]] * z / pz[ph1_idx[0]]
    y2 = py[ph2_idx[0]] * z / pz[ph2_idx[0]]
    sep = np.abs(y2-y1)
    return sep

h_thzmax_eff = TH1D("h_thzmax_eff","",16,0,80)

ee_num = np.zeros(17)
ee_denom = np.zeros(17)
eh_num = np.zeros(17)
eh_denom = np.zeros(17)
he_num = np.zeros(17)
he_denom = np.zeros(17)
hh_num = np.zeros(17)
hh_denom = np.zeros(17)

h_thzmax_veto_eff = TH1D("h_thzmax_eff","",16,0,80)

eh_num_v = np.zeros(17)
eh_denom_v = np.zeros(17)
he_num_v = np.zeros(17)
he_denom_v = np.zeros(17)
hh_num_v = np.zeros(17)
hh_denom_v = np.zeros(17)

h_ph_sep_dis = TH1D("h_ph_sep_dis","",25,0,1000)
h_ph_sep_dis_top = TH1D("h_ph_sep_dis_top","",15,0,600)
h_ph_sep_dis_bottom = TH1D("h_ph_sep_dis_bottom","",15,0,600)
h_ph_sep_dis_left = TH1D("h_ph_sep_dis_left","",15,0,600)
h_ph_sep_dis_right = TH1D("h_ph_sep_dis_right","",15,0,600)

h_ph_sep_dis_x = TH1D("h_ph_sep_dis_x","",25,0,1000)
h_ph_sep_dis_y = TH1D("h_ph_sep_dis_y","",25,0,1000)

h_bottom1 = TH1D("h_bottom1","",100,0,4000)
h_bottom2 = TH1D("h_bottom2","",100,0,4000)

denom1 = 0
denomv1 = 0

#for entryNum in range(0,50):
for entryNum in range(0,ana_tree.GetEntries()):
    ana_tree.GetEntry(entryNum)
    n_sim_p = getattr(ana_tree,"n_sim_p")
    pi0_idx = getattr(ana_tree,"pi0_idx")
    ph1_idx = getattr(ana_tree,"pi0_photon1_idx")
    ph2_idx = getattr(ana_tree,"pi0_photon2_idx")
    ph1_det = getattr(ana_tree,"pi0_photon1_det")
    ph2_det = getattr(ana_tree,"pi0_photon2_det")
    ph_e_true = getattr(ana_tree,"sim_p_e")
    ph_e_lw = getattr(ana_tree,"sim_p_ecal_e_lw")
    ph_h = getattr(ana_tree,"sim_p_hcal_e")
    ph_h_top = getattr(ana_tree,"hcal_top_e")
    ph_h_bottom = getattr(ana_tree,"hcal_bottom_e")
    ph_h_left = getattr(ana_tree,"hcal_left_e")
    ph_h_right = getattr(ana_tree,"hcal_right_e")
    ph_px = getattr(ana_tree,"sim_p_px")
    ph_py = getattr(ana_tree,"sim_p_py")
    ph_pz = getattr(ana_tree,"sim_p_pz")
    ph_p = getattr(ana_tree,"sim_p_p")
    ph_thetaz = getattr(ana_tree,"sim_p_thetaz")

    #print(f'entryNum:{entryNum}')
    #print(f'ph1idx:{ph1_idx[0]}')

#    if(entryNum<3):
#        print(f'event:{entryNum}')
#        print(len(ph1_idx))
#        print(ph1_idx[0])
    #if(len(ph_e_lw)>ph1_idx[0] and len(ph_h)>ph1_idx[0] and len(ph_e_lw)>ph2_idx[0] and len(ph_h)>ph2_idx[0]):
    if(len(ph1_det)>0 and len(ph2_det)>0 and n_sim_p>=ph1_idx[0] and n_sim_p>=ph2_idx[0] and len(ph_thetaz)>ph1_idx[0] and len(ph_thetaz)>ph2_idx[0]):
#    if(len(ph1_det)>0 and len(ph2_det)>0 and n_sim_p>=ph1_idx[0] and n_sim_p>=ph2_idx[0]):

        # don't require that both photons go into hcal?
        for i in range(1,17):
           #if(hcal_theta_max(ph_thetaz,ph1_idx,ph2_idx)>=(5.*(i-1))*180/3.1415 and hcal_theta_max(ph_thetaz,ph1_idx,ph2_idx)*180/3.1415<=(5.*i)):
            #if(ph_thetaz[ph1_idx[0]]*180/3.1415<=5.*i and ph_thetaz[ph2_idx[0]]*180/3.1415<=5.*i):
            if(hcal_theta_max(ph_thetaz,ph1_idx,ph2_idx)*180/3.1415<=(5.*i)):
                if(ph1_det[0] == 1 and ph2_det[0] ==1):
                    ee_denom[i] += 1
                    if(invarmass(ecal_energy(ph_e_lw,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>ecal_min and invarmass(ecal_energy(ph_e_lw,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<=ecal_max):
                        ee_num[i] += 1
                if(ph1_det[0] == 1 and ph2_det[0] == 2):
                    eh_denom[i] += 1
                    if(invarmass(ecal_energy(ph_e_lw,ph1_idx),hcal_energy_nominal(ph_h,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>hcal_min and invarmass(ecal_energy(ph_e_lw,ph1_idx),hcal_energy_nominal(ph_h,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<=hcal_max):
                        eh_num[i] += 1
                if(ph1_det[0] == 2 and ph2_det[0] == 1):
                    he_denom[i] += 1
                    if(invarmass(hcal_energy_nominal(ph_h,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>hcal_min and invarmass(hcal_energy_nominal(ph_h,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<=hcal_max):
                        he_num[i] += 1
                if(ph1_det[0] == 2 and ph2_det[0] == 2):
                    hh_denom[i] += 1
                    if(invarmass(hcal_energy_nominal(ph_h,ph1_idx),hcal_energy_nominal(ph_h,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>hcal_min and invarmass(hcal_energy_nominal(ph_h,ph1_idx),hcal_energy_nominal(ph_h,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<=hcal_max):
                        hh_num[i] += 1
        
        for i in range(1,17):
            if(hcal_theta_max(ph_thetaz,ph1_idx,ph2_idx)*180/3.1415<=(5*i)):
            #if(hcal_theta_max(ph_thetaz,ph1_idx,ph2_idx)*180/3.1415>=40 and hcal_theta_max(ph_thetaz,ph1_idx,ph2_idx)*180/3.1415<=(85-5.*i)):
                if(ph1_det[0] == 1 and ph2_det[0] == 2):
                    eh_denom_v[i] += 1
                    if(invarmass(ecal_energy(ph_e_lw,ph1_idx),hcal_energy_nominal(ph_h,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>hcal_min and invarmass(ecal_energy(ph_e_lw,ph1_idx),hcal_energy_nominal(ph_h,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<=hcal_max):
                        eh_num_v[i] += 1
                if(ph1_det[0] == 2 and ph2_det[0] == 1):
                    he_denom_v[i] += 1
                    if(invarmass(hcal_energy_nominal(ph_h,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>hcal_min and invarmass(hcal_energy_nominal(ph_h,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<=hcal_max):
                        he_num_v[i] += 1
                if(ph1_det[0] == 2 and ph2_det[0] == 2):
                    hh_denom_v[i] += 1
                    if(invarmass(hcal_energy_nominal(ph_h,ph1_idx),hcal_energy_nominal(ph_h,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)>hcal_min and invarmass(hcal_energy_nominal(ph_h,ph1_idx),hcal_energy_nominal(ph_h,ph2_idx),ph_px,ph_py,ph_pz,ph1_idx,ph2_idx)<=hcal_max):
                        hh_num_v[i] += 1


        if(ph1_det[0] == 2 and ph2_det[0] ==2):
            h_ph_sep_dis.Fill(ph_sep(ph1_idx,ph2_idx,ph_px,ph_py,ph_pz,ph_thetaz))
            h_ph_sep_dis_x.Fill(ph_sep_x(ph1_idx,ph2_idx,ph_px,ph_py,ph_pz,ph_thetaz))
            h_ph_sep_dis_y.Fill(ph_sep_y(ph1_idx,ph2_idx,ph_px,ph_py,ph_pz,ph_thetaz))
#        #if(ph1_det[0] == 2 and ph2_det[0] == 2 and ph_h_top[ph1_idx[0]] > 0 and ph_h_top[ph2_idx[0]] > 0 and (ph_h_bottom[ph1_idx[0]]+ph_h_bottom[ph2_idx[0]]+ph_h_left[ph1_idx[0]]+ph_h_left[ph2_idx[0]]+ph_h_right[ph1_idx[0]]+ph_h_right[ph2_idx[0]]) < 1000):
#        #if(ph1_det[0] == 2 and ph2_det[0] == 2 and ph_h_top[ph1_idx[0]] > 0 and ph_h_top[ph2_idx[0]] > 0):
#        if(ph_h_top[ph1_idx[0]] > 0 and ph_h_top[ph2_idx[0]] > 0 and (ph_h_left[ph1_idx[0]]+ph_h_left[ph2_idx[0]]+ph_h_right[ph1_idx[0]]+ph_h_right[ph2_idx[0]]+ph_h_bottom[ph1_idx[0]]+ph_h_bottom[ph2_idx[0]]) < 50):
#            h_ph_sep_dis_top.Fill(ph_sep(ph1_idx,ph2_idx,ph_px,ph_py,ph_pz,ph_thetaz))
#            #h_ph_sep_dis_top.Fill(ph_h_top[ph1_idx[0]])
#        #if(ph1_det[0] == 2 and ph2_det[0] == 2 and ph_h_bottom[ph1_idx[0]] > 0 and ph_h_bottom[ph2_idx[0]] > 0):
#        if(ph_h_bottom[ph1_idx[0]] > 0 and ph_h_bottom[ph2_idx[0]] > 0 and (ph_h_left[ph1_idx[0]]+ph_h_left[ph2_idx[0]]+ph_h_top[ph1_idx[0]]+ph_h_top[ph2_idx[0]]+ph_h_right[ph1_idx[0]]+ph_h_right[ph2_idx[0]]) < 50):
#            h_ph_sep_dis_bottom.Fill(ph_sep(ph1_idx,ph2_idx,ph_px,ph_py,ph_pz,ph_thetaz))
#            #h_ph_sep_dis_bottom.Fill(ph_h_bottom[ph1_idx[0]])
#        #if(ph1_det[0] == 2 and ph2_det[0] == 2 and ph_h_left[ph1_idx[0]] > 0 and ph_h_left[ph2_idx[0]] > 0):
#        if(ph_h_left[ph1_idx[0]] > 0 and ph_h_left[ph2_idx[0]] > 0 and (ph_h_right[ph1_idx[0]]+ph_h_right[ph2_idx[0]]+ph_h_top[ph1_idx[0]]+ph_h_top[ph2_idx[0]]+ph_h_bottom[ph1_idx[0]]+ph_h_bottom[ph2_idx[0]]) < 50):
#            h_ph_sep_dis_left.Fill(ph_sep(ph1_idx,ph2_idx,ph_px,ph_py,ph_pz,ph_thetaz))
#            #h_ph_sep_dis_left.Fill(ph_h_left[ph1_idx[0]])
#        #if(ph1_det[0] == 2 and ph2_det[0] == 2 and ph_h_right[ph1_idx[0]] > 0 and ph_h_right[ph2_idx[0]] > 0):
#        if(ph_h_right[ph1_idx[0]] > 0 and ph_h_right[ph2_idx[0]] > 0 and (ph_h_left[ph1_idx[0]]+ph_h_left[ph2_idx[0]]+ph_h_top[ph1_idx[0]]+ph_h_top[ph2_idx[0]]+ph_h_bottom[ph1_idx[0]]+ph_h_bottom[ph2_idx[0]]) < 50):
#            h_ph_sep_dis_right.Fill(ph_sep(ph1_idx,ph2_idx,ph_px,ph_py,ph_pz,ph_thetaz))
#            #h_ph_sep_dis_right.Fill(ph_h_right[ph1_idx[0]])

#        if(ph_h_top[ph1_idx[0]]>0 and ph_h_top[ph2_idx[0]]>0):
#            print(f'entrynum:{entryNum}')
#            print(f'ph1idx:{ph1_idx[0]}')
#            print(f'ph_h_top:{ph_h_top[ph1_idx[0]]}')
#            print(f'ph_h_bottom:{ph_h_bottom[ph1_idx[0]]}')
#            h_temp.Fill(ph_h_bottom[ph1_idx[0]])

ana_tree.Project("h_bottom1","hcal_bottom_e[pi0_photon1_idx]","hcal_top_e[pi0_photon1_idx]>1.5&&hcal_top_e[pi0_photon2_idx]>1.5")
ana_tree.Project("h_bottom2","hcal_bottom_e[pi0_photon1_idx]","hcal_top_e[pi0_photon1_idx]>100&&hcal_top_e[pi0_photon2_idx]>100")

#print(ee_num)
#print(hh_num)

tmp_num_1 = np.add(ee_num,hh_num)
#print(tmp_num_1)
tmp_num_2 = np.add(eh_num,he_num)
num_total = np.add(tmp_num_1,tmp_num_2)
tmp_denom_1 = np.add(ee_denom,hh_denom)
tmp_denom_2 = np.add(eh_denom,he_denom)
denom_total = np.add(tmp_denom_1,tmp_denom_2)
denom_total_integral = denom_total.sum()

#print(denom1)
#print(denom_total_integral)

tmp_num_v_1 = np.add(eh_num,he_num)
num_v_total = np.add(tmp_num_v_1,hh_num_v)
tmp_denom_v_1 = np.add(eh_denom,he_denom)
denom_v_total = np.add(tmp_denom_v_1,hh_num_v)
denom_v_total_integral = denom_v_total.sum()

for i in range(1,17):
    h_thzmax_eff.SetBinContent(i,num_total[i]/denom_total_integral)
    #h_thzmax_eff.SetBinContent(i,num_total[i])

for i in range(1,17):
    h_thzmax_veto_eff.SetBinContent(i,num_v_total[i]/denom_v_total_integral)

print(h_thzmax_eff.Integral())

c1 = TCanvas("c1")
h_thzmax_eff.SetStats(0)
h_thzmax_eff.SetLineColor(1)
h_thzmax_eff.GetXaxis().SetTitle("hcal #theta_{max} [degrees]")
h_thzmax_eff.GetYaxis().SetTitle("integrated efficiency")
#h_thzmax_eff.GetYaxis().SetRangeUser(0,1)
h_thzmax_eff.Draw("hist")
c1.SaveAs("hcalthzmax_eff.png")

c2 = TCanvas("c2")
h_thzmax_veto_eff.SetStats(0)
h_thzmax_veto_eff.SetLineColor(1)
h_thzmax_veto_eff.GetXaxis().SetTitle("hcal #theta_{max} [degrees]")
h_thzmax_veto_eff.GetYaxis().SetTitle("integrated efficiency")
#h_thzmax_veto_eff.GetYaxis().SetRangeUser(0,1)
h_thzmax_veto_eff.Draw("hist")
c2.SaveAs("hcalthzmax_veto_eff.png")

c8 = TCanvas("c8")
h_bottom1.Draw("hist")
c8.SaveAs("bottom1.png")

c9 = TCanvas("c9")
h_bottom2.Draw("hist")
c9.SaveAs("bottom2.png")

#c4 = TCanvas("c4")
#h_ph_sep_dis_top.Draw("hist")
#c4.SaveAs("ph_sep_top.png")

#c5 = TCanvas("c5")
#h_ph_sep_dis_bottom.Draw("hist")
#c5.SaveAs("ph_sep_bottom.png")

#c6 = TCanvas("c6")
#h_ph_sep_dis_left.Draw("hist")
#c6.SaveAs("ph_sep_left.png")

#c7 = TCanvas("c7")
#h_ph_sep_dis_right.Draw("hist")
#c7.SaveAs("ph_sep_right.png")

#h_ph_sep_dis = h_ph_sep_dis_top.Clone()
#h_ph_sep_dis.Add(h_ph_sep_dis_bottom)
#h_ph_sep_dis.Add(h_ph_sep_dis_left)
#h_ph_sep_dis.Add(h_ph_sep_dis_right)

h_ph_sep_dis.Scale(1./h_ph_sep_dis.Integral())

c3 = TCanvas("c3")
h_ph_sep_dis.SetStats(0)
h_ph_sep_dis.GetXaxis().SetTitle("#gamma#gamma separation [cm]")
h_ph_sep_dis.GetYaxis().SetTitle("fraction of events")
h_ph_sep_dis.SetLineColor(1)
#h_ph_sep_dis_top.SetLineColor(2)
#h_ph_sep_dis_bottom.SetLineColor(2)
#h_ph_sep_dis_bottom.SetLineStyle(2)
#h_ph_sep_dis_left.SetLineColor(2)
#h_ph_sep_dis_left.SetLineStyle(3)
#h_ph_sep_dis_right.SetLineColor(2)
#h_ph_sep_dis_right.SetLineStyle(4)
h_ph_sep_dis.Draw("hist")
#h_ph_sep_dis_top.Draw("hist same")
#h_ph_sep_dis_bottom.Draw("hist same")
#h_ph_sep_dis_left.Draw("hist same")
#h_ph_sep_dis_right.Draw("hist same")
c3.SaveAs("h_ph_sep_dis.png")

h_ph_sep_dis_x.Scale(1./h_ph_sep_dis_x.Integral())
h_ph_sep_dis_y.Scale(1./h_ph_sep_dis_y.Integral())

c4 = TCanvas("c4")
h_ph_sep_dis_y.SetStats(0)
h_ph_sep_dis_y.GetXaxis().SetTitle("|#Delta#gamma_{x/y}| [cm]")
h_ph_sep_dis_y.GetYaxis().SetTitle("fraction of events")
h_ph_sep_dis_y.SetLineColor(1)
h_ph_sep_dis_y.SetLineColor(2)
h_ph_sep_dis_y.Draw("hist")
h_ph_sep_dis_x.Draw("hist same")
l_xy = TLegend(0.67,0.68,0.9,0.88)
l_xy.SetBorderSize(0)
l_xy.AddEntry(h_ph_sep_dis_x,"Difference in x coords")
l_xy.AddEntry(h_ph_sep_dis_y,"Difference in y coords")
l_xy.Draw("hist same")
c4.SaveAs("h_ph_sep_dis_coords.png")

