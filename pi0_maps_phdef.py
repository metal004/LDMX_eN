from ROOT import *
import numpy as np
import random

myfile = TFile("pi0_maps_phdef.root","RECREATE")
f = TFile("1pi0_6May.root")
ana_tree = f.Get("ana_tree")

# calibration factors
#cf_ecal_mean = 125.
#cf_ecal_median = 134.7
#cf_hcal_mean = 10.45
#cf_hcal_mean = 9.83
#cf_hcal_median = 10.59
cf_lw = 1.60
#cf_byslices = [13.1,10.2,9.6,9.6,9.3,9.8,8.7,13.9]
#cf_byslices500keV = [9.1,9.9,9.7,9.7,9.3,9.8,9.2,12.5] # Previoius
#cf_byslices500keV = [9.1,9.9,9.7,9.7,9.3,9.8,9.2,10.8]
cf_byslices500keV = [9.25,10.27,9.81,9.83,9.74,9.57,8.99,9.91]
cf_byslices = [22., 11.1, 8.9, 8.4, 8.]

# Mean and cutoffs
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

def cos_theta(px, py, pz, ph1, ph2):
    return (px[ph1[0]]*px[ph2[0]]+py[ph1[0]]*py[ph2[0]]+pz[ph1[0]]*pz[ph2[0]])/(np.sqrt(px[ph1[0]]*px[ph1[0]]+py[ph1[0]]*py[ph1[0]]+pz[ph1[0]]*pz[ph1[0]])*np.sqrt(px[ph2[0]]*px[ph2[0]]+py[ph2[0]]*py[ph2[0]]+pz[ph2[0]]*pz[ph2[0]]))

def invarmass(e1,e2,px,py,pz,ph1,ph2):
    return np.sqrt(2.*e1*e2*(1-cos_theta(px,py,pz,ph1,ph2)))

def ecal_energy(e,ph):
    return cf_lw*e[ph[0]]

#def hcal_energy(e,ph):
#    return cf_hcal_mean*e[ph[0]]

#def hcal_energy(e,ph):
#    if(cf_byslices500keV[0]*e[ph[0]]>=0 and cf_byslices500keV[0]*e[ph[0]]<100):
#        return cf_byslices500keV[0]*e[ph[0]]
#    if(cf_byslices500keV[1]*e[ph[0]]>=100 and cf_byslices500keV[1]*e[ph[0]]<200):
#        return cf_byslices500keV[1]*e[ph[0]]
#    if(cf_byslices500keV[2]*e[ph[0]]>=200 and cf_byslices500keV[2]*e[ph[0]]<300):
#        return cf_byslices500keV[2]*e[ph[0]]
#    if(cf_byslices500keV[3]*e[ph[0]]>=300 and cf_byslices500keV[3]*e[ph[0]]<400):
#        return cf_byslices500keV[3]*e[ph[0]]
#    if(cf_byslices500keV[4]*e[ph[0]]>=400 and cf_byslices500keV[4]*e[ph[0]]<500):
#        return cf_byslices500keV[4]*e[ph[0]]
#    if(cf_byslices500keV[5]*e[ph[0]]>=500 and cf_byslices500keV[5]*e[ph[0]]<600):
#        return cf_byslices500keV[5]*e[ph[0]]
#    if(cf_byslices500keV[6]*e[ph[0]]>=600 and cf_byslices500keV[6]*e[ph[0]]<800):
#        return cf_byslices500keV[6]*e[ph[0]]
#    if(cf_byslices500keV[7]*e[ph[0]]>=800 and cf_byslices500keV[7]*e[ph[0]]<1400):
#        return cf_byslices500keV[7]*e[ph[0]]

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

# historgram definitions
# Energy
h_ee_e_ph_l = TH1F("h_ee_e_ph_l","#pi^{0} Efficiency (photon energy)",40,0,2500)
h_ee_e_ph_sl = TH1F("h_ee_e_ph_sl","#pi^{0} Sub-Leading Photon Energy Selected Events",40,0,2500)
h_hh_e_ph_l = TH1F("h_hh_e_ph_l","#pi^{0} Leading Photon Energy [hh]",40,0,2500)
h_hh_e_ph_sl = TH1F("h_hh_e_ph_sl","#pi^{0} Sub-Leading Photon Energy [hh]",40,0,2500)
h_eh_e_ph_l = TH1F("h_eh_e_ph_l","#pi^{0} Leading Photon Energy [eh]",40,0,2500)
h_eh_e_ph_sl = TH1F("h_eh_e_ph_sl","#pi^{0} Sub-Leading Photon Energy [eh]",40,0,2500)
h_e_total_ph_l = TH1F("h_e_total_ph_l","",40,0,2500)
h_e_total_ph_sl = TH1F("h_e_total_ph_sl","",40,0,2500)
h_e_sum_ph_l = TH1F("h_e_sum_ph_l","",40,0,2500)
h_e_sum_ph_l_dist = TH1F("h_e_sum_ph_l_dist","",40,0,2500)
h_e_sum_ph_sl = TH1F("h_e_sum_ph_sl","",40,0,2500)
h_e_sum_ph_sl_dist = TH1F("h_e_sum_ph_sl_dist","",40,0,2500)

h_ee_e_pi0 = TH1F("h_ee_e_pi0","#pi^{0} Efficiency (#pi^{0} energy)",40,0,2500)
h_hh_e_pi0 = TH1F("h_hh_e_pi0","#pi^{0} Energy [hh]",40,0,2500)
h_eh_e_pi0 = TH1F("h_eh_e_pi0","#pi^{0} Energy [eh]",40,0,2500)
h_e_total_pi0 = TH1F("h_e_total_pi0","",40,0,2500)
h_e_sum_pi0 = TH1F("h_e_sum_pi0","",40,0,2500)

h_ee_thz_ph_l = TH1F("h_ee_thz_ph_l","#pi^{0} Efficiency (photon angle)",50,0,100)
h_ee_thz_ph_sl = TH1F("h_ee_thz_ph_sl","#pi^{0} Sub-Leading Photon Angle Selected Events",50,0,100)
h_hh_thz_ph_l = TH1F("h_hh_thz_ph_l","#pi^{0} Leading Photon Angle [hh]",50,0,100)
h_hh_thz_ph_sl = TH1F("h_hh_thz_ph_sl","#pi^{0} Sub-Leading Photon Angle [hh]",50,0,100)
h_eh_thz_ph_l = TH1F("h_eh_thz_ph_l","#pi^{0} Leading Photon Angle [eh]",50,0,100)
h_eh_thz_ph_sl = TH1F("h_eh_thz_ph_sl","#pi^{0} Sub-Leading Photon Angle [eh]",50,0,100)
h_thz_total_ph_l = TH1F("h_thz_total_ph_l","",50,0,100)
h_thz_total_ph_sl = TH1F("h_thz_total_ph_sl","",50,0,100)
h_thz_sum_ph_l = TH1F("h_thz_sum_ph_l","",50,0,100)
h_thz_sum_ph_sl = TH1F("h_thz_sum_ph_sl","",50,0,100)
h_thz_sum_ph_l_dist = TH1F("h_thz_sum_ph_l_dist","",50,0,100)
h_thz_sum_ph_sl_dist = TH1F("h_thz_sum_ph_sl_dist","",50,0,100)

h_ee_thz_pi0 = TH1F("h_ee_thz_pi0","#pi^{0} Efficiency (#pi^{0} angle)",50,0,100)
h_hh_thz_pi0 = TH1F("h_hh_thz_pi0","#pi^{0} Angle [hh]",50,0,100)
h_eh_thz_pi0 = TH1F("h_eh_thz_pi0","#pi^{0} Angle [eh]",50,0,100)
h_thz_total_pi0 = TH1F("h_thz_total_pi0","",50,0,100)
h_thz_sum_pi0 = TH1F("h_thz_sum_pi0","",50,0,100)

# 2d histograms
# photon energy
h_ee_2d_ph_l = TH2F("h_ee_2d_ph_l","Efficiency Map [ecal-ecal] (leading photon kinematics)",40,0,4000,40,0,100)
h_hh_2d_ph_l = TH2F("h_hh_2d_ph_l","Efficiency Map [hcal-hcal] (leading photon kinematics)",40,0,4000,40,0,100)
h_eh_2d_ph_l = TH2F("h_eh_2d_ph_l","Efficiency Map [ecal-hcal] (leading photon kinematics)",40,0,4000,40,0,100)
h_ee_2d_ph_sl = TH2F("h_ee_2d_ph_sl","Efficiency Map [ecal-ecal] (sub-leading photon kinematics)",40,0,4000,40,0,100)
h_hh_2d_ph_sl = TH2F("h_hh_2d_ph_sl","Efficiency Map [hcal-hcal] (sub-leading photon kinematics)",40,0,4000,40,0,100)
h_eh_2d_ph_sl = TH2F("h_eh_2d_ph_sl","Efficiency Map [ecal-hcal] (sub-leading photon kinematics)",40,0,4000,40,0,100)
h_total_2d_ph_l = TH2F("h_total_2d_ph_l","",40,0,4000,40,0,100)
h_total_2d_ph_sl = TH2F("h_total_2d_ph_sl","",40,0,4000,40,0,100)
h_sum_2d_ph_l = TH2F("h_sum_2d_ph_l","#pi^{0} Photon Efficiency Map",40,0,4000,40,0,100)
h_sum_2d_ph_sl = TH2F("h_sum_2d_ph_sl","#pi^{0} Photon Efficiency Map",40,0,4000,40,0,100)

h_ee_2d_pi0 = TH2F("h_ee_2d_pi0","#pi^{0} Efficiency Map (#pi^{0} kinematics)",40,0,4000,40,0,100)
h_hh_2d_pi0 = TH2F("h_hh_2d_pi0","",40,0,4000,40,0,100)
h_eh_2d_pi0 = TH2F("h_eh_2d_pi0","",40,0,4000,40,0,100)
h_total_2d_pi0 = TH2F("h_total_2d_pi0","",40,0,4000,40,0,100)
h_sum_2d_pi0 = TH2F("h_sum_2d_pi0","#pi^{0} Efficiency Map",40,0,4000,40,0,100)

# separating 2d maps by topology
h_ee_2d_ph_l_clone = TH2F("h_ee_2d_ph_l_clone","Efficiency Map [ecal-ecal] (leading photon kinematics)",40,0,4000,40,0,100)
h_ee_2d_ph_sl_clone = TH2F("h_ee_2d_ph_sl_clone","Efficiency Map [ecal-ecal] (sub-leading photon kinematics)",40,0,4000,40,0,100)
h_ee_2d_ph_l_clone_clone = TH2F("h_ee_2d_ph_l_clone_clone","Efficiency Map [ecal-ecal] (leading photon kinematics)",40,0,4000,40,0,100)
h_ee_2d_ph_sl_clone_clone = TH2F("h_ee_2d_ph_sl_clone_clone","Efficiency Map [ecal-ecal] (sub-leading photon kinematics)",40,0,4000,40,0,100)
h_ee_2d_ph_l_tot = TH2F("h_ee_2d_ph_l_tot","Efficiency Map [ecal-ecal] (leading photon kinematics)",40,0,4000,40,0,100)
h_hh_2d_ph_l_tot = TH2F("h_hh_2d_ph_l_tot","Efficiency Map [hcal-hcal] (leading photon kinematics)",40,0,4000,40,0,100)
h_eh_2d_ph_l_tot = TH2F("h_eh_2d_ph_l_tot","Efficiency Map [ecal-hcal] (leading photon kinematics)",40,0,4000,40,0,100)
h_ee_2d_ph_sl_tot = TH2F("h_ee_2d_ph_sl_tot","Efficiency Map [ecal-ecal] (sub-leading photon kinematics)",40,0,4000,40,0,100)
h_hh_2d_ph_sl_tot = TH2F("h_hh_2d_ph_sl_tot","Efficiency Map [hcal-hcal] (sub-leading photon kinematics)",40,0,4000,40,0,100)
h_eh_2d_ph_sl_tot = TH2F("h_eh_2d_ph_sl_tot","Efficiency Map [ecal-hcal] (sub-leading photon kinematics)",40,0,4000,40,0,100)
h_hh_2d_ph_l_clone = TH2F("h_hh_2d_ph_l_clone","Efficiency Map [hcal-hcal] (leading photon kinematics)",40,0,4000,40,0,100)
h_hh_2d_ph_sl_clone = TH2F("h_hh_2d_ph_sl_clone","Efficiency Map [hcal-hcal] (leading photon kinematics)",40,0,4000,40,0,100)
h_eh_2d_ph_l_clone = TH2F("h_eh_2d_ph_l_clone","Efficiency Map [hcal-hcal] (leading photon kinematics)",40,0,4000,40,0,100)
h_eh_2d_ph_sl_clone = TH2F("h_eh_2d_ph_sl_clone","Efficiency Map [hcal-hcal] (leading photon kinematics)",40,0,4000,40,0,100)
h_2d_pi0_ee = TH2F("h_2d_pi0_ee","Efficiency Map [ecal-ecal]",40,0,4000,40,0,100)
h_2d_pi0_hh = TH2F("h_2d_pi0_hh","Efficiency Map [hcal-hcal]",40,0,4000,40,0,100)
h_2d_pi0_eh = TH2F("h_2d_pi0_eh","Efficiency Map [ecal-hcal]",40,0,4000,40,0,100)
h_2d_pi0_ee_clone = TH2F("h_2d_pi0_ee_clone","Efficiency Map [ecal-ecal]",40,0,4000,40,0,100)
h_2d_pi0_hh_clone = TH2F("h_2d_pi0_hh_clone","Efficiency Map [hcal-hcal]",40,0,4000,40,0,100)
h_2d_pi0_eh_clone = TH2F("h_2d_pi0_eh_clone","Efficiency Map [ecal-hcal]",40,0,4000,40,0,100)

for entryNum in range(0,ana_tree.GetEntries()):
    ana_tree.GetEntry(entryNum)
    pi0_idx = getattr(ana_tree,"pi0_idx")
    ph1_idx = getattr(ana_tree,"pi0_photon1_idx")
    ph2_idx = getattr(ana_tree,"pi0_photon2_idx")
    ph_e = getattr(ana_tree,"sim_p_ecal_e")
    ph_h = getattr(ana_tree,"sim_p_hcal_e")
    ph_e_lw = getattr(ana_tree,"sim_p_ecal_e_lw")
    p_px = getattr(ana_tree,"sim_p_px")
    p_py = getattr(ana_tree,"sim_p_py")
    p_pz = getattr(ana_tree,"sim_p_pz")
    p_e_true = getattr(ana_tree,"sim_p_e")
    p_thetaz = getattr(ana_tree,"sim_p_thetaz")
    ph1_det = getattr(ana_tree,"pi0_photon1_det")
    ph2_det = getattr(ana_tree,"pi0_photon2_det")
    if(len(ph_e)>ph1_idx[0] and len(ph_h)>ph1_idx[0] and len(ph_e)>ph2_idx[0] and len(ph_h)>ph2_idx[0]):
        h_e_total_pi0.Fill(p_e_true[pi0_idx[0]])
        h_thz_total_pi0.Fill(p_thetaz[pi0_idx[0]]*180/3.1415)
        h_total_2d_ph_l.Fill(p_e_true[ph1_idx[0]],p_thetaz[ph1_idx[0]]*180/3.1415)
        h_total_2d_ph_sl.Fill(p_e_true[ph2_idx[0]],p_thetaz[ph2_idx[0]]*180/3.1415)
        h_total_2d_pi0.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)
        h_e_total_ph_l.Fill(p_e_true[ph1_idx[0]])
        h_e_total_ph_sl.Fill(p_e_true[ph2_idx[0]])
        h_thz_total_ph_l.Fill(p_thetaz[ph1_idx[0]]*180/3.1415)
        h_thz_total_ph_sl.Fill(p_thetaz[ph2_idx[0]]*180/3.1415)

        #if(ph_e_lw[ph1_idx[0]]>0 and ph_h[ph1_idx[0]]==0 and ph_e_lw[ph2_idx[0]]>0 and ph_h[ph2_idx[0]]==0):
        if(ph1_det[0]==1 and ph2_det[0]==1):
            h_ee_2d_ph_l_tot.Fill(p_e_true[ph1_idx[0]],p_thetaz[ph1_idx[0]]*180/3.1415)
            h_ee_2d_ph_sl_tot.Fill(p_e_true[ph2_idx[0]],p_thetaz[ph2_idx[0]]*180/3.1415)
            # ee case
            h_2d_pi0_ee_clone.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)
            if(invarmass(ecal_energy(ph_e_lw,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),p_px,p_py,p_pz,ph1_idx,ph2_idx)<ecal_max and invarmass(ecal_energy(ph_e_lw,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),p_px,p_py,p_pz,ph1_idx,ph2_idx)>ecal_min):
                h_ee_e_pi0.Fill(p_e_true[pi0_idx[0]])
                h_ee_thz_pi0.Fill(p_thetaz[pi0_idx[0]]*180/3.1415)
                h_ee_2d_pi0.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)
                h_ee_e_ph_l.Fill(p_e_true[ph1_idx[0]])
                h_ee_thz_ph_l.Fill(p_thetaz[ph1_idx[0]]*180/3.1415)
                h_ee_2d_ph_l.Fill(p_e_true[ph1_idx[0]],p_thetaz[ph1_idx[0]]*180/3.1415)
                h_ee_e_ph_sl.Fill(p_e_true[ph2_idx[0]])
                h_ee_thz_ph_sl.Fill(p_thetaz[ph2_idx[0]]*180/3.1415)
                h_ee_2d_ph_sl.Fill(p_e_true[ph2_idx[0]],p_thetaz[ph2_idx[0]]*180/3.1415)
                h_2d_pi0_ee.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)

        # hh case
        #if(ph_h[ph1_idx[0]]>0.5 and ph_e_lw[ph1_idx[0]]==0 and ph_h[ph2_idx[0]]>0.5 and ph_e_lw[ph2_idx[0]]==0):
        if(ph1_det[0]==2 and ph2_det[0]==2):
            h_hh_2d_ph_l_tot.Fill(p_e_true[ph1_idx[0]],p_thetaz[ph1_idx[0]]*180/3.1415)
            h_hh_2d_ph_sl_tot.Fill(p_e_true[ph2_idx[0]],p_thetaz[ph2_idx[0]]*180/3.1415)
            h_2d_pi0_hh_clone.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)
            if(invarmass(hcal_energy(ph_h,ph1_idx),hcal_energy(ph_h,ph2_idx),p_px,p_py,p_pz,ph1_idx,ph2_idx)<hcal_max and invarmass(hcal_energy(ph_h,ph1_idx),hcal_energy(ph_h,ph2_idx),p_px,p_py,p_pz,ph1_idx,ph2_idx)>hcal_min):
                h_hh_e_pi0.Fill(p_e_true[pi0_idx[0]])
                h_hh_thz_pi0.Fill(p_thetaz[pi0_idx[0]]*180/3.1415)
                h_hh_2d_pi0.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)
                h_hh_e_ph_l.Fill(p_e_true[ph1_idx[0]])
                h_hh_thz_ph_l.Fill(p_thetaz[ph1_idx[0]]*180/3.1415)
                h_hh_2d_ph_l.Fill(p_e_true[ph1_idx[0]],p_thetaz[ph1_idx[0]]*180/3.1415)
                h_hh_e_ph_sl.Fill(p_e_true[ph2_idx[0]])
                h_hh_thz_ph_sl.Fill(p_thetaz[ph2_idx[0]]*180/3.1415)
                h_hh_2d_ph_sl.Fill(p_e_true[ph2_idx[0]],p_thetaz[ph2_idx[0]]*180/3.1415)
                h_2d_pi0_hh.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)

        # eh case
        #if((ph_e_lw[ph1_idx[0]]>0 and ph_h[ph1_idx[0]]==0 and ph_h[ph2_idx[0]]>0.5 and ph_e_lw[ph2_idx[0]]==0) or (ph_e_lw[ph2_idx[0]]>0 and ph_h[ph2_idx[0]]==0 and ph_h[ph1_idx[0]]>0.5 and ph_e_lw[ph1_idx[0]]==0)):
        if((ph1_det[0]==1 and ph2_det[0]==2) or (ph1_det[0]==2 and ph2_det[0]==1)):
            h_eh_2d_ph_l_tot.Fill(p_e_true[ph1_idx[0]],p_thetaz[ph1_idx[0]]*180/3.1415)
            h_eh_2d_ph_sl_tot.Fill(p_e_true[ph2_idx[0]],p_thetaz[ph2_idx[0]]*180/3.1415)
            h_2d_pi0_eh_clone.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)
            if((invarmass(ecal_energy(ph_e_lw,ph1_idx),hcal_energy(ph_h,ph2_idx),p_px,p_py,p_pz,ph1_idx,ph2_idx)<eh_max and invarmass(ecal_energy(ph_e_lw,ph1_idx),hcal_energy(ph_h,ph2_idx),p_px,p_py,p_pz,ph1_idx,ph2_idx)>eh_min) or (invarmass(hcal_energy(ph_h,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),p_px,p_py,p_pz,ph1_idx,ph2_idx)<eh_max and invarmass(hcal_energy(ph_h,ph1_idx),ecal_energy(ph_e_lw,ph2_idx),p_px,p_py,p_pz,ph1_idx,ph2_idx)<eh_min)):
                h_eh_thz_pi0.Fill(p_thetaz[pi0_idx[0]]*180/3.1415)
                h_eh_e_pi0.Fill(p_e_true[pi0_idx[0]])
                h_eh_2d_pi0.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)
                h_eh_e_ph_l.Fill(p_e_true[ph1_idx[0]])
                h_eh_thz_ph_l.Fill(p_thetaz[ph1_idx[0]]*180/3.1415)
                h_eh_2d_ph_l.Fill(p_e_true[ph1_idx[0]],p_thetaz[ph1_idx[0]]*180/3.1415)
                h_eh_e_ph_sl.Fill(p_e_true[ph2_idx[0]])
                h_eh_thz_ph_sl.Fill(p_thetaz[ph2_idx[0]]*180/3.1415)
                h_eh_2d_ph_sl.Fill(p_e_true[ph2_idx[0]],p_thetaz[ph2_idx[0]]*180/3.1415)
                h_2d_pi0_eh.Fill(p_e_true[pi0_idx[0]],p_thetaz[pi0_idx[0]]*180/3.1415)

# clones
h_hh_2d_ph_l_clone = h_hh_2d_ph_l.Clone()
h_hh_2d_ph_sl_clone = h_hh_2d_ph_sl.Clone()
h_eh_2d_ph_l_clone = h_eh_2d_ph_l.Clone()
h_eh_2d_ph_sl_clone = h_eh_2d_ph_sl.Clone()

# photon Energy
h_e_sum_ph_l = h_ee_e_ph_l.Clone()
h_e_sum_ph_l.Add(h_hh_e_ph_l)
h_e_sum_ph_l.Add(h_eh_e_ph_l)
h_e_sum_ph_sl = h_ee_e_ph_sl.Clone()
h_e_sum_ph_sl.Add(h_hh_e_ph_sl)
h_e_sum_ph_sl.Add(h_eh_e_ph_sl)

h_e_sum_ph_l.Divide(h_e_total_ph_l)
h_e_sum_ph_sl.Divide(h_e_total_ph_sl)

h_ee_e_ph_l.Divide(h_e_total_ph_l)
h_hh_e_ph_l.Divide(h_e_total_ph_l)
h_eh_e_ph_l.Divide(h_e_total_ph_l)

h_e_sum_ph_l_dist = h_e_total_ph_l.Clone()
h_e_sum_ph_sl_dist = h_e_total_ph_sl.Clone()
h_e_sum_ph_l_dist.Scale(1./h_e_sum_ph_l_dist.Integral())
h_e_sum_ph_sl_dist.Scale(1./h_e_sum_ph_sl_dist.Integral())

# ROOT line colors 9, 30, 46

c1 = TCanvas("c1")
h_e_sum_ph_l.SetLineColor(1)
h_e_sum_ph_l_dist.SetLineColor(1)
h_e_sum_ph_l_dist.SetLineStyle(2)
h_e_sum_ph_sl.SetLineColor(4)
h_e_sum_ph_sl_dist.SetLineColor(4)
h_e_sum_ph_sl_dist.SetLineStyle(2)
h_e_sum_ph_l.GetXaxis().SetTitle("Photon Energy")
h_e_sum_ph_l.GetYaxis().SetRangeUser(0,1)
h_e_sum_ph_l.SetStats(0)
h_e_sum_ph_l.Draw("hist")
h_e_sum_ph_sl.Draw("hist same")
h_e_sum_ph_l_dist.Draw("hist same")
h_e_sum_ph_sl_dist.Draw("hist same")
#h_ee_e_ph_l.SetLineColor(9)
#h_ee_e_ph_l.Draw("hist same")
#h_hh_e_ph_l.SetLineColor(30)
#h_hh_e_ph_l.Draw("hist same")
#h_eh_e_ph_l.SetLineColor(46)
#h_eh_e_ph_l.Draw("hist same")
l1 = TLegend(0.45,0.15,0.9,0.5)
l1.AddEntry(h_e_sum_ph_l,"Leading photon current efficiency")
l1.AddEntry(h_e_sum_ph_sl,"Sub-leading photon current efficiency")
l1.AddEntry(h_e_sum_ph_l_dist,"#splitline{Leading photon distribution}{(normalized to 1)}")
l1.AddEntry(h_e_sum_ph_sl_dist,"#splitline{Sub-leading photon distribution}{(normalized to 1)}")
#l1.AddEntry(h_ee_e_ph_l,"Leading photon EE efficiency")
#l1.AddEntry(h_hh_e_ph_l,"Leading photon HH efficiency")
#l1.AddEntry(h_eh_e_ph_l,"Leading photon EH efficiency")
#l1.Draw("hist same")
c1.SaveAs("efficiency_e_ph_kin.png")
myfile.WriteObject(c1,"phkin_e")

# pi0 Energy
h_e_sum_pi0 = h_ee_e_pi0.Clone()
h_e_sum_pi0.Add(h_hh_e_pi0)
h_e_sum_pi0.Add(h_eh_e_pi0)

h_e_sum_pi0.Divide(h_e_total_pi0)

h_ee_e_pi0.Divide(h_e_total_pi0)
h_hh_e_pi0.Divide(h_e_total_pi0)
h_eh_e_pi0.Divide(h_e_total_pi0)

c2 = TCanvas("c2")
h_e_sum_pi0.SetLineColor(1)
h_e_sum_pi0.GetXaxis().SetTitle("#pi^{0} Energy")
h_e_sum_pi0.GetYaxis().SetRangeUser(0,1)
h_e_sum_pi0.SetStats(0)
h_e_sum_pi0.Draw("hist")
h_ee_e_pi0.SetLineColor(9)
h_hh_e_pi0.SetLineColor(30)
h_eh_e_pi0.SetLineColor(46)
h_ee_e_pi0.Draw("hist same")
h_hh_e_pi0.Draw("hist same")
h_eh_e_pi0.Draw("hist same")
l3 = TLegend(0.1,0.6,0.36,0.85)
l3.AddEntry(h_e_sum_pi0,"Current efficiency")
l3.AddEntry(h_ee_e_pi0,"EE efficiency")
l3.AddEntry(h_hh_e_pi0,"HH efficiency")
l3.AddEntry(h_eh_e_pi0,"EH efficiency")
l3.Draw("hist same")
c2.SaveAs("efficiency_e_pi0_kin.png")
myfile.WriteObject(c2,"pi0kin_e")

### ph validation
##c4 = TCanvas("c4")
##h_ee_thz_ph.Draw("hist")
##c4.SaveAs("h_ee_thz_ph.png")
##
##c5 = TCanvas("c5")
##h_hh_thz_ph.Draw("hist")
##c5.SaveAs("h_hh_thz_ph.png")
##
##c6 = TCanvas("c6")
##h_eh_thz_ph.Draw("hist")
##c6.SaveAs("h_eh_thz_ph.png")

# photon Angle
h_thz_sum_ph_l = h_ee_thz_ph_l.Clone()
h_thz_sum_ph_l.Add(h_hh_thz_ph_l)
h_thz_sum_ph_l.Add(h_eh_thz_ph_l)
h_thz_sum_ph_sl = h_ee_thz_ph_sl.Clone()
h_thz_sum_ph_sl.Add(h_hh_thz_ph_sl)
h_thz_sum_ph_sl.Add(h_eh_thz_ph_sl)

h_thz_sum_ph_l.Divide(h_thz_total_ph_l)
h_thz_sum_ph_sl.Divide(h_thz_total_ph_sl)

h_thz_sum_ph_l_dist = h_thz_total_ph_l.Clone()
h_thz_sum_ph_sl_dist = h_thz_total_ph_sl.Clone()
h_thz_sum_ph_l_dist.Scale(1./h_thz_sum_ph_l_dist.Integral())
h_thz_sum_ph_sl_dist.Scale(1./h_thz_sum_ph_sl_dist.Integral())

c3 = TCanvas("c3")
h_thz_sum_ph_l.SetLineColor(1)
h_thz_sum_ph_l_dist.SetLineColor(1)
h_thz_sum_ph_l_dist.SetLineStyle(2)
h_thz_sum_ph_sl.SetLineColor(4)
h_thz_sum_ph_sl_dist.SetLineColor(4)
h_thz_sum_ph_sl_dist.SetLineStyle(2)
h_thz_sum_ph_l.GetXaxis().SetTitle("Photon Angle [degrees]")
h_thz_sum_ph_l.GetYaxis().SetRangeUser(0,1)
h_thz_sum_ph_l.SetStats(0)
h_thz_sum_ph_l.Draw("hist")
h_thz_sum_ph_sl.Draw("hist same")
h_thz_sum_ph_l_dist.Draw("hist same")
h_thz_sum_ph_sl_dist.Draw("hist same")
l2 = TLegend(0.45,0.5,0.9,0.85)
l2.AddEntry(h_thz_sum_ph_l,"Leading photon current Efficiency")
l2.AddEntry(h_thz_sum_ph_sl,"Sub-leading photon current efficiency")
l2.AddEntry(h_thz_sum_ph_l_dist,"#splitline{Leading photon distribution}{(normalized to 1)}")
l2.AddEntry(h_thz_sum_ph_sl_dist,"#splitline{Sub-leading photon distribution}{(normalized to 1)}")
#l2.Draw("hist same")
c3.SaveAs("efficiency_thz_ph_kin.png")
myfile.WriteObject(c3,"phkin_thz")

# pi0 Angle
h_thz_sum_pi0 = h_ee_thz_pi0.Clone()
h_thz_sum_pi0.Add(h_hh_thz_pi0)
h_thz_sum_pi0.Add(h_eh_thz_pi0)

h_thz_sum_pi0.Divide(h_thz_total_pi0)

h_ee_thz_pi0.Divide(h_thz_total_pi0)
h_hh_thz_pi0.Divide(h_thz_total_pi0)
h_eh_thz_pi0.Divide(h_thz_total_pi0)

c4 = TCanvas("c4")
h_thz_sum_pi0.SetLineColor(1)
h_thz_sum_pi0.GetXaxis().SetTitle("#pi^{0} Angle [degrees]")
h_thz_sum_pi0.GetYaxis().SetRangeUser(0,0.6)
h_thz_sum_pi0.GetXaxis().SetRangeUser(0,80)
h_thz_sum_pi0.SetStats(0)
h_thz_sum_pi0.Draw("hist")
h_ee_thz_pi0.SetLineColor(9)
h_hh_thz_pi0.SetLineColor(30)
h_eh_thz_pi0.SetLineColor(46)
h_ee_thz_pi0.Draw("hist same")
h_hh_thz_pi0.Draw("hist same")
h_eh_thz_pi0.Draw("hist same")
l4 = TLegend(0.45,0.55,0.75,0.7)
l4.SetBorderSize(0)
l4.AddEntry(h_e_sum_pi0,"Current efficiency")
l4.AddEntry(h_ee_e_pi0,"EE efficiency")
l4.AddEntry(h_hh_e_pi0,"HH efficiency")
l4.AddEntry(h_eh_e_pi0,"EH efficiency")
l4.Draw("hist same")
c4.SaveAs("efficiency_thz_pi0_kin.png")
myfile.WriteObject(c4,"pi0kin_thz")

#c5 = TCanvas("c5")
#h_total_2d_ph.Draw("COLZ")
#c5.SaveAs("temp.png")

# clone
h_ee_2d_ph_l_clone = h_ee_2d_ph_l.Clone()
h_ee_2d_ph_sl_clone = h_ee_2d_ph_sl.Clone()

h_ee_2d_ph_l_clone_clone = h_ee_2d_ph_l.Clone()
h_ee_2d_ph_sl_clone_clone = h_ee_2d_ph_sl.Clone()

# 2d e ph
h_sum_2d_ph_l = h_ee_2d_ph_l.Clone()
h_sum_2d_ph_l.Add(h_hh_2d_ph_l)
h_sum_2d_ph_l.Add(h_eh_2d_ph_l)

h_sum_2d_ph_l.Divide(h_total_2d_ph_l)

h_sum_2d_ph_sl = h_ee_2d_ph_sl.Clone()
h_sum_2d_ph_sl.Add(h_hh_2d_ph_sl)
h_sum_2d_ph_sl.Add(h_eh_2d_ph_sl)

h_sum_2d_ph_sl.Divide(h_total_2d_ph_sl)


c6 = TCanvas("c6")
h_sum_2d_ph_l.GetXaxis().SetTitle("Leading Photon Energy")
h_sum_2d_ph_l.GetYaxis().SetTitle("Leading Photon Angle [degrees]")
h_sum_2d_ph_l.SetStats(0)
h_sum_2d_ph_l.Draw("COLZ")
c6.SaveAs("efficiency_map_e_ph_l_kin.png")

c7 = TCanvas("c7")
h_sum_2d_ph_sl.GetXaxis().SetTitle("Sub-leading Photon Energy")
h_sum_2d_ph_sl.GetYaxis().SetTitle("Sub-leading Photon Angle [degrees]")
h_sum_2d_ph_sl.SetStats(0)
h_sum_2d_ph_sl.Draw("COLZ")
c7.SaveAs("efficiency_map_e_ph_sl_kin.png")

h_sum_2d_pi0 = h_ee_2d_pi0.Clone()
h_sum_2d_pi0.Add(h_hh_2d_pi0)
h_sum_2d_pi0.Add(h_eh_2d_pi0)

h_sum_2d_pi0.Divide(h_total_2d_pi0)

c8 = TCanvas("c8")
h_sum_2d_pi0.GetXaxis().SetTitle("#pi^{0} Energy [MeV]")
h_sum_2d_pi0.GetYaxis().SetTitle("#pi^{0} Angle [degrees]")
h_sum_2d_pi0.SetStats(0)
h_sum_2d_pi0.Draw("COLZ")
c8.SaveAs("efficiency_map_e_pi0_kin.png")
myfile.WriteObject(c8,"pi0kin_map")

h_ee_2d_ph_l_clone.Divide(h_ee_2d_ph_l_tot)

c9 = TCanvas("c9")
h_ee_2d_ph_l_clone.GetXaxis().SetTitle("Leading Photon Energy [MeV]")
h_ee_2d_ph_l_clone.GetYaxis().SetTitle("Leading Photon Angle [degrees]")
h_ee_2d_ph_l_clone.SetStats(0)
h_ee_2d_ph_l_clone.Draw("COLZ")
c9.SaveAs("efficiency_map_ee_ph_l.png")

h_ee_2d_ph_sl_clone.Divide(h_ee_2d_ph_sl_tot)

c10 = TCanvas("c10")
h_ee_2d_ph_sl_clone.GetXaxis().SetTitle("Sub-leading Photon Energy [MeV]")
h_ee_2d_ph_sl_clone.GetYaxis().SetTitle("Sub-leading Photon Angle [degrees]")
h_ee_2d_ph_sl_clone.SetStats(0)
h_ee_2d_ph_sl_clone.Draw("COLZ")
c10.SaveAs("efficiency_map_ee_ph_sl.png")

h_hh_2d_ph_l.Divide(h_hh_2d_ph_l_tot)

c11 = TCanvas("c11")
h_hh_2d_ph_l.GetXaxis().SetTitle("Leading Photon Energy [MeV]")
h_hh_2d_ph_l.GetYaxis().SetTitle("Leading Photon Angle [degrees]")
h_hh_2d_ph_l.SetStats(0)
h_hh_2d_ph_l.Draw("COLZ")
c11.SaveAs("efficiency_map_hh_ph_l.png")

h_hh_2d_ph_sl.Divide(h_hh_2d_ph_sl_tot)

c12 = TCanvas("c12")
h_hh_2d_ph_sl.GetXaxis().SetTitle("Sub-leading Photon Energy [MeV]")
h_hh_2d_ph_sl.GetYaxis().SetTitle("Sub-leading Photon Angle [degrees]")
h_hh_2d_ph_sl.SetStats(0)
h_hh_2d_ph_sl.Draw("COLZ")
c12.SaveAs("efficiency_map_hh_ph_sl.png")

h_eh_2d_ph_l.Divide(h_eh_2d_ph_l_tot)

c13 = TCanvas("c13")
h_eh_2d_ph_l.GetXaxis().SetTitle("Leading Photon Energy [MeV]")
h_eh_2d_ph_l.GetYaxis().SetTitle("Leading Photon Angle [degrees]")
h_eh_2d_ph_l.SetStats(0)
h_eh_2d_ph_l.Draw("COLZ")
c13.SaveAs("efficiency_map_eh_ph_l.png")

h_eh_2d_ph_sl.Divide(h_eh_2d_ph_sl_tot)

c14 = TCanvas("c14")
h_eh_2d_ph_sl.GetXaxis().SetTitle("Sub-leading Photon Energy [MeV]")
h_eh_2d_ph_sl.GetYaxis().SetTitle("Sub-leading Photon Angle [degrees]")
h_eh_2d_ph_sl.SetStats(0)
h_eh_2d_ph_sl.Draw("COLZ")
c14.SaveAs("efficiency_map_eh_ph_sl.png")

h_ee_2d_ph_l_clone_clone.Divide(h_total_2d_ph_l)

c15 = TCanvas("c15")
h_ee_2d_ph_l_clone_clone.GetXaxis().SetTitle("Leading Photon Energy [MeV]")
h_ee_2d_ph_l_clone_clone.GetYaxis().SetTitle("Leading Photon Angle [degrees]")
h_ee_2d_ph_l_clone_clone.SetStats(0)
h_ee_2d_ph_l_clone_clone.Draw("COLZ")
c15.SaveAs("efficiency_map_ee_ph_l_allph.png")

h_ee_2d_ph_sl_clone_clone.Divide(h_total_2d_ph_sl)

c16 = TCanvas("c16")
h_ee_2d_ph_sl_clone_clone.GetXaxis().SetTitle("Sub-leading Photon Energy [MeV]")
h_ee_2d_ph_sl_clone_clone.GetYaxis().SetTitle("Sub-leading Photon Angle [degrees]")
h_ee_2d_ph_sl_clone_clone.SetStats(0)
h_ee_2d_ph_sl_clone_clone.Draw("COLZ")
c16.SaveAs("efficiency_map_ee_ph_sl_allph.png")

h_hh_2d_ph_l_clone.Divide(h_total_2d_ph_l)

c17 = TCanvas("c17")
h_hh_2d_ph_l_clone.GetXaxis().SetTitle("Leading Photon Energy [MeV]")
h_hh_2d_ph_l_clone.GetYaxis().SetTitle("Leading Photon Angle [degrees]")
h_hh_2d_ph_l_clone.SetStats(0)
h_hh_2d_ph_l_clone.Draw("COLZ")
c17.SaveAs("efficiency_map_hh_ph_l_allph.png")

h_hh_2d_ph_sl_clone.Divide(h_total_2d_ph_sl)

c18 = TCanvas("c18")
h_hh_2d_ph_sl_clone.GetXaxis().SetTitle("Sub-leading Photon Energy [MeV]")
h_hh_2d_ph_sl_clone.GetYaxis().SetTitle("Sub-leading Photon Angle [degrees]")
h_hh_2d_ph_sl_clone.SetStats(0)
h_hh_2d_ph_sl_clone.Draw("COLZ")
c18.SaveAs("efficiency_map_hh_ph_sl_allph.png")

h_eh_2d_ph_l_clone.Divide(h_total_2d_ph_l)

c19 = TCanvas("c19")
h_eh_2d_ph_l_clone.GetXaxis().SetTitle("Leading Photon Energy [MeV]")
h_eh_2d_ph_l_clone.GetYaxis().SetTitle("Leading Photon Angle [degrees]")
h_eh_2d_ph_l_clone.SetStats(0)
h_eh_2d_ph_l_clone.Draw("COLZ")
c19.SaveAs("efficiency_map_eh_ph_l_allph.png")

h_eh_2d_ph_sl_clone.Divide(h_total_2d_ph_sl)

c20 = TCanvas("c20")
h_eh_2d_ph_sl_clone.GetXaxis().SetTitle("Sub-leading Photon Energy [MeV]")
h_eh_2d_ph_sl_clone.GetYaxis().SetTitle("Sub-leading Photon Angle [degrees]")
h_eh_2d_ph_sl_clone.SetStats(0)
h_eh_2d_ph_sl_clone.Draw("COLZ")
c20.SaveAs("efficiency_map_eh_ph_sl_allph.png")

h_2d_pi0_ee.Divide(h_2d_pi0_ee_clone)

c21 = TCanvas("c21")
h_2d_pi0_ee.GetXaxis().SetTitle("#pi^{0} Energy [MeV]")
h_2d_pi0_ee.GetYaxis().SetTitle("#pi^{0} Angle [degrees]")
h_2d_pi0_ee.SetStats(0)
h_2d_pi0_ee.Draw("COLZ")
c21.SaveAs("efficiency_map_pi0_ee.png")

h_2d_pi0_hh.Divide(h_2d_pi0_hh_clone)

c22 = TCanvas("c22")
h_2d_pi0_hh.GetXaxis().SetTitle("#pi^{0} Energy [MeV]")
h_2d_pi0_hh.GetYaxis().SetTitle("#pi^{0} Angle [degrees]")
h_2d_pi0_hh.SetStats(0)
h_2d_pi0_hh.Draw("COLZ")
c22.SaveAs("efficiency_map_pi0_hh.png")

h_2d_pi0_eh.Divide(h_2d_pi0_eh_clone)

c23 = TCanvas("c23")
h_2d_pi0_eh.GetXaxis().SetTitle("#pi^{0} Energy [MeV]")
h_2d_pi0_eh.GetYaxis().SetTitle("#pi^{0} Angle [degrees]")
h_2d_pi0_eh.SetStats(0)
h_2d_pi0_eh.Draw("COLZ")
c23.SaveAs("efficiency_map_pi0_eh.png")
