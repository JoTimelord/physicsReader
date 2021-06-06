#include "C2V4p5.h"

int main(int argc, char** argv)
{
    C2V4p5();
    return 0;
}

void C2V4p5()
{

    // Open the output file
    TFile* output_file = new TFile("/home/jotimelord/physicsReader/VBSWWH_4p5.root", "RECREATE");

    // Create histograms for invariant masses & rapidity
    enum Objects {
        kW0 = 0,  // leading W
        kW1,      // subleading W
        kHiggs,   // higgs
        kB0,      // leading bottom quark
        kB1,      // subleading bottom quark
        kL0,      // leading charged lepton
        kL1,      // subleading charged lepton
        kVBF0,    // leading VBF quark
        kVBF1,    // subleading VBF quark
        kMET,     // neutrino+neutrino
        NObjects, // Number of objects
    };

    // ordering the array
    std::vector<TString> hist_types = {
        "leading W",
        "subleading W",
        "higgs",
        "leading bottom quark",
        "subleading bottom quark",
        "leading charged lepton",
        "subleading charged lepton",
        "leading VBF quark",
        "subleading VBF quark",
        "neutrino+neutrino"
    };

    // ordering the array
    std::vector<TString> hist_shortname_types = {
        "W0",
        "W1",
        "Higgs",
        "B0",
        "B1",
        "L0",
        "L1",
        "VBF0",
        "VBF1",
        "MET"
    };

    // histogram range
    std::vector<std::pair<float, float>> hist_pt_ranges = {
        {0., 1500.}, // W0
        {0., 1500.}, // W1
        {0., 1500.}, // Higgs
        {0.,  750.}, // B0
        {0.,  550.}, // B1
        {0.,  750.}, // L0
        {0.,  550.}, // L1
        {0.,  250.}, // VBF0
        {0.,  150.}, // VBF1
        {0.,  650.}, // MET
    };

    std::vector<std::pair<float, float>> hist_eta_ranges = {
        {-5, 5.}, // W0
        {-5, 5.}, // W1
        {-5, 5.}, // Higgs
        {-5, 5.}, // B0
        {-5, 5.}, // B1
        {-5, 5.}, // L0
        {-5, 5.}, // L1
        {-5, 5.}, // VBF0
        {-5, 5.}, // VBF1
        {-5, 5.}, // MET
    };

    std::vector<std::pair<float, float>> hist_phi_ranges = {
        {-TMath::Pi(), TMath::Pi()}, // W0
        {-TMath::Pi(), TMath::Pi()}, // W1
        {-TMath::Pi(), TMath::Pi()}, // Higgs
        {-TMath::Pi(), TMath::Pi()}, // B0
        {-TMath::Pi(), TMath::Pi()}, // B1
        {-TMath::Pi(), TMath::Pi()}, // L0
        {-TMath::Pi(), TMath::Pi()}, // L1
        {-TMath::Pi(), TMath::Pi()}, // VBF0
        {-TMath::Pi(), TMath::Pi()}, // VBF1
        {-TMath::Pi(), TMath::Pi()}, // MET
    };
    //Cross section value
    Double_t cs = 0.624 * 137;

    //Create histograms for invariant masses & rapidity
    TH1F* wBoson = new TH1F("wBoson", "Mass of W Bosons", 20000, 0, 300);
    TH1F* higg = new TH1F("higg", "Mass of Higgs", 10000, 0, 300);
    TH1F* rapidWlead = new TH1F("rapidWlead", "Rapidity of leading W Bosons relative to the Z axis", 10000, -100, 100);
    TH1F* rapidWsub = new TH1F("rapidWsub", "Rapidity of subleading W Bosons relative to the Z axis", 10000, -100, 100);
    TH1F* rapidHigg = new TH1F("rapidHigg", "Rapidity of Higgs relative to the Z axis", 10000, -100, 100);

    // Create an array of histograms
    std::vector<TH1F*> histograms_pt;
    std::vector<TH1F*> histograms_eta;
    std::vector<TH1F*> histograms_phi;

    for (unsigned int i = 0; i < NObjects; ++i)
    {
        // Binning into 1080 bins because 1080 contains a lot of prime factors, and therefore it's easy to rebin in the future. (cf. highly composite number)
        // It's not perfect but it gets the job done
        int nbins = 1080;
        histograms_pt.push_back(new TH1F(TString::Format("pt%s", hist_shortname_types[i].Data()), TString::Format("pt of %s", hist_types[i].Data()), nbins, hist_pt_ranges[i].first, hist_pt_ranges[i].second));
        histograms_eta.push_back(new TH1F(TString::Format("eta%s", hist_shortname_types[i].Data()), TString::Format("eta of %s", hist_types[i].Data()), nbins, hist_eta_ranges[i].first, hist_eta_ranges[i].second));
        histograms_phi.push_back(new TH1F(TString::Format("phi%s", hist_shortname_types[i].Data()), TString::Format("phi of %s", hist_types[i].Data()), nbins, hist_phi_ranges[i].first, hist_phi_ranges[i].second));
    }
    // =================================================================================================================
    // Do your stuff
    // -----------------------------------------------------------------------------------------------------------------

    // Open the root file
    TFile* f = TFile::Open("/home/jotimelord/physicsReader/VBSWWH_4p5_cmsgrid_final.root");

    // Simple pointer t
    TTree* t = (TTree*)f->Get("Physics");

    // Declare pointers;
    vector<Double_t>* mass = nullptr;
    vector<Double_t>* energy = nullptr;
    vector<Int_t>* pid = nullptr;
    vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>>* p4 = nullptr;

    // Access pid branches
    t->SetBranchAddress("PID", &pid);
    t->SetBranchAddress("M", &mass);
    t->SetBranchAddress("E", &energy);
    t->SetBranchAddress("P4", &p4);

    // Get entry number
    Long64_t entries = t->GetEntries();
    if (entries != 10000)
    {
        cout << "Event size is not equal to 10000." << endl;
    }

    //loop over all of the events
    for (Long64_t i = 0; i < entries; i++) {
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> ordered[NObjects];
        // =================================================================================================================
    //find all the bottom quarks
        t->GetEntry(i);
        vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>> bQ, Neu, Muon, qVBF;
        for (Long64_t j = 0; j < 8; j++) {
            Int_t pidj = pid->at(j);
            //Sort out different particles
            if (TMath::Abs(pidj) == 5) {
                bQ.push_back(p4->at(j));
            }
            if (TMath::Abs(pidj) == 11 || TMath::Abs(pidj) == 13 || TMath::Abs(pidj) == 15) {
                Muon.push_back(p4->at(j));
            }
            if (TMath::Abs(pidj) == 12 || TMath::Abs(pidj) == 14 || TMath::Abs(pidj) == 16) {
                Neu.push_back(p4->at(j));
            }
            else {
                qVBF.push_back(p4->at(j));
            }
        }
        ordered[kB0] = leading(bQ);
        ordered[kB1] = subleading(bQ);
        ordered[kL0] = leading(Muon);
        ordered[kL1] = subleading(Muon);
        ordered[kVBF0] = leading(qVBF);
        ordered[kVBF1] = subleading(qVBF);
        ordered[kMET] = Neu[0] + Neu[1];
        ordered[kHiggs] = bQ[0] + bQ[1];
        //reconstruct w Bosons
        vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>> wB;
        const auto mW1 = (Neu[0] + Muon[0]).M();
        const auto mW2 = (Neu[0] + Muon[1]).M(); 
        // std::cout << mW1 << endl;
        if (TMath::Abs(mW1 - 80) < TMath::Abs(mW2 - 80)) {
            wB.push_back((Neu[0] + Muon[0]));
            wB.push_back((Neu[1] + Muon[1]));
        }
        else {
            wB.push_back((Neu[0] + Muon[1]));
            wB.push_back((Neu[1] + Muon[0]));
        }
        ordered[kW0] = leading(wB);
        ordered[kW1] = subleading(wB);
        //fill histograms
        
        for (unsigned int n = 0; n < 10; n++) {
            histograms_pt[n]->Fill(ordered[n].Pt());
            histograms_eta[n]->Fill(ordered[n].Eta());
            histograms_phi[n]->Fill(ordered[n].Phi());
        } 
        //std::cout << ordered[kW0].M() << endl;
        wBoson->Fill(ordered[kW0].M());
        wBoson->Fill(ordered[kW1].M());
        higg->Fill(ordered[kHiggs].M());
        rapidWlead->Fill(ordered[kW0].Rapidity());
        rapidWsub->Fill(ordered[kW1].Rapidity());
        rapidHigg->Fill(ordered[kHiggs].Rapidity()); 
    }
     // Go to the output file and write the file
    output_file->cd(); // (TFile has a weird "ownership" that can be confusing...)
    for (Int_t i = 0; i < NObjects; ++i) {
        histograms_eta[i]->Scale(cs/(histograms_eta[i]->Integral()));
        histograms_pt[i]->Scale(cs/(histograms_pt[i]->Integral()));
        histograms_phi[i]->Scale(cs/(histograms_phi[i]->Integral()));
        histograms_eta[i]->Write();
        histograms_pt[i]->Write();
        histograms_phi[i]->Write();
    }
    wBoson->Scale(cs / (wBoson->Integral()));
    wBoson->Write();
    higg->Scale(cs / (higg->Integral()));
    higg->Write();
    rapidWlead->Scale(cs / (rapidWlead->Integral()));
    rapidWlead->Write();
    rapidWsub->Scale(cs / (rapidWsub->Integral()));
    rapidWsub->Write();
    rapidHigg->Scale(cs / (rapidHigg->Integral()));
    rapidHigg->Write();

    output_file->Close(); 
}

ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> leading(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>> v1) {
    if (v1[0].Pt() > v1[1].Pt())
    {
        return v1[0];
    }
    else
    {
        return v1[1];
    }
}

ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> subleading(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>> v1) {
    if (v1[0].Pt() > v1[1].Pt())
    {
        return v1[1];
    }
    else
    {
        return v1[0];
    }
}
