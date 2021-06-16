#include "C2V1.h"

int main(int argc, char** argv)
{
    C2V1();
    return 0;
}

void C2V1()
{

    // Open the output file
    // TFile* output_file = new TFile("/home/users/joytzphysics/plots/VBSWWH_1.root", "RECREATE");
    TFile* output_file = new TFile("/home/jotimelord/physicsReader/VBSWWH_1.root", "RECREATE");

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
    Double_t cs = 0.018 * 137;

    //Create histograms for invariant masses & rapidity
    TH1F* massW = new TH1F("massW", "Invariant mass of W Bosons", 20000, 0, 300);
    TH1F* massH = new TH1F("massH", "Invariant mass of Higgs", 10000, 100, 150);
    TH1F* rapidWlead = new TH1F("rapidWlead", "Rapidity of leading W Bosons relative to the Z axis", 10000, -100, 100);
    TH1F* rapidWsub = new TH1F("rapidWsub", "Rapidity of subleading W Bosons relative to the Z axis", 10000, -100, 100);
    TH1F* rapidHigg = new TH1F("rapidHigg", "Rapidity of Higgs relative to the Z axis", 10000, -100, 100);
    TH1F* deltaEta = new TH1F("deltaEta", "Delta Eta of VBF quarks", 10000, -10, 10);
    TH1F* massVBF = new TH1F("massVBF", "Invariant mass of VBF quarks", 20000, -5, 5);

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
    // TFile* f = TFile::Open("/home/users/joytzphysics/plots/VBSWWH_1_cmsgrid_final.root");
    TFile* f = TFile::Open("/home/jotimelord/physicsReader/VBSWWH_1_cmsgrid_final.root");

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

        if (TMath::Abs(mW1 - 80) < TMath::Abs(mW2 - 80)) {
            wB.push_back((Neu[0] + Muon[0]));
            wB.push_back((Neu[1] + Muon[1]));
            wB.push_back((Neu[0] + Muon[1]));
            wB.push_back((Neu[1] + Muon[0]));
        }
        else {
            wB.push_back((Neu[0] + Muon[1]));
            wB.push_back((Neu[1] + Muon[0]));
        }
        ordered[kW0] = leading(wB);
        ordered[kW1] = subleading(wB);

        //fill histograms
        for (unsigned int n = 0; n < NObjects; n++) {
            histograms_pt[n]->Fill(ordered[n].Pt());
            histograms_eta[n]->Fill(ordered[n].Eta());
            histograms_phi[n]->Fill(ordered[n].Phi());
        }
        //std::cout << ordered[kW0].M() << endl;
	deltaEta->Fill(TMath::Abs(ordered[kVBF0].Eta()-ordered[kVBF1].Eta()));
        rapidWlead->Fill(ordered[kW0].Rapidity());
        rapidWsub->Fill(ordered[kW1].Rapidity());
        rapidHigg->Fill(ordered[kHiggs].Rapidity()); 
	massW->Fill(ordered[kW0].M());
	massW->Fill(ordered[kW1].M());
	massH->Fill(ordered[kHiggs].M());
	massVBF->Fill(ordered[kVBF0].M());
	massVBF->Fill(ordered[kVBF1].M());
    }
    
     // Go to the output file and write the file
    output_file->cd(); // (TFile has a weird "ownership" that can be confusing...)
    for (Int_t i = 0; i < NObjects; ++i) {
        histograms_eta[i]->Scale(cs/(histograms_eta[i]->Integral()));
        histograms_pt[i]->Scale(cs/(histograms_pt[i]->Integral()));
        histograms_phi[i]->Scale(cs/ (histograms_phi[i]->Integral()));
        histograms_eta[i]->Write();
        histograms_pt[i]->Write();
        histograms_phi[i]->Write();
    }
    massW->Scale(cs / (massW->Integral()));
    massW->Write();
    massH->Scale(cs / (massH->Integral()));
    massH->Write();
    massVBF->Scale(cs / (massVBF->Integral()));
    massVBF->Write();
    rapidWlead->Scale(cs / (rapidWlead->Integral()));
    rapidWlead->Write();
    rapidWsub->Scale(cs / (rapidWsub->Integral()));
    rapidWsub->Write();
    rapidHigg->Scale(cs / (rapidHigg->Integral()));
    rapidHigg->Write();
    deltaEta->Scale(cs / (deltaEta->Integral()));
    deltaEta->Write();
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
