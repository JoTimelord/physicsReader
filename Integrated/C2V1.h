#ifndef C2V1_h
#define C2V1_h

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "Math/LorentzVector.h"
#include <Math/GenVector/GenVector_exception.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzE4D.h>
#include <array>
#include <iostream>
#include <iterator>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

// Declare helper functions
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> leading(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>> v1);
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> subleading(vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>> v1);

int main(int argc, char** argv);
void C2V1();

#endif
