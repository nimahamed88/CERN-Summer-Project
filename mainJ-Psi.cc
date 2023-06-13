#include "Pythia8/Pythia.h"
#include "TH1F.h"
#include "TFile.h"

using namespace Pythia8;
int main() {

  // You can always read an plain LHE file,
  // but if you ran "./configure --with-gzip" before "make"
  // then you can also read a gzipped LHE file.
#ifdef GZIP
  bool useGzip = true;
#else
  bool useGzip = false;
#endif
  cout << " useGzip = " << useGzip << endl;

  // Generator. We here stick with default values, but changes
  // could be inserted with readString or readFile.
  Pythia pythia;

  // Initialize Les Houches Event File run. List initialization information.
  pythia.readString("Beams:frameType = 4");
  if (useGzip) pythia.readString("Beams:LHEF = pp_ttxbbx_NLO_MG5_PY8.lhe.gz");
  else         pythia.readString("Beams:LHEF = pp_ttxbbx_NLO_MG5_PY8.lhe");




const int B0_ID   = 511;
  const int B0s_ID  = 531;
  const int Jpsi_ID = 443;
  const int K1_1270_ID = 10313;
  const int Kstar892_ID = 323;
  const int phi_ID  = 333;
  const int K0s_ID  = 310;
  const int Kch_ID  = 321;
  const int pich_ID = 211;
  const int mu_ID = 13;
  const int Dch_ID = 411;
  const int D0_ID = 421;
  const int Ds_ch_ID = 431;



  // Switch-off MPI and others to keep it simple. From pythia team
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:FSR = off");
  pythia.readString("BeamRemnants:primordialKt = off");

  // pythia.readString("6:onMode = off"); // no decay for top quark


  pythia.readString("443:onMode = off");  // All J/psi decays off
  pythia.readString("443:onIfAny = 13"); // J/psi->mu+mu- allowed



  pythia.readString("511:onMode = off");  // All B0 decay off
  pythia.readString("511:onIfAny = 443"); // B0 -> J/psi+X allowed

  pythia.readString("521:onMode = off");  // All B+ decay off
  pythia.readString("521:onIfAny = 443"); // B+ -> J/psi+X allowed

  pythia.readString("-521:onMode = off");  // All B- decay off
  pythia.readString("-521:onIfAny = 443"); // B- -> J/psi+X allowed

  pythia.readString("531:onMode = off");  // All B_s decay off
  pythia.readString("531:onIfAny = 443"); // B_s -> J/psi+X allowed

  pythia.readString("541:onMode = off");  // All B_c decay off
  pythia.readString("541:onIfAny = 443"); // B_c -> J/psi+X allowed


  pythia.readString("411:onMode = off");  // All D+/- decays off
  pythia.readString("421:onMode = off");  // All D0 decays off
  pythia.readString("431:onMode = off");  // All Ds+/- decays off


  pythia.init();


// Histograms
  TH1F* hPtJPsi = new TH1F("hPtJPsi", "Transverse Momentum of J/Psi", 100, 0, 100);
  TH1F* hYJPsi = new TH1F("hYJPsi", "Rapidity of J/Psi", 100, -10, 10);
  TH1F* hEtaJPsi = new TH1F("hEtaJPsi", "Pseudorapidity of J/Psi", 100, -10, 10);
  TH1F* hTransverseVertex = new TH1F("hTransverseVertex", "Transverse Vertex Distribution; Transverse Vertex (µm); Counts", 100, 0, 1000);




 // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;
  int nBmesons = 0;        // Counter for B mesons
  int nBtoJPsi = 0;       // Counter for B mesons decaying to J/Psi
  int nBtoJPsiAfterCuts = 0; // Counter for B mesons decaying to J/Psi after cuts
  int nJPsi = 0;
  int nJPsiPassingCuts = 0;     // Counter for J/Psi particles passing cuts

  int nJPsiToMuMu = 0;                    // Counter for J/Psi decays to muon-muon
  int nJPsiToMuMuPassingCuts = 0;         // Counter for J/Psi decays to muon-muon passing cuts



    // Event loop
    for (int iEvent = 0;  ; ++iEvent) {
           // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }


    for (int i = 0; i < pythia.event.size(); ++i) {
     
      if (pythia.event[i].id() == 443) //&& pythia.event[i].isFinal()) {
                 {

                  ++nJPsi;
                  // Get the transverse momentum (Pt) of the J/Psi
                  double pt = pythia.event[i].pT();
                  hPtJPsi->Fill(pt);
                  double Y = pythia.event[i].y();
                  hYJPsi->Fill(Y);
                  double eta = pythia.event[i].eta();
                  hEtaJPsi->Fill(eta);
                  // Get vertex position of J/Psi in millimeters
                double x = pythia.event[i].xProd();
                double y = pythia.event[i].yProd();

                // Calculate transverse vertex
                double L = sqrt(x * x + y * y);

                // Fill histogram with transverse vertex values
                hTransverseVertex->Fill(L * 1000); // Convert from mm to µm
                if (pt > 6 && abs(Y) < 2.4 && abs(L) < 0.03) {
                ++nJPsiPassingCuts;
                int muon1 = -1;
            int muon2 = -1;
             // Loop over J/Psi daughters to find the muon-muon pair
            for (int daughter : pythia.event[i].daughterList()) {
                if (abs(pythia.event[daughter].id()) == 13) {
                    if (muon1 == -1)
                        muon1 = daughter;
                    else if (muon2 == -1)
                        muon2 = daughter;
                    else
                        break;  // Already found both muons
                    }
                }
                // Check if two muons were found
    if (muon1 != -1 && muon2 != -1) {
        // Calculate the invariant mass of the muon-muon pair
        double massMuMu = (pythia.event[muon1].p() + pythia.event[muon2].p()).mCalc();

        // Apply the specified cuts on the muon-muon system
        double ptMuon1 = pythia.event[muon1].pT();
        double ptMuon2 = pythia.event[muon2].pT();
        double etaMuon1 = pythia.event[muon1].eta();
        double etaMuon2 = pythia.event[muon2].eta();

        if (2.9 < massMuMu && massMuMu < 3.3 && ((abs(etaMuon1) < 1.2 && ptMuon1 > 3.5) || (1.2 < abs(etaMuon1) && abs(etaMuon1) < 2.4 && ptMuon1 > 2.5)) && ((abs(etaMuon2) < 1.2 && ptMuon2 > 3.5) || (1.2 < abs(etaMuon2) && abs(etaMuon2) < 2.4 && ptMuon2 > 2.5))) {
            ++nJPsiToMuMuPassingCuts;
        }

        ++nJPsiToMuMu;
    }
            }
            



      }



            // Check if the particle is a B meson
      if (abs(pythia.event[i].id()) == 511 || abs(pythia.event[i].id()) == 521 || abs(pythia.event[i].id()) == 531 || abs(pythia.event[i].id()) == 541) {
          ++nBmesons;

          // Check if the B meson decays to J/Psi
            bool decaysToJPsi = false;
            for (int daughter : pythia.event[i].daughterList()) {
                if (abs(pythia.event[daughter].id()) == 443) {
                    decaysToJPsi = true;
                    break;
                }
            }

            if (decaysToJPsi) {
                ++nBtoJPsi;
                // Fill histograms for non-prompt J/Psi properties
                int jpsiIndex = pythia.event[i].daughterList().at(0); // Assuming J/Psi is the first daughter
                if (!pythia.event[jpsiIndex].isFinal()) continue;      // Skip if J/Psi is not final

                // Get the properties of non-prompt J/Psi
                 


              }

                // // Apply cuts on J/Psi properties
                // int jpsiIndex = pythia.event[i].daughterList().at(0); // Assuming J/Psi is the first daughter
                // if (pythia.event[jpsiIndex].pT() > 6 && abs(pythia.event[jpsiIndex].eta()) < 2.4) {
                //     ++nBtoJPsiAfterCuts;
                // }
            
        }



}



// end loop
    }



cout << "Number of total J/Psi mesons: " << nJPsi<<endl;
cout << "Number of B mesons: " << nBmesons<<endl;
cout << "Number of B to J/Psi: " << nBtoJPsi<<endl;


// Calculate the fraction of B mesons decaying to J/Psi
double fractionBtoJPsi = static_cast<double>(nBtoJPsi) / nBmesons;
cout << "Fraction of B mesons decaying to J/Psi: " << fractionBtoJPsi << endl;


// Calculate the fraction of J/Psi particles passing the cuts
double fractionJPsiPassingCuts = static_cast<double>(nJPsiPassingCuts) / nJPsi;
cout << "Fraction of J/Psi particles passing the initial cuts (only for J/Psi): " << fractionJPsiPassingCuts << endl;


double fractionJPsiToMuMuPassingCuts = static_cast<double>(nJPsiToMuMuPassingCuts) / nJPsi;
cout << "Number of J/Psi (passed the cuts) to Mu Mu: " <<nJPsiToMuMu<<endl;
cout << "Number of J/Psi to Mu Mu passing all cuts (included muon cuts): " <<nJPsiToMuMuPassingCuts<<endl;
cout << "Fraction of J/Psi decays to muon-muon passing all cuts: " << fractionJPsiToMuMuPassingCuts << endl;




  // Save histograms to a ROOT file
  TFile outputFile("output.root", "RECREATE");
  hPtJPsi->Write();
  hEtaJPsi->Write();
  hYJPsi->Write();
 
  hTransverseVertex->Write();
  outputFile.Close();

    // // Delete histograms
    delete hPtJPsi;
    delete hYJPsi;
    // delete hetaJPsi;
    delete hTransverseVertex;

    // End Pythia
    pythia.stat();

    // Done
    return 0;
  
}