#include "TTree.h"
#include "TFile.h"

 #ifndef BooNeNtuples
 #define BooNeNtuples


namespace bsim
{
struct BooNENtuple
{
    float beamwgt;        /// Magic Weight: CrossSect * BooBeamNT
    int ntp;              /// 1,2,3,4: nue, nuebar, numu, numubar
    int npart;            /// number of particles in the chain
                          ///    npart-1 == proton
                          ///    0 == neutrino
    int id[20];           /// id of each particle in before chain, array length 'npart'
    float ini_pos[20][3]; /// 3-pos of particle in before chain, array length 'npart'
    float ini_mom[20][3]; /// 3-mom of particle in before chain, array length 'npart'
    float ini_eng[20];    /// E of particle in before chain, array length 'npart'
    float ini_t[20];      /// "decay" time of particle (starting from proton)
                          ///             in before chain, array length 'npart'
    float fin_mom[20][3]; /// final 3-mom of particle in before chain, array length 'npart'
    float fin_pol[20][3]; /// final 3-polarization of particle in before chain, array length 'npart'

    float run   = -9999.;
    float eventn= -9999.;
};

    class BooNe
    {
        TTree *myTree;
        TFile *myFile;
        bool debug = false;
        // bool debug = true;
        std::string filename;
        float POT   = 10000;
    public:
        BooNENtuple myNtuple;
        BooNe(){};
        BooNe(std::string Filename)
        {
            filename = Filename;
            std::cout << "loading file" << filename << std::endl;
            SetRootFile();
            std::cout << "Branches addres set" << std::endl;
            if (debug)
            {
                std::cout << "Printting first entry" << std::endl;
                PrintFirstEntry();
                std::cout << "POT:\t" << GetPOT() << std::endl;
            }
        };
        void GetEntry(int entry = 0)
        {
            myTree->GetEntry(entry);
        }

        float GetPOT()
        {//Returns the POT/entry (meant to be accumulated on each iteration)
            return POT / myTree->GetEntries();
        }
        
    private:
        void PrintFirstEntry()
        {
            std::cout << "beamwgt:\t" << myNtuple.beamwgt << std::endl;
            std::cout << "ntp:\t" << myNtuple.ntp << std::endl;
            std::cout << "npart:\t" << myNtuple.npart << std::endl;
            std::cout << "id:\t" << myNtuple.id[0] << std::endl;
            std::cout << "ini_pos:\t" << myNtuple.ini_pos[0][0] << std::endl;
            std::cout << "ini_mom:\t" << myNtuple.ini_mom[0][0] << std::endl;
            std::cout << "ini_eng:\t" << myNtuple.ini_eng[0] << std::endl;
            std::cout << "ini_t:\t" << myNtuple.ini_t[0] << std::endl;
            std::cout << "fin_mom:\t" << myNtuple.fin_mom[0][0] << std::endl;
            std::cout << "fin_pol:\t" << myNtuple.fin_pol[0][0] << std::endl;
        }

        void SetRootFile()
        {
            myFile = new TFile(filename.c_str());
            if (debug)
                myFile->ls();
            myTree = dynamic_cast<TTree *>(myFile->Get("h101"));

            myTree->SetBranchAddress("beamwgt", &myNtuple.beamwgt);
            myTree->SetBranchAddress("ntp", &myNtuple.ntp);
            myTree->SetBranchAddress("npart", &myNtuple.npart);
            myTree->SetBranchAddress("id", myNtuple.id);
            myTree->SetBranchAddress("ini_pos", &myNtuple.ini_pos[0][0]);
            myTree->SetBranchAddress("ini_mom", &myNtuple.ini_mom[0][0]);
            myTree->SetBranchAddress("ini_eng", myNtuple.ini_eng);
            myTree->SetBranchAddress("ini_t", myNtuple.ini_t);
            myTree->SetBranchAddress("fin_mom", &myNtuple.fin_mom[0][0]);
            myTree->SetBranchAddress("fin_pol", &myNtuple.fin_pol[0][0]);
            myTree->GetEntry(0);
        };
    };
}

 #endif