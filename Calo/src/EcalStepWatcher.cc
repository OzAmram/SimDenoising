//custom headers
#include "SimDenoising/Calo/interface/EcalStepWatcher.h"

//CMSSW headers
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "SimDataFormats/RandomEngine/interface/RandomEngineState.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"

//Geant4 headers
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "Randomize.hh"

//STL headers
#include <algorithm>

typedef math::XYZTLorentzVector LorentzVector;

namespace {
//modifiable class w/ same layout as edm::StreamID
class FakeStreamID {
public:
	FakeStreamID(unsigned int value) : value_(value) {}
private:
	unsigned int value_;
};
}

EcalStepWatcher::EcalStepWatcher(const edm::ParameterSet& iConfig)
{
	//store list of volumes to watch
	const auto& vols = iConfig.getParameter<std::vector<std::string>>("volumes");
	volumes_.insert(vols.begin(),vols.end());
	
	//get parameters
	reset_random = iConfig.getParameter<bool>("reset_random");
	image_only = iConfig.getParameter<bool>("image_only");
	xbins = iConfig.getParameter<int>("xbins");
	ybins = iConfig.getParameter<int>("ybins");
	xmin = iConfig.getParameter<int>("xmin");
	xmax = iConfig.getParameter<int>("xmax");
	ymin = iConfig.getParameter<int>("ymin");
	ymax = iConfig.getParameter<int>("ymax");

	//create output tree
	tree_ = fs_->make<TTree>("tree","tree");

	//assign branches
	tree_->Branch("prim_pt",&entry_.prim_pt,"prim_pt/D");
	tree_->Branch("prim_eta",&entry_.prim_eta,"prim_eta/D");
	tree_->Branch("prim_phi",&entry_.prim_phi,"prim_phi/D");
	tree_->Branch("prim_E",&entry_.prim_E,"prim_E/D");
	tree_->Branch("prim_id",&entry_.prim_id,"prim_id/I");
	if (!image_only){
		tree_->Branch("step_t" , "vector<double>", &entry_.step_t, 32000, 0);
		tree_->Branch("step_x" , "vector<double>", &entry_.step_x, 32000, 0);
		tree_->Branch("step_y" , "vector<double>", &entry_.step_y, 32000, 0);
		tree_->Branch("step_z" , "vector<double>", &entry_.step_z, 32000, 0);
		tree_->Branch("step_E" , "vector<double>", &entry_.step_E, 32000, 0);
		tree_->Branch("step_t" , "vector<double>", &entry_.step_t, 32000, 0);

        tree_->Branch("t_avg_bin_weights", "vector<double>", &entry_.t_avg_bin_weights, 32000, 0);
        tree_->Branch("t_Eavg_bin_weights", "vector<double>", &entry_.t_Eavg_bin_weights, 32000, 0);
        tree_->Branch("t_max_bin_weights", "vector<double>", &entry_.t_max_bin_weights, 32000, 0);
        tree_->Branch("n_bin_weights", "vector<double>", &entry_.n_bin_weights, 32000, 0);
	}
    tree_->Branch("bin_weights", "vector<double>", &entry_.bin_weights, 32000, 0);
    tree_->Branch("xbins",&xbins,"xbins/I");
    tree_->Branch("ybins",&ybins,"ybins/I");
    tree_->Branch("xmin",&xmin,"xmin/I");
    tree_->Branch("xmax",&xmax,"xmax/I");
    tree_->Branch("ymin",&ymin,"ymin/I");
    tree_->Branch("ymax",&ymax,"ymax/I");
    h_E = new TH2F("h_E", "hist", xbins, xmin, xmax, ybins, ymin, ymax);
    h_t_avg = new TH2F("h_t_avg", "hist", xbins, xmin, xmax, ybins, ymin, ymax);
    h_t_Eavg = new TH2F("h_t_Eavg", "hist", xbins, xmin, xmax, ybins, ymin, ymax);
    h_t_max = new TH2F("h_t_max", "hist", xbins, xmin, xmax, ybins, ymin, ymax);
    h_n = new TH2F("h_n", "hist", xbins, xmin, xmax, ybins, ymin, ymax);
}

void EcalStepWatcher::update(const BeginOfEvent* evt) {  
	//reset branches
	entry_ = SimNtuple();
	h_E->Reset("ICESM");
	if (!image_only){
        h_t_avg->Reset("ICESM");
        h_t_Eavg->Reset("ICESM");
        h_t_max->Reset("ICESM");
        h_n->Reset("ICESM");

        entry_.step_x.clear();
        entry_.step_y.clear();
        entry_.step_z.clear();
        entry_.step_E.clear();
        entry_.step_E.clear();
        entry_.step_t.clear();
    }


	//reset random number generator
	if(reset_random){
		edm::Service<edm::RandomNumberGenerator> rng;
		//mockup of a stream ID: assume single thread
		FakeStreamID fid(0);
		edm::StreamID* sid(reinterpret_cast<edm::StreamID*>(&fid));
		//make a copy of initial cache
		if(orig_seeds.empty()) {
			std::vector<RandomEngineState> cache = rng->getEventCache(*sid);
			for(auto& state : cache){
				if(state.getLabel() == "g4SimHits"){
					//store original state
					orig_seeds = state.getSeed();
					return; //don't need to increment for first event
				}
			}
		}
		else {
			//increment all seeds
			std::for_each(orig_seeds.begin(), orig_seeds.end(), [](auto& n){ n++; });
			//reset G4 seed explicitly (also resets state consistently)
			G4Random::setTheSeed(orig_seeds[0]);
		}
	}
}

void EcalStepWatcher::update(const G4Step* step) {
	G4StepPoint* preStepPoint = step->GetPreStepPoint();
	const G4ThreeVector& hitPoint = preStepPoint->GetPosition();
	G4VPhysicalVolume* currentPV = preStepPoint->GetPhysicalVolume();
	std::string name(currentPV->GetName());	
	std::string subname(name.substr(0,4));
	if(volumes_.find(subname)==volumes_.end()) return;
	if (!image_only){
        auto time = step->GetTrack()->GetGlobalTime();
		entry_.step_x.push_back(hitPoint.x());
		entry_.step_y.push_back(hitPoint.y());
		entry_.step_z.push_back(hitPoint.z());
		entry_.step_E.push_back(step->GetTotalEnergyDeposit()); 
		entry_.step_t.push_back(time);

        h_n->Fill(hitPoint.x(), hitPoint.y(), 1);
        h_t_avg->Fill(hitPoint.x(), hitPoint.y(), time);
        h_t_Eavg->Fill(hitPoint.x(), hitPoint.y(), step->GetTotalEnergyDeposit() * time);

        //fill max of t and bin content
        auto bin = h_t_max->FindBin(hitPoint.x(), hitPoint.y());
        float cont = h_t_max->GetBinContent(bin);
        if( time > cont) h_t_max->SetBinContent(bin, time);

	}
	h_E->Fill(hitPoint.x(), hitPoint.y(), step->GetTotalEnergyDeposit());
}

void EcalStepWatcher::update(const EndOfEvent* evt) {

	//assume single particle gun
	G4PrimaryParticle* prim = (*evt)()->GetPrimaryVertex(0)->GetPrimary(0);
	LorentzVector vprim(prim->GetPx(),prim->GetPy(),prim->GetPz(),prim->GetTotalEnergy());

	entry_.prim_pt = vprim.pt();
	entry_.prim_eta = vprim.eta();
	entry_.prim_phi = vprim.phi();
	entry_.prim_E = vprim.energy();
	entry_.prim_id = prim->GetPDGcode();

    // get bin weights from TH2 and store in tree
    Int_t x, y;
    h_E->ClearUnderflowAndOverflow();
    entry_.bin_weights.reserve(xbins*ybins);
    for (x=1; x <= xbins; x++){
        for (y=1; y <= ybins; y++){
            entry_.bin_weights.push_back(h_E->GetBinContent(x, y));
        }
    }

    if(!image_only){
        entry_.t_avg_bin_weights.reserve(xbins*ybins);
        entry_.t_Eavg_bin_weights.reserve(xbins*ybins);
        entry_.t_max_bin_weights.reserve(xbins*ybins);
        entry_.n_bin_weights.reserve(xbins*ybins);
        for (x=1; x <= xbins; x++){
            for (y=1; y <= ybins; y++){
                if( h_n->GetBinContent(x,y) > 0) entry_.t_avg_bin_weights.push_back(h_t_avg->GetBinContent(x, y) / h_n->GetBinContent(x,y));
                else  entry_.t_avg_bin_weights.push_back(0.);
                if( h_E->GetBinContent(x,y) > 0.)  entry_.t_Eavg_bin_weights.push_back(h_t_Eavg->GetBinContent(x, y) / h_E->GetBinContent(x,y));
                else entry_.t_Eavg_bin_weights.push_back(0.);
                entry_.t_max_bin_weights.push_back(h_t_max->GetBinContent(x, y));
                entry_.n_bin_weights.push_back(h_n->GetBinContent(x, y));
            }
        }
    }


	//fill tree
	tree_->Fill();	
}

DEFINE_SIMWATCHER(EcalStepWatcher);
