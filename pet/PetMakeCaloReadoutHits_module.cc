//
// An EDProducer Module that reads StepPointMC objects and turns them into
// CaloHit, CaloHitMCTruth, CaloHitMCPtr, CrystalOnlyHit, CrystalHitMCPtr
// objects.
//
// Original author Ivan Logashenko; major changes by Krzysztof Genser and Rob Kutschke.
//
// Notes
// 1) The CrystalOnlyHitCollection is a form of MC truth on a per crystal basis.
//    It represents an idea per crystal response if no readouts were hit directly.
//
// 2) We still have two times the crystal StepPointMC hits saved per hit in a crystal. Note sure we can 
//    keep only a single copy before simulating the whole readout chain, so keep them for now. 


// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/sort_functors.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"

#include "murat/pet/PetGeometryService.hh"
#include "murat/pet/PetGeomHandle.hh"
#include "murat/pet/BrainImager.hh"

#include "HitMakers/inc/CaloReadoutUtilities.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "SeedService/inc/SeedService.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"


using namespace std;

namespace mu2e {

  // Anonymous namespace to hold some helper classes.
  namespace {


       // A helper class to hold temporary information.
       class ROHit {
	   
	   public:

	     // Is this StepPointMC from a hit in the crystal or a hit in the readout device?
	     enum step_type {crystal, readout};

	     art::Ptr<StepPointMC> _step;    // Ptr back to the StepPointMC.
	     double    _edep;                // copy of energy from the StepPointMC.
	     double    _edep_corr;           // Energy corrected for absorption within the crystal.
	     step_type _type;                // Is this a hit in the crystal or the readout?
	     double    _time;                // copy of the time from the StepPointMC.

	     ROHit(art::Ptr<StepPointMC> const& step, double edep, double edep1, step_type charged, double time):
               _step(step), _edep(edep), _edep_corr(edep1),
               _type(charged), _time(time) { }

	     // This operator is overloaded in order to time-sort the hits
	     bool operator <(const ROHit& b) const { return (_time < b._time); }

       }; 


       // A helper class to add Ptr's to the appropriate PtrStepPointMCVector collection.
       class PtrAdder{

	   public:
	      PtrAdder( PtrStepPointMCVector& crystals, 
                	PtrStepPointMCVector& readouts ) : 
			_crystals(crystals),_readouts(readouts){}

	      // Choose the appropriate collection and add the Ptr for this hit.
	      void operator()( ROHit const& hit )
	      {
        	if ( hit._type == ROHit::crystal) _crystals.push_back( hit._step );
        	else                              _readouts.push_back( hit._step );
	      }


	   private:

	     PtrStepPointMCVector& _crystals;
	     PtrStepPointMCVector& _readouts;

       }; 


  } // end anonymous namespace





  //--------------------------------------------------------------------
  //
  //
  class PetMakeCaloReadoutHits : public art::EDProducer {
  
     public:

       // First vector is list of crystal steps, associated with particular readout element.
       // Second vector is list of readout steps, associated with particular readout element.

       explicit PetMakeCaloReadoutHits(fhicl::ParameterSet const& pset) :

	 // Parameters
	 _diagLevel(pset.get<int>("diagLevel",0)),
	 _maxFullPrint(pset.get<int>("maxFullPrint",5)),
	 _fillDetailedHit(pset.get<bool>("fillDetailedHit",0)),
	 _caloLRUcorrection(pset.get<bool>("caloLRUcorrection",0)),
	 _caloNonLinCorrection(pset.get<bool>("caloNonLinCorrection",0)),
	 _randGauss( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
	 _stepPoints(pset.get<string>("calorimeterStepPoints","calorimeter")),
	 _rostepPoints(pset.get<string>("calorimeterROStepPoints","calorimeterRO")),
	 _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
	 _messageCategory("CaloReadoutHitsMakerNew"){

	 // Tell the framework what we make.
	 produces<CaloHitCollection>();
	 produces<CaloHitMCTruthCollection>();
	 produces<CaloCrystalOnlyHitCollection>();
	 produces<PtrStepPointMCVectorCollection>("CaloHitMCCrystalPtr");
	 produces<PtrStepPointMCVectorCollection>("CaloHitMCReadoutPtr");
	 produces<CaloHitSimPartMCCollection>();

       }

       virtual ~PetMakeCaloReadoutHits() { }
       virtual void beginJob();
       void produce( art::Event& e);


     private:

	 typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
	 typedef art::Ptr<StepPointMC>   StepPtr;
	 typedef std::vector<StepPtr >   StepPtrs;
	 typedef art::Ptr<SimParticle>   SimPtr;
	 typedef std::vector<SimPtr >    SimPtrs;
	 typedef std::map<int,StepPtrs > HitMap;

	 
	 int _diagLevel;  
	 int _maxFullPrint;
         bool _fillDetailedHit;
         bool _caloLRUcorrection;
	 bool _caloNonLinCorrection;
         CLHEP::RandGaussQ _randGauss;

	 std::string _stepPoints;
	 std::string _rostepPoints;
	 string _g4ModuleLabel;  // Name of the module that made these hits.



	 const std::string _messageCategory;


	 void makeCalorimeterHits (HandleVector const& crystalStepsHandles,
                        	   HandleVector const& readoutStepsHandles,
                        	   CaloHitCollection &,
                        	   CaloHitMCTruthCollection&,
                        	   CaloCrystalOnlyHitCollection&,
                        	   PtrStepPointMCVectorCollection&,
                        	   PtrStepPointMCVectorCollection&,
				   CaloHitSimPartMCCollection& );


	 void nonLinearityCorrection(double kinetic_energy, double energy, int cryId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations );
	 void longitudinalResponseUniformityCorrection(double posZ, double cryhalflength, double energy, int crid,ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
	 void fillMapById(HitMap& map,HandleVector const& crystalStepsHandles);


	 // Print information about the data products found by the selector functions. 
	 void printDataProductInfo( HandleVector const& crystalStepsHandles,
                        	    HandleVector const& readoutStepsHandles );


  };

  //--------------------------------------------------------------------



  void PetMakeCaloReadoutHits::beginJob(){}


  void  PetMakeCaloReadoutHits::produce(art::Event& event) {

    if ( _diagLevel > 0 ) {
      cout << "PetMakeCaloReadoutHits: produce() begin" << endl;
      printf(" Run:Event: %i:%i\n",event.run(),event.event());
    }

    static int ncalls(0);
    ++ncalls;

    // Check that calorimeter geometry description exists
    art::ServiceHandle<PetGeometryService> geom;    
    if( !(geom->hasElement<BrainImager>()) ) return;
   

    
    // A container to hold the output hits.
    unique_ptr<CaloHitCollection>               caloHits         (new CaloHitCollection);
    unique_ptr<CaloHitMCTruthCollection>        caloMCHits       (new CaloHitMCTruthCollection);
    unique_ptr<CaloCrystalOnlyHitCollection>    caloCrystalMCHits(new CaloCrystalOnlyHitCollection);
    unique_ptr<PtrStepPointMCVectorCollection>  caloMCptrHits    (new PtrStepPointMCVectorCollection);
    unique_ptr<PtrStepPointMCVectorCollection>  caloMCroptrHits  (new PtrStepPointMCVectorCollection);
    unique_ptr<CaloHitSimPartMCCollection>      caloMCSimParts   (new CaloHitSimPartMCCollection);


    // These selectors will select data products with the given instance name, and ignore
    // all other fields of the product ID.
    art::ProductInstanceNameSelector getCrystalSteps(_stepPoints);
    art::ProductInstanceNameSelector getReadoutSteps(_rostepPoints);

    // Get the StepPointMCs from the event.
    HandleVector crystalStepsHandles, readoutStepsHandles;
    event.getMany( getCrystalSteps, crystalStepsHandles);
    event.getMany( getReadoutSteps, readoutStepsHandles);


    static bool firstEvent(true);
    if ( firstEvent ) {
      printDataProductInfo( crystalStepsHandles, readoutStepsHandles);
      firstEvent = false;
    }


    makeCalorimeterHits(crystalStepsHandles, readoutStepsHandles,
                        *caloHits, *caloMCHits, *caloCrystalMCHits,
                        *caloMCptrHits, *caloMCroptrHits, *caloMCSimParts);


    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      cout << "PetMakeCaloReadoutHits: Total number of calorimeter hits = "
           << caloHits->size()
           << endl;
      cout << "PetMakeCaloReadoutHits: Total number of crystal MC hits = "
           << caloCrystalMCHits->size()
           << endl;
    }

    
    
    // Add the output hit collection to the event
    event.put(std::move(caloHits));
    event.put(std::move(caloMCHits));
    event.put(std::move(caloCrystalMCHits));
    event.put(std::move(caloMCptrHits),"CaloHitMCCrystalPtr");
    event.put(std::move(caloMCroptrHits),"CaloHitMCReadoutPtr");
    event.put(std::move(caloMCSimParts));

    if ( _diagLevel > 0 ) cout << "PetMakeCaloReadoutHits: produce() end" << endl;

  } 













  void PetMakeCaloReadoutHits::makeCalorimeterHits (HandleVector const& crystalStepsHandles,
                                                 HandleVector const& readoutStepsHandles,
                                                 CaloHitCollection &caloHits,
                                                 CaloHitMCTruthCollection& caloHitsMCTruth,
                                                 CaloCrystalOnlyHitCollection& caloCrystalHitsMCTruth,
                                                 PtrStepPointMCVectorCollection& caloHitsMCCrystalPtr,
                                                 PtrStepPointMCVectorCollection& caloHitsMCReadoutPtr,
						 CaloHitSimPartMCCollection& caloMCSimParts ){
    
    
    
    GlobalConstantsHandle<ParticleDataTable> pdt;
    // Handle to the conditions service
    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
    
    BrainImager const & cal = *(PetGeomHandle<BrainImager>());
    double timeGap    = cal.getTimeGap();
    double addEdep    = cal.getElectronEdep();
    double cryhalflength = cal.crystalHalfLength();



    // Fill map: crystal id -> StepMCPtr and RO id -> StepMCPtr
    HitMap hitmapCrystal,hitmapRO;        
    fillMapById( hitmapCrystal, crystalStepsHandles);
    fillMapById( hitmapRO, readoutStepsHandles);

    CaloReadoutUtilities readoutUtil;


    //loop over each crystal, collect hits in the crystal, then create RO hits
    for(HitMap::const_iterator crIter = hitmapCrystal.begin(); crIter != hitmapCrystal.end(); ++crIter ) {
      
	
	vector<ROHit> cr_hits;
	int crid = crIter->first;
	StepPtrs const& isteps = crIter->second;

	// Loop over steps inside the crystal for a given crystal id
	for ( StepPtrs::const_iterator i=isteps.begin(), e=isteps.end(); i != e; ++i ){

            StepPointMC const& h = **i;

            if ( h.eDep()<=0.0 ) continue;
	    double edep_corr = h.eDep();

	    // Calculate correction for edep
	    CLHEP::Hep3Vector const& posInMu2e = h.position();
	    CLHEP::Hep3Vector posInCrystal     = cal.toCrystalFrame(crid,posInMu2e);

//must check if this works for disks	
double posZ = posInCrystal.x();

	    if (_caloNonLinCorrection && h.simParticle().isNonnull())
	    {

		 double edep_save(edep_corr);

		 art::Ptr<SimParticle> const& simptr = h.simParticle();
		 SimParticle const& sim = *simptr;

		 if( (std::abs(sim.pdgId()) == 11 || sim.pdgId() == 22))
		 {
		     const HepPDT::ParticleData& data = pdt->particle(sim.pdgId()).ref();
		     double mass = data.mass().value();

		     //add non-linearity effect
		     double trackKine = std::sqrt(h.momentum().mag2() + std::pow(mass, 2) ) - mass;

		     if (trackKine < 1.0) nonLinearityCorrection(trackKine, edep_corr, crid, calorimeterCalibrations);

		     if(_diagLevel > 0)		    
		       std::cout<<"************************** BEFORE ? AFTER NON-LINEARITY EFFECT-> edep_corr = "<< edep_save<<"  /  "<<edep_corr<<std::endl
		        	<<", energyKin = "<< trackKine<<", mass = "<< mass<< ", momentum.mag2() = "<< h.momentum().mag2()<<std::endl;
		 } 

	    }


	    if (_caloLRUcorrection)
	    {
		double edep_save(edep_corr);
		longitudinalResponseUniformityCorrection(posZ, cryhalflength, edep_corr,crid, calorimeterCalibrations);
		if (_diagLevel > 0) std::cout<<"***************BEFORE /  AFTER LRU EFFECT-> edep_corr = "<< edep_save<<"  /  "<<edep_corr<<std::endl;	  
	    }

            cr_hits.push_back(ROHit(*i,h.eDep(),edep_corr,ROHit::crystal,h.time()));
	}


        


	
	//loop over all RO for a given crystal id (Roid = Roidbase + j, j=0,nROPerCrystal-1)
	//for each readout, assign hits in the crystal + hits in the readout
	int ROidBase = cal.ROBaseByCrystal(crid);
	int ROidEnd  = ROidBase+cal.nROPerCrystal();

	for (int roid=ROidBase;roid<ROidEnd;++roid)
	{

	     vector<ROHit> ro_hits(cr_hits);

	     //find the entry in RO map and add the RO hits
	     HitMap::const_iterator irIter = hitmapRO.find(roid);

	     if (irIter != hitmapRO.end() ){

		StepPtrs const& irosteps = irIter->second;
		for( StepPtrs::const_iterator i=irosteps.begin(); i!=irosteps.end(); ++i ){

        	    StepPointMC const& h = **i;
        	    // There is no cut on energy deposition here - may be, we need to add one?
        	    ro_hits.push_back(ROHit(*i,0.,0.,ROHit::readout,h.time()));
		}
             }


	     // Sort hits by time
	     sort(ro_hits.begin(), ro_hits.end());


	     // A buffer to collect output.
	     PtrStepPointMCVector mcptr_crystal;
	     PtrStepPointMCVector mcptr_readout;
	     
	     // A tool to add the art::Ptr to the correct output collection.
	     PtrAdder addPtr( mcptr_crystal, mcptr_readout );


	     // Loop over sorted hits and form complete ro/calorimeter hits
	     // We will need to digitize the hits here, simulate APD response

	     double h_time    = ro_hits[0]._time;
	     double h_edep    = ro_hits[0]._edep;
	     double h_edepc   = ro_hits[0]._edep_corr;
	     ROHit::step_type h_type = ro_hits[0]._type;
	     addPtr( ro_hits[0] );

	     for( size_t i=1; i<ro_hits.size(); ++i ) {

        	 if( (ro_hits[i]._time- h_time) > timeGap ) {

		   // Save current hit
		   CaloHitSimPartMC  caloHitSimPartMC;
		   if (_fillDetailedHit) readoutUtil.fillSimMother(cal,mcptr_crystal,caloHitSimPartMC);
		   
		   caloHits.push_back(       CaloHit(       roid,h_time,h_edepc+h_type*addEdep));
        	   caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_type));
        	   caloHitsMCCrystalPtr.push_back(mcptr_crystal);
        	   caloHitsMCReadoutPtr.push_back(mcptr_readout);
                   caloMCSimParts.push_back(caloHitSimPartMC);	   
        	   
		   // ...and create new hit
        	   mcptr_crystal.clear();
        	   mcptr_readout.clear();
        	   addPtr( ro_hits[i] );

        	   h_time    = ro_hits[i]._time;
        	   h_edep    = ro_hits[i]._edep;
        	   h_edepc   = ro_hits[i]._edep_corr;
        	   h_type    = ro_hits[i]._type;

        	 } else {
        	   // Append data to hit
        	   h_edep  += ro_hits[i]._edep;
        	   h_edepc += ro_hits[i]._edep_corr;
        	   if( ro_hits[i]._type != ROHit::crystal ) h_type = ROHit::readout; // this does not count the charge...
        	   addPtr( ro_hits[i] );
        	 }
	     }

	     CaloHitSimPartMC  caloHitSimPartMC;
	     if (_fillDetailedHit) readoutUtil.fillSimMother(cal,mcptr_crystal,caloHitSimPartMC);
	     
	     caloHits.push_back(       CaloHit(roid,h_time,h_edepc+h_type*addEdep));
	     caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_type));
	     caloHitsMCCrystalPtr.push_back(mcptr_crystal);
	     caloHitsMCReadoutPtr.push_back(mcptr_readout);
             caloMCSimParts.push_back(caloHitSimPartMC);	   

	}

    }







    for(HitMap::const_iterator cr = hitmapCrystal.begin(); cr != hitmapCrystal.end(); ++cr) {

	CaloCrystalOnlyHitCollection cr_hits;

	int cid = cr->first;
	for( size_t i=0; i<cr->second.size(); i++ ) {
           StepPointMC const& h2 = *cr->second[i];
           cr_hits.push_back(CaloCrystalOnlyHit(cid,h2.time(),h2.eDep()));
	}

	sort(cr_hits.begin(), cr_hits.end(), lessByTime<CaloCrystalOnlyHitCollection::value_type>());


	// now form final hits if they are close enough in time
	CaloCrystalOnlyHitCollection::value_type cHitMCTruth  = cr_hits[0];

	for ( size_t i=1; i<cr_hits.size(); ++i ) {

          if ( (cr_hits[i].time()-cr_hits[i-1].time()) > timeGap ) {

            // Save current hit and create new onw
            caloCrystalHitsMCTruth.push_back(cHitMCTruth);
            cHitMCTruth  = cr_hits[i];

          } else {

            // Add energy to the old hit (keep the "earlier" time)
            cHitMCTruth.setEnergyDep(cHitMCTruth.energyDep()+cr_hits[i].energyDep());

          }

	}

	caloCrystalHitsMCTruth.push_back(cHitMCTruth);

      }





  } // end makeCalorimeterHits







  void PetMakeCaloReadoutHits::fillMapById( HitMap& hitmap,HandleVector const& crystalStepsHandles)
  {
      for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end(); i != e; ++i ){

	art::Handle<StepPointMCCollection> const& handle(*i);
	StepPointMCCollection const& steps(*handle);

	StepPointMCCollection::const_iterator j0=steps.begin();
	for ( StepPointMCCollection::const_iterator j=j0, je=steps.end(); j != je; ++j ){
          StepPointMC const& step(*j);
          hitmap[step.volumeId()].push_back( StepPtr(handle,j-j0) );
	}
      }

  } 



  void PetMakeCaloReadoutHits::printDataProductInfo( HandleVector const& crystalStepsHandles,
                                                  HandleVector const& readoutStepsHandles ){
      mf::LogInfo log(_messageCategory);
      log << "MakeCaloReadoutHit::produce will use StepPointMCs from: \n";
      for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end();
            i != e; ++i ){
	art::Provenance const& prov(*(i->provenance()));
	log  << "   " << prov.branchName() << "\n";
      }
      log << "\nMakeCaloReadoutHit::produce will use StepPointMCs from: \n";
      for ( HandleVector::const_iterator i=readoutStepsHandles.begin(), e=readoutStepsHandles.end();
            i != e; ++i ){
	art::Provenance const& prov(*(i->provenance()));
	log  << "   " << prov.branchName() << "\n";
      }
  } 
  
  
  void PetMakeCaloReadoutHits::nonLinearityCorrection(double kinetic_energy, double energy, int cryId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations ){
  
		double MeV2keV = 1000.0;
				
		//the formula used returns positive values if tmpEnergy>(approximately) 3keV, see Mu2e doc 1748-v1
		kinetic_energy = MeV2keV*kinetic_energy;
		double trackKine = kinetic_energy;
		kinetic_energy = std::log10(kinetic_energy);
		kinetic_energy *= _randGauss.fire(calorimeterCalibrations->LINpar2(cryId), calorimeterCalibrations->LINpar2Err(cryId) );
		kinetic_energy += _randGauss.fire(calorimeterCalibrations->LINpar1(cryId), calorimeterCalibrations->LINpar1Err(cryId) );
		kinetic_energy = std::log10(kinetic_energy);
		kinetic_energy *= _randGauss.fire(calorimeterCalibrations->LINpar0(cryId), calorimeterCalibrations->LINpar0Err(cryId) );
		kinetic_energy /= _randGauss.fire(calorimeterCalibrations->LINpar3(cryId), calorimeterCalibrations->LINpar3Err(cryId) );
		kinetic_energy /= std::log10(trackKine);
		
		if (kinetic_energy>0) energy *= kinetic_energy;
		
  }
  
  void PetMakeCaloReadoutHits::longitudinalResponseUniformityCorrection(double posZ, double cryhalflength, double energy, int crid, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations ){
     posZ *= -1.0;
     posZ += cryhalflength;
     posZ *= -_randGauss.fire(calorimeterCalibrations->LRUpar0(crid), calorimeterCalibrations->LRUpar0Err(crid) );
     posZ += 1.0;
     energy *= posZ;
  }


} // end namespace mu2e

using mu2e::PetMakeCaloReadoutHits;
DEFINE_ART_MODULE(PetMakeCaloReadoutHits);



