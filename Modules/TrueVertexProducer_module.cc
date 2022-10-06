////////////////////////////////////////////////////////////////////////
// Class:       TrueVertexProducer
// Plugin Type: producer (art v3_01_02)
// File:        TrueVertexProducer_module.cc
//
// Generated at Wed Aug 21 17:07:38 2019 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardata/Utilities/AssociationUtil.h"

#include <memory>

class TrueVertexProducer;


class TrueVertexProducer : public art::EDProducer {
public:
  explicit TrueVertexProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrueVertexProducer(TrueVertexProducer const&) = delete;
  TrueVertexProducer(TrueVertexProducer&&) = delete;
  TrueVertexProducer& operator=(TrueVertexProducer const&) = delete;
  TrueVertexProducer& operator=(TrueVertexProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

};


TrueVertexProducer::TrueVertexProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  produces< std::vector<recob::Vertex> >();
  produces< art::Assns<recob::Slice,recob::Vertex,void> >();
}

void TrueVertexProducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  std::unique_ptr< std::vector<recob::Vertex> > vertexcol(new std::vector<recob::Vertex>);
  std::unique_ptr< art::Assns<recob::Slice, recob::Vertex> > slicevtxassn(new art::Assns<recob::Slice, recob::Vertex>);

  //auto const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  //auto const* detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>(); 

  // generator, simb::MCTruth
  art::InputTag GenInputTag("generator");
  art::ValidHandle<std::vector<simb::MCTruth> > inputMCTruth = e.getValidHandle<std::vector<simb::MCTruth> >(GenInputTag);

  if (inputMCTruth->size()<1) return;

  const auto& mct = inputMCTruth->at(0);
  /*
     double nux = mct.GetNeutrino().Lepton().Vx();
     double nuy = mct.GetNeutrino().Lepton().Vy();
     double nuz = mct.GetNeutrino().Lepton().Vz();
     */
 
  if(mct.NParticles() < 1) return;
  simb::MCParticle Part = mct.GetParticle(0); 
  double nux = Part.Position().X();
  double nuy = Part.Position().Y();
  double nuz = Part.Position().Z();
     
  auto scecorr = SCE->GetPosOffsets({nux,nuy,nuz});
  //double g4Ticks = detClocks->TPCG4Time2Tick(mct.GetNeutrino().Lepton().T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
  double xOffset = /*theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-*/-scecorr.X()+0.6;
  double yOffset = scecorr.Y();
  double zOffset = scecorr.Z();
  recob::tracking::Point_t mcpos_unc(nux,nuy,nuz);
  recob::tracking::Point_t mcpos(nux+xOffset,nuy+yOffset,nuz+zOffset);

  std::cout << "True PV = " << nux << " " << nuy << " " << nuz << std::endl;
  std::cout << "Corrections = " << xOffset << " " << yOffset << " " << zOffset << std::endl;
/*
  std::cout << "mct.GetNeutrino().Lepton().T()=" << mct.GetNeutrino().Lepton().T() 
	    << " detClocks->TPCG4Time2Tick(mct.GetNeutrino().Lepton().T())=" << detClocks->TPCG4Time2Tick(mct.GetNeutrino().Lepton().T()) 
	    << " detClocks->TPCTime()=" << detClocks->TPCTime()
	    << " detClocks->TPCClock().TickPeriod()=" << detClocks->TPCClock().TickPeriod()
	    << " detClocks->TriggerTime()=" << detClocks->TriggerTime()
	    << " detClocks->TriggerOffsetTPC()=" << detClocks->TriggerOffsetTPC()
	    << " theDetector->GetXTicksOffset(0,0,0)=" << theDetector->GetXTicksOffset(0,0,0) 
	    << " theDetector->TriggerOffset()=" << theDetector->TriggerOffset() 
	    << std::endl;
  std::cout << "g4Ticks=" << g4Ticks << " theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)=" << theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0) << " scecorr.X()=" << scecorr.X() << std::endl;
  std::cout << "mcpos_unc=" << mcpos_unc << " mcpos=" << mcpos << std::endl;
*/
  double xyz[3] = { mcpos.X(), mcpos.Y(), mcpos.Z() };
  recob::Vertex trueVtx(xyz, -1);

  art::InputTag PfInputTag("pandora");
  art::ValidHandle<std::vector<recob::PFParticle> > inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle> >(PfInputTag);
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(inputPfParticle, e, PfInputTag);
  auto assocSlice = std::unique_ptr<art::FindManyP<recob::Slice> >(new art::FindManyP<recob::Slice>(inputPfParticle, e, PfInputTag));


  // first case: vertexcol->push_back(trueVtx);
  art::InputTag PatRecInputTag("pandora");
  art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(PatRecInputTag);
  art::InputTag CandVtxInputTag("pandoraPatRec","allcandidatevertices");
  art::FindManyP< recob::Vertex > sliceToVertexAssoc(inputSlice, e, CandVtxInputTag);
  // second case: need to get the collection of slices produced at the same time as the candidate vertices, then we map through the slice.ID() value
  // art::InputTag PatRecInputTag("pandoraPatRec");
  // art::ValidHandle<std::vector<recob::Slice> > inputSlice = e.getValidHandle<std::vector<recob::Slice> >(PatRecInputTag);
  // art::InputTag CandVtxInputTag("pandoraPatRec","allcandidatevertices");
  // art::FindManyP< recob::Vertex > sliceToVertexAssoc(inputSlice, e, CandVtxInputTag);

  for (unsigned int inpf=0; inpf<inputPfParticle->size(); ++inpf) {
    art::Ptr<recob::PFParticle> npfp(inputPfParticle,inpf);
    bool isTheNeutrino = false;
    const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(npfp.key()));
    if (!pfParticleMetadataList.empty()) {
      for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j) {
	const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
	const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
	for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
	  if (it->first=="IsNeutrino" && it->second==1) isTheNeutrino = true;
	}
      }
    }
    if (npfp->IsPrimary()==false || isTheNeutrino==false) continue; 
    auto slices = assocSlice->at(npfp.key());
    if (slices.size()!=1) {
      std::cout << "WRONG!!! n slices = " << slices.size() << std::endl;
    }

    if (1) {
      vertexcol->push_back(trueVtx);
    } else {
      for (size_t isl=0;isl<inputSlice->size();isl++) {
	if (inputSlice->at(isl).ID()!=slices[0]->ID()) continue;
	float minDist = 9999999.;
	int closestIdx = -1;
	auto candVertices = sliceToVertexAssoc.at(isl);
	for (size_t iv = 0; iv < candVertices.size(); iv++) {
	  float dist = (trueVtx.position() - candVertices.at(iv)->position()).R();
	  std::cout << "consider vertex at iv=" << iv << " pos=" << candVertices.at(iv)->position() << " dist=" << dist << std::endl;
	  if (dist<minDist) {
	    minDist = dist;
	    closestIdx = iv;
	  }
	}
	if (closestIdx>=0) {
	  std::cout << "closestIdx=" << closestIdx << " minDist=" << minDist << std::endl;
	  recob::Vertex closestVtx(*candVertices.at(closestIdx));
	  vertexcol->push_back(closestVtx);
	}
      }
    }
    if (vertexcol->size()!=0) {
      //the association is mapped wrt the "pandora" slice collection
      std::cout << "create assn slice ID=" << slices[0]->ID() << " id=" << slices[0].id().value()  << " key=" << slices[0].key() << std::endl;
      util::CreateAssn(*this, e, *vertexcol, slices[0], *slicevtxassn);
      // auto slicehit = assocSliceHit->at(slices[0].key());
    }
  }

  e.put(std::move(vertexcol));
  e.put(std::move(slicevtxassn));
  return;
}

DEFINE_ART_MODULE(TrueVertexProducer)
