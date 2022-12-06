////////////////////////////////////////////////////////////////////////
// Class:       KaonVertexProducer
// Plugin Type: producer (art v3_01_02)
// File:        KaonVertexProducer_module.cc
//
// Generated at Wed Aug 21 17:07:38 2019 by Giuseppe Cerati using cetskelgen
// Modified by C Thorpe Sept 2022.
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
#include "lardataobj/RecoBase/Track.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "ubana/HyperonProduction/Modules/SubModules/SubModuleReco.h"

#include <memory>

namespace hyperon {
class KaonVertexProducer;
}

class hyperon::KaonVertexProducer : public art::EDProducer {
   public:
      explicit KaonVertexProducer(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      KaonVertexProducer(KaonVertexProducer const&) = delete;
      KaonVertexProducer(KaonVertexProducer&&) = delete;
      KaonVertexProducer& operator=(KaonVertexProducer const&) = delete;
      KaonVertexProducer& operator=(KaonVertexProducer&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

   private:

      // Declare member data here.
      std::string f_TrackLabel;
      std::string f_PFParticleLabel;
      bool f_KaonCheat = false;
      fhicl::ParameterSet f_Reco;
      size_t f_PassNo;

      BDTHandle KaonBDT;
};


hyperon::KaonVertexProducer::KaonVertexProducer(fhicl::ParameterSet const& p)
   : EDProducer{p},
   f_TrackLabel(p.get<std::string>("TrackLabel","pandora")),
   f_PFParticleLabel(p.get<std::string>("PFParticleLabel","pandora")),
   f_KaonCheat(p.get<bool>("KaonCheat",false)),
   f_Reco(p.get<fhicl::ParameterSet>("Reco")),
   f_PassNo(p.get<size_t>("PassNo",0)),
   KaonBDT(p.get<std::string>("BDTWeightDir"))
   // More initializers here.
{
   // Call appropriate produces<>() functions here.
   // Call appropriate consumes<>() for any products to be retrieved by this module.
   produces< std::vector<recob::Vertex> >();
   produces< art::Assns<recob::Slice,recob::Vertex,void> >();
}

void hyperon::KaonVertexProducer::produce(art::Event& e)
{
   std::unique_ptr< std::vector<recob::Vertex> > vertexcol(new std::vector<recob::Vertex>);
   std::unique_ptr< art::Assns<recob::Slice, recob::Vertex> > slicevtxassn(new art::Assns<recob::Slice, recob::Vertex>);

   int tracktouseindex = -1;

   SubModuleReco* Reco_SM = new hyperon::SubModuleReco(e,false,f_Reco);
   Reco_SM->PrepareInfo();
   RecoData RecoD = Reco_SM->GetInfo(); 
   if(f_KaonCheat){
      for(size_t i_tr=0;i_tr<RecoD.TrackPrimaryDaughters.size();i_tr++){
         if(RecoD.TrackPrimaryDaughters.at(i_tr).TrackTrueOrigin == 1 && RecoD.TrackPrimaryDaughters.at(i_tr).TrackTruePDG == 321){
            tracktouseindex = RecoD.TrackPrimaryDaughters.at(i_tr).Index;
            std::cout << "Found true kaon at index " << tracktouseindex << std::endl;
            break;
         }
      }
   }
   else {
      std::sort(RecoD.TrackPrimaryDaughters.begin(),RecoD.TrackPrimaryDaughters.end(),
            [&](RecoParticle A,RecoParticle B) -> bool {
            //return A.Track_LLR_PID_Kaon > B.Track_LLR_PID_Kaon;
            return KaonBDT.GetScore(A.Track_LLR_PID_Kaon,A.Track_Bragg_PID_Kaon) > KaonBDT.GetScore(B.Track_LLR_PID_Kaon,B.Track_Bragg_PID_Kaon); 
          });

      if(RecoD.TrackPrimaryDaughters.size() > f_PassNo) tracktouseindex = RecoD.TrackPrimaryDaughters.at(f_PassNo).Index; 
   }    

   delete Reco_SM; 
  
   std::cout << "Looking for track with index " << tracktouseindex << std::endl;

   // Load the reco'd tracks from the event
   art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
   std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;
   art::Handle<std::vector<recob::Track>> Handle_Track;
   std::vector<art::Ptr<recob::Track>> Vect_Track;

   if(!e.getByLabel(f_PFParticleLabel,Handle_PFParticle)) 
      throw cet::exception("KaonVertexProducer") << "No PFParticle Data Products Found!" << std::endl;

   if(!e.getByLabel(f_TrackLabel,Handle_Track)) 
      throw cet::exception("KaonVertexProducer") << "No Track Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);
   art::fill_ptr_vector(Vect_Track,Handle_Track);
   art::FindManyP<recob::Track> Assoc_PFParticleTrack = art::FindManyP<recob::Track>(Vect_PFParticle,e,f_TrackLabel);

   // Get the ID of the neutrino candidate
   size_t neutrinoid = -9999;
   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle)
      if(pfp->IsPrimary() && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14))
         neutrinoid = pfp->Self();

   std::cout << "Found neutrino candidate with ID " << neutrinoid << std::endl;  

   double xyz[3] = {-1000,-1000,-1000};
   int index = -1;
   std::vector<size_t> PFP_IDs;

   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
      if(pfp->Parent() != neutrinoid && std::find(PFP_IDs.begin(),PFP_IDs.end(),pfp->Parent()) == PFP_IDs.end()) continue;
      std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack.at(pfp.key());
      PFP_IDs.push_back(pfp->Self());
      if(pfpTracks.size() != 1) continue;
      index++;
      if(index == tracktouseindex){
         std::cout << "Using track with index " << index << std::endl;
         art::Ptr<recob::Track> trk = pfpTracks.at(0);
         xyz[0] = trk->End().X();
         xyz[1] = trk->End().Y();
         xyz[2] = trk->End().Z();
         break;
      }
   }

   if(tracktouseindex != -1 && xyz[0] == -1000 && xyz[1] == -1000 && xyz[2] == -1000)
      throw cet::exception("KaonVertexProducer") << "Couldn't find a track with the index requested" << std::endl;
   
   recob::Vertex newVtx(xyz, -1);

   art::InputTag PfInputTag("pandora");
   art::ValidHandle<std::vector<recob::PFParticle> > inputPfParticle = e.getValidHandle<std::vector<recob::PFParticle> >(PfInputTag);
   art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(inputPfParticle, e, PfInputTag);
   auto assocSlice = std::unique_ptr<art::FindManyP<recob::Slice> >(new art::FindManyP<recob::Slice>(inputPfParticle, e, PfInputTag));

   // first case: vertexcol->push_back(newVtx);
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
         vertexcol->push_back(newVtx);
      } else {
         for (size_t isl=0;isl<inputSlice->size();isl++) {
            if (inputSlice->at(isl).ID()!=slices[0]->ID()) continue;
            float minDist = 9999999.;
            int closestIdx = -1;
            auto candVertices = sliceToVertexAssoc.at(isl);
            for (size_t iv = 0; iv < candVertices.size(); iv++) {
               float dist = (newVtx.position() - candVertices.at(iv)->position()).R();
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

DEFINE_ART_MODULE(hyperon::KaonVertexProducer)
