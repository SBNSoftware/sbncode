#include "sbncode/CAFMaker/FillQELambda.h"
  
void caf::FillFRANSLambda( const qelambda::QELambdaReco LAr_QELambdaRecoProducts, caf::SRQELambda& CAF_QELambdaReco)
{

  CAF_QELambdaReco.SliceId = LAr_QELambdaRecoProducts.SliceId;

  CAF_QELambdaReco.FransNUnassocOrigins       = LAr_QELambdaRecoProducts.FransNUnassocOrigins;
  CAF_QELambdaReco.FransNMaxDirtUnassocHits   = LAr_QELambdaRecoProducts.FransNMaxDirtUnassocHits;
  CAF_QELambdaReco.FransNAngles               = LAr_QELambdaRecoProducts.FransNAngles;
  CAF_QELambdaReco.FransContainedPhi          = LAr_QELambdaRecoProducts.FransContainedPhi;
  CAF_QELambdaReco.FransPzsign                = LAr_QELambdaRecoProducts.FransPzsign;

  CAF_QELambdaReco.FransScore                 = LAr_QELambdaRecoProducts.FransScore;
  CAF_QELambdaReco.FransIota                  = LAr_QELambdaRecoProducts.FransIota;
  CAF_QELambdaReco.FransDelta                 = LAr_QELambdaRecoProducts.FransDelta;
  CAF_QELambdaReco.FransEta                   = LAr_QELambdaRecoProducts.FransEta;
  CAF_QELambdaReco.FransFitScore              = LAr_QELambdaRecoProducts.FransFitScore;

  CAF_QELambdaReco.FransTrack1ID = LAr_QELambdaRecoProducts.FransTrack1ID;
  CAF_QELambdaReco.FransTrack2ID = LAr_QELambdaRecoProducts.FransTrack2ID;
  CAF_QELambdaReco.FransMainTrackID = LAr_QELambdaRecoProducts.FransMainTrackID;

  CAF_QELambdaReco.HyperonVtxX      = LAr_QELambdaRecoProducts.HyperonVtxX;
  CAF_QELambdaReco.HyperonVtxY      = LAr_QELambdaRecoProducts.HyperonVtxY;
  CAF_QELambdaReco.HyperonVtxZ      = LAr_QELambdaRecoProducts.HyperonVtxZ;
  CAF_QELambdaReco.PrimLeptonVtxX   = LAr_QELambdaRecoProducts.PrimLeptonVtxX;
  CAF_QELambdaReco.PrimLeptonVtxY   = LAr_QELambdaRecoProducts.PrimLeptonVtxY;
  CAF_QELambdaReco.PrimLeptonVtxZ   = LAr_QELambdaRecoProducts.PrimLeptonVtxZ;

  CAF_QELambdaReco.MichelNMuonTracks  = LAr_QELambdaRecoProducts.MichelNMuonTracks;
  CAF_QELambdaReco.MichelRssRatio     = LAr_QELambdaRecoProducts.MichelRssRatio;
  CAF_QELambdaReco.MichelAngleQ80     = LAr_QELambdaRecoProducts.MichelAngleQ80;
  CAF_QELambdaReco.MichelDropRadius   = LAr_QELambdaRecoProducts.MichelDropRadius;
  CAF_QELambdaReco.MichelNBlips       = LAr_QELambdaRecoProducts.MichelNBlips;
  CAF_QELambdaReco.MichelScore        = LAr_QELambdaRecoProducts.MichelScore;

}

void caf::FillLambdaPID( const qelambda::QELambdaPID LAr_QELambdaPID, caf::SRQELambda& CAF_QELambdaReco)
{
  CAF_QELambdaReco.SecondPandoraVtxX            = LAr_QELambdaPID.VtxX;
  CAF_QELambdaReco.SecondPandoraVtxY            = LAr_QELambdaPID.VtxY;
  CAF_QELambdaReco.SecondPandoraVtxZ            = LAr_QELambdaPID.VtxZ;
  CAF_QELambdaReco.SecondPandoraProtonKE        = LAr_QELambdaPID.ProtonKE;
  CAF_QELambdaReco.SecondPandoraPionKE          = LAr_QELambdaPID.PionKE;
  CAF_QELambdaReco.SecondPandoraProtonPX        = LAr_QELambdaPID.ProtonPX;
  CAF_QELambdaReco.SecondPandoraProtonPY        = LAr_QELambdaPID.ProtonPY;
  CAF_QELambdaReco.SecondPandoraProtonPZ        = LAr_QELambdaPID.ProtonPZ;
  CAF_QELambdaReco.SecondPandoraPionPX          = LAr_QELambdaPID.PionPX;
  CAF_QELambdaReco.SecondPandoraPionPY          = LAr_QELambdaPID.PionPY;
  CAF_QELambdaReco.SecondPandoraPionPZ          = LAr_QELambdaPID.PionPZ;
}
