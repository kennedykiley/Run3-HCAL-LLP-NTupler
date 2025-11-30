import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# MET Filter Recommendations: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Run_3_2022_and_2023_data_and_MC
# Also see: https://github.com/cms-sw/cmssw/blob/a31424ec6d9989e1390981466d09b137ad068318/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py#L8
"""
from RecoMET.METFilters.metFilters_cff import primaryVertexFilter
from RecoMET.METFilters.metFilters_cff import globalSuperTightHalo2016Filter
from RecoMET.METFilters.metFilters_cff import EcalDeadCellTriggerPrimitiveFilter
from RecoMET.METFilters.metFilters_cff import BadPFMuonFilter
from RecoMET.METFilters.metFilters_cff import BadPFMuonDzFilter
from RecoMET.METFilters.metFilters_cff import BadChargedCandidateFilter
from RecoMET.METFilters.metFilters_cff import hfNoisyHitsFilter
from RecoMET.METFilters.metFilters_cff import eeBadScFilter
from RecoMET.METFilters.metFilters_cff import ecalBadCalibFilter
"""

from RecoMET.METFilters.metFilters_cff import HBHENoiseFilterResultProducer, HBHENoiseFilter, HBHENoiseIsoFilter, hcalLaserEventFilter
from RecoMET.METFilters.metFilters_cff import EcalDeadCellTriggerPrimitiveFilter, eeBadScFilter, ecalLaserCorrFilter, EcalDeadCellBoundaryEnergyFilter, ecalBadCalibFilter
from RecoMET.METFilters.metFilters_cff import primaryVertexFilter, CSCTightHaloFilter, CSCTightHaloTrkMuUnvetoFilter, CSCTightHalo2015Filter, globalTightHalo2016Filter, globalSuperTightHalo2016Filter, HcalStripHaloFilter
from RecoMET.METFilters.metFilters_cff import goodVertices, trackingFailureFilter, trkPOGFilters, manystripclus53X, toomanystripclus53X, logErrorTooManyClusters
from RecoMET.METFilters.metFilters_cff import chargedHadronTrackResolutionFilter, muonBadTrackFilter
from RecoMET.METFilters.metFilters_cff import BadChargedCandidateFilter, BadPFMuonFilter, BadPFMuonDzFilter #2016 post-ICHEPversion
from RecoMET.METFilters.metFilters_cff import BadChargedCandidateSummer16Filter, BadPFMuonSummer16Filter #2016 ICHEP version
from RecoMET.METFilters.metFilters_cff import hfNoisyHitsFilter

# Recommended

primaryVertexFilter.src = cms.InputTag("offlineSlimmedPrimaryVertices")

# Tagging mode (true = events pass with boolean Flag_; false = events fail)
primaryVertexFilter.taggingMode = cms.bool(False)
globalSuperTightHalo2016Filter.taggingMode = cms.bool(True)
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)
BadPFMuonFilter.taggingMode = cms.bool(True)
BadPFMuonDzFilter.taggingMode = cms.bool(True)
BadChargedCandidateFilter.taggingMode = cms.bool(True)
hfNoisyHitsFilter.taggingMode = cms.bool(True)
eeBadScFilter.taggingMode = cms.bool(True)
ecalBadCalibFilter.taggingMode = cms.bool(False)

metFiltersRecommended = cms.Sequence(
primaryVertexFilter *
globalSuperTightHalo2016Filter *
EcalDeadCellTriggerPrimitiveFilter *
BadPFMuonFilter *
BadPFMuonDzFilter *
BadChargedCandidateFilter *
hfNoisyHitsFilter *
eeBadScFilter *
ecalBadCalibFilter
)

metFiltersRecommended_MINIAOD = cms.Sequence(
globalSuperTightHalo2016Filter *
EcalDeadCellTriggerPrimitiveFilter
#BadPFMuonFilter *
#BadPFMuonDzFilter *
#BadChargedCandidateFilter *
#hfNoisyHitsFilter *
#eeBadScFilter *
#ecalBadCalibFilter
)

# All

metFiltersAll = cms.Sequence(
HBHENoiseFilterResultProducer *
HBHENoiseFilter *
HBHENoiseIsoFilter *
hcalLaserEventFilter *
EcalDeadCellTriggerPrimitiveFilter *
eeBadScFilter *
ecalLaserCorrFilter *
EcalDeadCellBoundaryEnergyFilter *
ecalBadCalibFilter *
primaryVertexFilter *
CSCTightHaloFilter *
CSCTightHaloTrkMuUnvetoFilter *
CSCTightHalo2015Filter *
globalTightHalo2016Filter *
globalSuperTightHalo2016Filter *
HcalStripHaloFilter *
goodVertices *
trackingFailureFilter *
trkPOGFilters *
manystripclus53X *
toomanystripclus53X *
logErrorTooManyClusters *
chargedHadronTrackResolutionFilter *
muonBadTrackFilter *
BadChargedCandidateFilter *
BadPFMuonFilter *
BadPFMuonDzFilter *
BadChargedCandidateSummer16Filter *
BadPFMuonSummer16Filter *
hfNoisyHitsFilter
)

