import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from RecoMET.METFilters.metFilters_cff import HBHENoiseFilterResultProducer
from RecoMET.METFilters.metFilters_cff import HBHENoiseFilter
from RecoMET.METFilters.metFilters_cff import primaryVertexFilter
from RecoMET.METFilters.metFilters_cff import globalSuperTightHalo2016Filter
from RecoMET.METFilters.metFilters_cff import globalTightHalo2016Filter
from RecoMET.METFilters.metFilters_cff import EcalDeadCellTriggerPrimitiveFilter
from RecoMET.METFilters.metFilters_cff import BadPFMuonFilter
from RecoMET.METFilters.metFilters_cff import BadChargedCandidateFilter
from RecoMET.METFilters.metFilters_cff import eeBadScFilter
# from RecoMET.METFilters.metFilters_cff import ecalBadCalibFilter
# HBHENoiseFilter.taggingMode = cms.bool(True)
#primaryVertexFilter.taggingMode = cms.bool(True)
globalSuperTightHalo2016Filter.taggingMode = cms.bool(True)
globalTightHalo2016Filter.taggingMode = cms.bool(True)
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)
BadPFMuonFilter.taggingMode = cms.bool(True)
BadChargedCandidateFilter.taggingMode = cms.bool(True)
eeBadScFilter.taggingMode = cms.bool(True)
# ecalBadCalibFilter.taggingMode = cms.bool(True)
baddetEcallist2017 = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649])
baddetEcallist2018 = cms.vuint32(
    [872422436,872421950,872437185,872422564,
     872421566,872421695,872421955,872421567,
     872437184,872421951,872421694])
ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
"EcalBadCalibFilter",
EcalRecHitSource = cms.InputTag("reducedEcalRecHitsEE"),
ecalMinEt        = cms.double(50.),
baddetEcal    = baddetEcallist2017, #use baddetEcallist2018  for 2018 analysis
taggingMode = cms.bool(True),
debug = cms.bool(False)
)
metFilters = cms.Sequence(
HBHENoiseFilterResultProducer *
# HBHENoiseFilter *
primaryVertexFilter *
globalSuperTightHalo2016Filter *
globalTightHalo2016Filter *
EcalDeadCellTriggerPrimitiveFilter *
BadPFMuonFilter *
BadChargedCandidateFilter *
eeBadScFilter *
ecalBadCalibReducedMINIAODFilter
)

