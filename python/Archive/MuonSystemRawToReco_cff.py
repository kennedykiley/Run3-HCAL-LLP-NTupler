import FWCore.ParameterSet.Config as cms

import EventFilter.CSCRawToDigi.cscUnpacker_cfi
muonCSCDigis = EventFilter.CSCRawToDigi.cscUnpacker_cfi.muonCSCDigis.clone()
muonCSCDigis.InputObjects = 'rawDataCollector'

import EventFilter.DTRawToDigi.dtunpacker_cfi
muonDTDigis = EventFilter.DTRawToDigi.dtunpacker_cfi.muonDTDigis.clone()
muonDTDigis.inputLabel = 'rawDataCollector'

import EventFilter.RPCRawToDigi.rpcUnpacker_cfi
muonRPCDigis = EventFilter.RPCRawToDigi.rpcUnpacker_cfi.rpcunpacker.clone()
muonRPCDigis.InputLabel = 'rawDataCollector'

from RecoLocalMuon.CSCRecHitD.cscRecHitD_cfi import *
from RecoLocalMuon.DTRecHit.dt1DRecHits_LinearDriftFromDB_cfi import *
from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import *

muonSystemClusterSelSeq = cms.Sequence(
			   muonCSCDigis * muonDTDigis * muonRPCDigis * dt1DRecHits * csc2DRecHits * rpcRecHits
                           )

