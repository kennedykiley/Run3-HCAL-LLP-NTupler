// ======================================================================================
// TriggerSFTables.h
//
// Measured data/MC tagging efficiency tables (eData, eQCD) for each HLT filter leg,
// used by DisplacedHcalJetAnalyzer::GetEfficiencies() to derive per-jet scale factors.
//
// Each table is a set of three parallel arrays: bin edges, data efficiency, MC
// efficiency, plus a bin count. Lookup convention (see GetEfficiencies): the last
// bin acts as a catch-all for any value at or above its edge.
//
// ======================================================================================

#ifndef TRIGGERSFTABLES_H
#define TRIGGERSFTABLES_H

// DisplacedTrack_DTrk: 13 bins, [0, 12]
extern const double Bins_DisplacedTrack_DTrk[13];
extern const double Data_DisplacedTrack_DTrk[13];
extern const double MC_DisplacedTrack_DTrk[13];
extern const int N_DisplacedTrack_DTrk;

// DisplacedTrack_pT: 23 bins, [30, 118]
extern const double Bins_DisplacedTrack_pT[23];
extern const double Data_DisplacedTrack_pT[23];
extern const double MC_DisplacedTrack_pT[23];
extern const int N_DisplacedTrack_pT;

// DisplacedTrack_PTrk: 20 bins, [0, 19]
extern const double Bins_DisplacedTrack_PTrk[20];
extern const double Data_DisplacedTrack_PTrk[20];
extern const double MC_DisplacedTrack_PTrk[20];
extern const int N_DisplacedTrack_PTrk;

// Inclusive_pT: 23 bins, [30, 118]
extern const double Bins_Inclusive_pT[23];
extern const double Data_Inclusive_pT[23];
extern const double MC_Inclusive_pT[23];
extern const int N_Inclusive_pT;

// Inclusive_PTrk: 20 bins, [0, 19]
extern const double Bins_Inclusive_PTrk[20];
extern const double Data_Inclusive_PTrk[20];
extern const double MC_Inclusive_PTrk[20];
extern const int N_Inclusive_PTrk;

// PtrkShortSig5_pT_35GeV: 23 bins, [30, 118]
extern const double Bins_PtrkShortSig5_pT_35GeV[23];
extern const double Data_PtrkShortSig5_pT_35GeV[23];
extern const double MC_PtrkShortSig5_pT_35GeV[23];
extern const int N_PtrkShortSig5_pT_35GeV;

// PtrkShortSig5_pT_40GeV: 23 bins, [30, 118]
extern const double Bins_PtrkShortSig5_pT_40GeV[23];
extern const double Data_PtrkShortSig5_pT_40GeV[23];
extern const double MC_PtrkShortSig5_pT_40GeV[23];
extern const int N_PtrkShortSig5_pT_40GeV;

// PtrkShortSig5_PTrk: 20 bins, [0, 19]
extern const double Bins_PtrkShortSig5_PTrk[20];
extern const double Data_PtrkShortSig5_PTrk[20];
extern const double MC_PtrkShortSig5_PTrk[20];
extern const int N_PtrkShortSig5_PTrk;

// ---- HT tables (L1 and HLT-family), from earlier CSV uploads ----

// HT_L1: 11 bins, [75, 1252.5] -- SF-only table (no separate eData/eQCD; see note below)
extern const double HT_L1[11];
extern const double SF_L1[11];
extern const int N_L1;

// HT_Inclusive1PtrkShortSig5: SF-only table, measured at HT ~200-equivalent turn-on
extern const double HT_Inclusive1PtrkShortSig5[24];
extern const double SF_Inclusive1PtrkShortSig5[24];
extern const int N_Inclusive1PtrkShortSig5;

// HT_Inclusive: SF-only table, measured at HT320 (Inclusive family)
extern const double HT_Inclusive[30];
extern const double SF_Inclusive[30];
extern const int N_Inclusive;

#endif // TRIGGERSFTABLES_H
