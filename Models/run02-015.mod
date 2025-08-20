;; 1. Based on: run02-012
$PROBLEM 
;; 2. Description: MM + groups gow
;; x1. Author: munoz01r
;; 3. Label: Basic model + Michaelis Menten and 2 new compartments
;----------------------------------
;; 1 CMT infusion
;----------------------------------
$INPUT ID TIME TAD DV MDV CMT EVID AMT SS II VISIT TSFD WT ADA MTX ADA_G
$DATA fnm_adl_10.csv IGNORE=@
$SUBROUTINES ADVAN6 TRANS1 TOL=6

$MODEL
  COMP=(DEPOT) COMP=(CENTRAL) COMP=(PERI) COMP=(COMPLX) COMP=(ADA)

$PK

CL = THETA(1) * EXP(ETA(1))
CLX = THETA(2) * EXP(ETA(7))
CLA = CL ; THETA(3) 

V2 = THETA(3) * EXP(ETA(2))
V3 = THETA(4) * EXP(ETA(3))
V4 = V2
V5 = V2

KA = THETA(5) * EXP(ETA(4))
 ;KADA = THETA(6) * EXP(ETA(5))

Q23 = THETA(7) * EXP(ETA(6))

X = THETA(8)
VMAX = THETA(12)
KM   = THETA(13)

K20 = CL / V2
K23 = Q23 / V2
K32 = Q23 / V3
K40 = CLX / V4
K50 = CLA / V5

F1 = THETA(9)
S2 = V2
S4 = V4
CMTX = CMT
REP = IREP

; ADA group-specific maxADA (2 groups only)
IF (ADA_G.EQ.0) KADA = THETA(6)  * EXP(ETA(5))
IF (ADA_G.EQ.1) KADA = THETA(14)  * EXP(ETA(5))
IF (ADA_G.EQ.2) KADA = THETA(15)  * EXP(ETA(5))

$DES
  CONC_DRUG = A(2) / V2
  CONC_ADA  = A(5) / V5
  MM_BIND = VMAX * CONC_DRUG * CONC_ADA / (KM*KM + KM*(CONC_DRUG + CONC_ADA) + CONC_DRUG*CONC_ADA)

  DADT(1) = -KA * A(1)
  DADT(2) = KA * A(1) + K32 * A(3) - K23 * A(2) - K20 * A(2) - MM_BIND
  DADT(3) = K23 * A(2) - K32 * A(3)
  DADT(4) = - K40 * A(4) + MM_BIND
  DADT(5) = KADA - K50 * A(5) - X * MM_BIND

$ERROR
  IPRED = F
  W = SQRT(THETA(10)**2 * IPRED**2 + THETA(11)**2)
  Y = IPRED + W * EPS(1)
  IRES = DV - IPRED
  IWRES = IRES / W

$THETA 
  0.181 FIX ; T_CL (L/d)
  (0, 0.3)     ; T_CLX (L/d)
  ; 0.22     ; T_CLA (L/d)
  3.36 FIX ; T_VC (L)
  2.56 FIX ; T_VP (L)
  0.25 FIX ; T_KA (/d)
  (0, 0.25)  ; T_KADA0 (/d)
  0.57 FIX ; T_Qcp
  (0.01, 0.1)   ; T_X scaling factor
  0.67 FIX ; T_F (fraction) 
  (0, 0.2) ; Prop error
  (0, 1)   ; Add error
  (0.01, 0.5)      ; T_VMAX (mg/L/day)
  (0.01, 1.0)    ; T_KM (mg/L)
  (0.01, 0.25)  ; T_KADA1 (/d)
  (0.05, 0.5)  ; T_KADA2 (/d)

$OMEGA
  0.17  FIX ; IIV_CL
  0.075 FIX ; IIV_VC
  0.14  FIX ; IIV_VP
  0.2   FIX ; IIV_KA
  0.2    ; IIV_KADA
  0.37   FIX ; IIV_QCP
  0.3 ; IIV_CLX

$SIGMA
  1 FIX ; Residual error

$EST METHOD=1 INTER MAXEVAL=9999 NOABORT SIGL=6 NSIG=2 PRINT=1 POSTHOC
$COV PRINT=E UNCONDITIONAL SIGL=6

$TABLE ID TIME TAD DV MDV EVID PRED IPRED IWRES CWRES CMTX VISIT WT ADA MTX ADA_G ONEHEADER NOPRINT FILE=sdtab02-015
$TABLE ID CL CLX KM VMAX X V2 V3 V4 V5 Q23 KA KADA ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ONEHEADER NOPRINT FILE=patab02-015
