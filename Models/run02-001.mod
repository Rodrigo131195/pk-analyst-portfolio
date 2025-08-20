;; 1. Based on:
$PROBLEM 
;; 2. Description: Base
;; x1. Author: munoz01r
;; 3. Label: Basic model
;----------------------------------
;; 1 CMT infusion
;----------------------------------
$INPUT ID TIME TAD DV MDV CMT EVID AMT SS II VISIT TSFD WT ADA MTX
;----------------------------------
$DATA fnm_adl_base.csv IGNORE=@ ; The ignore is crucial
;----------------------------------
$SUBROUTINES ADVAN6 TRANS1 TOL=6 ;  TOL=9
;----------------------------------
$MODEL
  COMP=(DEPOT) COMP=(CENTRAL) COMP=(PERI)

$PK

CL = THETA(1) * EXP(ETA(1)) ; * (WT/70)^0.75 
V1  = THETA(2) * EXP(ETA(2))
V2  = THETA(3) * EXP(ETA(3))
Q = THETA(4) * EXP(ETA(4))
;S1 = V1
K20 = CL / V1
K23 = Q / V1
K32 = Q / V2
KA = THETA(5) * EXP(ETA(5))
F1 = THETA(6) 
S2 = V1
CMTX = CMT 
REP = IREP

$DES
  CONC    = A(2) / V1
  DADT(1) = -KA*A(1)
  DADT(2) = KA*A(1)  - K20*A(2) + K32*A(3) - K23*A(2) 
  DADT(3) = K23*A(2) - K32*A(3)

$ERROR
IPRED = F
    W = SQRT(THETA(7)**2*IPRED**2 + THETA(8)**2)   ; weighted error model
    Y = IPRED + W*EPS(1)
 IRES = DV-IPRED
IWRES = IRES/W

$THETA 
0.181 FIX ; T_CL (L/d)
3.36 FIX ; T_VC (L)
2.56 FIX ; T_VP (L)
0.57 FIX ; T_Q
0.25 FIX ; T_ka (/d) 
0.67 FIX ; T_F (fraction) 
(0, 0.2) ; T_Prop_RE (sd) 
(0, 1) ; T_Add_RE 

$OMEGA
0.17 FIX ; IIV_CL
0.075 FIX ; IIV_VC
0.14 FIX ; IIV_VP
0.37 FIX ; IIV_Q
0.2 FIX ; IIV_ka

$SIGMA
1 FIX ; Residual error

; $EST METHOD=1 INTER MAXEVAL=9999 NOABORT SIGL=9 NSIG=3 PRINT=1 POSTHOC
$EST METHOD=1 INTER MAXEVAL=9999 NOABORT SIGL=6 NSIG=2 PRINT=1 POSTHOC
$COV PRINT=E UNCONDITIONAL SIGL=6

$TABLE ID TIME TAD DV MDV EVID PRED IPRED IWRES CWRES CMTX VISIT WT ADA MTX ONEHEADER NOPRINT FILE=sdtab02-001
$TABLE ID CL V1 V2 Q KA ETA1 ETA2 ETA3 ETA4 ETA5 ONEHEADER NOPRINT FILE=patab02-001
