Lstrandstart=10969; # about 2/3 through the genome
#(* the fractions of A,C,T, and G on the heavy and light strands of mtDNA *)
fTH=0.309; fTL=0.247;
fCH=0.131; fCL=0.313;
fGH=0.313; fGL=0.131;
fAH=0.247; fAL=0.309;
#(* the Hill coefficient of TK2 for thymidine *) 
tk2hill=0.5;
#(* The length of both strands of mtDNA *)
DNAlength=33136;
#(*  the length of one strand of mtDNA *)
strandDNA=DNAlength/2;

#(* volume of a mitochondrion *)
volmito=2e-16;  # i.e. .2 fL
#(* conversion factor used to convert K's and concentrations from microMolar to molecules/mitochondria *) 
# (* conversion = 120.4; *)

conversion=1e-6 * 6.022e23 * volmito;               
secondsperminute=60;

#      (* factor used to decrease the vmax of the polymerase on double stranded templates with lower primer density *)
dsfact=1/2;

#      (*  Polymerase kinetic constants KA Johnson JBC 2001 *)
VmPolT= 25.0 *dsfact* secondsperminute;
VmPolC= 43.0 *dsfact* secondsperminute;
VmPolA= 45.0 *dsfact* secondsperminute; 
VmPolG= 37.0 *dsfact* secondsperminute;
KpolT= 0.63 * conversion;  
KpolC= 0.9 * conversion;    
KpolA= 0.8 * conversion;    
KpolG= 0.8 * conversion; 

#      (* 2.3 Ki (inhibition constant) (microMolar) of dTTP on thymidine kinase 2  Eriksson JBC 2003 278 6963 *)
KiTtptk2 = 2.3*conversion;
#      (* 40 Ki of dCytidine on thymidine kinase 2 Eriksson JBC 2003 278 6963 *)
KiCtk2 = 40*conversion;
#      (* 0.83 Ki of dCTP on thymidine kinase 2 Eriksson JBC 2003 278 6963 *)
KiCtptk2 = 0.83*conversion;
#      (* 4.9 Ki of dThymidine on thymidine kinase 2 Eriksson JBC 2003 278 6963 *)
KiTtk2 = 4.9*conversion;
#      (* Ki of dGuanosine on deoxyguanosine kinase Eriksson FEBS Lett 2003 554 319 *)
KiGdgk = 4*conversion;
#      (* Ki of dAdenosine on deoxyguanosine kinase Eriksson FEBS Lett 2003 554 319 *)
KiAdgk = 467*conversion;
#      (* Ki of dGMP on deoxyguanosine kinase Eriksson Antimicrob Agents Chem  45 739 2001 *)
KiGmpdgk=4*conversion;
#      (* 2.1 Ki of dGDP on deoxyguanosine kinase Yamada Y BBA 1982 *)
KiGdpdgk=2.1*conversion;
#      (* 28 Ki of dAMP on deoxyguanosine kinase Eriksson Antimicrob Agents Chem  45 739 2001 *)
KiAmpdgk=28*conversion;
#      (* 41 Ki of dATP on deoxyguanosine kinase Eriksson Antimicrob Agents Chem  45 739 2001 *)
KiAtpdgk = 41*conversion;
#      (* 0.4 Ki of dGTP on deoxyguanosine kinase Eriksson Antimicrob Agents Chem  45 739 2001 *)
KiGtpdgk = 0.4*conversion;        
#      (* Ki of AMP on general 5prime nucleotidase *)
KiAMPtidase=94*conversion;
#(* concentration of (non-deoxy) NTP in the mitochondrial matrix in microMolar 8 mM Traut 94, 
#4 times higher than ADP Kunz W BBA 593 196 1980, ADP higher Mathews CK PNAS 100 12159 2003   *)
rATPin =  8000   *conversion; #(* 4000 *) TR guessing these were commented alternatives
#(* concentration of (non-deoxy) NDP in the mitochondrial matrix in microMolar *)
rADPin =  6000  *conversion; #(* 3500 *) 
#(* concentration of (non-deoxy) NMP in the mitochondrial matrix in microMolar *)
rAMPin =  2900  *conversion; # (* 1500 *)


rCMPin = rCDPin= 50 * conversion;  #(* mitochondrial concentration of CMP and CDP *)
KiCMP= KiCDP= 23 * conversion;     #(* Ki of CMP and CDP on CMP kinase *)
UMPin=100 * conversion;            #(* mitochondrial concentration of UMP *)
UDPin=200 *conversion;             #(* mitochondrial concentration of UDP *)
KiUMP=KiUDP=95*conversion;         #(* Ki of UMP and UDP on CMP kinase *)
rGMPin=200*conversion;             #(* mitochondrial concentration of GMP *)
rGDPin=400*conversion;             #(* mitochondrial concentration of GDP *)
KiGMP=KiGDP=17*conversion;         #(* Ki of GMP and GDP on GMP kinase *)


#(* Ki ATP for nucleotidase Greger J Enzyme 25 26 1980 -- 660 uM Hassinen IE BBA 1099 238 1992  *)
KiATPtidase=660*conversion;
KiADPtidase=1400*conversion;
#(* 75 Ki for matrix NDPK for ribonucleoside di and tri phosphates Lambeth DO JBC 1996 *)
KiADPNDPK=KiATPNDPK=75*conversion;


#(* estimated nucleoside transporter molecular weight (kD) Lai Y Unadkat JD JBC 2003  *)
nuctransMW=50;
#(* nucleoside kinase (dGK and TK2) molecular weight in kD Eriksson Febs Lett 1999 Eriksson Mol Pharm 1998 *)
nuckinMW = 28;
#(* molecular weight of mitochondrial nucleotidase (in kD) Bianchi PNAS 2000 97 8239 *)
tidaseMW = 25;
#(* 24 dTMP kinase molecular weight in Kd Fujimura S BBA 995 28 1989 *)
tmpkMW=24;
#(* nucleoside monophosphate kinase (NMPK) (used for C and G) molecular weight in Kd BBRC Karlsson A 2003 *)
nmpkMW = 25;
#(* adenylate kinase (nmpk enzyme for A) molecular weight in Kd Yoshinobu K 2001 Biochem J *)
akMW = 26;

# (* 66 pigeon nucleoside diphosphokinase (NDPK) molecular weight in Kd, 25 in humans 72 % homology 
#? Lambeth JBC 1997  and Milon L JBC 2000 *)
ndpkMW = 66;
# (* deoxynucleotide di/triphosphate exchanger molecular weight in kiloDaltons  Walker PNAS 98 2284 2000 *)
exMW = 32;

# (* molecules of enzyme in each mitochondrion from Saada A Nat Genet 29, 342 2001  
#TK2 Saada 2003 Mol Gen Metabolism 79 1 and dGK (27 in muscle 188 heart 300 fibroblast 500 in liver)  *)
tk2culespermito = 100;
dgkculespermito = 200;
# (* molecules of A,C,G,T nucleotidase in each mitochondrion *)
tidaseTpermito=50;
# (* not currently used *)
tidasepermito=50;
#              (* molecules of deoxythymidinemonophosphate kinase (dTMPK) in each mitochondrion *)
culestmpkpermito=50
#              (* molecules nmpk in each mitochondrion *)
nmpkculespermito=50;
#              (* molecules ndpk in each mitochondrion *)
ndpkculespermito=300;
#              (* number of deoxynucleotide di/triphosphate exchangers in a single mitochondrion Walker PNAS 2000 *)
expermito = 500;

#              (* the factor that the reverse reaction is faster than the forward reaction for nucleoside monophosphokinase (NMPK) *)
factorMD = 1/2 #(* rAMPin/ rADPin *) ;
#              (* the factor that the reverse reaction is faster than the forward reaction for nucleoside diphosphokinase (NDPK) *)
factorDT = rADPin/rATPin #(* 1/1.33=0.75 *) ;


#              (* 38 molecules dtransporter per mitochondrion Life Sciences Camins A 58, 753 1996 *)
culestrans = 38;
#              (* 450 molecules adenylate kinase per mitochondrion Eur. J. Biochem. 93, 263 1979 Tomaselli *) 
culesAK = 450;
#              (* number of total proteins in a mitochondrion assuming average MW of 30 kD and 5x10^10 mito/mg mito protein   *)
proteinspermito=400000;

#              (* concentration of dT in cplasm Traut 1994 0.8 Dudman npb Anal Biochem 115 428 1981 *)
dTc = 0.5*conversion;
#              (* concentration of dC in cplasm Traut 1994  0.7 Dudman npb Anal Biochem 115 428 1981 *)
dCc=0.5*conversion;
#              (* 0.5 concentration of dA in cplasm 1 pmol/1000000 cells J Amer Soc Nephrol 12, 1721 2001  *)
dAc=0.5*conversion;
#              (* 0.5 concentration of dG in cplasm  *)
dGc= 0.5*conversion;

#              (* dtransporter Km for Adenosine, T, and C  2 E Escubedo Eur J Pharm 398 31 2000   *)
dnuctransK = 2*conversion;
#              (* dtransporter Km for dguanosine Lewis RA Mol Cell Biochem 77 71 1987 *)
GdnuctransK = 0.64*conversion;
#  (* dtransporter Vm converted from micromoles substrate/mg mito/min 
#   to molecule substrate/molecule enzyme/min using 2.1 picomoles enzyme/mg mito E Escubedo Eur J Pharm 398 31 2000   *)
dnuctransvmax = 0.000086/0.0000021*culestrans;
#  deoxyguanosine transporter Vm converted from micromoles substrate/mg mito/min to molecule 
#  substrate/molecule enzyme/min using 2.1 picomoles enzyme/mg mito Life Sciences 58 753 1996 E Escubedo Eur J Pharm 398 31 2000 *)
Gdnuctransvmax=0.000015*proteinspermito*nuctransMW*culestrans;

#              (* 1.288 Vm of the first P of dT in the forward direction converting from micromoles substrate/mg enzyme/minute to molecules substrate/molecule enzyme/minute Eriksson JBC 278 6963 2003 *)
Vm1PfT = 1.288*nuckinMW* tk2culespermito;
#              (* 0.789 Vm of the first P of dC in the forward direction converting from micromoles substrate/mg enzyme/minute to molecules substrate/molecule enzyme/minute Eriksson JBC 278 6963 2003 *)
Vm1PfC = 0.789*nuckinMW* tk2culespermito;
# 0.429 Vm of the first P of dA in the forward direction converting from micromoles substrate/mg enzyme/minute to molecules substrate/molecule enzyme/minute Eriksson FEBS Lett 2003 554 319 *)
Vm1PfA = 0.429*nuckinMW* dgkculespermito;
# 0.043 Vm of the first P of dG in the forward direction converting from micromoles substrate/mg enzyme/minute to molecules substrate/molecule enzyme/minute Eriksson FEBS Lett 2003 554 319 *)
Vm1PfG = 0.043*nuckinMW*dgkculespermito; 
# 74 Vm of the first P of dT in the reverse direction  Bianchi Biochem Pharm 66 471 2003 *)
Vm1PrT = 74*tidaseMW*tidaseTpermito;
# 121 Vm of the first P of dC in the reverse direction Greger J Enzyme 25 26 1980 *)
Vm1PrC = 0.031 * 6e23 / 1e6 /5e9;
# 120 Vm of the first P of dA in the reverse direction Greger J Enzyme 25 26 1980 *)
Vm1PrA = 0.031 * 6e23 / 1e6 /5e9;
# 71 Vm of the first P of dG in the reverse direction Greger J Enzyme 25 26 1980 *)
Vm1PrG = 0.031 * 6e23 / 1e6 /5e9;
# Vm of the second P (MP to DP) of T in the forward direction converting from micromoles 
# substrate/mg enzyme/minute to molecules substrate/molecule enzyme/minute*)
Vm2PfT = 4*tmpkMW*culestmpkpermito; 
# Vm of the second P (MP to DP) of C in the forward direction converting from micromoles 
# substrate/mg enzyme/minute to molecules substrate/molecule enzyme/minute  Karlsson A BBRC 311 440 2003 *)
Vm2PfC = 23*nmpkMW*nmpkculespermito; 
# Vm of the second P (MP to DP) of A in the forward direction converting from micromoles 
# substrate/mg enzyme/minute to molecules substrate/molecule enzyme/minute*)
Vm2PfA = 1320*akMW*culesAK; 
# Vm of the second P (MP to DP) in the forward direction converting from micromoles substrate/mg 
# enzyme/minute to molecules substrate/molecule enzyme/minute value for dC Karlsson A BBRC 311 440 2003 
# Kuhn Eur J. Biochem 161 551 1986 kcat 131/s K=14 uM *)
Vm2PfG = 131*60*nmpkculespermito;
#              (* Vm of the second P of dT in the reverse direction *)
Vm2PrT = Vm2PfT*factorMD;
#              (* Vm of the second P of dC in the reverse direction *)
Vm2PrC = Vm2PfC*factorMD;      
#              (* Vm of the second P of dA in the reverse direction *)
Vm2PrA = Vm2PfA*factorMD;        
#              (* Vm of the second P of dG in the reverse direction *)
Vm2PrG = Vm2PfG*factorMD;
#Vm of the third P (DP to TP) for T in the forward direction converting from micromoles substrate/mg enzyme/minute 
# to molecules substrate/molecule enzyme/minute   Lambeth DO J. Biol. Chem 272,24604, 1997 *)
Vm3PfT = 140*ndpkMW*ndpkculespermito;
# Vm of the third P (DP to TP) for T in the forward direction converting from micromoles substrate/mg enzyme/minute to 
# molecules substrate/molecule enzyme/minute Lambeth DO J. Biol. Chem 272,24604, 1997  *)
Vm3PfC = 140*ndpkMW*ndpkculespermito;
# Vm of the third P (DP to TP) for A in the forward direction converting from micromoles substrate/mg enzyme/minute 
# to molecules substrate/molecule enzyme/minute Lambeth DO J. Biol. Chem 272,24604, 1997  *)
Vm3PfA = 325*ndpkMW*ndpkculespermito;
# Vm of the third P (DP to TP) for G in the forward direction converting from micromoles substrate/mg enzyme/minute
# to molecules substrate/molecule enzyme/minute Lambeth DO J. Biol. Chem 272,24604, 1997  *)
Vm3PfG = 325*ndpkMW*ndpkculespermito; 
#              (* Vm of the third P of dT in the reverse direction *)
Vm3PrT = Vm3PfT*factorDT;
#              (* Vm of the third P of dC in the reverse direction *)
Vm3PrC = Vm3PfC*factorDT;              
#              (* Vm of the third P of dA in the reverse direction *)
Vm3PrA =  Vm3PfA* factorDT;
#              (* Vm of the third P of dG in the reverse direction *)
Vm3PrG = Vm3PfG* factorDT;

#              (* 16 Km of the first P of dT in the forward direction  JBC 266, 9032 1991 Eriksson  *)
K1PfT = 16*conversion; 
#              (* 11 Km of the first P of dC in the forward direction  JBC 278, 6963 2003 Eriksson *)
K1PfC = 11*conversion; 
#              (* 467 Km of the first P of dA in the forward direction  Eriksoon FEBS Lett 554, 319 2003 *)
K1PfA = 467*conversion;
#              (* 4 Km of the first P of dG in the forward direction Eriksoon FEBS Lett 554, 319 2003 *)
K1PfG = 4*conversion;
#              (* 200 Km of the first P of dT in the reverse direction Bianchi Biochem Pharm 66, 471 2003 *)
K1PrT = 200 *conversion;
#              (* 1900 Km of the first P of dC in the reverse direction Greger Enzyme 25, 26 1980 *)
K1PrC = 1900*conversion;
#              (* 1520 Km of the first P of dA in the reverse direction Greger Enzyme 25, 26 1980 *)
K1PrA = 1520*conversion;         
#              (* 560 Km of the first P of dG in the reverse direction Greger Enzyme 25, 26 1980 *)
K1PrG = 560*conversion;
# 4.9 Km of the second P (MP to DP) of T in the forward direction Tamiya N Bichim Biophys Acta 995, 28 1989  Furman PA PNAS 83, 8333 1986
K2PfT = 4.9*conversion;                       
#              (* 17 Km of the second P (MP to DP) of C in the forward direction Karlsson BBRC 311, 440 2003 JBC 259,12346 1984 *)
K2PfC = K2Pf = 17*conversion;
#              (* 139 Km of the second P (MP to DP) of A in the forward direction Noda LH Eur J. Biochem 93, 263 1979 *)
K2PfA = 139*conversion;
#              (* 17 Km of the second P (MP to DP) of G in the forward direction no reference used same as dCMP kinase *)
K2PfG = 17*conversion;
#              (* Km of the second P (MP to DP) of T in the reverse direction *)
K2PrT = K2PfT;
#              (* Km of the second P (MP to DP) of C in the reverse direction *)
K2PrC = K2PfC;
#              (* Km of the second P (MP to DP) of A in the reverse direction *)
K2PrA = K2PfA;
#              (* Km of the second P (MP to DP) of G in the reverse direction *)
K2PrG = K2PfG;
#              (* 300 Km of the third P (DP to TP) of T in the forward direction Lambeth DO J. Biol. Chem 272,24604, 1997 *)
K3PfT = 300*conversion;
#              (* 300 Km of the third P (DP to TP) of C in the forward direction Lambeth DO J. Biol. Chem 272,24604, 1997 *)
K3PfC = 300*conversion;
#              (* 75 Km of the third P (DP to TP) of A or G in the forward direction Lambeth DO J. Biol. Chem 272,24604, 1997 *)
K3PfA = 75*conversion;
#              (* 75 Km of the third P (DP to TP) of G in the forward direction Lambeth DO J. Biol. Chem 272,24604, 1997 *)
K3PfG = 75*conversion;
#              (* Km of the third P (DP to TP) of T in the reverse direction *)
K3PrT = K3PfT;
#              (* Km of the third P (DP to TP) of C in the reverse direction *)
K3PrC = K3PfC;
#              (* Km of the third P (DP to TP) of A in the reverse direction *)
K3PrA = K3PfA;
#              (* Km of the third P (DP to TP) of G in the reverse direction *)
K3PrG = K3PfG;
#              (* initial concentration of dnucleosides dXn and dnucleotides dNXP in the mitochondrial matrix set equal to cplasm  *)

#              (* Concentration of dnucleotide di and tri phosphates dNXP in the cplasm Traut Mol Cell Biochem 140, 1 1994 *)
#              (* cancer cytoplasmic dNDP and dNTP levels *)
if(celltype==1) dTTPc =  35 * conversion
if(celltype==1) dCTPc =  20 * conversion
if(celltype==1) dATPc =  25 * conversion
if(celltype==1) dGTPc =  20 * conversion
if(celltype==1) dTDPc =  dTTPc/5
if(celltype==1) dCDPc =  dCTPc/5
if(celltype==1) dADPc =  dATPc/5
if(celltype==1) dGDPc =  dGTPc/5
#              (* cplasmic dNDP and dNTP concentrations in rapidly dividing levels *) 
if(celltype==2) dTTPc =  20 * conversion
if(celltype==2) dCTPc =  10 * conversion
if(celltype==2) dATPc =  10 * conversion
if(celltype==2) dGTPc =  5  * conversion
if(celltype==2) dTDPc=  dTTPc * 1/5
if(celltype==2) dCDPc = dCTPc * 1/5
if(celltype==2) dADPc = dATPc * 1/5
if(celltype==2) dGDPc = dGTPc * 1/5
#  cplasmic dNDP and dNTP concentration in slowly dividing or resting cells  
#Sadee W J Chrom. 188)  149 1980 Gaillard RK ANtimicrob Ag. Chemother 46)  1005 2002 LePoivre Anal Biochem 269)  
#403 1999 Rabes BBA 992)  349 1989  *)
if(celltype==3) dTTPc =  2.5  * conversion
if(celltype==3) dCTPc =  1.0  * conversion
if(celltype==3) dATPc =  0.7  * conversion
if(celltype==3) dGTPc =  0.8  * conversion
if(celltype==3) dTDPc=  dTTPc
if(celltype==3) dCDPc = dCTPc
if(celltype==3) dADPc = dATPc
if(celltype==3) dGDPc = dGTPc
#              (* cplasmic dNDP and dNTP concentrations in postmitotic cells *)
if(celltype==4) dTTPc =  0.25  * conversion
if(celltype==4) dCTPc =  0.10  * conversion
if(celltype==4) dATPc =  0.07  * conversion
if(celltype==4) dGTPc =  0.08  * conversion
if(celltype==4) dTDPc=  dTTPc
if(celltype==4) dCDPc = dCTPc
if(celltype==4) dADPc = dATPc
if(celltype==4) dGDPc = dGTPc

dT0  =  dTc;
dTMP0 = dTDPc/2;
dTDP0 = dTDPc;
dTTP0 = dTTPc;
dC0  =  dCc;
dCMP0 = dCDPc/2;
dCDP0 = dCDPc;
dCTP0 = dCTPc;
dG0  =  dGc;  
dGMP0 = dGDPc/2;
dGDP0 = dGDPc;
dGTP0 = dGTPc;
dA0  =  dAc;
dAMP0 = dADPc/2;
dADP0 = dADPc;
dATP0 = dATPc;
DNA0  = 0;
LDNA0 = 0;
HDNA0 = 0;

# Ks for exchanger J Walker, F Palmieri PNAS 98, 2284 2000 *)
# deoxynucleotide exchanger Vm converted from micromoles substrate/mg enzyme/min 
# to molecules substrate/molecule enzyme/min Walker PNAS 2000 *)
VmDnc = 0.85*exMW*expermito;
#              (* deoxynucleotide exchanger K's J Walker, F Palmieri PNAS 98, 2284 2000 *)
#          (* deoxynucleotide exchanger Km for dTDP *)
exKdT = 117*conversion;
#(* deoxynucleotide exchanger Km for dCDP *)
exKdC = 99*conversion;
#(* deoxynucleotide exchanger Km for dADP *)
exKdA = 14*conversion;
#(* deoxynucleotide exchanger Km for dGDP *)
exKdG = 55*conversion;

#(* deoxynucleotide exchanger Km for dTTP *)
exKtT = 595*conversion;
#(* deoxynucleotide exchanger Km for dCTP *)
exKtC = 423*conversion;
#(* deoxynucleotide exchanger Km for dTDP *)
exKtA = 106*conversion;
#(* deoxynucleotide exchanger Km for dGTP *)
exKtG = 230*conversion;

#(* deoxynucleotide exchanger Km's for deoxynucleotide diphosphates on the matrix side of the inner membrane *)
exKdAin=exKdA;
exKdTin=exKdT;
exKdCin=exKdC;
exKdGin=exKdG;
exKtAin=exKtA;
exKtTin=exKtT;
exKtCin=exKtC;
exKtGin=exKtG;


# define ODE right hand side
fmito <- function(t, X, p)
{ # bring states back to normal names
  dT=X[1];dTMP=X[2];dTDP=X[3];dTTP=X[4];
  dC=X[5];dCMP=X[6];dCDP=X[7];dCTP=X[8];
  dG=X[9];dGMP=X[10];dGDP=X[11];dGTP=X[12];
  dA=X[13];dAMP=X[14];dADP=X[15];dATP=X[16];
  LDNA=X[17];HDNA=X[18];
  Lon=p["Lon"]
  Hon=p["Hon"]
# now define the fluxes
  Tnucsidetransport=dnuctransvmax*(dTc/dnuctransK/(1+dTc/dnuctransK+dCc/dnuctransK+dAc/dnuctransK) 
        -  dT/dnuctransK/(1+dT/dnuctransK+dC/dnuctransK+dA/dnuctransK))
  Tnucleosidekinase=Vm1PfT*(dT/K1PfT)^tk2hill/(1+(dT/K1PfT)^tk2hill+dTTP/KiTtptk2+dC/KiCtk2+dCTP/KiCtptk2)
  Tnucleotidase=Vm1PrT*dTMP/K1PrT/(1+dTMP/K1PrT)
  Tmpkforward=Vm2PfT*dTMP/K2PfT/(1+dTMP/K2PfT+dTDP/K2PfT)
  Tmpkreverse=Vm2PrT*dTDP/K2PrT/(1+dTDP/K2PrT+dTMP/K2PrT)
  Tdpkforward=Vm3PfT*dTDP/K3PfT/(1+dTDP/K3PfT+dCDP/K3PfC+dGDP/K3PfG+dADP/K3PfA+rADPin/KiADPNDPK+dTTP/K3PrT+dCTP/K3PrC+dGTP/K3PrG+dATP/K3PrA+rATPin/KiATPNDPK)
  Tdpkreverse=Vm3PrT*dTTP/K3PrT/(1+dTTP/K3PrT+dCTP/K3PrC+dGTP/K3PrG+dATP/K3PrA+rATPin/KiATPNDPK+dTDP/K3PfT+dCDP/K3PfC+dGDP/K3PfG+dADP/K3PfA+rADPin/KiADPNDPK) 
  dTDPexchangein=VmDnc*dTDPc/exKdT/(1+dTTPc/exKtT+dCTPc/exKtC+dGTPc/exKtG+dATPc/exKtA+dTDPc/exKdT+dCDPc/exKdC+dGDPc/exKdG+dADPc/exKdA)
  dTTPexchangein=VmDnc*dTTPc/exKtT/(1+dTTPc/exKtT+dCTPc/exKtC+dGTPc/exKtG+dATPc/exKtA+dTDPc/exKdT+dCDPc/exKdC+dGDPc/exKdG+dADPc/exKdA)
  dTDPexchangeout=VmDnc*dTDP/exKdTin/(1+dTTP/exKtTin+dCTP/exKtCin+dGTP/exKtGin+dATP/exKtAin+dTDP/exKdTin+dCDP/exKdCin+dGDP/exKdGin+dADP/exKdAin)
  dTTPexchangeout=VmDnc*dTTP/exKtTin/(1+dTTP/exKtTin+dCTP/exKtCin+dGTP/exKtGin+dATP/exKtAin+dTDP/exKdTin+dCDP/exKdCin+dGDP/exKdGin+dADP/exKdAin)
  Cnucsidetransport=dnuctransvmax*(dCc/dnuctransK/(1+dTc/dnuctransK+dCc/dnuctransK+dAc/dnuctransK) 
        -  dC/dnuctransK/(1+dT/dnuctransK+dC/dnuctransK+dA/dnuctransK))
  Cnucleosidekinase=Vm1PfC*dC/K1PfC/(1+dC/K1PfC+dT/KiTtk2+dTTP/KiTtptk2+dCTP/KiCtptk2)
  Cnucleotidase=Vm1PrC*dCMP/K1PrC/(1+dCMP/K1PrC+rAMPin/KiAMPtidase+rADPin/KiADPtidase+rATPin/KiATPtidase)
  Cmpkforward=Vm2PfC*dCMP/K2PfC/(1+dCMP/K2PfC+UMPin/KiUMP+rCMPin/KiCMP+dCDP/K2PrC+UDPin/KiUDP+rCDPin/KiCDP)
  Cmpkreverse=Vm2PrC*dCDP/K2PrC/(1+dCDP/K2PrC+UDPin/KiUDP+rCDPin/KiCDP+dCMP/K2PfC+UMPin/KiUMP+rCMPin/KiCMP)
  Cdpkforward=Vm3PfC*dCDP/K3PfC/(1+dTDP/K3PfT+dCDP/K3PfC+dGDP/K3PfG+dADP/K3PfA+rADPin/KiADPNDPK+dTTP/K3PrT+dCTP/K3PrC+dGTP/K3PrG+dATP/K3PrA+rATPin/KiATPNDPK)
  Cdpkreverse=Vm3PrC*dCTP/K3PrC/(1+dTTP/K3PrT+dCTP/K3PrC+dGTP/K3PrG+dATP/K3PrA+rATPin/KiATPNDPK+dTDP/K3PfT+dCDP/K3PfC+dGDP/K3PfG+dADP/K3PfA+rADPin/KiADPNDPK) 
  dCDPexchangein=VmDnc*dCDPc/exKdC/(1+dTTPc/exKtT+dCTPc/exKtC+dGTPc/exKtG+dATPc/exKtA+dTDPc/exKdT+dCDPc/exKdC+dGDPc/exKdG+dADPc/exKdA)
  dCTPexchangein=VmDnc*dCTPc/exKtC/(1+dTTPc/exKtT+dCTPc/exKtC+dGTPc/exKtG+dATPc/exKtA+dTDPc/exKdT+dCDPc/exKdC+dGDPc/exKdG+dADPc/exKdA)
  dCDPexchangeout=VmDnc*dCDP/exKdCin/(1+dTTP/exKtTin+dCTP/exKtCin+dGTP/exKtGin+dATP/exKtAin+dTDP/exKdTin+dCDP/exKdCin+dGDP/exKdGin+dADP/exKdAin)
  dCTPexchangeout=VmDnc*dCTP/exKtCin/(1+dTTP/exKtTin+dCTP/exKtCin+dGTP/exKtGin+dATP/exKtAin+dTDP/exKdTin+dCDP/exKdCin+dGDP/exKdGin+dADP/exKdAin)
  Gnucsidetransport=Gdnuctransvmax*(dGc/GdnuctransK/(1+dGc/GdnuctransK) - dG/GdnuctransK/(1+dG/GdnuctransK))
  Gnucleosidekinase=Vm1PfG*dG/K1PfG/(1+dG/K1PfG+dA/KiAdgk+dGTP/KiGtpdgk+dATP/KiAtpdgk+dGMP/KiGmpdgk+dAMP/KiAmpdgk) # I doubt the last 2 monos inhibit
  Gnucleotidase=Vm1PrG*dGMP/K1PrG/(1+dGMP/K1PrG+rAMPin/KiAMPtidase+rADPin/KiADPtidase+rATPin/KiATPtidase)
  Gmpkforward=Vm2PfG*dGMP/K2PfG/(1+dGMP/K2PfG+rGMPin/KiGMP+dGDP/K2PrG+rGDPin/KiGDP)
  Gmpkreverse=Vm2PrG*dGDP/K2PrG/(1+dGDP/K2PrG+rGDPin/KiGDP+dGMP/K2PfG+rGMPin/KiGMP)
  Gdpkforward=Vm3PfG*dGDP/K3PfG/(1+dTDP/K3PfT+dCDP/K3PfC+dGDP/K3PfG+dADP/K3PfA+rADPin/KiADPNDPK+dTTP/K3PrT+dCTP/K3PrC+dGTP/K3PrG+dATP/K3PrA+rATPin/KiATPNDPK)
  Gdpkreverse=Vm3PrG*dGTP/K3PrG/(1+dTTP/K3PrT+dCTP/K3PrC+dGTP/K3PrG+dATP/K3PrA+rATPin/KiATPNDPK+dTDP/K3PfT+dCDP/K3PfC+dGDP/K3PfG+dADP/K3PfA+rADPin/KiADPNDPK) 
  dGDPexchangein=VmDnc*dGDPc/exKdG/(1+dTTPc/exKtT+dCTPc/exKtC+dGTPc/exKtG+dATPc/exKtA+dTDPc/exKdT+dCDPc/exKdC+dGDPc/exKdG+dADPc/exKdA)
  dGTPexchangein=VmDnc*dGTPc/exKtG/(1+dTTPc/exKtT+dCTPc/exKtC+dGTPc/exKtG+dATPc/exKtA+dTDPc/exKdT+dCDPc/exKdC+dGDPc/exKdG+dADPc/exKdA)
  dGDPexchangeout=VmDnc*dGDP/exKdGin/(1+dTTP/exKtTin+dCTP/exKtCin+dGTP/exKtGin+dATP/exKtAin+dTDP/exKdTin+dCDP/exKdCin+dGDP/exKdGin+dADP/exKdAin)
  dGTPexchangeout=VmDnc*dGTP/exKtGin/(1+dTTP/exKtTin+dCTP/exKtCin+dGTP/exKtGin+dATP/exKtAin+dTDP/exKdTin+dCDP/exKdCin+dGDP/exKdGin+dADP/exKdAin)
  Anucsidetransport=dnuctransvmax*(dAc/dnuctransK/(1+dTc/dnuctransK+dCc/dnuctransK+dAc/dnuctransK) 
        -  dA/dnuctransK/(1+dT/dnuctransK+dC/dnuctransK+dA/dnuctransK))
  Anucleosidekinase=Vm1PfA*dA/K1PfA/(1+dG/KiGdgk+dA/K1PfA+dATP/KiAtpdgk+dGTP/KiGtpdgk+dGMP/KiGmpdgk+dAMP/KiAmpdgk)
  Anucleotidase=Vm1PrA*dAMP/K1PrA/(1+dAMP/K1PrA+rAMPin/KiAMPtidase+rADPin/KiADPtidase+rATPin/KiATPtidase)
  Ampkforward=Vm2PfA*dAMP/K2PfA/(1+dAMP/K2PfA+rAMPin/K2PfA+dADP/K2PrA+rADPin/K2PrA)
  Ampkreverse=Vm2PrA*dADP/K2PrA/(1+dADP/K2PrA+rADPin/K2PrA+dAMP/K2PfA+rAMPin/K2PfA)
  Adpkforward=Vm3PfA*dADP/K3PfA/(1+dTDP/K3PfT+dCDP/K3PfC+dGDP/K3PfG+dADP/K3PfA+rADPin/KiADPNDPK+dTTP/K3PrT+dCTP/K3PrC+dGTP/K3PrG+dATP/K3PrA+rATPin/KiATPNDPK)
  Adpkreverse=Vm3PrA*dATP/K3PrA/(1+dTTP/K3PrT+dCTP/K3PrC+dGTP/K3PrG+dATP/K3PrA+rATPin/KiATPNDPK+dTDP/K3PfT+dCDP/K3PfC+dGDP/K3PfG+dADP/K3PfA+rADPin/KiADPNDPK) 
  dADPexchangein=VmDnc*dADPc/exKdA/(1+dTTPc/exKtT+dCTPc/exKtC+dGTPc/exKtG+dATPc/exKtA+dTDPc/exKdT+dCDPc/exKdC+dGDPc/exKdG+dADPc/exKdA)
  dATPexchangein=VmDnc*dATPc/exKtA/(1+dTTPc/exKtT+dCTPc/exKtC+dGTPc/exKtG+dATPc/exKtA+dTDPc/exKdT+dCDPc/exKdC+dGDPc/exKdG+dADPc/exKdA)
  dADPexchangeout=VmDnc*dADP/exKdAin/(1+dTTP/exKtTin+dCTP/exKtCin+dGTP/exKtGin+dATP/exKtAin+dTDP/exKdTin+dCDP/exKdCin+dGDP/exKdGin+dADP/exKdAin)
  dATPexchangeout=VmDnc*dATP/exKtAin/(1+dTTP/exKtTin+dCTP/exKtCin+dGTP/exKtGin+dATP/exKtAin+dTDP/exKdTin+dCDP/exKdCin+dGDP/exKdGin+dADP/exKdAin)
  rT=VmPolT*dTTP/KpolT/(1+dTTP/KpolT)
  rC=VmPolC*dCTP/KpolC/(1+dCTP/KpolC)
  rG=VmPolG*dGTP/KpolG/(1+dGTP/KpolG)
  rA=VmPolA*dATP/KpolA/(1+dATP/KpolA)
  PolrateL=Lon*rT*rA*rG*rC/((fTL*rC*rA*rG)+(fCL*rT*rA*rG)+(fAL*rT*rC*rG)+ (fGL*rT*rC*rA)+0.1)
  PolrateH=Hon*rT*rA*rG*rC/((fTH*rC*rA*rG)+(fCH*rT*rA*rG)+(fAH*rT*rC*rG)+ (fGH*rT*rC*rA)+0.1)
  # now define dX/dt as fluxes into a node minus fluxes out
  dTp=Tnucsidetransport-Tnucleosidekinase+Tnucleotidase
  dTMPp=Tnucleosidekinase-Tnucleotidase-Tmpkforward+Tmpkreverse
  dTDPp=Tmpkforward-Tmpkreverse-Tdpkforward+Tdpkreverse+dTDPexchangein-dTDPexchangeout
  dTTPp=Tdpkforward-Tdpkreverse+dTTPexchangein-dTTPexchangeout- PolrateL*fTL- PolrateH*fTH
  dCp=Cnucsidetransport-Cnucleosidekinase+Cnucleotidase
  dCMPp=Cnucleosidekinase-Cnucleotidase-Cmpkforward+Cmpkreverse
  dCDPp=Cmpkforward-Cmpkreverse-Cdpkforward+Cdpkreverse+dCDPexchangein-dCDPexchangeout
  dCTPp=Cdpkforward-Cdpkreverse+dCTPexchangein-dCTPexchangeout- PolrateL*fCL-PolrateH*fCH
  dGp=Gnucsidetransport-Gnucleosidekinase+Gnucleotidase
  dGMPp=Gnucleosidekinase-Gnucleotidase-Gmpkforward+Gmpkreverse
  dGDPp=Gmpkforward-Gmpkreverse-Gdpkforward+Gdpkreverse+dGDPexchangein-dGDPexchangeout
  dGTPp=Gdpkforward-Gdpkreverse+dGTPexchangein-dGTPexchangeout-  PolrateL*fGL-  PolrateH*fGH
  dAp=Anucsidetransport-Anucleosidekinase+Anucleotidase
  dAMPp=Anucleosidekinase-Anucleotidase-Ampkforward+Ampkreverse
  dADPp=Ampkforward-Ampkreverse-Adpkforward+Adpkreverse+dADPexchangein-dADPexchangeout
  dATPp=Adpkforward-Adpkreverse+dATPexchangein-dATPexchangeout- PolrateL*fAL-PolrateH*fAH
  LDNAp=PolrateL 
  HDNAp=PolrateH   
  XP = c(dTp,dTMPp,dTDPp,dTTPp,dCp,dCMPp,dCDPp,dCTPp,dGp,dGMPp,dGDPp,dGTPp,dAp,dAMPp,dADPp,dATPp,LDNAp,HDNAp);
  V=c(Tnucsidetransport,Tnucleosidekinase,Tnucleotidase,Tmpkforward,Tmpkreverse,Tdpkforward,Tdpkreverse,
      dTDPexchangein,dTDPexchangeout,dTTPexchangein,dTTPexchangeout,
      Cnucsidetransport,Cnucleosidekinase,Cnucleotidase,Cmpkforward,Cmpkreverse,Cdpkforward,Cdpkreverse,
      dCDPexchangein,dCDPexchangeout,dCTPexchangein,dCTPexchangeout,
      Gnucsidetransport,Gnucleosidekinase,Gnucleotidase,Gmpkforward,Gmpkreverse,Gdpkforward,Gdpkreverse,
      dGDPexchangein,dGDPexchangeout,dGTPexchangein,dGTPexchangeout,
      Anucsidetransport,Anucleosidekinase,Anucleotidase,Ampkforward,Ampkreverse,Adpkforward,Adpkreverse,
      dADPexchangein,dADPexchangeout,dATPexchangein,dATPexchangeout,
      PolrateL,PolrateH)
  names(V)<-c("Tnucsidetransport","Tnucleosidekinase","Tnucleotidase","Tnmpkforward","Tnmpkreverse","Tndpkforward","Tndpkreverse",
      "dTDPexchangein","dTDPexchangeout","dTTPexchangein","dTTPexchangeout",
      "Cnucsidetransport","Cnucleosidekinase","Cnucleotidase","Cmpkforward","Cmpkreverse","Cdpkforward","Cdpkreverse",
      "dCDPexchangein","dCDPexchangeout","dCTPexchangein","dCTPexchangeout",
      "Gnucsidetransport","Gnucleosidekinase","Gnucleotidase","Gmpkforward","Gmpkreverse","Gdpkforward","Gdpkreverse",
      "dGDPexchangein","dGDPexchangeout","dGTPexchangein","dGTPexchangeout",
      "Anucsidetransport","Anucleosidekinase","Anucleotidase","Ampkforward","Ampkreverse","Adpkforward","Adpkreverse",
      "dADPexchangein","dADPexchangeout","dATPexchangein","dATPexchangeout",
      "PolrateL","PolrateH")
  list(XP,V)
}

