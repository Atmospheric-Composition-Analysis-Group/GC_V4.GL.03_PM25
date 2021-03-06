! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
! 
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : gckpp_Parameters.f90
! Time                 : Fri Aug 12 14:12:49 2016
! Working directory    : /misc/data9/chili/kpp-2.2.3/perl
! Equation file        : gckpp.kpp
! Output root filename : gckpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE gckpp_Parameters

  USE gckpp_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 204 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 193 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 127 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 11 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 514 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 194 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 1688 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 2014 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 194 
! CNEQN - (NREACT+1) Number stoicm elements in compressed col format
  INTEGER, PARAMETER :: CNEQN = 515 
! NHESS - Length of Sparse Hessian
  INTEGER, PARAMETER :: NHESS = 1435 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 204 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_MONX = 1 
  INTEGER, PARAMETER :: ind_DRYACET = 2 
  INTEGER, PARAMETER :: ind_DRYCH2O = 3 
  INTEGER, PARAMETER :: ind_DRYH2O2 = 4 
  INTEGER, PARAMETER :: ind_DRYHNO3 = 5 
  INTEGER, PARAMETER :: ind_DRYN2O5 = 6 
  INTEGER, PARAMETER :: ind_DRYNO2 = 7 
  INTEGER, PARAMETER :: ind_DRYO3 = 8 
  INTEGER, PARAMETER :: ind_DRYPAN = 9 
  INTEGER, PARAMETER :: ind_DRYIPMN = 10 
  INTEGER, PARAMETER :: ind_DRYNPMN = 11 
  INTEGER, PARAMETER :: ind_DRYPMNN = 12 
  INTEGER, PARAMETER :: ind_DRYPPN = 13 
  INTEGER, PARAMETER :: ind_DRYR4N2 = 14 
  INTEGER, PARAMETER :: ind_DRYMVK = 15 
  INTEGER, PARAMETER :: ind_DRYMACR = 16 
  INTEGER, PARAMETER :: ind_DRYHAC = 17 
  INTEGER, PARAMETER :: ind_DRYALD2 = 18 
  INTEGER, PARAMETER :: ind_SO4 = 19 
  INTEGER, PARAMETER :: ind_MSA = 20 
  INTEGER, PARAMETER :: ind_LISOPOH = 21 
  INTEGER, PARAMETER :: ind_LISOPNO3 = 22 
  INTEGER, PARAMETER :: ind_LBRO2H = 23 
  INTEGER, PARAMETER :: ind_LBRO2N = 24 
  INTEGER, PARAMETER :: ind_LTRO2H = 25 
  INTEGER, PARAMETER :: ind_LTRO2N = 26 
  INTEGER, PARAMETER :: ind_LXRO2H = 27 
  INTEGER, PARAMETER :: ind_LXRO2N = 28 
  INTEGER, PARAMETER :: ind_LNRO2H = 29 
  INTEGER, PARAMETER :: ind_LNRO2N = 30 
  INTEGER, PARAMETER :: ind_NRO2 = 31 
  INTEGER, PARAMETER :: ind_NAP = 32 
  INTEGER, PARAMETER :: ind_DRYHOBr = 33 
  INTEGER, PARAMETER :: ind_DRYHBr = 34 
  INTEGER, PARAMETER :: ind_DRYBrNO3 = 35 
  INTEGER, PARAMETER :: ind_DRYBr2 = 36 
  INTEGER, PARAMETER :: ind_PYAC = 37 
  INTEGER, PARAMETER :: ind_DRYISOPND = 38 
  INTEGER, PARAMETER :: ind_DRYISOPNB = 39 
  INTEGER, PARAMETER :: ind_DRYRIP = 40 
  INTEGER, PARAMETER :: ind_DRYIEPOX = 41 
  INTEGER, PARAMETER :: ind_DRYMACRN = 42 
  INTEGER, PARAMETER :: ind_DRYMVKN = 43 
  INTEGER, PARAMETER :: ind_DRYPROPNN = 44 
  INTEGER, PARAMETER :: ind_DRYHCOOH = 45 
  INTEGER, PARAMETER :: ind_DRYACTA = 46 
  INTEGER, PARAMETER :: ind_DRYLVOC = 47 
  INTEGER, PARAMETER :: ind_DRYGLYX = 48 
  INTEGER, PARAMETER :: ind_DRYMGLY = 49 
  INTEGER, PARAMETER :: ind_DRYGLYC = 50 
  INTEGER, PARAMETER :: ind_DRYISN1 = 51 
  INTEGER, PARAMETER :: ind_INDIOL = 52 
  INTEGER, PARAMETER :: ind_ISN1OA = 53 
  INTEGER, PARAMETER :: ind_ISN1OG = 54 
  INTEGER, PARAMETER :: ind_LVOCOA = 55 
  INTEGER, PARAMETER :: ind_SOAGX = 56 
  INTEGER, PARAMETER :: ind_SOAIE = 57 
  INTEGER, PARAMETER :: ind_SOAMG = 58 
  INTEGER, PARAMETER :: ind_SOAME = 59 
  INTEGER, PARAMETER :: ind_DRYIMAE = 60 
  INTEGER, PARAMETER :: ind_ISOPOH = 61 
  INTEGER, PARAMETER :: ind_ISOPO3 = 62 
  INTEGER, PARAMETER :: ind_ISOPNO3 = 63 
  INTEGER, PARAMETER :: ind_RIO2NO = 64 
  INTEGER, PARAMETER :: ind_RIO2HO2 = 65 
  INTEGER, PARAMETER :: ind_RIO2OTH = 66 
  INTEGER, PARAMETER :: ind_PSOA1 = 67 
  INTEGER, PARAMETER :: ind_PSOA2 = 68 
  INTEGER, PARAMETER :: ind_CO2 = 69 
  INTEGER, PARAMETER :: ind_DRYDEP = 70 
  INTEGER, PARAMETER :: ind_LVOC = 71 
  INTEGER, PARAMETER :: ind_IMAE = 72 
  INTEGER, PARAMETER :: ind_TOLU = 73 
  INTEGER, PARAMETER :: ind_XYLE = 74 
  INTEGER, PARAMETER :: ind_CHBr3 = 75 
  INTEGER, PARAMETER :: ind_CH2Br2 = 76 
  INTEGER, PARAMETER :: ind_CH3Br = 77 
  INTEGER, PARAMETER :: ind_BENZ = 78 
  INTEGER, PARAMETER :: ind_DHDN = 79 
  INTEGER, PARAMETER :: ind_BrNO2 = 80 
  INTEGER, PARAMETER :: ind_PMNN = 81 
  INTEGER, PARAMETER :: ind_PPN = 82 
  INTEGER, PARAMETER :: ind_BRO2 = 83 
  INTEGER, PARAMETER :: ind_TRO2 = 84 
  INTEGER, PARAMETER :: ind_N2O5 = 85 
  INTEGER, PARAMETER :: ind_XRO2 = 86 
  INTEGER, PARAMETER :: ind_IEPOX = 87 
  INTEGER, PARAMETER :: ind_ALK4 = 88 
  INTEGER, PARAMETER :: ind_HNO2 = 89 
  INTEGER, PARAMETER :: ind_MAP = 90 
  INTEGER, PARAMETER :: ind_MPN = 91 
  INTEGER, PARAMETER :: ind_IMAO3 = 92 
  INTEGER, PARAMETER :: ind_DMS = 93 
  INTEGER, PARAMETER :: ind_C3H8 = 94 
  INTEGER, PARAMETER :: ind_ETP = 95 
  INTEGER, PARAMETER :: ind_HNO4 = 96 
  INTEGER, PARAMETER :: ind_MP = 97 
  INTEGER, PARAMETER :: ind_RA3P = 98 
  INTEGER, PARAMETER :: ind_RB3P = 99 
  INTEGER, PARAMETER :: ind_HOBr = 100 
  INTEGER, PARAMETER :: ind_Br2 = 101 
  INTEGER, PARAMETER :: ind_HBr = 102 
  INTEGER, PARAMETER :: ind_RP = 103 
  INTEGER, PARAMETER :: ind_PAN = 104 
  INTEGER, PARAMETER :: ind_C2H6 = 105 
  INTEGER, PARAMETER :: ind_H2O2 = 106 
  INTEGER, PARAMETER :: ind_BrNO3 = 107 
  INTEGER, PARAMETER :: ind_PP = 108 
  INTEGER, PARAMETER :: ind_PRPN = 109 
  INTEGER, PARAMETER :: ind_ATOOH = 110 
  INTEGER, PARAMETER :: ind_R4P = 111 
  INTEGER, PARAMETER :: ind_HC187 = 112 
  INTEGER, PARAMETER :: ind_RIP = 113 
  INTEGER, PARAMETER :: ind_VRP = 114 
  INTEGER, PARAMETER :: ind_IAP = 115 
  INTEGER, PARAMETER :: ind_MRP = 116 
  INTEGER, PARAMETER :: ind_MOBA = 117 
  INTEGER, PARAMETER :: ind_MAOP = 118 
  INTEGER, PARAMETER :: ind_DHMOB = 119 
  INTEGER, PARAMETER :: ind_NPMN = 120 
  INTEGER, PARAMETER :: ind_INPN = 121 
  INTEGER, PARAMETER :: ind_ISNP = 122 
  INTEGER, PARAMETER :: ind_ETHLN = 123 
  INTEGER, PARAMETER :: ind_MACRNO2 = 124 
  INTEGER, PARAMETER :: ind_ROH = 125 
  INTEGER, PARAMETER :: ind_MOBAOO = 126 
  INTEGER, PARAMETER :: ind_DIBOO = 127 
  INTEGER, PARAMETER :: ind_IPMN = 128 
  INTEGER, PARAMETER :: ind_MVKOO = 129 
  INTEGER, PARAMETER :: ind_CH3CHOO = 130 
  INTEGER, PARAMETER :: ind_ACET = 131 
  INTEGER, PARAMETER :: ind_GAOO = 132 
  INTEGER, PARAMETER :: ind_ISNOHOO = 133 
  INTEGER, PARAMETER :: ind_MGLYOO = 134 
  INTEGER, PARAMETER :: ind_MVKN = 135 
  INTEGER, PARAMETER :: ind_BrO = 136 
  INTEGER, PARAMETER :: ind_ISOP = 137 
  INTEGER, PARAMETER :: ind_PRPE = 138 
  INTEGER, PARAMETER :: ind_ISNOOB = 139 
  INTEGER, PARAMETER :: ind_MGLOO = 140 
  INTEGER, PARAMETER :: ind_HNO3 = 141 
  INTEGER, PARAMETER :: ind_ISOPNB = 142 
  INTEGER, PARAMETER :: ind_A3O2 = 143 
  INTEGER, PARAMETER :: ind_IEPOXOO = 144 
  INTEGER, PARAMETER :: ind_CH2OO = 145 
  INTEGER, PARAMETER :: ind_GLYX = 146 
  INTEGER, PARAMETER :: ind_MACRN = 147 
  INTEGER, PARAMETER :: ind_PROPNN = 148 
  INTEGER, PARAMETER :: ind_MACROO = 149 
  INTEGER, PARAMETER :: ind_PO2 = 150 
  INTEGER, PARAMETER :: ind_ISNOOA = 151 
  INTEGER, PARAMETER :: ind_MAOPO2 = 152 
  INTEGER, PARAMETER :: ind_B3O2 = 153 
  INTEGER, PARAMETER :: ind_MAN2 = 154 
  INTEGER, PARAMETER :: ind_KO2 = 155 
  INTEGER, PARAMETER :: ind_GLYC = 156 
  INTEGER, PARAMETER :: ind_ISOPND = 157 
  INTEGER, PARAMETER :: ind_HC5OO = 158 
  INTEGER, PARAMETER :: ind_VRO2 = 159 
  INTEGER, PARAMETER :: ind_PRN1 = 160 
  INTEGER, PARAMETER :: ind_ETO2 = 161 
  INTEGER, PARAMETER :: ind_RCO3 = 162 
  INTEGER, PARAMETER :: ind_ATO2 = 163 
  INTEGER, PARAMETER :: ind_R4N1 = 164 
  INTEGER, PARAMETER :: ind_ISN1 = 165 
  INTEGER, PARAMETER :: ind_NMAO3 = 166 
  INTEGER, PARAMETER :: ind_MGLY = 167 
  INTEGER, PARAMETER :: ind_MRO2 = 168 
  INTEGER, PARAMETER :: ind_HC5 = 169 
  INTEGER, PARAMETER :: ind_RIO2 = 170 
  INTEGER, PARAMETER :: ind_CH2O = 171 
  INTEGER, PARAMETER :: ind_ISOPNBO2 = 172 
  INTEGER, PARAMETER :: ind_ISOPNDO2 = 173 
  INTEGER, PARAMETER :: ind_INO2 = 174 
  INTEGER, PARAMETER :: ind_R4O2 = 175 
  INTEGER, PARAMETER :: ind_HAC = 176 
  INTEGER, PARAMETER :: ind_ALD2 = 177 
  INTEGER, PARAMETER :: ind_R4N2 = 178 
  INTEGER, PARAMETER :: ind_MACR = 179 
  INTEGER, PARAMETER :: ind_MVK = 180 
  INTEGER, PARAMETER :: ind_SO2 = 181 
  INTEGER, PARAMETER :: ind_RCHO = 182 
  INTEGER, PARAMETER :: ind_MCO3 = 183 
  INTEGER, PARAMETER :: ind_OH = 184 
  INTEGER, PARAMETER :: ind_MO2 = 185 
  INTEGER, PARAMETER :: ind_NO = 186 
  INTEGER, PARAMETER :: ind_O3 = 187 
  INTEGER, PARAMETER :: ind_NO2 = 188 
  INTEGER, PARAMETER :: ind_NO3 = 189 
  INTEGER, PARAMETER :: ind_CO = 190 
  INTEGER, PARAMETER :: ind_HO2 = 191 
  INTEGER, PARAMETER :: ind_MEK = 192 
  INTEGER, PARAMETER :: ind_Br = 193 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_ACTA = 194 
  INTEGER, PARAMETER :: ind_CH4 = 195 
  INTEGER, PARAMETER :: ind_EMISSION = 196 
  INTEGER, PARAMETER :: ind_EOH = 197 
  INTEGER, PARAMETER :: ind_H2 = 198 
  INTEGER, PARAMETER :: ind_H2O = 199 
  INTEGER, PARAMETER :: ind_HCOOH = 200 
  INTEGER, PARAMETER :: ind_MOH = 201 
  INTEGER, PARAMETER :: ind_O2 = 202 
  INTEGER, PARAMETER :: ind_RCOOH = 203 
  INTEGER, PARAMETER :: ind_DUMMY = 204 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_ACTA = 1 
  INTEGER, PARAMETER :: indf_CH4 = 2 
  INTEGER, PARAMETER :: indf_EMISSION = 3 
  INTEGER, PARAMETER :: indf_EOH = 4 
  INTEGER, PARAMETER :: indf_H2 = 5 
  INTEGER, PARAMETER :: indf_H2O = 6 
  INTEGER, PARAMETER :: indf_HCOOH = 7 
  INTEGER, PARAMETER :: indf_MOH = 8 
  INTEGER, PARAMETER :: indf_O2 = 9 
  INTEGER, PARAMETER :: indf_RCOOH = 10 
  INTEGER, PARAMETER :: indf_DUMMY = 11 

! NJVRP - Length of sparse Jacobian JVRP
  INTEGER, PARAMETER :: NJVRP = 804 

! NSTOICM - Length of Sparse Stoichiometric Matrix
  INTEGER, PARAMETER :: NSTOICM = 2150 

END MODULE gckpp_Parameters

