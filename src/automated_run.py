import config
import os
import pandas as pd
from pathlib import Path
import subprocess


def run_commands(*commands):
	os.system(' ; '.join(commands))


# for P3 grch38 b154
def makeCommands(names, data):
	cmds = []

	for name in names:
		cmd = 'python3 main.py -g ' + name + '-d ' + data
		cmds.append(cmd)

	return cmds
	
def main():
	# these are the 832 genes on chr 1 long arm
	# names = ['H3-2', 'FAM72C', 'PPIAL4E', 'NBPF15', 'PPIAL4F', 'SRGAP2B', 'LOC101929805', 'PPIAL4D', 'NBPF20', 'GPR89A', 'ANKRD35', 'ITGA10', 'PEX11B', 'RBM8A', 'LIX1L', 'POLR3GL', 'TXNIP', 'HJV', 'NBPF10', 'NOTCH2NLA', 'PPIAL4H', 'NBPF12', 'CHD1L', 'GJA5', 'GJA8', 'GPR89B', 'NBPF11', 'PPIAL4G', 'NBPF14', 'NOTCH2NLB', 'NUDT4B', 'PDE4DIP', 'NBPF9', 'NOTCH2NLC', 'NBPF19', 'PPIAL4C', 'H2BC18', 'H3C13', 'H4C14', 'H3C14', 'H2AC18', 'H2AC19', 'H3C15', 'H4C15', 'H2BC21', 'H2AC20', 'H2AC21', 'BOLA1', 'SV2A', 'SF3B4', 'OTUD7B', 'VPS45', 'PLEKHO1', 'ANP32E', 'APH1A', 'CIART', 'MRPS21', 'PRPF3', 'RPRD2', 'TARS2', 'ECM1', 'ADAMTSL4', 'MCL1', 'ENSA', 'GOLPH3L', 'HORMAD1', 'CTSS', 'CTSK', 'ARNT', 'CTXND2', 'SETDB1', 'CERS2', 'ANXA9', 'MINDY1', 'PRUNE1', 'BNIPL', 'CDC42SE1', 'MLLT11', 'GABPB2', 'SEMA6C', 'LYSMD1', 'VPS72', 'PIP5K1A', 'PSMD4', 'PI4KB', 'RFX5', 'SELENBP1', 'PSMB4', 'POGZ', 'CGN', 'TUFT1', 'SNX27', 'CELF3', 'RIIAD1', 'TDRKH', 'LINGO4', 'RORC', 'C2CD4D', 'THEM5', 'THEM4', 'S100A10', 'S100A11', 'LOC100131107', 'TCHHL1', 'TCHH', 'RPTN', 'HRNR', 'FLG', 'FLG2', 'CRNN', 'LCE5A', 'CRCT1', 'LCE3E', 'LCE3D', 'LCE3C', 'LCE3B', 'LCE3A', 'LCE2D', 'LCE2C', 'LCE2B', 'LCE2A', 'LCE4A', 'C1orf68', 'KPRP', 'LCE1F', 'LCE1E', 'LCE1D', 'LCE1C', 'LCE1B', 'LCE1A', 'LCE6A', 'SMCP', 'IVL', 'SPRR4', 'SPRR1A', 'SPRR3', 'SPRR1B', 'SPRR2D', 'SPRR2A', 'SPRR2B', 'SPRR2E', 'SPRR2F', 'SPRR2G', 'LELP1', 'PRR9', 'LORICRIN', 'PGLYRP3', 'PGLYRP4', 'S100A9', 'S100A12', 'S100A8', 'S100A7', 'S100A6', 'S100A5', 'S100A4', 'S100A3', 'S100A2', 'S100A16', 'S100A14', 'S100A13', 'NPR1', 'INTS3', 'SLC27A3', 'GATAD2B', 'DENND4B', 'CRTC2', 'SLC39A1', 'RAB13', 'RPS27', 'NUP210L', 'TPM3', 'C1orf189', 'UBAP2L', 'HAX1', 'AQP10', 'ATP8B2', 'IL6R', 'UBE2Q1', 'ADAR', 'KCNN3', 'PMVK', 'PBXIP1', 'PYGO2', 'SHC1', 'CKS1B', 'FLAD1', 'LENEP', 'ZBTB7B', 'DCST1', 'ADAM15', 'EFNA4', 'EFNA3', 'EFNA1', 'SLC50A1', 'DPM3', 'KRTCAP2', 'TRIM46', 'MUC1', 'THBS3', 'GBA', 'FAM189B', 'SCAMP3', 'CLK2', 'HCN3', 'FDPS', 'RUSC1', 'ASH1L', 'MSTO1', 'DAP3', 'GON4L', 'SYT11', 'RIT1', 'KHDC4', 'RXFP4', 'ARHGEF2', 'SSR2', 'UBQLN4', 'LAMTOR2', 'RAB25', 'MEX3A', 'LMNA', 'SEMA4A', 'SLC25A44', 'PMF1-BGLAP', 'GLMP', 'VHLL', 'CCT3', 'RHBG', 'MEF2D', 'IQGAP3', 'TTC24', 'NAXE', 'HAPLN2', 'BCAN', 'NES', 'CRABP2', 'METTL25B', 'MRPL24', 'HDGF', 'PRCC', 'NTRK1', 'PEAR1', 'ARHGEF11', 'ETV3L', 'ETV3', 'FCRL5', 'FCRL4', 'FCRL3', 'FCRL2', 'FCRL1', 'CD5L', 'KIRREL1', 'CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 'ACKR1', 'FCER1A', 'OR10J1', 'OR10J5', 'APCS', 'CRP', 'DUSP23', 'FCRL6', 'SLAMF8', 'VSIG8', 'CFAP45', 'TAGLN2', 'IGSF9', 'SLAMF9', 'PIGM', 'KCNJ10', 'KCNJ9', 'IGSF8', 'ATP1A2', 'ATP1A4', 'CASQ1', 'PEA15', 'DCAF8', 'PEX19', 'COPA', 'NCSTN', 'NHLH1', 'VANGL2', 'SLAMF6', 'CD84', 'SLAMF1', 'CD48', 'SLAMF7', 'LY9', 'CD244', 'ITLN1', 'ITLN2', 'F11R', 'TSTD1', 'USF1', 'ARHGAP30', 'NECTIN4', 'KLHDC9', 'PFDN2', 'UFC1', 'USP21', 'ADAMTS4', 'FCER1G', 'APOA2', 'NR1I3', 'PCP4L1', 'MPZ', 'SDHC', 'FCGR2A', 'HSPA6', 'FCGR3A', 'FCGR3B', 'FCGR2B', 'FCRLA', 'FCRLB', 'DUSP12', 'ATF6', 'OLFML2B', 'NOS1AP', 'SPATA46', 'C1orf226', 'SH2D1B', 'UHMK1', 'UAP1', 'DDR2', 'HSD17B7', 'CCDC190', 'RGS4', 'RGS5', 'NUF2', 'PBX1', 'LMX1A', 'RXRG', 'LRRC52', 'MGST3', 'ALDH9A1', 'TMCO1', 'UCK2', 'FAM78B', 'POGK', 'TADA1', 'MAEL', 'GPA33', 'STYXL2', 'POU2F1', 'CD247', 'CREG1', 'RCSD1', 'MPZL1', 'ADCY10', 'DCAF6', 'GPR161', 'TIPRL', 'SFT2D2', 'TBX19', 'XCL2', 'XCL1', 'DPT', 'NME7', 'SLC19A2', 'F5', 'SELP', 'SELL', 'SELE', 'C1orf112', 'KIFAP3', 'NTMT2', 'GORAB', 'PRRX1', 'MROH9', 'FMO3', 'FMO2', 'FMO1', 'FMO4', 'PRRC2C', 'MYOCOS', 'MYOC', 'VAMP4', 'METTL13', 'DNM3', 'C1orf105', 'SUCO', 'FASLG', 'TNFSF18', 'TNFSF4', 'PRDX6', 'SLC9C2', 'KLHL20', 'DARS2', 'ZBTB37', 'SERPINC1', 'RC3H1', 'RABGAP1L', 'CACYBP', 'MRPS14', 'TNN', 'KIAA0040', 'TNR', 'COP1', 'PAPPA2', 'ASTN1', 'BRINP2', 'CRYZL2P-SEC16B', 'TEX35', 'RALGPS2', 'FAM20B', 'TOR3A', 'ABL2', 'SOAT1', 'AXDND1', 'TDRD5', 'FAM163A', 'TOR1AIP2', 'TOR1AIP1', 'CEP350', 'QSOX1', 'ACBD6', 'XPR1', 'KIAA1614', 'STX6', 'MR1', 'IER5', 'CACNA1E', 'ZNF648', 'GLUL', 'TEDDM1', 'RGSL1', 'RNASEL', 'RGS16', 'RGS8', 'NPL', 'DHX9', 'SHCBP1L', 'LAMC1', 'LAMC2', 'NMNAT2', 'SMG7', 'NCF2', 'ARPC5', 'RGL1', 'COLGALT2', 'TSEN15', 'C1orf21', 'EDEM3', 'NIBAN1', 'RNF2', 'SWT1', 'IVNS1ABP', 'HMCN1', 'TPR', 'ODR4', 'PDC', 'PTGS2', 'PLA2G4A', 'BRINP3', 'RGS18', 'RGS21', 'RGS1', 'RGS13', 'RGS2', 'UCHL5', 'GLRX2', 'CDC73', 'KCNT2', 'CFH', 'CFHR3', 'CFHR1', 'CFHR4', 'CFHR2', 'CFHR5', 'F13B', 'ASPM', 'ZBTB41', 'CRB1', 'DENND1B', 'C1orf53', 'LHX9', 'NEK7', 'ATP6V1G3', 'PTPRC', 'NR5A2', 'ZNF281', 'KIF14', 'DDX59', 'CAMSAP2', 'GPR25', 'INAVA', 'KIF21B', 'CACNA1S', 'ASCL5', 'TMEM9', 'IGFN1', 'PKP1', 'TNNT2', 'LAD1', 'TNNI1', 'PHLDA3', 'CSRP1', 'NAV1', 'IPO9', 'SHISA4', 'LMOD1', 'TIMM17A', 'RNPEP', 'ELF3', 'GPR37L1', 'PTPN7', 'LGR6', 'UBE2T', 'PPP1R12B', 'KDM5B', 'RABIF', 'KLHL12', 'ADIPOR1', 'CYB5R1', 'TMEM183A', 'PPFIA4', 'MYOG', 'ADORA1', 'MYBPH', 'CHI3L1', 'CHIT1', 'BTG2', 'FMOD', 'PRELP', 'OPTC', 'ATP2B4', 'LAX1', 'ZC3H11A', 'SNRPE', 'SOX13', 'ETNK2', 'REN', 'KISS1', 'GOLT1A', 'PLEKHA6', 'PPP1R15B', 'PIK3C2B', 'MDM4', 'LRRN2', 'NFASC', 'CNTN2', 'TMEM81', 'RBBP5', 'DSTYK', 'TMCC2', 'NUAK2', 'KLHDC8A', 'LEMD1', 'CDK18', 'MFSD4A', 'ELK4', 'SLC45A3', 'NUCKS1', 'RAB29', 'SLC41A1', 'PM20D1', 'SLC26A9', 'RAB7B', 'CTSE', 'RHEX', 'AVPR1B', 'SRGAP2', 'IKBKE', 'RASSF5', 'DYRK3', 'MAPKAPK2', 'IL19', 'IL20', 'FCMR', 'PIGR', 'FCAMR', 'C1orf116', 'PFKFB2', 'C4BPB', 'C4BPA', 'CD55', 'CR2', 'CR1', 'CR1L', 'CD46', 'CD34', 'PLXNA2', 'CAMK1G', 'LAMB3', 'G0S2', 'HSD11B1', 'TRAF3IP3', 'IRF6', 'UTP25', 'SYT14', 'SERTAD4', 'HHAT', 'KCNH1', 'RCOR3', 'TRAF5', 'RD3', 'SLC30A1', 'NEK2', 'LPGAT1', 'INTS7', 'DTL', 'PPP2R5A', 'PACC1', 'NENF', 'ATF3', 'FAM71A', 'BATF3', 'NSL1', 'SPATA45', 'FLVCR1', 'VASH2', 'ANGEL2', 'RPS6KC1', 'PROX1', 'SMYD2', 'PTPN14', 'CENPF', 'KCNK2', 'KCTD3', 'USH2A', 'ESRRG', 'GPATCH2', 'SPATA17', 'RRP15', 'TGFB2', 'LYPLAL1', 'ZC3H11B', 'SLC30A10', 'EPRS1', 'BPNT1', 'IARS2', 'RAB3GAP2', 'MARK1', 'C1orf115', 'MTARC2', 'MTARC1', 'HLX', 'DUSP10', 'HHIPL2', 'TAF1A', 'MIA3', 'AIDA', 'FAM177B', 'DISP1', 'TLR5', 'SUSD4', 'CCDC185', 'CAPN8', 'CAPN2', 'TP53BP2', 'FBXO28', 'DEGS1', 'NVL', 'CNIH4', 'CNIH3', 'DNAH14', 'LBR', 'ENAH', 'SRP9', 'TMEM63A', 'LEFTY1', 'PYCR2', 'LEFTY2', 'SDE2', 'H3-3A', 'ACBD3', 'MIXL1', 'LIN9', 'PARP1', 'STUM', 'ITPKB', 'PSEN2', 'COQ8A', 'CDC42BPA', 'ZNF678', 'SNAP47', 'PRSS38', 'WNT9A', 'WNT3A', 'ARF1', 'C1orf35', 'MRPL55', 'GUK1', 'GJC2', 'IBA57', 'OBSCN', 'TRIM11', 'TRIM17', 'H3-4', 'H2AW', 'H2BU1', 'RNF187', 'RHOU', 'RAB4A', 'CCSAP', 'ACTA1', 'NUP133', 'ABCB10', 'TAF5L', 'URB2', 'GALNT2', 'PGBD5', 'COG2', 'AGT', 'CAPN9', 'C1orf198', 'TTC13', 'ARV1', 'FAM89A', 'TRIM67', 'C1orf131', 'GNPAT', 'EXOC8', 'SPRTN', 'EGLN1', 'TSNAX', 'DISC1', 'SIPA1L2', 'MAP10', 'PCNX2', 'MAP3K21', 'KCNK1', 'SLC35F3', 'COA6', 'TARBP1', 'IRF2BP2', 'TOMM20', 'RBM34', 'ARID4B', 'GNG4', 'LYST', 'NID1', 'GPR137B', 'EDARADD', 'HEATR1', 'ACTN2', 'MTR', 'MT1HL1', 'RYR2', 'ZP4', 'CHRM3', 'FMN2', 'GREM2', 'RGS7', 'FH', 'KMO', 'WDR64', 'EXO1', 'BECN2', 'MAP1LC3C', 'PLD5', 'AKT3', 'ZBTB18', 'C1orf100', 'CATSPERE', 'DESI2', 'COX20', 'HNRNPU', 'EFCAB2', 'KIF26B', 'SMYD3', 'TFB2M', 'CNST', 'SCCPDH', 'AHCTF1', 'ZNF695', 'ZNF670', 'ZNF669', 'ZNF124', 'ZNF496', 'NLRP3', 'OR2B11', 'GCSAML', 'OR2G2', 'OR2G3', 'OR13G1', 'OR6F1', 'OR14A2', 'OR1C1', 'OR14A16', 'OR11L1', 'TRIM58', 'OR2W3', 'OR2T8', 'OR2AJ1', 'OR2L13', 'OR2M5', 'OR2M2', 'OR2M3', 'OR2M4', 'OR2T33', 'OR2T12', 'OR2M7', 'OR14C36', 'OR2T4', 'OR2T6', 'OR2T1', 'OR2T7', 'OR2T2', 'OR2T3', 'OR2T5', 'OR2G6', 'OR2T29', 'OR2T34', 'OR2T10', 'OR2T11', 'OR2T35', 'OR2T27', 'OR14I1', 'LYPD8', 'SH3BP5L', 'ZNF672', 'ZNF692', 'PGBD2']
	names = ['ACKR1']
	data = '/projects/molonc/roth_lab/afratz/projects/HapPy/data/'
	print(len(names))
	cmd = makeCommands(names, data)
	run_commands(*cmd)

if __name__ == '__main__':
	main()
	