import numpy as np
#from scipy.integrate import quad
from matplotlib import pyplot as plt
import os,inspect
import matplotlib.ticker as ticker
from matplotlib import gridspec
import ROOT
import rootpy
import rootpy.plotting as rplot
import rootpy.plotting.root2matplotlib as rplt
from ROOT import TMath as tm

def parameters_kin(sam,var,par):
	topmassup=185
	topmassdown=165
	zmassup=110
	zmassdown=70
	bmassup=5
	bmassdown=3
	umassup=2
	umassdown=0
	wmassup=95
	wmassdown=70

	default_range_up=250
	default_range_down=0
	default_nbins = 64

	mass_nbins =40

	nbin 		= 	default_nbins
	lower_range	=	default_range_down
	upper_range	=	default_range_up
	yup			=	1.2
	ydown		=	1e-4
	ylog		=	1
	up_l 		= 	0.001
	if var=="PT":
		if par == "Z":
			up_l=0.1
			upper_range=250
			lower_range=0
			nbin=20
		if par == "UQuark":
			up_l=0.1
			upper_range=250
			lower_range=0
			nbin=10
		if par == "BQuark":
			upper_range=100
			lower_range=0
			nbin=10
		if par == "WBoson":
			upper_range=250
			lower_range=0
			nbin=10
		if par == "TopQuark":
			upper_range=250
			lower_range=0
			nbin=10
	if var == "Phi":
		lower_range=-4.2
		upper_range=4.2
		nbin = 10
		up_l = 0
	if var == "Eta":
		lower_range=-2.5
		upper_range=2.5
		nbin=10
		up_l = 0

	if(var=="M"):
		nbin = mass_nbins
		up_l = 0
		ylog = 0
		yup = 0
		ydown = 0
		if("W" in par):
			lower_range=wmassdown
			upper_range=wmassup
		if("Z" in par):
			lower_range=zmassdown
			upper_range=zmassup
		if("Top" in par):
			lower_range=topmassdown
			upper_range=topmassup
		if("BQuark" in par):
			lower_range=bmassdown
			upper_range=bmassup
		if("UQuark" in par):
			lower_range=umassdown
			upper_range=umassup
	return [lower_range,upper_range,nbin,up_l,yup,ydown,ylog]

def parameters_RM(p1,p2,va):
	up_l=0
	if va == "R":
		lower_range = 0
		upper_range = 7
		nbin = 64
	if va == "M":
		lower_range = 0
		nbin = 64
		if p1 =="Z":
			lower_range = 0
			upper_range = 1400
			if p2 == "TopQuark":
				lower_range = 150
				upper_range = 800
			if p2 == "UQuark":
				lower_range = 125
				upper_range = 220
				nbin=80
		if p1 =="TopQuark":
			lower_range = 150
			upper_range = 500
			if p2 == "WBoson":
				lower_range = 250
				upper_range = 300
			if p2 == "Z":
				lower_range = 150
				upper_range = 800
	return [lower_range,upper_range,nbin,up_l]

def labelkin(va,par):
	label = r"$"
	if va == "PT":
		label += r"p_T"
	if va == "Eta":
		label += r"\eta"
	if va == "Phi":
		label += r"\Phi"
	if va == "M":
		label += r"m"
	if par == "Z":
		label += r"(Z)$"
	if par == "TopQuark":
		label += r"(t)$"
	if par == "BQuark":
		label += r"(b)$"
	if par == "WBoson":
		label += r"(W)$"
	if par == "UQuark":
		label += r"(u)$"
	if va == "PT" or va == "M":
		label += "/GeV"
	plt.xlabel(label)
	plt.ylabel("number of entries (normalized)",rotation=90)
	ax = plt.gca()
	ax.yaxis.set_label_coords(-0.085,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	plt.grid(alpha=0.5)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

	if va == "M":
		plt.legend(loc="upper right")
	else:
		plt.legend(loc="best")
	if va == "Eta":
		plt.legend(loc="best",fontsize="small")		

def labelkinax(ax,ax2,va,par):
	label = r"$"
	if va == "PT":
		label += r"p_T"
	if va == "Eta":
		label += r"\eta"
	if va == "Phi":
		label += r"\Phi"
	if va == "M":
		label += r"m"
	if par == "Z":
		label += r"(Z)$"
	if par == "TopQuark":
		label += r"(t)$"
	if par == "BQuark":
		label += r"(b)$"
	if par == "WBoson":
		label += r"(W)$"
	if par == "UQuark":
		label += r"(u)$"
	if va == "PT" or va == "M":
		label += "/GeV"
	#ax.set_xlabel(label)
	ax.set_ylabel("number of entries (normalized)",rotation=90)
	ax.yaxis.set_label_coords(-0.085,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax.grid(alpha=0.5)
	#ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax2.set_xlabel(label)
	#ax.yaxis.set_label_coords(-0.085,0.5)
	ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax2.grid(alpha=0.5)
	ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

	if va == "M":
		ax.legend(loc="upper right")
	else:
		ax.legend(loc="best")
	if va == "Eta":
		ax.legend(loc="best",fontsize="small")		


def labelRM(va,p1,p2):
	vari = r"$"
	if va == "R":
		vari += r"\Delta R"
	if va == "M":
		vari += r"m"

	parts = ""

	if p1 == "Z":
		parts += r"Z,"
	if p1 == "TopQuark":
		parts += r"t,"

	if p2 == "Z":
		parts += r"Z"
	if p2 == "TopQuark":
		parts += r"t"
	if p2 == "BQuark":
		parts += r"b"
	if p2 == "WBoson":
		parts += r"W"
	if p2 == "UQuark":
		parts += r"u"

	einh = ""
	if va == "M":
		einh += "/GeV"

	plt.xlabel(vari+"("+parts+")"+r"$"+einh)
	plt.ylabel("number of entries (normalized)",rotation=90)
	ax = plt.gca()
	ax.yaxis.set_label_coords(-0.085,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	plt.grid(alpha=0.5)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.legend(loc="best")
	if p1 == "Z" and p2=="UQuark" and va=="M":
		plt.yscale("log")
		plt.ylim(0.0005,0.37)

def labelRMax(ax,ax2,va,p1,p2):
	vari = r"$"
	if va == "R":
		vari += r"\Delta R"
	if va == "M":
		vari += r"m"

	parts = ""

	if p1 == "Z":
		parts += r"Z,"
	if p1 == "TopQuark":
		parts += r"t,"

	if p2 == "Z":
		parts += r"Z"
	if p2 == "TopQuark":
		parts += r"t"
	if p2 == "BQuark":
		parts += r"b"
	if p2 == "WBoson":
		parts += r"W"
	if p2 == "UQuark":
		parts += r"u"

	einh = ""
	if va == "M":
		einh += "/GeV"

	#ax.set_xlabel(vari+"("+parts+")"+r"$"+einh)
	ax.set_ylabel("number of entries (normalized)",rotation=90)
	ax.yaxis.set_label_coords(-0.085,0.5)
	ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax.grid(alpha=0.5)
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax.legend(loc="best")
	ax2.set_xlabel(vari+"("+parts+")"+r"$"+einh)
	#ax.yaxis.set_label_coords(-0.085,0.5)
	ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
	ax2.grid(alpha=0.5)
	ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	if p1 == "Z" and p2=="UQuark" and va=="M":
		ax.set_yscale("log")
		ax.set_ylim(0.0005,0.37)

def set_dyn_binning(va,lo,up,n,err=0.05):
	#eps=0.0001
	#va = va[(np.abs(va)>up_l) & (va!=999.9)]
	#va = np.clip(va,lo+eps,up-eps)
	binning = np.linspace(lo,up,n)
	h1 = rplot.Hist(binning)
	map(h1.Fill,va)
	h1.Scale(1/(h1.Integral(0,h1.GetNbinsX()+1)))
	#print(binning)

	for i in range(n-1,1,-1):
		#print(i," ",binning[i])
		#if (h1.GetBinContent(i)==0):
			#print(0)
		#	continue
		#print(h1.GetBinError(i)/h1.GetBinContent(i))

		if (h1.GetBinContent(i)==0):
			binning = np.delete(binning, i-1)
			h1 = rplot.Hist(binning)
			map(h1.Fill,va)
			h1.Scale(1/(h1.Integral(0,h1.GetNbinsX()+1)))
			continue

		if (h1.GetBinError(i)/h1.GetBinContent(i)) > err:
			binning = np.delete(binning, i-1)
			h1 = rplot.Hist(binning)
			map(h1.Fill,va)
			h1.Scale(1/(h1.Integral(0,h1.GetNbinsX()+1)))
			continue

	return binning


def KMS(a,b):
	na = len(a)
	nb = len(b)
	a = np.sort(a)
	b = np.sort(b)
	h = tm.KolmogorovTest(na,a,nb,b,"")
	return h