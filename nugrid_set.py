''' Tool for analysis of mesa and mppnp data of paper2 calculations
    v0.1, 15JAN2013: Christian Ritter,Falk Herwig
    v0.2, 01APRIL2013: Christian Ritter,Falk Herwig,Raphael Hirschi,Marco Pignatari + Collaboration
    v0.3, July2016: CR: added to nugridpy

	##########################
	#### Many functions still under development  ####
	##########################



	This tools contains two classes which can deal with set data
	of multsple mesa runs and the class mppnp_set the analysis
	of multiple MPPNP dirs. Check out the documentation in
	the DOC dir.
	

	Both classes allow during their initialisation the set of 
	one dir path - containing the set calculations or the set of 
	multiple run dirs.


	Python 2.7.3 (default, Jul 24 2012, 10:05:38) 
	Type "copyright", "credits" or "license" for more information.

	IPython 0.12 -- An enhanced Interactive Python.
	?         -> Introduction and overview of IPython's features.
	%quickref -> Quick reference.
	help      -> Python's own help system.
	object?   -> Details about 'object', use 'object??' for extra details.

	IPython profile: numpy

	Welcome to pylab, a matplotlib-based Python environment [backend: GTKAgg].
	For more information, type 'help(pylab)'.



	################################################################
	#for working with mesa set - e.g. from inside the set directory

	In [1]: import nugrid_set as set
	In [2]: setresults=set.mesa_set(".")
	In [3]: setresults.multi_DUP()

	e.g for calculation of DUP parameter evolution of all set runs.



	################################################################
	#for working with mppnp set - e.g. from inside the set directory
 

        In [2]: setresults=set.mppnp_set(".")
        In [3]: setresults.

	#for calculation of yields of the set 

'''
##need to be changed, matplotlib import
import matplotlib.pyplot as plt
from matplotlib.pyplot import * 
import numpy as np
import os
from nugridse import *
from mesa import *
#from paper2_analysis import *
import glob
from mpl_toolkits.mplot3d import Axes3D
import utils


symbs=utils.symbol_list('lines1')



class mppnp_set(se):

	'''
		This class allows access to multiple runs in the path dir.
		Instead of defining one path, an array of paths , mult_dir can be
		defined. In this case each path in multi_dir represents a run.
		

                rundir - dir with run directory
                multi_dir - specify dirs at different location
                extra_label - set labels of runs, else run dir name will be chosen      

                e.g.:
                        
                        or
                                

                Using the set_plot* functions you can instantly plot diagrams
                for all runs defined in MPPNP set


		Functions:


		Abundance extraction:

		def set_pocket: plots composition at certain cycle and mass cell.
				uses def plot_abu_atmasscoord

		def set_surface_abundances: basically the same as def set_pocket
						but for surface

		def set_surface_plots: Plot surface abundance evolution of multiple runs
					needs to be checked, plots [hs/ls', [Rb/Sr]..

		def plot_abu_atmasscoord: plots isotope distribution at cycle,masscell, exists already above...

		def plot_surf_element: plots surface evolution of certain element

		def plot_surf_isotope: plots surface evolution of certain isotope

		def def set_plot_prof: abundance profile of certain cycle for different isotopes

		def surface_plot: Plots [hs/Fe] vs [ls/Fe] and [Pb/hs] vs [ls/Fe] , used in thesis?!

		def pocket_composition: uses read_iso_abund_marco, useful, but already there?

		def get_elem_abu: calculates elemental massfraction from isotopic massfractions

		def final_bottom_envelope_set1

		Calculates remnant mass for set1

		def final_bottom_envelope
		Calculates remnant mass for Set1 extension


		Yields related:
			
		    To write element tables: (example)	
			write_gce_input(file_name=filename,elem=True,exp_dir='/nfs/rpod3/critter/Results/DATA/set1.2/ppd_exp',delay=True)
		    To write isotopic tables: (example)
			write_gce_input(file_name=filename,elem=False,GCE_input=False,GCE_input_starkenburg=True,exp_dir='/nfs/rpod3/critter/Results/DATA/set1.2/ppd_exp',delay=True)

			To add certain parameter to yield tables:
			write_gce_input_lifetimes(file_name=filename) : add lifetimes
			write_gce_input_exp_energy(file_name=filename) :  add kinetic energies
		

		######Calculate production factors:		

		#for data tables of elements or isotopes:

			#stellar winds + pre-exp +  exp (delay and rapid)?
			def write_prod_fact_stellarwinds:

		#for figures (prodfac vs. ini mass) of elements or isotopes:

			#for 1 set only:
				#stellar winds
				prod_fac_massgrid_single_wind	
				#pre-exp +exp
				prod_fac_massgrid
			#To compare more than 1 set:
				prod_fac_massgrid

		#######


		def prod_fac_massgrid: Plots behaviour star mass dependent yields of different isotopes

		def yields_massgrid: not in use, what is the meaning

		def write_prod_fact_stellarwinds:       calculates yields or production factors, and input for GCE functios
							uses self.weighted_yields
							uses de write_latex_table_set1:

		def set1_extension_weight: calculates the IMF weights with weights of paper1 with mass extension

		def set1_weight_marco: calculates IMF as done in paper1

		def weighted_yields_pre_explosion: Marcos routine adapted, calc pre SN yields

		def def prodfac_explosion: Marcos routine adapted, calc SN yields

		def weighted_yields_explosion: really weighted???

		def weighted_yields: important function, used, returns prod_facs,isotopes,yields,iniabus,remn_mass


		Other:

		def plot_saga_data: reads saga platform plots
					uses:forum.astro.keele.ac.uk/utils/spectroscopy_data/ saga_platform

		def set_triple_isotope_plot: grain data plot, uses def plot_isoratios
						may be adaption necessary?

		def def write_GCE_input_fxt2: Reads input_file of input for GCE code of fxt and replaces
						uses def write_prod_fact_stellarwinds
		
		def write_GCE_input_fxt: Reads input_file of input for GCE code of fxt and replaces
		
		def write_GCE_input_fxt_OLD

		def write_latex_table_set1: writes latex table as in paper1

		def write_ascii_table: Writes ascii table - needs to be replaced by paper1 function
	


	'''

	def __init__(self,rundir='.',multi_dir=[],extra_label=[],exp_type=""):



		self.explosion_type=exp_type    #delay or rapid

                if len(multi_dir)==0:
			slist = os.listdir(rundir)
			if len(exp_type)>0:
				expr = re.compile(exp_type)
				slist=(filter(expr.search,slist))
			
		else:
                        slist=multi_dir
		pwd=os.getcwd()
		pattern='.h5'
		expr = re.compile(pattern)
		expr = re.compile(pattern)
		self.runs_H5_out=[]
		self.runs_H5_surf=[]
		self.runs_H5_restart=[]
		self.run_dirs_name=[]
		self.extra_label=[]
                runs_H5_out=[]
                runs_H5_surf=[]
		runs_H5_restart=[]
                run_dirs_name=[]
		extra_label_1=[]
		i=0
		for element in slist:
                        
			if len(multi_dir)==0:
				if rundir[0]=="/":
					run_path=rundir+"/"+element
				else:
                                	run_path=pwd+"/"+rundir+"/"+element
                        else:
                                if multi_dir[i][0] == "/":
                                        run_path=multi_dir[i]
                                else:
					if len(rundir)>0:
						run_path=rundir+"/"+multi_dir[i]
					else:
                                        	run_path=pwd+"/"+multi_dir[i]

			if os.path.isdir(run_path+"/H5_out") and os.path.isdir(run_path+"/H5_surf"):
				sefiles = os.listdir(run_path+"/H5_out")
				if (filter(expr.search,sefiles)) <1:
					print "Warning: No hdf5 out files found in "+run_path+"/H5_out"
				else:
					runs_H5_out.append(run_path+"/H5_out")	
				sefiles = os.listdir(run_path+"/H5_surf")
				if (filter(expr.search,sefiles)) <1:
					print "Warning: No hdf5 surf files found in "+run_path+"/H5_surf"
				else:
					runs_H5_surf.append(run_path+"/H5_surf")
				sefiles = os.listdir(run_path+"/H5_restart")
				if (filter(expr.search,sefiles)) <1:
                                        print "Warning: No hdf5 restart files found in "+run_path+"/H5_surf"
                                else:
					runs_H5_restart.append(run_path+"/H5_restart")
				if len(multi_dir)>0:
					element_1=element.split('/')[-1]	
						
				if len(extra_label)>0:
					extra_label_1.append(extra_label[i])
				else:
					extra_label_1.append(element)
					element_1=element
				run_dirs_name.append(element_1)
				print "Read "+run_path		
			i+=1
		###order lists###
		to_sort_names=[]
		for i in range(len(run_dirs_name)):
			test=float(run_dirs_name[i][1:4])*100 
			to_sort_names.append(test)
			print test
		
		sorted_indices=sorted(range(len(to_sort_names)),key=lambda k: to_sort_names[k])	
		print sorted_indices
		print run_dirs_name
		k=0
		for i in sorted_indices:
			print i
			self.run_dirs_name.append(run_dirs_name[i])
			self.extra_label.append(extra_label_1[i])
			self.runs_H5_surf.append(runs_H5_surf[i])
			self.runs_H5_out.append(runs_H5_out[i])
			self.runs_H5_restart.append(runs_H5_restart[i])
			k+=1
		print self.run_dirs_name


        def initial_finall_mass_relation(self,marker='o',linestyle='--'):
		'''
			INtiial to final mass relation
		'''

		final_m=[]
                ini_m=[]
                for i in range(len(self.runs_H5_surf)):
			sefiles=se(self.runs_H5_out[i])
                        ini_m.append(sefiles.get("mini")) 
			h1=sefiles.get(int(sefiles.se.cycles[-2]),'H-1')
			mass=sefiles.get(int(sefiles.se.cycles[-2]),'mass')
			idx=-1
			for k in range(len(h1)):
				if h1[k]>0.1:
					idx=k
					break
			final_m.append(mass[idx])
		label='Z='+str(sefiles.get('zini'))
		plt.plot(ini_m,final_m,label=label,marker=marker,linestyle=linestyle)
		plt.xlabel('$M_{Initial} [M_{\odot}]$',size=23)
		plt.ylabel('$M_{Final} [M_{\odot}]$',size=23)


	def final_bottom_envelope_set1(self):
		'''
			For paper1 marco routine:
			 Numbers of remnant mass shell masses, exists also in mesa_set!
		'''
		inim=[]
		remnm=[]
		for i in range(len(self.runs_H5_surf)):
			m1p65_last=se(self.runs_H5_out[i])	
			mass_dummy=m1p65_last.se.get(m1p65_last.se.cycles[len(m1p65_last.se.cycles)-1],'mass')
			top_of_envelope=mass_dummy[len(mass_dummy)-1]
			h_dummy=m1p65_last.se.get(m1p65_last.se.cycles[len(m1p65_last.se.cycles)-1],'iso_massf','H-1')

			for j in range(len(mass_dummy)):
        			if h_dummy[j] > 0.05:
                			bottom_of_envelope = mass_dummy[j]
                			break
			inim.append(m1p65_last.get("mini"))
			remnm.append(bottom_of_envelope)
		print "M_initial | M_remn/bottom of envelope"
		for i in range(len(inim)):
			print inim[i],"|",remnm[i]	

        def remnant_lifetime_agb(self):
                '''
                        For paper1 extension:
			bottom_envelope
                         Numbers of remnant mass shell masses, exists also in mesa_set + star age!
                '''
                inim=[]
                remnm=[]
		time11=[]
		tottime=[]
		c_core=[]
		o_core=[]
		small_co_core=[]
		c_core_center=[]	
                for i in range(len(self.runs_H5_surf)):
                        m1p65_last=se(self.runs_H5_out[i])
                        mass_dummy=m1p65_last.se.get(m1p65_last.se.cycles[len(m1p65_last.se.cycles)-2],'mass')
                        top_of_envelope=mass_dummy[len(mass_dummy)-1]
                        h_dummy=m1p65_last.se.get(m1p65_last.se.cycles[len(m1p65_last.se.cycles)-2],'iso_massf','H-1')
			c_dummy=m1p65_last.se.get(m1p65_last.se.cycles[len(m1p65_last.se.cycles)-2],'iso_massf','C-12')
                        o_dummy=m1p65_last.se.get(m1p65_last.se.cycles[len(m1p65_last.se.cycles)-2],'iso_massf','O-16')
			for j in range(len(mass_dummy)):
                                if h_dummy[j] > 1e-1:
                                        bottom_of_envelope = mass_dummy[j]
                                        break
                        inim.append(m1p65_last.get("mini"))
                        remnm.append(bottom_of_envelope)

			###Calculate the lifetime (MS)
                        sefiles=m1p65_last
                        cycs=[]
                        for k in range(5,len(sefiles.se.cycles),5):
                                cycs.append(int(sefiles.se.cycles[k]))
                        w=0
                        for cyc in cycs:
                                c12_center=sefiles.get(cyc,'C-12')[0]
                                #c12_center=c12[w][0]
                                w+=1
                                if c12_center>1e-1:
                                        time1=(sefiles.get(cyc,'age')*sefiles.get('age_unit'))/31557600.
                                        time11.append(time1)
                                        break
			tottime.append(sefiles.get(int(sefiles.se.cycles[-1]),'age')/31557600.)



                print "M_initial | M_remn/bottom of envelope | total lifetime"
                for i in range(len(inim)):
                        print inim[i],"|",'{:.3E}'.format(remnm[i]),"|",'{:.3E}'.format(tottime[i])


	def plot_surf_composition(self,hsfe_hsls = True,lsfe_hsls = False,rbsr_hsls = False,rbzr_hsls = False,sry_srzr = False,rbfe_sfe = False,iso_ratio= False,isotopes_x = ['Mg-25','Mg-24'],isotopes_y = ['Mg-26','Mg-24'],sparsity = 5000,marker1=['o'],symbols1=['k-'],color=['k']):
		'''
			Set only one of the vriables True to plot quantitie4s.
			In case you set iso_ratio = True:
				isotooes_x and isototopes_y are used.  
		'''

		import nugridse as mp
		import matplotlib.pyplot as pl
		import sys

		plotting=True

		### SET THIS IF PLOTTING = TRUE ####################
		###>>>>!!!! ONLY ONE OPTION CAN BE TRUE !!!!<<<< ###
		#test=[hsfe_hsls,lsfe_hsls,
		#hsfe_hsls = False
		#lsfe_hsls = False
		#rbsr_hsls = False
		#rbzr_hsls = False
		#sry_srzr = False
		#rbfe_sfe = False
		#iso_ratio= True
		###################################

                for i in range(len(self.runs_H5_surf)):
                       	m3z2m2=se(self.runs_H5_surf[i])
			run_mass=[m3z2m2]

			markers=[marker1[i]]#['ks']#,'r^','gh','bd','mv','cH','kD']
			symbol=[symbols1[i]]#['k-']#,'r--','g-','b:','m-','c--','k-']
			color1=color[i]
			ls_element = ['Sr-84','Sr-86','Sr-87','Sr-88','Y-89','Zr-90','Zr-91','Zr-92','Zr-94','Zr-96']
			hs_element = ['Ba-130','Ba-132','Ba-134','Ba-135','Ba-136','Ba-137','Ba-138','La-138','La-139','Nd-142','Nd-143','Nd-144','Nd-145','Nd-146','Nd-148','Nd-150','Sm-144','Sm-147','Sm-148','Sm-149','Sm-150','Sm-152','Sm-154']
			# notice that Pb (3rd s-process index) is not included here.
			s_element =  ['Sr-84','Sr-86','Sr-87','Sr-88','Y-89','Zr-90','Zr-91','Zr-92','Zr-94','Zr-96','Ba-130','Ba-132','Ba-134','Ba-135','Ba-136','Ba-137','Ba-138','La-138','La-139','Nd-142','Nd-143','Nd-144','Nd-145','Nd-146','Nd-148','Nd-150','Sm-144','Sm-147','Sm-148','Sm-149','Sm-150','Sm-152','Sm-154']	

			ls_element = ['Sr-84','Sr-86','Sr-87','Sr-88','Y-89','Zr-90','Zr-91','Zr-92','Zr-94','Zr-96']
			hs_element = ['Ba-130','Ba-132','Ba-134','Ba-135','Ba-136','Ba-137','Ba-138','La-138','La-139','Nd-142','Nd-143','Nd-144','Nd-145','Nd-146','Nd-148','Nd-150','Sm-144','Sm-147','Sm-148','Sm-149','Sm-150','Sm-152','Sm-154']
			# notice that Pb (3rd s-process index) is not included here.
			s_element =  ['Sr-84','Sr-86','Sr-87','Sr-88','Y-89','Zr-90','Zr-91','Zr-92','Zr-94','Zr-96','Ba-130','Ba-132','Ba-134','Ba-135','Ba-136','Ba-137','Ba-138','La-138','La-139','Nd-142','Nd-143','Nd-144','Nd-145','Nd-146','Nd-148','Nd-150','Sm-144','Sm-147','Sm-148','Sm-149','Sm-150','Sm-152','Sm-154']	

			s_ini = []
			ls_ini = []
			hs_ini = []
			pb_ini = []
			fe_ini = []
			rb_ini = []
			sr_ini = []
			zr_ini = []
			y_ini = []
			zra_ini = []
			zrb_ini = []
			gda_ini = []
			gdb_ini = []

			if plotting:

				if iso_ratio:

					isotopes_x = ['Mg-25','Mg-24']
					isotopes_y = ['Mg-26','Mg-24']

					#isotopes_x = ['Zr-96','Zr-94']
					#isotopes_y = ['Gd-152','Gd-154']
					#isotopes_y = ['Zr-90','Zr-94']
					#isotopes_y = ['Zr-91','Zr-94']
					#isotopes_y = ['Zr-92','Zr-94']

					#isotopes_x = ['Ba-135','Ba-136']
					#isotopes_y = ['Ba-138','Ba-136']
					#isotopes_y = ['Ba-134','Ba-136']
					#isotopes_y = ['Ba-137','Ba-136']

					#isotopes_x = ['Sr-84','Sr-86']
					#isotopes_y = ['Sr-87','Sr-86']
					#isotopes_y = ['Sr-88','Sr-86']

					# get initial value
					initial_isotopic_ratio_x = float(run_mass[0].se.get(min(run_mass[0].se.cycles),'iso_massf',isotopes_x[0]))/float(run_mass[0].se.get(min(run_mass[0].se.cycles),'iso_massf',isotopes_x[1]))
					initial_isotopic_ratio_y = float(run_mass[0].se.get(min(run_mass[0].se.cycles),'iso_massf',isotopes_y[0]))/float(run_mass[0].se.get(min(run_mass[0].se.cycles),'iso_massf',isotopes_y[1]))

					# sparcity for cycles I am looking at.
					#sparsity = 5000

					I_want_iso_ratio = False
					I_want_delta = True

					isotopic_ratio_x  = []
					isotopic_ratio_y = []
					isotopic_ratio_x_tps = []
					isotopic_ratio_y_tps = []
					isotopic_ratio_x_tps_co = []
					isotopic_ratio_y_tps_co = []

					k = 0

					for i in run_mass:
						dum_isotopic_ratio_x = []
						dum_isotopic_ratio_y = []
						dumdum_isotopic_ratio_x = []
						dumdum_isotopic_ratio_y = []
						dum_isotopic_ratio_x_tps = []
						dum_isotopic_ratio_y_tps = []
						dum_isotopic_ratio_x_tps_co = []
						dum_isotopic_ratio_y_tps_co = []
						co_ratio=[]
						for j in i.se.cycles[0::sparsity]:
							print 'j= ', j
							dum_isotopic_ratio_x.append(float(i.se.get(j,'iso_massf',isotopes_x[0]))/float(i.se.get(j,'iso_massf',isotopes_x[1])))
							dum_isotopic_ratio_y.append(float(i.se.get(j,'iso_massf',isotopes_y[0]))/float(i.se.get(j,'iso_massf',isotopes_y[1])))
							co_ratio.append((float(i.se.get(j,'iso_massf','C-12')*4.))/(float(i.se.get(j,'iso_massf','O-16'))*3))

							if (len(co_ratio)>1):
								if (co_ratio[len(co_ratio)-1]>(co_ratio[len(co_ratio)-2]+0.02)):
									dum_isotopic_ratio_x_tps.append(float(i.se.get(j,'iso_massf',isotopes_x[0]))/float(i.se.get(j,'iso_massf',isotopes_x[1])))
									dum_isotopic_ratio_y_tps.append(float(i.se.get(j,'iso_massf',isotopes_y[0]))/float(i.se.get(j,'iso_massf',isotopes_y[1])))

									if (co_ratio[len(co_ratio)-1]>1.):
										dum_isotopic_ratio_x_tps_co.append(float(i.se.get(j,'iso_massf',isotopes_x[0]))/float(i.se.get(j,'iso_massf',isotopes_x[1])))
										dum_isotopic_ratio_y_tps_co.append(float(i.se.get(j,'iso_massf',isotopes_y[0]))/float(i.se.get(j,'iso_massf',isotopes_y[1])))

							if I_want_delta:
								dumdum_isotopic_ratio_x=(np.array(dum_isotopic_ratio_x)/initial_isotopic_ratio_x-1.)*1000.
								dumdum_isotopic_ratio_y=(np.array(dum_isotopic_ratio_y)/initial_isotopic_ratio_y-1.)*1000.
								dumdum_isotopic_ratio_x_tps=(np.array(dum_isotopic_ratio_x_tps)/initial_isotopic_ratio_x-1.)*1000.
								dumdum_isotopic_ratio_y_tps=(np.array(dum_isotopic_ratio_y_tps)/initial_isotopic_ratio_y-1.)*1000.
								dumdum_isotopic_ratio_x_tps_co=(np.array(dum_isotopic_ratio_x_tps_co)/initial_isotopic_ratio_x-1.)*1000.
								dumdum_isotopic_ratio_y_tps_co=(np.array(dum_isotopic_ratio_y_tps_co)/initial_isotopic_ratio_y-1.)*1000.

						if I_want_iso_ratio:
							isotopic_ratio_x.append(dum_isotopic_ratio_x)
							isotopic_ratio_y.append(dum_isotopic_ratio_y)
						else:
							isotopic_ratio_x.append(dumdum_isotopic_ratio_x)
							isotopic_ratio_y.append(dumdum_isotopic_ratio_y)
							isotopic_ratio_x_tps.append(dumdum_isotopic_ratio_x_tps)
							isotopic_ratio_y_tps.append(dumdum_isotopic_ratio_y_tps)
							isotopic_ratio_x_tps_co.append(dumdum_isotopic_ratio_x_tps_co)
							isotopic_ratio_y_tps_co.append(dumdum_isotopic_ratio_y_tps_co)

						k = k+1

					markersss=['ko','b^','rh','gd','mv','cH','kD']
					symbol=['k-','b--','r-.','g:','m-','c--','k-']

					mass_label =[]
					for i in run_mass:
						mini=float(i.se.get('mini'))
						zini=float(i.se.get('zini'))
						label=str(mini)+'$M_{\odot}$, Z='+str(zini)
						mass_label.append(label)

					metallicity_label=['M3z2m2','M3z2m2_nomol','M3z1m2','3 Msun, Z=0.02, PI13']
					metallicity_label=[label]

					array_to_plot_x = isotopic_ratio_x
					array_to_plot_y = isotopic_ratio_y
					array_to_plot_x_tps = isotopic_ratio_x_tps
					array_to_plot_y_tps = isotopic_ratio_y_tps
					array_to_plot_x_tps_co = isotopic_ratio_x_tps_co
					array_to_plot_y_tps_co = isotopic_ratio_y_tps_co

					if I_want_delta:
						pl.axhline(y=0.,linewidth=2, color='k')
						pl.axvline(x=0.,linewidth=2, color='k')
					else:
						#	pl.axhline(y=initial_isotopic_ratio_y/initial_isotopic_ratio_y,linewidth=2, color='k')
						pl.axhline(y=initial_isotopic_ratio_y,linewidth=2, color='k')
						#	pl.axvline(x=initial_isotopic_ratio_x/initial_isotopic_ratio_x,linewidth=2, color='k')
						pl.axvline(x=initial_isotopic_ratio_x,linewidth=2, color='k')

					for k in range(0,len(array_to_plot_x)):
						if I_want_delta:
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=4.,label=str(metallicity_label[k]))
							pl.plot(array_to_plot_x_tps[k],array_to_plot_y_tps[k],markersss[k],markersize=10.,linewidth=4.)
							pl.plot(array_to_plot_x_tps_co[k],array_to_plot_y_tps_co[k],markersss[k],markersize=20.,linewidth=4.)

						else:
							#		pl.loglog(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=2.,label=str(mass_label[k])+'Msun ')
							pl.loglog(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=2.,label=str(mass_label[k])+'Msun ')

					pl.legend(numpoints=1,loc=2,prop={'size':20})
					size=20
					pl.xticks(size=size)
					pl.yticks(size=size)
					#pl.axis([0.08,0.5,0.1,1.])
					#pl.grid()

					if I_want_delta:
						#pl.xlabel(r'$\delta$($^{25}$Mg/$^{24}$Mg)', fontsize=30)
						#pl.ylabel(r'$\delta$($^{26}$Mg/$^{24}$Mg)', fontsize=30)

						pl.xlabel(r'$\delta$($^{96}$Zr/$^{94}$Zr)', fontsize=20)
						pl.ylabel(r'$\delta$($^{152}$Gd/$^{154}$Gd)', fontsize=20)
						#pl.ylabel(r'$\delta$($^{90}$Zr/$^{94}$Zr)', fontsize=30)
						#pl.ylabel(r'$\delta$($^{91}$Zr/$^{94}$Zr)', fontsize=30)
						#pl.ylabel(r'$\delta$($^{92}$Zr/$^{94}$Zr)', fontsize=30)
					
						#pl.xlabel('d($^{135}$Ba/$^{136}$Ba)', fontsize=20)
						#pl.ylabel('d($^{138}$Ba/$^{136}$Ba)', fontsize=20)
						#pl.ylabel('d($^{134}$Ba/$^{136}$Ba)', fontsize=20)
						#pl.ylabel('d($^{137}$Ba/$^{136}$Ba)', fontsize=20)

						#pl.xlabel(r'$\delta$($^{84}$Sr/$^{86}$Sr)', fontsize=20)
						#pl.ylabel(r'$\delta$($^{87}$Sr/$^{86}$Sr)', fontsize=20)
						#pl.ylabel(r'$\delta$($^{88}$Sr/$^{86}$Sr)', fontsize=20)
					else:
						pl.xlabel('$^{25}$Mg/$^{24}$Mg', fontsize=20)
						pl.ylabel('$^{26}$Mg/$^{24}$Mg', fontsize=20)

						#pl.xlabel('$^{96}$Zr/$^{94}$Zr', fontsize=30)
						#pl.ylabel('$^{152}$Gd/$^{154}$Gd', fontsize=30)
						#pl.ylabel('$^{90}$Zr/$^{94}$Zr', fontsize=30)
						#pl.ylabel('$^{91}$Zr/$^{94}$Zr', fontsize=30)
						#pl.ylabel('$^{92}$Zr/$^{94}$Zr', fontsize=30)
					
						#pl.xlabel('d($^{135}$Ba/$^{136}$Ba)', fontsize=20)
						#pl.ylabel('d($^{138}$Ba/$^{136}$Ba)', fontsize=20)
						#pl.ylabel('d($^{134}$Ba/$^{136}$Ba)', fontsize=20)
						#pl.ylabel('d($^{137}$Ba/$^{136}$Ba)', fontsize=20)

						#pl.xlabel('$^{84}$Sr/$^{86}$Sr', fontsize=20)
						#pl.ylabel('$^{87}$Sr/$^{86}$Sr', fontsize=20)
						#pl.ylabel('$^{88}$Sr/$^{86}$Sr', fontsize=20)
						
					sys.exit()             

			for i in run_mass:
				dum_s_ini = 0.
				dum_ls_ini = 0.
				dum_hs_ini = 0.
				dum_pb_ini = 0.
				dum_fe_ini = 0.
				dum_rb_ini = 0.
				dum_sr_ini = 0.
				dum_zr_ini = 0.	
				dum_pb_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Pb-204'))+float(i.se.get(min(i.se.cycles),'iso_massf','Pb-206'))+float(i.se.get(min(i.se.cycles),'iso_massf','Pb-207'))+float(i.se.get(min(i.se.cycles),'iso_massf','Pb-208'))
				#
				dum_fe_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Fe-54'))+float(i.se.get(min(i.se.cycles),'iso_massf','Fe-56'))+float(i.se.get(min(i.se.cycles),'iso_massf','Fe-57'))+float(i.se.get(min(i.se.cycles),'iso_massf','Fe-58'))	
				dum_rb_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Rb-85'))+float(i.se.get(min(i.se.cycles),'iso_massf','Rb-87'))
				dum_sr_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Sr-84'))+float(i.se.get(min(i.se.cycles),'iso_massf','Sr-86'))+float(i.se.get(min(i.se.cycles),'iso_massf','Sr-87'))+float(i.se.get(min(i.se.cycles),'iso_massf','Sr-88'))
				dum_zr_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Zr-90'))+float(i.se.get(min(i.se.cycles),'iso_massf','Zr-91'))+float(i.se.get(min(i.se.cycles),'iso_massf','Zr-92'))+float(i.se.get(min(i.se.cycles),'iso_massf','Zr-94'))+float(i.se.get(min(i.se.cycles),'iso_massf','Zr-96'))
				dum_y_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Y-89'))
				dum_zra_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Zr-96'))
				dum_gda_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Gd-152'))
				dum_zrb_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Zr-94'))
				dum_gdb_ini = float(i.se.get(min(i.se.cycles),'iso_massf','Gd-154'))


				for j in ls_element:
					dum_ls_ini = dum_ls_ini + float(i.se.get(min(i.se.cycles),'iso_massf',j))
				for j in hs_element:
					dum_hs_ini = dum_hs_ini + float(i.se.get(min(i.se.cycles),'iso_massf',j))
				for j in s_element:
					dum_s_ini = dum_s_ini + float(i.se.get(min(i.se.cycles),'iso_massf',j))
				fe_ini.append(dum_fe_ini)
				pb_ini.append(dum_pb_ini)
				ls_ini.append(dum_ls_ini)    #  /float(len(ls_element)))
				hs_ini.append(dum_hs_ini)     #/float(len(hs_element)))
				s_ini.append(dum_s_ini)  #/float(len(s_element)))				
				rb_ini.append(dum_rb_ini)
				sr_ini.append(dum_sr_ini)
				zr_ini.append(dum_zr_ini)
				y_ini.append(dum_y_ini)
				zra_ini.append(dum_zra_ini)
				gda_ini.append(dum_gda_ini)
				zrb_ini.append(dum_zrb_ini)
				gdb_ini.append(dum_gdb_ini)



			# sparcity for cycles I am looking at.
			sparsity = 5000
			s_fe  = []
			ls_fe = []
			hs_fe = []
			pb_fe = []
			hs_ls = []
			rb_sr = []
			rb_zr = []
			rb_fe = []
			sr_y = []
			sr_zr = []
			sr_zr_tps=[]
			sr_y_tps=[]
			sr_zr_tps_co=[]
			sr_y_tps_co=[]
			s_fe_tps  = []
			ls_fe_tps = []
			hs_fe_tps = []
			pb_fe_tps = []
			hs_ls_tps = []
			rb_sr_tps = []
			rb_zr_tps = []
			rb_fe_tps = []
			s_fe_tps_co  = []
			ls_fe_tps_co = []
			hs_fe_tps_co = []
			pb_fe_tps_co = []
			hs_ls_tps_co = []
			rb_sr_tps_co = []
			rb_zr_tps_co = []
			rb_fe_tps_co = []
			gd=[]
			zr=[]
			gd_tps=[]
			zr_tps=[]
			gd_tps_co=[]
			zr_tps_co=[]
			co=[]
			k = 0

			for i in run_mass:

				jjjj=0
				dum_s_fe = []
				dum_ls_fe = []
				dum_hs_fe = []
				dum_pb_fe = []
				dum_hs_ls = []
				dum_rb_sr = []
				dum_rb_zr = []
				dum_rb_fe = []
				dum_sr_y = []
				dum_sr_zr = []
				dum_sr_zr_tps=[]
				dum_sr_y_tps=[]
				dum_sr_zr_tps_co=[]
				dum_sr_y_tps_co=[]
				dum_s_fe_tps = []
				dum_ls_fe_tps = []
				dum_hs_fe_tps = []
				dum_pb_fe_tps = []
				dum_hs_ls_tps = []
				dum_rb_sr_tps = []
				dum_rb_fe_tps = []
				dum_rb_zr_tps = []
				dum_sr_y_tps = []
				dum_sr_zr_tps = []
				dum_s_fe_tps_co = []
				dum_ls_fe_tps_co = []
				dum_hs_fe_tps_co = []
				dum_pb_fe_tps_co = []
				dum_hs_ls_tps_co = []
				dum_rb_sr_tps_co = []
				dum_rb_zr_tps_co = []
				dum_rb_fe_tps_co = []
				dum_rb_sr_tps_co = []
				dum_sr_y_tps_co = []
				dum_sr_zr_tps_co = []
				dum_gd=[]
				dum_zrr=[]
				dum_gd_tps=[]
				dum_zr_tps=[]
				dum_gd_tps_co=[]
				dum_zr_tps_co=[]
				dum_co=[]

				for j in i.se.cycles[0::sparsity]:
					dum_s = 0.
					dum_ls = 0.
					dum_hs = 0.
					dum_fe = 0.
					dum_pb = 0.
					dum_rb = 0.
					dum_sr = 0.
					dum_zr = 0.
					dum_y = 0.
					dum_c = 0.
					dum_o = 0.
					dum_gda = 0.
					dum_gdb = 0.
					dum_zra = 0.
					dum_zrb = 0.

					print j
					dum_fe = float(i.se.get(j,'iso_massf','Fe-54'))+float(i.se.get(j,'iso_massf','Fe-56'))+float(i.se.get(j,'iso_massf','Fe-57'))+float(i.se.get(j,'iso_massf','Fe-58')) 
					dum_pb = float(i.se.get(j,'iso_massf','Pb-204'))+float(i.se.get(j,'iso_massf','Pb-206'))+float(i.se.get(j,'iso_massf','Pb-207'))+float(i.se.get(j,'iso_massf','Pb-208'))
					dum_rb = float(i.se.get(j,'iso_massf','Rb-85'))+float(i.se.get(j,'iso_massf','Rb-87'))
					dum_sr = float(i.se.get(j,'iso_massf','Sr-84'))+float(i.se.get(j,'iso_massf','Sr-86'))+float(i.se.get(j,'iso_massf','Sr-87'))+float(i.se.get(j,'iso_massf','Sr-88'))
					dum_zr = float(i.se.get(j,'iso_massf','Zr-90'))+float(i.se.get(j,'iso_massf','Zr-91'))+float(i.se.get(j,'iso_massf','Zr-92'))+float(i.se.get(j,'iso_massf','Zr-94'))+float(i.se.get(j,'iso_massf','Zr-96'))
					dum_y = float(i.se.get(j,'iso_massf','Y-89'))
					dum_c = float(i.se.get(j,'iso_massf','C-12'))
					dum_o = float(i.se.get(j,'iso_massf','O-16'))
					dum_gda = float(i.se.get(j,'iso_massf','Gd-152'))
					dum_gdb = float(i.se.get(j,'iso_massf','Gd-154'))
					dum_zra = float(i.se.get(j,'iso_massf','Zr-96'))
					dum_zrb = float(i.se.get(j,'iso_massf','Zr-94'))

					dum_c=(float((i.se.get((int(j)),'iso_massf','C-12')))+float((i.se.get(j,'iso_massf','C-13'))))
					dum_o=(float((i.se.get((int(j)),'iso_massf','O-16')))+float((i.se.get(j,'iso_massf','O-17'))))

					dum_c12=(float((i.se.get((int(j)),'iso_massf','C-12'))))
					dum_c13=(float((i.se.get((int(j)),'iso_massf','C-13'))))

					dum_co.append(((dum_c/dum_o)*(16./12.)))

					for jj in ls_element:
						dum_ls = dum_ls + float(i.se.get(j,'iso_massf',jj)) #/float(len(ls_element)))
					for jj in hs_element:
						dum_hs = dum_hs + float(i.se.get(j,'iso_massf',jj))  #/float(len(hs_element)))
					for jj in s_element:
						dum_s = dum_s + float(i.se.get(j,'iso_massf',jj)) #/float(len(s_element)))
					
					dum_s_fe.append(log10((dum_s/dum_fe)/(s_ini[k]/fe_ini[k])))
					dum_ls_fe.append(log10((dum_ls/dum_fe)/(ls_ini[k]/fe_ini[k])))
					dum_hs_fe.append(log10((dum_hs/dum_fe)/(hs_ini[k]/fe_ini[k])))
					dum_pb_fe.append(log10((dum_pb/dum_fe)/(pb_ini[k]/fe_ini[k])))
					dum_hs_ls.append(log10((dum_hs/dum_ls)/(hs_ini[k]/ls_ini[k])))
					dum_rb_sr.append(log10((dum_rb/dum_sr)/(rb_ini[k]/sr_ini[k])))
					dum_rb_zr.append(log10((dum_rb/dum_zr)/(rb_ini[k]/zr_ini[k])))
					dum_sr_y.append(log10((dum_sr/dum_y)/(sr_ini[k]/y_ini[k])))
					dum_sr_zr.append(log10((dum_sr/dum_zr)/(sr_ini[k]/zr_ini[k])))
					dum_gd.append((dum_gda/dum_gdb)/(gda_ini[k]/gdb_ini[k]))
					dum_zrr.append((dum_zra/dum_zrb)/(zra_ini[k]/zrb_ini[k]))
			#                dum_rb_fe.append(log10((dum_rb/dum_fe)/(rb_ini[k]/fe_ini[k]))-0.2) ## correction for Lambert 1995 data
					dum_rb_fe.append(log10((dum_rb/dum_fe)/(rb_ini[k]/fe_ini[k]))-0.068) ## correction for Zamora 2009 data

					if (len(dum_co)>1):
						if (dum_co[len(dum_co)-1]>(dum_co[len(dum_co)-2]+0.02)):
							dum_sr_zr_tps.append(log10((dum_sr/dum_zr)/(sr_ini[k]/zr_ini[k])))
							dum_sr_y_tps.append(log10((dum_sr/dum_y)/(sr_ini[k]/y_ini[k])))
							dum_s_fe_tps.append(log10((dum_s/dum_fe)/(s_ini[k]/fe_ini[k])))
							dum_ls_fe_tps.append(log10((dum_ls/dum_fe)/(ls_ini[k]/fe_ini[k])))
							dum_hs_fe_tps.append(log10((dum_hs/dum_fe)/(hs_ini[k]/fe_ini[k])))
							dum_pb_fe_tps.append(log10((dum_pb/dum_fe)/(pb_ini[k]/fe_ini[k])))
							dum_hs_ls_tps.append(log10((dum_hs/dum_ls)/(hs_ini[k]/ls_ini[k])))
							dum_rb_sr_tps.append(log10((dum_rb/dum_sr)/(rb_ini[k]/sr_ini[k])))
							dum_rb_zr_tps.append(log10((dum_rb/dum_zr)/(rb_ini[k]/zr_ini[k])))
			#                                dum_rb_fe_tps.append(log10((dum_rb/dum_fe)/(rb_ini[k]/fe_ini[k]))-0.2) ## correction for Lambert 1995 data
							dum_rb_fe_tps.append(log10((dum_rb/dum_fe)/(rb_ini[k]/fe_ini[k]))-0.068) ## correction for Zamora 2009 data
							dum_gd_tps.append((((dum_gda/dum_gdb)/(gda_ini[k]/gdb_ini[k]))-1)*1000.)
							dum_zr_tps.append((((dum_zra/dum_zrb)/(zra_ini[k]/zrb_ini[k]))-1)*1000.)

							if (dum_co[len(dum_co)-1]>1.):
								dum_sr_zr_tps_co.append(log10((dum_sr/dum_zr)/(sr_ini[k]/zr_ini[k])))
								dum_sr_y_tps_co.append(log10((dum_sr/dum_y)/(sr_ini[k]/y_ini[k])))
								dum_s_fe_tps_co.append(log10((dum_s/dum_fe)/(s_ini[k]/fe_ini[k])))
								dum_ls_fe_tps_co.append(log10((dum_ls/dum_fe)/(ls_ini[k]/fe_ini[k])))
								dum_hs_fe_tps_co.append(log10((dum_hs/dum_fe)/(hs_ini[k]/fe_ini[k])))
								dum_pb_fe_tps_co.append(log10((dum_pb/dum_fe)/(pb_ini[k]/fe_ini[k])))
								dum_hs_ls_tps_co.append(log10((dum_hs/dum_ls)/(hs_ini[k]/ls_ini[k])))
								dum_rb_sr_tps_co.append(log10((dum_rb/dum_sr)/(rb_ini[k]/sr_ini[k])))
								dum_rb_zr_tps_co.append(log10((dum_rb/dum_zr)/(rb_ini[k]/zr_ini[k])))
			###                                        dum_rb_fe_tps_co.append(log10((dum_rb/dum_fe)/(rb_ini[k]/fe_ini[k]))-0.2)
								dum_rb_zr_tps_co.append(log10((dum_rb/dum_zr)/(rb_ini[k]/zr_ini[k])))
			#                                        dum_rb_fe_tps_co.append(log10((dum_rb/dum_fe)/(rb_ini[k]/fe_ini[k]))-0.2) ## correction for Lambert 1995 data
								dum_rb_fe_tps_co.append(log10((dum_rb/dum_fe)/(rb_ini[k]/fe_ini[k]))-0.068) ## correction for Zamora 2009 data
								dum_gd_tps_co.append((((dum_gda/dum_gdb)/(gda_ini[k]/gdb_ini[k]))-1)*1000.)
								dum_zr_tps_co.append((((dum_zra/dum_zrb)/(zra_ini[k]/zrb_ini[k]))-1)*1000.)
							jjjj=jjjj+1
					
					#
				s_fe.append(dum_s_fe)
				ls_fe.append(dum_ls_fe)
				hs_fe.append(dum_hs_fe)
				pb_fe.append(dum_pb_fe)
				hs_ls.append(dum_hs_ls)
				rb_sr.append(dum_rb_sr)
				rb_zr.append(dum_rb_zr)
				rb_fe.append(dum_rb_fe)
				sr_y.append(dum_sr_y)
				sr_zr.append(dum_sr_zr)
				sr_zr_tps.append(dum_sr_zr_tps)
				sr_y_tps.append(dum_sr_y_tps)
				sr_zr_tps_co.append(dum_sr_zr_tps_co)
				sr_y_tps_co.append(dum_sr_y_tps_co)
				s_fe_tps.append(dum_s_fe_tps)
				ls_fe_tps.append(dum_ls_fe_tps)
				hs_fe_tps.append(dum_hs_fe_tps)
				pb_fe_tps.append(dum_pb_fe_tps)
				hs_ls_tps.append(dum_hs_ls_tps)
				rb_sr_tps.append(dum_rb_sr_tps)
				rb_zr_tps.append(dum_rb_zr_tps)
				rb_fe_tps.append(dum_rb_fe_tps)
				s_fe_tps_co.append(dum_s_fe_tps_co)
				ls_fe_tps_co.append(dum_ls_fe_tps_co)
				hs_fe_tps_co.append(dum_hs_fe_tps_co)
				pb_fe_tps_co.append(dum_pb_fe_tps_co)
				hs_ls_tps_co.append(dum_hs_ls_tps_co)
				rb_sr_tps_co.append(dum_rb_sr_tps_co)
				rb_zr_tps_co.append(dum_rb_zr_tps_co)
				rb_fe_tps_co.append(dum_rb_fe_tps_co)
				gd.append(dum_gd)
				zr.append(dum_zrr)
				gd_tps.append(dum_gd_tps)
				zr_tps.append(dum_zr_tps)
				gd_tps_co.append(dum_gd_tps_co)
				zr_tps_co.append(dum_zr_tps_co)
				co.append(dum_co)
				k = k+1

			if plotting:

				mass_label =[]
				for i in run_mass:
					mini=float(i.se.get('mini'))
					zini=float(i.se.get('zini'))
					label=str(mini)+'$M_{\odot}$, Z='+str(zini)
					mass_label.append(label)

				metallicity_label=['3 Msun, Z=0.02','3 Msun, Z=0.02, no mol diff','3 Msun, Z=0.01','M3 set1.2']

				fig = plt.figure(0)            # Figure object
				ax = fig.add_subplot(1,1,1)     # Axes object: one row, one column, first plot (one plot!)

				mpl.rcParams['xtick.major.size'] = 20
				#mpl.rcParams['xtick.major.width'] = 4
				mpl.rcParams['xtick.minor.size'] = 10
				#mpl.rcParams['xtick.minor.width'] = 2
				
				mpl.rcParams['ytick.major.size'] = 20
				#mpl.rcParams['ytick.major.width'] = 4
				mpl.rcParams['ytick.minor.size'] = 10
				#mpl.rcParams['ytick.minor.width'] = 2
				
				if hsfe_hsls:
					array_to_plot_x = hs_ls
					array_to_plot_x_tps = hs_ls_tps
					array_to_plot_x_co = hs_ls_tps_co
					
					array_to_plot_y = hs_fe
					array_to_plot_y_tps = hs_fe_tps
					array_to_plot_y_co = hs_fe_tps_co

					pl.axhline(y=0,linewidth=2, color='k')
					pl.axvline(x=0,linewidth=2, color='k')
					for k in range(0,len(array_to_plot_x)):
						#if k > 0:
						pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=3.,color=color1)
						#pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=3.,label=str(metallicity_label[k]))
						#if k >= 1:
						pl.plot(array_to_plot_x_tps[k],array_to_plot_y_tps[k],markers[k],markersize=12.,color=color1)
						pl.plot(array_to_plot_x_co[k],array_to_plot_y_co[k],markers[k],markersize=20.,color=color1,label=str(mass_label[k]))

					pl.legend(numpoints=1,loc=2,prop={'size':20})

					pl.xlabel('[hs/ls]', fontsize=20)
					pl.ylabel('[hs/Fe]', fontsize=20)
					y_min = -0.1
					y_max = 1.0
					pl.ylim(y_min,y_max)
					pl.xlim(-0.1,0.5)

					size=20
					pl.xticks(size=size)
					pl.yticks(size=size)

					ax = pl.gca()

					for line in ax.xaxis.get_ticklines():
						line.set_markersize(20)
						line.set_markeredgewidth(3)

					for line in ax.yaxis.get_ticklines():
						line.set_markersize(20)
						line.set_markeredgewidth(3)

				if lsfe_hsls:

					array_to_plot_x = hs_ls
					array_to_plot_x_tps = hs_ls_tps
					array_to_plot_x_co = hs_ls_tps_co
					
					array_to_plot_y = ls_fe
					array_to_plot_y_tps = ls_fe_tps
					array_to_plot_y_co = ls_fe_tps_co

					pl.axhline(y=0,linewidth=2, color='k')
					pl.axvline(x=0,linewidth=2, color='k')
					for k in range(0,len(array_to_plot_x)):
						#                if k > 0:
						pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=4.,label=str(mass_label[k]))
						pl.plot(array_to_plot_x_tps[k],array_to_plot_y_tps[k],markers[k],markersize=12.)
						pl.plot(array_to_plot_x_co[k],array_to_plot_y_co[k],markers[k],markersize=20.)

					pl.legend(numpoints=1,loc='upper right',prop={'size':20})

					pl.xlabel('[hs/ls]', fontsize=20)
					pl.ylabel('[ls/Fe]', fontsize=20)
					y_min = -0.1
					y_max = 0.8
					pl.ylim(y_min,y_max)
					pl.xlim(-0.2,0.5)

					size=20
					pl.xticks(size=size)
					pl.yticks(size=size)

					ax = pl.gca()

					for line in ax.xaxis.get_ticklines():
						line.set_markersize(20)
						line.set_markeredgewidth(3)

					for line in ax.yaxis.get_ticklines():
						line.set_markersize(20)
						line.set_markeredgewidth(3)

				if rbsr_hsls:
					array_to_plot_x = hs_ls
					array_to_plot_y = rb_sr

					pl.axhline(y=0,linewidth=2, color='k')
					pl.axvline(x=0,linewidth=2, color='k')
					for k in range(0,len(array_to_plot_x)):
						#                if k > 0:
						if k >= 1:
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=4.,label=str(mass_label[k]))
						else:
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=4.,label=str(mass_label[k]))

					pl.legend(numpoints=1,loc='upper right',prop={'size':20})

					pl.xlabel('[hs/ls]', fontsize=20)
					pl.ylabel('[Rb/Sr]', fontsize=20)
					y_min = -0.4
					y_max = 1.5
					x = ax.get_position()

					pl.ylim(y_min,y_max)
					pl.xlim(-0.7,0.4)

					size=30
					pl.xticks(size=size)
					pl.yticks(size=size)

					ax = pl.gca()

					for line in ax.xaxis.get_ticklines():
						line.set_markersize(20)
						line.set_markeredgewidth(3)
					for line in ax.yaxis.get_ticklines():
						line.set_markersize(20)
						line.set_markeredgewidth(3)

				if rbzr_hsls:
					array_to_plot_x = hs_ls
					array_to_plot_y = rb_zr

					pl.axhline(y=0,linewidth=2, color='k')
					pl.axvline(x=0,linewidth=2, color='k')
					for k in range(0,len(array_to_plot_x)):
			#                if k > 0:
						if k >= 1:
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=2.,label=str(mass_label[k]))
							pl.plot(hs_ls_tps[k],rb_zr_tps[k],markers[k],markersize=12.)
							pl.plot(hs_ls_tps_co[k],rb_zr_tps_co[k],markers[k],markersize=20.)
						else:
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=2.,label=str(mass_label[k]))
							pl.plot(hs_ls_tps[k],rb_zr_tps[k],markers[k],markersize=12.)
							pl.plot(hs_ls_tps_co[k],rb_zr_tps_co[k],markers[k],markersize=20.)

					pl.legend(numpoints=1,loc='upper right',prop={'size':20})

					pl.xlabel('[hs/ls]', fontsize=20)
					pl.ylabel('[Rb/Zr]', fontsize=20)
					y_min = -0.4
					y_max = 0.4
					pl.ylim(y_min,y_max)
					pl.xlim(-0.5,0.5)

					size=20
					pl.xticks(size=size)
					pl.yticks(size=size)

					ax = pl.gca()

					for line in ax.xaxis.get_ticklines():
						line.set_markersize(20)
						line.set_markeredgewidth(3)

					for line in ax.yaxis.get_ticklines():
						line.set_markersize(20)
						line.set_markeredgewidth(3)

				if rbfe_sfe:

					array_to_plot_x = s_fe
					array_to_plot_y = rb_fe

					pl.axhline(y=0,linewidth=3, color='k')
					pl.axvline(x=0,linewidth=3, color='k')
					for k in range(0,len(array_to_plot_x)):
			#                if k > 0:
						if k >= 1:
			#                        pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=2.,label=str(mass_label[k])+'Msun '+str(metallicity_label[k]))
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=2.,label=str(metallicity_label[k]))
							pl.plot(s_fe_tps[k],rb_fe_tps[k],markers[k],markersize=12.)
							pl.plot(s_fe_tps_co[k],rb_fe_tps_co[k],markers[k],markersize=20.)
						else:
			#                        pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=2.,label=str(mass_label[k])+'Msun '+str(metallicity_label[k]))
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=2.,label=str(metallicity_label[k]))
							pl.plot(s_fe_tps[k],rb_fe_tps[k],markers[k],markersize=12.)
							pl.plot(s_fe_tps_co[k],rb_fe_tps_co[k],markers[k],markersize=20.)

					pl.legend(numpoints=1,loc=4,prop={'size':20})

					pl.xlabel('[s/Fe]', fontsize=20)
					pl.ylabel('[Rb/Fe]', fontsize=20)
					y_min = -0.4
					y_max = 0.4
					pl.ylim(y_min,y_max)
					pl.xlim(-0.1,2.0)

					size=20
					pl.xticks(size=size)
					pl.yticks(size=size)

					ax = pl.gca()

					for line in ax.xaxis.get_ticklines():
						line.set_markersize(25)
						line.set_markeredgewidth(3)

					for line in ax.yaxis.get_ticklines():
						line.set_markersize(25)
						line.set_markeredgewidth(3)

				if sry_srzr:
					array_to_plot_x = sr_zr
					array_to_plot_y = sr_y

					pl.axhline(y=0,linewidth=2, color='k')
					pl.axvline(x=0,linewidth=2, color='k')
					for k in range(0,len(array_to_plot_x)):
						if k >= 1:
							#                if k > 0:
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=3.5,label=str(metallicity_label[k]))
							pl.plot(sr_zr_tps[k],sr_y_tps[k],markers[k],markersize=12.)
							pl.plot(sr_zr_tps_co[k],sr_y_tps_co[k],markers[k],markersize=20.)
						else:
							pl.plot(array_to_plot_x[k],array_to_plot_y[k],symbol[k],markersize=10.,linewidth=3.5,label=str(metallicity_label[k]))
							pl.plot(sr_zr_tps[k],sr_y_tps[k],markers[k],markersize=12.,label='C/O < 1')
							pl.plot(sr_zr_tps_co[k],sr_y_tps_co[k],markers[k],markersize=20.,label='C/O > 1')

					pl.xlabel('[Sr/Zr]', fontsize=20)
					pl.ylabel('[Sr/Y]', fontsize=20)
					y_min = -0.1
					y_max = 0.1
					pl.ylim(y_min,y_max)
					pl.xlim(-0.3,0.3)
					
					size=20
					pl.xticks(size=size)
					pl.yticks(size=size)

					ax = pl.gca()

					for line in ax.xaxis.get_ticklines():
						line.set_markersize(25)
						line.set_markeredgewidth(3)

					for line in ax.yaxis.get_ticklines():
						line.set_markersize(25)
						line.set_markeredgewidth(3)

					box = ax.get_position()

					ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

					pl.legend(loc=2, borderaxespad=1.0,prop={'size':25})

			box = ax.get_position()
	

	def set_plot_CO_mass(self,fig=3123,xaxis='mass',linestyle=['-'],marker=['o'],color=['r'],age_years=True,sparsity=500,markersparsity=200,withoutZlabel=False,t0_model=[]):
		
		'''
			PLots C/O surface number fraction
		'''

		if len(t0_model)==0:
			t0_model = len(self.runs_H5_surf)*[0]

		plt.figure(fig)
                for i in range(len(self.runs_H5_surf)):
                        sefiles=se(self.runs_H5_surf[i])
			cycles=range(int(sefiles.se.cycles[0]),int(sefiles.se.cycles[-1]),sparsity)
			mini=sefiles.get("mini")
			zini=sefiles.get("zini")
			label=str(mini)+'$M_{\odot}$, Z='+str(zini)
			if xaxis=='cycles':
				x=cycles
			if xaxis=='age':
				x=sefiles.get(cycles,'age')
                        	if age_years==True:
                                	x=np.array(x)*sefiles.get('age_unit')/(365*24*3600)
				x = x - x[t0_model[i]]
			if xaxis=='mass':
				x=sefiles.get(cycles,'mass')
			x=x[t0_model[i]:]
			c12=sefiles.get(cycles,'C-12')[t0_model[i]:]
			o16=sefiles.get(cycles,'O-16')[t0_model[i]:]
			if withoutZlabel==True:
				plt.plot(x,4./3.*np.array(c12)/np.array(o16),label=label.split(',')[0],marker=marker[i],linestyle=linestyle[i],markevery=markersparsity,color=color[i])

			else:
				plt.plot(x,4./3.*np.array(c12)/np.array(o16),label=label,marker=marker[i],linestyle=linestyle[i],markevery=markersparsity,color=color[i])
			if xaxis=='mass':
				plt.xlim(7,0.5)
				#plt.gca().invert_xaxis()	
	                        plt.xlabel('$M/M_{\odot}$',fontsize=18)
	
                        plt.ylabel('C/O Ratio', fontsize=18)
			plt.legend(loc=1)

	def set_plot_surface_abu_ratios(self,fig=2,x_species=['C-12','C-13'],y_species=['Zr-96','Zr-94'], samefigure=True,marker=['o'],linestyle=['-'],markevery=500,samefigureall=True,sparsity=100):
		import nugridse as mp	
                plt.figure(fig)
                for i in range(len(self.runs_H5_surf)):
			idx=0
                        sefiles=mp.se(self.runs_H5_surf[i])
			if samefigure==True:
				plt.figure(fig)
			else:
				plt.figure(i)
			cycles=range(int(sefiles.se.cycles[0]),int(sefiles.se.cycles[-1]),sparsity)
	
			x_species1=sefiles.get(cycles,x_species[0])
			x_species2=sefiles.get(cycles,x_species[1])
			ratio_x=np.array(x_species1)/np.array(x_species2)
			
                        y_species1=sefiles.get(cycles,y_species[0])
                        y_species2=sefiles.get(cycles,y_species[1])
                        ratio_y=np.array(y_species1)/np.array(y_species2)
			mini=sefiles.get('mini')
			zini=sefiles.get('zini')
                        label=str(mini)+'$M_{\odot}$, Z='+str(zini)
			plt.plot(ratio_x,ratio_y,marker=marker[i],linestyle=linestyle[i],markevery=markevery,label=label)
			
		xlabel1=x_species[0]+'/'+x_species[1]
		ylabel1=y_species[0]+'/'+y_species[1]
		plt.xlabel(xlabel1)
		plt.ylabel(ylabel1)	
	

	def set_plot_surface_abu(self,fig=2,species=['Sr-88','Ba-136'],decay=False,number_frac=False,xaxis='cycles',age_years=False,ratio=False,sumiso=False,eps=False,samefigure=False,samefigureall=False,withkip=False,sparsity=200,linestyle=['--'],marker=['o'],color=['r'],label=[],markevery=100,t0_model=-1,savefig=''):
		
		'''	
			Simply plots surface abundance versus model number or time

		'''
		extralabel=False
		if len(label)>0:
			extralabel=True
                import nugridse as mp
                import utils as u
		idx=0
		if eps==True:
			species= species + ['H-1']
		if samefigureall==True and ratio==False:
			plt.figure(fig)
                for i in range(len(self.runs_H5_surf)):
			idx=0
                        sefiles=mp.se(self.runs_H5_surf[i])
			if samefigure==True:
				plt.figure(i)
			cycles=range(int(sefiles.se.cycles[0]),int(sefiles.se.cycles[-1]),sparsity)
			mini=sefiles.get("mini")
			zini=sefiles.get("zini")
			if not extralabel:
				label1=str(mini)+'$M_{\odot}$, Z='+str(zini)
			if xaxis=='cycles':
				x=cycles
			if xaxis=='age':
				x=sefiles.get(cycles,'age')
				if age_years==True:
					x=np.array(x)*sefiles.get('age_unit')/(365*24*3600)
					print 'x is age'
				if t0_model>0:
					#print cycles
					idxt0=0
					for kk in range(len(cycles)):
						print int(cycles[kk]),t0_model
						if int(cycles[kk]) == t0_model:
							idxt0=kk
							break
					print 'found t0_model idx',idxt0
					#idxt0=cycles.index(t0_model)
					cycles=cycles[idxt0:]
					if idxt0==0:
						print 'Warning, t0modle not found'
					x=x[idxt0:] - x[idxt0]
				else:
					idxt0=0
			if xaxis=='mass':
				x=sefiles.get(cycles,'mass')
			if decay==False:
				
				species_abu1=sefiles.get(cycles,species)
			else:
				species_abu11=sefiles.get(cycles,'iso_massf_decay')
				species_abu1=[]
				for jj in range(len(cycles)):
				    species_abu1.append([])	
				    for j in range(len(species)):
					species_abu1[-1].append(species_abu11[jj][sefiles.se.isotopes.index(species[j])])

			if len(species)==1:
				species_abu11=[]
				for kk in range(len(species_abu1)):
					species_abu11.append([species_abu1[kk]])
				species_abu1=species_abu11	
			species_abu=[]
			for k in range(len(species)):
				print 'species ',k
				species_abu.append([])
			for k in range(len(species)):
				for h in range(len(cycles)):
					species_abu[k].append(species_abu1[h][k])
			#print species_abu
			#if t0_model>0:
		#		species_abu=species_abu[t0_model:]
			for k in range(len(species)):
				if samefigure==False and ratio==False:
					fig=plt.figure(species[k])
				if xaxis=='cycles':
					plt.xlabel('model number')
				if xaxis=='age':
					plt.xlabel('Age [yrs]')
				if xaxis=='mass':
					plt.xlabel('$M/M_{\odot}$')
				plt.ylabel('X$_i$')
				if ratio==True:
					continue
				if extralabel:
					label=label[k]
				else:
					label=label1
				if samefigure==True:
					if sumiso == True:
						sumiso_massfrac=np.array(species_abu[0])
						for hh in range(1,len(species_abu)):
							sumiso_massfrac = sumiso_massfrac + np.array(species_abu[hh])
						plt.plot(x,sumiso_massfrac,linestyle=linestyle[idx],marker=marker[idx],label=species[k]+', '+label,color=color[idx],markevery=markevery)
						break #leave iso looop 
					
					else:
                                                if eps==True:
                                                        species_abu[0]=np.log10(np.array(species_abu[0])/(np.array(species_abu[1])*7))   + 12.

						plt.plot(x,species_abu[k],linestyle=linestyle[idx],marker=marker[idx],label=species[k]+', '+label,color=color[idx],markevery=markevery)
					idx+=1
					if eps==True:
						break

				else:
					if withkip==True:
						print 'test'
					else:
						plt.ylabel('X('+species[k]+')')
                                                if eps==True:
                                                        species_abu[0]=np.log10(np.array(species_abu[0])/(np.array(species_abu[1])*7))   + 12.

						plt.plot(x,species_abu[k],linestyle=linestyle[i],marker=marker[i],label=label,color=color[i],markevery=markevery)
                                        if eps==True:
                                                break
						
				plt.legend(loc=2)
				plt.yscale('log')
			if ratio==True:
				if number_frac==True:
					print 'plot number frac'
					plt.plot(x,4./3.*np.array(species_abu[1])/np.array(species_abu[0]),linestyle=linestyle[i],marker=marker[i],label=label,color=color[i],markevery=markevery)
				else:
					plt.plot(x,np.array(species_abu[1])/np.array(species_abu[0]),linestyle=linestyle[i],marker=marker[i],label=label,color=color[i],markevery=markevery)
                                plt.legend(loc=2)
                                plt.yscale('log')
			name='M'+str(mini)+'Z'+str(zini)
			plt.legend(loc=4)
			plt.savefig(savefig+'/surf_'+name+'.png')



	def set_write_chris_prof(self,table='profiles_chris.txt'):


		import nugridse as mp
		out = '' #9fields
		out = out + ' Zone '+ ' & ' + ' R[Rsun] '+' & ' +' M[cgs] '+' & ' +' P[cgs]  '+' & ' +'Rho[cgs] '+' & ' +'  T[K]   '+r'\n'
                for i in range(len(self.runs_H5)):
                        sefiles=mp.se(self.runs_H5_restart[i])
                        mini=sefiles.get("mini")
                        zini=sefiles.get("zini")
			name='M='+str(mini)+'Z='+str(zini)
			out = out + (name+r'\n')
			lastcyc=sefiles.se.cycles[-1]
			print 'last cycle',lastcyc
			T=sefiles.get(lastcyc,'temperature')*sefiles.get('temperature_unit')
			R=sefiles.get(lastcyc,'radius')*sefiles.get('radius_unit')
			rho=sefiles.get(lastcyc,'rho')*sefiles.get('rho_unit')
			P=sefiles.get(lastcyc,'pressure')
			M=sefiles.get(lastcyc,'mass')*sefiles.get('mass_unit')
			out1=''
			for k in range(len(M)):
				print 'idx ',k			
				print 'R ',R[k]
				print 'M ',M[k]
				print 'P ',
				print 'P ',P[k]
				print 'rho ',rho[k]
				print 'T ',T[k]	
				out1= out1 + ( format(k, '05d')+ ' & '+'{:.3E}'.format(R[k]) + ' & '+'{:.3E}'.format(M[k]) + ' & '+'{:.3E}'.format(P[k]) + ' & '+'{:.3E}'.format(rho[k]) + ' & '+'{:.3E}'.format(T[k]) +r'\n' )			
			out = out + out1
		f1=open(table,'w')
		f1.write(out)



	def set_plot_profile_decay(self,cycles=20*[-1],mass_range=20*[[0,0]],ylim=20*[[0,0]],isotopes=[],linestyle=[],save_dir=''): 
		'''
			Plots HRDs


			end_model - array, control how far in models a run is plottet, if -1 till end
			symbs_1  - set symbols of runs
		'''
		if len(linestyle)==0:
			linestyle=200*['-']
		import nugridse as mp
		import utils as u
		#print self.runs_H5_restart
                for i in range(len(self.runs_H5_restart)):
        		sefiles=mp.se(self.runs_H5_restart[i])
			cycle=cycles[i]
                        if cycle==-1:
                                cycle=int(sefiles.se.cycles[-1])			
			if mass_range[i][0] ==0 and mass_range[i][1]==0:
				mass_range[i][1]=sefiles.get(cycle,'mass')[-1]				
			sefiles.read_iso_abund_marco(mass_range[i],cycle)
		    	u.stable_specie()
		    	sefiles.decay(sefiles.mass_frac)
			idx_species=[]
			for k in range(len(isotopes)):
				other_name_scheme=isotopes[k].split("-")[0].upper()+(5-len(isotopes[k])+1)*" "+isotopes[k].split("-")[1]
				#other_name_scheme=other_name_scheme.capitalize()

				idx_specie=u.back_ind[other_name_scheme]
				idx_species.append(idx_specie)
			mass_abu_array=[]
			for idx_specie in idx_species:
				mass_abu_array.append([])
				for idx_mass in range(len(mp.decayed_multi_d)):
					mass_abu_array[-1].append(mp.decayed_multi_d[idx_mass][idx_specie])
					
			#plotting
			plt.figure(self.run_dirs_name[i])
			#print len(mp.used_masses),len(mass_abu_array[0])
			#print mass_abu_array[0]
			for k in range(len(isotopes)):
				plt.plot(mp.used_masses,mass_abu_array[k],linestyle=linestyle[k],label=isotopes[k])
			plt.legend()
			plt.yscale('log')
			#print sefiles.get(cycle,'mass')[-1]
			plt.xlabel('M/Msun')
			plt.ylabel('$X_i$')
			plt.xlim(mass_range[i][0],mass_range[i][1])
			if (ylim[i][0]>0 or ylim[i][1]>0) or (ylim[i][0]>0 and ylim[i][1]>0):
				plt.ylim(ylim[i][0],ylim[i][1])
			if len(save_dir)>0:
				star_mass=sefiles.get("mini")
				star_z=sefiles.get("zini")
				plt.savefig(save_dir+'/'+self.run_dirs_name[i]+'_decay_profiles.png')

	def set_get_abu_distr_decay(self,cycles=20*[-1],mass_range=20*[[0,0]],ylim=20*[[0,0]],isotopes=['all']):
                import nugridse as mp
                import utils as u
                print self.runs_H5_restart
                massfrac_all=[]
                iso_all=[]
                for i in range(len(self.runs_H5_restart)):		
                        sefiles=mp.se(self.runs_H5_restart[i])
                        cycle=cycles[i]
                        if cycle==-1:
                                cycle=int(sefiles.se.cycles[-1])
                        if mass_range[i][0] ==0 and mass_range[i][1]==0:
                                mass_range[i][1]=sefiles.get(cycle,'mass')[-1]
			print 'use cycle',cycle
			sefiles.average_iso_abund_marco(mass_range=mass_range[i],cycle=cycle,stable=True,i_decay=2)			
			massfrac1=mp.average_mass_frac_decay
			other_name_scheme=u.back_ind.keys()
			isos=[]
			massfrac=[]
			print len(massfrac),len(other_name_scheme)
			for kk in range(len(other_name_scheme)):
				list1=re.split('(\d+)',other_name_scheme[kk])
				newname=list1[0].capitalize().strip()+'-'+list1[1].strip()
				#print other_name_scheme[kk],newname,massfrac[kk]
				isos.append(newname)
				massfrac.append(massfrac1[u.back_ind[other_name_scheme[kk]]])
			massfrac_all.append(massfrac)
			iso_all.append(isos)
		return iso_all,massfrac_all
		

	def set_get_abu_distr_decay_old(self,cycles=20*[-1],mass_range=20*[[0,0]],ylim=20*[[0,0]],isotopes=['all'],linestyle=[],save_dir=''): 
		'''
			Plots HRDs


			end_model - array, control how far in models a run is plottet, if -1 till end
			symbs_1  - set symbols of runs
		'''
		if len(linestyle)==0:
			linestyle=200*['-']
		import nugridse as mp
		import utils as u
		print self.runs_H5_restart
		massfrac_all=[]
		iso_all=[]
                for i in range(len(self.runs_H5_restart)):
        		sefiles=mp.se(self.runs_H5_restart[i])
			cycle=cycles[i]
                        if cycle==-1:
                                cycle=int(sefiles.se.cycles[-1])			
			if mass_range[i][0] ==0 and mass_range[i][1]==0:
				mass_range[i][1]=sefiles.get(cycle,'mass')[-1]				
			sefiles.read_iso_abund_marco(mass_range[i],cycle)
		    	u.stable_specie()
		    	sefiles.decay(sefiles.mass_frac)
			idx_species=[]
			massfrac=[]
			iso=[]
			if not isotopes[0]=='all':
				for k in range(len(isotopes)):
					other_name_scheme=isotopes[k].split("-")[0].upper()+(5-len(isotopes[k])+1)*" "+isotopes[k].split("-")[1]
					#other_name_scheme=other_name_scheme.capitalize()

					idx_specie=u.back_ind[other_name_scheme]
					idx_species.append(idx_specie)
					massfrac.append(average_massfrac_decay[idx_specie])
				iso=isotopes
			else:
				massfrac=mp.average_mass_frac_decay
				other_name_scheme=u.back_ind
				iso=[]
				import re
				for kk in range(len(other_name_scheme)):
					list1=re.split('(\d+)',other_name_scheme[kk])
					newname=list1[0].capitalize()+'-'+list1[1]
					iso.append(newname)
			massfrac_all.append(massfrac)
			iso_all.append(iso)

		return iso_all,massfrac_all

	def set_cores_massive(self,filename='core_masses_massive.txt'):

		'''
			Uesse function cores in nugridse.py
		'''
		
		core_info=[]
		minis=[]
                for i in range(len(self.runs_H5_surf)):
                        sefiles=se(self.runs_H5_out[i])
			mini=sefiles.get('mini')
			minis.append(mini)
			incycle=int(sefiles.se.cycles[-1])
			core_info.append(sefiles.cores(incycle=incycle))
		print_info=''
                for i in range(len(self.runs_H5_surf)):
			if i ==0:
				print 'Following returned for each initial mass'
				print core_info[i][1]
                        #print '----Mini: ',minis[i],'------'
			print_info+=(str(minis[i])+' & ')
			info=core_info[i][0]
			for k in range(len(info)):
				print_info+=('{:.3E}'.format(float(core_info[i][0][k]))+' & ')
			print_info=(print_info+'\n')
			#print core_info[i][2]
		f1=open(filename,'a')
		f1.write(print_info)
		f1.close()

	def get_fallback_coord(self,isotope='Ni-56',masslimit=0.1,masscutlim=False,delay=True):

		'''
			Returns fallback mass coordinate so that the amount of masslimit
			of the isotope isotope is ejected. Explosion type chosen with delay option.
			
			masscutlim: If true, new fallback coordinate can only as small as the original fallback
				  prescription by C. Fryer. Useful for more massive stars which would not eject any metals
		 		with Freyer's prescription.
		'''

		def getmasscut(m_ini,z_ini,delay): 
			if int(m_ini)==12:
				m_ini=15
			z_metal=z_ini/0.02

			print 'MINI',m_ini,z_metal
			if ((m_ini>=11.) and (m_ini<30.)):
				if delay==True:
					mass_cut = 1.1 + 0.2*np.exp((m_ini-11.0)/4.) - (2.0 + z_metal)*np.exp(0.4*(m_ini -26.0))
				####rapid cc
				else:
					if m_ini<22.:
						mass_cut= 1.1 +0.2*np.exp((m_ini-11.0)/7.5) + 10*(1.0+z_metal)*np.exp(-(m_ini-23.5)**2/(1.0+z_metal)**2)
					elif m_ini<30 :
						mass_cut=  1.1 + 0.2*np.exp((m_ini-11.0)/4.) - (2.0 + z_metal)*np.exp(0.4*(m_ini -26.0)) - 1.85 + 0.25*z_metal +10.0*(1.0+z_metal)*np.exp(-(m_ini-23.5)**2/(1.0+z_metal)**2)
			##at higher mass difference
			elif ((m_ini>30) and (m_ini<50)):
				#delay
				if delay==True:
					mass_cut= min( 33.35 + (4.75 + 1.25*z_metal)*(m_ini-34.),m_ini-z_metal**0.5 *(1.3*m_ini - 18.35))
				else:
					mass_cut = min( 33.35 + (4.75 + 1.25*z_metal)*(m_ini-34.),m_ini-z_metal**0.5 *(1.3*m_ini - 18.35)) - 1.85 + z_metal*(75. -m_ini)/20.
			elif m_ini>50:
				#Page 7, Fryer12, only at solar Z valid
				if z_metal==1:
					if m_ini<90.:
						mass_cut = 1.8 + 0.04*(90. - m_ini)
					else:
						mass_cut = 1.8 + np.log10(m_ini - 89.)
				#The part below will probably never be used
				if z_metal <1:
					if m_ini<90.:
						mass_cut = max(min( 33.35 + (4.75 + 1.25*z_metal)*(m_ini-34.),m_ini-z_metal**0.5 *(1.3*m_ini - 18.35)),1.8 + 0.04*(90. - m_ini))
					else:
						mass_cut = max(min( 33.35 + (4.75 + 1.25*z_metal)*(m_ini-34.),m_ini-z_metal**0.5 *(1.3*m_ini - 18.35)),1.8 + np.log10(m_ini - 89.))

			mass_cut=round(mass_cut,2)
						   
			return mass_cut	
	
		fallback_coords=[]
		orig_fallback=[]
		minis=[]
		ni56_mass_all=[]
		o16_mass_all=[]
                for i in range(len(self.runs_H5_surf)):
                        sefiles=se(self.runs_H5_restart[i])
			m_ini=sefiles.get("mini")
			z_ini=sefiles.get("zini")
			minis.append(m_ini)
			mass_cut=getmasscut(m_ini,z_ini,delay)
			print 'mass cut',mass_cut		
			cycle=int(sefiles.se.cycles[-1])
			mass_cycle=sefiles.get(cycle,"mass")
			mass_limit=mass_cycle[-1] #maximum mass_test can be
			idx_start=min(range(len(mass_cycle)), key=lambda i: abs(mass_cycle[i]-mass_cut))
			ni56_frac=sefiles.get(cycle,'Ni-56')
			o16_frac=sefiles.get(cycle,'O-16')
			deltam=mass_cycle[1:]-mass_cycle[:-1]
			ni56_mass=0
			o16_mass=0
			ni56_mass_orig=0
			newremnant=mass_cut
			#go fromm outside to the inside
			for k in range(len(mass_cycle)-1)[::-1]:
				cellm=deltam[k]*ni56_frac[k]
				#in case fallback coordinate should not be smaller then mass_cut (fryer)
				if masscutlim==True:
					if mass_cycle[k]<mass_cut:
						break
				if ni56_mass>masslimit:
					newremnant=mass_cycle[k]
					print 'found new  remnant',newremnant,'ni56:',ni56_mass
					break
				ni56_mass+=cellm
				o16_mass+= (deltam[k]*o16_frac[k])
			if newremnant == mass_limit:
				print 'Ni-56 does not reach 0.1Msun, take old remnant',newremnant
			fallback_coords.append(newremnant)
			orig_fallback.append(mass_cut)
			ni56_mass_all.append(ni56_mass)
			o16_mass_all.append(o16_mass)
		print '########Results:######'
		for k in range(len(minis)):
			print 'Initial mass: '+str(minis[k])+'Original fallback coord (fryer): '+str(orig_fallback[k])+',New fallback coord: '+str(fallback_coords[k])+'Ni-56 ejected: '+str(ni56_mass_all[k])+'O16: '+str(o16_mass_all[k])


		return minis, fallback_coords





	def set_burnstages_upgrade_massive(self):
		'''
			Outputs burnign stages as done in burningstages_upgrade (nugridse)
		'''
		burn_info=[]
		burn_mini=[]
                for i in range(len(self.runs_H5_surf)):
                        sefiles=se(self.runs_H5_out[i])
			burn_info.append(sefiles.burnstage_upgrade())
                        mini=sefiles.get('mini')
                        #zini=sefiles.get('zini')
			burn_mini.append(mini)
		for i in range(len(self.runs_H5_surf)):
			print 'Following returned for each initial mass'
			print '[burn_cycles,burn_ages, burn_abun, burn_type,burn_lifetime]'
			print '----Mini: ',burn_mini[i],'------'
			print burn_info[i]
			
	def set_plot_CC_T_rho_new(self,fig='CC evol',linestyle=['-'],burn_limit=0.997,color=['r'],marker=['o'],markevery=500): 
		'''
			Plots 
			end_model - array, control how far in models a run is plottet, if -1 till end
			symbs_1  - set symbols of runs
		'''
		if len(linestyle)==0:
			linestyle=200*['-']
		plt.figure(fig)
                for i in range(len(self.runs_H5_surf)):
                        sefiles=se(self.runs_H5_out[i])
			t1_model=-1
                        mini=sefiles.get('mini')
			zini=sefiles.get('zini')
			label=str(mini)+'$M_{\odot}$, Z='+str(zini)
                        model=sefiles.se.cycles
                        model_list=[]
                        for k in range(0,len(model),1):
                                model_list.append(model[k])
			print 'REad Rho,T, this might take a while...'
                        rho1=sefiles.get(model_list,'rho')   #[:(t1_model-t0_model)]
                        T1=sefiles.get(model_list,'temperature')#[:(t1_model-t0_model)]
			print 'finished'
                        rho=[]
                        T=[]
			T_unit=sefiles.get('temperature_unit')
			labeldone=False
			for k in range(len(model_list)):
				#print 'test model ',model_list[k]
                                t9max=max(np.array(T1[k])*T_unit/1.e9)
                                #T.append(max(t9))
				rho1max=max(rho1[k])
				print 'model',model_list[k]
				print 'maxT, maxrho'
				print t9max,rho1max
				if k==0:
					t9_prev=t9max
					rho1_prev=rho1max
					idx_T=0
					idx_rho=0
					continue
				if t9max>t9_prev:
					idx_T=k
					t9_prev=t9max
				if rho1max>rho1_prev:
					idx_rho=k
					rho1_prev=rho1max
			print 'found highest rho',idx_rho,max(np.array(T1[idx_rho])*T_unit/1.0e9),max(rho1[idx_rho]),model_list[idx_rho]
			print 'found highest T',idx_T,max(np.array(T1[idx_T])*T_unit/1.0e9),max(rho1[idx_T]),model_list[idx_T]
			if idx_T==idx_rho:
				x=np.array(T1[idx_T])*T_unit/1e9
				y=rho1[idx_T]
				rho1=[]	
				T1=[]
				for k in range(len(x)):
				    if not y[k]==1.0:
			        	rho1.append(y[k])
        				T1.append(x[k])
				x=T1
				y=rho1
				plt.plot(x,y,label=label,color=color[i],marker=marker[i],linestyle=linestyle[i],markevery=markevery)
                                #rhoi.append(max(rho1[k]))	       
			else:
				#for max T
				x=np.array(T1[idx_T])*T_unit/1e9
				y=rho1[idx_T]
                                rho_temp=[]
                                T_temp=[]
                                for k in range(len(x)):
                                    if not y[k]==1.0:
                                        rho_temp.append(y[k])
                                        T_temp.append(x[k])
                                x=T_temp
                                y=rho_temp
				plt.plot(x,y,label=label,color=color[i],marker=marker[i],linestyle=linestyle[i],markevery=markevery)
				#for max rho
        			x=np.array(T1[idx_rho])*T_unit/1e9
				y=rho1[idx_rho]
                                rho_temp=[]
                                T_temp=[]
                                for k in range(len(x)):
                                    if not y[k]==1.0:
                                        rho_temp.append(y[k])
                                        T_temp.append(x[k])
                                x=T_temp
                                y=rho_temp
				plt.plot(x,y,color=color[i],marker=marker[i],linestyle=linestyle[i],markevery=markevery)
 

	def set_plot_CC_T_rho_max(self,linestyle=[],burn_limit=0.997,color=['r'],marker=['o'],markevery=500): 
		'''
			Plots 
			end_model - array, control how far in models a run is plottet, if -1 till end
			symbs_1  - set symbols of runs
		'''
		if len(linestyle)==0:
			linestyle=200*['-']
		plt.figure('CC evol')
                for i in range(len(self.runs_H5_surf)):
                        sefiles=se(self.runs_H5_out[i])
			t1_model=-1
			sefiles.get('temperature')
			sefiles.get('density')
                        mini=sefiles.get('mini')
			zini=sefiles.get('zini')
                        model=sefiles.se.cycles
                        model_list=[]
                        for k in range(0,len(model),1):
                                model_list.append(model[k])
                        rho1=sefiles.get(model_list,'rho')   #[:(t1_model-t0_model)]
                        T1=sefiles.get(model_list,'temperature')#[:(t1_model-t0_model)]
                        rho=[]
                        T=[]
			T_unit=sefiles.get('temperature_unit')
			labeldone=False
                        for k in range(len(model_list)):
				t9=np.array(T1[k])*T_unit/1e9
				T.append(max(t9))
				rho.append(max(rho1[k]))
			label=str(mini)+'$M_{\odot}$, Z='+str(zini)
			plt.plot(T,rho,label=label,color=color[i],marker=marker[i],markevery=markevery)					
		plt.xlabel('$T_{9,max} (GK)$')
		plt.ylabel(r'$\rho [cm^{-3}]$')
		plt.yscale('log')
		plt.xscale('log')
		plt.legend(loc=2)



	def set_plot_CC_T_rho(self,linestyle=[],burn_limit=0.997,color=['r'],marker=['o'],nolabelZ=False,markevery=500): 
		'''
			Plots HRDs
			end_model - array, control how far in models a run is plottet, if -1 till end
			symbs_1  - set symbols of runs
		'''
		if len(linestyle)==0:
			linestyle=200*['-']
		plt.figure('CC evol')
                for i in range(len(self.runs_H5_surf)):
                        sefiles=se(self.runs_H5_out[i])
			t1_model=-1
			sefiles.get('temperature')
			sefiles.get('density')
                        mini=sefiles.get('mini')
			zini=sefiles.get('zini')
                        model=sefiles.se.cycles
                        model_list=[]
                        for k in range(0,len(model),1):
                               model_list.append(model[k])
                        rho1=sefiles.get(model_list,'rho')   #[:(t1_model-t0_model)]
                        T1=sefiles.get(model_list,'temperature')#[:(t1_model-t0_model)]
                        rho=[]
                        T=[]
			T_unit=sefiles.get('temperature_unit')
			labeldone=False
                        for k in range(len(model_list)):
				t9=np.array(T1[k])*T_unit/1e9
				T_delrem=[]
				rho_delrem=[]
				if k ==0:
					T.append(max(t9))
					rho.append(max(rho1[k]))
					for h in range(len(rho1[k])):
						if rho1[k][h] < 1e3:
							T_delrem.append(t9[h])
							rho_delrem.append(rho1[k][h])
					if nolabelZ==True:
						plt.plot(T_delrem,rho_delrem,label=self.extra_label[i].split(',')[0],color=color[i],marker=marker[i],markevery=markevery)
					else:
						plt.plot(T_delrem,rho_delrem,label=self.extra_label[i],color=color[i],marker=marker[i],markevery=markevery)
				else:
					if (max(rho)<max(rho1[k]) or max(T)<max(t9)):
                                        	for h in range(len(rho1[k])):
							if rho1[k][h] > 1e3:
                                                		T_delrem.append(t9[h])
                                        	      		rho_delrem.append(rho1[k][h])
						if labeldone==True:	
							plt.plot(T_delrem,rho_delrem,color=color[i],marker=marker[i],markevery=markevery)
						else:
							label=str(mini)+'$M_{\odot}$, Z='+str(zini)
							if nolabelZ==True:
								plt.plot(T_delrem,rho_delrem,label=label.split(',')[0],color=color[i],marker=marker[i],markevery=markevery)
							else:
								plt.plot(T_delrem,rho_delrem,label=label,color=color[i],marker=marker[i],markevery=markevery)
							labeldone=True
						T.append(max(t9))
						rho.append(max(rho1[k]))
					else:
						break
					#else:
						
		plt.xlabel('$T_9 [GK]$',size=22)
		plt.ylabel(r'$\rho [cm^{-3}]$',size=22)
		plt.yscale('log')
		plt.xscale('log')
		plt.legend(loc=2)
	def set_plot_tcrhoc(self,linestyle=['-'],burn_limit=0.997,marker=['o'],markevery=500,end_model=[-1],deg_line=True): 
		'''
			Plots HRDs
			end_model - array, control how far in models a run is plottet, if -1 till end
			symbs_1  - set symbols of runs
		'''
		if len(linestyle)==0:
			linestyle=200*['-']

                for i in range(len(self.runs_H5_surf)):
                        m1p65_last=se(self.runs_H5_out[i])
			t1_model=-1
			if end_model[i] != -1:
				t1_model=end_model[i]

			model=m1p65_last.se.cycles
			model_list=[]
			for k in range(0,len(model),5):
				model_list.append(model[k])
                        rho1=m1p65_last.get(model_list,'rho')   #[:(t1_model-t0_model)]
                        T1=m1p65_last.get(model_list,'temperature')#[:(t1_model-t0_model)]
			T_unit=m1p65_last.get('temperature_unit')
			mass=m1p65_last.get(model_list,'mass')
			mini=m1p65_last.get('mini')
			zini=m1p65_last.get('zini')
			#info=m1p65_last.burnstage_upgrade()
			#burn_info=[]
			#burn_cycle=[]
			#for k in range(len(info[0])):
			#    if 'start' in info[3][k]:
			#	burn_info.append( info[3][k])
			#	burn_cycle.append(info[0][k])
			#print burn_info
			#print burn_cycle
			'''
			#H
			#Maybe use get_elemental_abunds(), too slow if individually taken
			h=m1p65_last.get(model_list,'H-1')
			h_ini=mini*h[0][-1]
			h_limit=h_ini*burn_limit
			print h_limit
			#He
			he=m1p65_last.get(model_list,'He-4')
			he+=m1p65_last.get(model_list,'He-3')
			he_ini=mini*he[0][-1]
			he_limit=he_ini*burn_limit
			print he_limit
			if mini>5:
				c=m1p65_last.get(model_list,'C-12')
				c_ini=mini*c[0][-1]
				c_limit=c_ini*burn_limit
			if mini>10:
				ne=m1p65_last.get(model_list,'Ne-20')
				ne_ini=	mini*ne[0][-1]
				ne_limit=ne_ini*burn_limit
				o=m1p65_last.get(model_list,'O-16')
				o_ini=mini*o[0][-1]
				o_limit=o_ini*burn_limit
				si=m1p65_last.get(model_list,'Si-28')
				si_ini=mini*si[0][-1]
				si_limit=si_ini*burn_limit
			print '#############################3'
			print h_ini,he_ini
			rho=[]
			T=[]
			he_depl_idx=-1
			h_depl_idx=-1
			c_depl_idx=-1
			ne_depl_idx=-1
			o_depl_idx=-1
			si_depl_idx=-1
			'''
			rho=[]
			T=[]
			for k in range(len(model_list)):
				rho_center=rho1[k][0]
				T_center=T1[k][0]
				rho.append(rho_center)
				T.append(T_center)
				'''
				mass11=np.array([0]+list(mass[k]))
				delta_mass=mass11[1:]-mass11[:-1]
				#for the low-mass + massive AGB
				if (sum(np.array(he[k])*np.array(delta_mass))<he_limit) and  (he_depl_idx==-1):
					he_depl_idx=k			
                                if (sum(np.array(h[k])*np.array(delta_mass))<h_limit) and  (h_depl_idx==-1):
                                        h_depl_idx=k
				#for the SAGB + massive
				if mini>5:
					if (sum(np.array(c[k])*np.array(delta_mass))<c_limit) and  (c_depl_idx==-1):
                                        	c_depl_idx=k
				#for the massive stars
				if mini>10:
					if (sum(np.array(ne[k])*np.array(delta_mass))<ne_limit) and  (ne_depl_idx==-1):
                                        	ne_depl_idx=k
					if (sum(np.array(o[k])*np.array(delta_mass))<o_limit) and  (o_depl_idx==-1):
                                        	o_depl_idx=k
					if (sum(np.array(si[k])*np.array(delta_mass))<si_limit) and  (si_depl_idx==-1):
                                        	si_depl_idx=k
				'''

			T=np.log10(np.array(T)*T_unit)
			#T_degeneracy=np.log10(np.array(rho)**(2./3.))
			rho=np.log10(np.array(rho))
                        figure(1)
                        #plot(logTeff,logL,marker=symbs[i],label=self.run_label[i],linestyle=linestyle[i],markevery=markevery)
			#pl.plot(rho,T,marker=symbs[i],label=self.run_label[i],linestyle=linestyle[i],markevery=markevery)
			label=str(mini)+'$M_{\odot}$, Z='+str(zini)
			plt.plot(rho,T,label=label,linestyle=linestyle[i],marker=marker[i],markevery=markevery)
			'''
			#plt.plot(rho,T_degeneracy,color='b',linestyle='-',label='$P_e = P_{e,deg}$')
			print burn_info
			print burn_cycle
			burn_cycle=np.array(burn_cycle)/100.
			print 'changed cycles:', burn_cycle
			print len(rho)
			#For different burning stages:
			if 'H_start' in burn_info:
				plt.plot(rho[burn_cycle[0]],T[burn_cycle[0]],marker='o',color='b')
			if 'He_start' in burn_info:
				plt.plot(rho[burn_cycle[1]],T[burn_cycle[1]],marker='o',color='r')
			if 'C_start' in burn_info:
				plt.plot(rho[burn_cycle[2]],T[burn_cycle[2]],marker='o',color='g')
			if 'Ne_start' in burn_info:
				plt.plot(rho[burn_cycle[3]],T[burn_cycle[3]],marker='D',color='b')
			if 'O_start' in burn_info:
				plt.plot(rho[burn_cycle[4]],T[burn_cycle[4]],marker='D',color='r')
			if 'Si_start' in burn_info:
				plt.plot(rho[burn_cycle[5]],T[burn_cycle[5]],marker='D',color='g')
			'''
		
			#ax = plt.gca()
			#ax.invert_xaxis()
			plt.rcParams.update({'font.size': 16})
			plt.rc('xtick', labelsize=16)
			plt.rc('ytick', labelsize=16)
			legend(loc=4)
		        plt.xlabel('log $\\rho_{\\rm c}$',fontsize=18)
        		plt.ylabel('log $T_{\\rm c}$',fontsize=18)
			'''
			#plt.gca().invert_xaxis()	
			figure(0)
			#pl.plot(rho,T,marker=symbs[i],label=self.run_label[i],linestyle=linestyle[i],markevery=markevery)
                        plt.plot(rho,T,label=self.extra_label[i],linestyle=linestyle[i])
			plt.plot(rho,T_degeneracy,color='b',linestyle='-',label='$P_e = P_{e,deg}$')
			#For different burning stages:
			if not h_depl_idx ==-1:
				plt.plot(rho[h_depl_idx],T[h_depl_idx],marker='o',color='b')
			if not he_depl_idx ==-1:
				plt.plot(rho[he_depl_idx],T[he_depl_idx],marker='o',color='r')
			if not c_depl_idx ==-1:
				plt.plot(rho[c_depl_idx],T[c_depl_idx],marker='o',color='g')
			if not ne_depl_idx ==-1:
				plt.plot(rho[ne_depl_idx],T[ne_depl_idx],marker='D',color='b')
			if not o_depl_idx ==-1:
				plt.plot(rho[o_depl_idx],T[o_depl_idx],marker='D',color='r')
			if not si_depl_idx ==-1:
				plt.plot(rho[si_depl_idx],T[si_depl_idx],marker='D',color='g')
			

			#ax = plt.gca()
			#ax.invert_xaxis()
			plt.rcParams.update({'font.size': 16})
                	plt.rc('xtick', labelsize=16)
                	plt.rc('ytick', labelsize=16)
                	legend(loc=4)
		        plt.xlabel('log $\\rho_{\\rm c}$',fontsize=18)
        		plt.ylabel('log $T_{\\rm c}$',fontsize=18)
			#plt.gca().invert_xaxis()
			'''
			i+=1	
		if deg_line==True:
			rho=np.arange(0,9,0.01)
			T_degeneracy=2./3. *rho +np.log10(1.207e5 * 1.8/(2.**(5./3.)))
			#T_degeneracy=np.log10( (10**np.array(rho))**(2./3.))
			plt.plot(rho,T_degeneracy,color='b',linestyle='-',label='$P_e = P_{e,deg}$')
		plt.legend(loc=2)



	def prod_fac_massgrid_A_1(self,weighted=True,log=True,runs=[],isotopes=[],elements=[],cycles=[],plot_set_diagram=True,color=['r','b','g','k'],marker_type=['o','p','s','D'],line_style=['--','-','-.',':'],markersize=[6,6,6,6],line_width=[14,14,14,14],title='',withlabel=True,label='',plot_lines=True,exp_only=False,pre_exp=False,delay=True,exp_dir='',fontsizelabel='x-small',iniabupath='/astro/critter/critter/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn',xlim=[0,0]):
		'''
			Plots behaviour star mass dependent yields of different isotopes - specify dirs at the beginning
			Beware  of different z,

			runs : If array is empty function uses all available directories.					
				
			!! If len of elements longer than zero, than isotopes will be ignored and elements used!!

		'''
		runs=[]

		HDF5_surf=[]
		HDF5_out=[]
		HDF5_restart=[]
		for i in range(len(self.run_dirs_name)):
			HDF5_surf.append(self.runs_H5_surf[i])
			HDF5_out.append(self.runs_H5_out[i])
			HDF5_restart.append(self.runs_H5_restart[i])
		legend_k=0
		sefiles=[]
		sefiles_hout=[]
		sefiles_restart=[]
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		
				sefiles_hout.append(se(HDF5_out[i]))
				sefiles_restart.append(se(HDF5_restart[i]))


		import utils as u
		iniabu=u.iniabu(iniabupath)
		x_iniabu=[]

		
		isotopes11=sefiles[0].se.isotopes
		if len(isotopes)>0:
			if isotopes[0]=='allstable':
				isotopes=get_stable(isotopes11,get_elements=False)
				color=len(isotopes)*[color[0]]
				marker_type=len(isotopes)*[marker_type[0]]
				line_style=len(isotopes)*[line_style[0]]
				markersize=len(isotopes)*[markersize[0]]
				line_width=len(isotopes)*[line_width[0]]
				print isotopes
			x_iniabu=iniabu.iso_abundance(isotopes)
			mass_numbers=iniabu.a
			charge=iniabu.z

			


		if len(elements)>0:
			if elements[0]=='allstable':
				elements=get_stable(isotopes11,get_elements=True)
				color=len(elements)*[color[0]]
				marker_type=len(elements)*[marker_type[0]]
				line_style=len(elements)*[line_style[0]]
				markersize=len(elements)*[markersize[0]]
				line_width=len(elements)*[line_width[0]]
				print elements
			isotopes=get_stable(isotopes11,get_elements=False)
			import re
			x_iniabu=[0]*len(elements)
			x_iniabu_names=[]
			for k in range(len(iniabu.names)):
				iso=iniabu.names[k].replace(' ','')
				charge1=iniabu.z
				mass_numbers1=iniabu.a
				match = re.match(r"([a-z]+)([0-9]+)",iso, re.I)
				ele1=match.groups()[0].upper()
				iso=ele1+'-'+match.groups()[1]
				for r in range(len(elements)):
					if elements[r] == ele1:
						x_iniabu_names.append(iso)
						x_iniabu[r]+=iniabu.iso_abundance(iso)
			print 'test ouptput',x_iniabu_names	

		#import utils as u

		#iniabu=u.iniabu(iniabupath)

		#iniabu.iso_abundance(isotopes)


		####dealign with explosion pinput
                #test to take right dir:
                #
                if delay:
                        exp_type='delay'
                else:
                        exp_type='rapid'
                slist = os.listdir(exp_dir)
                expr = re.compile(exp_type)
                slist=(filter(expr.search,slist))


		exp_runs_H5_restart=[]
		sefiles_exp=[]


		#to dinstiungish between pre exp and exp sources
		if pre_exp==False:

			slist = os.listdir(exp_dir)
			expr = re.compile(exp_type)
			slist=(filter(expr.search,slist))

			for element in slist:
				run_path=exp_dir+'/'+element
				if not os.path.isdir(run_path):
					continue
				if os.path.isdir(run_path+"/H5_restart"):
					sefiles1 = os.listdir(run_path+"/H5_restart")
					if (filter(expr.search,sefiles1)) <1:
						print "Warning: No hdf5 restart files found in "+run_path+"/H5_restart"
					else:
						exp_runs_H5_restart.append(run_path+"/H5_restart")
						sefiles_exp.append(se(run_path+"/H5_restart"))
		else:
	
			exp_runs_H5_restart=HDF5_restart
			for k in range(len(HDF5_restart)):
				sefiles_exp.append(se(HDF5_restart[k]))


		z_index_files=[]
		z_values=[]
		j=-1
		for i in range(len(HDF5_surf)):
			j+=1
			star_z=sefiles[i].get("zini")
			if star_z not in z_values:
				z_values.append(star_z)
				z_index_files.append([])
			z_index_files[ z_values.index(star_z)].append(i)
		max_yield=[]
		color_iso=-1
		yields_1=[]
		t=0
		if len(elements)>0:
			isotopes=elements
		legend_k=0
		for w in range(len(z_index_files)):
			star_mass_array=[]
			yields=[]
			legend_k+=1
			#iso_yields=[]
			#iniabu_yields_folded=[]
			production_factor=[]
			for i in range(len(isotopes)):
				#iso_yields.append(np.zeros(len( z_index_files[w]   )))
				#iniabu_yields_folded.append(np.zeros(len( z_index_files[w]   )))
				if exp_only==False:	
					production_factor.append(np.zeros(len( z_index_files[w]   )))
				else:
					production_factor.append([])
			ttt=0
			for k in z_index_files[w]:
				if type(sefiles[k].get("mini")) == np.ndarray:
					star_mass=sefiles[k].get("mini")[0]
					star_z=sefiles[k].get("zini")[0]
				else:
                                        star_mass=sefiles[k].get("mini")
                                        star_z=sefiles[k].get("zini")
					
					
				#star_mass=sefiles[k].get("mini")[0]
				#star_z=sefiles[k].get("zini")[0]				
				if cycles[k][1]==-1:
					endcycle=int(sefiles[k].se.cycles[-1]) #+ cycles[k][2]  #1000
					endcycle= int( round( endcycle,-3) -2000 )
				else:
					endcycle=cycles[k][1]
				if exp_only==False:
					star_mass_array.append(star_mass)
					prod_factor,isotopes_prod_fac,yields,iniabu,mass_frac_ini,remn_mass =self.weighted_yields(sefiles[k],sefiles_hout[k],sefiles_restart[k],isotopes,elements,cycles[k][0],endcycle,cycles[k][2])
					#print 'yield output###################################:'
					#print 'wind: ',yields,prod_factor,iniabu
					for tt in range(len(prod_factor)):
                                        	prod_factor[tt]=yields[tt]/(star_mass-remn_mass)    #iniabu[tt]


				else:
					prod_factor=[]
					if star_mass<=8:
						continue
					for pp in range(len(isotopes)): 
						production_factor[pp].append([0])
							
					star_mass_array.append(star_mass)
					
				mass=star_mass
				metallicity=star_z
				if mass>8:

					#in case no ppn_exp run dir is available, skip explosion contribution
					#this is because mass cut is then larger then actual mass of star...no mass lost

					for t in range(len(exp_runs_H5_restart)):
						if type(sefiles_exp[t].get("mini")) == np.ndarray:
							mass1=sefiles_exp[t].get("mini")[0]
							metallicity1=sefiles_exp[t].get("zini")[0]
						else:
							mass1=sefiles_exp[t].get("mini")
							metallicity1=sefiles_exp[t].get("zini")

						###for currupted files

						corrupt=False
						if type(mass1) == list:
							corrupt=True
						if type(metallicity1) ==list:
							corrupt=True
						if ((float(mass1) < 1) or  (float(mass1) >70 )) or type(mass1) == type(None):
							corrupt=True
						if ((float(metallicity1) < 0.00001) or  (float(metallicity1) >1. )) or type(metallicity1) == type(None):
							corrupt=True


						if corrupt == True:
							metallicity1=0.02
							rundir11=exp_runs_H5_restart[t].split('/')[-2]
							mass1=float(rundir11.split('.')[0][1:])
							#print 'corruption, assinge new values',mass1,metallicity1
						if mass1 == mass and metallicity1 == metallicity:
							#sefiles_re2=sefiles_exp[t]
							#calculate exp part
							import nugridse as mp
							reload(mp)
							#Have to find the right cycle in restart file
							#.se.cycles does not work when reading all files	
							#therefore identify file with last cycle and read it:
							refiles=sefiles_exp[t].se.files
							tosort=[]
							refiles1=[]
							for i in range(len(refiles)):
								cycle_file=refiles[i][-18:-11]
								if cycle_file.isdigit():
									tosort.append(int(cycle_file))
									refiles1.append(refiles[i])
							idx_sorted=sorted(range(len(tosort)), key=lambda k: tosort[k])
							idx_file = idx_sorted[-1]
							sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_file])
							if len(sefiles_re_cycle.se.cycles)==0:
								#in the case there is not cycle in the last file...set1.2 M2-
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_sorted[-2]])
							print 'tewststest',exp_cycle_inp
							if exp_cycle_inp>0:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],exp_cycle_inp)
							print 'Use custom exp cycle ', exp_cycle_inp
							#print 'identifies as ',exp_runs_H5_restart[t],'with file',refiles1[idx_file]
							isotopes_exp,yields_exp,iniabu_exp,mass_cut_exp,first_cycle_exp=self.weighted_yields_explosion(mp,sefiles_re_cycle,sefiles[k],isotopes=isotopes,elements=elements,delay=delay)

							if exp_only==False:	
								#add pre-exp yields
								for tt in range(len(prod_factor)):
									#print 'star above 8M'
									#print isotopes_prod_fac[tt],yields[tt],yields_exp[tt]
									#the iniabu from the explosion is not added
									#because it is already in the iniabu from wind
									prod_factor[tt]=(yields[tt]+yields_exp[tt])/(star_mass-mass_cut_exp)     #(iniabu[tt])
									#print 'yields wind:',yields[tt]
									#print 'prodfac wind:',(yields[tt]/iniabu[tt])
									#print 'wind+exp yields : ',yields[tt]+yields_exp[tt]	
							else:
								for tt in range(len(isotopes_exp)):		
									prod_factor.append((yields_exp[tt])/( (star_mass-mass_cut_exp) - (star_mass-remn_mass) ))   #(iniabu_exp[tt]))
							remn_mass=mass_cut_exp
							#print 'exp yield',yields_exp
							#print 'prodfac exp:',(yields_exp[tt]/iniabu_exp[tt])
							#print 'wind+exp prodfac:',prod_factor[tt]
							break		


		

				#print 'iso exp',len(isotopes_exp)
				#prod_factor=np.log10(prod_factor)
				print len(isotopes),len(prod_factor),len(production_factor)
				#Normalization to solar abundance
				for i in range(len(isotopes)):
                                        production_factor[i][ttt]=prod_factor[i]/x_iniabu[i]
				ttt+=1
			
			###plotting
			#if plot_set_diagram==True:
			#if plot_set_diagram==True:
			#	fig_2=plt.figure(isotopes[h])
			#	#fig_2=plt.figure()
			#	ax_2 = fig_2.add_subplot(1,1,1)
			#	plt.rcParams.update({'font.size': 16})
			#	plt.rc('xtick', labelsize=16)
			#	plt.rc('ytick', labelsize=16)
			#	if log==True:
			#		ax_2.set_yscale('log')
                        color_iso=0
                        legend_k=0


			for h in range(len(isotopes)):
				#yield_1=iso_yields[i]
				#mass_1=star_mass_array
				yield_1=[]
				mass_1=[]			
				iniabu_folded=[]
				prod_fac_sorted=[]	
               			indices_array=sorted(range(len(star_mass_array)),key=lambda x:star_mass_array[x])
                		for i in indices_array:
				#	yield_1.append(iso_yields[h][i])
					mass_1.append(star_mass_array[i])
				#	iniabu_folded.append(iniabu_yields_folded[h][i])
					prod_fac_sorted.append(production_factor[h][i])
				####Plotting prodfactor
				#plt.figure(fig_2.number)

				fig_2=plt.figure(isotopes[h])
				ax_2 = fig_2.add_subplot(1,1,1)
				plt.rcParams.update({'font.size': 16})
				plt.rc('xtick', labelsize=16)
				plt.rc('ytick', labelsize=16)
				if log==True:
					ax_2.set_yscale('log')
				ax_2.legend()
				ax_2.set_xlabel("M/M$_{\odot}$",fontsize=16)
				ax_2.minorticks_on()
				ax_2.set_ylabel("Overproduction factor",fontsize=16)
				ax_2.set_title(title)
	
				if len(label)>0:
					label=", "+label
				if len(elements)==0:
					plot_quantity=isotopes[h]
					plot_quantity="$^{"+plot_quantity.split("-")[1]+"}$"+plot_quantity.split("-")[0]
				else:
					plot_quantity=elements[h]
				if withlabel==True:
                                	plt.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso],label=plot_quantity+" , Z="+str(star_z)+label  )				
				else:
                                	plt.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso])				
				

				m_max=max(star_mass_array)+2
				
				if plot_lines==True:
					plt.plot([0,m_max],[1,1],"k--",linewidth=3)
					plt.plot([0,m_max],[2,2],"k--",linewidth=1)
					plt.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)				
	
				color_iso+=1
				legend_k+=1
		#####
		x_imf=[]
		y_imf=[]
		m_max=max(star_mass_array)+2.
		#if plot_set_diagram==True:
		#	plt.figure(fig_2.number)
		##	###plot figure 1
		#	plt.legend()
		##	plt.xlabel("M/M$_{\odot}$",fontsize=20)
		#	plt.minorticks_on()
		#	plt.ylabel("Overproduction factor",fontsize=20)
		
	def prod_fac_massgrid_A(self,fig=0,weighted=True,log=True,runs=[],isotopes=[],elements=[],cycles=[],plot_set_diagram=True,color=['r','b','g','k'],marker_type=['o','p','s','D'],line_style=['--','-','-.',':'],markersize=[6,6,6,6],line_width=[14,14,14,14],title='',withlabel=True,label='',labelzonly=False,plot_lines=True,exp_only=False,pre_exp=False,delay=True,exp_dir='',fontsizelabel='x-small',save_dir='.',arange=[0,0],plot_o16_lines=False,fallback_coords=[],fallback_masses=[],fixed_labels=True):
		'''
			Plots production factor of elements (!!) vs. charge/mass number for different masses.

			runs : If array is empty function uses all available directories.					
				
		'''
		#if len(isotopes)>0 or len(elements)<1:
		#	print 'Specify only elements as input'
		#	return 0
	

                self.fallback_coords=fallback_coords
                self.fallback_masses=fallback_masses
                print 'fallback coords',self.fallback_coords
		
		import utils as u
		
		runs=[]

		HDF5_surf=[]
		HDF5_out=[]
		HDF5_restart=[]
		for i in range(len(self.run_dirs_name)):
			HDF5_surf.append(self.runs_H5_surf[i])
			HDF5_out.append(self.runs_H5_out[i])
			HDF5_restart.append(self.runs_H5_restart[i])
		legend_k=0
		sefiles=[]
		sefiles_hout=[]
		sefiles_restart=[]
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		
				sefiles_hout.append(se(HDF5_out[i]))
				sefiles_restart.append(se(HDF5_restart[i]))
		
		isotopes11=sefiles[0].se.isotopes
		if len(isotopes)>0:
			if isotopes[0]=='allstable':
				isotopes=get_stable(isotopes11,get_elements=False)
				color=len(isotopes)*[color[0]]
				marker_type=len(isotopes)*[marker_type[0]]
				line_style=len(isotopes)*[line_style[0]]
				markersize=len(isotopes)*[markersize[0]]
				line_width=len(isotopes)*[line_width[0]]
				print isotopes
		if len(elements)>0:
			if elements[0]=='allstable':
				elements=get_stable(isotopes11,get_elements=True)
				color=len(elements)*[color[0]]
				marker_type=len(elements)*[marker_type[0]]
				line_style=len(elements)*[line_style[0]]
				markersize=len(elements)*[markersize[0]]
				line_width=len(elements)*[line_width[0]]
				print elements
                a_isotopes=[]
                for k in range(len(isotopes)):
                        a_isotopes.append(int(isotopes[k].split('-')[1]))
	

		####dealign with explosion pinput
                #test to take right dir:
                #
                if delay:
                        exp_type='delay'
                else:
                        exp_type='rapid'
                slist = os.listdir(exp_dir)
                expr = re.compile(exp_type)
                slist=(filter(expr.search,slist))


		exp_runs_H5_restart=[]
		sefiles_exp=[]


		#to dinstiungish between pre exp and exp sources
		if pre_exp==False:

			slist = os.listdir(exp_dir)
			expr = re.compile(exp_type)
			slist=(filter(expr.search,slist))

			for element in slist:
				run_path=exp_dir+'/'+element
				if not os.path.isdir(run_path):
					continue
				if os.path.isdir(run_path+"/H5_restart"):
					sefiles1 = os.listdir(run_path+"/H5_restart")
					if (filter(expr.search,sefiles1)) <1:
						print "Warning: No hdf5 restart files found in "+run_path+"/H5_restart"
					else:
						exp_runs_H5_restart.append(run_path+"/H5_restart")
						sefiles_exp.append(se(run_path+"/H5_restart"))
		else:
	
			exp_runs_H5_restart=HDF5_restart
			for k in range(len(HDF5_restart)):
				sefiles_exp.append(se(HDF5_restart[k]))


		z_index_files=[]
		z_values=[]
		j=-1
		for i in range(len(HDF5_surf)):
			j+=1
			star_z=sefiles[i].get("zini")
			if star_z not in z_values:
				z_values.append(star_z)
				z_index_files.append([])
			z_index_files[ z_values.index(star_z)].append(i)
		max_yield=[]
		color_iso=-1
		yields_1=[]
		t=0
		if len(elements)>0:
			isotopes=elements
		legend_k=0
		for w in range(len(z_index_files)):
			star_mass_array=[]
			yields=[]
			legend_k+=1
			#iso_yields=[]
			#iniabu_yields_folded=[]
			production_factor=[]
			for i in range(len(isotopes)):
				#iso_yields.append(np.zeros(len( z_index_files[w]   )))
				#iniabu_yields_folded.append(np.zeros(len( z_index_files[w]   )))
				if exp_only==False:	
					production_factor.append(np.zeros(len( z_index_files[w]   )))
				else:
					production_factor.append([])
			ttt=0
			for k in z_index_files[w]:
				if type(sefiles[k].get("mini")) == np.ndarray:
					star_mass=sefiles[k].get("mini")[0]
					star_z=sefiles[k].get("zini")[0]
				else:
                                        star_mass=sefiles[k].get("mini")
                                        star_z=sefiles[k].get("zini")
				print 'mini',star_mass
					
				#star_mass=sefiles[k].get("mini")[0]
				#star_z=sefiles[k].get("zini")[0]


                                if cycles[k][1]==-1:
					print 'cycles'
					#print sefiles_restart[k].se.cycles
                                        endcycle=int(sefiles_restart[k].se.cycles[-2])
                                else:
                                        endcycle=cycles[k][1]
                                cycs=range(cycles[k][0],endcycle,cycles[k][2])
                                if endcycle not in cycs:
                                        cycs=cycs+[endcycle]
                                star_mass11=sefiles[k].get(cycs,'mass')
                                if not len(star_mass11)==len(cycs):
                                        cycs1=[int(cycs[0])]
                                        for cyc in cycs[1:]:
                                                if  cyc <= max(cycs1):
                                                        continue
                                                a=sefiles[k].get(cyc,'mass')
                                                cyctest=int(cyc)
                                                #if cycle 0 or a non-existing cycle, or cycle was alrady added, 
                                                #go one cycle further   
                                                while ((type(a)==list) or (a==0)):
                                                        if   cyctest<int(cycs[-1]):
                                                                print int(cycs[-1])
                                                                print 'found bad cycle',cyctest
                                                                #if last cycle does not exist, go backwards
                                                                #if cyctest == cycs[-1]:
                                                                #       cyctest-=1
                                                                #else:
                                                                cyctest+=1
                                                                a=sefiles[k].get(cyctest,'mass')
                                                        elif cyctest==int(cycs[-1]):
                                                                cyctest=int(sefiles_restart[k].se.cycles[-3])
                                                                a=sefiles[k].get(cyctest,'mass')
                                                                h=-4
                                                                while ((type(a)==list) or (a==0)):
                                                                        cyctest=int(sefiles_restart[k].se.cycles[h])
                                                                        a=sefiles[k].get(cyctest,'mass')
                                                                        h-=1
                                                                for kk in cycs1:
                                                                        if kk>=cyctest:
                                                                                cycs1=cycs1[:kk]
                                                                break
                                                        else:
                                                                break
                                                print 'Use cycle ',cyctest
                                                if cyctest>max(cycs1):
                                                        cycs1.append(cyctest)
                                        print len(cycs),len(cycs1)
                                        cycs=cycs1
                                print cycs[-5:]
				if exp_only==False:
					star_mass_array.append(star_mass)
					prod_factor,isotopes_prod_fac,yields,iniabu,mass_frac_ini,remn_mass = self.weighted_yields(sefiles[k],sefiles_restart[k],isotopes,elements,cycs)

					#print 'yield output###################################:'
					#print 'wind: ',yields,prod_factor,iniabu
					for tt in range(len(prod_factor)):
                                        	prod_factor[tt]=yields[tt]/iniabu[tt]


				else:
					prod_factor=[]
					if star_mass<=8:
						continue
					for pp in range(len(isotopes)): 
						production_factor[pp].append([0])
							
					star_mass_array.append(star_mass)
					
				mass=star_mass
				metallicity=star_z
				if mass>8:

					#in case no ppn_exp run dir is available, skip explosion contribution
					#this is because mass cut is then larger then actual mass of star...no mass lost

					for t in range(len(exp_runs_H5_restart)):
						if type(sefiles_exp[t].get("mini")) == np.ndarray:
							mass1=sefiles_exp[t].get("mini")[0]
							metallicity1=sefiles_exp[t].get("zini")[0]
						else:
							mass1=sefiles_exp[t].get("mini")
							metallicity1=sefiles_exp[t].get("zini")

						###for currupted files

						corrupt=False
						if type(mass1) == list:
							corrupt=True
						if type(metallicity1) ==list:
							corrupt=True
						if ((float(mass1) < 1) or  (float(mass1) >70 )) or type(mass1) == type(None):
							corrupt=True
						if ((float(metallicity1) < 0.00001) or  (float(metallicity1) >1. )) or type(metallicity1) == type(None):
							corrupt=True


						if corrupt == True:
							metallicity1=0.02
							rundir11=exp_runs_H5_restart[t].split('/')[-2]
							mass1=float(rundir11.split('.')[0][1:])
							#print 'corruption, assinge new values',mass1,metallicity1
						if mass1 == mass and metallicity1 == metallicity:
							#sefiles_re2=sefiles_exp[t]
							#calculate exp part
							import nugridse as mp
							reload(mp)
							#Have to find the right cycle in restart file
							#.se.cycles does not work when reading all files	
							#therefore identify file with last cycle and read it:
							refiles=sefiles_exp[t].se.files
							tosort=[]
							refiles1=[]
							for i in range(len(refiles)):
								cycle_file=refiles[i][-18:-11]
								if cycle_file.isdigit():
									tosort.append(int(cycle_file))
									refiles1.append(refiles[i])
							idx_sorted=sorted(range(len(tosort)), key=lambda k: tosort[k])
							idx_file = idx_sorted[-1]
							if pre_exp==True:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t])
							else:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_file])
							if len(sefiles_re_cycle.se.cycles)==0:
								#in the case there is not cycle in the last file...set1.2 M2-
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_sorted[-2]])
                                                        if pre_exp==True:
                                                                cycleend=cycs[-1]
                                                        else:
                                                                cycleend=sefiles_re_cycle.se.cycles[-1]
							#print 'identifies as ',exp_runs_H5_restart[t],'with file',refiles1[idx_file]
							###prodfac A
							isotopes_exp,yields_exp,iniabu_exp,mass_cut_exp,first_cycle_exp= self.weighted_yields_explosion(mp,sefiles_re_cycle,sefiles[k],isotopes=isotopes,elements=elements,cycleend=cycleend,delay=delay)


							if exp_only==False:	
								#add pre-exp yields
								for tt in range(len(prod_factor)):
									#print 'star above 8M'
									#print isotopes_prod_fac[tt],yields[tt],yields_exp[tt]
									#the iniabu from the explosion is not added
									#because it is already in the iniabu from wind
									prod_factor[tt]=(yields[tt]+yields_exp[tt])/(iniabu[tt] + iniabu_exp[tt])#(iniabu[tt])
									#print 'yields wind:',yields[tt]
									#print 'prodfac wind:',(yields[tt]/iniabu[tt])
									#print 'wind+exp yields : ',yields[tt]+yields_exp[tt]	
							else:
								for tt in range(len(isotopes_exp)):		
									prod_factor.append((yields_exp[tt])/(iniabu_exp[tt]))
							remn_mass=mass_cut_exp
							#print 'exp yield',yields_exp
							#print 'prodfac exp:',(yields_exp[tt]/iniabu_exp[tt])
							#print 'wind+exp prodfac:',prod_factor[tt]
							break		
							
				#plot for every star prodfac vs A.
                                if fig ==0:
                                        figname=str(mass)+'_'+str(star_z)
                                else:
                                        figname=fig
                                plt.figure(figname)
				if labelzonly==True:
					label='Z='+str(star_z)
				else:
					label=str(mass)+' M$_{\odot}$'
				#plt.plot(a_isotopes,prod_factor,label=label,marker='o',color='k')
				plt.xlabel('Mass number A')
				plt.ylabel('Overproduction factor')
				ax=plt.gca()
				label_done=[]
				idx_min=0
				hh=0
                                idx_sorted=sorted(range(len(a_isotopes)),key=lambda x:a_isotopes[x])
                                a_isotopes_sorted=[]
				isotopes_sorted=[]
                                prod_factor_sorted=[]

				isotopes_sorted_group=[]
				prod_factor_sorted_group=[]
				element_group=[]
				a_isotopes_sorted_group=[]

                                for idxx in idx_sorted:
                                        a_isotopes_sorted.append(a_isotopes[idxx])
                                        prod_factor_sorted.append(prod_factor[idxx])
					isotopes_sorted.append(isotopes[idxx])
					ele1=isotopes[idxx].split('-')[0]		
					if ele1 in element_group:
						idx1=element_group.index(ele1)
						prod_factor_sorted_group[idx1].append(prod_factor[idxx])
						a_isotopes_sorted_group[idx1].append(a_isotopes[idxx])
						isotopes_sorted_group[idx1].append(isotopes[idxx])
					else:
						element_group.append(ele1)
						prod_factor_sorted_group.append([])
						a_isotopes_sorted_group.append([])
						isotopes_sorted_group.append([])
						prod_factor_sorted_group[-1].append(prod_factor[idxx])
						a_isotopes_sorted_group[-1].append(a_isotopes[idxx])
						isotopes_sorted_group[-1].append(isotopes[idxx])
				marker='o'
				color1='k'
				linestyle='-'
				if len(marker_type)>0:
					marker=marker_type[k]
				if len(color)>0:
					color1=color[k]
				if len(line_style)>0:
					linestyle=line_style[k]


				for z in range(len(element_group)):
					if (z == (len(element_group)-1)):
						plt.plot(a_isotopes_sorted_group[z],prod_factor_sorted_group[z],marker=marker,color=color1,linestyle=linestyle,label=label)
					else:	
                                        	plt.plot(a_isotopes_sorted_group[z],prod_factor_sorted_group[z],marker=marker,color=color1,linestyle=linestyle)
					print isotopes_sorted_group[z] 
					print a_isotopes_sorted_group[z]

				for z in range(len(prod_factor)):
					if isotopes[z] == 'O-16':
						if plot_o16_lines==True:
							plt.plot([0,208],[prod_factor[z],prod_factor[z]],"k-",linewidth=3)
							plt.plot([0,208],[prod_factor[z]*2,prod_factor[z]*2],"k--",linewidth=1)
							plt.plot([0,208],[prod_factor[z]/2.,prod_factor[z]/2.],"k--",linewidth=1)
							ax.text(208,prod_factor[z]*1.1, 'O-16',horizontalalignment='right',verticalalignment='bottom',clip_on=True,fontsize=9)
							ax.text(208, prod_factor[z]*2*1.1, 'O-16*2',horizontalalignment='right',verticalalignment='bottom',clip_on=True,fontsize=9)
							ax.text(208, (prod_factor[z]/2)*1.05, 'O-16/2',horizontalalignment='right',verticalalignment='bottom',clip_on=True,fontsize=9)

					#for last isotope, if its element appears (not in label_done) with it
					#if (z == (len(prod_factor)-1)) and (not (isotopes_sorted[z].split('-')[0]  in label_done)):
						#plt.plot(a_isotopes_sorted[z],prod_factor_sorted[z],label=label,marker=marker,color=color1,linestyle=linestyle)
					#if its not the first time the element appears but reaches end of array
					elif (z == (len(prod_factor)-1)):
						idx_max=z
                                        	#plt.plot(a_isotopes_sorted[idx_min:idx_max],prod_factor_sorted[idx_min:idx_max],label=label,marker=marker,color=color1,linestyle=linestyle)
						break

					if (isotopes_sorted[z].split('-')[0]  in label_done):			
						continue
					idx_max=z
					label_done.append(isotopes_sorted[z].split('-')[0])
					#plt.plot(a_isotopes_sorted[idx_min:idx_max],prod_factor_sorted[idx_min:idx_max],marker=marker,color=color1,linestyle=linestyle)
					idx_min=z
					if fixed_labels==True:
						if hh==0:
							nmin=0.1
						if hh==1:
							nmin=0.15
						if hh==2:
							nmin=0.2
						if hh==3:
							nmin=0.3
					else:
						nmin=prod_factor_sorted[z]*1.3
					#label below belongs to new element
					print a_isotopes_sorted[z],prod_factor_sorted[z]
					ax.text(a_isotopes_sorted[z],nmin,isotopes_sorted[z],horizontalalignment='center',verticalalignment='center',\
                              		fontsize=fontsizelabel,clip_on=True) #x-msall




                                        if hh==0:
                                                hh=1
                                                continue
                                        if hh==1:
                                                hh=2
                                                continue
                                        if hh==2:
                                                hh=3
                                                continue
					if hh==3:
						hh=0
						continue
					
				m_max=208
                                if plot_lines==True:
                                        plt.plot([0,m_max],[1,1],"k-",linewidth=3)
                                        plt.plot([0,m_max],[2,2],"k--",linewidth=1)
                                        plt.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)
				#set the plot x axis range
				if arange[1]>0:
					plt.xlim(arange[0],arange[1])
				else:
					plt.xlim(0,208)



				#plt.ylim(ymax=(max(prod_factor)+0.1))
				plt.yscale('log')
                                plt.minorticks_on()
                                plt.ylabel("Overproduction factor",fontsize=16)
				if len(save_dir)>0:
					plt.savefig(save_dir+'/'+figname+'.prodfac_vsZ.png')
				#print 'iso exp',len(isotopes_exp)
				#prod_factor=np.log10(prod_factor)
				print len(isotopes),len(prod_factor),len(production_factor)
				for i in range(len(isotopes)):
                                        production_factor[i][ttt]=prod_factor[i]
				ttt+=1
			
			###plotting
			#if plot_set_diagram==True:
			#if plot_set_diagram==True:
			#	fig_2=plt.figure(isotopes[h])
			#	#fig_2=plt.figure()
			#	ax_2 = fig_2.add_subplot(1,1,1)
			#	plt.rcParams.update({'font.size': 16})
			#	plt.rc('xtick', labelsize=16)
			#	plt.rc('ytick', labelsize=16)
			#	if log==True:
			#		ax_2.set_yscale('log')
                        color_iso=0
                        legend_k=0

			#here loop over star masses needed
			break
			for h in range(len(isotopes)):
				#yield_1=iso_yields[i]
				#mass_1=star_mass_array
				yield_1=[]
				mass_1=[]			
				iniabu_folded=[]
				prod_fac_sorted=[]	
               			indices_array=sorted(range(len(star_mass_array)),key=lambda x:star_mass_array[x])
                		for i in indices_array:
				#	yield_1.append(iso_yields[h][i])
					mass_1.append(star_mass_array[i])
				#	iniabu_folded.append(iniabu_yields_folded[h][i])
					prod_fac_sorted.append(production_factor[h][i])
				####Plotting prodfactor
				#plt.figure(fig_2.number)

				fig_2=plt.figure(isotopes[h])
				ax_2 = fig_2.add_subplot(1,1,1)
				plt.rcParams.update({'font.size': 16})
				plt.rc('xtick', labelsize=16)
				plt.rc('ytick', labelsize=16)
				if log==True:
					ax_2.set_yscale('log')
				ax_2.legend()
				ax_2.set_xlabel("M/M$_{\odot}$",fontsize=16)
				ax_2.minorticks_on()
				ax_2.set_ylabel("Overproduction factor",fontsize=16)
				ax_2.set_title(title)
	
				if len(label)>0:
					label=", "+label
				if len(elements)==0:
					plot_quantity=isotopes[h]
					plot_quantity="$^{"+plot_quantity.split("-")[1]+"}$"+plot_quantity.split("-")[0]
				else:
					plot_quantity=elements[h]
				if withlabel==True:
                                	plt.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso],label=plot_quantity+" , Z="+str(star_z)+label  )				
				else:
                                	plt.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso])				
				

				m_max=max(star_mass_array)+2
				
				if plot_lines==True:
					plt.plot([0,m_max],[1,1],"k--",linewidth=3)
					plt.plot([0,m_max],[2,2],"k--",linewidth=1)
					plt.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)
									
	
				color_iso+=1
				legend_k+=1
		#####
		x_imf=[]
		y_imf=[]
		m_max=max(star_mass_array)+2.
		#if plot_set_diagram==True:
		#	plt.figure(fig_2.number)
		##	###plot figure 1
		#	plt.legend()
		##	plt.xlabel("M/M$_{\odot}$",fontsize=20)
		#	plt.minorticks_on()
		#	plt.ylabel("Overproduction factor",fontsize=20)
		#	plt.title(title)






	def prod_fac_massgrid_Z(self,fig=0,weighted=True,log=True,runs=[],isotopes=[],elements=[],cycles=[],plot_set_diagram=True,color=['r','b','g','k'],marker_type=['o','p','s','D'],line_style=['--','-','-.',':'],markersize=[6,6,6,6],line_width=[14,14,14,14],title='',withlabel=True,label='',labelZonly=False,plot_lines=True,exp_only=False,pre_exp=False,delay=True,exp_dir='',fontsizelabel='x-small',save_dir='.',arange=[0,0],yields_output=False,fallback_coords=[],fallback_masses=[]):
		'''
			Plots production factor of elements (!!) vs. charge/mass number for different masses.

			runs : If array is empty function uses all available directories.					
				
		'''

                self.fallback_coords=fallback_coords
                self.fallback_masses=fallback_masses
                print 'fallback coords',self.fallback_coords

		#if len(isotopes)>0 or len(elements)<1:
		#	print 'Specify only elements as input'
		#	return 0
			
		import utils as u
		
		runs=[]

		HDF5_surf=[]
		HDF5_out=[]
		HDF5_restart=[]
		for i in range(len(self.run_dirs_name)):
			HDF5_surf.append(self.runs_H5_surf[i])
			HDF5_out.append(self.runs_H5_out[i])
			HDF5_restart.append(self.runs_H5_restart[i])
		legend_k=0
		sefiles=[]
		sefiles_hout=[]
		sefiles_restart=[]
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		
				sefiles_hout.append(se(HDF5_out[i]))
				sefiles_restart.append(se(HDF5_restart[i]))
		
		isotopes11=sefiles[0].se.isotopes
		if len(isotopes)>0:
			if isotopes[0]=='allstable':
				isotopes=get_stable(isotopes11,get_elements=False)
				#color=len(isotopes)*[color[0]]
				#marker_type=len(isotopes)*[marker_type[0]]
				#line_style=len(isotopes)*[line_style[0]]
				#markersize=len(isotopes)*[markersize[0]]
				#line_width=len(isotopes)*[line_width[0]]
				print isotopes
		if len(elements)>0:
			if elements[0]=='allstable':
				elements=get_stable(isotopes11,get_elements=True)
				#color=len(elements)*[color[0]]
				#marker_type=len(elements)*[marker_type[0]]
				#line_style=len(elements)*[line_style[0]]
				#markersize=len(elements)*[markersize[0]]
				#line_width=len(elements)*[line_width[0]]
				print elements
                z_elements=[]
                for k in range(len(elements)):
                        z_elements.append(u.get_z_from_el(elements[k]))
	

		####dealign with explosion pinput
                #test to take right dir:
                #
                if delay:
                        exp_type='delay'
                else:
                        exp_type='rapid'
                slist = os.listdir(exp_dir)
                expr = re.compile(exp_type)
                slist=(filter(expr.search,slist))


		exp_runs_H5_restart=[]
		sefiles_exp=[]


		#to dinstiungish between pre exp and exp sources
		if pre_exp==False:

			slist = os.listdir(exp_dir)
			expr = re.compile(exp_type)
			slist=(filter(expr.search,slist))

			for element in slist:
				run_path=exp_dir+'/'+element
				if not os.path.isdir(run_path):
					continue
				if os.path.isdir(run_path+"/H5_restart"):
					sefiles1 = os.listdir(run_path+"/H5_restart")
					if (filter(expr.search,sefiles1)) <1:
						print "Warning: No hdf5 restart files found in "+run_path+"/H5_restart"
					else:
						exp_runs_H5_restart.append(run_path+"/H5_restart")
						sefiles_exp.append(se(run_path+"/H5_restart"))
		else:
	
			exp_runs_H5_restart=HDF5_restart
			for k in range(len(HDF5_restart)):
				sefiles_exp.append(se(HDF5_restart[k]))


		z_index_files=[]
		z_values=[]
		j=-1
		for i in range(len(HDF5_surf)):
			j+=1
			star_z=sefiles[i].get("zini")
			if star_z not in z_values:
				z_values.append(star_z)
				z_index_files.append([])
			z_index_files[ z_values.index(star_z)].append(i)
		max_yield=[]
		color_iso=-1
		yields_1=[]
		t=0
		if len(elements)>0:
			isotopes=elements
		legend_k=0
		for w in range(len(z_index_files)):
			star_mass_array=[]
			yields=[]
			legend_k+=1
			#iso_yields=[]
			#iniabu_yields_folded=[]
			production_factor=[]
			for i in range(len(isotopes)):
				#iso_yields.append(np.zeros(len( z_index_files[w]   )))
				#iniabu_yields_folded.append(np.zeros(len( z_index_files[w]   )))
				if exp_only==False:	
					production_factor.append(np.zeros(len( z_index_files[w]   )))
				else:
					production_factor.append([])
			ttt=0
			for k in z_index_files[w]:
				if type(sefiles[k].get("mini")) == np.ndarray:
					star_mass=sefiles[k].get("mini")[0]
					star_z=sefiles[k].get("zini")[0]
				else:
                                        star_mass=sefiles[k].get("mini")
                                        star_z=sefiles[k].get("zini")
				print 'mini',star_mass
					
				#star_mass=sefiles[k].get("mini")[0]
				#star_z=sefiles[k].get("zini")[0]


                                if cycles[k][1]==-1:
					print 'cycles'
					#print sefiles_restart[k].se.cycles
                                        endcycle=int(sefiles_restart[k].se.cycles[-2])
                                else:
                                        endcycle=cycles[k][1]
                                cycs=range(cycles[k][0],endcycle,cycles[k][2])
                                if endcycle not in cycs:
                                        cycs=cycs+[endcycle]
                                star_mass11=sefiles[k].get(cycs,'mass')
                                if not len(star_mass11)==len(cycs):
                                        cycs1=[int(cycs[0])]
                                        for cyc in cycs[1:]:
                                                if  cyc <= max(cycs1):
                                                        continue
                                                a=sefiles[k].get(cyc,'mass')
                                                cyctest=int(cyc)
                                                #if cycle 0 or a non-existing cycle, or cycle was alrady added, 
                                                #go one cycle further   
                                                while ((type(a)==list) or (a==0)):
                                                        if   cyctest<int(cycs[-1]):
                                                                print int(cycs[-1])
                                                                print 'found bad cycle',cyctest
                                                                #if last cycle does not exist, go backwards
                                                                #if cyctest == cycs[-1]:
                                                                #       cyctest-=1
                                                                #else:
                                                                cyctest+=1
                                                                a=sefiles[k].get(cyctest,'mass')
                                                        elif cyctest==int(cycs[-1]):
                                                                cyctest=int(sefiles_restart[k].se.cycles[-3])
                                                                a=sefiles[k].get(cyctest,'mass')
                                                                h=-4
                                                                while ((type(a)==list) or (a==0)):
                                                                        cyctest=int(sefiles_restart[k].se.cycles[h])
                                                                        a=sefiles[k].get(cyctest,'mass')
                                                                        h-=1
                                                                for kk in cycs1:
                                                                        if kk>=cyctest:
                                                                                cycs1=cycs1[:kk]
                                                                break
                                                        else:
                                                                break
                                                print 'Use cycle ',cyctest
                                                if cyctest>max(cycs1):
                                                        cycs1.append(cyctest)
                                        print len(cycs),len(cycs1)
                                        cycs=cycs1
                                print cycs[-5:]
				if exp_only==False:
					star_mass_array.append(star_mass)
					prod_factor,isotopes_prod_fac,yields,iniabu,mass_frac_ini,remn_mass = self.weighted_yields(sefiles[k],sefiles_restart[k],isotopes,elements,cycs)
					#print 'yield output###################################:'
					print 'wind: ',yields,prod_factor,iniabu,isotopes_prod_fac,isotopes
					for tt in range(len(prod_factor)):
                                                if yields_output==False:
                                                        prod_factor[tt]=yields[tt]/iniabu[tt]
                                                else:
                                                        prod_factor[tt]=yields[tt]


				else:
					prod_factor=[]
					if star_mass<=8:
						continue
					for pp in range(len(isotopes)): 
						production_factor[pp].append([0])
							
					star_mass_array.append(star_mass)
					
				mass=star_mass
				metallicity=star_z
				if mass>8:

					#in case no ppn_exp run dir is available, skip explosion contribution
					#this is because mass cut is then larger then actual mass of star...no mass lost

					for t in range(len(exp_runs_H5_restart)):
						if type(sefiles_exp[t].get("mini")) == np.ndarray:
							mass1=sefiles_exp[t].get("mini")[0]
							metallicity1=sefiles_exp[t].get("zini")[0]
						else:
							mass1=sefiles_exp[t].get("mini")
							metallicity1=sefiles_exp[t].get("zini")

						###for currupted files

						corrupt=False
						if type(mass1) == list:
							corrupt=True
						if type(metallicity1) ==list:
							corrupt=True
						if ((float(mass1) < 1) or  (float(mass1) >70 )) or type(mass1) == type(None):
							corrupt=True
						if ((float(metallicity1) < 0.00001) or  (float(metallicity1) >1. )) or type(metallicity1) == type(None):
							corrupt=True


						if corrupt == True:
							metallicity1=0.02
							rundir11=exp_runs_H5_restart[t].split('/')[-2]
							mass1=float(rundir11.split('.')[0][1:])
							#print 'corruption, assinge new values',mass1,metallicity1
						if mass1 == mass and metallicity1 == metallicity:
							#sefiles_re2=sefiles_exp[t]
							#calculate exp part
							import nugridse as mp
							reload(mp)
							#Have to find the right cycle in restart file
							#.se.cycles does not work when reading all files	
							#therefore identify file with last cycle and read it:
							refiles=sefiles_exp[t].se.files
							tosort=[]
							refiles1=[]
							for i in range(len(refiles)):
								cycle_file=refiles[i][-18:-11]
								if cycle_file.isdigit():
									tosort.append(int(cycle_file))
									refiles1.append(refiles[i])
							idx_sorted=sorted(range(len(tosort)), key=lambda k: tosort[k])
							idx_file = idx_sorted[-1]
							if pre_exp==True:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t])
							else:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_file])
							if len(sefiles_re_cycle.se.cycles)==0:
								#in the case there is not cycle in the last file...set1.2 M2-
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_sorted[-2]])
                                                        if pre_exp==True:
                                                                cycleend=cycs[-1]
                                                        else:
                                                                cycleend=sefiles_re_cycle.se.cycles[-1]
							#print 'identifies as ',exp_runs_H5_restart[t],'with file',refiles1[idx_file]
							isotopes_exp,yields_exp,iniabu_exp,mass_cut_exp,first_cycle_exp= self.weighted_yields_explosion(mp,sefiles_re_cycle,sefiles[k],isotopes=isotopes,elements=elements,cycleend=cycleend,delay=delay)


							if exp_only==False:	
								#add pre-exp yields
								for tt in range(len(prod_factor)):
									#print 'star above 8M'
									#print isotopes_prod_fac[tt],yields[tt],yields_exp[tt]
									#the iniabu from the explosion is not added
									#because it is already in the iniabu from wind

                                                                        if yields_output==False:
                                                                                prod_factor[tt]=(yields[tt]+yields_exp[tt])/(iniabu[tt] + iniabu_exp[tt])#(iniabu[tt])
                                                                        else:
                                                                                prod_factor[tt]=yields[tt]+yields_exp[tt]

									#print 'yields wind:',yields[tt]
									#print 'prodfac wind:',(yields[tt]/iniabu[tt])
									#print 'wind+exp yields : ',yields[tt]+yields_exp[tt]	
							else:
								for tt in range(len(isotopes_exp)):		
                                                                        if yields_output==False:
                                                                                prod_factor.append((yields_exp[tt])/(iniabu_exp[tt]))
                                                                        else:
                                                                                prod_factor.append(yields_exp[tt])
							remn_mass=mass_cut_exp
							#print 'exp yield',yields_exp
							#print 'prodfac exp:',(yields_exp[tt]/iniabu_exp[tt])
							#print 'wind+exp prodfac:',prod_factor[tt]
							break		
							
				#plot for every star prodfac vs A.
				if fig ==0:
					figname=str(mass)+'_'+str(star_z)
				else:
					figname=fig
				print 'test output##########'
				print prod_factor
				print elements
				plt.figure(figname)
				label=str(mass)+' M$_{\odot}$'
				if labelZonly==True:
					label='Z='+str(star_z)
				idx_sorted=sorted(range(len(z_elements)),key=lambda x:z_elements[x])
				z_elements_sorted=[]
				prod_factor_sorted=[]
				elements_sorted=[]
				for idxx in idx_sorted:
					z_elements_sorted.append(z_elements[idxx])
					prod_factor_sorted.append(prod_factor[idxx])
					elements_sorted.append(elements[idxx])
				marker='o'
				color1='k'
				linestyle='-'
				if len(marker_type)>0:
					marker=marker_type[k]
				if len(color)>0:
					color1=color[k]
				if len(line_style)>0:
					linestyle=line_style[k]
				plt.plot(z_elements_sorted,prod_factor_sorted,label=label,marker=marker,color=color1,linestyle=linestyle)
				plt.xlabel('Charge Z')
				plt.ylabel('Produ..')
				ax=plt.gca()
				hh=0
				import matplotlib.transforms as transforms
				trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
				print z_elements_sorted
				print prod_factor_sorted
				print elements_sorted
				for z in range(len(prod_factor)):
					pos=0.90
					if hh==0:
						nmin=pos#max(prod_factor)#0.2
					if hh==1:
						nmin=pos+0.02#max(prod_factor)+0.05#0.25
					if hh==2:
						nmin=pos+0.04#max(prod_factor)+#0.3
					fontsizelabel='x-small'
					ax.text(z_elements_sorted[z],nmin,elements_sorted[z],horizontalalignment='center',verticalalignment='center',\
                              		fontsize=fontsizelabel,clip_on=True,transform=trans) #x-msall
					if hh==0:
						hh=1
						continue
					if hh==1:
						hh=2
						continue
					if hh==2:
						hh=0
						continue
				m_max=86
                                if plot_lines==True:
                                        plt.plot([0,m_max],[1,1],"k-",linewidth=3)
                                        #plt.plot([0,m_max],[2,2],"k--",linewidth=1)
                                        #plt.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)
				
                               #set the plot x axis range
                                if arange[1]>0:
                                        plt.xlim(arange[0],arange[1])
                                else:
					plt.xlim(0,84)
				plt.ylim(1,(max(prod_factor)+0.1))
				plt.yscale('log')
                                plt.minorticks_on()
                                if yields_output==False:
                                        plt.ylabel("Overproduction factor",fontsize=16)
                                else:
                                        plt.ylabel("Yields/M$_{\odot}$",fontsize=16)
				plt.legend()
				#print 'iso exp',len(isotopes_exp)
				#prod_factor=np.log10(prod_factor)
				if yields_output==False:
					plt.savefig(save_dir+'/'+figname+'.prodfac_vsZ.png')
				else:
					plt.savefig(save_dir+'/'+figname+'.yields_vsZ.png')
				print len(isotopes),len(prod_factor),len(production_factor)
				for i in range(len(isotopes)):
                                        production_factor[i][ttt]=prod_factor[i]
				ttt+=1
			
			###plotting
			#if plot_set_diagram==True:
			#if plot_set_diagram==True:
			#	fig_2=plt.figure(isotopes[h])
			#	#fig_2=plt.figure()
			#	ax_2 = fig_2.add_subplot(1,1,1)
			#	plt.rcParams.update({'font.size': 16})
			#	plt.rc('xtick', labelsize=16)
			#	plt.rc('ytick', labelsize=16)
			#	if log==True:
			#		ax_2.set_yscale('log')
                        color_iso=0
                        legend_k=0

			#here loop over star masses needed
			break
			for h in range(len(isotopes)):
				#yield_1=iso_yields[i]
				#mass_1=star_mass_array
				yield_1=[]
				mass_1=[]			
				iniabu_folded=[]
				prod_fac_sorted=[]	
               			indices_array=sorted(range(len(star_mass_array)),key=lambda x:star_mass_array[x])
                		for i in indices_array:
				#	yield_1.append(iso_yields[h][i])
					mass_1.append(star_mass_array[i])
				#	iniabu_folded.append(iniabu_yields_folded[h][i])
					prod_fac_sorted.append(production_factor[h][i])
				####Plotting prodfactor
				#plt.figure(fig_2.number)




				fig_2=plt.figure(isotopes[h])
				ax_2 = fig_2.add_subplot(1,1,1)
				plt.rcParams.update({'font.size': 16})
				plt.rc('xtick', labelsize=16)
				plt.rc('ytick', labelsize=16)
				if log==True:
					ax_2.set_yscale('log')
				ax_2.legend()
				ax_2.set_xlabel("M/M$_{\odot}$",fontsize=16)
				ax_2.minorticks_on()
                                if yields_output==False:
                                        ax_2.set_ylabel("Overproduction factor",fontsize=16)
                                else:
                                        ax_2.set_ylabel("Yields/M$_{\odot}$",fontsize=16)
				ax_2.set_title(title)
	
				if len(label)>0:
					label=", "+label
				if len(elements)==0:
					plot_quantity=isotopes[h]
					plot_quantity="$^{"+plot_quantity.split("-")[1]+"}$"+plot_quantity.split("-")[0]
				else:
					plot_quantity=elements[h]
				if withlabel==True:
                                	plt.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso],label=plot_quantity+" , Z="+str(star_z)+label  )				
				else:
                                	plt.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso])				
				

				m_max=max(star_mass_array)+2
				
				if plot_lines==True:
					plt.plot([0,m_max],[1,1],"k--",linewidth=3)
					plt.plot([0,m_max],[2,2],"k--",linewidth=1)
					plt.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)				
	
				color_iso+=1
				legend_k+=1
		#####
		x_imf=[]
		y_imf=[]
		m_max=max(star_mass_array)+2.
		#if plot_set_diagram==True:
		#	plt.figure(fig_2.number)
		##	###plot figure 1
		#	plt.legend()
		##	plt.xlabel("M/M$_{\odot}$",fontsize=20)
		#	plt.minorticks_on()
		#	plt.ylabel("Overproduction factor",fontsize=20)
		#	plt.title(title)

	def prod_fac_massgrid_metallicity(self,weighted=True,log=True,runs=[],isotopes=[],elements=[],cycles=[],plot_set_diagram=True,color=['r','b','g','k'],marker_type=['o','p','s','D'],line_style=['--','-','-.',':'],markersize=[6,6,6,6],line_width=[14,14,14,14],title='',withlabel=True,withlabelelement=False,label='',plot_lines=True,exp_only=False,pre_exp=False,delay=True,exp_dir='',yields_output=False):
		'''
			Plots behaviour star mass dependent yields of different isotopes - specify dirs at the beginning
			Beware  of different z,

			runs : If array is empty function uses all available directories.					
				
			!! If len of elements longer than zero, than isotopes will be ignored and elements used!!

		'''
		runs=[]

		HDF5_surf=[]
		HDF5_out=[]
		HDF5_restart=[]
		for i in range(len(self.run_dirs_name)):
			HDF5_surf.append(self.runs_H5_surf[i])
			HDF5_out.append(self.runs_H5_out[i])
			HDF5_restart.append(self.runs_H5_restart[i])
		legend_k=0
		sefiles=[]
		sefiles_hout=[]
		sefiles_restart=[]
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		
				sefiles_hout.append(se(HDF5_out[i]))
				sefiles_restart.append(se(HDF5_restart[i]))
		
		isotopes11=sefiles[0].se.isotopes
		if len(isotopes)>0:
			if isotopes[0]=='allstable':
				isotopes=get_stable(isotopes11,get_elements=False)
				color=len(isotopes)*[color[0]]
				marker_type=len(isotopes)*[marker_type[0]]
				line_style=len(isotopes)*[line_style[0]]
				markersize=len(isotopes)*[markersize[0]]
				line_width=len(isotopes)*[line_width[0]]
				print isotopes
		if len(elements)>0:
			if elements[0]=='allstable':
				elements=get_stable(isotopes11,get_elements=True)
				color=len(elements)*[color[0]]
				marker_type=len(elements)*[marker_type[0]]
				line_style=len(elements)*[line_style[0]]
				markersize=len(elements)*[markersize[0]]
				line_width=len(elements)*[line_width[0]]
				print elements
	

		####dealign with explosion pinput
                #test to take right dir:
                #
                if delay:
                        exp_type='delay'
                else:
                        exp_type='rapid'
                slist = os.listdir(exp_dir)
                expr = re.compile(exp_type)
                slist=(filter(expr.search,slist))


		exp_runs_H5_restart=[]
		sefiles_exp=[]


		#to dinstiungish between pre exp and exp sources
		if pre_exp==False:

			slist = os.listdir(exp_dir)
			expr = re.compile(exp_type)
			slist=(filter(expr.search,slist))

			for element in slist:
				run_path=exp_dir+'/'+element
				if not os.path.isdir(run_path):
					continue
				if os.path.isdir(run_path+"/H5_restart"):
					sefiles1 = os.listdir(run_path+"/H5_restart")
					if (filter(expr.search,sefiles1)) <1:
						print "Warning: No hdf5 restart files found in "+run_path+"/H5_restart"
					else:
						exp_runs_H5_restart.append(run_path+"/H5_restart")
						sefiles_exp.append(se(run_path+"/H5_restart"))
		else:
	
			exp_runs_H5_restart=HDF5_restart
			for k in range(len(HDF5_restart)):
				sefiles_exp.append(se(HDF5_restart[k]))


		z_index_files=[]
		z_values=[]
		j=-1
		for i in range(len(HDF5_surf)):
			j+=1
			star_z=sefiles[i].get("zini")
			if star_z not in z_values:
				z_values.append(star_z)
				z_index_files.append([])
			z_index_files[ z_values.index(star_z)].append(i)
		max_yield=[]
		color_iso=-1
		yields_1=[]
		t=0
		if len(elements)>0:
			isotopes=elements
		legend_k=0
		for w in range(len(z_index_files)):
			star_mass_array=[]
			yields=[]
			legend_k+=1
			#iso_yields=[]
			#iniabu_yields_folded=[]
			production_factor=[]
			for i in range(len(isotopes)):
				#iso_yields.append(np.zeros(len( z_index_files[w]   )))
				#iniabu_yields_folded.append(np.zeros(len( z_index_files[w]   )))
				if exp_only==False:	
					production_factor.append(np.zeros(len( z_index_files[w]   )))
				else:
					production_factor.append([])
			ttt=0
			for k in z_index_files[w]:
				if type(sefiles[k].get("mini")) == np.ndarray:
					star_mass=sefiles[k].get("mini")[0]
					star_z=sefiles[k].get("zini")[0]
				else:
                                        star_mass=sefiles[k].get("mini")
                                        star_z=sefiles[k].get("zini")
				print 'mini',star_mass
					
				#star_mass=sefiles[k].get("mini")[0]
				#star_z=sefiles[k].get("zini")[0]


                                if cycles[k][1]==-1:
					print 'cycles'
					#print sefiles_restart[k].se.cycles
                                        endcycle=int(sefiles_restart[k].se.cycles[-2])
                                else:
                                        endcycle=cycles[k][1]
                                cycs=range(cycles[k][0],endcycle,cycles[k][2])
                                if endcycle not in cycs:
                                        cycs=cycs+[endcycle]
                                star_mass11=sefiles[k].get(cycs,'mass')
                                if not len(star_mass11)==len(cycs):
                                        cycs1=[int(cycs[0])]
                                        for cyc in cycs[1:]:
                                                if  cyc <= max(cycs1):
                                                        continue
                                                a=sefiles[k].get(cyc,'mass')
                                                cyctest=int(cyc)
                                                #if cycle 0 or a non-existing cycle, or cycle was alrady added, 
                                                #go one cycle further   
                                                while ((type(a)==list) or (a==0)):
                                                        if   cyctest<int(cycs[-1]):
                                                                print int(cycs[-1])
                                                                print 'found bad cycle',cyctest
                                                                #if last cycle does not exist, go backwards
                                                                #if cyctest == cycs[-1]:
                                                                #       cyctest-=1
                                                                #else:
                                                                cyctest+=1
                                                                a=sefiles[k].get(cyctest,'mass')
                                                        elif cyctest==int(cycs[-1]):
                                                                cyctest=int(sefiles_restart[k].se.cycles[-3])
                                                                a=sefiles[k].get(cyctest,'mass')
                                                                h=-4
                                                                while ((type(a)==list) or (a==0)):
                                                                        cyctest=int(sefiles_restart[k].se.cycles[h])
                                                                        a=sefiles[k].get(cyctest,'mass')
                                                                        h-=1
                                                                for kk in cycs1:
                                                                        if kk>=cyctest:
                                                                                cycs1=cycs1[:kk]
                                                                break
                                                        else:
                                                                break
                                                print 'Use cycle ',cyctest
                                                if cyctest>max(cycs1):
                                                        cycs1.append(cyctest)
                                        print len(cycs),len(cycs1)
                                        cycs=cycs1
                                print cycs[-5:]
				if exp_only==False:
					star_mass_array.append(star_mass)
					prod_factor,isotopes_prod_fac,yields,iniabu,mass_frac_ini,remn_mass = self.weighted_yields(sefiles[k],sefiles_restart[k],isotopes,elements,cycs)

					#print 'yield output###################################:'
					#print 'wind: ',yields,prod_factor,iniabu
					for tt in range(len(prod_factor)):
						if yields_output==False:
                                        		prod_factor[tt]=yields[tt]/iniabu[tt]
						else:
							prod_factor[tt]=yields[tt]

				else:
					prod_factor=[]
					if star_mass<=8:
						continue
					for pp in range(len(isotopes)): 
						production_factor[pp].append([0])
							
					star_mass_array.append(star_mass)
					
				mass=star_mass
				metallicity=star_z
				if mass>8:

					#in case no ppn_exp run dir is available, skip explosion contribution
					#this is because mass cut is then larger then actual mass of star...no mass lost

					for t in range(len(exp_runs_H5_restart)):
						if type(sefiles_exp[t].get("mini")) == np.ndarray:
							mass1=sefiles_exp[t].get("mini")[0]
							metallicity1=sefiles_exp[t].get("zini")[0]
						else:
							mass1=sefiles_exp[t].get("mini")
							metallicity1=sefiles_exp[t].get("zini")

						###for currupted files

						corrupt=False
						if type(mass1) == list:
							corrupt=True
						if type(metallicity1) ==list:
							corrupt=True
						if ((float(mass1) < 1) or  (float(mass1) >70 )) or type(mass1) == type(None):
							corrupt=True
						if ((float(metallicity1) < 0.00001) or  (float(metallicity1) >1. )) or type(metallicity1) == type(None):
							corrupt=True


						if corrupt == True:
							metallicity1=0.02
							rundir11=exp_runs_H5_restart[t].split('/')[-2]
							mass1=float(rundir11.split('.')[0][1:])
							#print 'corruption, assinge new values',mass1,metallicity1
						if mass1 == mass and metallicity1 == metallicity:
							#sefiles_re2=sefiles_exp[t]
							#calculate exp part
							import nugridse as mp
							reload(mp)
							#Have to find the right cycle in restart file
							#.se.cycles does not work when reading all files	
							#therefore identify file with last cycle and read it:
							refiles=sefiles_exp[t].se.files
							tosort=[]
							refiles1=[]
							for i in range(len(refiles)):
								cycle_file=refiles[i][-18:-11]
								if cycle_file.isdigit():
									tosort.append(int(cycle_file))
									refiles1.append(refiles[i])
							idx_sorted=sorted(range(len(tosort)), key=lambda k: tosort[k])
							idx_file = idx_sorted[-1]
							if pre_exp==True:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t])
							else:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_file])
							if len(sefiles_re_cycle.se.cycles)==0:
								#in the case there is not cycle in the last file...set1.2 M2-
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_sorted[-2]])
                                                        if pre_exp==True:
                                                                cycleend=cycs[-1]
                                                        else:
                                                                cycleend=sefiles_re_cycle.se.cycles[-1]
							#print 'identifies as ',exp_runs_H5_restart[t],'with file',refiles1[idx_file]
							isotopes_exp,yields_exp,iniabu_exp,mass_cut_exp,first_cycle_exp= self.weighted_yields_explosion(mp,sefiles_re_cycle,sefiles[k],isotopes=isotopes,elements=elements,cycleend=cycleend,delay=delay)


							if exp_only==False:	
								#add pre-exp yields
								for tt in range(len(prod_factor)):
									#print 'star above 8M'
									#print isotopes_prod_fac[tt],yields[tt],yields_exp[tt]
									#the iniabu from the explosion is not added
									#because it is already in the iniabu from wind
									if yields_output==False:
										prod_factor[tt]=(yields[tt]+yields_exp[tt])/(iniabu[tt] + iniabu_exp[tt])#(iniabu[tt])
									else:
										prod_factor[tt]=yields[tt]+yields_exp[tt]
									#print 'yields wind:',yields[tt]
									#print 'prodfac wind:',(yields[tt]/iniabu[tt])
									#print 'wind+exp yields : ',yields[tt]+yields_exp[tt]	
							else:
								for tt in range(len(isotopes_exp)):	
									if yields_output==False:	
										prod_factor.append((yields_exp[tt])/(iniabu_exp[tt]))
									else:
										prod_factor.append(yields_exp[tt])
							remn_mass=mass_cut_exp
							#print 'exp yield',yields_exp
							#print 'prodfac exp:',(yields_exp[tt]/iniabu_exp[tt])
							#print 'wind+exp prodfac:',prod_factor[tt]
							break		


		

				#print 'iso exp',len(isotopes_exp)
				#prod_factor=np.log10(prod_factor)
				print len(isotopes),len(prod_factor),len(production_factor)
				for i in range(len(isotopes)):
                                        production_factor[i][ttt]=prod_factor[i]
				ttt+=1
			
			###plotting
			#if plot_set_diagram==True:
			#if plot_set_diagram==True:
			#	fig_2=plt.figure(isotopes[h])
			#	#fig_2=plt.figure()
			#	ax_2 = fig_2.add_subplot(1,1,1)
			#	plt.rcParams.update({'font.size': 16})
			#	plt.rc('xtick', labelsize=16)
			#	plt.rc('ytick', labelsize=16)
			#	if log==True:
			#		ax_2.set_yscale('log')
                        color_iso=0
                        legend_k=0

			# prod_factor [i][h] with i for isotopes
			return isotopes,star_mass_array,production_factor


	def prod_fac_massgrid(self,weighted=True,log=True,runs=[],isotopes=[],elements=[],cycles=[],plot_set_diagram=True,color=['r','b','g','k'],marker_type=['o','p','s','D'],line_style=['--','-','-.',':'],markersize=[6,6,6,6],line_width=[14,14,14,14],title='',withlabel=True,withlabelelement=False,label='',plot_lines=True,exp_only=False,pre_exp=False,delay=True,exp_dir='',yields_output=False,fallback_coords=[],fallback_masses=[]):
		'''
			Plots behaviour star mass dependent yields of different isotopes - specify dirs at the beginning
			Beware  of different z,

			runs : If array is empty function uses all available directories.					
				
			!! If len of elements longer than zero, than isotopes will be ignored and elements used!!

		'''
	
                self.fallback_coords=fallback_coords
                self.fallback_masses=fallback_masses


		runs=[]

		HDF5_surf=[]
		HDF5_out=[]
		HDF5_restart=[]
		for i in range(len(self.run_dirs_name)):
			HDF5_surf.append(self.runs_H5_surf[i])
			HDF5_out.append(self.runs_H5_out[i])
			HDF5_restart.append(self.runs_H5_restart[i])
		legend_k=0
		sefiles=[]
		sefiles_hout=[]
		sefiles_restart=[]
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		
				sefiles_hout.append(se(HDF5_out[i]))
				sefiles_restart.append(se(HDF5_restart[i]))
		
		isotopes11=sefiles[0].se.isotopes
		if len(isotopes)>0:
			if isotopes[0]=='allstable':
				isotopes=get_stable(isotopes11,get_elements=False)
				color=len(isotopes)*[color[0]]
				marker_type=len(isotopes)*[marker_type[0]]
				line_style=len(isotopes)*[line_style[0]]
				markersize=len(isotopes)*[markersize[0]]
				line_width=len(isotopes)*[line_width[0]]
				print isotopes
		if len(elements)>0:
			if elements[0]=='allstable':
				elements=get_stable(isotopes11,get_elements=True)
				color=len(elements)*[color[0]]
				marker_type=len(elements)*[marker_type[0]]
				line_style=len(elements)*[line_style[0]]
				markersize=len(elements)*[markersize[0]]
				line_width=len(elements)*[line_width[0]]
				print elements
	

		####dealign with explosion pinput
                #test to take right dir:
                #
                if delay:
                        exp_type='delay'
                else:
                        exp_type='rapid'
                slist = os.listdir(exp_dir)
                expr = re.compile(exp_type)
                slist=(filter(expr.search,slist))


		exp_runs_H5_restart=[]
		sefiles_exp=[]


		#to dinstiungish between pre exp and exp sources
		if pre_exp==False:

			slist = os.listdir(exp_dir)
			expr = re.compile(exp_type)
			slist=(filter(expr.search,slist))

			for element in slist:
				run_path=exp_dir+'/'+element
				if not os.path.isdir(run_path):
					continue
				if os.path.isdir(run_path+"/H5_restart"):
					sefiles1 = os.listdir(run_path+"/H5_restart")
					if (filter(expr.search,sefiles1)) <1:
						print "Warning: No hdf5 restart files found in "+run_path+"/H5_restart"
					else:
						exp_runs_H5_restart.append(run_path+"/H5_restart")
						sefiles_exp.append(se(run_path+"/H5_restart"))
		else:
	
			exp_runs_H5_restart=HDF5_restart
			for k in range(len(HDF5_restart)):
				sefiles_exp.append(se(HDF5_restart[k]))


		z_index_files=[]
		z_values=[]
		j=-1
		for i in range(len(HDF5_surf)):
			j+=1
			star_z=sefiles[i].get("zini")
			if star_z not in z_values:
				z_values.append(star_z)
				z_index_files.append([])
			z_index_files[ z_values.index(star_z)].append(i)
		max_yield=[]
		color_iso=-1
		yields_1=[]
		t=0
		if len(elements)>0:
			isotopes=elements
		legend_k=0
		for w in range(len(z_index_files)):
			star_mass_array=[]
			yields=[]
			legend_k+=1
			#iso_yields=[]
			#iniabu_yields_folded=[]
			production_factor=[]
			for i in range(len(isotopes)):
				#iso_yields.append(np.zeros(len( z_index_files[w]   )))
				#iniabu_yields_folded.append(np.zeros(len( z_index_files[w]   )))
				if exp_only==False:	
					production_factor.append(np.zeros(len( z_index_files[w]   )))
				else:
					production_factor.append([])
			ttt=0
			for k in z_index_files[w]:
				if type(sefiles[k].get("mini")) == np.ndarray:
					star_mass=sefiles[k].get("mini")[0]
					star_z=sefiles[k].get("zini")[0]
				else:
                                        star_mass=sefiles[k].get("mini")
                                        star_z=sefiles[k].get("zini")
				print 'mini',star_mass
					
				#star_mass=sefiles[k].get("mini")[0]
				#star_z=sefiles[k].get("zini")[0]


                                if cycles[k][1]==-1:
					print 'cycles'
					#print sefiles_restart[k].se.cycles
                                        endcycle=int(sefiles_restart[k].se.cycles[-2])
                                else:
                                        endcycle=cycles[k][1]
                                cycs=range(cycles[k][0],endcycle,cycles[k][2])
                                if endcycle not in cycs:
                                        cycs=cycs+[endcycle]
                                star_mass11=sefiles[k].get(cycs,'mass')
                                if not len(star_mass11)==len(cycs):
                                        cycs1=[int(cycs[0])]
                                        for cyc in cycs[1:]:
                                                if  cyc <= max(cycs1):
                                                        continue
                                                a=sefiles[k].get(cyc,'mass')
                                                cyctest=int(cyc)
                                                #if cycle 0 or a non-existing cycle, or cycle was alrady added, 
                                                #go one cycle further   
                                                while ((type(a)==list) or (a==0)):
                                                        if   cyctest<int(cycs[-1]):
                                                                print int(cycs[-1])
                                                                print 'found bad cycle',cyctest
                                                                #if last cycle does not exist, go backwards
                                                                #if cyctest == cycs[-1]:
                                                                #       cyctest-=1
                                                                #else:
                                                                cyctest+=1
                                                                a=sefiles[k].get(cyctest,'mass')
                                                        elif cyctest==int(cycs[-1]):
                                                                cyctest=int(sefiles_restart[k].se.cycles[-3])
                                                                a=sefiles[k].get(cyctest,'mass')
                                                                h=-4
                                                                while ((type(a)==list) or (a==0)):
                                                                        cyctest=int(sefiles_restart[k].se.cycles[h])
                                                                        a=sefiles[k].get(cyctest,'mass')
                                                                        h-=1
                                                                for kk in cycs1:
                                                                        if kk>=cyctest:
                                                                                cycs1=cycs1[:kk]
                                                                break
                                                        else:
                                                                break
                                                print 'Use cycle ',cyctest
                                                if cyctest>max(cycs1):
                                                        cycs1.append(cyctest)
                                        print len(cycs),len(cycs1)
                                        cycs=cycs1
                                print cycs[-5:]
				if exp_only==False:
					star_mass_array.append(star_mass)
					prod_factor,isotopes_prod_fac,yields,iniabu,mass_frac_ini,remn_mass = self.weighted_yields(sefiles[k],sefiles_restart[k],isotopes,elements,cycs)

					#print 'yield output###################################:'
					#print 'wind: ',yields,prod_factor,iniabu
					for tt in range(len(prod_factor)):
						if yields_output==False:
                                        		prod_factor[tt]=yields[tt]/iniabu[tt]
						else:
							prod_factor[tt]=yields[tt]

				else:
					prod_factor=[]
					if star_mass<=8:
						continue
					for pp in range(len(isotopes)): 
						production_factor[pp].append([0])
							
					star_mass_array.append(star_mass)
					
				mass=star_mass
				metallicity=star_z
				if mass>8:

					#in case no ppn_exp run dir is available, skip explosion contribution
					#this is because mass cut is then larger then actual mass of star...no mass lost

					for t in range(len(exp_runs_H5_restart)):
						if type(sefiles_exp[t].get("mini")) == np.ndarray:
							mass1=sefiles_exp[t].get("mini")[0]
							metallicity1=sefiles_exp[t].get("zini")[0]
						else:
							mass1=sefiles_exp[t].get("mini")
							metallicity1=sefiles_exp[t].get("zini")

						###for currupted files

						corrupt=False
						if type(mass1) == list:
							corrupt=True
						if type(metallicity1) ==list:
							corrupt=True
						if ((float(mass1) < 1) or  (float(mass1) >70 )) or type(mass1) == type(None):
							corrupt=True
						if ((float(metallicity1) < 0.00001) or  (float(metallicity1) >1. )) or type(metallicity1) == type(None):
							corrupt=True


						if corrupt == True:
							metallicity1=0.02
							rundir11=exp_runs_H5_restart[t].split('/')[-2]
							mass1=float(rundir11.split('.')[0][1:])
							#print 'corruption, assinge new values',mass1,metallicity1
						if mass1 == mass and metallicity1 == metallicity:
							#sefiles_re2=sefiles_exp[t]
							#calculate exp part
							import nugridse as mp
							reload(mp)
							#Have to find the right cycle in restart file
							#.se.cycles does not work when reading all files	
							#therefore identify file with last cycle and read it:
							refiles=sefiles_exp[t].se.files
							tosort=[]
							refiles1=[]
							for i in range(len(refiles)):
								cycle_file=refiles[i][-18:-11]
								if cycle_file.isdigit():
									tosort.append(int(cycle_file))
									refiles1.append(refiles[i])
							idx_sorted=sorted(range(len(tosort)), key=lambda k: tosort[k])
							idx_file = idx_sorted[-1]
							if pre_exp==True:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t])
							else:
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_file])
							if len(sefiles_re_cycle.se.cycles)==0:
								#in the case there is not cycle in the last file...set1.2 M2-
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_sorted[-2]])
                                                        if pre_exp==True:
                                                                cycleend=cycs[-1]
                                                        else:
                                                                cycleend=sefiles_re_cycle.se.cycles[-1]
							print 'files exp from ',sefiles_re_cycle.se.filename
							#print 'identifies as ',exp_runs_H5_restart[t],'with file',refiles1[idx_file]
							isotopes_exp,yields_exp,iniabu_exp,mass_cut_exp,first_cycle_exp= self.weighted_yields_explosion(mp,sefiles_re_cycle,sefiles[k],isotopes=isotopes,elements=elements,cycleend=cycleend,delay=delay)


							if exp_only==False:	
								#add pre-exp yields
								for tt in range(len(prod_factor)):
									#print 'star above 8M'
									#print isotopes_prod_fac[tt],yields[tt],yields_exp[tt]
									#the iniabu from the explosion is not added
									#because it is already in the iniabu from wind
									if yields_output==False:
                                                                                print '##################'
                                                                                print 'pre exp and exp: ',yields[tt],yields_exp[tt]
										prod_factor[tt]=(yields[tt]+yields_exp[tt])/(iniabu[tt] + iniabu_exp[tt])#(iniabu[tt])
									else:
										#print '##################'
										#print 'pre exp and exp: ',yields[tt],yields_exp[tt]
										prod_factor[tt]=yields[tt]+yields_exp[tt]
									#print 'yields wind:',yields[tt]
									#print 'prodfac wind:',(yields[tt]/iniabu[tt])
									#print 'wind+exp yields : ',yields[tt]+yields_exp[tt]	
							else:
								for tt in range(len(isotopes_exp)):	
									if yields_output==False:	
										prod_factor.append((yields_exp[tt])/(iniabu_exp[tt]))
									else:
										prod_factor.append(yields_exp[tt])
							remn_mass=mass_cut_exp
							#print 'exp yield',yields_exp
							#print 'prodfac exp:',(yields_exp[tt]/iniabu_exp[tt])
							#print 'wind+exp prodfac:',prod_factor[tt]
							break		


		

				#print 'iso exp',len(isotopes_exp)
				#prod_factor=np.log10(prod_factor)
				print len(isotopes),len(prod_factor),len(production_factor)
				for i in range(len(isotopes)):
                                        production_factor[i][ttt]=prod_factor[i]
				ttt+=1
			
			###plotting
			#if plot_set_diagram==True:
			#if plot_set_diagram==True:
			#	fig_2=plt.figure(isotopes[h])
			#	#fig_2=plt.figure()
			#	ax_2 = fig_2.add_subplot(1,1,1)
			#	plt.rcParams.update({'font.size': 16})
			#	plt.rc('xtick', labelsize=16)
			#	plt.rc('ytick', labelsize=16)
			#	if log==True:
			#		ax_2.set_yscale('log')
                        color_iso=0
                        legend_k=0


			for h in range(len(isotopes)):
				#yield_1=iso_yields[i]
				#mass_1=star_mass_array
				yield_1=[]
				mass_1=[]			
				iniabu_folded=[]
				prod_fac_sorted=[]	
               			indices_array=sorted(range(len(star_mass_array)),key=lambda x:star_mass_array[x])
                		for i in indices_array:
				#	yield_1.append(iso_yields[h][i])
					mass_1.append(star_mass_array[i])
				#	iniabu_folded.append(iniabu_yields_folded[h][i])
					prod_fac_sorted.append(production_factor[h][i])
				####Plotting prodfactor
				#plt.figure(fig_2.number)

				fig_2=plt.figure(isotopes[h])
				#ax_2 = fig_2.add_subplot(1,1,1)
				ax_2=fig_2.gca()
				#plt.rcParams.update({'font.size': 16})
				#plt.rc('xtick', labelsize=16)
				#plt.rc('ytick', labelsize=16)
				if log==True:
					ax_2.set_yscale('log')
				ax_2.legend()
				ax_2.set_xlabel("M/M$_{\odot}$",fontsize=16)
				ax_2.minorticks_on()
				if yields_output==False:
					ax_2.set_ylabel("Overproduction factor",fontsize=16)
				else:
					ax_2.set_ylabel("Yields/M$_{\odot}$",fontsize=16)
				ax_2.set_title(title)
	
				if len(label)>0:
					label=", "+label
				if len(elements)==0:
					plot_quantity=isotopes[h]
					plot_quantity="$^{"+plot_quantity.split("-")[1]+"}$"+plot_quantity.split("-")[0]
				else:
					plot_quantity=elements[h]
				if withlabel==True:
					if withlabelelement==True:
                                		ax_2.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso],label=plot_quantity+" , Z="+str(star_z)+label  )
					else:
                                        	ax_2.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso],label="Z="+str(star_z)+label  )
					#plt.legend(loc=2)
					print 'legend created'
				else:
                                	ax_2.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso])				
				

				m_max=max(star_mass_array)+2
				
				if plot_lines==True:
					ax_2.plot([0,m_max],[1,1],"k--",linewidth=3)
					ax_2.plot([0,m_max],[2,2],"k--",linewidth=1)
					ax_2.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)				
	
				color_iso+=1
				legend_k+=1
		#####
		x_imf=[]
		y_imf=[]
		m_max=max(star_mass_array)+2.
		#if plot_set_diagram==True:
		#	plt.figure(fig_2.number)
		##	###plot figure 1
		#	plt.legend()
		##	plt.xlabel("M/M$_{\odot}$",fontsize=20)
		#	plt.minorticks_on()
		#	plt.ylabel("Overproduction factor",fontsize=20)
		#	plt.title(title)


	
	def prod_fac_massgrid_single_wind(self,weighted=True,log=True,runs=[],isotopes=[],elements=[],cycles=[],plot_set_diagram=True,color=['r','b','g','k'],marker_type=['o','p','s','D'],line_style=['--','-','-.',':'],markersize=[6,6,6,6],line_width=12*[6],title='',withlabel=False,label='wind',plot_lines=True):
		'''
			Plots behaviour star mass dependent yields of different isotopes - specify dirs at the beginning
			Beware  of different z
			The isotopes should be being produced because of the logscale...
			At the moment 4 isotopes supported

			runs : If array is empty function uses all available directories.					
				
			!! If len of elements longer than zero, than isotopes will be ignored and elements used!!

			" -Final wind yields - isotopic composition folded with IMF" 
			dirs=['M1.650Z0.0001', 'M2.000Z0.0001','M3.000Z0.0001','M5.000Z0.0001']
			dirs=['M1.650Z0.0010', 'M2.000Z0.0010','M3.000Z0.0010','M5.000Z0.0010']
			dirs=['M1.650Z0.0060', 'M2.000Z0.0060','M3.000Z0.0060','M5.000Z0.0060']		
		
			cycles=[[0,80000,1000],[0,80000,1000],[0,80000,1000],[0,120000,1000]]	
		
			pap.weighted_yields_massgrid(runs=['M1.650Z0.0001', 'M2.000Z0.0001','M3.000Z0.0001','M5.000Z0.0001'],isotopes=["He-4","C-12","O-16"],cycles=[[0,10000,1000],[0,10000,1000],[0,10000,1000],[0,10000,1000]])
	
		'''
		
		runs=[]
		HDF5_surf=[]
		HDF5_out=[]
		HDF5_restart=[]
		for i in range(len(self.run_dirs_name)):
				HDF5_surf.append(self.runs_H5_surf[i])
				HDF5_out.append(self.runs_H5_out[i])
				HDF5_restart.append(self.runs_H5_restart[i])
		legend_k=0
		sefiles=[]
		sefiles_hout=[]
		sefiles_restart=[]
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		
				sefiles_hout.append(se(HDF5_out[i]))
				sefiles_restart.append(se(HDF5_restart[i]))


		isotopes11=sefiles[0].se.isotopes

		if len(isotopes)>0:
			if isotopes[0]=='allstable':
				isotopes=get_stable(isotopes11,get_elements=False)
				color=len(isotopes)*[color[0]]
				marker_type=len(isotopes)*[marker_type[0]]
				line_style=len(isotopes)*[line_style[0]]
				markersize=len(isotopes)*[markersize[0]]
				line_width=len(isotopes)*[line_width[0]]
				print isotopes
		if len(elements)>0:
			if elements[0]=='allstable':
				elements=get_stable(isotopes11,get_elements=True)
				color=len(elements)*[color[0]]
				marker_type=len(elements)*[marker_type[0]]
				line_style=len(elements)*[line_style[0]]
				markersize=len(elements)*[markersize[0]]
				line_width=len(elements)*[line_width[0]]
				print elements
	


		z_index_files=[]
		z_values=[]
		j=-1
		for i in range(len(HDF5_surf)):
			j+=1
			star_z=sefiles[i].get("zini")
			if star_z not in z_values:
				z_values.append(star_z)
				z_index_files.append([])
			z_index_files[ z_values.index(star_z)].append(i)
		max_yield=[]
		color_iso=-1
		yields_1=[]
		t=0

		if len(elements)>0:
			isotopes=elements

		legend_k=0
		for w in range(len(z_index_files)):
			star_mass_array=[]
			yields=[]
			legend_k+=1
			iso_yields=[]
			iniabu_yields_folded=[]
			production_factor=[]
			for i in range(len(isotopes)):
				iso_yields.append(np.zeros(len( z_index_files[w]   )))
				iniabu_yields_folded.append(np.zeros(len( z_index_files[w]   )))
				production_factor.append(np.zeros(len( z_index_files[w]   )))
			ttt=0
			for k in z_index_files[w]:
				if type(sefiles[k].get("mini")) == np.ndarray:
					star_mass=sefiles[k].get("mini")[0]
					star_z=sefiles[k].get("zini")[0]
				else:
                                        star_mass=sefiles[k].get("mini")
                                        star_z=sefiles[k].get("zini")
				print 'mini',star_mass				
				#star_mass=sefiles[k].get("mini")[0]
				#star_z=sefiles[k].get("zini")[0]				
				star_mass_array.append(star_mass)
                                if cycles[k][1]==-1:
                                        endcycle=int(sefiles_restart[k].se.cycles[-2])
                                else:
                                        endcycle=cycles[k][1]
                                cycs=range(cycles[k][0],endcycle,cycles[k][2])
                                if endcycle not in cycs:
                                        cycs=cycs+[endcycle]
                                star_mass11=sefiles[k].get(cycs,'mass')
                                #do a check for missing cycles
                                if not len(star_mass11)==len(cycs):
                                        cycs1=[int(cycs[0])]
                                        for cyc in cycs[1:]:
                                                if  cyc <= max(cycs1):
                                                        continue
                                                a=sefiles[k].get(cyc,'mass')
                                                cyctest=int(cyc)
                                                #if cycle 0 or a non-existing cycle, or cycle was alrady added, 
                                                #go one cycle further   
                                                while ((type(a)==list) or (a==0)):
                                                        if   cyctest<int(cycs[-1]):
                                                                print int(cycs[-1])
                                                                print 'found bad cycle',cyctest
                                                                #if last cycle does not exist, go backwards
                                                                #if cyctest == cycs[-1]:
                                                                #       cyctest-=1
                                                                #else:
                                                                cyctest+=1
                                                                a=sefiles[k].get(cyctest,'mass')
                                                        elif cyctest==int(cycs[-1]):
                                                                cyctest=int(sefiles_restart[k].se.cycles[-3])
                                                                a=sefiles[k].get(cyctest,'mass')
                                                                h=-4
                                                                while ((type(a)==list) or (a==0)):
                                                                        cyctest=int(sefiles_restart[k].se.cycles[h])
                                                                        a=sefiles[k].get(cyctest,'mass')
                                                                        h-=1
                                                                for kk in cycs1:
                                                                        if kk>=cyctest:
                                                                                cycs1=cycs1[:kk]
                                                                break
                                                        else:
                                                                break
                                                print 'Use cycle ',cyctest
                                                if cyctest>max(cycs1):
                                                        cycs1.append(cyctest)
                                        print len(cycs),len(cycs1)
					cycs=cycs1
				print cycs[-5:]

				prod_factor,isotopes_prod_fac,yields,iniabu,mass_frac_ini,remn_mass =self.weighted_yields(sefiles[k],sefiles_restart[k],isotopes,elements,cycs)

				#print 'yield output###################################:'
				#print  'yields: ',yields
				#print 'prod fact: ',prod_factor
				#print 'iniabu: ',iniabu 
				#print 'isotope: ',isotopes_prod_fac
				mass=star_mass
				metallicity=star_z

				#prod_factor=np.log10(prod_factor)
				for i in range(len(isotopes)):
                                        production_factor[i][ttt]=yields[i]/iniabu[i]
				ttt+=1
			
			###plotting
			color_iso=0

			if plot_set_diagram==True:
				
				fig_2=plt.figure(isotope)
				ax_2 = fig_2.add_subplot(1,1,1)
				plt.rcParams.update({'font.size': 16})
				plt.rc('xtick', labelsize=16)
				plt.rc('ytick', labelsize=16)
				if log==True:
					ax_2.set_yscale('log')
                        color_iso=0
                        legend_k=0


			for h in range(len(isotopes)):
				#yield_1=iso_yields[i]
				#mass_1=star_mass_array
				yield_1=[]
				mass_1=[]			
				iniabu_folded=[]
				prod_fac_sorted=[]	
               			indices_array=sorted(range(len(star_mass_array)),key=lambda x:star_mass_array[x])
                		for i in indices_array:
				#	yield_1.append(iso_yields[h][i])
					mass_1.append(star_mass_array[i])
				#	iniabu_folded.append(iniabu_yields_folded[h][i])
					prod_fac_sorted.append(production_factor[h][i])
				####Plotting prodfactor
				#plt.figure(fig_2.number)
				if plot_set_diagram==False:
					fig_2=plt.figure(isotopes[h])
					ax_2 = fig_2.gca() #fig_2.add_subplot(1,1,1)
					plt.rcParams.update({'font.size': 16})
					plt.rc('xtick', labelsize=16)
					plt.rc('ytick', labelsize=16)
					if log==True:
						ax_2.set_yscale('log')
					ax_2.legend()
					ax_2.set_xlabel("M/M$_{\odot}$",fontsize=16)
					ax_2.minorticks_on()
					ax_2.set_ylabel("Overproduction factor",fontsize=16)
					ax_2.set_title(title)
				if len(label)>0:
					label1=", "+label
				if len(elements)==0:
					plot_quantity=isotopes[h]
					plot_quantity="$^{"+plot_quantity.split("-")[1]+"}$"+plot_quantity.split("-")[0]
				else:
					plot_quantity=elements[h]
				if withlabel==True:
                                	plt.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso],label=plot_quantity+label1  )				
				else:
					plt.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso])                     

				m_max=max(star_mass_array)+2
				
				if plot_lines==True:
					plt.plot([0,m_max],[1,1],"k--",linewidth=3)
					plt.plot([0,m_max],[2,2],"k--",linewidth=1)
					plt.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)				
	
				color_iso+=1
				legend_k+=1
		#####
		x_imf=[]
		y_imf=[]
		m_max=max(star_mass_array)+2.
		if plot_set_diagram==True:
			fig_2=plt.figure(isotopes[h])
			plt.figure(fig_2.number)
			###plot figure 1
			plt.legend()
			plt.xlabel("M/M$_{\odot}$",fontsize=16)
			plt.minorticks_on()
			plt.ylabel("Overproduction factor",fontsize=16)
			plt.title(title)

	


	def set_pocket(self,runs=[],cycle=[45532,47566],mass_cell=[0.641981,0.641981],isotopes=[],x_charge=False,mass_number_range=[],decayed=False,iso_label=False,title="PDCZ",colors=["red","blue","green","black"],yax_log=True,filename_norm="iniab2.0E-02GN93.ppn",elem_norm=''):

		'''
			plots isotopic composition for different runs, each with specific cycle number cycle[i] and specific mass cell mass_cell[i].
			Initial abundace file for normalization can be chosen as filename_norm, but MUST be in run directory

			decayed - plot only decayed isotopes
			mass_number_range - array of min mass number and max mass number, isotopes inbetween will be plottet, if set, isotope array will be ignored
			iso_label - set labels of isotopes True or False
			title - 
			colors -
			filename_norm - not yet included, normalization with initial abundance file
			....
			e.g.
			
			setse.set_pocket(cycle=[32840,29390],mass_cell=[0.6300,0.6076],isotopes=["He-4","C-12","O-16","Ne-22"],mass_number_range=[80,140],decayed=True)
			

	
		'''
	
		import nugridse as mp	
                sefiles=[]
                legend=[]
                HDF5_out=[]
                extra_label=[]

		if len(runs) ==0:
			HDF5_out=self.runs_H5_out
			runs=self.run_dirs_name
			extra_label=self.extra_label
		else:
			for i in range(len(self.run_dirs_name)):
				if self.run_dirs_name[i] in runs:
					HDF5_out.append(self.runs_H5_out[i])
                                        extra_label.append(self.extra_label[i])
		for i in range(len(HDF5_out)):
			reload(mp)
			file_norm=HDF5_out[i][:-6]+filename_norm
			sefiles=se(HDF5_out[i])
                        mass=sefiles.get("mini")
                        z=sefiles.get("zini")
                        legend=str(mass)+"M$_{\odot}$ Z= "+str(z)+", "+extra_label[i]
			if len(isotopes)==0:
                                self.plot_abu_atmasscoord(mp,sefiles,cycle[i],mass_cell[i],["He-4","C-12","O-16","Ne-22"],x_charge,mass_number_range,decayed,file_norm,elem_norm,yax_log,label=iso_label,legend=legend,color=colors[i],title=title)

				#plt.figure(111);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=["He-4","C-12","O-16","Ne-22"],label=True,legend=legend,color=colors[i],title=title)
				#plt.figure(222);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=['Fe-56','Co-60','Ni-61','Cu-65','Zn-67','Ga-71','Ge-73','As-75','Se-77','Br-82','Kr-84','Rb-87','Sr-88','Y-89','Zr-92','Nb-94','Zr-96','Mo-96','Ba-137','Ba-138','La-139','Ce-140','Nd-142','Sm-144','Pb-206'],label=True,legend=legend,color=colors[i],title=title)
			else:
				#plt.figure(111);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=isotopes,label=True,legend=legend,color=colors[i],title=title)
				self.plot_abu_atmasscoord(mp,sefiles,cycle[i],mass_cell[i],isotopes,x_charge,mass_number_range,decayed,file_norm,elem_norm,yax_log,label=iso_label,legend=legend,color=colors[i],title=title)


        def set_surface_abundances(self,runs=[],cycle=[-1,-1],isotopes=[],x_charge=True,mass_number_range=[],decayed=False,iso_label=False,title="",colors=["red","blue","green","black"],yax_log=True,filename_norm="iniab2.0E-02GN93.ppn",elem_norm=''):

                '''
                        plots isotopic composition for different runs, each with specific cycle number cycle[i] and specific mass cell mass_cell[i].
                        Initial abundace file for normalization can be chosen as filename_norm, but MUST be in run directory

                        decayed - plot only decayed isotopes
                        mass_number_range - array of min mass number and max mass number, isotopes inbetween will be plottet, if set, isotope array will be ignored
                        iso_label - set labels of isotopes True or False
                        title - 
                        colors -
                        filename_norm - not yet included, normalization with initial abundance file
                        ....
                        e.g.
                        
                        setse.set_pocket(cycle=[32840,29390],mass_cell=[0.6300,0.6076],isotopes=["He-4","C-12","O-16","Ne-22"],mass_number_range=[80,140],decayed=True)
                        

        
                '''


                import nugridse as mp
                sefiles=[]
                legend=[]
                HDF5_out=[]
                extra_label=[]


                if len(runs) ==0:
                        HDF5_out=self.runs_H5_out
                        runs=self.run_dirs_name
                        extra_label=self.extra_label
                else:
                        for i in range(len(self.run_dirs_name)):
                                if self.run_dirs_name[i] in runs:
                                        HDF5_out.append(self.runs_H5_out[i])
                                        extra_label.append(self.extra_label[i])
                for i in range(len(HDF5_out)):
                        reload(mp)
                        file_norm=HDF5_out[i][:-6]+filename_norm
                        sefiles=se(HDF5_out[i])
			if cycle[i] ==-1:
				cycle_1=sefiles.se.cycles[-1]
			else:
				cycle_1=cycle[i]
			masscell=sefiles.get(cycle_1,"mass")[0]
                        mass=sefiles.get("mini")
                        z=sefiles.get("zini")
                        legend=str(mass)+"M$_{\odot}$ Z= "+str(z)+", "+extra_label[i]
                        if len(isotopes)==0:
                                self.plot_abu_atmasscoord(mp,sefiles,cycle_1,masscell,isotopes,x_charge,mass_number_range,decayed,file_norm,elem_norm,yax_log,label=iso_label,legend=legend,color=colors[i],title=title)

                                #plt.figure(111);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=["He-4","C-12","O-16","Ne-22"],label=True,legend=legend,color=colors[i],title=title)
                                #plt.figure(222);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=['Fe-56','Co-60','Ni-61','Cu-65','Zn-67','Ga-71','Ge-73','As-75','Se-77','Br-82','Kr-84','Rb-87','Sr-88','Y-89','Zr-92','Nb-94','Zr-96','Mo-96','Ba-137','Ba-138','La-139','Ce-140','Nd-142','Sm-144','Pb-206'],label=True,legend=legend,color=colors[i],title=title)
                        else:
                                #plt.figure(111);self.pocket_composition(mp,sefiles=sefiles,cycle=cycle[i], massbot=massrange[i][0],masstop=massrange[i][1],isotopes=isotopes,label=True,legend=legend,color=colors[i],title=title)
                                self.plot_abu_atmasscoord(mp,sefiles,cycle_1,masscell,isotopes,x_charge,mass_number_range,decayed,file_norm,elem_norm,yax_log,label=iso_label,legend=legend,color=colors[i],title=title)


	def plot_saga_data(self,fig=7,label_plot=["C-rich RGB","CEMP RGB","EMP RGB"],path="/nfs/rpod3/critter/DPG/saga_platform"):
	
		'''
			Plots data from the SAGA database. It uses the scripts in 
			forum.astro.keele.ac.uk/utils/spectroscopy_data/saga_platform
			Adapt them for your purpose and use this function to combine 
			SAGA data and your own simulation data in one plot

	
		'''

		figure_data=plt.figure(fig)
		pwd=os.getcwd()
		os.chdir(path)
		execfile(path+"/"+"read_files_images.py")
		print pwd
		os.chdir(pwd)
		print "done"	
		


	def set_surface_plots(self,runs=[],cycles=[[10000,1000],[10000,1000]],t0_model=[],decayed=True,mesarunpath="",mesa_multi_dir=[],ini_abu_file="/astro/critter/critter/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn",x_range=[-1.8,0.1],withmarker=True,y_range=[]):
		'''
		Plot surface abundance evolution of multiple runs

		plots [hs/ls] vs star age
		plots  [ls/Fe] vs [hs/ls] 
		[Rb/Sr] [hs/ls]
		plots Mg25/Mg24 vs Mg26/Mg24
          	plots plt.xlabel("$^{96}$Zr/$^{94}$Zr") plt.ylabel("$^{152}$Gd/$^{154}$Gd")
		plots [ls/Fe] vs star age
		cycles - 2d array, containing the start and stop cycle for each run and the sparsity
		t0_model - the starting model to plot, for the time dependent plot, t0_model will be..
		x_range - define range of 
	
		!!Check/c setse.run_dirs_name array for order of runs for setting the cycles parameter!!	
	
		# creating prefix variable, we have agreed that this is an 11 character
		# name with the initial mass and metallicity encoded: M5.00Z0.010 
		
		e.g. 
		
		set_surface_plots(self,runs=[],cycles=[[3113,50000,1000],[3031,60000,1000]],t0_model=[3113,3031],x_range=[],y_range=[]	
		test2_template: [3620, 5145, 3030, 4757, 11997, 5288, 9470, 5123, 2996, 3312, 5828, 4749]

		['./M3.000Z0.0060/LOGS',
		 './M5.000Z0.0010/LOGS',
		 './M2.000Z0.0001/LOGS',
		 './M3.000Z0.0001/LOGS',
		 './M1.650Z0.0060/LOGS',
		 './M3.000Z0.0010/LOGS',
		 './M1.650Z0.0010/LOGS',
		 './M5.000Z0.0060/LOGS',
		 './M2.000Z0.0010/LOGS',
		 './M2.000Z0.0060/LOGS',
	 	'./M1.650Z0.0001/LOGS',
		 './M5.000Z0.0001/LOGS
	
		'''	
		

                color=['r','b','k','g']	#'r','r','r','b','b','k','g','k','k','g','b','g']
                marker_type=['o','o','o','p','p','p','D','D','D','s','s','s']
		line_style=3*['--','-','-.',':']
		#IDEA PASS PARAMETE BETWEEN BOTH CLASSES
	        #if len(runs)==0:
                #	t0_model=self.set_find_first_TP()
		#cycles=[]
		#for i in range(len(cyc)):
		#	cycles.append([t0_model[i],cyc[i][0],cyc[i][1]])
		
		print "Do you want a TP-dependece on time axis?"
		print "Read corresponding MESA runs?"
		#answer=raw_input("Read corresponding MPPNP runs?")
		#if answer ==True:
		if len(mesarunpath) >0:
		#mesa_path=raw_input("Read corresponding MPPNP runs?")
			dir_1=mesarunpath
			multi_dir=[]
			if len(mesa_multi_dir)>0:
                                for i in range(len(self.run_dirs_name)):
                                        multi_dir.append(mesarunpath+"/"+mesa_multi_dir[i])
			else:
				for i in range(len(self.run_dirs_name)):
					multi_dir.append(mesarunpath+"/"+self.run_dirs_name[i])
			mesaset=mesa_set(multi_dir=multi_dir)
			peak_lum_model_array_1,h1_mass_min_DUP_model_array = mesaset.multi_DUP(t0_model=[],plot_fig=False)
			t0_model=[]
			peak_lum_model_array=[]
			cycles_1=[]
			for i in range(len(peak_lum_model_array_1)):
				cycles_1.append(map(int,np.array(peak_lum_model_array_1[i])))
				t0_model.append(int(peak_lum_model_array_1[i][0]))
			if len(cycles)==0:
				cycles=cycles_1
		print "Cycles: ",cycles

		sefiles=[]
		legend=[]
		HDF5_surf=[]
		extra_label=[]
		if len(runs) ==0:
			HDF5_surf=self.runs_H5_surf
			runs=self.run_dirs_name
			extra_label=self.extra_label
		else:
			for i in range(len(self.run_dirs_name)):
				if self.run_dirs_name[i] in runs:
					HDF5_surf.append(self.runs_H5_surf[i])
					extra_label.append(self.extra_label[i])
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		

		for i in range(len(HDF5_surf)):							
			mass=sefiles[i].get("mini")[0]	
			z=sefiles[i].get("zini")
			if z <0.001:
				color_1=color[0]
			elif z <0.006:
				color_1=color[1]
			else:
				color_1=color[2]	
			legend=str(mass)+"M$_{\odot}$ Z= "+str(z)+", "+extra_label[i]
			self.surface_plots(HDF5_surf[i],cycles[i],t0_model[i],decayed,legend,ini_abu_file,withmarker,marker_type[i],color_1,line_style[i],title="")
	
		#HDF5_dirs=["/nfs/rpod3/critter/Results/test2_template_mppnp/M1.650Z0.0001","/nfs/rpod3/critter/Results/test2_template_mppnp/M5.000Z0.0001"],	
			
	
		plt.figure(2)	
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()		
		plt.plot([-1000,1000],[0,0],"-",c="k",linewidth=2)
		plt.plot([0,0],[-1000,1000],"-",c="k",linewidth=2)		
		if xmax<0:
			xmax=0.1
		if xmin >0:
			xmin=-0.1	
		if ymax<0:
			ymax=0.1	
		if ymin >0:
			ymin=-0.1		
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)					
		plt.figure(3)	
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()	
		plt.plot([-1000,1000],[0,0],"-",c="k",linewidth=2)
		plt.plot([0,0],[-1000,1000],"-",c="k",linewidth=2)	
		if xmax<0:
			xmax=0.1
		if xmin >0:
			xmin=-0.1	
		if ymax<0:
			ymax=0.1	
		if ymin >0:
			ymin=-0.1
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)			
			
		plt.figure(4)	
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()	
		plt.plot([-1000,1000],[0,0],"-",c="k",linewidth=2)
		plt.plot([0,0],[-1000,1000],"-",c="k",linewidth=2)	
		if xmax<0:
			xmax=0.1
		if xmin >0:
			xmin=-0.1	
		if ymax<0:
			ymax=0.1	
		if ymin >0:
			ymin=-0.1
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)		
	
		plt.figure(5)	
		xmin, xmax = plt.xlim()
		ymin, ymax = plt.ylim()	
		plt.plot([-1000,1000],[0,0],"-",c="k",linewidth=2)
		plt.plot([0,0],[-1000,1000],"-",c="k",linewidth=2)	
		if xmax<0:
			xmax=0.1
		if xmin >0:
			xmin=-0.1	
		if ymax<0:
			ymax=0.1	
		if ymin >0:
			ymin=-0.1
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)		
	
		if len(x_range)>0:
			plt.figure(1);plt.xlim(x_range[0],x_range[1])
			plt.figure(2);plt.xlim(x_range[0],x_range[1])
			plt.figure(3);plt.xlim(x_range[0],x_range[1])
			plt.figure(4);plt.xlim(x_range[0],x_range[1])
			plt.figure(5);plt.xlim(x_range[0],x_range[1])


		plt.figure(6)



	def set_triple_isotope_plot(self,runs=[],xiso=['Zr',96,'Zr',94],yiso=['Zr',92,'Zr',94],graintype=['sic'],C_star_only=True,dcycle=100,USEEPP_path='USEEPP',extra_label=[],sens_runs=True):
		'''
        		xiso:       give isotopes as ['Fe',57,'Fe',56]. x axis, yaxis
	


		'''

		import os
		pwd=os.getcwd()
		iniabufiles=[]			
		
                sefiles=[]
                legend=[]
                HDF5_surf=[]
                if len(runs) ==0:
                        HDF5_surf=self.runs_H5_surf
                        runs=self.run_dirs_name
			if len(extra_label)==0:
                        	extra_label=self.extra_label
                else:
                        for i in range(len(self.run_dirs_name)):
                                if self.run_dirs_name[i] in runs:
                                        HDF5_surf.append(self.runs_H5_surf[i])
                                        extra_label.append(self.extra_label[i])
		if sens_runs==True:
			import glob
			import os
			color=['r','r','b','b','k','g','k','g','k','b','r','g','g','b']
			marker_type=['o','D','s','p','o','D','s','p','p','o','p','o','s','p']
			line_style=2*['--','-','-.',':']
			first_plot=True
			for i in range(len(HDF5_surf)):
				sefiles=se(HDF5_surf[i])
                        	mass=sefiles.get("mini")
                        	z=sefiles.get("zini")
	                        new_notation='{:.1E}'.format(float(z))
        	                iniabufile="iniab"+new_notation+"GN93.ppn"     #pwd+"/"+USEEPP_path+'/'+"iniab"+new_notation+"GN93.ppn"
                        	legend=str(mass)+"M$_{\odot}$, "+str(z)+"Z$_{\odot}$, "+extra_label[i]
				if first_plot==True:
					sefiles.plot_isoratios(xiso=xiso,yiso=yiso,graintype=graintype,tp_finding=False,C_star_only=C_star_only,deltax=True,deltay=True,logx=False,logy=False,title=None,legend=True,extra_legend=legend,dcycle=dcycle,errbar=True,iniabufile=iniabufile,plt_symb=marker_type[i],plt_col=color[i])		
					first_plot=False			
				else:
					print ""				
					sefiles.plot_isoratios(xiso=xiso,yiso=yiso,tp_finding=False,C_star_only=C_star_only,deltax=True,deltay=True,logx=False,logy=False,title=None,legend=True,extra_legend=legend,dcycle=dcycle,errbar=True,iniabufile=iniabufile,plt_symb=marker_type[i],plt_col=color[i])		


		#plot_isoratios(self,xiso,yiso,graintype=None,tp_finding=False,deltax=True,deltay=True,logx=False,logy=False,title=None,legend=True,dcycle=500,errbar=True,iniabufile='iniab2.0E-02GN93.ppn'):




	def chem(self,imf=-2.35,Z=0.001,dt=5e9,table='isotope_yield_table.txt',sn1a_table='sn1a_t86.txt'):
		'''
			Method to perform simple chemical evolution
			Uses:

		'''

		#### Using the IMF

		yields,masses=self.read_yield_tables(table,Z)

		salp_fac_wind,salp_fac_exp,sum_mass_term=self.imf(imf,masses)
		#mass_ini_wind=[1.0,1.65,2.0,3.0,5.0,6.0,7.0,12.0,15.0,20.0,25.0]
		
		weighted_yields=[]
		isotopes=[]
		lifetimes=[]
		#getting the remnants
		remn_masses=[]
		remn_number=[]
		for  star_mass in masses:
			print "Use set 1 extension weightening for",star_mass
			indx=masses.index(star_mass)
			weight_factor=salp_fac_wind[indx]
			weighted_yields.append(np.array(yields.get(star_mass,Z,'Yields'))*weight_factor*sum_mass_term)		
			remn_masses.append(yields.get(star_mass,Z,'Mfinal'))
			remn_number.append(weight_factor*sum_mass_term)
			lifetimes.append(yields.get(star_mass,Z,'age'))
			#iniabu_yields=np.array(iniabu)*weight_factor*sum_mass_term

		#### Got ini masses, weighed yield and lifetimes

		isotopes=yields.get(star_mass,Z,'Isotopes')
		
		#Read SN1a contribution

		sn1a_yields=self.read_yield_sn1a_tables(sn1a_table,isotopes)
		weighted_sn1a_yields=number1a*np.array(sn1a_yields)


		#storage
		yields_t=[[0]*len(weighted_yields[0])]  #initial yields to zero at t=0
                t_contri=[0] #time variable
		dt_steps=[0] #timestep sizes

		remn_t=[[0]] #evolution of remnant distribution

                ##Contribution times/Ejection
                contr_idx=sorted(range(len(masses)), key=lambda k: lifetimes[k])
                contribution_times=[]
                for idx in contr_idx:
                        contribution_times.append(lifetimes[idx])
                print 'Contribution at ',contribution_times
		
		#Do the evolution steps in sizes of contribution times

		i=-1
		t=0
		t_last=0
		for time in contribution_times:
			#time
			i+=1
			t=time
			#full stop if time > dt
			if t > dt:
				break
			print 'time: '+str(time)
			t_contri.append(t)
			dt_steps.append(t - t_last)
			t_last=time
			#increase yields due to contribution
			yields_new=np.array(yields_t[-1]) + np.array(weighted_yields[contr_idx[i]])
			yields_t.append(yields_new)

			#dealing with remnants
			#Assume only the current generation which goes of produces SN1a
			#
			#mass[contr_idx[i]]
			if masses[contr_idx[i]] < 8:
				print 'SN1a contribution'
				number1a=remn_number[contr_idx[i]]
				#assume half of the AGB stars become SN1a
				yields_new=np.array(yields_t[-1]) + 0.5*weighted_sn1a_yields
			remn_new=np.array(remn_t[-1]) + remn_number[contr_idx[i]]
			remn_t.append(remn_new)
		
		#save the history
		self.gce_yields_t=yields_t
		self.gce_t_contri=t_contri
		self.gce_isotopes=isotopes		
	
	def chem_plot(self,xaxis,yaxis):

		'''
			Plotting tool for the chem chemical evolution function
		'''	
	
		#Assume isotope input
		if '-' in xaxis:
			yields=self.gce_yields_t
			iso_idx=self.gce_isotopes.index(xaxis)
			x=[]
			for k in range(len(yields)):
				x.append(yields[k][iso_idx])
		if '-' in yaxis:
                        yields=self.gce_yields_t
                        iso_idx=self.gce_isotopes.index(yaxis)
                        y=[]
                        for k in range(len(yields)):
                                y.append(yields[k][iso_idx])
			
		if 'time' in xaxis:
			times=self.gce_t_contri
			x=times
		plt.plot(x,y)
		plt.xlabel(xaxis)
		plt.ylabel(yaxis)


	def imf(self,imf,masses):
		
		'''
			
			Calculates the IMF by using an adaption of the mass grid of paper1 (based on Marcos scripts).
			INcludes masses till 25 since set1.1 is the only one with 32 u. 60
			and plotting sets together give a more complete picture without them
		'''


		mass_ini_exp=[12.,15.,20.,25.]
		#mass_ini_wind=[1.0,1.65,2.0,3.0,5.0,6.0,7.0,12.,15.,20.,25.]
		mass_ini_wind=masses
		bounds = [0.65,1.45, 1.85, 2.5, 3.5,4.5,5.5,6.5, 8., 13.5, 17.5, 22.5, 28.5]
		salp_fac = imf #-2.35
		# the crazyness below is just because the mass between 7. and 11.5 Msun 
		# is missing, super AGB stars. Marco
		# set 1.2
		facs=[bounds[0]**salp_fac-bounds[1]**salp_fac,bounds[1]**salp_fac-bounds[2]**salp_fac,bounds[2]**salp_fac-bounds[3]**salp_fac,bounds[3]**salp_fac-bounds[4]**salp_fac,bounds[4]**salp_fac-bounds[5]**salp_fac,bounds[5]**salp_fac-bounds[6]**salp_fac,bounds[6]**salp_fac-bounds[7]**salp_fac,bounds[7]**salp_fac-bounds[8]**salp_fac,bounds[8]**salp_fac-bounds[9]**salp_fac,bounds[9]**salp_fac-bounds[10]**salp_fac,bounds[10]**salp_fac-bounds[11]**salp_fac,bounds[11]**salp_fac-bounds[12]**salp_fac]
		# set 1.1
		#facs=[bounds[0]**salp_fac-bounds[1]**salp_fac,bounds[1]**salp_fac-bounds[2]**salp_fac,bounds[2]**salp_fac-bounds[3]**salp_fac,bounds[3]**salp_fac-bounds[4]**salp_fac,bounds[5]**salp_fac-bounds[6]**salp_fac,bounds[6]**salp_fac-bounds[7]**salp_fac,bounds[7]**salp_fac-bounds[8]**salp_fac]
		weighted_facs = [x/sum(facs) for x in facs]

		sum_mass_term = 0.
		for i in range(len(mass_ini_wind)):
			sum_mass_term = weighted_facs[i]*mass_ini_wind[i] + sum_mass_term
		sum_mass_term = 1./(sum_mass_term)

		# weighted_facs = weighted_facs/sum_mass_term

		# -2 is because we not consider m=32 and m=60.
		salp_fac_wind = weighted_facs[0:len(weighted_facs)]
		# add the -2 term if 32 Msun and 60 Msun are included for salp_fac_exp. Marco
		# add for set1.2
		#salp_fac_exp = weighted_facs[len(weighted_facs)-3-2:len(weighted_facs)-2]
		# add for set1.1
		salp_fac_exp = weighted_facs[len(weighted_facs)-4:len(weighted_facs)]
		
		return salp_fac_wind,salp_fac_exp,sum_mass_term
	







	def read_yield_sn1a_tables(self,sn1a_table,isotopes):
		f1=open(sn1a_table)
		lines=f1.readlines()
		f1.close()
		iso_1a=[]
		yield_1a=[]
		for line in lines:
			#for header
			if '#' in line:
				continue
			iso_1a.append(line.split()[0])
			yield_1a.append(float(line.split()[1]))				


		yields=[]
		#fill up the missing isotope yields with zero
		for iso in isotopes:
			if iso in iso_1a:
				idx=iso_1a.index(iso)
				yields.append(yield_1a[idx])
			else:
				yields.append(0.)		
		return yields

		
	def read_yield_tables(self,table,Z):
	
		'''
			Reads out yield tables in table format
		'''
	
		print 'reading yield tables'
		import read_yield_tables as y
		y1=y.yields(table)
		masses=[]
		for header in y1.table_header:
			if str(Z) in header:
				mass=header.split(',')[0].split('=')[1]
				masses.append(float(mass))
		#print y1.get(1.65,0.01,'age')	
		#print 'Available masses: ',masses	

		return y1,masses 
	

	def write_gce_input_wind_ejection_rate(self,file_name="isotopic_table.txt",final_models=[]):

		'''
			Calculates total mass ejection rate from stellar wind phase.
		'''

		e_kin_array=[]
		mz=[]
		m_final=[]
		mz=[]
		ini_m=[]
		p=-1
		for i in range(len(self.runs_H5_out)):
			p+=1
			sefiles=se(self.runs_H5_out[i])
			sefiles_restart=se(self.runs_H5_restart[i])
			mass=sefiles.get("mini")
			metallicity=sefiles.get("zini")	
			mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
			endcycle=int(sefiles.se.cycles[-1])
			if len(final_models)>0:
				cycs=self.get_last_cycle([0,final_models[p]+500,500],sefiles,sefiles_restart)
			else:
				cycs=self.get_last_cycle([0,-1,500],sefiles,sefiles_restart)
			endcycle=cycs[-1]
			h_mass_frac=sefiles.get(endcycle,'H-1')
			mass_array=sefiles.get(endcycle,'mass')
			#earlier: 1.e-1 but now as in set1 scripts 0.05
			h_free_core=mass_array[np.where(h_mass_frac<5.e-2)[0][-1]]
			m_final.append(h_free_core)		
			mz.append(mz1)
			ini_m.append(mass)

                f1=open(file_name,'r')
                lines=f1.readlines()
                f1.close()
                i=-1
                line1=''
                while (True):
                        i+=1
                        if i>len(lines)-1:
                                break
                        line=lines[i]
                        line1+=lines[i]
                        for k in range(len(mz)):
                                if mz[k] in lines[i]:
					## check ahead for H Mfinal:
					#for h in range(i,i+10):
					#	if 'H Lifetime:' in lines[h]:
					#		lifetime=float(lines[h].split(':')[1])
					#		break
					#	if h==i+9:
					#		print 'ERROR: no lifetime in table found!!'
					#		return
					ejection=(ini_m[k]-m_final[k]) #/lifetime
                                        line1+=('H Wind ejection: '+'{:.3E}'.format(ejection)+'\n')
                                        break

                f1=open(file_name,'w')
                f1.write(line1)
                f1.close()

	def write_gce_input_exp_energy(self,file_name="isotopic_table.txt",velocities_file='velocities.txt',sparsity=500,final_models=[]):

		'''
			Calculates and writes total kinetic energy [erg]
			fromo the Wind (no exp included!) 
		'''

		e_kin_array=[]
		mz=[]
		p=-1
		for i in range(len(self.runs_H5_surf)):
			p +=1
			sefiles=se(self.runs_H5_surf[i])
			mass=sefiles.get("mini")
			metallicity=sefiles.get("zini")
			mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
			if True:	
				msol = 1.9892e33  # solar mass (g)
				standard_cgrav = 6.67428e-8 #(g^-1 cm^3 s^-2)
				rsol = 6.9598e10 # solar radius (cm)
		

				if len(final_models)>0:
					endccycle=final_models[p]
				else:	
					endcycle=int(sefiles.se.cycles[-1])
					endcycle = int(endcycle) // 100 * 100
					endcycle = int(endcycle-20)
				startcycle=int(sefiles.se.cycles[0])
				#choose in sparsity steps			
				cycs=range(startcycle,endcycle+sparsity,sparsity)	
				print 'from ',startcycle,' to ',endcycle
				star_mass=sefiles.get(cycs,'mass')
				star_age=sefiles.get(cycs,'age')
				print 'got star_mass'
				delta_star_mass=np.array(star_mass[1:])-np.array(star_mass[:-1])
				print 'got delta mass'
				star_radius=sefiles.get(cycs,'radius')
				delta_star_radius=np.array(star_radius[1:])-np.array(star_radius[:-1])
				#between two models, where I assume mass lost of delta_star_mass
				star_radius_mid=np.array(star_radius[:-1]) + (delta_star_radius/2.)
				star_mass_mid=np.array(star_mass[:-1]) + (delta_star_mass/2.)


				#Escape if Ekin=Epot
				#Ekin = 1/2*m*v**2 = GMm/r
				e_kin=0.
				v_esc=[]
				out= '(M='+str(mass)+',Z='+str(metallicity)+')\n'
				for p in range(len( star_mass_mid)):
					if star_radius_mid[p]>0:
						e_kin = e_kin + standard_cgrav*(star_mass_mid[p]*msol)*(abs(delta_star_mass[p])*msol)/(star_radius_mid[p]*rsol) #ergs
						v_esc.append( 1e-5*np.sqrt(e_kin*2./(star_mass_mid[p]*msol))) #cm/s to km/s 1e-5
						out= out + (str(p)+' & '+'{:.10E}'.format(star_age[p]/3.1558149984e7)+' & '+str(v_esc[-1])+'\n') #cm/s
				f1=open(velocities_file,'a')
				f1.write(out)
				f1.close()
				#print 'E _kin wind:',e_kin	
				#for massive stars, add explosion energy
			'''
			else:
				###Get EKin###
				#if metallicity==0.02 or metallicity==0.01:
					#print 'add exp energy'
				             #15,20,25 in foe /10**51 ergs

				#following assumptions here
				#1st: exp M12 is equal M15
				#2nd: exp M32 and M60 is equal M25
                		exp_energy={12:1.5*10**51 ,15:1.5*10**51,20:50*10**51,25:5*10**51,32:5*10**51,60:5*10**51 } #ergs
				#else:
				#	return 'Error, inpout missing'	
				print 'E _kin exp',exp_energy[int(mass)]
				e_kin=exp_energy[int(mass)]
			'''
			e_kin_array.append(e_kin)
			mz.append(mz1)	
		f1=open(file_name,'r')
		lines=f1.readlines()
		f1.close()
		i=-1
		line1=''
		while (True):
			i+=1
			if i>len(lines)-1:
				break
			line=lines[i]	
			line1+=lines[i]
			for k in range(len(mz)):
				if mz[k] in lines[i]:
					line1+=('H Wind kinetic energy: '+'{:.3E}'.format(e_kin_array[k])+'\n')
					break			
			
		f1=open(file_name,'w')
                f1.write(line1)
                f1.close()

	def write_gce_input_lifetimes(self,file_name="isotopic_table.txt",final_cycles=[]):

		'''
			Calculates and writes lifetime in file_name
			X_carbon_center >  0.4
			
		'''

		###Get lifetimes###

		print '##############GETTING THE LIFETIMES##############'
		time=[]
		model=[]
		mz=[]
		for i in range(len(self.runs_H5_surf)):
			sefiles=se(self.runs_H5_surf[i])
			mass=sefiles.get("mini")
			metallicity=sefiles.get("zini")
			mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
			cycs=[]
			if final_cycles[i]>-1:
				final_cycle=final_cycles[i]
			else:
				final_cycle=int(sefiles.se.cycles[-1])
			time1=(sefiles.get(final_cycle,'age')*sefiles.get('age_unit'))/31557600.
			time.append(time1)
			mz.append(mz1)

		f1=open(file_name,'r')
		lines=f1.readlines()
		f1.close()
		i=-1
		line1=''
		print 'Adding for ',mz,'lifetime : ',time
		while (True):
			i+=1
			if i>len(lines)-1:
				break
			line=lines[i]	
			line1+=lines[i]
			for k in range(len(mz)):
				if mz[k] in lines[i]:
					line1+=('H Lifetime: '+'{:.3E}'.format(time[k])+'\n')
					break			
			
		f1=open(file_name,'w')
                f1.write(line1)
                f1.close()



	def write_gce_input_lifetimes_old(self,file_name="isotopic_table.txt"):

		'''
			Calculates and writes lifetime in file_name
			X_carbon_center >  0.4
			
		'''

		###Get lifetimes###

		print '##############GETTING THE LIFETIMES##############'
		time=[]
		model=[]
		mz=[]
		for i in range(len(self.runs_H5_out)):
			sefiles=se(self.runs_H5_out[i])
			mass=sefiles.get("mini")
			metallicity=sefiles.get("zini")
			mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
			cycs=[]
			for k in range(5,len(sefiles.se.cycles),5):
				cycs.append(int(sefiles.se.cycles[k]))
			w=0
			for cyc in cycs:
				c12_center=sefiles.get(cyc,'C-12')[0]
				#c12_center=c12[w][0]
				w+=1
				#print c12_center, cyc
				if c12_center>1e-1:
					model.append(cyc)
					#print sefiles.get('age_unit')
					time1=(sefiles.get(cyc,'age')*sefiles.get('age_unit'))/31557600.
					time.append(time1)
					mz.append(mz1)
					break

		f1=open(file_name,'r')
		lines=f1.readlines()
		f1.close()
		i=-1
		line1=''
		print 'Adding for ',mz,'lifetime : ',time
		while (True):
			i+=1
			if i>len(lines)-1:
				break
			line=lines[i]	
			line1+=lines[i]
			for k in range(len(mz)):
				if mz[k] in lines[i]:
					line1+=('H Lifetime: '+'{:.3E}'.format(time[k])+'\n')
					break			
			
		f1=open(file_name,'w')
                f1.write(line1)
                f1.close()



	def write_gce_input(self,elem=True,cycles=20*[[0,-1,500]],file_name="isotopic_table.txt_elem",GCE_tables=True,exp_dir='',exp_type=None,yields_output=True,pre_exp=False,wind_only=False,exp_only=False,fallback_coords=[],fallback_masses=[],exp_cycle_inp=-1):

		'''
			Write yield input
		'''

                #self.fallback_coords=fallback_coords
                #self.fallback_masses=fallback_masses
                #print 'fallback coords',self.fallback_coords


		import getpass
		user=getpass.getuser()
		import time
		date=time.strftime("%d %b %Y", time.localtime())
		sefiles=se(self.runs_H5_surf[0])
		z=sefiles.get("zini")

                sefiles=se(self.runs_H5_surf[0])
                if elem==True:
			elem_write=''
			elements_2=get_stable(sefiles.se.isotopes,get_elements=True)
			for elem1 in elements_2:
				elem_write+=elem1
				if not elem1 == elements_2[-1]:
					elem_write+=', '
			print elements_2
			#elements_2=['H','C']
			isotopes_2=[]
		else:
			isotopes_1=sefiles.se.isotopes
                	isotopes_2=[]
                	available_iso=[]
			iso_write=''
                	for j in range(len(isotopes_1)):
                	        if is_stable(isotopes_1[j])=='t':
                	                isotopes_2.append(isotopes_1[j])
					iso_write+=isotopes_1[j]
					if not j ==(len(isotopes_1)-2):
						iso_write+=', '
					
                	elements_2=[]
			print 'All isotopes for GCE available'
 
                if os.path.exists(file_name):
                        print 'file already exists, update header and add new tables'
			f1=open(file_name)
			lines=f1.readlines()
			f1.close()
			old_z=[]
			new_z=[]
			new_line=''
			for line in lines:
				if 'Table' in line:
					old_z.append(line.split(',')[1].split(')')[0].split('=')[1])	
			print old_z
			if str(z) in old_z:
				new_z=old_z
			else:
				new_z=old_z+[str(z)]
			print new_z
			for line in lines:
				if 'Data prepared date' in line:
					l='H Data prepared date: '+date+'\n'
				elif 'Data prepared by' in line:
					l='H Data prepared by: '+user+'\n'
				elif 'Number of metallicities' in line:
					if not len(old_z) == len(new_z):
						l='H Number of metallicities: '+str(int(line.split(':')[1])+1)+'\n'
					else:
						l=line
				elif 'H NuGrid yields Set1' in line:
					l=line
					if not len(old_z) == len(new_z):
						sets={'0.02':'set1.2','0.01':'set1.1','0.006':'set1.3a','0.001':'set1.4a','0.0001':'set1.5a'}
						l=line[:-1]+', '+sets[new_z[-1]]+'\n'
				else:
					l=line
				new_line+=l

			f1=open(file_name,'w')			
			f1.write(new_line)
			f1.close()

			#import fileinput
			#for line in fileinput.input(file_name, inplace = 1): 
      			#	print line.replace("bar", "barrrr").resplace,


	
		else:
			f1=open(file_name,'w')
			sets={'0.02':'set1.2','0.01':'set1.1','0.006':'set1.3a','0.001':'set1.4a','0.0001':'set1.5a'}
			zz=sets[str(z)]
			headerline='H NuGrid yields Set1: '+zz+'\n' 
			headerline+='H Data prepared by: '+user+'\n'
			headerline+='H Data prepared date: '+date+'\n'
			if len(isotopes_2)>0:
				headerline+='H Isotopes: '+iso_write+'\n'
			else:
				headerline+='H Elements: '+elem_write+'\n'
			headerline+='H Number of metallicities: 1\n'
			headerline+='H Units: Msun, year, erg\n'
			f1.write(headerline)
			f1.close()
			print 'file with header created'

               ##############read new simulation input###########

                #cycles=20*[[0,-1,500]]
		#extract and write data
                #self.write_prod_fact_stellarwinds(cycles=cycles,isotopes=isotopes_2,elements=elements_2,ascii_1=False,latex=False,GCE_input_fxt=False,GCE_input=GCE_input,GCE_input_starkenburg=GCE_input_starkenburg,file_name=file_name,exp_dir=exp_dir,exp_type=exp_type,yields_output=yields_output,pre_exp=pre_exp,exp_only=exp_only,wind_only=wind_only,exp_cycle_inp=exp_cycle_inp)
		self.get_stellar_ejecta(cycles=cycles,isotopes=isotopes_2,elements=elements_2,GCE_tables=True,file_name=file_name,exp_dir=exp_dir,exp_type=exp_type,fallback_masses=fallback_masses,fallback_coords=fallback_coords,yields_output=yields_output,pre_exp=pre_exp,exp_only=exp_only,wind_only=wind_only,exp_cycle_inp=exp_cycle_inp)

		elem_number=1
		elem=['H']
		for i in range(1,len(isotopes_2)):
			if not isotopes_2[i].split('-')[0] == isotopes_2[i-1].split('-')[0]:
				elem_number+=1
				elem.append(isotopes_2[i].split('-')[0])
		print elem
		print 'Number of elements:',elem_number


        def write_gce_input_fxt3(self,input_file='maltov1.new',input_file2='maltov1.orig',output_file='maltov1.new2'):

		'''
			Function allows to add massive stars of fxt to already filled maltov file with own data
		'''

		###read in input_file and input_file2

		isotopes,bb_massfrac,type1a_yields,nova_yields,metallicities,masses,yields,label=self.read_gce_fxt_input(input_file=input_file)	
		
		isotopes2, bb_massfrac2,type1a_yields2,nova_yields2,metallicities2,masses2,yields2,label2=self.read_gce_fxt_input(input_file=input_file2)

		if not isotopes == isotopes2:
			print 'must be same number of isotopes'
			
		isotopes_out=isotopes	
		metallicities_out=metallicities
                isotopes_out=isotopes
                bb_massfrac_out=bb_massfrac
                type1a_yields_out=type1a_yields
		nova_yields_out=nova_yields	

		#file1
		idx_z=[2,3,4]
		

		#file2
		idx_z=[2,3,4]
		idx_mass_to_add=range(12,25)		
		#z_to_add=[float(1.0000E-04),float(6.0000E-03),float(1.0000E-02),float(2.0000E-02)]	
		z_to_add=[float(1.0000E-04),float(6.0000E-03),float(2.0000E-02)]
		masses_out=len(isotopes)*[0]
		yields_out=len(isotopes)*[0]
		label_out=len(isotopes)*[0]
		print metallicities_out[0]
		for t in range(len(isotopes)):
			masses_out[t]=len(metallicities_out)*[0]
			yields_out[t]=len(metallicities_out)*[0]
			label_out[t]=len(metallicities_out)*[0]
			w=0	
			for wt in range(len(metallicities_out[0])):
				#merge file1 data and file2 data
				if metallicities_out[0][wt] in z_to_add:
					yield_sum=[]
					mass_sum=[]
					label_sum=[]
					#first the non-massive stars,
					for mass in range(len(masses[t][wt])):
						yield_sum.append(yields[t][wt][mass])
						mass_sum.append(masses[t][wt][mass])
						label_sum.append(label[t][wt][mass])
					#than add the massive ones
					for mass in range(len(idx_mass_to_add)):
						yield_sum.append(yields2[t][idx_z[w]][idx_mass_to_add[mass]])
						mass_sum.append(masses2[t][idx_z[w]][idx_mass_to_add[mass]])
                                                label_sum.append(label2[t][idx_z[w]][idx_mass_to_add[mass]])
							

					masses_out[t][wt]=mass_sum
					yields_out[t][wt]=yield_sum
					label_out[t][wt]=label_sum
					w+=1
					#print wt,w
				#file1 data only
				else:
                                        masses_out[t][wt]=masses[t][wt]     
                                        yields_out[t][wt]=yields[t][wt]
                                        label_out[t][wt]=label[t][wt]
		

		###Write output

		self.write_gce_fxt_input(isotopes_out,bb_massfrac_out,type1a_yields_out,nova_yields_out,metallicities_out,masses_out,yields_out,label_out,output_file)


	def write_GCE_input_fxt4(self,input_file='maltov1.orig',output_file='maltov1.new',first_change=True,label_masses='set1.1',cycles=20*[[0,-1,500]],exp_dir='',delay=False,bb_inp='/astro/critter/critter/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab1.0E-04GN93_alpha.ppn',isotopes_inp=[]):

		'''
			Function allows to create maltov files only from simulation data
		'''

		#to get sn1a and nova data
                isotopesfxtscheme,bb_massfrac,type1a_yields,nova_yields,metallicities,masses,yields,label=self.read_gce_fxt_input(input_file=input_file)


		##############read new simulation input###########

                isotopes=["H-1","H-2","He-3","He-4","Li-6","Li-7","Be-9","B-10","B-11","C-12","C-13","N-14","N-15","O-16","O-17","O-18","F-19","Ne-20","Ne-21","Ne-22","Na-23","Mg-24","Mg-25","Mg-26","Al-26","Al-27","Si-28","Si-29","Si-30","P-31","S-32","S-33","S-34","S-36","Cl-35","Cl-37","Ar-36","Ar-38","Ar-40","K-39","K-40","K-41","Ca-40","Ca-42","Ca-43","Ca-44","Ca-46","Ca-48","Sc-45","Ti-46","Ti-47","Ti-48","Ti-49","Ti-50","V-50","V-51","Cr-50","Cr-52","Cr-53","Cr-54","Mn-55","Fe-54","Fe-56","Fe-57","Fe-58","Fe-60","Co-59","Ni-58","Ni-60","Ni-61","Ni-62","Ni-64","Cu-63","Cu-65","Zn-64","Zn-66","Zn-67","Zn-68","Zn-70","rem"]
		if len(isotopes_inp)>0:
			#custom isotope number
			isotopesfxtscheme1=[]
			#bb_massfrac
			type1a_yields1=[]
			nova_yields1=[]
			for i in range(len(isotopes_inp)):
				#if data from fxt is available
				if isotopes_inp[i] in isotopes:
					idx=isotopes.index(isotopes_inp[i])
					type1a_yields1.append(type1a_yields[idx])
					nova_yields1.append(nova_yields[idx])
					isotopesfxtscheme1.append(isotopesfxtscheme[idx])
				else:
				#else set to zero till availabl
					type1a_yields1.append( len(type1a_yields[0])*[0.] )
					nova_yields1.append( len(nova_yields[0])*[0.] )
					isotopesfxtscheme1.append(isotopes_inp[i].replace('-','').lower())
			type1a_yields=type1a_yields1
			nova_yields=nova_yields1
			isotopesfxtscheme=isotopesfxtscheme1
			isotopes=isotopes_inp



		print isotopes		

                #check if GCE isotopes are all available in simulation input
                sefiles=se(self.runs_H5_surf[0])
                isotopes_1=sefiles.se.isotopes
                isotopes_2=[]
                no_iso=[]
                available_iso=[]
                for j in range(len(isotopes_1)):
                        if is_stable(isotopes_1[j])=='t':
                                isotopes_2.append(isotopes_1[j])
                for j in range(len(isotopes)):
                        #if isotope is not available stop       
                        if isotopes[j] not in isotopes_2:
                                print 'Isotope ',isotopes[j],'is not available'
                                no_iso.append(j)
                                #print 'STOP'
                        else:
                                available_iso.append(isotopes[j])
                #cycles=20*[[0,-1,500]]
                labelout,cols,remn_masses_1 = self.write_prod_fact_stellarwinds(cycles=cycles,isotopes=available_iso,ascii_1=False,latex=False,GCE_input_fxt=True,exp_dir=exp_dir,delay=delay)

                #to include remnants
                remn_masses=[]
                #remn_masses_1=20*[3]
                for i in range(len(remn_masses_1)):
                        #remnant_1='{:.4E}'.format(float(remn_masses_1[i]))
                        remn_masses.append(remn_masses_1[i])

                #add values for missing isotopes, value=99
		#print 'testteest'
		#print no_iso
                for i in no_iso:
                        for j in range(len(cols)):
				if i == no_iso[-1]:
				# to add the remnant masses
					print 'adding normal value for rem'
					cols[j].insert(i,remn_masses[j])
				else:
                                	cols[j].insert(i,0)

                #because the set runs have the same Z:
                #but first item in label is 'specie', so skip this
		#trick here, measure z only from first input file
		#massive stars are therefore not checked and treated as other z
                metallicity_new=float(labelout[1].split('=')[1])
                #metallicity1='{:.4E}'.format(float(metallicity))
                masses_new1=[]
                for i in range(1,len(labelout)):
                        mass=labelout[i].split('Msun')[0]
                        #mass='{:.4E}'.format(float(mass))
                        masses_new1.append(float(mass))
                #print 'Available input masses',masses

	
		yields_new=[0]*len(cols[0])
		masses_new=[0]*len(cols[0])
		label_new=[0]*len(cols[0])
		
		#for every mass
		for k in range(len(yields_new)):	
			yields_new[k]=[[]]
			masses_new[k]=[masses_new1]
			label_new[k]=[[label_masses]*len(masses_new1)]
			for w in range(len(cols)):
				yields_new[k][0].append(cols[w][k])

		#print metallicity_new
                #print yields_new
		#print label_new
		#print masses   

	
                isotopes_out=isotopesfxtscheme
		print isotopes_out
		print isotopes
		if len(bb_inp)<1:
                	bb_massfrac_out=bb_massfrac
		else:
			###get ini abu of stars and feed in as bb abu				
			bb_massfrac_out=[]
  	       		#import utils as u
        	        initial_abu=iniabu(bb_inp)
			names=initial_abu.names
			names1=[]
			for j in range(len(names)):
				names1.append(names[j].replace(' ',''))
			for i in range(len(isotopes_out)):
				if isotopes_out[i] in names1:
					iniabu_1=initial_abu.abu[ initial_abu.names.index(isotopes[i].split("-")[0].lower()+(5-len(isotopes[i])+1)*" "+isotopes[i].split("-")[1])]
					bb_massfrac_out.append(iniabu_1)	
				else:
					bb_massfrac_out.append(0.)

                type1a_yields_out=type1a_yields
		nova_yields_out=nova_yields	

		print 'bb massfrac'
		print bb_massfrac_out
		print masses_new
		print yields_new
		print metallicity_new

		masses_out=masses_new
		yields_out=yields_new
		label_out=label_new
		metallicities_out=len(isotopes_out)*[[metallicity_new]]
		





		###Write output

		self.write_gce_fxt_input(isotopes_out,bb_massfrac_out,type1a_yields_out,nova_yields_out,metallicities_out,masses_out,yields_out,label_out,output_file)





	def read_gce_fxt_input(self,input_file='maltov1.new4'):
	
                '''
			For reading gce fxt input
                '''

	
	
                import re
                f1=open(input_file)
                lines=f1.readlines()
                f1.close()
                isotopes=[]
		isotopes_out=[]
                iso_line_idx=[]
                i=0

		##########read input file##############
                #gce isotopes
                for line in lines:
                        if line[0] != ' ' and line[0:3] != 'rem':
                                if line[0].isalpha():
                                        a= line[:6].split()[0]
					isotopes_out.append(a)
                                        match = re.match(r"([a-z]+)([0-9]+)",a, re.I)
                                        items = match.groups()
                                        isotopes.append(items[0].capitalize()+'-'+items[1])
                                        iso_line_idx.append(i)
                        if line[0:3] == 'rem':
				isotopes_out.append('rem')
                                isotopes.append('rem')
                                iso_line_idx.append(i)
                        i+=1

		bb_massfrac=[0]*len(isotopes)
		type1a_yields=[0]*len(isotopes)
		nova_yields=[0]*len(isotopes)
		metallicities=[0]*len(isotopes)
		masses=[0]*len(isotopes)
		yields=[0]*len(isotopes)
		label=[0]*len(isotopes)

		i=-1
		j=-1
		while (True):
			i+=1
                        if i>len(lines)-1:
                                break
			line=lines[i]
			if i in iso_line_idx:
				j+=1
			#read data
                        if 'big bang' in line:
				
                                bb_massfrac[j]=float(line[:13])
                                continue
                        #type 1a input ignore
                        if 'number of type1a' in line:
                                number=int(line[0])
				type1a_yields[j]=[0]*number
				k=0
				for sn in range(i+1,i+number+1):
					type1a_yields[j][k]=(float(lines[sn][:13]))				
					k+=1
				i=i+number	
				continue
                        #nova input ignore
                        if 'number of nova' in line:
                                number=int(line[0])
				k=0
				nova_yields[j]=[0]*number
                                for nova in range(i+1,i+number+1):
                                         nova_yields[j][k]=(float(lines[nova][:13])) 
                                	 k+=1
				i=i+number   
                                continue
			if 'number of metallicity values' in line:
				number_metallicities=int(line[0])
				metallicities[j]=[0]*(number_metallicities)
				masses[j]=[0]*(number_metallicities)
				yields[j]=[0]*(number_metallicities)
				label[j]=[0]*(number_metallicities)	
				z_counter=-1
				continue
			if '* metallicity' in line:
				z_counter+=1
				metallicities[j][z_counter]=float(line[:13])				
				continue
			if '* number of masses @ this metallicity' in line:
				number_masses=int(line[:2])
				masses[j][z_counter]=[0]*(number_masses)
				yields[j][z_counter]=[0]*(number_masses)
				label[j][z_counter]=[0]*(number_masses)
				k=0
				for mass in range(i+1,i+number_masses+1):
					masses[j][z_counter][k]=float(lines[mass][:13])
					yields[j][z_counter][k]=float(lines[mass][14:24])
					#check for label
					if len(lines[mass])>27:
						label1=lines[mass][27:].strip()
									
					label[j][z_counter][k] = label1
					k+=1
				i=i+number_masses
                                continue
                #print bb_massfrac
                #print type1a_yields
                #print nova_yields
                #print metallicities
                #print masses
                #print yields
		##print label	
		
		return isotopes_out,bb_massfrac,type1a_yields,nova_yields,metallicities,masses,yields,label


	def write_gce_fxt_input(self,isotopes,bb_massfrac,type1a_yields,nova_yields,metallicities,masses,yields,label,output_file):

		'''
			For writing gce fxt output
		'''

		isotopes_out=isotopes
		bb_massfrac_out=bb_massfrac
                type1a_yields_out=type1a_yields
                nova_yields_out=nova_yields
                metallicities_out=metallicities
                masses_out=masses
                yields_out=yields
                label_out=label 		
		newline=''
		for i in range(len(isotopes)):
			newline+= (isotopes_out[i].ljust(29)+'* symbol name'+'\n') 	
			newline+= ('  '+'{:.4E}'.format(bb_massfrac_out[i]).ljust(27)+'* big bang mass fraction'+'\n')
			newline+= (str(len(type1a_yields_out[i])).ljust(29)+'* number of type1a'+'\n')
			ia_label=['w7 tny86','sub chandra 1','sub chandra 2','sub chandra 7','sub chandra 8','nse 1a']
			for k in range(len(type1a_yields_out[i])):
				newline+= ('  '+'{:.4E}'.format(type1a_yields_out[i][k]).ljust(29)+ia_label[k]+'\n')
			newline+= (str(len(nova_yields_out[i])).ljust(29)+'* number of nova'+'\n')
			nova_label='* nova yield 1.0 1.25 & 1.35 msun'
			for k in range(len(nova_yields_out[i])):
				if k==0:
					newline+= ('  '+'{:.4E}'.format(nova_yields_out[i][k]).ljust(27)+nova_label+'\n')
				else:
					newline+= ('  '+'{:.4E}'.format(nova_yields_out[i][k]).ljust(27)+'\n')
			newline+= (str(len(metallicities_out[i])).ljust(29)+'* number of metallicity values'+'\n')
			for k in range(len(metallicities_out[i])):
				newline+= ('  '+'{:.4E}'.format(metallicities_out[i][k]).ljust(27)+'* metallicity'+'\n') 
				newline+= (str(len(masses_out[i][k])).ljust(29)+'* number of masses @ this metallicity'+'\n')
				for w in range(len(masses_out[i][k])):
					#the following is needed for the right label setting
					if w==0:
						newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+(7*' ')+label_out[i][k][w]+'\n')
						label_old=label_out[i][k][w]	
					elif label_out[i][k][w] != label_old:
						newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+(7*' ')+label_out[i][k][w]+'\n')	
						label_old=label_out[i][k][w]
					else:
						newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+'\n')
			


		###write out data
		
		#to fix bug where fxt GCe cannot deal with 99 yields
		newline=newline.replace('9.9000E+01','0.0000E+00')
		f2=open(output_file,'w')
		f2.write(newline)
		f2.close()

		


	def write_gce_input_fxt2(self,input_file='maltov1.new4',input_file2='isotopic_table_set1.1_yields_exp.txt',output_file='maltov1.new5',exp_type='delay',label_masses='set1.2'):
	
                '''
                        reads input_file for gce code of fxt and input for in format of set1 tables and add yields
			sum up yields. usefull when combining wind and exp data. output file is input for gce code of fxt.
			define explosion type, delay or rapdid

                '''

	
	
                import re
                f1=open(input_file)
                lines=f1.readlines()
                f1.close()
                isotopes=[]
		isotopes_out=[]
                iso_line_idx=[]
                i=0
		first_change=True

		##########read input file##############
                #gce isotopes
                for line in lines:
                        if line[0] != ' ' and line[0:3] != 'rem':
                                if line[0].isalpha():
                                        a= line[:6].split()[0]
					isotopes_out.append(a)
                                        match = re.match(r"([a-z]+)([0-9]+)",a, re.I)
                                        items = match.groups()
                                        isotopes.append(items[0].capitalize()+'-'+items[1])
                                        iso_line_idx.append(i)
                        if line[0:3] == 'rem':
				isotopes_out.append('rem')
                                isotopes.append('rem')
                                iso_line_idx.append(i)
                        i+=1

		bb_massfrac=[0]*len(isotopes)
		type1a_yields=[0]*len(isotopes)
		nova_yields=[0]*len(isotopes)
		metallicities=[0]*len(isotopes)
		masses=[0]*len(isotopes)
		yields=[0]*len(isotopes)
		label=[0]*len(isotopes)

		i=-1
		j=-1
		while (True):
			i+=1
                        if i>len(lines)-1:
                                break
			line=lines[i]
			if i in iso_line_idx:
				j+=1
			#read data
                        if 'big bang' in line:
				
                                bb_massfrac[j]=float(line[:13])
                                continue
                        #type 1a input ignore
                        if 'number of type1a' in line:
                                number=int(line[0])
				type1a_yields[j]=[0]*number
				k=0
				for sn in range(i+1,i+number+1):
					type1a_yields[j][k]=(float(lines[sn][:13]))				
					k+=1
				i=i+number	
				continue
                        #nova input ignore
                        if 'number of nova' in line:
                                number=int(line[0])
				k=0
				nova_yields[j]=[0]*number
                                for nova in range(i+1,i+number+1):
                                         nova_yields[j][k]=(float(lines[nova][:13])) 
                                	 k+=1
				i=i+number   
                                continue
			if 'number of metallicity values' in line:
				number_metallicities=int(line[0])
				metallicities[j]=[0]*(number_metallicities)
				masses[j]=[0]*(number_metallicities)
				yields[j]=[0]*(number_metallicities)
				label[j]=[0]*(number_metallicities)	
				z_counter=-1
				continue
			if '* metallicity' in line:
				z_counter+=1
				metallicities[j][z_counter]=float(line[:13])				
				continue
			if '* number of masses @ this metallicity' in line:
				number_masses=int(line[:2])
				masses[j][z_counter]=[0]*(number_masses)
				yields[j][z_counter]=[0]*(number_masses)
				label[j][z_counter]=[0]*(number_masses)
				k=0
				for mass in range(i+1,i+number_masses+1):
					masses[j][z_counter][k]=float(lines[mass][:13])
					yields[j][z_counter][k]=float(lines[mass][14:24])
					#check for label
					if len(lines[mass])>27:
						label1=lines[mass][27:].strip()
									
					label[j][z_counter][k] = label1
					k+=1
				i=i+number_masses
                                continue
                #print bb_massfrac
                #print type1a_yields
                #print nova_yields
                #print metallicities
                #print masses
                #print yields
		##print label	
		
		##############read new simulation input###########

		f1=open(input_file2)
		lines_txt=f1.readlines()
		f1.close()
		isotopes_1=[]
		for i in range(len(lines_txt)):
			if lines_txt[i][0]=='H':
				continue
			if 'specie' in lines_txt[i]:
				continue
			name1=lines_txt[i].split('&')[1].strip()
			name1= name1.replace(' ','')
			match = re.match(r"([a-z]+)([0-9]+)",name1, re.I)
			items = match.groups()
			isotopes_1.append(items[0].capitalize()+'-'+items[1])

		#print isotopes_1

                #check if GCE isotopes are all available in simulation input
                #sefiles=se(self.runs_H5_surf[0])
                isotopes_2=[]
                no_iso=[]
                available_iso=[]
                for j in range(len(isotopes_1)):
                        if is_stable(isotopes_1[j])=='t':
                                isotopes_2.append(isotopes_1[j])
                for j in range(len(isotopes)):
                        #if isotope is not available stop       
                        if isotopes[j] not in isotopes_2:
                                print 'Isotope ',isotopes[j],'is not available'
                                no_iso.append(j)
                                #print 'STOP'
                        else:
                                available_iso.append(isotopes[j])
                print 'All isotopes for GCE available'
		
		label_done=False
		for i in range(len(lines_txt)):
			#header
			if 'specie' in lines_txt[i]:
				lins=lines_txt[i].split('&')[2:-1]
				print lins
				label1=[]
				idx_runs=[]
				for k in range(len(lins)):
					if exp_type in lins[k]:
						label1.append(lins[k])
						idx_runs.append(k)
				#print idx_runs
				label_done=True
				cols=[0]*len(label1)
				for k in range(len(cols)):
					cols[k]=[]
				continue
			if label_done==False:
				continue

			#if 'hline' in lines_txt[i] and i<15:
			#	continue
			#if 'hline' in lines_txt[i]:
			#	break

			#print lines_txt
			lins=lines_txt[i].split('&')[2:-1]
			j=0
			#print cols[j]
			for k in idx_runs:
				cols[j].append(float(lins[k].strip('\n\t').replace("\\", "")))
                		j+=1
		#print cols
		

                #labelout,cols,remn_masses_1 = self.write_prod_fact_stellarwinds(cycles=cycles,isotopes=available_iso,ascii_1=False,latex=False,GCE_input_fxt=True)


		#from paper1
                if label_masses == 'set1.1':
                        remn_masses=[12.349,14.120,14.232]
			metallicity_new=0.01
                if label_masses == 'set1.2':
                        metallicity_new=0.02
			remn_masses=[12.132,13.974,13.738,12.495,13.428]
                if label_masses == 'set1.3a':
                        metallicity_new=0.006
			remn_masses=[12.349,14.120,14.232]
                if label_masses== 'set1.5a':
                        metallicity_new=0.0001
			remn_masses=[12.349,14.120,14.232]

                #add values for missing isotopes, value=99
		print 'testteest'
		print no_iso
                for i in no_iso:
			print no_iso
                        for j in range(len(cols)):
				if i == no_iso[-1]:
				# here come the right remnant masses for massive stars
				#not just pre-collapse masses
					print 'adding normal value for rem'
					cols[j].insert(i,remn_masses[j])
				else:
                                	cols[j].insert(i,0)

                #metallicity1='{:.4E}'.format(float(metallicity))
                masses_new1=[]
		labelout=label1
                for i in range(len(labelout)):
                        mass=labelout[i].split('Msun')[0]
                        #mass='{:.4E}'.format(float(mass))
                        masses_new1.append(float(mass))
                #print 'Available input masses',masses

	
		yields_new=[0]*len(cols[0])
		masses_new=[0]*len(cols[0])
		label_new=[0]*len(cols[0])
		
		#for every mass
		for k in range(len(yields_new)):	
			yields_new[k]=[]
			masses_new[k]=masses_new1
			label_new[k]=[label_masses]*len(masses_new1)
			for w in range(len(cols)):
				yields_new[k].append(cols[w][k])

		#print metallicity_new
                #print yields_new
		#print label_new
		#print masses   

		##### Merge data###############################


		


		metallicities_out=len(isotopes)*[0]
		for t in range(len(isotopes)):
			pos_idx=0
			#assuming that if this value is still in it, fresh file from fxt without any new z
			#print metallicities,float('1.9000E-04')	
			if float('1.9000E-04') in metallicities[0]:
				metallicities_out[t]=metallicities[t][:2]+[metallicity_new]
				pos_idx=2
			else:
				metal_sum=[]
				st=0	
				for metal in metallicities[t]:
					#metallicities_out[t]=metallicities[t]
					if (metal>metallicity_new) and (metallicity_new not in metal_sum):
						metal_sum.append(metallicity_new)
						metal_sum.append(metallicities[t][st])
					else:
						metal_sum.append(metallicities[t][st])
					st+=1
				#in the case the added Z is larger than all available Z
				if (metallicity_new not in metal_sum):
					if len(metal_sum) == len(metallicities[t]):
						metal_sum.append(metallicity_new)
				metallicities_out[t]=metal_sum
				pos_idx=metallicities_out[t].index(metallicity_new)		
		
		masses_out=len(isotopes)*[0]
		yields_out=len(isotopes)*[0]
		label_out=len(isotopes)*[0]
		print metallicities_out[0]
		for t in range(len(isotopes)):
			masses_out[t]=len(metallicities_out)*[0]
			yields_out[t]=len(metallicities_out)*[0]
			label_out[t]=len(metallicities_out)*[0]
			w=0	
			for wt in range(len(metallicities_out[0])):
				#if case that it is a fresh file 
				if float('1.9000E-04') in metallicities[0] and wt <2:
					masses_out[t][wt]=masses[t][w]
					yields_out[t][wt]=yields[t][w]
					label_out[t][wt]=label[t][w]
					w+=1				
				else:
					#if we are a t the position to add new metallicity
					#or in primary case of this function to add yield values of exp
					#print isotopes[t]
					if wt==pos_idx:
						#print 'add new z or primary add yields to old Z'
					#	print t
						#number of to added masses different from already existing data
						yield_sum=[]
						for idx1 in range(len(masses[t][w])):
							if masses[t][w][idx1] in masses_new[t]:
								idx2=masses_new[t].index(masses[t][w][idx1])
								if not isotopes[t]=='rem':
									yield_sum.append( float(yields[t][w][idx1] + yields_new[t][idx2] ))
								else:
									yield_sum.append( float(yields_new[t][idx2] ))
							else:
									yield_sum.append(float(yields[t][w][idx1]))
						masses_out[t][wt]=masses[t][w]
						#print yield_sum
						yields_out[t][wt]=yield_sum
						label_out[t][wt]=label[t][w]
						#print wt,w
					else:
						#print metallicities_out
						#print 'else',wt,w
						#print masses_out[t],masses[t]
						masses_out[t][wt]=masses[t][w]
						yields_out[t][wt]=yields[t][w]
						label_out[t][wt]=label[t][w]
					w+=1


		###### Prepare output in right format######################
		#isotopes_out defined already above
		bb_massfrac_out=bb_massfrac
		type1a_yields_out=type1a_yields
		nova_yields_out=nova_yields
		#metallicities_out=metallicities
		#masses_out=masses
		#yields_out=yields
		#label_out=label 		
		newline=''
		for i in range(len(isotopes)):
			newline+= (isotopes_out[i].ljust(29)+'* symbol name'+'\n') 	
			newline+= ('  '+'{:.4E}'.format(bb_massfrac_out[i]).ljust(27)+'* big bang mass fraction'+'\n')
			newline+= (str(len(type1a_yields_out[i])).ljust(29)+'* number of type1a'+'\n')
			ia_label=['w7 tny86','sub chandra 1','sub chandra 2','sub chandra 7','sub chandra 8','nse 1a']
			for k in range(len(type1a_yields_out[i])):
				newline+= ('  '+'{:.4E}'.format(type1a_yields_out[i][k]).ljust(29)+ia_label[k]+'\n')
				newline+= (str(len(nova_yields_out[i])).ljust(29)+'* number of nova'+'\n')
				nova_label='* nova yield 1.0 1.25 & 1.35 msun'
				for k in range(len(nova_yields_out[i])):
					if k==0:
						newline+= ('  '+'{:.4E}'.format(nova_yields_out[i][k]).ljust(27)+nova_label+'\n')
					else:
						newline+= ('  '+'{:.4E}'.format(nova_yields_out[i][k]).ljust(27)+'\n')
						newline+= (str(len(metallicities_out[i])).ljust(29)+'* number of metallicity values'+'\n')
					for k in range(len(metallicities_out[i])):
						newline+= ('  '+'{:.4E}'.format(metallicities_out[i][k]).ljust(27)+'* metallicity'+'\n') 
						newline+= (str(len(masses_out[i][k])).ljust(29)+'* number of masses @ this metallicity'+'\n')
						for w in range(len(masses_out[i][k])):
							#the following is needed for the right label setting
							if w==0:
								newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+(7*' ')+label_out[i][k][w]+'\n')
								label_old=label_out[i][k][w]	
							elif label_out[i][k][w] != label_old:
								newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+(7*' ')+label_out[i][k][w]+'\n')	
								label_old=label_out[i][k][w]
							else:
								newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+'\n')



		###write out data

		#to fix bug where fxt GCe cannot deal with 99 yields
		newline=newline.replace('9.9000E+01','0.0000E+00')
		f2=open(output_file,'w')
		f2.write(newline)
		f2.close()



	def write_GCE_input_fxt(self,input_file='maltov1.orig',output_file='maltov1.new',first_change=True,label_masses='set1.1',cycles=20*[[0,-1,500]],exp_dir='',delay=False):

		'''
		Reads input_file of input for GCE code of fxt and replaces
		data with available data from the simulations. Also creates
		a file *output_file*_differences to keep track of the replacements. latter not yet
		'''



		import re
		f1=open(input_file)
		lines=f1.readlines()
		f1.close()
		isotopes=[]
		isotopes_out=[]
		iso_line_idx=[]
		i=0


		##########read input file##############
		#GCE isotopes
		for line in lines:
			if line[0] != ' ' and line[0:3] != 'rem':
				if line[0].isalpha():
					a= line[:6].split()[0]
					isotopes_out.append(a)
					match = re.match(r"([a-z]+)([0-9]+)",a, re.I)
					items = match.groups()
					isotopes.append(items[0].capitalize()+'-'+items[1])
					iso_line_idx.append(i)
			if line[0:3] == 'rem':
				isotopes_out.append('rem')
				isotopes.append('rem')
				iso_line_idx.append(i)
			i+=1

		bb_massfrac=[0]*len(isotopes)
		type1a_yields=[0]*len(isotopes)
		nova_yields=[0]*len(isotopes)
		metallicities=[0]*len(isotopes)
		masses=[0]*len(isotopes)
		yields=[0]*len(isotopes)
		label=[0]*len(isotopes)


		print isotopes,isotopes_out

		i=-1
		j=-1
		while (True):
			i+=1
			if i>len(lines)-1:
				break
			line=lines[i]
			print line
			if i in iso_line_idx:
				j+=1
			#read data
			if 'big bang' in line:
				bb_massfrac[j]=float(line[:13])
				continue
			#type 1a input ignore
			if 'number of type1a' in line:
				number=int(line[0])
				type1a_yields[j]=[0]*number
				k=0
				for sn in range(i+1,i+number+1):
					type1a_yields[j][k]=(float(lines[sn][:13]))				
					k+=1
				i=i+number	
				continue
			#nova input ignore
			if 'number of nova' in line:
				number=int(line[0])
				k=0
				nova_yields[j]=[0]*number
				for nova in range(i+1,i+number+1):
					nova_yields[j][k]=(float(lines[nova][:13])) 
					k+=1
				i=i+number   
				continue
			if 'number of metallicity values' in line:
				number_metallicities=int(line[0])
				metallicities[j]=[0]*(number_metallicities)
				masses[j]=[0]*(number_metallicities)
				yields[j]=[0]*(number_metallicities)
				label[j]=[0]*(number_metallicities)	
		 		z_counter=-1
		 		continue
		 	if '* metallicity' in line:
				z_counter+=1
		 		metallicities[j][z_counter]=float(line[:13])				
		 		continue
		 	if '* number of masses @ this metallicity' in line:
				number_masses=int(line[:2])
				masses[j][z_counter]=[0]*(number_masses)
				yields[j][z_counter]=[0]*(number_masses)
				label[j][z_counter]=[0]*(number_masses)
				k=0
				for mass in range(i+1,i+number_masses+1):
					masses[j][z_counter][k]=float(lines[mass][:13])
					yields[j][z_counter][k]=float(lines[mass][14:24])
					#check for label
					if len(lines[mass])>27:
						label1=lines[mass][27:].strip()
					label[j][z_counter][k] = label1
					k+=1
				i=i+number_masses
				continue
		#print bb_massfrac
		#print type1a_yields
		#print nova_yields
		#print '##########################################3'
		#print isotopes
		#print metallicities
		#print masses
		#print yields
		#print label	
		##############read new simulation input###########

		#check if GCE isotopes are all available in simulation input
		sefiles=se(self.runs_H5_surf[0])
		isotopes_1=sefiles.se.isotopes
		print isotopes_1
		isotopes_2=[]
		no_iso=[]
		available_iso=[]
		for j in range(len(isotopes_1)):
			if is_stable(isotopes_1[j])=='t':
				isotopes_2.append(isotopes_1[j])
		for j in range(len(isotopes)):
			#if isotope is not available stop       
			if isotopes[j] not in isotopes_2:
				print 'Isotope ',isotopes[j],'is not available'
				no_iso.append(j)
				#print 'STOP'
			else:
				available_iso.append(isotopes[j])
		#cycles=20*[[0,-1,500]]
		labelout,cols,remn_masses_1 = self.write_prod_fact_stellarwinds(cycles=cycles,isotopes=available_iso,ascii_1=False,latex=False,GCE_input_fxt=True,exp_dir=exp_dir,delay=delay)

                #to include remnants
                remn_masses=[]
                #remn_masses_1=20*[3]
                for i in range(len(remn_masses_1)):
                        #remnant_1='{:.4E}'.format(float(remn_masses_1[i]))
                        remn_masses.append(remn_masses_1[i])

                #add values for missing isotopes, value=99
		#print 'testteest'
		#print no_iso
                for i in no_iso:
                        for j in range(len(cols)):
				if i == no_iso[-1]:
				# to add the remnant masses
					print 'adding normal value for rem'
					cols[j].insert(i,remn_masses[j])
				else:
                                	cols[j].insert(i,0)

                #because the set runs have the same Z:
                #but first item in label is 'specie', so skip this
		#trick here, measure z only from first input file
		#massive stars are therefore not checked and treated as other z
                metallicity_new=float(labelout[1].split('=')[1])
                #metallicity1='{:.4E}'.format(float(metallicity))
                masses_new1=[]
                for i in range(1,len(labelout)):
                        mass=labelout[i].split('Msun')[0]
                        #mass='{:.4E}'.format(float(mass))
                        masses_new1.append(float(mass))
                #print 'Available input masses',masses

	
		yields_new=[0]*len(cols[0])
		masses_new=[0]*len(cols[0])
		label_new=[0]*len(cols[0])
		
		#for every mass
		for k in range(len(yields_new)):	
			yields_new[k]=[]
			masses_new[k]=masses_new1
			label_new[k]=[label_masses]*len(masses_new1)
			for w in range(len(cols)):
				yields_new[k].append(cols[w][k])

		#print metallicity_new
                #print yields_new
		#print label_new
		#print masses   

		##### Merge data###############################
		metallicities_out=len(isotopes)*[0]
		for t in range(len(isotopes)):
			pos_idx=0
			#assuming that if this value is still in it, fresh file from fxt without any new z
			#print metallicities,float('1.9000E-04')	
			if float('1.9000E-04') in metallicities[0]:
				metallicities_out[t]=metallicities[t][:2]+[metallicity_new]
				pos_idx=2
			else:
				metal_sum=[]
				st=0	
				print metallicities[t]
				print isotopes[t]
				for metal in metallicities[t]:
					#metallicities_out[t]=metallicities[t]
					if (metal>metallicity_new) and (metallicity_new not in metal_sum):
						metal_sum.append(metallicity_new)
						metal_sum.append(metallicities[t][st])
					else:
						metal_sum.append(metallicities[t][st])
					st+=1
				#in the case the added Z is larger than all available Z
				if len(metal_sum) == len(metallicities[t]):
					metal_sum.append(metallicity_new)
				metallicities_out[t]=metal_sum
				pos_idx=metallicities_out[t].index(metallicity_new)		
		
		masses_out=len(isotopes)*[0]
		yields_out=len(isotopes)*[0]
		label_out=len(isotopes)*[0]
		print metallicities_out[0]
		for t in range(len(isotopes)):
			#print isotopes[t]
			masses_out[t]=len(metallicities_out)*[0]
			yields_out[t]=len(metallicities_out)*[0]
			label_out[t]=len(metallicities_out)*[0]
			w=0	
			for wt in range(len(metallicities_out[0])):
				#if case that it is a fresh file 
				if float('1.9000E-04') in metallicities[0] and wt <2:
					masses_out[t][wt]=masses[t][w]
					yields_out[t][wt]=yields[t][w]
					label_out[t][wt]=label[t][w]
					w+=1				
				else:
					#if we are a t the position to add new metallicity
					if wt==pos_idx:
						#print 'add new z'
						masses_out[t][wt]=masses_new[t]
						yields_out[t][wt]=yields_new[t]
						label_out[t][wt]=label_new[t]
						#print wt,w
					else:
						#print metallicities_out
						#print 'else',wt,w
						#print masses_out[t],masses[t]
						masses_out[t][wt]=masses[t][w]
						yields_out[t][wt]=yields[t][w]
						label_out[t][wt]=label[t][w]
						w+=1


		###### Prepare output in right format######################
		#isotopes_out defined already above
		bb_massfrac_out=bb_massfrac
                type1a_yields_out=type1a_yields
                nova_yields_out=nova_yields
                #metallicities_out=metallicities
                #masses_out=masses
                #yields_out=yields
                #label_out=label 		
		newline=''
		for i in range(len(isotopes)):
			newline+= (isotopes_out[i].ljust(29)+'* symbol name'+'\n') 	
			newline+= ('  '+'{:.4E}'.format(bb_massfrac_out[i]).ljust(27)+'* big bang mass fraction'+'\n')
			newline+= (str(len(type1a_yields_out[i])).ljust(29)+'* number of type1a'+'\n')
			ia_label=['w7 tny86','sub chandra 1','sub chandra 2','sub chandra 7','sub chandra 8','nse 1a']
			for k in range(len(type1a_yields_out[i])):
				newline+= ('  '+'{:.4E}'.format(type1a_yields_out[i][k]).ljust(29)+ia_label[k]+'\n')
			newline+= (str(len(nova_yields_out[i])).ljust(29)+'* number of nova'+'\n')
			nova_label='* nova yield 1.0 1.25 & 1.35 msun'
			for k in range(len(nova_yields_out[i])):
				if k==0:
					newline+= ('  '+'{:.4E}'.format(nova_yields_out[i][k]).ljust(27)+nova_label+'\n')
				else:
					newline+= ('  '+'{:.4E}'.format(nova_yields_out[i][k]).ljust(27)+'\n')
			newline+= (str(len(metallicities_out[i])).ljust(29)+'* number of metallicity values'+'\n')
			for k in range(len(metallicities_out[i])):
				newline+= ('  '+'{:.4E}'.format(metallicities_out[i][k]).ljust(27)+'* metallicity'+'\n') 
				newline+= (str(len(masses_out[i][k])).ljust(29)+'* number of masses @ this metallicity'+'\n')
				for w in range(len(masses_out[i][k])):
					#the following is needed for the right label setting
					if w==0:
						newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+(7*' ')+label_out[i][k][w]+'\n')
						label_old=label_out[i][k][w]	
					elif label_out[i][k][w] != label_old:
						newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+(7*' ')+label_out[i][k][w]+'\n')	
						label_old=label_out[i][k][w]
					else:
						newline+= ('  '+'{:.4E}'.format(masses_out[i][k][w])+'  '+'{:.4E}'.format(yields_out[i][k][w])+'\n')
			


		###write out data

		f2=open(output_file,'w')
		f2.write(newline)
		f2.close()




        def write_GCE_input_fxt_OLD(self,input_file='maltov1.orig',output_file='maltov1.new',first_change=True,label_masses='set1.2',cycles=20*[[0,-1,500]]):

                '''
                        Reads input_file of input for GCE code of fxt and replaces
                        data with available data from the simulations. Also creates
                        a file *output_file*_differences to keep track of the replacements.
                '''

		import re
		f1=open(input_file)
		lines=f1.readlines()
		f1.close()
		isotopes=[]
		iso_line_idx=[]
		i=0
		#GCE isotopes
		for line in lines:
			if line[0] != ' ' and line[0:3] != 'rem':
				if line[0].isalpha():
					a= line[:6].split()[0]
					match = re.match(r"([a-z]+)([0-9]+)",a, re.I)
					items = match.groups()
					isotopes.append(items[0].capitalize()+'-'+items[1])
					iso_line_idx.append(i)
			if line[0:3] == 'rem':
                		isotopes.append('rem')
                		iso_line_idx.append(i)
			i+=1
		
		print isotopes
		print iso_line_idx
		#isotopes=['H-1', 'H-2', 'He-3', 'He-4', 'Li-6', 'Li-7']

		#check if GCE isotopes are all available in simulation input
		sefiles=se(self.runs_H5_surf[0])
		isotopes_1=sefiles.se.isotopes
		isotopes_2=[]
		no_iso=[]
		available_iso=[]
		for j in range(len(isotopes_1)):
			if is_stable(isotopes_1[j])=='t':
				isotopes_2.append(isotopes_1[j])
		print isotopes_2
		for j in range(len(isotopes)):
			#if isotope is not available stop	
			if isotopes[j] not in isotopes_2:
				#print 'Isotope ',isotopes[j],'is not available'
				no_iso.append(j)
				#print 'STOP'
			else:
				available_iso.append(isotopes[j])
		print 'All isotopes for GCE available'
		#cycles=20*[[0,-1,500]]
		label,cols,remn_masses_1 = self.write_prod_fact_stellarwinds(cycles=cycles,isotopes=available_iso,ascii_1=False,latex=False,GCE_input_fxt=True)
		

		#to include remnants
		remn_masses=[]
		#remn_masses_1=20*[3]
		for i in range(len(remn_masses_1)):
		        remnant_1='{:.4E}'.format(float(remn_masses_1[i]))
        		remn_masses.append(remnant_1)


		
		#add values for missing isotopes, value=99
		for i in no_iso:
			for j in range(len(cols)):
				cols[j].insert(i,0)		
		print isotopes
		print cols
		#because the set runs have the same Z:
		#but first item in label is 'specie', so skip this
		metallicity=label[-1].split('=')[1]
		metallicity='{:.4E}'.format(float(metallicity))
		masses=[]
		for i in range(1,len(label)):
			mass=label[i].split('Msun')[0]
			mass='{:.4E}'.format(float(mass))
			masses.append(mass)
		print 'Available input masses',masses

		#Write them into output_file

		print 'Writing in file ',output_file
				
		i=-1
		j=-1
		new_lines=''
		new_data_added=False	
		while (True):
			i+=1
			if i>len(lines)-1:
				break
			#if 200<i:
			#	break
			line=lines[i]
			#for each isotope
			if i in iso_line_idx:
				#if new data was not added for last isotope
				#jump back to last metallicity line and add new metallicity
				if new_data_added==False and i>iso_line_idx[0]:
					i=last_metallicity_line
					need_to_add_z=True
					continue
				j+=1
				#if j>1:
				 #      break
				#if models do not have isotope needed, take one from fxt
				if j in no_iso and iso_line_idx[-1]>i:
					print 'Isotope data for isotope',isotopes[j],'taken from fxt'
					for old_input_line in range(i,149+i):
						new_lines+=lines[old_input_line]
					print 'jump to next isotope'
					print 'j:',j
					i=iso_line_idx[j+1]-1
					#break
					continue
				print 'isotope',isotopes[j]
				print 'set new data false'
				new_data_added=False
				need_to_add_z=False

				to_skip_1=0
				to_skip_2=0
				to_skip_3=0
				label_masses_done=False
				new_lines+=line
				continue
			if 'big bang' in line:
				new_lines+=line
				continue
			#type 1a input ignore
			if 'type1a' in line:
				to_skip_1=int(line[0])

				new_lines+=line
				continue
			if to_skip_1>0:
				to_skip_1-=1
				new_lines+=line
				continue
			#nova input ignore
			if 'number of nova' in line:
				to_skip_2=int(line[0])
				new_lines+=line
				continue
			if to_skip_2>0:
				to_skip_2-=1
				new_lines+=line
				continue
			#input from models
			#how many metallicities? 2 standard values from fxt + 
			if 'number of metallicity values' in line:
				#z_value=0
				#check if z is already available:
				following_z_array=[]
				#the following if statement because iso_line_idx does not contain
				#iso after 'remn', no further index j
				if iso_line_idx[-1]<i:
					max1=len(lines)
				else:
					max1=iso_line_idx[j+1]
				for t in range(i,max1):
					if '* metallicity' in lines[t]: 
						following_z_array.append(lines[t][2:12])
				#if z is already there
				print metallicity,following_z_array
				if metallicity in following_z_array:
					new_lines+=line
					print 'Z found'
					z_values=0
				#else add new z	
				else:
					print 'add +1 to metal number for element'
					new_metal_number=int(line[:5]) +1
					new_lines+=str(new_metal_number)+line[len(str(new_metal_number)):]
					z_values=int(line[0])

                                #if first_change == True:
                                z_values=3
				continue
		       #go over metallicities
			if '* metallicity' in line:
				last_metallicity_line=i
				z_values-=1
				skip_masses=0
				add_own_masses=False
				if z_values>0:
					new_lines+=line
				else:
					#check if available Z is larger then Z from input simulation
					print float(metallicity),float(line[2:12]),'teststest'
					print line[2:12].strip()
					if (float(metallicity)< float(line[2:12]) and (new_data_added==False)) or need_to_add_z==True:
						       #add own Z
							print 'adding new metallicity'
							new_lines+='  '+metallicity+line[12:]
							new_data_added=True
							line_replaced=i
					#if z is already there, compare strings
					elif metallicity == line[2:12].strip() and (new_data_added==False):
						print 'location of z found'
						add_own_masses=True
						new_lines+=line
						new_data_added=True
					else:
						#add massive stars, M>=8
						new_lines+=line
						if first_change==False:
							skip_masses=11
						else:
							skip_masses=11
						#in remnant case
						print 'skipping masses'
				continue
			#go over masses
			if 'number of masses @ this metallicity' in line:
				#for available data
				masses_available=0
				if z_values>0:
					to_skip_3=int(line[0:5])
					new_lines+=line
					continue
				#for new input
				else:
					if skip_masses>0:
						once_done=False
						masses_available=int(line[:2])
						#when reaching the last 'isotope' remn:
						print masses_available,skip_masses
						if iso_line_idx[-1]<i:
							new_lines+=line
						else:
							if skip_masses==15:
								new_lines+=line
							else:
								if masses_available<=15:
									new_lines+=line
									once_done=True
								else:
									new_lines+=str( masses_available-skip_masses+1 )+line[len(str(masses_available-skip_masses+1)):]
						continue
					else:
						#if z is already there
						if add_own_masses==True:
							new_sum_masses=len(masses)+int(line[:2])
							new_lines+=str(new_sum_masses)+'  '+line[len(str(new_sum_masses))+2:]
							mass_counter=-int(line[:2])
						else:
							print 'add new Z'
							new_lines+=str(len(masses))+'  '+line[len(str(len(masses)))+2:]
							mass_counter=0
						continue
			#to skip due to available data
			if to_skip_3>0:
				to_skip_3-=1
				new_lines+=line
				continue
			if masses_available>0:
				#print 'this is executed too many times'
				masses_available-=1
				skip_masses-=1
				if iso_line_idx[-1]<i:
					new_lines+=line
					continue
				if once_done==True:
					new_lines+=line
					continue
				if skip_masses<1:
					new_lines+=line
				if first_change==False:
					new_lines+=line
				continue
			print 'this line is not executed...'
			#add new masses
			if mass_counter < (len(masses)):
				#if z is already there and available masses needs to be skipped to add new later on
                                if add_own_masses==True and mass_counter<0:
					mass_counter+=1
					new_lines+=line
					i-=1
					continue
				#if create new z or add new masses to available z
				#when reaching the last 'isotope' remn:
				if iso_line_idx[-1]<i:
					yield_1=remn_masses[mass_counter]
				#for all real isotopes:
				else:
					yield_1 = cols[mass_counter][j]

				yield_line='  '+masses[mass_counter]+'  '+'{:.4E}'.format(float(yield_1))
				if label_masses_done==False:
					new_lines+=yield_line+7*' '+label_masses+'\n'
					label_masses_done=True
				else:
					new_lines+=yield_line+'\n'
				mass_counter+=1
				continue
			# jump to line with last metallicity which was replaced by new data
			elif mass_counter == (len(masses)):
				print 'jump to original'
				if need_to_add_z==False:
					i=line_replaced-1
					mass_counter+=1
				continue
			#all input for isotope j done, jump to line where next isotope starts

			print 'jump to next isotope'
			print 'j:',j
			#i=iso_line_idx[j+1]-1
			#for rem_line in range(11771,11920):
			#       new_lines+=lines[rem_line]
			#       print new_lines
        		break

#write output file
		f2=open(output_file,'w')
		f2.write(new_lines)
		f2.close()		


	
	
	#the followign routine needs to be adapted..
	def write_prod_fact_exp(self,cycles=[],isotopes=['H-1','C-12'],elements=[],delay=True,yields_output=True,headers=['this is a test','to write strings in the first col'],file_name="isotopic_table_set1.2_prodfac.txt",ascii_1=False,latex=False,GCE_input_fxt=True,GCE_input=False):

		'''
	
			Calculates the yields or production factors of isotopes or elements and writes them as
			- an ascii output
			- latex output
			- GCE input file , this is used by the write_GCE_input.... functions above

		
		'''
		import ascii_table as ascii1

		production_factors=[]
		case_label_and_first_column = ['specie']
		yields_all=[]
		iso_names=[]
		remn_masses=[]
		for i in range(len(self.runs_H5_surf)):
                	sefiles=se(self.runs_H5_surf[i])
			sefiles_hout=se(self.runs_H5_out[i])
			indx_iso=[]
			yields_1=[]
			iso_names=[]
			mass=sefiles.get("mini")
			metallicity=sefiles.get("zini")
			case_label_and_first_column.append(str(mass)+" Msun, Z="+str(metallicity))
			##if all stable ones:
			if isotopes[0]=="all":
				isotopes_1=sefiles.se.isotopes
				isotopes=[]
				for j in range(len(isotopes_1)):
					if is_stable(isotopes_1[j])=='t':
						isotopes.append(isotopes_1[j])			
			print isotopes		
			if len(elements)>0:
				if elements[0]=="all":
					elements=sefiles.se.elements
					elements=get_stable(sefiles.se.isotopes,get_elements=True)[:12]		
					#exclude neutrons
					print elements
                        if cycles[i][1]==-1:
                        	endcycle=int(sefiles.se.cycles[-1]) #+ cycles[k][2]  #1000
                                print endcycle
                               	endcycle = int(endcycle) // 100 * 100
				endcycle = int(endcycle-20)
				#endcycle= int( round( endcycle,-3) -20 )
                                print "endcycle",endcycle
                        else:
                                endcycle=cycles[i][1]
			
			import nugridse as mp
			reload(mp)
			prod_facs,isotopes,yields,iniabus,mass_cut,first_cycle=self.weighted_yields_explosion(mp,sefiles,sefiles_hout,isotopes=isotopes,delay=delay)
			#print 'isotopes',isotopes_prod_fac
			#print 'prod factor:',prod_factor
			#print 'yield:',yields
			
		for k in range(len(yields)):
			yields[k]='{:.3E}'.format(yields[k])			
	
		return 	isotopes,yields,mass_cut
	
	#new function to combine ability to produce NUPYCEE tables and create plots	
	def get_stellar_ejecta(self,cycles=[],isotopes=["C-12","N-14","O-16","Ne-22"],elements=[],yields_output=True,file_name="isotopic_table_set1.2_prodfac.txt",GCE_tables=True,plots=False,plots_singleZ=False,plots_singleA=False,exp_dir='',exp_type=None,fallback_masses=[],fallback_coords=[],pre_exp=False,exp_only=False,wind_only=False,plot_masses_only=[],exp_cycle_inp=-1,plot_log=True,color=['r'],marker_type=['o'],line_style=['--'],markersize=[6],line_width=[14],title='',withlabel=True,withlabelelement=False,label='',plot_lines=True,fig=0,labelZonly=False,label_isotopes=[],plot_o16_lines=False):

		'''
	
			Calculates the yields or overproduction factors of isotopes or elements and writes them as
			1)  NUPYCEE tables if GCE_tables=False			
			2)  plots if plots=True
			
			yields_output: if true, plots or writes yields, else writes overproduction factors
			species chosen: default uses isotopes array, if elements set uses elements array
			###Currently works not for set1.2 for exp_only=True and m_remn too large (ppd_exp run dir not available)
		
		'''

		self.fallback_masses=fallback_masses
		self.fallback_coords=fallback_coords
		

		if (GCE_tables==True and  plots==True) or (GCE_tables==True and (plots_singleA==True or  plots_singleZ==True)):
			print 'Warning: GCE_tables==True and  plots==True'
			print 'Choose one option only'
			return 0

		import ascii_table as ascii1


                case_label_and_first_column = ['specie']
		case_label_and_first_column_tex = ['specie']
                yields_all=[]
                remn_masses=[]

		HDF5_surf=[]
		HDF5_out=[]
		HDF5_restart=[]
		HDF5_surf=self.runs_H5_surf
		HDF5_out=self.runs_H5_out
		HDF5_restart=self.runs_H5_restart
		
		legend_k=0
		sefiles=[]
		sefiles_hout=[]
		sefiles_restart=[]
		for i in range(len(HDF5_surf)):
				sefiles.append(se(HDF5_surf[i]))		
				sefiles_hout.append(se(HDF5_out[i]))
				sefiles_restart.append(se(HDF5_restart[i]))	
		isotopes11=sefiles[0].se.isotopes
		if (len(isotopes)>0) and (len(elements)==0):
			if isotopes[0]=='allstable':
				isotopes=get_stable(isotopes11,get_elements=False)
				#color=len(isotopes)*[color[0]]
				#marker_type=len(isotopes)*[marker_type[0]]
				#line_style=len(isotopes)*[line_style[0]]
				#markersize=len(isotopes)*[markersize[0]]
				#line_width=len(isotopes)*[line_width[0]]
				print isotopes
		else: #len(elements)>0:
			if elements[0]=='allstable':
				elements=get_stable(isotopes11,get_elements=True)
				#color=len(elements)*[color[0]]
				#marker_type=len(elements)*[marker_type[0]]
				#line_style=len(elements)*[line_style[0]]
				#markersize=len(elements)*[markersize[0]]
				#line_width=len(elements)*[line_width[0]]
				print elements
				isotopes=get_stable(isotopes11,get_elements=False)
			#for plotting prodfac vs Z
                	z_elements=[]
                	for k in range(len(elements)):
                	        z_elements.append(u.get_z_from_el(elements[k]))
		a_isotopes=[]
                for k in range(len(isotopes)):
                        a_isotopes.append(int(isotopes[k].split('-')[1]))
	

		####dealign with explosion pinput
                #test to take right dir:
                #
                # use exp_type variable from now on
                #if exp_type == None:    #to be backwards compatible
                #        if delay:
                #                exp_type='delay'
                #        else:
                #                exp_type='rapid'
		#relevant for weighted  yields expl funciton below
		delay=True
		if exp_type == 'delay':
			delay=True
		if exp_type == 'rapid':
			delay=False
		#else exp_type=='' is ye mcut/user input
                slist = os.listdir(exp_dir)
                expr = re.compile(exp_type)
                slist=(filter(expr.search,slist))


		exp_runs_H5_restart=[]
		sefiles_exp=[]


		#to dinstiungish between pre exp and exp sources
		if pre_exp==False:

			slist = os.listdir(exp_dir)
			expr = re.compile(exp_type)
			slist=(filter(expr.search,slist))

			for element in slist:
				run_path=exp_dir+'/'+element
				if not os.path.isdir(run_path):
					continue
				if os.path.isdir(run_path+"/H5_restart"):
					sefiles1 = os.listdir(run_path+"/H5_restart")
					if (filter(expr.search,sefiles1)) <1:
						print "Warning: No hdf5 restart files found in "+run_path+"/H5_restart"
					else:
						exp_runs_H5_restart.append(run_path+"/H5_restart")
						sefiles_exp.append(se(run_path+"/H5_restart"))
		else:
	
			exp_runs_H5_restart=HDF5_restart #H5_out if no restart available
			for k in range(len(HDF5_restart)):
				sefiles_exp.append(se(HDF5_restart[k]))


		z_index_files=[]
		z_values=[]
		j=-1
		print 'HDF5_surf'
		print HDF5_surf
		for i in range(len(HDF5_surf)):
			j+=1
			if type(sefiles[i].get("mini")) == np.ndarray:
				star_z=sefiles[i].get("zini")[0]
			else:
				star_z=sefiles[i].get("zini")
			if star_z not in z_values:
				z_values.append(star_z)
				z_index_files.append([])
			z_index_files[ z_values.index(star_z)].append(i)
		max_yield=[]
		color_iso=-1
		yields_1=[]
		t=0
		if len(elements)>0:
			isotopes=elements
		legend_k=0
		for w in range(len(z_index_files)):
			star_mass_array=[]
			#yields=[]
			legend_k+=1
			iso_yields=[]
			iniabu_yields_folded=[]
			yields_all_masses=[]
			exclude_masses=[]
			for i in range(len(isotopes)):
				if exp_only==False:
					yields_all_masses.append(np.zeros(len( z_index_files[w]   )))
				#	yields.append(np.zeros(len( z_index_files[w]   )))
				else:
					#yields_all_masses.append([])
					yields_all_masses.append(np.zeros(len( z_index_files[w]   )))
				#	yields.append([])
			ttt=0
			### loop over masses
			plot_k=0
			for k in z_index_files[w]:
				if type(sefiles[k].get("mini")) == np.ndarray:
					star_mass=sefiles[k].get("mini")[0]
					star_z=sefiles[k].get("zini")[0]
				else:
                                        star_mass=sefiles[k].get("mini")
                                        star_z=sefiles[k].get("zini")
					
				#star_mass=sefiles[k].get("mini")[0]
				#star_z=sefiles[k].get("zini")[0]
				cycs=self.get_last_cycle(cycles[k],sefiles[k],sefiles_restart[k])

				'''
				if cycles[k][1]==-1:
					#if star_mass>9:
					endcycle=int(sefiles_restart[k].se.cycles[-2])
					#else:
					#	endcycle=int(sefiles[k].se.cycles[-2])
				else:
					endcycle=cycles[k][1]
				cycs=range(cycles[k][0],endcycle,cycles[k][2])
				if endcycle not in cycs:
					cycs=cycs+[endcycle]
				print 'first, last cycle: ',cycles[k][0],endcycle
				star_mass11=sefiles[k].get(cycs,'mass')				
				#do a check for missing cycles
				if not len(star_mass11)==len(cycs):
					cycs1=[int(cycs[0])]
					for cyc in cycs[1:]:
						if  cyc <= max(cycs1):
							continue
						a=sefiles[k].get(cyc,'mass')
						cyctest=int(cyc)
						#if cycle 0 or a non-existing cycle, or cycle was alrady added, 
						#go one cycle further   
						while ((type(a)==list) or (a==0)):
							if   cyctest<int(cycs[-1]):
								print int(cycs[-1])
								print 'found bad cycle',cyctest
								cyctest+=1
								a=sefiles[k].get(cyctest,'mass')
							elif cyctest==int(cycs[-1]):
								cyctest=int(sefiles_restart[k].se.cycles[-3])
								a=sefiles[k].get(cyctest,'mass')	
								h=-4
								while ((type(a)==list) or (a==0)):
									cyctest=int(sefiles_restart[k].se.cycles[h])
									a=sefiles[k].get(cyctest,'mass')
									h-=1
								for kk in cycs1:
									if kk>=cyctest:
										cycs1=cycs1[:kk]
								break
							else:
								break
						print 'Use cycle ',cyctest
						if cyctest>max(cycs1):
							cycs1.append(cyctest)
					print len(cycs),len(cycs1)
					#test if last cycle is also present in 
					#outlastcyc=int(sefiles_hout[k].se.cycles[-1])
					#if not cyctest == outlastcyc:
						#testout = sefiles_hout[k].get(cyctest,'mass')
						#if ((type(testout)==list) or (testout==0))
					cycs=cycs1
				print cycs[-5:]
				'''

				#decide if wind contribution needs to be taken into account	
				if exp_only==False:
					case_label_and_first_column_tex.append(str(star_mass)+" Msun")
					case_label_and_first_column.append(str(star_mass)+" Msun, Z="+str(star_z))
					star_mass_array.append(star_mass)
					prod_factor,isotopes_prod_fac,yields,iniabu,mass_frac_ini,remn_mass =self.weighted_yields(sefiles[k],sefiles_restart[k],isotopes,elements,cycs)
					for tt in range(len(prod_factor)):
						prod_factor[tt]=yields[tt]/iniabu[tt]
					#print 'yield output###################################:'
					#print yields,prod_factor,iniabu
				else:
					prod_factor=[]
					yields=[]
                                        if star_mass<=8:
                                                continue
					#just to get mass_frac_ini...
                                        prod_factor_dum,isotopes_prod_fac_dum,yields_dum,iniabu_dum,mass_frac_ini,remn_mass_dum =self.weighted_yields(sefiles[k],sefiles_restart[k],isotopes,elements,cycs)
					
                                        #for pp in range(len(isotopes)):
                                        #        yields_all_masses[pp].append([0])
						#yields[pp].append([0])
					case_label_and_first_column_tex.append(str(star_mass)+" Msun")
					case_label_and_first_column.append(str(star_mass)+" Msun, Z="+str(star_z))

                                        star_mass_array.append(star_mass)

					
				mass=star_mass
				metallicity=star_z
				print 'test for stellar mass ',mass,wind_only
				if (mass>8) and (wind_only==False):

					#in case no ppn_exp run dir is available, skip explosion contribution
					#this is because mass cut is then larger then actual mass of star...no mass lost
					print 'too bad length is zero ',len(exp_runs_H5_restart)
					for t in range(len(exp_runs_H5_restart)):
						if type(sefiles_exp[t].get("mini")) == np.ndarray:
							mass1=sefiles_exp[t].get("mini")[0]
							metallicity1=sefiles_exp[t].get("zini")[0]
						else:
							mass1=sefiles_exp[t].get("mini")
							metallicity1=sefiles_exp[t].get("zini")

						###for currupted files

						corrupt=False
						if type(mass1) == list:
							corrupt=True
						if type(metallicity1) ==list:
							corrupt=True
						if ((float(mass1) < 1) or  (float(mass1) >70 )) or type(mass1) == type(None):
							corrupt=True
						if ((float(metallicity1) < 0.00001) or  (float(metallicity1) >1. )) or type(metallicity1) == type(None):
							corrupt=True


						if corrupt == True:
							metallicity1=0.02
							rundir11=exp_runs_H5_restart[t].split('/')[-2]
							mass1=float(rundir11.split('.')[0][1:])
							#print 'corruption, assinge new values',mass1,metallicity1
						if mass1 == mass and metallicity1 == metallicity:
							print 'found Z-M pair'
							#sefiles_re2=sefiles_exp[t]
							#calculate exp part
							import nugridse as mp
							reload(mp)
							#Have to find the right cycle in restart file
							#.se.cycles does not work when reading all files	
							#therefore identify file with last cycle and read it:
							refiles=sefiles_exp[t].se.files
							tosort=[]
							refiles1=[]
							for i in range(len(refiles)):
								cycle_file=refiles[i][-18:-11]
								if cycle_file.isdigit():
									tosort.append(int(cycle_file))
									refiles1.append(refiles[i])
							idx_sorted=sorted(range(len(tosort)), key=lambda k: tosort[k])
							idx_file = idx_sorted[-1]
							sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_file])
							if len(sefiles_re_cycle.se.cycles)==0:
								#in the case there is not cycle in the last file...set1.2 M2-
								sefiles_re_cycle=mp.se(exp_runs_H5_restart[t],refiles1[idx_sorted[-2]])
							if pre_exp==True:
								cycleend=cycs[-1]
							else:
								cycleend=sefiles_re_cycle.se.cycles[-1]	
                                                        print 'tewststest',exp_cycle_inp
                                                        if exp_cycle_inp>0:
                                                                sefiles_re_cycle=mp.se(exp_runs_H5_restart[t])		
								cycleend=exp_cycle_inp
							#cycleend=-1
                                                        print 'Use custom exp cycle cycleend',cycleend
							print 'for explosion profile use ',sefiles_re_cycle.se.filename
							print 'make sure that cycles you request are avilable'
							#print 'use file ',sefiles[k].se.filename
							print 'mode ',delay
							#print 'identifies as ',exp_runs_H5_restart[t],'with file',refiles1[idx_file]
							isotopes_exp,yields_exp,iniabu_exp,mass_cut_exp,first_cycle_exp=self.weighted_yields_explosion(mp,sefiles_re_cycle,sefiles[k],isotopes=isotopes,elements=elements,cycleend=cycleend,delay=delay)

							#if pre-exp prescription (H5_surf from ppd_wind) but remnant larger than total mass.
							if len(isotopes_exp)==0:
								print 'For M=',mass,'Z=',metallicity,' skip because mtotal<_remn'
								continue

							if exp_only==False:	
								#add pre-exp yields
								for tt in range(len(prod_factor)):
									#print 'star above 8M'
									#print isotopes_prod_fac[tt],yields[tt],yields_exp[tt]
									#the iniabu from the explosion is not added
									#because it is already in the iniabu from wind
									prod_factor[tt]=(yields[tt]+yields_exp[tt])/(iniabu[tt]+iniabu_exp[tt])
									#print 'for isotope',isotopes[tt],yields[tt],yields_exp[tt]
									yields[tt]=yields[tt]+yields_exp[tt]
									#print 'yields wind:',yields[tt]
									#print 'prodfac wind:',(yields[tt]/iniabu[tt])
									#print 'wind+exp yields : ',yields[tt]+yields_exp[tt]	
							else:
								for tt in range(len(isotopes_exp)):		
									prod_factor.append((yields_exp[tt])/(iniabu_exp[tt]))
									yields.append(yields_exp[tt])
							remn_mass=mass_cut_exp
							#print 'exp yield',yields_exp
							#print 'prodfac exp:',(yields_exp[tt]/iniabu_exp[tt])
							#print 'wind+exp prodfac:',prod_factor[tt]
							break		
				#### Specific  star ####

				#when exp_only without wind but exp run does not exist (BH)
				if len(prod_factor)==0:
					print 'zero production factors'	
					print 'when exp_only without wind but exp run does not exist (BH)'
					exclude_masses.append(star_mass)
					continue
				#yields_all carrying all info
				remn_masses.append(remn_mass)	
				if yields_output==False:	
					yields_all.append(list(prod_factor))
					#if plots vs mass are requested
					if plots==True:
						for i in range(len(isotopes)):
							yields_all_masses[i][ttt]=prod_factor[i]
						ttt+=1
				else:
					yields_all.append(list(yields))
					#if plots vs mass are requested
                                        if plots==True:
                                                for i in range(len(isotopes)):
                                                        yields_all_masses[i][ttt]=yields[i]
						ttt+=1

				####to write out table for specific star for  SYGMA and OMEGA
				if GCE_tables==True:
					#needs mass_frac_ini,m_final
					#m_final=round(remn_mass,2)
					#print '                       ',mass_frac_ini
					#mass_frac_ini=mass_frac_ini
					dcols=['Isotopes','Yields','X0']
					finalmheader='Mfinal: '+'{:.3E}'.format(remn_mass)
					special_header='Table: (M='+str(mass)+',Z='+str(metallicity)+')'
					data=[isotopes,yields_all[-1],mass_frac_ini]
					ascii1.writeGCE_table(filename=file_name,headers=[special_header,finalmheader],data=data,dcols=dcols)	
				####for each star plot overprodfac vs Z
				if plots_singleZ:
					if (not mass in plot_masses_only) and len(plot_masses_only)>0:
						continue
					#plot for every star prodfac vs A.
					if fig ==0:
						figname=str(mass)+'_'+str(star_z)
					else:
						figname=fig
					#print 'test output##########'
					#print prod_factor
					#print elements
					plt.figure(figname)
					label=str(mass)+' M$_{\odot}$'
					if labelZonly==True:
						label='Z='+str(star_z)
					idx_sorted=sorted(range(len(z_elements)),key=lambda x:z_elements[x])
					z_elements_sorted=[]
					prod_factor_sorted=[]
					elements_sorted=[]
					for idxx in idx_sorted:
						z_elements_sorted.append(z_elements[idxx])
						prod_factor_sorted.append(prod_factor[idxx])
						elements_sorted.append(elements[idxx])
					marker='o'
					color1='k'
					linestyle='-'
					if len(marker_type)>0:
						marker=marker_type[plot_k]
					if len(color)>0:
						color1=color[plot_k]
					if len(line_style)>0:
						linestyle=line_style[plot_k]
					plt.plot(z_elements_sorted,prod_factor_sorted,label=label,marker=marker,color=color1,linestyle=linestyle)			 
					plot_k = plot_k+1
					plt.xlabel('Charge Z')
					plt.ylabel('Produ..')
					ax=plt.gca()
					hh=0
					import matplotlib.transforms as transforms
					trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
					#print z_elements_sorted
					#print prod_factor_sorted
					#print elements_sorted
					for z in range(len(prod_factor)):
						pos=0.90
						if hh==0:
							nmin=pos#max(prod_factor)#0.2
						if hh==1:
							nmin=pos+0.02#max(prod_factor)+0.05#0.25
						if hh==2:
							nmin=pos+0.04#max(prod_factor)+#0.3
						fontsizelabel='x-small'
						ax.text(z_elements_sorted[z],nmin,elements_sorted[z],horizontalalignment='center',verticalalignment='center',\
						fontsize=fontsizelabel,clip_on=True,transform=trans) #x-msall
						if hh==0:
							hh=1
							continue
						if hh==1:
							hh=2
							continue
						if hh==2:
							hh=0
							continue
					m_max=86
					#if plot_lines==True:
					#	plt.plot([0,m_max],[1,1],"k-",linewidth=3)
						#plt.plot([0,m_max],[2,2],"k--",linewidth=1)
						#plt.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)
					
				       #set the plot x axis range
					#if arange[1]>0:
					#	plt.xlim(arange[0],arange[1])
					#else:
					plt.xlim(0,84)
					plt.ylim(1,(max(prod_factor)+0.1))
					plt.yscale('log')
					plt.minorticks_on()
					if yields_output==False:
						plt.ylabel("Overproduction factor",fontsize=16)
					else:
						plt.ylabel("Yields/M$_{\odot}$",fontsize=16)
					plt.legend()

                                ####for each star plot overprodfac vs A
				if plots_singleA:
                                        if (not mass in plot_masses_only) and len(plot_masses_only)>0:
                                                continue
					#plot for every star prodfac vs A.
					if fig ==0:
						figname=str(mass)+'_'+str(star_z)
					else:
						figname=fig
					plt.figure(figname)
					if labelZonly==True:
						label='Z='+str(star_z)
					else:
						label=str(mass)+' M$_{\odot}$'
					#plt.plot(a_isotopes,prod_factor,label=label,marker='o',color='k')
					plt.xlabel('Mass number A')
					ax=plt.gca()
					label_done=[]
					idx_min=0
					hh=0
					idx_sorted=sorted(range(len(a_isotopes)),key=lambda x:a_isotopes[x])
					a_isotopes_sorted=[]
					isotopes_sorted=[]
					prod_factor_sorted=[]

					isotopes_sorted_group=[]
					prod_factor_sorted_group=[]
					element_group=[]
					a_isotopes_sorted_group=[]

					for idxx in idx_sorted:
						a_isotopes_sorted.append(a_isotopes[idxx])
						prod_factor_sorted.append(prod_factor[idxx])
						isotopes_sorted.append(isotopes[idxx])
						ele1=isotopes[idxx].split('-')[0]		
						if ele1 in element_group:
							idx1=element_group.index(ele1)
							prod_factor_sorted_group[idx1].append(prod_factor[idxx])
							a_isotopes_sorted_group[idx1].append(a_isotopes[idxx])
							isotopes_sorted_group[idx1].append(isotopes[idxx])
						else:
							element_group.append(ele1)
							prod_factor_sorted_group.append([])
							a_isotopes_sorted_group.append([])
							isotopes_sorted_group.append([])
							prod_factor_sorted_group[-1].append(prod_factor[idxx])
							a_isotopes_sorted_group[-1].append(a_isotopes[idxx])
							isotopes_sorted_group[-1].append(isotopes[idxx])
					marker='o'
					color1='k'
					linestyle='-'
					if len(marker_type)>0:
						marker=marker_type[plot_k]
					if len(color)>0:
						color1=color[k]
					if len(line_style)>0:
						linestyle=line_style[plot_k]
					plot_k = plot_k+1
					if labelZonly ==True:
                                                label='Z='+str(star_z)
					for z in range(len(element_group)):
						if (z == (len(element_group)-1)):
							plt.plot(a_isotopes_sorted_group[z],prod_factor_sorted_group[z],marker=marker,color=color1,linestyle=linestyle,label=label)
						else:	
							plt.plot(a_isotopes_sorted_group[z],prod_factor_sorted_group[z],marker=marker,color=color1,linestyle=linestyle)
						print isotopes_sorted_group[z] 
						print a_isotopes_sorted_group[z]
					if withlabel==False:
						continue
					print '######## I Plot O LINES111',plot_o16_lines
					for z in range(len(prod_factor)):
						if isotopes[z] == 'O-16':
							#plot_o16_lines=False
							if plot_o16_lines==True:
								print '######## I Plot O LINES',prod_factor[z]
								plt.plot([0,208],[prod_factor[z],prod_factor[z]],"k-",linewidth=3)
								plt.plot([0,208],[prod_factor[z]*2,prod_factor[z]*2],"k--",linewidth=1)
								plt.plot([0,208],[prod_factor[z]/2.,prod_factor[z]/2.],"k--",linewidth=1)
								ax.text(208,prod_factor[z]*1.1, 'O-16',horizontalalignment='right',verticalalignment='bottom',clip_on=True,fontsize=9)
								ax.text(208, prod_factor[z]*2*1.1, 'O-16*2',horizontalalignment='right',verticalalignment='bottom',clip_on=True,fontsize=9)
								ax.text(208, (prod_factor[z]/2)*1.05, 'O-16/2',horizontalalignment='right',verticalalignment='bottom',clip_on=True,fontsize=9)

						#for last isotope, if its element appears (not in label_done) with it
						#if (z == (len(prod_factor)-1)) and (not (isotopes_sorted[z].split('-')[0]  in label_done)):
							#plt.plot(a_isotopes_sorted[z],prod_factor_sorted[z],label=label,marker=marker,color=color1,linestyle=linestyle)
						#if its not the first time the element appears but reaches end of array
						#elif (z == (len(prod_factor)-1)):
						#	idx_max=z
							#plt.plot(a_isotopes_sorted[idx_min:idx_max],prod_factor_sorted[idx_min:idx_max],label=label,marker=marker,color=color1,linestyle=linestyle)
						#	break

						if (isotopes_sorted[z].split('-')[0]  in label_done):			
							continue
						idx_max=z
						label_done.append(isotopes_sorted[z].split('-')[0])
						#plt.plot(a_isotopes_sorted[idx_min:idx_max],prod_factor_sorted[idx_min:idx_max],marker=marker,color=color1,linestyle=linestyle)
						idx_min=z
						fixed_labels=True
						if fixed_labels==True:
							if hh==0:
								nmin=0.1
							if hh==1:
								nmin=0.15
							if hh==2:
								nmin=0.2
							if hh==3:
								nmin=0.3
						else:
							nmin=prod_factor_sorted[z]*1.3
						#label below belongs to new element
						fontsizelabel='x-small'
						print a_isotopes_sorted[z],prod_factor_sorted[z]
                                                if (isotopes_sorted[z] in label_isotopes) and (len(label_isotopes)>0):
							ax.text(a_isotopes_sorted[z],nmin,isotopes_sorted[z],horizontalalignment='center',verticalalignment='center',fontsize=fontsizelabel,clip_on=True) #x-msall

						if hh==0:
							hh=1
							continue
						if hh==1:
							hh=2
							continue
						if hh==2:
							hh=3
							continue
						if hh==3:
							hh=0
							continue

                                                if (z == (len(prod_factor)-1)):
                                                        idx_max=z
                                                        break



					plt.legend()
                                        plt.yscale('log')
                                        plt.minorticks_on()
                                        if yields_output==False:
                                                plt.ylabel("Overproduction factor",fontsize=16)
                                        else:
                                                plt.ylabel("Yields/M$_{\odot}$",fontsize=16)

			if (GCE_tables) or ((plots_singleA) or (plots_singleZ)):
				return

			#after loop over all stars do plotting

                        color_iso=0
                        legend_k=0

			for h in range(len(isotopes)):
				yield_1=[]
				mass_1=[]			
				iniabu_folded=[]
				prod_fac_sorted=[]	
               			indices_array=sorted(range(len(star_mass_array)),key=lambda x:star_mass_array[x])
                		for i in indices_array:
					#if explosive contribution does not exist but requested
					if star_mass_array[i] in exclude_masses:
						continue
					if ((not star_mass_array[i] in plot_masses_only) and len(plot_masses_only)>0):
						continue
					mass_1.append(star_mass_array[i])
					prod_fac_sorted.append(yields_all_masses[h][i])
				####Plotting prodfactor
				fig_2=plt.figure(isotopes[h])
				ax_2=fig_2.gca()
				if plot_log==True:
					ax_2.set_yscale('log')
				ax_2.legend()
				ax_2.set_xlabel("M/M$_{\odot}$",fontsize=16)
				ax_2.minorticks_on()
				if yields_output==False:
					ax_2.set_ylabel("Overproduction factor",fontsize=16)
				else:
					ax_2.set_ylabel("Yields/M$_{\odot}$",fontsize=16)
				ax_2.set_title(title)
	
				if len(label)>0:
					label=", "+label
				if len(elements)==0:
					plot_quantity=isotopes[h]
					plot_quantity="$^{"+plot_quantity.split("-")[1]+"}$"+plot_quantity.split("-")[0]
				else:
					plot_quantity=elements[h]
				if withlabel==True:
					if withlabelelement==True:
                                		ax_2.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso],label=plot_quantity+" , Z="+str(star_z)+label  )
					else:
                                        	ax_2.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso],label="Z="+str(star_z)+label  )
				else:
                                	ax_2.plot(mass_1,prod_fac_sorted,marker=marker_type[legend_k],color=color[color_iso],markersize=markersize[legend_k],mfc=color[color_iso],linewidth=line_width[legend_k],linestyle=line_style[color_iso])				
				

				m_max=max(star_mass_array)+2
				
				if plot_lines==True:
					ax_2.plot([0,m_max],[1,1],"k--",linewidth=3)
					ax_2.plot([0,m_max],[2,2],"k--",linewidth=1)
					ax_2.plot([0,m_max],[0.5,0.5],"k--",linewidth=1)				
	
				color_iso+=1
				legend_k+=1
	

	def get_last_cycle(self,cycles,sefiles,sefiles_restart):

		'''
			Because there are a lot of issues with cycles we need this routine :(
		'''

		if cycles[1]==-1:
			#if star_mass>9:
			endcycle=int(sefiles_restart.se.cycles[-2])
			#else:
			#       endcycle=int(sefiles[k].se.cycles[-2])
		else:
			endcycle=cycles[1]
		cycs=range(cycles[0],endcycle,cycles[2])
		if endcycle not in cycs:
			cycs=cycs+[endcycle]
		print 'first, last cycle: ',cycles[0],endcycle
		star_mass11=sefiles.get(cycs,'mass')				
		#do a check for missing cycles
		if not len(star_mass11)==len(cycs):
			cycs1=[int(cycs[0])]
			for cyc in cycs[1:]:
				if  cyc <= max(cycs1):
					continue
				a=sefiles.get(cyc,'mass')
				cyctest=int(cyc)
				#if cycle 0 or a non-existing cycle, or cycle was alrady added, 
				#go one cycle further   
				while ((type(a)==list) or (a==0)):
					if   cyctest<int(cycs[-1]):
						print int(cycs[-1])
						print 'found bad cycle',cyctest
						cyctest+=1
						a=sefiles.get(cyctest,'mass')
					elif cyctest==int(cycs[-1]):
						cyctest=int(sefiles_restart.se.cycles[-3])
						a=sefiles.get(cyctest,'mass')	
						h=-4
						while ((type(a)==list) or (a==0)):
							cyctest=int(sefiles_restart.se.cycles[h])
							a=sefiles.get(cyctest,'mass')
							h-=1
						for kk in cycs1:
							if kk>=cyctest:
								cycs1=cycs1[:kk]
						break
					else:
						break
				print 'Use cycle ',cyctest
				if cyctest>max(cycs1):
					cycs1.append(cyctest)
			print len(cycs),len(cycs1)
			cycs=cycs1

		return cycs



	def write_prod_fact_stellarwinds_old(self,cycles=[],isotopes=["C-12","N-14","O-16","Ne-22"],elements=[],yields_output=True,headers=['this is a test','to write strings in the first col'],file_name="isotopic_table_set1.2_prodfac.txt",ascii_1=False,latex=False,GCE_input_fxt=True,GCE_input=False,GCE_input_starkenburg=False,exp_dir='',delay=True,exp_type=None,plot=False):

		'''
	
			Calculates the yields or production factors of isotopes or elements and writes them as
			- an ascii output
			- latex output
			- GCE input file , this is used by the write_GCE_input.... functions above

		
		'''
		import ascii_table as ascii1

		case_label_and_first_column = ['specie']
		yields_all=[]
		iso_names=[]
		remn_masses=[]

		#test to take right dir:
		# use exp_type variable from now on
		if exp_type == None:	#to be backwards compatible
			if delay:
				exp_type='delay'
			else:
				exp_type='rapid'
                slist = os.listdir(exp_dir)
                expr = re.compile(exp_type)
                slist=(filter(expr.search,slist))
		exp_runs_H5_out=[]
		exp_runs_H5_surf=[]
		
		for element in slist:
			run_path=exp_dir+'/'+element
			if not os.path.isdir(run_path):
				continue
			if os.path.isdir(run_path+"/H5_out") and os.path.isdir(run_path+"/H5_surf"):
                                sefiles = os.listdir(run_path+"/H5_out")
                                if (filter(expr.search,sefiles)) <1:
                                        print "Warning: No hdf5 out files found in "+run_path+"/H5_out"
                                else:
                                        exp_runs_H5_out.append(run_path+"/H5_out")
                                sefiles = os.listdir(run_path+"/H5_surf")
                                if (filter(expr.search,sefiles)) <1:
                                        print "Warning: No hdf5 surf files found in "+run_path+"/H5_surf"
                                else:
                                        exp_runs_H5_surf.append(se(run_path+"/H5_surf"))



		for i in range(len(self.runs_H5_surf)):
                	sefiles=se(self.runs_H5_surf[i])
			sefiles_hout=se(self.runs_H5_out[i])
			indx_iso=[]
			yields_1=[]
			iso_names=[]
			mass=sefiles.get("mini")
			metallicity=sefiles.get("zini")

			case_label_and_first_column.append(str(mass)+" Msun, Z="+str(metallicity))
			##if all stable ones:
			if len(isotopes)>0:
				if isotopes[0]=="all":
					isotopes_1=sefiles.se.isotopes
					isotopes=[]
					for j in range(len(isotopes_1)):
						if is_stable(isotopes_1[j])=='t':
							isotopes.append(isotopes_1[j])			
			#print isotopes		
			if len(elements)>0:
				if elements[0]=="all":
					elements=sefiles.se.elements
					elements=get_stable(sefiles.se.isotopes,get_elements=True)[:12]		
					#exclude neutrons
					#print elements
                        if cycles[i][1]==-1:
                        	endcycle=int(sefiles.se.cycles[-1]) #+ cycles[k][2]  #1000
                               	endcycle = int(endcycle) // 100 * 100
				endcycle = int(endcycle-20)
				#endcycle= int( round( endcycle,-3) -20 )
                                #print "endcycle",endcycle
                        else:
                                endcycle=cycles[i][1]
			#np.array(prod_facs),isotopes,yields,iniabus,remn_mass
			prod_factor,isotopes_prod_fac,yields,iniabu,mass_frac_ini,remn_mass =self.weighted_yields(sefiles,sefiles_hout,isotopes,elements,cycles[i][0],endcycle,cycles[i][2])
			#print '                                               '
			#print mass_frac_ini
			#for massive stars add pre-sn component
			if mass>8:

				#in case no ppn_exp run dir is available, skip explosion contribution
				#this is because mass cut is then larger then actual mass of star...no mass lost
				for t in range(len(exp_runs_H5_surf)):
					if type(exp_runs_H5_surf[t].get("mini")) == np.ndarray:
						mass1=exp_runs_H5_surf[t].get("mini")[0]
						metallicity1=exp_runs_H5_surf[t].get("zini")[0]
					else:
						mass1=exp_runs_H5_surf[t].get("mini")
						metallicity1=exp_runs_H5_surf[t].get("zini")

					###for currupted files

					corrupt=False
					if type(mass1) == list:
						corrupt=True
					if type(metallicity1) ==list:
						corrupt=True
					if ((float(mass1) < 1) or  (float(mass1) >70 )) or type(mass1) == type(None):
						corrupt=True
					if ((float(metallicity1) < 0.00001) or  (float(metallicity1) >1. )) or type(metallicity) == type(None):
						corrupt=True

					if corrupt == True:
						metallicity1=0.02
						rundir11=exp_runs_H5_out[t].split('/')[-2]
						mass1=float(rundir11.split('.')[0][1:])
						print 'corruption, assinge new values',mass1,metallicity1
					if mass1 == mass and metallicity1 == metallicity:
						sefiles=exp_runs_H5_surf[t]
                                        	sefiles_hout=se(exp_runs_H5_out[t])
						
						#calculate exp part
						import nugridse as mp
						reload(mp)
						prod_facs_exp,isotopes_exp,yields_exp,iniabus_exp,mass_cut_exp,first_cycle_exp=self.weighted_yields_explosion(mp,sefiles,sefiles_hout,isotopes=isotopes,elements=elements,delay=delay)
						#add pre-exp yields
						for tt in range(len(yields)):
							#print 'star above 8M'
							#print isotopes_prod_fac[tt],yields[tt],yields_exp[tt]
							yields[tt]=yields[tt]+yields_exp[tt]
						remn_mass=mass_cut_exp
						print 'exp masses added'
						break		

			#print 'isotopes yield',isotopes_prod_fac
			#print 'yields len',len(yields)
			#print 'prod factor:',prod_factor
			#print 'yield:',yields

			
			remn_masses.append(remn_mass)	
			if yields_output==False:	
				yields_all.append(list(prod_factor))
			else:
				yields_all.append(list(yields))
			if GCE_input==True:

				if len(elements)>0:
					dcols=['Elements','Yields']
				else:
					dcols=['Isotopes','Yields']
				special_header='Table: (M='+str(mass)+',Z='+str(metallicity)+')'
				data=[isotopes_prod_fac,yields_all[-1]]
                                ascii1.writeGCE_table(filename=file_name,headers=[special_header],data=data,dcols=dcols)
			if GCE_input_starkenburg==True:
				#needs mass_frac_ini,m_final
				#m_final=round(remn_mass,2)
				#print '                       ',mass_frac_ini
				#mass_frac_ini=mass_frac_ini
				dcols=['Isotopes','Yields','X0']
				finalmheader='Mfinal: '+'{:.3E}'.format(remn_mass)
				special_header='Table: (M='+str(mass)+',Z='+str(metallicity)+')'
				data=[isotopes_prod_fac,yields_all[-1],mass_frac_ini]
				ascii1.writeGCE_table(filename=file_name,headers=[special_header,finalmheader],data=data,dcols=dcols)	
			
		cols=yields_all
		if plot==True:
			return cols,isotopes_prod_fac
		#print cols
		#print case_label_and_first_column
		if latex==True:
			print 'latex is true'
			self.write_latex_table_set1(isotopes_prod_fac,cols,headers,file_name,case_label_and_first_column)

			#self.write_latex_table(header_file="Weighted stellar winds",cols=cols,label=label,caption=caption,table_header=table_header,latexfile=file_1,table_size=table_size)
		if ascii_1==True:
			self.write_ascii_table(header_file="weighted stellar winds",cols=cols,txtfile=file_1)
	
		#for GCE model of fxt	
		if GCE_input_fxt==True:
			return case_label_and_first_column,cols,remn_masses	








	


	def write_latex_table_set1(self,isotopes_for_table,data,headers=['this is a test','to write strings in the first col'],file_name_table_species_wind = 'isotopic_table_set1.2_prodfac.txt',case_label_and_first_column = ['specie','1.65 Msun','2 Msun','3 Msun','5 Msun','15 Msun','20 Msun','25 Msun',]):

		'''
			Table write method from Paper1 scripts (from Marco)

			isotopes_for_table - isotopes for table (first column), preferable stable
			data - 2d array e.g. [tab_w_1p65,tab_w_2,tab_w_3,tab_w_5,tab_w_15,tab_w_20,tab_w_25]


		'''


		# headers and format
		form_str='%7.3E'
		sep_string=' & '

		all_data=[]
		all_data.append(isotopes_for_table)
		for i in range(len(data)):
        		datal=list(form_str%data[i][j] for j in range(len(data[i])))
        		all_data.append(datal)


		### attempt to add the trailing '\\' to a latex table
		final_col=len(isotopes_for_table)*[r'\\ ']
		case_label_and_first_column.append(r'\\ ')
		all_data.append(final_col)
		att.write(file_name_table_species_wind,headers,case_label_and_first_column,all_data,sep=sep_string)
				



        def write_ascii_table(self,row_type="isotopes",header_file="Table element",cols=[[],[]],table_header=[],txtfile="test.txt",table_size="normal"):


                '''

                        cols - 2d array containing rows and cols as shown below
                        header - header of table in latex
                        label - label of the array, for reference in latex
                        caption - latex-typical description of plot
                        header_latexfile - header in latexfile 
                        table_size - if "small" , set latex table size to small, else else normall
                        
                        writes tables of the form:
			row_type	3Msun Z=...	4Msun Z=...
                        cols[0][0]      cols[1][0]      ... 
                        cols[0][1]      cols[1][1]      ...
                        ..              ...     
                        ..              ...             ...


                '''
                if txtfile[-4:] != ".txt":
                        txtfile=txtfile+".txt"


                file_1=open(txtfile,"w")
                #write header
                file_1.write("#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"+"\n")
                file_1.write("#   "+header_file+"\n")
                file_1.write("#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"+"\n")
                file_1.write("\n")
                file_1.write("\n")
                #write content
                if len(table_header)==0:
                        for i in range(len(self.runs_H5_surf)):
                                sefiles=se(self.runs_H5_surf[i])
                                mass=sefiles.get("mini")
                                z=sefiles.get("zini")
                                new_notation='{:.1E}'.format(float(z))
                                legend=str(mass)+"Msun Z= "+str(z)
                                table_header.append(legend)
                #write table header
                header_write=row_type
                for i in range(0,len(table_header)):
                        header_write+="  "+table_header[i]
                file_1.write(header_write+"\n")
                file_1.write("\n")
                #write columns
                for i in range(len(cols[0])):
                        line=cols[0][i]
                        for j in range(1,len(cols)):
                                a=cols[j][i]
				a="%.3E" % a
                                line+= "	"+ a		#str(round(float(str(a).split("e")[0]),3))+"E"+str(a).split("e")[1]
                        file_1.write(line+"\n")




	#########################################
        #####Internal or single run functions for MPPNP ####
	#########################################




	def set1_extension_weight(self):
		'''
			
			Calculates the IMF by using an adaption of the mass grid of paper1 (based on Marcos scripts).
			INcludes masses till 25 since set1.1 is the only one with 32 u. 60
			and plotting sets together give a more complete picture without them
		'''


		mass_ini_exp=[12.,15.,20.,25.]
		mass_ini_wind=[1.0,1.65,2.0,3.0,5.0,6.0,7.0,12.,15.,20.,25.]
		bounds = [0.65,1.45, 1.85, 2.5, 3.5,5.5,6.5, 8., 13.5, 17.5, 22.5, 28.5]
		salp_fac = -2.35
		# the crazyness below is just because the mass between 7. and 11.5 Msun 
		# is missing, super AGB stars. Marco
		# set 1.2
		facs=[bounds[0]**salp_fac-bounds[1]**salp_fac,bounds[1]**salp_fac-bounds[2]**salp_fac,bounds[2]**salp_fac-bounds[3]**salp_fac,bounds[3]**salp_fac-bounds[4]**salp_fac,bounds[4]**salp_fac-bounds[5]**salp_fac,bounds[5]**salp_fac-bounds[6]**salp_fac,bounds[6]**salp_fac-bounds[7]**salp_fac,bounds[7]**salp_fac-bounds[8]**salp_fac,bounds[8]**salp_fac-bounds[9]**salp_fac,bounds[9]**salp_fac-bounds[10]**salp_fac,bounds[10]**salp_fac-bounds[11]**salp_fac]
		# set 1.1
		#facs=[bounds[0]**salp_fac-bounds[1]**salp_fac,bounds[1]**salp_fac-bounds[2]**salp_fac,bounds[2]**salp_fac-bounds[3]**salp_fac,bounds[3]**salp_fac-bounds[4]**salp_fac,bounds[5]**salp_fac-bounds[6]**salp_fac,bounds[6]**salp_fac-bounds[7]**salp_fac,bounds[7]**salp_fac-bounds[8]**salp_fac]


		weighted_facs = [x/sum(facs) for x in facs]

		sum_mass_term = 0.
		for i in range(len(mass_ini_wind)):
			sum_mass_term = weighted_facs[i]*mass_ini_wind[i] + sum_mass_term
		sum_mass_term = 1./(sum_mass_term)

		# weighted_facs = weighted_facs/sum_mass_term

		# -2 is because we not consider m=32 and m=60.
		salp_fac_wind = weighted_facs[0:len(weighted_facs)]
		# add the -2 term if 32 Msun and 60 Msun are included for salp_fac_exp. Marco
		# add for set1.2
		#salp_fac_exp = weighted_facs[len(weighted_facs)-3-2:len(weighted_facs)-2]
		# add for set1.1
		salp_fac_exp = weighted_facs[len(weighted_facs)-4:len(weighted_facs)]
		 
		return salp_fac_wind,salp_fac_exp,sum_mass_term
	



	def set1_weight_marco(self,set1_2=True):

		'''
			Calculates the IMF as in paper1 (based on Marcos scripts)
			

		'''

		mass_ini_exp=[15.,20.,25.]
		salp_fac = -2.35
		if set1_2==True:
			mass_ini_wind=[1.65,2.0,3.0,5.0,15.,20.,25.,32.,60.]
			bounds = [1.45, 1.85, 2.5, 3.5, 7., 11.5, 17.5, 22.5, 28.5, 42, 80]
			facs=[bounds[0]**salp_fac-bounds[1]**salp_fac,bounds[1]**salp_fac-bounds[2]**salp_fac,bounds[2]**salp_fac-bounds[3]**salp_fac,bounds[3]**salp_fac-bounds[4]**salp_fac,bounds[5]**salp_fac-bounds[6]**salp_fac,bounds[6]**salp_fac-bounds[7]**salp_fac,bounds[7]**salp_fac-bounds[8]**salp_fac,bounds[8]**salp_fac-bounds[9]**salp_fac,bounds[9]**salp_fac-bounds[10]**salp_fac]
		else:
			mass_ini_wind=[1.65,2.0,3.0,5.0,15.,20.,25.]
			bounds = [1.45, 1.85, 2.5, 3.5, 7., 11.5, 17.5, 22.5, 30.]
			facs=[bounds[0]**salp_fac-bounds[1]**salp_fac,bounds[1]**salp_fac-bounds[2]**salp_fac,bounds[2]**salp_fac-bounds[3]**salp_fac,bounds[3]**salp_fac-bounds[4]**salp_fac,bounds[5]**salp_fac-bounds[6]**salp_fac,bounds[6]**salp_fac-bounds[7]**salp_fac,bounds[7]**salp_fac-bounds[8]**salp_fac]

		# the crazyness above is just because the mass between 7. and 11.5 Msun 
		# is missing, super AGB stars. Marco

		weighted_facs = [x/sum(facs) for x in facs]

		sum_mass_term = 0.
		for i in range(len(mass_ini_wind)):
			sum_mass_term = weighted_facs[i]*mass_ini_wind[i] + sum_mass_term
		sum_mass_term = 1./(sum_mass_term)

		# weighted_facs = weighted_facs/sum_mass_term

		# -2 is because we not consider m=32 and m=60.
		salp_fac_wind = weighted_facs[0:len(weighted_facs)]
		# add the -2 term if 32 Msun and 60 Msun are included for salp_fac_exp. Marco
		if set1_2==True:
			# add for set1.2
			salp_fac_exp = weighted_facs[len(weighted_facs)-3-2:len(weighted_facs)-2]
		else:
			 #add for set1.1
			salp_fac_exp = weighted_facs[len(weighted_facs)-3:len(weighted_facs)]
		 
		return salp_fac_wind,salp_fac_exp,sum_mass_term

	def weighted_yields_pre_explosion(self,mp,sefiles,sefiles_hout,endcycle,mass_cut,isotopes):
		
		'''
			Marcos routine adapted, calc pre SN yields
			explosion cycle is endcycle, at which pre SN yield is taken from
		'''


		import utils as u
		

		cycle=endcycle
		first_cycle=sefiles.se.cycles[0]	
		#mass_cut=1.60
                mass_last_cycle=sefiles_hout.get(cycle,"mass")
                idx_mcut=int(np.where(mass_last_cycle>mass_cut)[0][0])
                

		# I take the max mass coord
		mass_limit=mass_last_cycle[-1]
		# I calculate average_mass_frac_decay
		sefiles_hout.average_iso_abund_marco([mass_cut,mass_limit],cycle,True,2)

		E_i_pre_15=array(mp.average_mass_frac_decay)*(mass_limit-mass_cut)

		#follwoing needed for isotope identification,must be stable
		# define back_ind presn
		back_ind_pre=u.back_ind

		yields=[]
		iniabus=[]
		prod_facs=[]
		for i in range(len(isotopes)):
			other_name_scheme=isotopes[i].split("-")[0].upper()+(5-len(isotopes[i])+1)*" "+isotopes[i].split("-")[1] 
			if other_name_scheme not in u.stable:
				print "error: Chosen isotope is not stable:",other_name_scheme
				return 1
			idx_iso=back_ind_pre[other_name_scheme]
			yield_decayed=E_i_pre_15[idx_iso]
			yields.append(yield_decayed)
			iniabu_1=sefiles.get(0,isotopes[i])	
                        print "iniabu for exp",iniabu_1
                        iniabus.append(iniabu_1)
                        total_yield= yield_decayed
                        total_prod_factor=  total_yield/     ((mass_limit-mass_cut)*iniabu_1)
                        prod_facs.append(total_prod_factor)
                        print total_yield, total_prod_factor
                
		ini_star_mass=sefiles.get("mini")
		iniabus=np.array(iniabus)*(ini_star_mass)#-end_star_mass)
                
		return np.array(prod_facs),isotopes,yields,iniabus

	def prodfac_explosion(self,mp,sefiles,sefiles_hout,isotopes,ini_abu=""):
		'''
			Marcos routine adapted, calc SN yields
		'''

		import utils as u

                initial_abu=iniabu(ini_abu)                     #,iniabufile='../../frames/mppnp/USEEPP/"

		cycle=sefiles_hout.se.cycles[-1]
		first_cycle=sefiles_hout.se.cycles[0]	
		#mass_cut=1.60
                ###identify mcut at mass cell where rho >1 at first cycle,from core to surface
                rho_first_cycle=sefiles_hout.get(first_cycle,"rho")
                mass_first_cycle=sefiles_hout.get(first_cycle,"mass")
		print rho_first_cycle
		
		print cycle
		print sefiles_hout.get("mini")
		print sefiles_hout.get("zini")
                idx_mcut=int(np.where(rho_first_cycle!=1)[0][0])
                mass_cut=mass_first_cycle[idx_mcut]

		# I take the max mass coord
		mass_limit=mass_first_cycle[-1]
		# I calculate average_mass_frac_decay
		sefiles_hout.average_iso_abund_marco([mass_cut,mass_limit],cycle,True,2)

		E_i_pre_15=array(mp.average_mass_frac_decay)*(mass_limit-mass_cut)

		#follwoing needed for isotope identification,must be stable
		# define back_ind presn
		back_ind_pre=u.back_ind

		yields=[]
		iniabus=[]
		prod_facs=[]
		for i in range(len(isotopes)):
			other_name_scheme=isotopes[i].split("-")[0].upper()+(5-len(isotopes[i])+1)*" "+isotopes[i].split("-")[1] 
			if other_name_scheme not in u.stable:
				print "error: Chosen isotope is not stable: ",other_name_scheme
				return 1
			idx_iso=back_ind_pre[other_name_scheme]
			yield_decayed=E_i_pre_15[idx_iso]
			yields.append(yield_decayed)	
			iniabu_1=initial_abu.abu[ initial_abu.names.index(isotopes[i].split("-")[0].lower()+(5-len(isotopes[i])+1)*" "+isotopes[i].split("-")[1])]	
                        print "iniabu for exp",iniabu_1
                        iniabus.append(iniabu_1)
                        total_yield= yield_decayed
                        total_prod_factor=  total_yield/     ((mass_limit-mass_cut)*iniabu_1)
                        prod_facs.append(total_prod_factor)
                        print total_yield, total_prod_factor
                
		ini_star_mass=sefiles.get("mini")
		iniabus=np.array(iniabus)*(ini_star_mass)#-end_star_mass)
                
		return np.array(prod_facs),isotopes,yields,iniabus,mass_cut,first_cycle


	def weighted_yields_explosion(self,mp,sefiles_re,sefiles,isotopes=[],elements=[],cycleend=1,delay=True):
		'''
			Marcos routine adapted, calc SN yields
			meeds tje correct path ini_abu
		'''

		import utils as u
		import nugridse as mp

		isotopes1=sefiles.se.isotopes	
		print 'Entering weighted_yields_explosion'
		#get all isotopes related to element input
                if len(elements)>0:
			ini_all=sefiles.get(1,'iso_massf') 
                        isotopes=[]
			ini_isotopes=[]
                        for i in range(len(isotopes1)):
				elem1=isotopes1[i].split('-')[0]
				if elem1 in elements:
					#stable condition mustbe due to u.stable below
					if is_stable(isotopes1[i])=='t':
                        			isotopes.append(isotopes1[i])
						ini_isotopes.append(ini_all[i])

		#get isotopes foriso inputinrightorder
		if len(elements)<1:
			print 'try iniall'
                        ini_all=sefiles.get(1,'iso_massf')
                        #isotopes=[]
			ini_isotopes=[0]*len(isotopes)
                        for i in range(len(isotopes1)):
                                if isotopes1[i] in isotopes:
					idx=isotopes.index(isotopes1[i])
                                        ini_isotopes[idx]=ini_all[i]
					

		cycle=cycleend#sefiles_re.se.cycles[-1]
		print 'last cycle ',cycle
		print 'last cycle spec',sefiles_re.se.cycles[-1]
		#mass_cut=1.60i
                ###identify mcut at mass cell where rho >1 at first cycle,from core to surface
                mass_cycle=sefiles_re.get(cycle,"mass")
		#print rho_first_cycle
		
		#print cycle

                if type(sefiles.get("mini")) == np.ndarray:
                        m_ini=sefiles.get("mini")[0]
			z_ini=sefiles.get("zini")[0]
                else:
                        m_ini=sefiles.get("mini")
			z_ini=sefiles.get("zini")
		#IMPORTANT CHANGE FOR M12 - SEE NUGRID MEETING 2014
		m_ini_save=m_ini
		if int(m_ini)==12:
			m_ini=15
		z_metal=z_ini/0.02

		print 'MINI',m_ini,z_metal
		if ((m_ini>=11.) and (m_ini<30.)):				
			if delay==True:
				mass_cut = 1.1 + 0.2*np.exp((m_ini-11.0)/4.) - (2.0 + z_metal)*np.exp(0.4*(m_ini -26.0))
			####rapid cc
			else:
				if m_ini<22.:
					mass_cut= 1.1 +0.2*np.exp((m_ini-11.0)/7.5) + 10*(1.0+z_metal)*np.exp(-(m_ini-23.5)**2/(1.0+z_metal)**2)
				elif m_ini<30 :
					mass_cut=  1.1 + 0.2*np.exp((m_ini-11.0)/4.) - (2.0 + z_metal)*np.exp(0.4*(m_ini -26.0)) - 1.85 + 0.25*z_metal +10.0*(1.0+z_metal)*np.exp(-(m_ini-23.5)**2/(1.0+z_metal)**2)
		##at higher mass difference
		elif ((m_ini>30) and (m_ini<50)):
			#delay
			if delay==True:
				mass_cut= min( 33.35 + (4.75 + 1.25*z_metal)*(m_ini-34.),m_ini-z_metal**0.5 *(1.3*m_ini - 18.35))		
			else:
				mass_cut = min( 33.35 + (4.75 + 1.25*z_metal)*(m_ini-34.),m_ini-z_metal**0.5 *(1.3*m_ini - 18.35)) - 1.85 + z_metal*(75. -m_ini)/20.
		elif m_ini>50:
			#Page 7, Fryer12, only at solar Z valid
			if z_metal==1:
				if m_ini<90.:
					mass_cut = 1.8 + 0.04*(90. - m_ini)		
				else:
					mass_cut = 1.8 + np.log10(m_ini - 89.)
			#The part below will probably never be used
			if z_metal <1:
				if m_ini<90.:
					mass_cut = max(min( 33.35 + (4.75 + 1.25*z_metal)*(m_ini-34.),m_ini-z_metal**0.5 *(1.3*m_ini - 18.35)),1.8 + 0.04*(90. - m_ini))	
				else:
					mass_cut = max(min( 33.35 + (4.75 + 1.25*z_metal)*(m_ini-34.),m_ini-z_metal**0.5 *(1.3*m_ini - 18.35)),1.8 + np.log10(m_ini - 89.))

		mass_cut=round(mass_cut,2)
		#mass_cut=1.6
	
		#Adjust mass cut if chosen as input, variable delay becomes irrelevant
                print 'self.fallback_coords'
                print self.fallback_coords
                if len(self.fallback_coords)>0:
			if m_ini_save in self.fallback_masses:
                        	idx=self.fallback_masses.index(m_ini_save)
                        	mass_cut=self.fallback_coords[idx]
			else:
				print 'input fallback prescription for ',m_ini_save,' not available'
		# I take the max mass coord
		print 'Use mass cut ',mass_cut
		mass_limit=mass_cycle[-1]
		#case only if pre-exp and remnant larger than total mass
		if mass_limit<= mass_cut:
			#species,yields,iniabus,mass_cut,cycle
			return [],[],[],0,0
			#species,yields,iniabus,mass_cut,cycle 
		# I calculate average_mass_frac_decay
		reload(u);reload(mp)
		print 'choice of mass cut:',m_ini,z_ini,mass_cut,mass_limit,cycle
		sefiles_re.average_iso_abund_marco([mass_cut,mass_limit],cycle,True,2)
		#print 'Calculated mass_limit-mass_cut',mass_limit-mass_cut
		E_i_pre_15=array(mp.average_mass_frac_decay)*(mass_limit-mass_cut)
		#follwoing needed for isotope identification,must be stable
		# define back_ind presn
		back_ind_pre=u.back_ind
  	

		idx = back_ind_pre['K  39']
		#print 'test K39',E_i_pre_15[idx],idx



		yields=[]
		iniabus=[]
		#prod_facs=[]
		for i in range(len(isotopes)):
			other_name_scheme=isotopes[i].split("-")[0].upper()+(5-len(isotopes[i])+1)*" "+isotopes[i].split("-")[1] 
			if other_name_scheme not in u.stable:
				print "Chosen isotope is not stable: ",other_name_scheme
				continue
			idx_iso=back_ind_pre[other_name_scheme]
			yield_decayed=E_i_pre_15[idx_iso]
			#print isotopes[i],other_name_scheme,yield_decayed,idx_iso
			yields.append(yield_decayed)
			#iniabu_1=ini_isotopes[i] *m_ini
			iniabu_1=ini_isotopes[i] * (mass_limit-mass_cut)
			if isotopes[i] == 'H-1':
				print ini_isotopes[i]
				print mass_limit
				print mass_cut
			#iniabu_1=initial_abu.abu[ initial_abu.names.index(isotopes[i].split("-")[0].lower()+(5-len(isotopes[i])+1)*" "+isotopes[i].split("-")[1])]	
                        #print "iniabu for exp",iniabu_1
                        iniabus.append(iniabu_1)
                        total_yield= yield_decayed
                        #total_prod_factor=  total_yield /     ((mass_limit-mass_cut)*iniabu_1)
                        #prod_facs.append(total_prod_factor)
                        #print total_yield, total_prod_factor


		#print 'elements:',elements
		if len(elements)>0:
			elem_data=[0]*len(isotopes)
			elem_name=[]
			iniabus1=[0]*len(isotopes)
			for ww in range(len(isotopes)):
                                elem1=isotopes[ww].split('-')[0]
				idx=elements.index(elem1)
                                if not elem1 in elem_name:
                                	elem_name.append(elem1)
                                        elem_data[idx]=yields[ww]
					iniabus1[idx]=iniabus[ww]
				else: 
					elem_data[idx] = elem_data[idx] + yields[ww]	
					iniabus1[idx]=iniabus1[idx] + iniabus[ww]
			iniabus=iniabus1
			species=elements
			yields=elem_data
			#prod_facs=np.array(elem_data)/np.array(iniabus) 
		else:
			species=isotopes

		print 'weighted_yields_explosion successfull' 
		return species,yields,iniabus,mass_cut,cycle


	def weighted_yields(self,sefiles,sefiles_restart,isotopes,elements=[],cycs=[]):
		'''
	  	    This function returns the production factors, isotope names,yields and initial
		    abundances of isotopes or elements
		    If elements>0 than ignore isotope input.

			!!!Isotopes and elements input must be available in hdf5 files!!!
	 
		'''
	

                import utils as u
                import re
                import nugridse as mp
		import gc
		import os

                if type(sefiles.get("mini")) == np.ndarray:
			ini_star_mass=sefiles.get("mini")[0]
		else:
			ini_star_mass=sefiles.get("mini")

		print "Initial mass:",ini_star_mass
		
	        #cycs=range(cyclestart,cycleend+sparse,sparse)#cycleend,sparse)
		#print cycs[0],cycs[-1],cycleend
		print 
		print 'h5 surf'
		star_mass=sefiles.get(cycs,'mass')
		cycleend=cycs[-1]
		print 'final cycle',cycleend
		#test for missing cycle
		print 'h5 restart'
                end_star_mass=sefiles_restart.get(cycleend,"mass")[-1]
                h_mass_frac=sefiles_restart.get(cycleend,"H-1")
                mass_array=sefiles_restart.get(cycleend,"mass")
                #earlier: 1.e-1 but now as in set1 scripts 0.05
		h_free_core=mass_array[np.where(h_mass_frac<5.e-2)[0][-1]]
		envelope_mass=end_star_mass-h_free_core

                #print "Final star mass = "+str(end_star_mass)
                #print "Remnant mass = "+str(h_free_core)

		#print 'star data read'

		delta_star_mass=np.array(star_mass[1:])-np.array(star_mass[0:-1])

		#if len(elements)>0:
		#	isotopes=elements
		iso_data1=sefiles.get(cycs,'iso_massf')
		isotopes1=sefiles.se.isotopes

		# in the case of isotope input
		#choose the right isotopes

		if len(elements)<1:
			iso_data=[0]*len(isotopes)		
			for w in range(len(iso_data1[0])):
				if isotopes1[w] in isotopes:
					idx=isotopes.index(isotopes1[w])
					iso_data[idx]=[]
					for k in range(len(cycs)):
						iso_data[idx].append(iso_data1[k][w])
			species_data=iso_data
			species=isotopes
		#in the case of elements
		if len(elements)>0:

			elem_data=[]
			elem_name=[]
	
			#the part doewn commented because all wind isotopes are already decayed
			#for isotope in isotopes1: 
			#	if is_stable(isotope)=='t':
			#		isotopes.append(isotope)
			isotopes=[]
			iso_data=[]	
			for w in range(len(iso_data1[0])):
				isotopes.append(isotopes1[w])
				iso_data.append([])
				for k in range(len(cycs)):
					iso_data[-1].append(iso_data1[k][w])

			elem_data=[0]*len(elements)
			for ii in range(len(isotopes)):
				elem1=isotopes[ii].split('-')[0]
				#Assume that 
				if elem1 in elements:
					idx=elements.index(elem1)
					if not elem1 in elem_name:
						elem_name.append(elem1)		
						elem_data[idx]=iso_data[ii]
					else:
						elem_data[idx]= np.array(list(elem_data[idx])) + np.array(list(iso_data[ii]))
			#elements to have elements input order
			species=elements
			species_data=elem_data

		prod_facs=[]
		yields=[]
		iniabus=[]
		envelope_yield=[]
        
		#specie out
		mass_frac_ini=[]
		for i in range(len(species)):
			specie=species_data[i]
			specie_ave=0.5*(np.array(specie[0:-1])+np.array(specie[1:]))
        		yield_specie_wind=np.sum(specie_ave*-1.*delta_star_mass)
			iniabu_1=specie[0]
			mass_frac_ini.append(iniabu_1)  
			envelope_yield=envelope_mass*specie[-1]
			###assume M>=8 stars do collapse
			if ( ini_star_mass < 8. ):
				total_yield= yield_specie_wind + envelope_yield
				#print  species[i],total_yield,ini_star_mass,iniabu_1
				total_prod_factor=  total_yield/ (ini_star_mass*iniabu_1)
				#iniabus.append(ini_star_mass*iniabu_1 )
				iniabus.append(	(ini_star_mass-h_free_core)*iniabu_1 )
				#print species
				#total_prod_factor=  total_yield/((ini_star_mass-end_star_mass)*iniabu_1) 
				#iniabus.append( (ini_star_mass-end_star_mass)*iniabu_1 )

			else:
				total_yield= yield_specie_wind
				#print  total_yield,ini_star_mass,iniabu_1
				total_prod_factor=  total_yield/(ini_star_mass*iniabu_1)     
				#total_prod_factor=  total_yield/((ini_star_mass-end_star_mass)*iniabu_1)
				#iniabus.append( ini_star_mass*iniabu_1)
				iniabus.append( (ini_star_mass-end_star_mass)*iniabu_1 )
			yields.append(total_yield)
			prod_facs.append(total_prod_factor)
			#print total_yield, yield_c12_wind,envelope_yield,(ini_star_mass-end_star_mass)*iniabu_1
			#print "massloss",ini_star_mass,end_star_mass,iniabu_1
	
		#from masss frac iniabu to mass unit
		#iniabus=np.array(iniabus)*(ini_star_mass)#-end_star_mass)
		remn_mass=h_free_core
		#print mass_frac_ini
		#return memory/garbage
		gc.collect()
		print 'weihted yield sucuceesfull'
		#output: elements,isotopes in elements/isotopes input order

		return np.array(prod_facs),species,yields,iniabus,mass_frac_ini,remn_mass
	


	def plot_abu_atmasscoord(self,mp,sefiles,cycle,masscell,isotopes_input,x_charge,mass_number_range,decayed,filename_solar_norm,elem_norm,yax_log,label,legend,color,title):
		import utils as u

		'''
			if mass_number_range greater 0, ignore isotopes input
	
		'''
		
		####decide which isotopes to process
		isotopes_stable=[]	
		isotopes=[]
		#if mass number range is needed
		if len(mass_number_range)>0:
			#get all isotopes in mass number range
			if mass_number_range[0] == 0:
				mass_number_min=sefiles.se.isotopes[0]
			else:
				mass_number_min=mass_number_range[0]
			if mass_number_range[1] == 0:
				mass_number_max=sefiles.se.isotopes[-1]
			else:
				mass_number_max=mass_number_range[1]	
			for i in range(len(sefiles.se.isotopes)):
				iso=sefiles.se.isotopes[i]
				#print iso.split("-")[1]
				#print mass_number_min 
				#print mass_number_max
				iso_number= iso.split("-")[1]
				iso_number= re.findall ('[\d ]+',iso_number)[0]
				if ( (int(iso_number)>=mass_number_min) and (int(iso_number)<=mass_number_max) ):
					if (is_stable(iso) == 't'):
						isotopes_stable.append(iso)
					isotopes.append(iso)
		else:
			isotopes=isotopes_input
			for i in range(len(sefiles.se.isotopes)):
				iso=sefiles.se.isotopes[i]
                                #print iso.split("-")[1]
                                #print mass_number_min 
                                #print mass_number_max
                                iso_number= iso.split("-")[1]
                                iso_number= re.findall ('[\d ]+',iso_number)[0]
                                if (is_stable(iso) == 't'):
                                	isotopes_stable.append(iso)
		
		print isotopes

		#############get isotopic distribution
		u.convert_specie_naming_from_h5_to_ppn(sefiles.se.isotopes)
		masses_for_this_cycle = sefiles.se.get(cycle,'mass') #<- find mass index you want
		idx=0	
		#print sefiles.isotopic_production_factors
		idx=np.abs(np.array(masses_for_this_cycle)-masscell).argmin()
		print "closest value: ", masses_for_this_cycle[idx]	
		plotiso=[]
		mass=[]
		plotiso_norm=[]
		iso_abu=[]
		iniabu_array=[i]
		name=[]
		
		mass=sefiles.se.A
		charge=sefiles.se.Z
		#sefiles.stable_isotope_list[isotope] <> get isotopes you want
		sefiles.iso_index=[]
		iniabu_norm=0	
		for j in range(len(sefiles.se.isotopes)):
			if sefiles.se.isotopes[j] in isotopes:
				name.append(sefiles.se.isotopes[j])
				if decayed == False:
					iso_abu.append(sefiles.get(cycle,sefiles.se.isotopes[j])[idx])
				sefiles.iso_index.append(j)
			if len(elem_norm)==1:
				if elem_norm in sefiles.se.isotopes[j]:
					iniabu_norm+=sefiles.get(cycle,sefiles.se.isotopes[j])[0]
					print "iniabu_norm....",iniabu_norm

		if decayed ==True:
               		########isotope decay, new abundances#######
			#isomers are excluded for now, in convert_specie naming...last ones 5
			#to decay whole network use all isotopes, change iso_abu
			iso_abu=sefiles.get(cycle,'iso_massf')[idx]
                	isotope_names = sefiles.se.isotopes
                	sefiles.mass_frac = [iso_abu] #double array important for decay

                	u.convert_specie_naming_from_h5_to_ppn(isotope_names) #-define u.spe
                	names_ppn_world_sefiles = u.spe
                	number_names_ppn_world = u.n_array
                	u.define_zip_index_for_species(names_ppn_world_sefiles,number_names_ppn_world) #--define u.cl!

                	#let species decay
                	u.stable_specie()
			print name
			print u.stable
			print "lentest:",len(u.stable),len(mp.decayed_multi_d[0])
			#print iso_abu
                	sefiles.decay([iso_abu]) #double array is important for decay()
                	iso_abu_decayed=mp.decayed_multi_d[0]
			print mp.decayed_multi_d


		#convert specified isotopes in correct name scheme
		u.convert_specie_naming_from_h5_to_ppn(isotopes) #-define u.spe
		names_ppn_world_chosen_isotopes = u.spe
	
		for i in range(len(names_ppn_world_chosen_isotopes)):
			names_ppn_world_chosen_isotopes[i]=names_ppn_world_chosen_isotopes[i].upper()
		
		#for normalization
		#u.solar(filename_solar_norm,solar_factor=1)


		mass_specific_isotopes=[]
		charge_specific_isotopes=[]
		label_specific_isotopes=[]
		###########Calculate production factor
		if decayed == True:
			for j in range(len(u.stable)):
				if u.stable[j] in names_ppn_world_chosen_isotopes:
					#if sefiles isotope in both name schemes have same order, than following is correct
					iso_idx=names_ppn_world_sefiles.index(u.stable[j].capitalize())
					#print u.stable[j].capitalize(), sefiles.se.isotopes[iso_idx]
					if len(elem_norm)==0:
						iniabu_norm=sefiles.get(0,sefiles.se.isotopes[iso_idx])[0]
						print "this is not "
                                        print iniabu_norm
					if iniabu_norm <1e-20:
                                                	print "Warning: Isotope "+sefiles.se.isotopes[i]+" not in initial abundance file, no plotting"
                                        else: 
						#norm=iniabu_array[iso_idx]		#u.solar_abundance[isotopes_ppn_world]
						plotiso.append(iso_abu_decayed[j]/iniabu_norm)
						mass_specific_isotopes.append(mass[iso_idx])
						label_specific_isotopes.append(sefiles.se.isotopes[iso_idx])
		else:			
			j=0
			for i in range(len(sefiles.se.isotopes)):
				if sefiles.se.isotopes[i] in isotopes:
					print sefiles.se.isotopes[i],isotopes
					if len(elem_norm)==0:
						print "count"+str(i)
						iniabu_norm=1
						print "print info:",sefiles.get(0,sefiles.se.isotopes[i])[0]
						print "aaaa",iniabu_norm
					print "aaaaaaaaaa",iniabu_norm
					#if len(iniabu_norm)>1:
					#	iniabu_norm=iniabu_norm[0]
					if iniabu_norm <1e-20:
						print "Warning: Isotope "+sefiles.se.isotopes[i]+" not in initial abundance file, no plotting"	
					else:	
						plotiso.append(iso_abu[j]/iniabu_norm)
						mass_specific_isotopes.append(mass[i])
						charge_specific_isotopes.append(charge[i])
						label_specific_isotopes.append(sefiles.se.isotopes[i])
					j+=1
		
		####Plotting


		fig=figure(0)
		if yax_log==True:
                	ax = fig.add_subplot(1,1,1)
                        ax.set_yscale('log')

		if x_charge==False:
			plt.plot(mass_specific_isotopes,plotiso,'*',label=legend,marker='*',markersize=8,mfc=color,linestyle='None')

	                if label ==True:
        	                for j in range(len(plotiso)):
                	                plt.annotate(label_specific_isotopes[j], xytext = (0, 10),textcoords = 'offset points' ,xy=(mass_specific_isotopes[j], plotiso[j]))


		else:
			plt.plot(charge_specific_isotopes,plotiso,'*',label=legend,marker='*',markersize=8,mfc=color,linestyle='None')
			if label ==True:
				for j in range(len(plotiso)):
					plt.annotate(label_specific_isotopes[j], xytext = (0, 10),textcoords = 'offset points' ,xy=(charge_specific_isotopes[j], plotiso[j]))

                plt.rcParams.update({'font.size': 16})
                plt.rc('xtick', labelsize=16)
                plt.rc('ytick', labelsize=16)
		plt.minorticks_on()
		plt.xlim(mass_specific_isotopes[0],mass_specific_isotopes[-1])
		plt.legend()
		plt.title(title)
 		plt.ylabel("$X_{i}/X_{ini}$",fontsize=20)
		if x_charge==False:
			plt.xlabel("Mass Number (A)",fontsize=20)
		else:
			plt.xlabel("Charge (Z)",fontsize=20)


	def plot_final_surf_abu(self,yaxis='[Ba/Eu]',color='r',marker='o',linestyle='-',iniabupath='/astro/critter/critter/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn'):


		elem1=yaxis.split('/')[0].split('[')[1]
		elem2=yaxis.split('/')[1].split(']')[0]

		elem1_ini=get_ini_elem_abu(element=elem1,iniabupath=iniabupath)
		elem2_ini=get_ini_elem_abu(element=elem2,iniabupath=iniabupath)

		spec_values=[]
		masses=[]
		for i in range(len(self.runs_H5_surf)):
			
			sefiles=se(self.runs_H5_surf[i])
			abu=[]
			abu1=self.get_elem_abu(cycles=-1,element=elem1,sefiles=sefiles)[0]	
			abu2=self.get_elem_abu(cycles=-1,element=elem2,sefiles=sefiles)[0]
			spec_values.append(np.log10(abu1/elem1_ini * elem2_ini/abu2))
                        print abu1,abu2
                        print elem1_ini,elem2_ini

			mass=sefiles.get("mini")
                        z=sefiles.get("zini")
			masses.append(mass)
                        legend=str(mass)+"M$_{\odot}$ Z= "+str(z)
		print spec_values
		print 'for ',elem1,elem2
		plt.figure(elem1+'_'+elem2)#+', Z='+str(z))
		plt.plot(masses,spec_values,marker=marker,color=color,linestyle=linestyle,label='Z='+str(z))
		plt.ylabel(yaxis)
		plt.xlabel('M/M$_{\odot}$')
		plt.legend(loc=1)


	def plot_surf_element(self,element='Ba'):
		
		'''
			Plots surface evolution of certain element

		'''

		for i in range(len(self.runs_H5_surf)):
			sefiles=se(self.runs_H5_surf[i])
			y=[]
			mass=sefiles.get("mini")
                        z=sefiles.get("zini")
                        legend=str(mass)+"M$_{\odot}$ Z= "+str(z)
			for j in arange(0,len(sefiles.se.cycles)-1000,500):
    				y.append(sefiles.get(j,"elem_massf",element))
			plt.plot(arange(0,len(sefiles.se.cycles)-1000,500),y,label=legend)
			plt.xlabel('Cycles')
			plt.ylabel('X('+element+')')

        def plot_surf_isotope(self,isotope='C-12'):
               
                '''
                        Plots surface evolution of certain isotope

                '''

                for i in range(len(self.runs_H5_surf)):
                        sefiles=se(self.runs_H5_surf[i])
                        y=[]
                        mass=sefiles.get("mini")
                        z=sefiles.get("zini")
                        legend=str(mass)+"M$_{\odot}$ Z= "+str(z)
                        for j in arange(0,len(sefiles.se.cycles)-1000,500):
                                y.append(sefiles.get(j,"iso_massf",isotope))
                        plt.plot(arange(0,len(sefiles.se.cycles)-1000,500),y,label=legend)
                        plt.xlabel('Cycles')
                        plt.ylabel('X$_{'+isotope+'}$')


	def find_most_abundant(self,coords=[],cycle=[]):
		
		minis=[]
		isotopes=[]
		abus=[]
                for i in range(len(self.runs_H5_out)):
                        sefiles=se(self.runs_H5_out[i])
                        mass_ini=sefiles.get("mini")
                        z=sefiles.get("zini")
			if cycle[i]==-1:
				cyc=sefiles.se.cycles[-1]
			else:
				cyc=cycle[i]
			coord=coords[i]
			mass=sefiles.get(cyc,'mass')
			mass_idx=min(range(len(mass)), key=lambda i: abs(mass[i]-coord))
			abu=sefiles.get(cyc,'iso_massf')[mass_idx]
			idx=list(abu).index(max(abu))
			print 'at mass coord',coord,': ',sefiles.se.isotopes[idx],'with ',max(abu)
			minis.append(mass_ini)
			isotopes.append(sefiles.se.isotopes[idx])
			abus.append(max(abu))

		for k in range(len(minis)):
			print 'M=',minis[k],'Z=',z,'at mass coord',coords[k],': ',isotopes[k],'with ',abus[k]

	def set_plot_prof_ratio(self,cycles=[],isotopes=['C-12','O-16'],savepath='/nfs/rpod3/critter/Results/DATA/DataTest/M15.0Z2.0e-02_SamTest/MESA_MPPNP_comparison',color='r',line_style='--',marker='x',markevery=200,label='',ylim=[0,0],write_label=''):
		''' 
			Plots abundance profiles of cycle cycle
			If cycle=-1 take last cycle.
		'''
		fig_num=20
		for i in range(len(self.runs_H5_out)):
			plt.figure(fig_num)
			fig_num = fig_num + 1
			sefiles=se(self.runs_H5_out[i])
			mass_ini=sefiles.get("mini")
			z=sefiles.get("zini")
			cyc=cycles[i]
			print cyc,isotopes[0]
			print cyc,isotopes[1]
			legend="M="+str(mass_ini)+", Z= "+str(z)
			mass=sefiles.get(cyc,"mass")
			iso1=sefiles.get(cyc,isotopes[0])
			iso2=sefiles.get(cyc,isotopes[1])
			print 'iso1',iso1
			print 'iso2',iso2
			plt.plot(mass,iso1/iso2,color=color,linestyle=line_style,label=isotopes[0]+'/'+isotopes[1]+label,marker=marker,markevery=markevery)
			plt.xlabel("M/M$_{\odot}$",fontsize=20)
			plt.ylabel('Ratio mass fraction', fontsize=20)
                        plt.legend()
			if ylim[1]>0:
				plt.ylim(ylim[0],ylim[1])
				
			if len(savepath)>0:
				plt.savefig(savepath+'/'+legend+write_label+'.png')
				#plt.close('all')
        
	def set_plot_prof(self,fig=-1,cycle=0,isotopes=['He-4','C-12','O-16','Ne-20','Mg-24','Si-28'],savepath='/astro/critter/Documents/PhD/Physics/na_low_zmass_paper/c13pocket_profiles',color=['r','k','b','g','y','m','c','0.35'],line_style=[],mass_range=[],marker=['s','d','o','p','^','>','*'],markevery=200,withdensity=False,firstplot=True,logplot=False,ylimits=[]):
		''' 
			Plots abundance profiles of cycle cycle
			of isotopes isotopes.
			If cycle=-1 take last cycle.
		'''
		for i in range(len(self.runs_H5_out)):
			sefiles=se(self.runs_H5_out[i])
			mass_ini=sefiles.get("mini")
			z=sefiles.get("zini")
			if type(cycle)==int:
                                if cycle==-1:
                                        cyc=int(sefiles.se.cycles[-1])-20
				else:
					cyc=cycle
				legend="M="+str(mass_ini)+", Z= "+str(z)+", cycle="+str(cyc)
				if fig==-1:
                                	fig1=plt.figure(legend)
                        	else:
                                	fig1=plt.figure(fig)
				mass=sefiles.get(cyc,"mass")
				rho=sefiles.get(cyc,'rho')
				print  rho
				for j in range(len(isotopes)):
					iso=sefiles.get(cyc,isotopes[j])
					print iso,mass
					if len(line_style)==0:
						line_style=50*['-']
					
					if not withdensity==True:
						if logplot==True:
							y=np.log10(iso)
						else:
							y=iso
						if not fig==-1:		
							plt.plot(mass,y,marker=marker[i],markevery=markevery,linestyle=line_style[i],color=color[i],label=isotopes[j]+', '+legend)
						else:
							plt.plot(mass,y,marker=marker[j],markevery=markevery,linestyle=line_style[j],color=color[j],label=isotopes[j])
				plt.legend()
				if withdensity==True:
					if firstplot==True:
						if j==0:
							ax=plt.gca()
							ax2=ax.twinx()
						else:
							ax=fig1.axes[0]
							ax2=fig1.axes[1]
					else:
						ax=fig1.axes[0]
						ax2=fig1.axes[1]
					if logplot==True:
						y=np.log10(iso)
					else:
						y=iso
					ax.plot(mass,y,marker=marker[i],markevery=markevery,linestyle=line_style[i],color=color[j],label=isotopes[j]+', M='+str(mass_ini)+'Z='+str(z))
					ax2.plot(mass,rho,marker=marker[i],markevery=markevery,linestyle=line_style[i],color=color[j],label='rho, M='+str(mass_ini)+'Z='+str(z))
					ax2.set_yscale('log')
					ax2.set_ylim(1,1e9)
					ax2.legend(loc=4)
				else:
					ax=plt.gca()
				if not logplot==True:
					ax.set_yscale('log')	
					ax.set_ylim(-4,1)
				else:
					ax.set_ylim(-4,0)
				plt.xlabel("M/M$_{\odot}$",fontsize=20)
				plt.ylabel('Log($X_{i})$', fontsize=20)	
				ax.legend()
				if len(mass_range)>0:
					plt.xlim(mass_range[0],mass_range[1])
				if len(savepath)>0:
					if not fig==-1:
						plt.savefig(savepath+'/'+str(cyc)+'.png')
					else:   
						plt.savefig(savepath+'/'+legend+str(cyc)+'.png')
						
			if type(cycle)==list:
				print 'dealwith list'
				for cyc in cycle:
					legend="M="+str(mass_ini)+", Z= "+str(z)+", cycle="+str(cyc)
		                        if fig==-1:
        		                        fig1=plt.figure(legend)
        		                else:
                	                	fig1=plt.figure(fig)

					if cyc==-1:
                                        	cyc=int(sefiles.se.cycles[-1])
                                	mass=sefiles.get(cyc,"mass")
                                	for j in range(len(isotopes)):
                                        	iso=sefiles.get(cyc,isotopes[j])
                                        	print iso,mass
					        if len(line_style)==0:
                                                	line_style='-'
                                        	plt.plot(mass,np.log10(iso),color=color[j],linestyle=line_style,label=isotopes[j])
                                	plt.xlabel("M/M$_{\odot}$",fontsize=20)
                                	plt.ylabel('Log($X_{i})$', fontsize=20)
                                	plt.ylim(-8,0)
                                	plt.legend()
					if len(ylimits)>0:			
						plt.ylim(ylimits[0],ylimits[1])
					if len(mass_range)>0:
                                        	plt.xlim(mass_range[0],mass_range[1])
					if len(savepath)>0:
                                        	plt.savefig(savepath+'/'+str(cyc)+'.png')
						plt.close('all')
			if len(ylimits)>0:
				plt.ylim(ylimits[0],ylimits[1])
                               
	def surface_plot_2(self,runs=[],extra_label=[],x_axis='C+N',y_axis='Na',steps=1000,timestep=0,color=["r","b","k","g","r","b","k","g","r"],line_style=["-","--","-.","-.","--","-"]):

		'''

			x_axis : Either element .e.g. 'N', or sum of elements 'C+N'.
			y_axis : same
			
			Plots in the following form
			esilon(element/sumelement) = log10(N(element/sumelement)) +12
			
			if x_axis='age' : plot evolution over time
				model : model

		'''	
		import utils as u
                sefiles=[]
                legend=[]
                HDF5_surf=[]

		setting='iso_massf'
                if len(runs) ==0:
                        runs=self.run_dirs_name
                        if len(extra_label)==0:
                                extra_label=self.extra_label
		for i in range(len(self.runs_H5_surf)):
			sefiles.append(se(self.runs_H5_surf[i]))

		if 'age' in x_axis or 'model' in x_axis:
			species_inp=[y_axis.split('/')[0][1:],y_axis.split('/')[1][:-1]]
		else:
			if '+' in x_axis:
				x_values1=x_axis.split('+')
			else:
				x_values1=[x_axis]
			if '+' in y_axis:
				y_values1=y_axis.split('+')
			else:
				y_values1=[y_axis]
			species_inp=x_values1+y_values1

		
		for i in range(len(sefiles)):
			print "open",sefiles[i]

			x_abu=[]
			y_abu=[]
			star_mass=sefiles[i].get("mini")[0]
                        star_z=sefiles[i].get("zini")
			legend=str(star_mass)+"$M_{\odot}$, Z="+str(star_z)
			cycles=[]
			age=[]
			model=[]
			cycs_tot=sefiles[i].se.cycles
			if timestep>0:
				cycles.append(int(sefiles[i].se.cycles[0]))
				timesteps=timestep
				print 'retrieve time'
				time1=sefiles[i].get(range(0,len(cycs_tot),1),'age')
				for j in range(0,len(sefiles[i].se.cycles),1):
					if time1[j]>timesteps:
						cycles.append(int(cycs_tot[j]))
						age.append(time1[j])			
						timesteps+=timestep
						#print 'take cycle',cycles[j]
			else:
				#time1=sefiles[i].get(range(0,len(cycs_tot),1),'star_age')	
				for j in range(0,len(cycs_tot),steps):
					cycles.append(int(cycs_tot[j]))	
					#age.append(time1[j])
					age.append(sefiles[i].get(int(cycs_tot[j]),'age'))
				#print sefiles[i].get(cycs_tot,'age')
				print 'aaaaaaaaaaaage',age
				print cycles
			cycles11=cycles
			model=cycles
			print species_inp
			###get the data
			print 'retrive data'
			species_abu=[]
			for k in range(len(species_inp)):
				species_abu.append([0]*len(cycles))	
			print species_abu
			data=sefiles[i].get(cycles,setting)			
			print 'done'
			cycles=cycles11
			for ww in range(len(cycles11)):
				#print 'cycle',cycles11[ww]
				for k in range(len(sefiles[i].se.isotopes)):
					for t in range(len(species_inp)):
						if species_inp[t] in sefiles[i].se.isotopes[k]:
							#print sefiles[i].se.isotopes[k]
							#print 't',t,ww,k
							#print species_inp[t]
							#print species_abu
							#print 'before'
							#print species_abu
							#print species_abu[t][ww]
							species_abu[t][ww]+= float(data[ww][k])			
							#print 'after'
							#print species_abu
			print species_abu
			print len(species_abu)	
			###calculate the quotients


                        #x_axis_out= np.array(species_abu[0])/np.array(species_abu[1])
                        #y_axis_out= np.array(species_abu[2])/np.array(species_abu[3])

			if x_axis=='age':
				x_axis_out=age
				y_axis_out= np.log10( (species_abu[0]/species_solar[0]) ) - np.log10 ( (species_abu[1]/species_solar[1]) )			
		
			elif x_axis=='model':
				x_axis_out=model
				y_axis_out= np.log10( (species_abu[0]/species_solar[0]) ) - np.log10 ( (species_abu[1]/species_solar[1]) )
	

			else:
				fake_iso=[]
				for ele in species_inp:
					fake_iso.append(ele+'-99')
				
				u.convert_specie_naming_from_h5_to_ppn(fake_iso)
                        	Z=u.znum_int
				print 'zzz',Z
				print fake_iso
				n_x=0
				n_y=0
				#sum up in case its epsilon((C+N)
				for tt in range(len(x_values1)):
					n_x+=species_abu[tt]/(Z[tt]*2.)
				#Assume A ~ 2Z for element
                                for tt in range(len(x_values1),len(species_abu)):
                                        n_y+=species_abu[tt]/(Z[tt]*2.)
				

				x_axis_out= np.log10(n_x) +12
				y_axis_out= np.log10(n_y) +12 

			print x_axis_out
			print y_axis_out

			print "------------------"
			print "plot",legend
			fig1=plt.figure(333)
			plt.rcParams.update({'font.size': 16})
                	plt.rc('xtick', labelsize=16)
                	plt.rc('ytick', labelsize=16)
			ax=fig1.add_subplot(1,1,1)
			ax.plot(x_axis_out,y_axis_out,marker='<',label=legend)
			plt.xlabel('log $\epsilon$('+x_axis+')')
			plt.ylabel('log $\epsilon$('+y_axis+')')
			plt.legend()
	


	def surface_plot_1(self,runs=[],extra_label=[],x_axis='[Mg-25/Mg-24]',y_axis='[Mg-26/Mg-24]',steps=1000,timestep=0,color=["r","b","k","g","r","b","k","g","r"],line_style=["-","--","-.","-.","--","-"],marker=['o','D','s','v'],markersize=8*[6],iniabupath='/rpod2/critter/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn'):

		'''
			Plots [hs/Fe] vs [ls/Fe] and [Pb/hs] vs [ls/Fe] 
			iniabupath - set path of solar normalization file
			setting - do not change this mode

			if x_axis='age' : plot evolution over time
				model : model

		'''	
		import utils as u
                sefiles=[]
                legend=[]
                HDF5_surf=[]

		setting='iso_massf'
                if len(runs) ==0:
                        runs=self.run_dirs_name
                        if len(extra_label)==0:
                                extra_label=self.extra_label
		for i in range(len(self.runs_H5_surf)):
			sefiles.append(se(self.runs_H5_surf[i]))

		if 'age' in x_axis or 'model' in x_axis:
			species_inp=[y_axis.split('/')[0][1:],y_axis.split('/')[1][:-1]]
		else:
			species_inp=[x_axis.split('/')[0][1:],x_axis.split('/')[1][:-1],y_axis.split('/')[0][1:],y_axis.split('/')[1][:-1]]

		if species_inp[0][-1].isdigit():
			iso_inp=True
		else:
			iso_inp=False
		
		#get initial abu	
		species_solar=[0]*len(species_inp)
		a=u.iniabu(iniabupath)
		for iso in a.names:
			iso_namescheme=iso.capitalize()[:2].split()[0]+"-"+iso[3:].split()[-1]
			for k in range(len(species_inp)):
				if species_inp[k] in iso_namescheme:
					species_solar[k]+=a.habu[iso]

		print 'iniabu solar',species_solar	
		
		for i in range(len(sefiles)):
			print "open",sefiles[i]

			x_abu=[]
			y_abu=[]
			star_mass=sefiles[i].get("mini")[0]
                        star_z=sefiles[i].get("zini")
			legend=str(star_mass)+"$M_{\odot}$, Z="+str(star_z)
			cycles=[]
			age=[]
			model=[]
			cycs_tot=sefiles[i].se.cycles
			if timestep>0:
				cycles.append(int(sefiles[i].se.cycles[0]))
				timesteps=timestep
				print 'retrieve time'
				time1=sefiles[i].get(range(0,len(cycs_tot),1),'age')
				for j in range(0,len(sefiles[i].se.cycles),1):
					if time1[j]>timesteps:
						cycles.append(int(cycs_tot[j]))
						age.append(time1[j])			
						timesteps+=timestep
						#print 'take cycle',cycles[j]
			else:
				#time1=sefiles[i].get(range(0,len(cycs_tot),1),'star_age')	
				for j in range(0,len(cycs_tot),steps):
					print j
					cycles.append(int(cycs_tot[j]))	
					#age.append(time1[j])
					age.append(sefiles[i].get(int(cycs_tot[j]),'age'))
				#print sefiles[i].get(cycs_tot,'age')
				print 'aaaaaaaaaaaage',age
				print cycles
			cycles11=cycles
			model=cycles
			print species_inp
			###get the data
			print 'retrive data'
			species_abu=[]
			for k in range(len(species_inp)):
				species_abu.append([0]*len(cycles))	
			print species_abu
			data=sefiles[i].get(cycles,setting)			
			print 'done'
			cycles=cycles11
			for ww in range(len(cycles11)):
				#print 'cycle',cycles11[ww]
				for k in range(len(sefiles[i].se.isotopes)):
					for t in range(len(species_inp)):
						if species_inp[t] in sefiles[i].se.isotopes[k]:
							#print sefiles[i].se.isotopes[k]
							#print 't',t,ww,k
							#print species_inp[t]
							#print species_abu
							#print 'before'
							#print species_abu
							#print species_abu[t][ww]
							species_abu[t][ww]+= float(data[ww][k])			
							#print 'after'
							#print species_abu
			print species_abu
			print len(species_abu)	
			###calculate the quotients


                        #x_axis_out= np.array(species_abu[0])/np.array(species_abu[1])
                        #y_axis_out= np.array(species_abu[2])/np.array(species_abu[3])

			if x_axis=='age':
				print 'plot xasxis out'
				x_axis_out=age
				y_axis_out= np.log10( (species_abu[0]/species_solar[0]) ) - np.log10 ( (species_abu[1]/species_solar[1]) )			
		
			elif x_axis=='model':
				x_axis_out=model
				y_axis_out= np.log10( (species_abu[0]/species_solar[0]) ) - np.log10 ( (species_abu[1]/species_solar[1]) )
	

			else:
				x_axis_out= np.log10( (species_abu[0]/species_solar[0]) ) - np.log10 ( (species_abu[1]/species_solar[1]) )
				y_axis_out= np.log10( (species_abu[2]/species_solar[2]) ) - np.log10 ( (species_abu[3]/species_solar[3]) )

			print x_axis_out
			print y_axis_out

			print "------------------"
			print "plot",legend
			fig1=plt.figure(333)
			plt.rcParams.update({'font.size': 16})
                	plt.rc('xtick', labelsize=16)
                	plt.rc('ytick', labelsize=16)
                	if not x_axis=='age' and not x_axis=='model':
				plt.plot([0,0],[-10,10],linewidth=2,color="K")
                	plt.plot([-10,10],[0,0],linewidth=2,color="K")
			ax=fig1.add_subplot(1,1,1)
			ax.plot(x_axis_out,y_axis_out,color=color[i],marker=marker[i],markersize=markersize[i],linestyle=line_style[i],label=self.extra_label[i])
			plt.xlabel(x_axis)
			plt.ylabel(y_axis)
			plt.legend()
			


	def surface_plot(self,runs=[],extra_label=[],steps=1000,timestep=[],color=["r","b","k","g","r","b","k","g","r"],line_style=["-","--","-.","-.","--","-"],setting='iso_massf',iniabupath='/home/christian/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn'):

		'''
			Plots [hs/Fe] vs [ls/Fe] and [Pb/hs] vs [ls/Fe] 
			iniabupath - set path of solar normalization file
			setting - do not change this mode

		'''	
		import utils as u
                sefiles=[]
                legend=[]
                HDF5_surf=[]
                if len(runs) ==0:
                        runs=self.run_dirs_name
                        if len(extra_label)==0:
                                extra_label=self.extra_label
		for i in range(len(self.runs_H5_surf)):
			sefiles.append(se(self.runs_H5_surf[i]))

		hs=["Ba","La","Nd","Sm"]
		ls=["Sr","Y-","Zr"]
		#get initial abu	
		ba_solar=0
		la_solar=0
		nd_solar=0
		sm_solar=0
		sr_solar=0
		y_solar=0
		zr_solar=0
		fe_solar=0
		pb_solar=0
		a=u.iniabu(iniabupath)
		for iso in a.names:
			iso_namescheme=iso.capitalize()[:2].split()[0]+"-"+iso[3:].split()[-1]	
			if "Ba-" in iso_namescheme:
				ba_solar+=a.habu[iso]
			if "La-" in iso_namescheme:
				la_solar+=a.habu[iso]
			if "Nd-" in iso_namescheme:
				nd_solar+=a.habu[iso]
			if "Sm-" in iso_namescheme:
				sm_solar+=a.habu[iso]
			if "Sr-" in iso_namescheme:
				sr_solar+=a.habu[iso]
			if "Y-" in iso_namescheme:
				y_solar+=a.habu[iso]
			if "Zr-" in iso_namescheme:
				zr_solar+=a.habu[iso]
			if "Fe-" in iso_namescheme:
				fe_solar+=a.habu[iso]
			if "Pb-" in iso_namescheme:
				pb_solar+=a.habu[iso]

		for i in range(len(sefiles)):
			print "open",sefiles[i]

			ls_fe_abu=[]
			hs_ls_abu=[]
			pb_hs_abu=[]
			hs_fe_abu=[]
			star_mass=sefiles[i].get("mini")
                        star_z=sefiles[i].get("zini")
			legend=str(star_mass)+"$M_{\odot}$, Z="+str(star_z)
			cycles=[]
			if timestep>0:
				starage=sefiles[i].get('age')
				cycles.append(sefiles[i].se.cycles[0])
				timesteps=timestep
				for j in range(0,len(sefiles[i].se.cycles),1):
					print starage[j]
					if starage[j]>timesteps:
						cycles.append(sefiles[i].se.cycles[j])				
						timesteps+=timestep
			else:
				for j in range(0,len(sefiles[i].se.cycles),steps):
					cycles.append(sefiles[i].se.cycles[j])					
			

			ba_abu=np.zeros(len(cycles))
			la_abu=np.zeros(len(cycles))
			nd_abu=np.zeros(len(cycles))
                        sm_abu=np.zeros(len(cycles))
                        sr_abu=np.zeros(len(cycles))
                        y_abu=np.zeros(len(cycles))
                        zr_abu=np.zeros(len(cycles))
                        fe_abu=np.zeros(len(cycles))
			pb_abu=np.zeros(len(cycles))	
			
			k=-1
			iso_abu=sefiles[i].get(cycles,setting)
			for cyc in cycles:
				k+=1
				print "read cycle",cyc
				if "iso" in setting:
					hh=-1
					for iso_namescheme in sefiles[i].se.isotopes:
						hh+=1
						if "Ba-" in iso_namescheme:
							ba_abu[k]+=iso_abu[k][hh]
						if "La-" in iso_namescheme:
							la_abu[k]+=iso_abu[k][hh]
						if "Nd-" in iso_namescheme:
							nd_abu[k]+=iso_abu[k][hh]
						if "Sm-" in iso_namescheme:
							sm_abu[k]+=iso_abu[k][hh]
						if "Sr-" in iso_namescheme:
							sr_abu[k]+=iso_abu[k][hh]
						if "Y-" in iso_namescheme:
							y_abu[k]+=iso_abu[k][hh]
						if "Zr-" in iso_namescheme:
							zr_abu[k]+=iso_abu[k][hh]
						if "Fe-" in iso_namescheme:					
							fe_abu[k]+=iso_abu[k][hh]
						if "Pb-" in iso_namescheme:
							pb_abu[k]+=iso_abu[k][hh]
				else:
					for iso_namescheme in sefiles[i].se.elements:
                                                if "Ba" in iso_namescheme:
                                                        ba_abu[k]+=iso_abu[k][hh]
                                                if "La" in iso_namescheme:
                                                        la_abu[k]+=iso_abu[k][hh]
                                                if "Nd" in iso_namescheme:
                                                        nd_abu[k]+=iso_abu[k][hh]
                                                if "Sm" in iso_namescheme:
                                                        sm_abu[k]+=iso_abu[k][hh]
                                                if "Sr" in iso_namescheme:
                                                        sr_abu[k]+=iso_abu[k][hh]
                                                if "Y" == iso_namescheme:
                                                        y_abu[k]+=iso_abu[k][hh]
                                                if "Zr" in iso_namescheme:
                                                        zr_abu[k]+=iso_abu[k][hh]
                                                if "Fe" in iso_namescheme:
							fe_abu[k]+=iso_abu[k][hh]
						#pb needs to be added here		
	
				#print fe_abu[k],sr_abu[k],zr_abu[k],y_abu[k]
				#print sr_solar,zr_solar,y_solar,fe_solar
				print "PB:",pb_solar,pb_abu[k]
				ls_fe_abu.append( np.log10((fe_solar/fe_abu[k])*((sr_abu[k] * zr_abu[k] * y_abu[k]/(sr_solar * zr_solar * y_solar))**(1./3.) )))
				#print ((fe_solar/fe_abu[k])*((sr_abu[k] * zr_abu[k] * y_abu[k]/(sr_solar * zr_solar * y_solar))**(1./3.) ))
				ls_ba=np.log10((ba_solar/ba_abu[k])*((sr_abu[k] * zr_abu[k] * y_abu[k]/(sr_solar * zr_solar * y_solar))**(1./3.) ))
				ls_la=np.log10((la_solar/la_abu[k])*((sr_abu[k] * zr_abu[k] * y_abu[k]/(sr_solar * zr_solar * y_solar))**(1./3.) ))
				ls_nd=np.log10((nd_solar/nd_abu[k])*((sr_abu[k] * zr_abu[k] * y_abu[k]/(sr_solar * zr_solar * y_solar))**(1./3.) ))
				ls_sm=np.log10((sm_solar/sm_abu[k])*((sr_abu[k] * zr_abu[k] * y_abu[k]/(sr_solar * zr_solar * y_solar))**(1./3.) ))
				hs_fe=np.log10((fe_solar/fe_abu[k])*((ba_abu[k] * la_abu[k] * nd_abu[k] * sm_abu[k]/( ba_solar * la_solar * nd_solar*sm_solar))**(1./4.) ))	
				hs_pb=np.log10((pb_solar/pb_abu[k])*((ba_abu[k] * la_abu[k] * nd_abu[k] * sm_abu[k]/( ba_solar * la_solar * nd_solar*sm_solar))**(1./4.) ))
				hs_fe_abu.append(hs_fe)
				hs_ls_abu.append(hs_fe - ls_fe_abu[-1])   #( -1*(ls_ba) + -1*(ls_la) + -1*(ls_nd) + -1*(ls_sm))/4.)
				pb_hs_abu.append( -1* hs_pb )
				#ls_fe_abu.append( np.log10( 1./3. *( sr_abu[k]/sr_solar + y_abu[k]/y_solar + zr_abu[k]/zr_solar)/(fe_abu[k]/fe_solar))) 
				#hs_ls_abu.append( np.log10( 1./4. *( ba_abu[k]/ba_solar + la_abu[k]/la_solar + nd_abu[k]/nd_solar + sm_abu[k]/sm_solar) / ( 1./3. *( sr_abu[k]/sr_solar + y_abu[k]/y_solar + zr_abu[k]/zr_solar)))) 
				


			print hs_ls_abu
			print "------------------"
			print ls_fe_abu
			print "plot",legend
			fig1=plt.figure(333)
			plt.rcParams.update({'font.size': 16})
                	plt.rc('xtick', labelsize=16)
                	plt.rc('ytick', labelsize=16)
                	plt.plot([0,0],[-10,10],linewidth=2,color="K")
                	plt.plot([-10,10],[0,0],linewidth=2,color="K")
			ax=fig1.add_subplot(1,1,1)
			ax.plot(ls_fe_abu,hs_fe_abu,marker='<',label=legend,markevery=500)
			plt.xlabel("[ls/Fe]")
			plt.ylabel("[hs/Fe]")
			plt.legend()
			
			fig2=plt.figure(555)
                        plt.rcParams.update({'font.size': 16})
                        plt.rc('xtick', labelsize=16)
                        plt.rc('ytick', labelsize=16)
                        plt.plot([0,0],[-10,10],linewidth=2,color="K")
                        plt.plot([-10,10],[0,0],linewidth=2,color="K")
			ax2=fig2.add_subplot(1,1,1)
			ax2.plot(ls_fe_abu,pb_hs_abu,marker='<',label=legend,markevery=500)
			plt.xlabel("[ls/Fe]")
			plt.ylabel("[Pb/hs]")
			plt.legend()

			fig2=plt.figure(888)
                        plt.rcParams.update({'font.size': 16})
                        plt.rc('xtick', labelsize=16)
                        plt.rc('ytick', labelsize=16)
                        plt.plot([0,0],[-10,10],linewidth=2,color="K")
                        plt.plot([-10,10],[0,0],linewidth=2,color="K")
			ax2=fig2.add_subplot(1,1,1)
			ax2.plot(hs_ls_abu,ls_fe_abu,marker='<',label=legend,markevery=500)
			plt.xlabel("[hs/ls]")
			plt.ylabel("[ls/Fe]")
			plt.legend()




	def pocket_composition(self,mp,sefiles,cycle,massbot,masstop,isotopes,label,legend,color,title):  ###hdf5out
        	'''

		mass_range    - required to plot data in a certain mass range. Needed for read_iso_abund_marco
        	cycle         - which cycle from the h5 file?. Needed for read_iso_abund_marco
        	stable        - logic if want to plot only stable or not.  
        	i_decay       - if = 1 I plot not decayed, if = 2 I plot decayed. Make sense only if stable is true.


		'''


		mass_range=[massbot,masstop]
		sefiles.average_iso_abund_marco(mass_range,cycle,stable=False,i_decay=1)
		mass=[]
		plotiso=[]
		startyields=[]
		plotiso_massfrac=[]
		for i in range(len(isotopes)):
			startyields.append(sefiles.get(sefiles.se.cycles[0],isotopes[i])[0])
		for j in range(len(sefiles.se.isotopes)):
			if sefiles.se.isotopes[j] in isotopes:
				plotiso.append(mp.average_mass_frac[j])
				mass.append(sefiles.se.A[j])
		for i in range(len(isotopes)):
			plotiso_massfrac.append(plotiso[i]/startyields[i])
		plt.plot(mass,plotiso_massfrac,marker='*',markersize=8,mfc=color,linestyle='None',label=legend)
		#plt.plot(mass,plotiso_massfrac,marker='*')
		if label ==True:
			for j in range(len(isotopes)):
				plt.annotate(isotopes[j], xytext = (0, 10),textcoords = 'offset points' ,xy=(mass[j], plotiso_massfrac[j]))
		mass=np.array(mass)
		plt.xlim(mass.min()-4,mass.max()+4)
		plt.legend()
		plt.title(title)
		plt.xlabel("mass number")
		plt.ylabel("Isotopic production factors")

	def get_elem_abu(self,cycles,element,sefiles):
		'''
			Get abundance of element without using elem_massf since
			its not working correctly...
			reads also hdf5 out files
		'''
		isotopes=sefiles.se.isotopes
		if type(cycles)==int:
			if cycles==-1:
				cycles=int(sefiles.se.cycles[-1])
			len_cycles=1
		else:
			len_cycles=len(cycles)
		type_hdf5= type(sefiles.get(sefiles.se.cycles[0],"H-1"))
		if type_hdf5 == np.float64:
		 	abu=np.zeros(len_cycles)
			for i in range(len(isotopes)):
				if element == isotopes[i].split('-')[0]: 	
					abu+=sefiles.get(cycles,isotopes[i])
		## hdf5 out
		#else:
			#abu

		return abu	





########################################################################################
############################ START OF MESA SPECIFIC CONTENT#############################
########################################################################################








class mesa_set(history_data):

	'''
		This class allows access to multiple MESA runs in the path dir.
		Instead of defining one path, an array of paths , mult_dir can be
		defined. In this case each path in multi_dir represents a run.

		rundir - dir with run directory, means this dir contains MESA workdirs
			 Only the LOGS dir inside each MESA workdir is processed.
		multi_dir - specify either run dirs (MESA workdir structure) or LOGS dirs at different locations
		extra_label - set labels of runs, else run dir name will be chosen	
	
		e.g. inside a dir containing MESA workdirs:
	
			mesaset=set.mesa_set(".",extra_label=["double-f","reference"])		
			
		or read certain MESA workdirs at certain locations: 
				
			mesaset=set.mesa_set(multi_dir=["M3.00Z2.0e-02_1_re_double_f","M3.00Z2.0e-02_re_referencerun"],extra_label=["double-f","reference]])

		or to specifiy certain LOGS dirs:

			mesaset=set.mesa_set(multi_dir=["M3.00Z2.0e-02/LOGS","M3.00Z2.0e-02/LOGS1"],extra_label=["double-f","reference]])



		Using the set_plot* functions you can instantly plot diagrams
		for all runs defined in mesa set
		e.g.

                	mesaset.set_plot_kipp_CO()


		List of available functions:


		Plotting:

		def set_plot_prof: Plots profile of isotopes for cycles.

		def set_plot_hrd: plots HRDS

		def set_plot_t_mass: plot star mass evolution over time

		def set_plot_model_mass: plot star mass evolution vs model number

		def set_plot_CO: Plots C/O ratio vs model number or age

		def set_plot_CO_mass: Plots C/O ratio vs star mass

		def set_plot_vsurf: Plots the surface velocity vs model number

		def set_plot_timesteps: Plots the timestep size vs model number

		def set_plot_mdot: Plots the mass loss vs time, model or mass

		def set_plot_R(self): Plots the stellar logarithmic radius vs model number

		def set_plot_kip_special: Kippenhahn which plots only h and he free bndry,

		def set_plot_kipp_cont: Plots kipp_cont diagrams

		def set_plot_kipp_CO: plots kippenhahn diagrams with c/o ratio

		def set_plot_surfabu(self): Plot surface abundances of t_surfabu in nugridse

		def set_plot_lumi: Plot luminosity vs. t

		def set_thing: Plots of what the user specifies

		def TPAGB_core_growth: Creates diagram of the core growth during AGB phase.,
                         Uses function multi_DUP to identify first TP cycle

		
		Tables:
		
		def lifetime: Calculate stellar lifetime till first TP dependent of initial mass
                		uses self.set_find_first_TP()

		def final_shells: Return table of different core masses (h1,he4,final).



		Other:



		def set_plot_infall_v: calculates the infall velocity vs model number
                        		uses forum.astro.keele.ac.uk_update/utils/explosion


		def explosion: Calculates synthetic explosions, delayed and rapid for
                        	typs "dummy" and "inter" for massive stars (as in set1)


		def multi_DUP: Plot (!) DUP parameter lambda but also returns peak lum and corresponding
                        model number

		def HDUP_profiles: Profil plot for analyzing the diffusion coefficients of the HDUP

		def HDUP_profiles1: Profile plot to analyze the HDUP H burning...



	'''

	def __init__(self,rundir='.',multi_dir=[],extra_label=99*[' '],newstarlog=False):

		import os.path
		if len(multi_dir)==0:
			#slist = os.listdir(rundir)
			slist = [d for d in os.listdir(rundir) if os.path.isdir(os.path.join(rundir, d))]
		else:
			slist = multi_dir
		self.runs_H5=[]
		self.run_LOGS=[]
		self.run_historydata=[]
		self.run_label=[]
		self.multi_dir=multi_dir
		i=0
		if newstarlog==None:
			clean_1=raw_input("Create new star.logsa/historydata.sa for all runs? [y/n]")
			if clean_1 == 'yes' or clean_1 == 'y':
				clean_starlog=True
			else:
				clean_starlog=False
		else:
			if newstarlog:
				clean_starlog=True
			else:
				clean_starlog=False	

		#sort 
		print 'slist',slist
		#structure? M1.00Z2.0e-02
		if len(multi_dir) == 0:
			struc=True
			masses=[]
			for element in slist:
				if (not 'M' in element) or (not 'Z' in element):
					struc=False
					break
				else:
					masses.append(float(element.split('Z')[0][1:]))		
			if struc == True:
				slist1=[]
				idx=sorted(range(len(masses)), key=lambda k: masses[k])
				for i in range(len(slist)):
					slist1.append(slist[idx[i]])
				slist=slist1
				print 'Sorted: ',slist
		i=0
		print 'test',slist
		print 'multi_dir',multi_dir
		for element in slist:
			print element
			#i = i + 1
			print i
			if len(multi_dir)==0:
				run_path=rundir+"/"+element
			else:
				if multi_dir[i][0] == "/":
					run_path=multi_dir[i]
				else:
					if len(rundir)==0:
						run_path=os.getcwd()+"/"+multi_dir[i]
					else:
						run_path=rundir+'/'+multi_dir[i]
			print 'runpath',run_path
			if os.path.isdir(run_path):
				run_path_1=glob.glob(run_path+"/*/*.h5")
				if len(run_path_1)>0:
					h5_dir=run_path_1[0].split("/")[-2]	
					self.runs_H5.append(run_path+"/"+h5_dir)
				else:
					print "Warning: h5 files are not available in "+run_path
				if (len(glob.glob(run_path+"/*/*.data"))>0) or (len(glob.glob(run_path+"/*/*.log"))>0):	
					self.run_dir=run_path[:-len(run_path.split('/')[-1])]
					
					self.run_LOGS.append(run_path+"/LOGS")
					print "Read "+run_path
					self.run_historydata.append(history_data(run_path+"/LOGS",clean_starlog=clean_starlog))
					if len(extra_label)>0:
						self.run_label.append(self.create_label(self.run_historydata[-1],extra_label[i]))
					else:
						self.run_label.append(self.create_label(self.run_historydata[-1],element))
					i+=1
				#if rundir is already LOGS dir
				elif (len(glob.glob(run_path+"/*.data"))>0) or (len(glob.glob(run_path+"/*.log"))>0):
					 self.run_LOGS.append(run_path)	
					 print "Read "+run_path
 					 self.run_historydata.append(history_data(run_path,clean_starlog=clean_starlog))
					 if len(extra_label)>0:
						self.run_label.append(self.create_label(self.run_historydata[-1],extra_label[i]))
					 else:
						self.run_label.append(self.create_label(self.run_historydata[-1],element))


				else:
					#if len(multi_dir)>=0:
					print "Error: not history.data or star.log file found in "+run_path		

				self.symbs=utils.symbol_list('lines1')


        def set_plot_infall_v(self,star_mass_1=7.,plots_dir='.'):
                '''
			Plots infaall velocity over model number
                        plots_dir - dir for plots
			Uses scripts from
			forum.astro.keele.ac.uk/utils/explosion

                '''
                for i in range(len(self.run_historydata)):

                        #star_mass=mesa_profile(self.run_LOGS[i],1,num_type="profile_num").header_attr['initial_mass']
			#star_z=mesa_profile(self.run_LOGS[i],1,num_type="profile_num").header_attr['initial_z']
                        star_mass=self.run_historydata[i].header_attr['initial_mass']
			star_z=self.run_historydata[i].header_attr['initial_z']
			if star_mass != star_mass_1  :
                                continue

                        #:wos.chdir(rundir_path)

                        #calculate collapse cycle
                        #1. ncycle,2.nramp,3.vsinput,4.coolingcycles,(5.cooolingcyclleswrite)

                        ncycle=21973
                        nramp=10
                        vsinput=2.e9
                        cooling_cycles=290

                        ##calculate collapse
                        import imp
                        import explosive_main as exp
                        #pylib_path= str(__file__)[:-(len("nugrid_set.py")+1)]
                        #exp=imp.load_source("explosive_main",pylib_path+'/../explosion/explosive_main.py')
                        calc=exp.explosion(run_H5=self.runs_H5[i],LOGS_path=self.run_LOGS[i],nramp=nramp,vsinput=vsinput,se_load='n')
			
			plot_name="M"+str(star_mass)+"Z"+str(star_z)
			plt.figure(plot_name+".png")
			cycle_collapse=calc.collapse_cycle_calc()
			plt.savefig(plot_name+".png")

	def explosion_v_traj(self,traj_dir='.'):
		'''

			Writes out shock (velocity) profiles.
			Calculates synthetic explosions, delayed and rapid for 
			typs "dummy" and "inter" for massive stars (as in set1)

			Use of explosion_main tool in explosion dir
			and from this calc_v_shock_prof

			see_exp_dir - dir for calculated h5 files
			traj_dir - dir for plots

                        Uses scripts from
                        forum.astro.keele.ac.uk/utils/explosion


		'''
		for i in range(len(self.run_historydata)):
			
			####Massie stars defined ast M<=8M
			star_mass=mesa_profile(self.run_LOGS[i],1).header_attr['initial_mass']
			if star_mass <= 8.  :
				continue
			rundir_path=traj_dir+'/'+self.run_LOGS[i].split("/")[-2]
			
			#:wos.chdir(rundir_path)
				
			#calculate collapse cycle
			#1. ncycle,2.nramp,3.vsinput,4.coolingcycles,(5.cooolingcyclleswrite)
			
			ncycle=21973
			nramp=10
			vsinput=2.e9
			cooling_cycles=290

			#Find the massive stars from the set
			##
			
			##calculate collapse
			import imp
			import explosive_main as exp
			#pylib_path= str(__file__)[:-(len("nugrid_set.py")+1)]
			#exp=imp.load_source("explosive_main",pylib_path+'/../explosion/explosive_main.py')
			calc=exp.explosion(run_H5=self.runs_H5[i],LOGS_path=self.run_LOGS[i],nramp=nramp,vsinput=vsinput)			
			for delay in [True,False]:
				if delay:
					run=rundir_path+"_prof_delay.txt"
				else:
					run=rundir_path+"_prof_rapid.txt"
				calc.calc_v_shock_prof(delay=True,output_file=run)

	
	def explosion(self,see_exp_dir='.',plots_dir='.'):
		'''
			Calculates synthetic explosions, delayed and rapid for 
			typs "dummy" and "inter" for massive stars (as in set1)

			Use of explosion_main tool in explosion dir

			see_exp_dir - dir for calculated h5 files
			plots_dir - dir for plots

                        Uses scripts from
                        forum.astro.keele.ac.uk/utils/explosion


		'''
		for i in range(len(self.run_historydata)):
			
			####Massie stars defined ast M<=8M
			star_mass=mesa_profile(self.run_LOGS[i],1).header_attr['initial_mass']
			if star_mass <= 8.  :
				continue
			rundir_path=see_exp_dir+'/'+self.run_LOGS[i].split("/")[-2]
			
			#:wos.chdir(rundir_path)
				
			#calculate collapse cycle
			#1. ncycle,2.nramp,3.vsinput,4.coolingcycles,(5.cooolingcyclleswrite)
			
			ncycle=21973
			nramp=10
			vsinput=2.e9
			cooling_cycles=290

			#Find the massive stars from the set
			##
			
			##calculate collapse
			import imp
			import explosive_main as exp
			#pylib_path= str(__file__)[:-(len("nugrid_set.py")+1)]
			#exp=imp.load_source("explosive_main",pylib_path+'/../explosion/explosive_main.py')
			calc=exp.explosion(run_H5=self.runs_H5[i],LOGS_path=self.run_LOGS[i],nramp=nramp,vsinput=vsinput)			
			for delay in [True,False]:
				if delay:
					run=rundir_path+".delay"
				else:
					run=rundir_path+".rapid"
				if not os.path.isdir(run):
                                	os.mkdir(run)
                                print "create exp run dir ",run

				print '############Calculate dummy is',delay
				calc.exp_dummy(delay=delay,write_dir=run,extra_name='')
				calc.exp_inter(delay=delay,write_dir=run,extra_name='')
			#create diagrams in plots_dir dir
				calc.plot_rho_T(delay=delay,plots_dir=plots_dir,energy_exp=2.0e51)			
				calc.plot_exp_details(delay=delay,plots_dir=plots_dir,model_label='')

	def set_plot_core_lum_relation(self,symb='o',linestyle='--',sparsity=500,label=''):

		'''
			Plots the core luminosity relation. Mean core mass
			from first TP till last TP in model steps of sparsity
			Mean lum derived first TP to last TP in model steps of sparsity
			time-weighted averages
		'''


                m=self.run_historydata
		core_mean=[]
		lum=[]
                for case in m:
	        	peak_lum_model,h1_mass_min_DUP_model=case.find_TP_attributes(t0_model=case.find_first_TP(),fig=3, color='r', marker_type='o')
			modelrange=range(int(peak_lum_model[0]),int(peak_lum_model[-1]+sparsity),int(sparsity))
			core_mean1=[]
			lum1=[]
			#weight for mean value
			ages=case.get('star_age')
			weight_sum=0.
			for h in range(len(modelrange)):
				mod=modelrange[h]
				if h==0:
					weight=(ages[modelrange[h+1]] - ages[mod])/2.
				elif h == len(modelrange)-1:
					weight=(ages[mod]-ages[modelrange[h-1]])/2.
				else:
					weight=(ages[modelrange[h+1]] - ages[mod])/2. + (ages[mod]-ages[modelrange[h-1]])/2.
				weight_sum+=weight

				core_mean1.append(case.get('h1_boundary_mass')[mod]*weight)
				lum1.append(case.get('log_L')[mod]*weight)
			
			core_mean.append(sum(core_mean1)/weight_sum)
			lum.append(sum(lum1)/weight_sum)
			
		plt.plot(core_mean,lum,marker=symb,linestyle=linestyle,label=label)
		print core_mean,lum
		plt.xlabel("M/M$_{\odot}$")
		plt.ylabel('log L/L$_{\odot}$')

	def set_plot_DUP(self,symbs_1=[],linestyle=[]):

                m=self.run_historydata
                i=0
                if len(symbs_1)>0:
                        symbs=symbs_1
                else:
                        symbs=self.symbs
                if len(linestyle)==0:
                        linestyle=200*['-']
                for case in m:
			TP_models,DUP_models,TPend,lambds = case.find_TPs_and_DUPs()
			plt.plot(range(len(lambds)),lambds,label=self.run_label[i],marker=symbs[i],linestyle=linestyle[i])	
			plt.legend()
			plt.figure(0)
                        plt.plot(range(len(lambds)),lambds,label=self.run_label[i],marker=symbs[i],linestyle=linestyle[i])
                        plt.legend()			
			i+=1


	def multi_DUP(self,t0_model=[],h_core_mass=False,plot_fig=True,linestyle=[],marker=[],color=[]):
		'''
			Plot (!) DUP parameter lambda versus mass.
			There is a algorithm finding the first TP (and set t0_model) for reach run dir 
			which worked fine for set1. But there could be a problem with SAGB stars
			so use with caution! You can use a manual input (the array t0_model) to
			skip this step.

			It also returns peak lum and corresponding
			model number of each TP.


			h_core_mass - If True: plot dependence from h free core , else star mass 
			t0_model -  Array of first TPs of models.			
		
			examples of set1:
			t0_model paramter:		
			z1e-2:1.65-5: [13690,3120,3163,5306]
			z2e-2:1.65-5: [16033,6214,3388,5368]

			dirs=["M1.65Z2.0e-02/LOGS","M2.00Z2.0e-02/LOGS","M3.00Z2.0e-02/LOGS","M5.00Z2.0e-02/LOGS"]
			dirs=["M1.65Z1.0e-02/LOGS","M2.00Z1.0e-02/LOGS","M3.00Z1.0e-02/LOGS","M5.00Z1.0e-02/LOGS"]
			note: M3_1e-2 shows decrease in core mass in the end
			set.multi_DUP(t0_model=[13690,3120,3163,5306])
			set.multi_DUP(t0_model=[16033,6214,3388,5368])
		
		'''
		if len(t0_model)==0:
			t0_model=self.set_find_first_TP()
		dirs=self.run_LOGS
		historydata=self.run_historydata
	
		marker_type=marker
		#line_style=['--','-','-.',':']
		peak_lum_model_array=[]
		h1_mass_min_DUP_model_array=[]
		for i in range(len(dirs)):
			#color,marker_type)
			peak_lum_model,h1_mass_min_DUP_model = historydata[i].find_TP_attributes(0,t0_model[i],color[i],marker_type[i],h_core_mass,no_fig=plot_fig)			
			peak_lum_model_array.append(peak_lum_model)
			h1_mass_min_DUP_model.append(h1_mass_min_DUP_model)

		return peak_lum_model_array,h1_mass_min_DUP_model_array

				
	###the following methods allow are parts of Falks vis3.py file, thanks Falk


	def HDUP_profiles(self,fig=0,cycles=[],conv_env_mass_range=False,markevery=50,xlim=[-1e-3,1e-3]):
		
		'''
			Profil plot for analyzing the diffusion coefficients of the HDUP

		'''	

                m=self.run_historydata
                i=0
                for case in m:
			
			fig = plt.figure(fig) #figsize=(24.0, 10.0),)
			fig.subplots_adjust(right=0.75)
			ax = fig.add_subplot(111)	
			profile=mesa_profile(self.run_LOGS[i],cycles[i])	
			modelnumber=profile.header_attr["model_number"]
			idx=modelnumber	
			modelnumbers=case.get("model_number")
			for t in range(len(modelnumbers)):		
				if modelnumbers[t]==modelnumber:
					idx=t
			print 'model',modelnumbers[idx],'vs cycle ',cycles[i]
			conv_envelope_mass=case.get("mx1_bot")[idx]*case.get("star_mass")[idx]
			print 'conv env',conv_envelope_mass
			gradr=profile.get("gradr")
			grada=profile.get("grada")
			#Find sb boudnasry
			m1=profile.get("mass")
			idx=min(range(len(m1)), key=lambda i: abs(m1[i]-conv_envelope_mass))
			print 'idx',idx
			if conv_env_mass_range ==True:
				for k in range(idx,0,-1):
					if abs(gradr[k] - grada[k])<0.0001:
						sb_mass=m1[k]
						break
				print 'short sb',sb_mass
				mass=profile.get("mass") - sb_mass
			else:
				mass=m1
			D_mlt=profile.get("log_mlt_D_mix")
			D = profile.get("log_D_mix")
			
			ax2 = plt.twinx()
			ax2.set_ylabel("grad(T)")	
			ax2.plot(mass,gradr,label="$grad(T)_{rad}$",color="yellow")
			ax2.plot(mass,grada,label="$grad(T)_{ad}$")
			ax2.set_ylim(0, 1)
			ax2.legend(loc=4)
			#ax2.grid()
			if conv_env_mass_range == True:
				ax.plot([conv_envelope_mass-sb_mass,conv_envelope_mass-sb_mass],[-25,45],color='b',linestyle=":")
				ax.plot([0,0],[-25,45],linestyle='--',color='r')

			ax.plot(mass,D_mlt,label="$D_{mlt}$",linestyle='-.',color='blue',marker='o',markevery=markevery)
			ax.plot(mass,D,label="D",color="red",marker='s',linestyle='-',markevery=markevery)
			#plt.title(title)
			ax.legend()
			ax.set_ylim(0,20)
			if conv_env_mass_range ==True:
				ax.set_xlabel("M - m$_{SB}$ [$M_{\odot}$]")
			else:
	 	                ax.set_xlabel("M [$M_{\odot}$]")
	
			ax.set_ylabel(r"Diffusion coefficient D [$\frac{cm^2}{s}$]")
		if conv_env_mass_range == True:
			#ax2.set_xlim(-5e-05+conv_envelope_mass,conv_envelope_mass+5.e-5)
			ax2.set_xlim(xlim[0],xlim[1])
		#ax.set_ylim(-5, 20)
			print conv_envelope_mass


	def HDUP_profiles_1(self,cycles=[],conv_env_mass_range=False,diff=False,sb_mass_input=0.975,xlim=[-1e-3,1e-3],const_log=-12):

		'''
			Profile plot to analyze the HDUP H burning...
		'''


		from mpl_toolkits.axes_grid1 import host_subplot
		import mpl_toolkits.axisartist as AA
		


                m=self.run_historydata
                i=0
                for case in m:
			if diff==True:
				fig = host_subplot(111, axes_class=AA.Axes)
				plt.subplots_adjust(right=0.75)
				ax=plt.gca()
				plt.tick_params(axis='both', which='major', labelsize=22)
				plt.tick_params(axis='both', which='minor', labelsize=22)

			else:
				fig = plt.figure(i+44)
				fig.subplots_adjust(right=0.75)
                        	ax = fig.add_subplot(111)
                        profile=mesa_profile(self.run_LOGS[i],cycles[i])
			modelnumber=profile.header_attr["model_number"]
			idx=modelnumber	
			modelnumbers=case.get("model_number")
			for k in range(len(modelnumbers)):		
				if modelnumbers[k]==modelnumber:
					idx=k
			print 'model',modelnumbers[idx],'vs cycle ',cycles[i]
			conv_envelope_mass=case.get("mx1_bot")[idx]*case.get("star_mass")[idx]
			print 'conv env',conv_envelope_mass
                        gradr=profile.get("gradr")
                        grada=profile.get("grada")
			m1=profile.get("mass")
			idx=min(range(len(m1)), key=lambda i: abs(m1[i]-conv_envelope_mass))
			print 'idx',idx
			if conv_env_mass_range ==True:
				for k in range(idx,0,-1):
					if abs(gradr[k] - grada[k])<0.0001:
						sb_mass=m1[k]
						break
				sb_mass=sb_mass_input
				print 'short sb',sb_mass
				mass=profile.get("mass") - sb_mass
			else:
				mass=m1	
			h1=profile.get("h1")
			c12=profile.get("c12")
			logT=profile.get("logT")
			energygen=profile.get("eps_nuc")
			entropy=profile.get("entropy")
			temp=profile.get('temperature')	
			logtemp=np.log10(np.array(temp))
			ax.plot(mass,np.log10(h1+1.e-99),label="X($^1H$)",color="blue",marker='o',linestyle='-',markevery=30)
			ax.plot(mass,np.log10(c12+1.e-99),label="X($^{12}C$)",color="red",marker='s',linestyle='--',markevery=30)
			ax.set_ylim(-4,0)
			plt.ylim(-4,0)
			ax2 = plt.twinx()
			ax2.set_ylabel(r'$log(\epsilon_{nuc} [\frac{erg}{g*s}])$',size=20)
			ax2.plot(mass,np.array(np.log10(energygen+1.e-99)),label="$\epsilon_{nuc}$",color="g",marker='*',linestyle='-',markevery=30)
			ax2.set_ylim(5,15)
			if diff==True:
				print 'for ax3'
				profile1=mesa_profile(self.run_LOGS[i],cycles[i]-1)
	                        D_mlt=profile1.get("log_mlt_D_mix")
        	                D = profile.get("log_D_mix")
				ax3=ax.twinx()
				offset = 60
				new_fixed_axis = ax3.get_grid_helper().new_fixed_axis
				ax3.axis["right"] = new_fixed_axis(loc="right",
								axes=ax3,
							offset=(offset, 0))
				#ax3.plot(mass,D_mlt,label='$D_{mlt}$')
				ax3.plot(mass,D,label='$D$',marker='^',markevery=30,color='k')
				ax3.set_ylabel('$log (D [cm^2/s])$')
				ax3.set_ylim(0,20)
				#ax=ax3
			if conv_env_mass_range ==True:
			
                       		#ax.plot([conv_envelope_mass-sb_mass,conv_envelope_mass-sb_mass],[-25,45],color='b',linestyle=":")
                        	ax.plot([0,0],[-25,45],linestyle='--',color='r')
				ax.set_xlabel("M - m$_{SB}$ [$M_{\odot}$]",size=20)
			else:
				ax.set_xlabel("M  [$M_{\odot}$]",size=20)
			ax.set_ylabel(r"$log (X_{i})$",size=20)
			ax.legend(loc=4)
			
			#ax2 = plt.twinx()
			#ax2.set_ylabel(r"Specific entropy [$\frac{erg}{mol*K}/(N_AK_B)$]")
			#ax2.plot(mass,entropy,label="entropy",color="green")
			#ax2.plot(mass,logtemp,label='T',marker='D',markevery=30,linestyle='-.')
			#ax2.set_ylabel('log T [K]')
			#ax2.set_ylim(7.77, 8.2)
			#ax2.legend(loc=3)
			if conv_env_mass_range == True:
				#ax2.set_xlim(-mass_range[0]+conv_envelope_mass,conv_envelope_mass+mass_range[1])
				ax.set_xlim(xlim[0],xlim[1])
		#plt.title(title)	
				print sb_mass
			plt.ylim(-4,1)
	def final_lifetime(self):
		''' 
			Returns the final time of the evolution'
		'''
		age=[]
		mass=[]
		m=self.run_historydata
                for case in m:
			##lifetime till first TP
			age.append(np.log10(case.get("star_age")[-1]))
			mass.append(case.get("star_mass")[0])
		#plt.plot(mass,age,"*",markersize=10,linestyle="-",label=label)
		#plt.xlabel('star mass $[M_{\odot}]$',fontsize=20)
		#plt.ylabel('Logarithmic stellar lifetime',fontsize=20)
		print 'Mini','|','Age'	
		for k in range(len(age)):
			print mass[k],'|',age[k]
		return mass,age

	def lifetime(self,label=""):
		'''
			Calculate stellar lifetime till first TP 
			dependent of initial mass
		'''
                plt.rcParams.update({'font.size': 20})
                plt.rc('xtick', labelsize=20)
                plt.rc('ytick', labelsize=20)
                t0_model=self.set_find_first_TP()
                m=self.run_historydata
                i=0
		age=[]
		mass=[]
                for case in m:
			##lifetime till first TP
			age.append(np.log10(case.get("star_age")[t0_model[i]]))
			mass.append(case.get("star_mass")[0])
			i+=1
		plt.plot(mass,age,"*",markersize=10,linestyle="-",label=label)
		plt.xlabel('star mass $[M_{\odot}]$',fontsize=20)
		plt.ylabel('Logarithmic stellar lifetime',fontsize=20)	
		
		
	def final_shells(self):

		'''
			Return table of different core masses (h1,he4,final).
		'''

                m=self.run_historydata
                i=0
		print "M_ini: | M_h1 | M_he4  | M_final (withoutloss!)"
                for case in m:
                        star_mass=case.get("star_mass")
			h1bndy=case.get('h1_boundary_mass')
			he4bndy=case.get('he4_boundary_mass')
			#if case.header_attr["initial_mass"]>=5: 
			#	c12bndy=case.get('c12_boundary_mass')
			#else:
			#	c12bndy=0
			print star_mass[0],"|",h1bndy[-1],"|",he4bndy[-1],"|","|",star_mass[-1]

		
	def TPAGB_core_growth(self,fig, label="",color='k',marker_type='<',linestyle='-'):
		'''
			Creates diagram of the core growth during AGB phase.
			Uses function mult_DUP to identify first TP cycle
			and subtract from last cycle of simulation.
			(Tested with 2 cases)
		'''
                m=self.run_historydata
                i=0
                t0_model=self.set_find_first_TP()
               	core_growth=[]
		ini_m=[]
		plt.figure(fig)
		for case in m:
			mH_first_TP=case.get('h1_boundary_mass')[t0_model[i]]
			mH_last_TP=case.get('h1_boundary_mass')[-1]                
			core_growth.append(mH_last_TP-mH_first_TP)
			ini_m.append(case.header_attr["initial_mass"])
			i += 1
                plt.plot(ini_m,core_growth,marker=marker_type,color=color,markersize=10,mfc=color,linewidth=2,linestyle=linestyle,label=label)
		xlabel('Initial star  mass $[M_{\odot}]$',fontsize=18)
                ylabel('Core growth $[M_{\odot}]$',fontsize=18)
                legend()


	def set_plot_prof_ratio(self,isotopes=['c12','o16'],cycles=10*[20000],label=""):

		'''

			Plots profile of isotopes for cycles.
			Each element/cycle is related to one rundir
			if cycle is -1 than take last cycle

		'''

		fig_num=20
		label=", "+label
		m=self.run_historydata
		for i in range(len(self.run_LOGS)):		
			figure(fig_num)
			if cycles[i] == -1:
				cyc=int(m[i].get('model_number')[-1])-100
			else:
				cyc=cycles[i]
			print 'cycle ',cyc
			pr=mesa_profile(self.run_LOGS[i],cyc)
			mass=pr.get("mass")
			abu1=pr.get(isotopes[0])
			abu2=pr.get(isotopes[1])
			plot(mass,abu1/abu2,label=isotopes[0]+'/'+isotopes[1]+label)
			fig_num+=1
                        xlabel('star mass $[M_{\odot}]$',fontsize=18)
                        ylabel('Ratio Massfraction',fontsize=18)			
			legend()
	



	def set_plot_prof(self,isotopes=['he4','c12','o16','ne20','mg24','si28'],cycles=10*[20000],label=""):

		'''

			Plots profile of isotopes for cycles.
			Each element/cycle is related to one rundir
			if cycle is -1 than take last cycle

		'''

		fig_num=20
		label=", "+label
		m=self.run_historydata
		for i in range(len(self.run_LOGS)):		
			figure(fig_num)
			if cycles[i] == -1:
				cyc=int(m[i].get('model_number')[-1])-100
			else:
				cyc=cycles[i]
			print 'cycle ',cyc
			pr=mesa_profile(self.run_LOGS[i],cyc)
			for j in range(len(isotopes)):
				mass=pr.get("mass")
				abu=pr.get(isotopes[j])
				plot(mass,np.log10(abu),label=isotopes[j]+label)
			fig_num+=1
                        xlabel('star mass $[M_{\odot}]$',fontsize=18)
                        ylabel('Log(X$_i)$',fontsize=18)			
			legend()
				

	def set_plot_hrd(self,fig,symbs_1=[],linestyle=[],markevery=500,end_model=[],single_plot=True,labelmassonly=False):

		'''
			Plots HRDs
			end_model - array, control how far in models a run is plottet, if -1 till end
			symbs_1  - set symbols of runs
		'''
		m=self.run_historydata
    		i=0
		if len(symbs_1)>0:
			symbs=symbs_1
		else:
			symbs=self.symbs
		if len(linestyle)==0:
			linestyle=200*['-']
    		for case in m:
			t1_model=-1
			print end_model[i]
			if not end_model[i] == -1:
				t1_model=end_model[i]
				print self.run_label[i],t1_model
			t0_model=case.get("model_number")[0]
        		logTeff=case.get('log_Teff')[:(t1_model-t0_model)]
        		logL=case.get('log_L')[:(t1_model-t0_model)]
			
			label=self.run_label[i]
			if labelmassonly == True:
				label=label.split('Z')[0][:-2]
			print 'label',label
			if single_plot==False:
                        	figure(i+1)
                        	plot(logTeff,logL,marker=symbs[i],label=label,linestyle=linestyle[i],markevery=markevery)
				case.xlimrev()
				ax = plt.gca()
				plt.rcParams.update({'font.size': 16})
				plt.rc('xtick', labelsize=16)
				plt.rc('ytick', labelsize=16)
				legend(loc=4)
				xlabel('log Teff',fontsize=18)
				ylabel('log L',fontsize=18)
			#plt.gca().invert_xaxis()	
			figure(fig)
			plot(logTeff,logL,marker=symbs[i],label=label,linestyle=linestyle[i],markevery=markevery)
			#case.xlimrev()
			ax = plt.gca()
			plt.rcParams.update({'font.size': 16})
                	plt.rc('xtick', labelsize=16)
                	plt.rc('ytick', labelsize=16)
                	legend(loc=4)
                	xlabel('log Teff',fontsize=18)
                	ylabel('log L',fontsize=18)
			#case.xlimrev()
			#plt.gca().invert_xaxis()
			i+=1	
		figure(0)
		plt.gca().invert_xaxis()

    	def intershell_abu_agb(self,plot=True,shifts=[0],shifts_t0=[0],models_end=[-1],label=[],symbs_1=[],linestyle=[],markevery=200,savepath=''):

		'''
		Either plots or returns the intershell abundance
		for each TP.
		t0 : float
			Model marking the first TP of the TP_AGB stage
		shift  : float
		        For starting not at the first TP but at TP of model number shift
		shift_t0 : int
			For starting not at the first TP but at TP number shift_t0
		model_end: int
			Last model to be taken into accoutn. If -1: Use all models
		'''
		import mesa as ms
                m=self.run_historydata
                i=0
                if len(symbs_1)>0:
                        symbs=symbs_1
                else:
                        symbs=self.symbs
                if len(linestyle)==0:
                        linestyle=200*['-']
		if len(label)==0:
			label=self.run_label
		kk=-1
		for case in m:
			kk+=1
			logspath=self.run_LOGS[kk]

			#######!!EDIT THIS!!###################
			#mod_max=100 ### Number of log files ###
			#######################################
			if plot==True:
				plot=1 ### 0=NO; 1=YES
			else:
				plot=0
			shift=shifts[kk]
			shift_t0=shifts_t0[kk]
			model_end=models_end[kk]
			if model_end==-1:
				cutting=999999999
			else:
				cutting=model_end
			last=True

			st1=ms.star_log(logspath)

			labell=label[kk]#['M2.z1m2','M2.z2m2','M3.z1m2','M3.z2m2']
			symbol=symbs[kk]#['r:o','g--d','b-*','r-p','g--^','b-h','r:v','g-s','b--H','r-D','g--<','b-.>']

			#metallicity_label=['M3z2m2.st','M3z2m2.dfA','M3z2m2.dfB.md','M3z2m2.dfA.mdd10','M3z2m2.dfA.mdd50']

			t0_model=self.run_historydata[kk].find_first_TP()+shift_t0
						
			t0_idx=(t0_model-st1.get("model_number")[0])	
			#t0_idx=(st1.get("model_number")[0])

			first_TP_he_lum=10**(st1.get("log_LHe")[t0_idx])
			he_lum=10**(st1.get("log_LHe")[t0_idx:])
			h_lum=10**(st1.get("log_LH")[t0_idx:])
			model=st1.get("model_number")[t0_idx:]	
			h1_bndry=st1.get("h1_boundary_mass")[t0_idx:]
					#define label	
			z=st1.header_attr["initial_z"]
			mass=st1.header_attr["initial_mass"]
			print 'Case ',z,mass
			leg=str(mass)+"M$_{\odot}$ Z= "+str(z)	
			peak_lum_model=[]
			h1_mass_tp=[]
			h1_mass_min_DUP_model=[]	
					##TP identification with he lum if within 1% of first TP luminosity
			perc=0.01
			min_TP_distance=300 #model
			lum_1=[]
			model_1=[]
			h1_mass_model=[]
			TP_counter=0
			new_TP=True
			TP_interpulse=False
			interpulse_counter=0
			for i in range(len(he_lum)):
				if (h_lum[i]>he_lum[i]):
					TP_interpulse=True
					interpulse_counter+=1		
				if (h_lum[i]<he_lum[i]):
					interpulse_counter=0	
					new_TP=True
					TP_interpulse=False
					if i > 0:	
						h1_mass_1.append(h1_bndry[i])
						h1_mass_model.append(model[i])
						#print i
				if i ==0:
							#peak_lum_model=[t0_model]
							#TP_counter=1 #for first TP
					lum_1.append(first_TP_he_lum)
					model_1.append(t0_model)
					h1_mass_1=[h1_bndry[0]]
					h1_mass_model=[t0_model]
							#peak_lum_model=[t0_model]
				else:
					lum_1.append(he_lum[i])										
					model_1.append(model[i])
									
				if (new_TP == True and TP_interpulse==True and interpulse_counter> 80):#>200): #and (model[i] - t0_model	>min_TP_distance):
													#if (he_lum[i]> (perc*first_TP_he_lum)) or (i == len(he_lum)-1):		
								#if (model[i] - model_1[-1]	>min_TP_distance):			
								#calculate maximum of peak lum of certain TP
					max_value=np.array(lum_1).max()					
					max_index = lum_1.index(max_value)
								#print max_index,i
					peak_lum_model.append(model_1[max_index])
								#for DUP calc
					max_lum_idx=h1_mass_model.index(model_1[max_index])												
					min_value=np.array(h1_mass_1[max_lum_idx:]).min()					
					min_index = h1_mass_1.index(min_value)										
					h1_mass_min_DUP_model.append(h1_mass_model[min_index])					
					TP_counter+=1
					lum_1=[]
					model_1=[]
					h1_mass_1=[]		
					h1_mass_model=[]
					new_TP=False
								#TP_interpulse=False
							
					print peak_lum_model
					#print h1_mass_min_DUP_model
					#print h1_mass_tp
			modeln=[]
			interp_md=[]

			for i in range(len(peak_lum_model)):
			#	modeln.append(peak_lum_model[i])
			#	modeln.append(h1_mass_min_DUP_model[i])
				if i < (len(peak_lum_model)-1):
					if ((peak_lum_model[i+1]+h1_mass_min_DUP_model[i])/2. < cutting):
						interp_md.append((peak_lum_model[i+1]+h1_mass_min_DUP_model[i])/2.)
			#	if no_fig==True:	
			#		st1.calc_DUP_parameter(fig,modeln,leg,color,marker_type,h_core_mass)


			H1=[]
			C12=[]
			He4=[]
			O16=[] 
			zone=[]
			qqqq=[]

			### LOGS_agb2 m3z2m2 ########
			#He4.append(0.404564258324)
			#C12.append(0.383074609403)
			#O16.append(0.170402091484)
			#qqqq.append(len(qqqq)+1+shift)

			not_found=[]
			for qqq in range(0,len(interp_md)):

				C12ck=[]
				He4ck=[]
				O16ck=[] 
				N14=[]

				print '   '
				print qqq
				print '   '

				at1=ms.mesa_profile(logspath,(interp_md[qqq]),num_type='nearest_model')


				jjj=0
				#lin=1
				n=0
				
				for lin in range(len(at1.get('h1'))):
					H1.append(at1.get('h1')[lin])
					He4ck.append(at1.get('he4')[lin])
					C12ck.append(at1.get('c12')[lin])
					O16ck.append(at1.get('o16')[lin])
					N14.append(at1.get('n14')[lin])
					
					if N14[lin] < 0.0000000000001:
						if He4ck[lin] > 0.10:
							if C12ck[lin] > 0.10:
								#if  He4ck[lin] > (1.1*(O16ck[lin])):
									#print 'OK '
								if (O16ck[lin]+C12ck[lin]+He4ck[lin]) > 0.90:
									jjj=jjj+1
									if jjj== 1.:
										He4.append(at1.get('he4')[lin])
										C12.append(at1.get('c12')[lin])
										O16.append(at1.get('o16')[lin])
										zone.append(at1.get('zone')[lin])
										qqqq.append(qqq+1+shift)
										print lin
										print lin, 'He4: ', He4[0] 
										print lin, 'C12: ', C12[0] 
										print lin, 'O16: ', O16[0] 


			############################

			print 'He4: ', He4
			print '   '
			print 'C12: ', C12
			print '   '
			print 'O16: ', O16
			print 'TP num: ', qqqq
			print ''
			print '***************'
			print 'TPs: ', peak_lum_model
			#print 'TPs: ', interp_md
			print '***************'
			print ''

			#	n=n+1

			#He4.append(0.628441740854)
			#C12.append(0.225860939201)
			#O16.append(0.086891156412)

			################ !!! F O R   P L O T T I N G !!! #####################################################################
			if plot == 1:

				plt.figure(kk)
			#	plt.semilogy(qqqq, O16, 'r:', lw=4.)
				plt.semilogy(qqqq, O16, marker='o', lw=3,markersize=10, label='O16; '+str(label),linestyle='--')
			#        plt.semilogy(qqqq, O16, symbol[0], lw=3,markersize=15)

			#	plt.semilogy(qqqq, O16, 'r:o', label='O16; '+str(labell[0]), markersize=10)

			#	plt.semilogy(qqqq, C12, 'g--', lw=4.)
				plt.semilogy(qqqq, C12, marker='s', lw=3,markersize=10, label='C12; '+str(labell),linestyle='-')
			#        plt.semilogy(qqqq, C12, symbol[1], lw=3,markersize=15)

			#	plt.semilogy(qqqq, C12, 'g--d', label='C12; '+str(labell[0]),  markersize=10)

			#	plt.semilogy(qqqq, He4, 'b-', lw=4.)
				plt.semilogy(qqqq, He4, marker='D', lw=3,markersize=10, label='He4; '+str(labell),linestyle=':')
			#        plt.semilogy(qqqq, He4, symbol[2], lw=3,markersize=15)

			#	plt.semilogy(qqqq, He4, 'b-*', label='He4; '+str(labell[0]),  markersize=10)

			##			plt.axis([0.6329, 0.634, 0.00001, 1])
			#	plt.legend(loc=4,prop={'size':40})
				plt.ylabel('Mass Fraction (X)', fontsize=20)
				plt.xlabel('TP number', fontsize=20)
				ax = plt.gca()
				plt.xlim(min(qqqq),max(qqqq))
				plt.xticks(qqqq)
				box = ax.get_position()
				ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
				if len(savepath)>0:
					plt.savefig(savepath+'/M'+str(mass)+'Z'+str(z)+'.png')




	def set_plot_surf_parameter(self,fig=1,xax='model',symbs_1=[],linestyle=[],markevery=500):

		'''
			xax either 'model' or 'time'
		'''

                m=self.run_historydata
                i=0
		plt.figure(fig)
                if len(symbs_1)>0:
                        symbs=symbs_1
                else:
                        symbs=self.symbs
                if len(linestyle)==0:
                        linestyle=200*['-']
		k=-1
                for case in m:
			k+=1
			if xax == 'time':
			    xaxisarray = case.get('star_age')
			elif xax == 'model':
			    xaxisarray = case.get('model_number')
			else:
			    print 'kippenhahn_error: invalid string for x-axis selction. needs to be "time" or "model"'


			logL    = case.get('log_L')
			logTeff    = case.get('log_Teff')

			pyl.plot(xaxisarray,logL,label=self.run_label[k]+'(L)',marker=symbs_1[k],linestyle=linestyle[k],markevery=markevery)
			pyl.plot(xaxisarray,logTeff,label=self.run_label[k]+'(T)',marker=symbs_1[k],linestyle=linestyle[k],markevery=markevery)
			pyl.ylabel('log L, log T$_{eff}$')		
			pyl.legend(loc=2)


			if xax == 'time':
			    pyl.xlabel('t / yrs')
			elif xax == 'model':
			    pyl.xlabel('model number')

	def set_plot_tcrhoc(self,symbs_1=[],linestyle=[],markevery=500,end_model=[]): 
		'''
			Plots HRDs
			end_model - array, control how far in models a run is plottet, if -1 till end
			symbs_1  - set symbols of runs
		'''
		m=self.run_historydata
    		i=0
		if len(symbs_1)>0:
			symbs=symbs_1
		else:
			symbs=self.symbs
		if len(linestyle)==0:
			linestyle=200*['-']
    		for case in m:
			t1_model=-1
			if end_model[i] != -1:
				t1_model=end_model[i]
			t0_model=case.get("model_number")[0]
			rho=case.get('log_center_Rho')[:(t1_model-t0_model)]

			T=case.get('log_center_T')[:(t1_model-t0_model)]

			h1=case.get('H-1')

                        figure(i+1)
                        #plot(logTeff,logL,marker=symbs[i],label=self.run_label[i],linestyle=linestyle[i],markevery=markevery)
			pl.plot(rho,T,marker=symbs[i],label=self.run_label[i],linestyle=linestyle[i],markevery=markevery)
			case.get			
			plt.gca().invert_xaxis()
			ax = plt.gca()
			plt.rcParams.update({'font.size': 16})
			plt.rc('xtick', labelsize=16)
			plt.rc('ytick', labelsize=16)
			legend(loc=4)
		        plt.xlabel('log $\\rho_{\\rm c}$',fontsize=18)
        		plt.ylabel('log $T_{\\rm c}$',fontsize=18)

			#plt.gca().invert_xaxis()	
			figure(0)
			pl.plot(rho,T,marker=symbs[i],label=self.run_label[i],linestyle=linestyle[i],markevery=markevery)
			plt.gca().invert_xaxis()
			ax = plt.gca()
			plt.rcParams.update({'font.size': 16})
                	plt.rc('xtick', labelsize=16)
                	plt.rc('ytick', labelsize=16)
                	legend(loc=4)
		        plt.xlabel('log $\\rho_{\\rm c}$',fontsize=18)
        		plt.ylabel('log $T_{\\rm c}$',fontsize=18)
			plt.gca().invert_xaxis()
			#plt.gca().invert_xaxis()
			i+=1	

 
	def set_plot_t_mass(self):
    		
		'''
			plot star mass evolution over time
		'''
		
		m=self.run_historydata
		figure(2)
    		i=0
    		for case in m:
        		star_age=case.get('star_age')
       			h1_boundary_mass=case.get('h1_boundary_mass')
        		star_mass=case.get('star_mass')
        		plot(star_age[noffset:],h1_boundary_mass[noffset:],self.symbs[i],label=self.run_label[i])
        		plot(star_age[noffset:],star_mass[noffset:],self.symbs[i],label=self.run_label[i])
        	i += 1
    		legend(loc=6)
    		xlabel('time / yr')
    		ylabel('star mass, core mass')


	def set_plot_model_mass(self):

		'''
			def set_plot_model_mass: plot star mass evolution vs model number
		'''

    		m=self.run_historydata
		figure(3)
    		i=0
    		for case in m:
        		model_number=case.get('model_number')
        		h1_boundary_mass=case.get('h1_boundary_mass')
      			star_mass=case.get('star_mass')
			plot(model_number[noffset:],h1_boundary_mass[noffset:],self.symbs[i],label=self.run_label[i])
        		plot(model_number[noffset:],star_mass[noffset:],self.symbs[i],label='')
        		i += 1
    		legend(loc=1)
    		xlabel('model number')
    		ylabel('star mass, core mass')

	def set_plot_CO_old(self,startfirstTP=True,xtime=True,label=[],symbs_1=[],linestyle=[],markevery=500,end_model=[]):

		'''
			Plots C/O ratio vs model number or age
	
			-startfirstTP is true than start at first TP

		'''
	
		#color=['r','b','k','g']

		m=self.run_historydata    
                plt.rcParams.update({'font.size': 20})
                plt.rc('xtick', labelsize=20)
                plt.rc('ytick', labelsize=20)

                m=self.run_historydata
                if startfirstTP==True:
                        t0_model=self.set_find_first_TP()

                else:
                        t0_model=len(self.run_historydata)*[0]
		figure(4)
		if len(label)==0:
			label=self.run_label
                if len(symbs_1)>0:
                        symbs=symbs_1
                else:
                        symbs=self.symbs
                if len(linestyle)==0:
                        linestyle=200*['-']
		i=0
                for case in m:
                        t1_model=-1
                        if end_model[i] != -1:
                                t1_model=end_model[i]
			if not startfirstTP==True:
                        	t0_model=case.get("model_number")[0]
			model_number=case.get('model_number')[t0_model[i]:t1_model]
        		C=case.get('surface_c12')[t0_model[i]:t1_model]
        		O=case.get('surface_o16')[t0_model[i]:t1_model]
			CO=C*4./(O*3.)
			age=case.get("star_age")[t0_model[i]:t1_model]
                        t0_age=age[0]
                        age=age-t0_age
			if xtime==True:
				model_number=age
        		plot(model_number,CO,label=label[i],marker=symbs[i],linestyle=linestyle[i])
			legend(loc=2)
			if xtime==True:
				plt.xlabel("Star age [yr]")
			else:
				xlabel('model number')
			#plt.xlabel('Star mass $[M_{\odot}]$',fontsize=18)
			plt.ylabel('C/O Ratio', fontsize=18)
			figure(0)
			#plt.xlabel('Star mass $[M_{\odot}]$',fontsize=18)
			
			plt.ylabel('C/O Ratio', fontsize=18)
			plot(model_number,CO,label=label[i],marker=symbs[i],linestyle=linestyle[i])
			legend(loc=2)
			i += 1

	def set_plot_CO(self,fig=1321,xaxis='mass',startfirstTP=True,t0_model=[],xtime=True,label=[],symbs_1=[],linestyle=[],markevery=500,color=[],end_model=[],withoutZlabel=False,savefig=''):

		'''
			Plots C/O ratio vs model number or age
	
			-startfirstTP is true than start at first TP

		'''
	
		#color=['r','b','k','g']

		m=self.run_historydata    
                plt.rcParams.update({'font.size': 20})
                plt.rc('xtick', labelsize=20)
                plt.rc('ytick', labelsize=20)

                m=self.run_historydata
                if startfirstTP==True:
                        t0_model=self.set_find_first_TP()
                elif t0_model==-1:
                        t0_model=len(self.run_historydata)*[0]
		print 'use t0_Model,',t0_model
		if len(label)==0:
			label=self.run_label
                if len(symbs_1)>0:
                        symbs=symbs_1
                else:
                        symbs=self.symbs
                if len(linestyle)==0:
                        linestyle=200*['-']
		i=0
                for case in m:
                        t1_model=-1
                        if end_model[i] != -1:
                                t1_model=end_model[i]
                        #t0_model=case.get("model_number")[0]
                        if xaxis=='cycles':
                                x=case.get('model_number')[t0_model[i]:t1_model]
                        elif xaxis=='age':
				print 'use age'
                                x=case.get('star_age')[t0_model[i]:t1_model] - case.get('star_age')[t0_model[i]]
				print 't0: ',case.get('star_age')[t0_model[i]],x[-1],t0_model[i]
                        elif xaxis=='mass':
                                x=case.get('star_mass')[t0_model[i]:t1_model]
			model_number=case.get('model_number')[t0_model[i]:t1_model]
        		C=case.get('surface_c12')[t0_model[i]:t1_model]
        		O=case.get('surface_o16')[t0_model[i]:t1_model]
			CO=C*4./(O*3.)
			#age=case.get("star_age")[t0_model[i]:t1_model]
			figure(i+int(fig)+1)
        		plot(x,CO,label=label[i],marker=symbs[i],color=color[i],linestyle=linestyle[i],markevery=markevery)
			legend(loc=2)
			plt.xlabel('$M/M_{\odot}$',fontsize=18)
			plt.ylabel('C/O Ratio', fontsize=18)
			legend(loc=2)
                        if len(savefig)>0:
                                plt.savefig(savefig+'/co_'+label[i]+'.png')
			plt.close(i+1)
			figure(fig)
			print 'x test',x[-1]
			if xaxis=='cycle':
				plt.xlabel('Model number')
			elif xaxis=='age':
				plt.xlabel('Age/yr')
			else:
				plt.xlabel('Star mass $[M_{\odot}]$',fontsize=23)
			plt.ylabel('Surface C/O ratio', fontsize=23)
			if withoutZlabel==True:
                                plot(x,CO,label=label[i].split(',')[0],marker=symbs[i],color=color[i],linestyle=linestyle[i],markevery=markevery)
				
			else:
				plot(x,CO,label=label[i],marker=symbs[i],color=color[i],linestyle=linestyle[i],markevery=markevery)
			legend(loc=2)
			if xaxis=='mass':
				plt.gca().invert_xaxis()
			i += 1

		return t0_model



	def set_plot_vsurf(self):

		'''
			Plots the surface velocity vs model number
		'''
		
    		m=self.run_historydata
		figure(55)
   		i=0
    		for case in m:
			plt.figure(55)
			model=case.get("model_number")
			vdiv=case.get("v_div_csound_surf")
			plt.plot(model,vdiv,label=self.run_label[i])
			plt.legend()
			plt.figure(i)
			plt.plot(model,vdiv,label=self.run_label[i])
			plt.legend()
        		i += 1
    		legend(loc=2)
    		xlabel('model number')
    		ylabel('v_div_csound_surf')
    		ylim(-10.,10.)

        def set_plot_timesteps(self):

		'''
			Plots the timestep size vs model number
		'''

                m=self.run_historydata
                figure(55)
                i=0
                for case in m:
                        plt.figure(55)
                        model=case.get("model_number")
                        vdiv=case.get("log_dt")
                        plt.plot(model,vdiv,label=self.run_label[i])
                        plt.legend()
                        plt.figure(i)
                        plt.plot(model,vdiv,label=self.run_label[i])
                        plt.legend()
                        i += 1
                legend(loc=2)
                xlabel('model number')
                ylabel('Log(dt)')

	def set_plot_mdot(self,fig,xaxis="model",masslabelonly=False,marker=[],markevery=500):
		
		'''
			Plots the mass loss vs time, model or mass
			xaxis: "model","time","mass" possible
		'''

                plt.rcParams.update({'font.size': 18})
                plt.rc('xtick', labelsize=18)
                plt.rc('ytick', labelsize=18)
                

		if xaxis=="time":
                        t0_model=self.set_find_first_TP()

		elif xaxis =="model" or xaxis=="mass":
			t0_model=len(self.run_historydata)*[0]
    		m=self.run_historydata
		figure(fig)
   		i=0
    		for case in m:
			if len(marker)>0:
				marker1=marker[i]
			else:
				marker1=None	
			if xaxis=="time":

	                        t0_time=case.get('star_age')[t0_model[i]]
				print t0_time
				star_mass=case.get('star_mass')[t0_model[i]]
        	                star_age=case.get('star_age')[t0_model[i]:]  -t0_time
				mdot=case.get('log_abs_mdot')[t0_model[i]:]
				plt.plot(star_age,mdot,self.symbs[i],label=self.run_label[i],marker=marker1,markevery=markevery)
			elif xaxis=="model":
				mdot=case.get('log_abs_mdot')
				model=case.get('model_number')				
        			plt.plot(model,mdot,self.symbs[i],label=self.symbs[i],marker=marker1,markevery=markevery)
			elif xaxis=="mass":
				star_mass=case.get('star_mass')
				mdot=case.get('log_abs_mdot')
				if masslabelonly==True:
					plt.plot(star_mass,mdot,self.symbs[i],label=self.run_label[i].split('Z')[0][:-2],marker=marker1,markevery=markevery)
				else:
					plt.plot(star_mass,mdot,self.symbs[i],label=self.run_label[i],marker=marker1,markevery=markevery)
				#case.plot('star_mass','log_abs_mdot',legend=self.run_label[i],shape=self.symbs[i])
       			i += 1
    		legend(loc=2)
		if xaxis=="time":
			plt.xlabel('star age')
		elif xaxis=="model":
			plt.xlabel('model number')
		elif xaxis=="mass":
			xlabel('M/M$_{\odot}$',fontsize=18)
    		ylabel('log($|\dot{M}|$)',fontsize=18)
    		#ylim(-7,-3.5)

	def set_plot_R(self):
	
		'''
			Plots the stellar logarithmic radius vs model number
		'''

		m=self.run_historydata
	    	figure(7)
    		i=0
    		for case in m:
        		case.plot('star_age','log_R',legend=self.run_label[i],shape=self.symbs[i])
        	i += 1
    		legend(loc=2)
    		xlabel('model number')
    		ylabel('log_R')
    		if xlim_mod_min >= 0:
        		xlim(xlim_mod_min,xlim_mod_max)

	def set_plot_kip_special(self,startfirstTP=True,xtime=True,label=[],color=[]):

		'''
			Kippenhahn which plots only h and he free bndry,
			label and color can be chosen.
			if label>0 then color must be set too!
			color=["r","b","g","k"]
			

		'''
		plt.rcParams.update({'font.size': 24})
		plt.rc('xtick', labelsize=24)
		plt.rc('ytick', labelsize=24)

                m=self.run_historydata
                i=0
                if startfirstTP==True:
                        t0_model=self.set_find_first_TP()

                else:
                	t0_model=len(self.run_historydata)*[0]

		for case in m:
			h1_boundary_mass  = case.get('h1_boundary_mass')[t0_model[i]:]
			he4_boundary_mass = case.get('he4_boundary_mass')[t0_model[i]:]
			star_mass         = case.get('star_mass')[t0_model[i]:]
			mx1_bot           = case.get('mx1_bot')[t0_model[i]:]*star_mass
			model = case.get("model_number")[t0_model[i]:]
			age=case.get("star_age")[t0_model[i]:]
			t0_age=age[0]
			age=age-t0_age
			if xtime==True:
				model=age
			if len(label)>0:
				plt.plot(model,h1_boundary_mass,color=color[i],label="$^1$H bndry, "+label[i])
				plt.plot(model,he4_boundary_mass,"--",color=color[i],label="$^4$He bndry, "+ label[i])
				#pltplot(model,star_mass,color=color[i],label="Total mass")
			else:
				plt.plot(model,h1_boundary_mass,label="h1_boundary_mass")
				plt.plot(model,he4_boundary_mass,"--",label="he4_boundary_mass")
			#plt.plot(model,mx1_bot,label="convective boundary")
			#title(self.run_label[i])
			i += 1
		if xtime==True:
			plt.xlabel('stellar age',size=28)
			if startfirstTP==True:
				plt.xlabel('t - t$_0$ $\mathrm{[yr]}$',size=28)
		else:	
			plt.xlabel('model number',size=28)
		plt.ylabel("M/M$_{\odot}$",size=28)		
		plt.legend()	
	def set_plot_kipp(self,fig=0,xaxis='model',startfirstTP=True,savepath=''):
		'''
			fig: id of first figure; following ones are fig+1,fig+2..
		'''
    		m=self.run_historydata
		i=0
                if startfirstTP==True:
                        t0_model=self.set_find_first_TP()

                else:
                        t0_model=len(self.run_historydata)*[0]
		
		if xaxis=='time':
    			for case in m:
        			case.kippenhahn(int(fig)+i,'time',t0_model=t0_model[i],c12_bm=False)
        			title(self.run_label[i])
        			i += 1
				mass=case.header_attr['initial_mass']
				z=case.header_attr['initial_z']
				if (len(savepath)>0) and (not '.png' in savepath):
                                	outfile1=savepath+'/M'+str(mass)+'_Z'+str(z)+'_kipp.png'
					plt.savefig(outfile1)
		elif xaxis=='logtimerev':
    			for case in m:
        			case.kippenhahn(int(fig)+i,'logtimerev',t0_model=0,c12_bm=False)
        			title(self.run_label[i])
        			i += 1
				mass=case.header_attr['initial_mass']
				z=case.header_attr['initial_z']
				if (len(savepath)>0) and (not '.png' in savepath):
                                	outfile1=savepath+'/M'+str(mass)+'_Z'+str(z)+'_kipp.png'
					plt.gca().invert_yaxis()
					plt.savefig(outfile1)


		else:
                	for case in m:
                                case.kippenhahn(int(fig)+i,'model',t0_model=t0_model[i],c12_bm=False)
                                title(self.run_label[i])
                                i += 1
				mass=case.header_attr['initial_mass']
				z=case.header_attr['initial_z']
				if (len(savepath)>0) and (not '.png' in savepath):
                                	outfile1=savepath+'/M'+str(mass)+'_Z'+str(z)+'_kipp.png'
					plt.savefig(outfile1)
		if len(savepath)>0:
			plt.close('all')


	def set_plot_kipp_cont(self,fig=32,mix_burn_zones=20,modstart=0,modstop=-1,ixaxis='model_number',xres=10000,yres=10000,xlims=[0,0],ylims=[0,0],plot_radius=False,boundaries=True,rad_lines=False,engenPlus=True,engenMinus=False,savepath='',filetype='png'):
		'''
				Plots kipp_cont diagrams.			
				
			        xres,yres         x and y resolution
        			xlims, ylims      x and y plot limits
        			ixaxis            'age', 'model_number' or 'log_time_left'
        			outfile           name of output file including extension, whic
		'''

                m=self.run_historydata
                i=fig
		for case in m:
			if ylims[1] ==0:
				ylims1=[ylims[0],case.header_attr['initial_mass']]
			else:
				ylims1=ylims
			print ylims1
			mass=case.header_attr['initial_mass']
			z=case.header_attr['initial_z']
			if (len(savepath)>0) and (not '.png' in savepath):
				outfile1=savepath+'/M'+str(mass)+'_Z'+str(z)+'_kippcont.png'
			else:
				outfile1=''
			case.kip_cont(ifig=i,modstart=modstart,modstop=modstop,xlims=xlims,ylims=ylims1,xres=xres,yres=yres,ixaxis=ixaxis,mix_zones=mix_burn_zones,burn_zones=mix_burn_zones,plot_radius=plot_radius,engenPlus=engenPlus,engenMinus=engenMinus,landscape_plot=False,rad_lines=False,profiles=[],showfig=True,outlines=True,boundaries=boundaries,outfile=outfile1,rasterise=False)
			i +=1



        def set_plot_hfree_core(self,fig=1,xaxis='model',marker=['o'],linestyle=['-']):

    		m1=self.run_historydata
		i=-1
		plt.figure(fig)
    		for m in m1:
			i+=1
			h1bdy=m.get('h1_boundary_mass')
			model=m.get('model_number')
			plt.plot(model,h1bdy,label=self.run_label[i])

		plt.xlabel('model_number')
		plt.ylabel('H-free core')


	def set_plot_kipp_CO(self, startfirstTP=False,savefig=''):

		'''

		plots kippenhahn diagrams with c/o ratio

		'''

		if startfirstTP==True:
			t0_model=self.set_find_first_TP()
		else:
			t0_model=len(self.run_historydata)*[0]
		m=self.run_historydata
    		i=0
    		for case in m:
        		case.kippenhahn_CO(i,'model',t0_model=t0_model[i])
        		title(self.run_label[i])
			if len(savefig)>0:
				plt.savefig(savefig+self.run_label[i]+'.png')
				plt.close()
        		i += 1
				

	def set_plot_surface_abu(self,num_frame=1,species=['c13','c12'],ratio=False,age_years=False, xax='model', t0_model=0,log=False,label=['Mesa','Mesa'],marker=['o'],linestyle=['-'],markevery=10000):

		'''
			Plot surface abundances of t_surfabu in nugridse
		'''

    		m1=self.run_historydata
		i=-1
    		for m in m1:
			i+=1
			if num_frame >= 0:
			    plt.figure(num_frame)
			if xax == 'time':
			    xaxisarray = m.get('star_age')[t0_model:]
			elif xax == 'model':
			    xaxisarray = m.get('model_number')[t0_model:]
			elif xax == 'logrevtime':
			    xaxisarray = m.get('star_age')
			    xaxisarray=np.log10(max(xaxisarray[t0_model:])+t_eps-xaxisarray[t0_model:])
                            if age_years==True:
				xaxisarray=xaxisarray*(365.*24.*3600.)
			elif xax == 'mass':
			    xaxisarray = m.get('star_mass')[t0_model:]
			else:
			    print 't-surfabu error: invalid string for x-axis selction.'+ \
				  ' needs to be "time" or "model"'

			star_mass         = m.get('star_mass')
			abu=[]

			t0_mod=xaxisarray[t0_model]
			
			for k in range(len(species)):
				abu.append(m.get('surface_'+species[k])[t0_model:])
			#surface_c12       = m.get('surface_c12')
			#surface_o16       = m.get('surface_o16')

			target_n14 = -3.5

			#symbs=['k:','-','--','-.','b:','-','--','k-.',':','-','--','-.']
			if ratio==False:
				if log==True:
					for k in range(len(species)):
						plt.plot(xaxisarray,np.log10(abu[k]),marker=marker[k],label=species[k]+label[0],linestyle=linestyle[k],markevery=markevery)
					plt.ylabel('mass fraction $\log X$')
				else:
					for k in range(len(species)):
						plt.plot(xaxisarray,abu[k],marker=marker[k],label=species[k]+label[0],linestyle=linestyle[k],markevery=markevery)
	
					plt.ylabel('mass fraction X$')

			else:
				pyl.plot(xaxisarray,np.array(abu[1])/np.array(abu[0]),\
                                                     marker=marker[i],linestyle=linestyle[i],markevery=markevery,label=label[i])	

				plt.ylabel(species[0]+'/'+species[1])
			pyl.legend(loc=2)

			if xax == 'time':
			    pyl.xlabel('t / yrs')
			elif xax == 'model':
			    pyl.xlabel('model number')
			elif xax == 'logrevtime':
			    pyl.xlabel('$\\log t-tfinal$')
			#if xax == 'logrevtime':
			 #   self._xlimrev()
			if xax=='mass':
			    plt.xlabel('M/Msun')


	def set_plot_lumi(self,fig=13,xaxis='model',totlum=False,startfirstTP=False,marker=[],linestyle=[],markevery=500,lastmodel=[-1]):              
		

		'''
			Plot luminosity vs. t; fig integer 
		'''

		m=self.run_historydata
                i=0
		#plt.figure(fig)
                if startfirstTP==True:
                        t0_model=self.set_find_first_TP()
                else:
                        t0_model=len(self.run_historydata)*[0]
                for case in m:
			#plt.figure(int(fig)+i)
			if lastmodel[i]==-1:
				t1_model=-1
			else:
			   	t1_model=lastmodel[i]
			#figure(12345)
			if totlum==False:
				if xaxis=='model':
					case.t_lumi(int(fig)+i,'model')
				elif xaxis=='time':
					case.t_lumi(int(fig)+i,'time')
				elif xaxis=='mass':
					plt.figure(int(fig)+i)
					logLH   = case.get('log_LH')[t0_model[i]:]
					logLHe  = case.get('log_LHe')[t0_model[i]:]
					mass=case.get('star_mass')[t0_model[i]:]
					plt.plot(mass,logLH,label='$L_H$')
					plt.plot(mass,logLHe,label='$L_{He}$')
					plt.xlabel('$M/M_{\odot}$')
			else:
				plt.figure(int(fig)+i)
				logL=np.log10(np.array(case.get('luminosity')))[t0_model[i]:t1_model]
				if xaxis=='model':
					x=case.get('model_number')[t0_model[i]:t1_model]
					plt.xlabel('Model number')
				elif xaxis=='time':
					x=case.get('star_age')[t0_model[i]:t1_model]
					plt.xlabel('age/years')
				elif xaxis=='mass':
					x=case.get('star_mass')[t0_model[i]:t1_model]
					plt.xlabel('$M/M_{\odot}$')
				plt.plot(x,logL,label=self.run_label[i].split(',')[0],marker=marker[i],linestyle=linestyle[i],markevery=markevery)
                        plt.ylabel('log $L/L_{\odot}$')
                        plt.legend(loc=2)
                        legend(loc=4)
                        i += 1

	

	def set_thing(self,things=['log_LH','log_LHe']):
		'''
			Plots of what the user specifies
			Plots list of things
		'''
    		m=self.run_historydata
		i=70
   		for case in m:
        		figure(i)
			j = 0
			for thing in things:
            			case.plot('model_number',thing,legend=thing,shape=self.symbs[j])
            			j += 1
        		title(self.run_label[i-70])
        		i += 1


	def coregrowth_time(self):
		'''
			Core grwoth + time in each interpulse phase, for Jeff
		'''		
	
    		m=self.run_historydata
		deltat_all=[]
		h1_growth_all=[]
		model_start_all=[]
		model_stop_all=[]
		corelum_relation_mc_all=[]
		corelum_relation_L_all=[]
		age_start_all=[]
   		for case in m:
			deltat=[]
			h1_growth=[]
			#TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,c13_pocket_min,c13_pocket_max= case.TPAGB_properties()
			peak_lum_model,h1_mass_min_DUP_model=case.find_TP_attributes( 3, t0_model=case.find_first_TP(), color='r', marker_type='o')
			h1_bdy=case.get('h1_boundary_mass')
			age=case.get('star_age')
			LH=case.get('log_LH')
			L=case.get('log_L')
			LHe=case.get('log_LHe')	
			model_number=case.get('model_number')
			interpulse=False
			corelum_relation_mc=[]
			corelum_relation_L=[]
			model_start=[]
			model_stop=[]
			age_start=[]
			for k in range(len(model_number)):
				mod=model_number[k]
				if mod<peak_lum_model[0]:
					continue
				if (LH[k]>LHe[k]) and not interpulse==True:
					interpulse=True
					idx_start=k #age[k+100]	
					#h1_start=h1_bdy[k+100]
					print 'start ',	mod,k
				if (LH[k]<LHe[k]) and interpulse==True:
					interpulse=False
					idx_end=k
					idx_diff=idx_end-idx_start
					phase=idx_diff*0.1
					idx_end=int(idx_end-phase)
					idx_start=int(idx_start+phase)
					if (age[idx_end] - age[idx_start])<0  or (h1_bdy[idx_end] -h1_bdy[idx_start])<0:
						continue
					deltat.append(age[idx_end] - age[idx_start])
					h1_growth.append(h1_bdy[idx_end] -h1_bdy[idx_start])
					print 'start ', idx_start
					print 'end ',idx_end
					model_start.append(idx_start)
					model_stop.append(idx_end)
					age_start.append(age[idx_start])
					idx_max=list(LH).index(max(LH[idx_start:idx_end]))
					corelum_relation_mc.append(h1_bdy[idx_max])
					corelum_relation_L.append( L[idx_max])

			'''
			for k in range(len(TPstart)-1):		
				if DUPmods[k] == -1:
					print 'Skip ', k
					

					continue
				h1_growth.append(h1_bdy[TPstart[k+1]]-h1_bdy[DUPmods[k]])
				deltat.append(age[TPstart[k+1]]-age[DUPmods[k]])
			'''
			deltat_all.append(deltat)
			h1_growth_all.append(h1_growth)	
			age_start_all.append(age_start)
			model_start_all.append(model_start)
			model_stop_all.append(model_stop)
			corelum_relation_mc_all.append(corelum_relation_mc)
			corelum_relation_L_all.append(corelum_relation_L)
		return h1_growth_all,deltat_all,model_start_all,model_stop_all,age_start_all,corelum_relation_mc_all,corelum_relation_L_all

	def set_find_first_TP_age(self,dirs=[]):

		'''
			Time at first TP
		'''
	        m=self.run_historydata
		firstTPs=[]
		ages=[]
		for case in m:
			firstTPs.append(case.find_first_TP())
			age=case.get('star_age')[int(case.find_first_TP())]
			ages.append(age)
		return ages

	def set_find_first_TP(self,dirs=[]):

		'''
			Find first TP of all runs and return array
		'''


                historydata=[]
                if (len(dirs)) == 0:
                        dirs=self.run_LOGS
                        historydata=self.run_historydata
                else:
                        for i in range(len(dirs)):
                                historydata.append(history_data(dirs[i]+"/LOGS"))
		t0_models=[]
		for j in range(len(dirs)):
			h1_boundary_mass  = historydata[j].get('h1_boundary_mass')
                	he4_boundary_mass = historydata[j].get('he4_boundary_mass')
                	star_mass         = historydata[j].get('star_mass')
                	#mx1_bot           = historyda.get('mx1_bot')*star_mass
                	#mx1_top           = historydata.get('mx1_top')*star_mass
                	mx2_bot           = historydata[j].get('mx2_bot')*star_mass
                	#mx2_top           = historydata.get('mx2_top')*star_mass
			he_lumi		  = historydata[j].get('log_LHe')
			h_lumi 		  = historydata[j].get('log_LH')
			#model_number            = historydata[j].get('model_number')
			lum_array=[]
			activate=False
			models=[]
			for i in range(len(h1_boundary_mass)):
				if (h1_boundary_mass[i]-he4_boundary_mass[i] <0.1) and (he4_boundary_mass[i]>0.2):
					if (mx2_bot[i]>he4_boundary_mass[i]) and (he_lumi[i]>h_lumi[i]):
						activate=True
						lum_array.append(he_lumi[i])
						models.append(i)
					if (activate == True) and (he_lumi[i]<h_lumi[i]):
						break	
			t0_models.append(models[np.argmax(lum_array)]  )
		return t0_models		

	def write_gce_input_total_lum(self,file_name="isotopic_table.txt",sparsity=100,final_models=[]):

		'''
			Calculates average total luminosity.

			sparsity: for averaging: split in intervals  of total time divided by sparsity
			final_models: list of final models until which evolution is considered
		'''

		#mesa model cases here
                m=self.run_historydata
		mz=[]
		lum=[]
		p=-1
                for case in m:
			p+=1
                        mass=case.header_attr['initial_mass']
                        metallicity=case.header_attr['initial_z']
                        mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
			if len(final_models)>0:
				models_all=case.get('model_number')
				idx_f=list(models_all).index(final_models[p])
                                L1=10**case.get('log_L')[:idx_f] #Lsu
                                T1=10**case.get('log_Teff')[:idx_f]  #K
                                age_all=case.get('star_age')[:idx_f]
			else:
				L1=10**case.get('log_L') #Lsu
				T1=10**case.get('log_Teff')  #K
				age_all=case.get('star_age')

			times=[]  #case.get('star_age')

			#mean luminosity taken during H burning phase, assume similiar
			final_age=age_all[-1]
			time_step_size=final_age/sparsity
			models=[] #case.get('model_number')
			for tt in arange(time_step_size,final_age+final_age,time_step_size):
				#find closest model
				idx=min(range(len(age_all)), key=lambda i: abs(age_all[i]-tt))
				models.append(idx)
				times.append(age_all[idx])
	
                        mz.append(mz1)

			kk=-1
			lum1=[]
			for i in models:
				kk+=1
				L_tot=L1[i]
				T=T1[i]
				if i==models[0]:
					delta_time= times[kk]+(times[kk+1]-times[kk])/2.
				elif i==models[-1]:
					delta_time= (times[kk]-times[kk-1])/2.+(final_age-times[kk])
				else:
					delta_time= (times[kk]-times[kk-1])/2.+(times[kk+1]-times[kk])/2.
				#weighted values
				lum1.append(L_tot* delta_time )
	
			#create weighted_values
			lum.append( sum(lum1)/final_age )

		#convert luminosity from Lsun to erg
		lsol = 3.8418e33 # from MESA ! solar luminosity (erg s^-1)
		lum = np.array(lum)*lsol	
		##write in file	
		f1=open(file_name,'r')
		lines=f1.readlines()
		f1.close()
		i=-1
		line1=''
		while (True):
			i+=1
			if i>len(lines)-1:
				break
			line=lines[i]	
			line1+=lines[i]
			for k in range(len(mz)):
				if mz[k] in lines[i]:
					line1+=('H Average Luminosity: '+'{:.3E}'.format(lum[k])+'\n')
			
		f1=open(file_name,'w')
                f1.write(line1)
                f1.close()

	def write_gce_input_spectra_input(self,file_name="param_for_spectra.txt",sparsity=100,genec_input=False,final_models=[]):


		'''
			Derives mean log g and Teff and produces table. Useful to adopt for model spectra.

			Table format:
				# run name & Teff & log g  & log R
			sparsity: for averaging: split in intervals  of total time divided by sparsity

		'''

		#mesa model cases here
                m=self.run_historydata
		mz=[]
		Teff_means=[]
		logg_means=[]
		logR_means=[]
		p=-1
                for case in m:
			p +=1
                        mass=case.header_attr['initial_mass']
                        metallicity=case.header_attr['initial_z']
                        mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
                        if len(final_models)>0:
                                models_all=case.get('model_number')
				if final_models[p]<0:
					final_models[p]=models_all[-1]
				print mz1,final_models[p]
                                idx_f=list(models_all).index(final_models[p])
				Teff_all=10**case.get('log_Teff')[:idx_f]  #K
				logg_all=10**case.get('log_g')[:idx_f]
				age_all=case.get('star_age')[:idx_f]
				logR_all=case.get('log_R')[:idx_f]
                        else:

				Teff_all=10**case.get('log_Teff')  #K
				logg_all=case.get('log_g')
				age_all=case.get('star_age')
				logR_all=case.get('log_R')

			models=[]   #case.get('model_number')
			times=[]  #case.get('star_age')
			#mean luminosity taken during H burning phase, assume similiar
			final_age=age_all[-1]
			time_step_size=final_age/sparsity
			models=[] #case.get('model_number')
			for tt in arange(time_step_size,final_age+final_age,time_step_size):
				#find closest model
				idx=min(range(len(age_all)), key=lambda i: abs(age_all[i]-tt))
				models.append(idx)
				times.append(age_all[idx])
	
                        mz.append(mz1)

			kk=-1
			logg11=[]
			Teff11=[]
			logR11=[]
			for i in models:
				kk+=1
				logg=logg_all[i]
				Teff=Teff_all[i]
				logR=logR_all[i]
				if i==models[0]:
					delta_time= times[kk]+(times[kk+1]-times[kk])/2.
				elif i==models[-1]:
					delta_time= (times[kk]-times[kk-1])/2.+(final_age-times[kk])
				else:
					delta_time= (times[kk]-times[kk-1])/2.+(times[kk+1]-times[kk])/2.
				#weighted values
				logg11.append(logg *delta_time)
				Teff11.append(Teff * delta_time)	
				logR11.append(logR *  delta_time)
			#create weighted_values
	                Teff_means.append( sum(Teff11)/final_age )
        	        logg_means.append( sum(logg11)/final_age )
			logR_means.append(sum(logR11)/final_age )

		print 'Teff_means',Teff_means
		print 'logg_means',logg_means
		print 'logR_means: ',logR_means
		##write in file
		import os.path
		line1=''
		if os.path.isfile(file_name):
			f1=open(file_name,'r')
			lines=f1.readlines()
			f1.close()
		else:
			lines=[]
		i=-1
		while (True):
			i+=1
			if i>len(lines)-1:
				for k in range(len(Teff_means)):
					line1+= (mz[k]+' & '+str(Teff_means[k]) + ' & '+str(logg_means[k])+' & '+str(logR_means[k])+'\n')
				break
			line1+=lines[i]
			
		f1=open(file_name,'w')
                f1.write(line1)
                f1.close()

	def write_gce_input_lum_bands(self,file_name="isotopic_table.txt",sparsity=100,final_models=[],genec_input=False):

		'''
			Calculates and writes average (!) luminosity in three bands:
			Lyman-Werner (11.18-13.6 eV) star formation;
			hydrogen-ionizing band (13.6 - 24.6 eV) 
			helium-ionizing/high energy band (> 24.6 eV)

			sparsity: divide total lifetime in sparsity equal time intervals for averaging
	
		'''
		import scipy.integrate as int
		
		lum_bands=[]
		for i in range(3):
			lum_bands.append([])
		mz=[]

		#for case of Genec models, need hdf5 file for extraction
		if genec_input: #self.run_historydata[0].header_attr['initial_z'] ==0.01 or self.run_historydata[0].header_attr['initial_z'] ==0.02:
			rundir=self.run_dir
			runs=[]
			slist = os.listdir(rundir)
			for element in slist:
				if len(self.multi_dir)>0:
					if not element in self.multi_dir:
						continue
				run_path=rundir+"/"+element
				run_path_1=glob.glob(run_path+"/*.h5")
				if len(run_path_1)>3:
					runs.append(run_path)	
			print 'GENEC models found:',runs	

			for run in runs:
				sefiles=se(run)
				mass=round(sefiles.get(0,'mass')[0])
				metallicity=self.run_historydata[0].header_attr['initial_z']
				mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
				print 'add GENEC',mz1
				models=[]
				times=[]
				a=sefiles.get(0,'age')
				for mod in range(0,np.int(sefiles.se.cycles[-1]),sparsity):
					#print type(mod)
					models.append(np.int(mod))
					#print mod
					times.append(sefiles.get(np.int(mod),'age'))
				mz.append(mz1)
				lum_bands1=[[],[],[]]
				test_bands=[[],[],[]]
				#print times
				#print 'go through models'
				kk=0
			
				c_light=2.99792458e10 #speed of light in vacuum (cm s^-1), as in MESA
				pi = 3.1415926535897932384626433832795029e0 #from MESA
				planck_h  = 6.62606896e-27 # Planck's constant (erg s), as in MESA
				ev2erg = 1.602176487e-12 # electron volt (erg), as in MESA
				rsol = 6.9598e10 # solar radius (cm) as in MESA
				k_B= 1.3806504e-16 #erg/K 
				for i in models[1:]:
					kk+=1
					L_tot=10**(sefiles.get(i,'logL'))
					T=10**(sefiles.get(i,'logTeff')) 
					R=10**(sefiles.get(i,'logR'))*rsol
					#T=T1[i] #K
					#L_tot=L1[i] #Lsun

					x1=[11.18,13.6,24.6] #eV
					x2=[13.6,24.6,inf] #eV
	
					for t in range(len(x1)):
				

						'''	
						u1= x1[t]/ (T * 8.6173324e-5) # calculation nom and denom in eV #no unit
						u2= x2[t]/ (T * 8.6173324e-5) # calculation nom and denom in eV #no unit

						int_u=int.quad(lambda x: (x**3./(np.exp(x) -1)) ,u1,u2)[0]
						L_delta_x = ( int_u/(np.pi**4/15.)) * L_tot * 3.8418e33 * 31536000 #Lsun/>erg/s/>year
						delta_time=times[kk]-times[kk-1]
						#weighted values
						'''
						
						lambda1= planck_h * c_light/(x2[t]*ev2erg) #cm
						lambda2= planck_h * c_light/(x1[t]*energy_band[0]*ev2erg) #cm
						#1e7 for erg, L in erg
						L=1e7 * 4*pi*R**2*4*pi*2*planck_h*c_light*int.quad(lambda x: (1./x**5 /(np.exp(planck_h*c_light/(x*T*k_B)) -1.)) ,lambda1,lambda2)[0]

						lum_bands1[t].append(L_delta_x * delta_time )
						test_bands[t].append(L_delta_x)	

				#create weighted_values
				for k in range(len(lum_bands1)):
					#plt.figure(0)
					lum_bands[k].append( sum(lum_bands1[k])/times[-1] )
					#plt.plot(times[1:],test_bands[k],label='band '+str(k))
					#plt.yscale('log')
					print 'mean '+str(k)+' '+str(sum(lum_bands1[k])/times[-1])
		######
		#mesa model cases here
                m=self.run_historydata
		p=-1
                for case in m:
			p+=1
                        mass=case.header_attr['initial_mass']
                        metallicity=case.header_attr['initial_z']
                        mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
                        if len(final_models)>0:
                                models_all=case.get('model_number')
                                idx_f=list(models_all).index(final_models[p])
				L1=10**case.get('log_L')[:idx_f] #Lsu
				T1=10**case.get('log_Teff')[:idx_f]  #K
				age_all=case.get('star_age')[:idx_f]
                        else:
				L1=10**case.get('log_L') #Lsu
				T1=10**case.get('log_Teff')  #K
				age_all=case.get('star_age')

			models=[]   #case.get('model_number')
			times=[]  #case.get('star_age')

			#mean luminosity taken during H burning phase, assume similiar

			final_age=age_all[-1]
			time_step_size=final_age/sparsity
			models=[] #case.get('model_number')
			for tt in arange(time_step_size,final_age+final_age,time_step_size):
				#find closest model
				idx=min(range(len(age_all)), key=lambda i: abs(age_all[i]-tt))
				#print type(mod)
				#models.append(int(models[idx]))
				models.append(idx)
				#print mod
				times.append(age_all[idx])
	
                        mz.append(mz1)
			lum_bands1=[[],[],[]]
			test_bands=[[],[],[]]

			kk=-1
			for i in models:
				kk+=1
				L_tot=L1[i]
				T=T1[i]
	
				x1=[11.18,13.6,24.6] #eV
				x2=[13.6,24.6,inf] #eV
	
				for t in range(len(x1)):
				
					u1= x1[t]/ (T * 8.6173324e-5) # calculation nom and denom in eV #no unit
					u2= x2[t]/ (T * 8.6173324e-5) # calculation nom and denom in eV #no unit

					int_u=int.quad(lambda x: (x**3./(np.exp(x) -1)) ,u1,u2)[0]
					L_delta_x = ( int_u/(np.pi**4/15.)) * L_tot * 3.8418e33 * 31536000 #Lsun/>erg/s/>year

					if i==models[0]:
						delta_time= times[kk]+(times[kk+1]-times[kk])/2.
					elif i==models[-1]:
						delta_time= (times[kk]-times[kk-1])/2.+(final_age-times[kk])
					else:
						delta_time= (times[kk]-times[kk-1])/2.+(times[kk+1]-times[kk])/2.
			
					#weighted values
					lum_bands1[t].append(L_delta_x * delta_time )
					test_bands[t].append(L_delta_x)
	
			#create weighted_values
			for k in range(len(lum_bands1)):
				lum_bands[k].append( sum(lum_bands1[k])/final_age )
				#plt.figure(1)
				#plt.plot(times[1:],test_bands[k],label='band '+str(k),marker='o')
				#plt.yscale('log')
				print 'mean '+str(k)+' '+str(sum(lum_bands1[k])/final_age)
		print mz
		print lum_bands
		##write in file	
		f1=open(file_name,'r')
		lines=f1.readlines()
		f1.close()
		i=-1
		line1=''
		while (True):
			i+=1
			if i>len(lines)-1:
				break
			line=lines[i]	
			line1+=lines[i]
			for k in range(len(mz)):
				if mz[k] in lines[i]:
					lum_bands_title=['H Average Luminosity, Lyman-Werner band','H Average  Luminosity, Hydrogen-ionizing band','H Average  Luminosity, High-energy band']
					for w in range(len(lum_bands)):
						line1+=(lum_bands_title[w]+': '+'{:.3E}'.format(lum_bands[w][k])+'\n')
					break			
			
		f1=open(file_name,'w')
                f1.write(line1)
                f1.close()



	def write_gce_input_lum_bands1(self,file_name="isotopic_table.txt",timesteps=500):

		'''
			Calculates and writes average (!) luminosity in three bands:
			Lyman-Werner (11.18-13.6 eV)
			hydrogen-ionizing band (13.6 - 24.6 eV)
			helium-ionizing/high energy band (> 24.6 eV)
			
		'''
		import scipy.integrate as int

		
		lum_bands=[]
		for i in range(3):
			lum_bands.append([])
		mz=[]
                m=self.run_historydata
                for case in m:
                        mass=case.header_attr['initial_mass']
                        metallicity=case.header_attr['initial_z']
                        mz1='(M='+str(round(mass,2))+',Z='+str(metallicity)+')'
			#average luminosity, 
			L1=10**case.get('log_L') #Lsu
			T1=10**case.get('log_Teff')  #K
			models=case.get('model_number')
			times=case.get('star_age')
			#h1=case.get('center_h1')
			#mean luminosity taken during H burning phase, assume similiar 
                        #for i in range(len(models)):
			#	if h1[i]<0.3:
                         #         def create_label(self,historydata,extra_label):
                        T=T1[i] #K
                        L_tot=L1[i] #Lsun

                        x1=[11.18,13.6,24.6] #eV
                        x2=[13.6,24.6,inf] #eV

                        for t in range(len(x1)):

                                u1= x1[t]/ (T * 8.6173324e-5) # calculation nom and denom in eV #no unit
                                u2= x2[t]/ (T * 8.6173324e-5) # calculation nom and denom in eV #no unit

                                int_u=int.quad(lambda x: (x**3./(np.exp(x) -1)) ,u1,u2)[0]
                                L_delta_x = ( int_u/(np.pi**4/15.)) * L_tot * 3.8418e33 * 31536000 #Lsun/>erg/s/>year
                                lum_bands[t].append(L_delta_x)

                print mz
                print lum_bands
                ##write in file 
                f1=open(file_name,'r')
                lines=f1.readlines()
                f1.close()
                i=-1
                line1=''
                while (True):
                        i+=1
                        if i>len(lines)-1:
                                break
                        line=lines[i]
                        line1+=lines[i]
                        for k in range(len(mz)):
                                if mz[k] in lines[i]:
                                        lum_bands_title=['H Lyman-Werner band','H Hydrogen-ionizing band','H High-energy band']
                                        for w in range(len(lum_bands)):
                                                line1+=(lum_bands_title[w]+': '+'{:.3E}'.format(lum_bands[w][k])+'\n')
                                        break

                f1=open(file_name,'w')
                f1.write(line1)
                f1.close()



        #####################Internal use

        def create_label(self,historydata,extra_label):
		mass=historydata.header_attr['initial_mass']		
		z=historydata.header_attr['initial_z']
		if extra_label==" ":
			return str(mass)+"$M_{\odot}$, Z="+str(z)
		else:
                        return str(mass)+"$M_{\odot}$, Z="+str(z)+" , "+extra_label


class model_set(se,history_data):

        '''
	    Class using mesa and mppnp runs
	'''
		
	def __init__(self,mesarundir='.',mppnprundir='.',mesa_dir=[],mppnp_dir=[],extra_label=99*[' ']):
	
		mesainit=mesa_set(rundir=mesarundir,multi_dir=mesa_dir,extra_label=99*[' '],newstarlog=False)
		mppnpini=mppnp_set(rundir=mppnprundir,multi_dir=mppnp_dir,extra_label=99*[' '])	
		self.run_historydata=mesainit.run_historydata		
                self.run_LOGS=mesainit.run_LOGS
		self.runs_H5_surf=mppnpini.runs_H5_surf
                self.runs_H5_out=mppnpini.runs_H5_out
                self.runs_H5_restart=mppnpini.runs_H5_restart


	def plot_surface_abu(self,fig=1,xaxis='model',species=[''],ratio=False,log=True,sparsity=200,markevery=10,savefig='',maxmod=[],closeall=False):
	
		'''
			Combine kip and surface abundance of isotopes
		'''


		for hh in range(len(self.run_historydata)):
                                historydata=self.run_historydata[hh]
				model=historydata.get('model_number')
                                sefiles=se(self.runs_H5_surf[hh])
                                mini=float(sefiles.se.get('mini'))
                                zini=float(sefiles.se.get('zini'))
                                label=str(mini)+'$M_{\odot}$, Z='+str(zini)
				cycs=[]
				for k in range(len(sefiles.se.cycles)):
					cycs.append(int(sefiles.se.cycles[k]))
				cycs=cycs[::sparsity]
				x1=[[],[]]
				for k in range(len(species)):
					if ratio==True:
						x1[k]=sefiles.get(cycs,species[k])
						continue
					historydata.kippenhahn(fig,xaxis)	
					specie=sefiles.get(cycs,species[k])
					fig+=1
					ax=plt.gca()
					ax2=ax.twinx()	
					ax2.plot(cycs,specie,marker='o',linestyle='--',label=species[k],markevery=markevery)
					ax2.legend(loc=4)
                                        plt.xlabel('Model number')
                                        ax2.set_ylabel('X('+species[k]+')')
					if log==True:
						ax2.set_yscale('log')
					min1=min(specie)
					max1=max(specie)
					powmin1=float('1E'+'{:.0E}'.format(min1).split('E')[1])
					powmax1=float('1E'+'{:.0E}'.format(max1).split('E')[1])
					if not powmin1==powmax1:
						ax2.set_ylim(min1,max1)
					else:
						if  abs(powmin1-min1)>abs(powmax1*10-max1):
							ax2.set_ylim(min1,powmax1*10)
						else:	
							ax2.set_ylim(powmin1,max1)
							
                                        if (len(maxmod)>0) and (maxmod[hh])>0:
                                                plt.xlim(0,maxmod[hh])
					if len(savefig)>0:
						plt.savefig(savefig+'/'+species[k]+'_'+'M'+str(mini)+'Z'+str(zini)+'.png')
						if closeall==True:
							plt.close()
				if ratio==True:
					quot=np.array(x1[0])/np.array(x1[1])
					historydata.kippenhahn(fig,xaxis)
                                        ax=plt.gca()
                                        ax2=ax.twinx()
					label_ratio='X('+species[0]+')/X('+species[1]+')'
                                        ax2.plot(cycs,quot,marker='o',linestyle='--',label=label_ratio,markevery=markevery)
                                        ax2.legend(loc=4)
                                        plt.xlabel('Model number')
                                        ax2.set_ylabel(label_ratio)
                                        if log==True:
                                                ax2.set_yscale('log')
                                        ax2.set_ylim(min(quot),max(quot))
                                        if (len(maxmod)>0) and (maxmod[hh])>0:
                                                plt.xlim(0,maxmod[hh])
                                        if len(savefig)>0:
                                                plt.savefig(savefig+'/'+species[0]+'_'+species[1]+'_'+'M'+str(mini)+'Z'+str(zini)+'.png')
                                                if closeall==True:
                                                        plt.close()


					
	  #      set.set_plot_surface_abu(self,fig=2,species=['Sr-88','Ba-136'],xaxis='cycles',ratio=False,samefigure=False,withkip=False,sparsity=200,linestyle=['--'],marker=['o'],color=['r'],label=[],markevery=100):


	def plot_intershell_massfrac(self,fig=-1,isotopes=['C-12'],abu_source='mppnp',coord_choice=1,kip=True,marker=['o'],shape=['-']):

		'''
		Plots the intershell abundance at the deepest extend of the TDUP
		during each interpulse phase. Take mass point close

		Uses only the find_TP_attributes
		For Falk

		isotopes: isotopes of intershell, if abu_source is mppnp choose mppnp notation e.g. 'C-12'
			  if abu_source is mesa choose mesa notation e.g. h1
		cood_choice: choice of intershell coordinate, see code, feel free to add your own choice

		abu_source: if 'mppnp' choose abundance from mppnp files; if 'mesa' choose abundance from mesa profile files

		kip: if true, plot also kippenhahn diagram with points where abundances are extracted

		'''
		import mesa as ms

		minis=[]	
		#plt.figure(fig)
		if not fig==-1:
			plt.figure(fig)
		for hh in range(len(self.run_historydata)):
				historydata=self.run_historydata[hh]
				sefiles_out=se(self.runs_H5_out[hh])
				LOGS=self.run_LOGS[hh]
				peak_lum_model,h1_mass_min_DUP_model=historydata.find_TP_attributes(t0_model=historydata.find_first_TP(),fig=3, color='r', marker_type='o')
				mini=float(sefiles_out.se.get('mini'))
				zini=float(sefiles_out.se.get('zini'))
				label=str(mini)+'$M_{\odot}$, Z='+str(zini)
				minis.append(mini)

				TP=[]
				models=[]
				iso_abu=[]
				for k in range(len(isotopes)):
					iso_abu.append([])
				mass_coords=[]
				for k in range(len(h1_mass_min_DUP_model)):
					cycle=int(h1_mass_min_DUP_model[k])
					if abu_source == 'mppnp':
						mass=sefiles_out.get(cycle,'mass')[::-1]	
                                		if coord_choice == 1: #this is Umberto's criteria for intershell coordinate
                                        		c12=sefiles_out.get(cycle,'C-12')[::-1]
                                        		o16=sefiles_out.get(cycle,'O-16')[::-1]
                                        		n14=sefiles_out.get(cycle,'N-14')[::-1]
                                        		he4=sefiles_out.get(cycle,'He-4')[::-1]
                                        		sumO_C_He=np.array(c12)+np.array(o16)+np.array(he4)
                                        		for hh in range(len(mass)):
                                        	        	if n14[hh]<1e-13:
                                        	        	        if he4[hh]>0.1:
                                       		         	            if c12[hh]>0.1:
                                       	                	                if sumO_C_He[hh]>0.9:
                                                	       	                    mass_coord=mass[hh]
                                                	                            break
	
						idx_max=min(range(len(mass)),key=lambda i: abs(mass[i]-mass_coord))
						for hh in range(len(isotopes)):		
							iso_abu[hh].append(sefiles_out.get(cycle,isotopes[hh])[::-1][idx_max])
						TP.append(k)
						mass_coords.append(mass_coord)
						models.append(cycle)
					if abu_source=='mesa':
						#find closest cycle
						import sys
						import cStringIO
						save_stdout = sys.stdout = cStringIO.StringIO()
						pr=ms.mesa_profile(LOGS,cycle)
						sys.stdout = save_stdout
						mass=pr.get('mass')
						if coord_choice == 1:
							c12=pr.get('c12')
							o16=pr.get('o16')
							n14=pr.get('n14')
							he4=pr.get('he4')
							sumO_C_He=np.array(c12)+np.array(o16)+np.array(he4)
	                                 		for hh in range(len(mass)):
                                        	        	if n14[hh]<1e-13:
                                        	        	        if he4[hh]>0.1:
                                       		         	            if c12[hh]>0.1:
                                       	                	                if sumO_C_He[hh]>0.9:
                                                	       	                    mass_coord=mass[hh]
                                                	                            break
						idx_max=min(range(len(mass)),key=lambda i: abs(mass[i]-mass_coord))
						for hh in range(len(isotopes)):
							iso_abu[hh].append(pr.get(isotopes[hh])[idx_max])
						TP.append(k)
						mass_coords.append(mass_coord)	
						models.append(pr.header_attr['model_number'])	
				### plotting
				title=abu_source+' '+label
				for k in range(len(isotopes)):
					plt.plot(models,iso_abu[k],label=isotopes[k],marker=marker[k],linestyle=shape[k])	
					plt.xlabel('Model number')
					plt.title(title)
				plt.ylabel('X$_i$')
				plt.yscale('log')
				plt.legend()
				if kip:
					historydata.kippenhahn(fig+1,'model')
					plt.plot(models,mass_coords,marker='o',label='Extraction point')		
					plt.legend()
					plt.title(title)
	def plot_intershell_massfrac_old(self,fig=-1,isotopes=['C-12'],ratio=False,norm=False,withsurfaceratio=False,marker=['o'],linestyle=['-']):
		''' after TP in intershell

			if ratio==True: need 2 elements in isotopes to create ratio

			norm: Normalize to initial abundance (first model)

		'''
		minis=[]	
		#plt.figure(fig)
		if not fig==-1:
			plt.figure(fig)
		for hh in range(len(self.run_historydata)):
				historydata=self.run_historydata[hh]
				sefiles_out=se(self.runs_H5_out[hh])
				sefiles=sefiles_out
				mini=float(sefiles_out.se.get('mini'))
				zini=float(sefiles_out.se.get('zini'))
				label=str(mini)+'$M_{\odot}$, Z='+str(zini)
				minis.append(mini)
				TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,c13_pocket_min,c13_pocket_max=self.TPAGB_properties(historydata,sefiles_out)
				massfrac_isos_noDUP=[]
				massfrac_isos=[]
				beginAGB_massfrac=[]
				beginAGB_massfrac_ratio=[]
				iso_idx=-1
				massfrac_surf=[]
				for isotope in isotopes:
					#for abu at begin of TPAGB
					mod=TPstart[0]
					mass=sefiles.get(mod,'mass')
					ini_iso=sefiles.get(1,isotope)[-1]
					if norm==False:
						massfrac=sefiles.get(mod,isotope)
					else:
						massfrac=sefiles.get(mod,isotope)/ini_iso

					##Important: choose mass coordinate where to extract abundances
					mass_coord=1./4.*min_m_TP[0]+3./4.*max_m_TP[0]

					idx_max=min(range(len(mass)),key=lambda i: abs(mass[i]-mass_coord)) #   max_m_TP[k-1]))
					beginAGB_massfrac=massfrac[idx_max]

					iso_idx+=1
					#pulse_mDUP=[]
					pulse_number=[]
					TP_prodfac=[]
					c13_pocket_prodfac=[]
					massfrac_all=[]
					massfrac_noDUP=[]
					pulse_noDUP=[]
					c13_pocket_yields=[]
					TP_begin=[]
					TP_end=[]
					c13_pocket_begin=[]
					c13_pocket_end=[]
					massfrac_surf1=[]
					print 'DUPmods',DUPmods
					for k in range(1,len(TPmods)):
						#model=int(float(TPmods[k]))
						#DUP_m, min_m_TP[k] max_m_TP
						print 'Pulse ',k-1,'model ',TPmods[k-1]
						if DUPmods[k-1]==-1:
							print 'no TDUP'
							mod=TPend[k-1]
						        mass=sefiles.get(mod,'mass')
							if norm==False:
								massfrac=sefiles.get(mod,isotope)
							else:
								massfrac=sefiles.get(mod,isotope)/ini_iso
							mass_coord=1./4.*min_m_TP[k-1]+3./4.*max_m_TP[k-1]
							idx_max=min(range(len(mass)),key=lambda i: abs(mass[i]-mass_coord)) #   max_m_TP[k-1]))
							massfrac_noDUP.append(massfrac[idx_max])	
							pulse_noDUP.append(k)
							continue

						###############For PDCZ
						if min_m_TP[k-1]<=min_m_TP[k]:
							lower_m=min_m_TP[k]
						else:
							lower_m=min_m_TP[k-1]
						#if DUP_m[k-1]>=min_m_TP[k] and DUP_m[k-1]<=max_m_TP[k]:
						#        upper_m=DUP_m[k-1]
						#elif DUP_m[k-1]>max_m_TP[k]:
						#        upper_m=max_m_TP[k]
						#else:
						#        upper_m=lower_m
					
						#mod=TPstart[k-1]
						#upper_m_TP=max_m_TP[k-1]
						#mass=sefiles.get(mod,'mass')
						#delta_star_mass=np.array(mass[1:])-np.array(mass[0:-1])
						#massfrac=sefiles.get(mod,isotope)
						#idx_min=min(range(len(mass)),key=lambda i: abs(mass[i]-lower_m))
						#idx_max=min(range(len(mass)),key=lambda i: abs(mass[i]-max_m_TP[k-1]))
						#print 'idx',idx_min,idx_max
						#abu_lastTPstart=sum(np.array(delta_star_mass[idx_min:idx_max])*np.array(massfrac[idx_min:idx_max]))
						mod=TPend[k-1]
						mass=sefiles.get(mod,'mass')
						if norm==False:
							massfrac=sefiles.get(mod,isotope)
						else:
							massfrac=sefiles.get(mod,isotope)/ini_iso
						idx_max=min(range(len(mass)),key=lambda i: abs(mass[i]-DUPm_min_h[k-1])) #   max_m_TP[k-1]))
						massfrac_all.append(massfrac[idx_max])
						pulse_number.append(k)
					        if withsurfaceratio==True:
                                                        massfrac_surf1.append(sefiles.get(mod,isotope)[-1])
						print 'After TP',k-1,'at model',mod
						print 'mass coord',DUPm_min_h[k-1] 
						massfrac_noDUP.append(massfrac[idx_max])
					#massfrac_ratio.append(massfrac[idx_max])
					if fig==-1:
						plt.figure(str(mini)+' 1, ratio='+str(ratio))
					if ratio==False:
						#plt.plot(0,beginAGB_massfrac,marker=marker[iso_idx],linestyle=linestyle[hh],markerfacecolor='w')#,label=isotope+', '+label)
						#pulse_noDUP,
						plt.plot([0]+list(np.array(range(len(massfrac_noDUP)))+1),[beginAGB_massfrac]+massfrac_noDUP,marker=marker[iso_idx],linestyle=linestyle[hh],markerfacecolor='w',label=isotope+', '+label)
						plt.plot(pulse_number,massfrac_all,marker=marker[iso_idx],linestyle=linestyle[hh])
					else:
						massfrac_isos.append(massfrac_all)
						massfrac_isos_noDUP.append([beginAGB_massfrac]+massfrac_noDUP)
						beginAGB_massfrac_ratio.append(beginAGB_massfrac)
						if withsurfaceratio==True:
							massfrac_surf.append(massfrac_surf1)
					plt.legend()
					plt.ylabel('Intershell massfraction')
					plt.xlabel('TP number')
					plt.yscale('log')
					#plt.ylabel('$X/X_{ini}$ in TP')
				if ratio==True:
					print 'ratio ',4./3*np.array(beginAGB_massfrac_ratio[0])/np.array(beginAGB_massfrac_ratio[1])
					#plt.plot(0,4./3*np.array(beginAGB_massfrac_ratio[0])/np.array(beginAGB_massfrac_ratio[1]),marker=marker[iso_idx],linestyle=linestyle[hh],markerfacecolor='w')

					if withsurfaceratio==True:
						plt.plot(pulse_number,4./3.*np.array(massfrac_surf[0])/np.array(massfrac_surf[1]),marker=marker[iso_idx-1],linestyle=linestyle[hh-1],label='Surface, '+label)  #(isotopes[0]+'/'+isotopes[1]+', surface, '+label))
					plt.plot([0]+list(np.array(range(len(massfrac_noDUP)))+1),4./3*np.array(massfrac_isos_noDUP[0])/np.array(massfrac_isos_noDUP[1]),marker=marker[iso_idx],linestyle=linestyle[hh],markerfacecolor='w',label='Intershell, '+label) #(isotopes[0]+'/'+isotopes[1]+', '+label),markerfacecolor='w')

					plt.plot(pulse_number,4./3.*np.array(massfrac_isos[0])/np.array(massfrac_isos[1]),marker=marker[iso_idx],linestyle='None')#linestyle[hh])
                                	
				        plt.legend()
                                        plt.ylabel('Intershell C/O ratio')
                                        plt.xlabel('TP number')
                                        plt.yscale('linear')
					



	def plot_TCEB(self,fig=13,marker=[],linestyle=[],sparsity=10):
		'''Plots DUP vs mass'''
		minis=[]	
		plt.figure(fig)
		for hh in range(len(self.run_historydata)):
				historydata=self.run_historydata[hh]
				sefiles_out=se(self.runs_H5_out[hh])
				sefiles=sefiles_out
				mini=float(sefiles_out.se.get('mini'))
				zini=float(sefiles_out.se.get('zini'))
				label=str(mini)+'$M_{\odot}$, Z='+str(zini)
				minis.append(mini)
				#TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,c13_pocket_min,c13_pocket_max=self.TPAGB_properties(historydata,sefiles_out)
				#get sefiles models
				models=[]
				for k in sefiles.se.cycles:
					models.append(int(float(k)))
				
				modelmin=historydata.find_first_TP()
				h_freecore1=[]
				lum1=[]
				radius1=[]
				maxtemp_CEB1=[]
				idx=0
				massenv1=np.array(historydata.get('conv_mx1_bot'))*np.array(historydata.get('star_mass'))
				
				#mesa_models=historydata.get('model_number')
				model_range=[]
				massenv=[]
				for k in range(0,len(models),sparsity):
					if models[k]<modelmin:
						continue
					#this is for the set1.2 M5 star only
					#if models[k]>104200 and models[k]<105001:
					#	continue
					model_range.append(int(float(models[k])))
					massenv.append(massenv1[k])

				#get sefiles profiles
				temp=np.array(sefiles.get(model_range,'temperature'))			
				mass=sefiles.get(model_range,'mass')
				model_range1=[]
				totmass=[]
				#get the maximum temperature in the AGB stage
				for k in range(len(model_range)):
					#to match sefiles model and mesa model
					idx_mass=min(range(len(mass[k])), key=lambda i: abs(mass[k][i]-(massenv[k]) ))
					if np.log10( temp[k][idx_mass]*sefiles.get('temperature_unit') ) <4:
						continue
					maxtemp_CEB1.append(temp[k][idx_mass]*sefiles.get('temperature_unit')/1e6   )
					totmass.append( mass[k][-1]  )
					model_range1.append(model_range[k])
				print 'finsihed max temp CEB calc'
				plt.figure()
				plt.plot(totmass,maxtemp_CEB1,label=(str(mini)+'$M_{\odot}$'),linestyle=linestyle[hh],marker=marker[hh],markevery=10)
				plt.legend()
				plt.xlabel('$M/M_{\odot}$');plt.ylabel('$T_{CEB}/K$')
				plt.legend(loc=1)
				plt.gca().invert_xaxis()
				plt.show()
				plt.figure()
				plt.plot(model_range1,maxtemp_CEB1,label=(str(mini)+'$M_{\odot}$'),linestyle=linestyle[hh],marker=marker[hh],markevery=10)
				plt.legend()
				plt.xlabel('Model number');plt.ylabel('$T_{CEB}/K$')
				plt.legend(loc=1)
	


			#plt.savefig('T_CEB_'+'1Z'+str(zini)+'.png')



	def plot_TPDCZ(self,fig=13,marker=[]):

		'''Plots TPDCZ vs TP'''
		
		minis=[]
		for hh in range(len(self.run_historydata)):
				historydata=self.run_historydata[hh]
				sefiles_out=se(self.runs_H5_out[hh])
				sefiles=sefiles_out
				mini=float(sefiles_out.se.get('mini'))
				zini=float(sefiles_out.se.get('zini'))
				label=str(mini)+'$M_{\odot}$, Z='+str(zini)
				minis.append(mini)
				TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,c13_pocket_min,c13_pocket_max=self.TPAGB_properties(historydata,sefiles_out)

				###Temperature

				pulse_peak_T=[]
				pulse_number=[]
				plt.figure(fig)
				for k in range(len(TPmods)):
					model=int(float(TPmods[k]))
					T_profile=sefiles_out.get(model,'temperature')*sefiles_out.get('temperature_unit')
					mass=sefiles_out.get(model,'mass')
					#plt.figure(2)
					print 'testse'
					print TPmods[k]
					#plt.plot(mass,T_profile,label='cycle '+str(TPmods[k]),marker=marker[k],markevery=100)
					#plt.figure(11)
					#models=sefiles_out.se.cycles
					#models=map(int,models)
					#print 'model taken: ',models[min(range(len(models)), key=lambda i: abs(models[i]-model))]
					idx_mass=min(range(len(mass)), key=lambda i: abs(mass[i]-min_m_TP[k]))	
					pulse_peak_T.append(max(T_profile[idx_mass:]))
					print k+1
					print '{:.3E}'.format(T_profile[idx_mass])
					pulse_number.append(k+1)

				plt.plot(pulse_number,pulse_peak_T,marker=marker[hh],linestyle='--',label=label)

				plt.legend()
				plt.ylabel('$T_{max,PDCZ}$')
				plt.xlabel('TP number')


	def plot_DUP_mass(self,fig=11,marker=[]):
		'''Plots DUP vs mass'''
		minis=[]
		for hh in range(len(self.run_historydata)):
				historydata=self.run_historydata[hh]
				sefiles_out=se(self.runs_H5_out[hh])
				sefiles=sefiles_out
				mini=float(sefiles_out.se.get('mini'))
				zini=float(sefiles_out.se.get('zini'))
				label=str(mini)+'$M_{\odot}$, Z='+str(zini)
				minis.append(mini)

				TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,c13_pocket_min,c13_pocket_max=self.TPAGB_properties(historydata,sefiles_out)
				pulse_mDUP=[]
				pulse_number=[]
				plt.figure(fig)
				for k in range(len(TPmods)):
					model=int(float(TPmods[k]))
					T_profile=sefiles_out.get(model,'temperature')*sefiles_out.get('temperature_unit')
					#models=sefiles_out.se.cycles
					#models=map(int,models)
					#print 'model taken: ',models[min(range(len(models)), key=lambda i: abs(models[i]-model))]
					pulse_mDUP.append(max_env[k]-DUP_m[k])
					print k+1
					pulse_number.append(k+1)

				plt.plot(pulse_number,pulse_mDUP,marker='o',linestyle='--',label=label)

				plt.legend()
				plt.ylabel('$M_{DUP}$[Msun]')
				plt.xlabel('TP number')


				
	def table_TPAGB_properties(self,filename,sparsity=50):
		'''AGB properties table 1'''


		f1=open(filename,'w')

		aver_interpulse_time=[]
		h_freecore=[]
		lum=[]
		radius=[]
		timefirstTP=[]
		pulse_mDUP=[]
		tot_mDUP=[]
		minis=[]
		num_TP=[]
		num_TP_wDUP=[]
		totmasslost=[]
		final_env_mass=[]
		final_core_mass=[]
		maxtemp_CEB=[]
		T_maxPDCZ=[]
		pulse_size_max=[]
		LHe_max=[]
		L_max=[]

		print 'historydatas'
		print self.run_historydata
		for hh in range(len(self.run_historydata)):
				historydata=self.run_historydata[hh]
				sefiles_out=se(self.runs_H5_out[hh])
				sefiles=sefiles_out
				mini=float(sefiles_out.se.get('mini'))
				zini=float(sefiles_out.se.get('zini'))
				label=str(mini)+'$M_{\odot}$, Z='+str(zini)
				minis.append(mini)

				#TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,DUP_m = self.TPAGB_properties(historydata,sefiles_out)
				TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,c13_pocket_min,c13_pocket_max=self.TPAGB_properties(historydata,sefiles_out)

				logL=historydata.get('log_L')
				logR=historydata.get('log_R')
				starage=historydata.get('star_age')
				modelmin=TPmods[0]
				modelmax=TPmods[-1]
				h_freecore1=[]
				lum1=[]
				radius1=[]
				maxtemp_CEB1=[]
				idx=0
				h1_bdy=historydata.get('h1_boundary_mass')
				massenv=np.array(historydata.get('conv_mx1_bot'))*np.array(historydata.get('star_mass'))

				##########get model range of interest/sparsity
				models=historydata.get('model_number')
				model_range=[]
				for k in range(0,len(models),sparsity):
					if models[k]>modelmax:
						break
					if k<modelmin:
						continue
					model_range.append(int(float(models[k])))
				#for corrupted set1.2 M5 data
                                if ((mini-5.)<0.1) and ((zini-0.02)<0.001):
                                        model_range11=[]
                                        for tt in range(len(model_range)):
                                #this is for the set1.2 M5 star only
                                                if model_range[tt]>104200 and model_range[tt]<105001:
                                                        continue
                                                model_range11.append(model_range[tt])
                                        model_range=model_range11
                                        #       continue


				temp=np.array(sefiles.get(model_range,'temperature'))			
				mass=sefiles.get(model_range,'mass')

				############get the maximum temperature in the AGB stage - measured from first TP to last TP
				for k in range(len(model_range)):
					#print model_range[k]
					#print int(model_range[k])
					#print mass[k]
					idx_mass=min(range(len(mass[k])), key=lambda i: abs(mass[k][i]-(massenv[int(model_range[k])]) ))
					maxtemp_CEB1.append( np.log10( temp[k][idx_mass]*sefiles.get('temperature_unit') )   )
				print 'finsihed max temp CEB calc'
				#plt.figure();plt.plot(model_range,maxtemp_CEB1,label=('M='+str(mini)));plt.legend()
				#plt.xlabel('Model number');plt.ylabel('$T_{CEB}/K$')
				#plt.savefig('T_CEB_M'+str(mini)+'Z'+str(zini)+'.png')
				maxtemp_CEB.append(max(maxtemp_CEB1))
				#print mini,zini,' Found maximum TCEB: ',maxtemp_CEB[-1]
				#############

				###### get mean luminosities between first and last TP
				timesteps=10**historydata.get('log_dt')[modelmin:modelmax]
				lum_mean_weight=  sum( np.array(logL[modelmin:modelmax])*np.array(timesteps)   ) / sum(timesteps)
                                lum.append(lum_mean_weight)

				###### get mean radius between first and last TP
                                radius_mean_weight= sum(np.array( 10**(logR[modelmin:modelmax]) ) *np.array(timesteps)  ) /sum(timesteps)
				radius.append(radius_mean_weight)

				######get h-free core at first TP
				h_freecore.append(h1_bdy[TPmods[0]])     #max(h1_bdy[modelmin:modelmax] ))
				
				#### get time at first TP
				timefirstTP.append(starage[TPmods[0]]/1e6) #in 1e6 yrs

				#### get number of TPs with TDUP & derive DUP masses	
				pulse_mDUP1=[0]	
				firstDUP_k=0
				for k in range(1,len(TPmods)):
					model=int(float(TPmods[k]))
					#print 'model taken: ',models[min(range(len(models)), key=lambda i: abs(models[i]-model))]
					if k>(len(DUPmods)-1):
						#DUP_parameter1.append(0)
						continue
					if (DUPmods[k]==-1):
						#DUP_parameter1.append(0)
						continue
						#pulse_mDUP1.append(max_env[k]-h1_bdy[k])
					m_dup=h1_bdy[TPmods[k]] - h1_bdy[DUPmods[k]]
					pulse_mDUP1.append(m_dup)
					if pulse_mDUP1[-1]<0:
						firstDUP_k=k+1
				print 'first TDUP at model', TPmods[firstDUP_k]
				num_TP_wDUP.append(len(pulse_mDUP1)-1)
			
				#### get maximum temperature in TPs 	
				pulse_peak_T=[]
				#get the peak T at each pulse first
				for k in range(len(TPmods)):
					model=int(float(TPmods[k]))
					T_profile=sefiles_out.get(model,'temperature')*sefiles_out.get('temperature_unit')
					mass=sefiles_out.get(model,'mass')
					idx_mass=min(range(len(mass)), key=lambda i: abs(mass[i]-min_m_TP[k]))
					pulse_peak_T.append(np.log10(max(T_profile[idx_mass:])))
				T_maxPDCZ.append(max(pulse_peak_T))

				#### get maximum helium luminosity between first and last TP
				LHe_max.append(max(historydata.get('log_LHe')[modelmin:modelmax]))

				#### get maximum hydrogen luminosity between first and last TP	
				L_max.append(max(historydata.get('log_L')[modelmin:modelmax]))

				#### get maximum TP size (mass) in 1e-2Msun
				pulse_size=[]
				for ww in range(len(TPmods)):
					pulse_size.append(max_m_TP[ww]-min_m_TP[ww])
					print max_m_TP[ww],min_m_TP[ww]
				pulse_size_max.append(max(pulse_size)*1e2)
				print 'pulse size max',pulse_size_max[-1]

				#### get number of TPs
				num_TP.append(len(TPmods))

				#### get maximum mass dredged-up in 1e-2Msun
				pulse_mDUP.append(max(pulse_mDUP1)/1e-2)
				print 'pulse_mDUP',pulse_mDUP[-1]

				#### get total mass dredged up in 1e-2Msun
				tot_mDUP.append(sum(pulse_mDUP1)/1e-2)
				print 'tot_mDUP',tot_mDUP[-1]

				#### get average interpulse time
				aver_interpulse_time1=[]
				#for k in range(firstDUP_k,len(TPmods)-1):
					#aver_interpulse_time1.append(starage[TPmods[k+1]]-starage[TPmods[k]])
				for k in range(len(TPstart)-1):
					aver_interpulse_time1.append(starage[TPend[k+1]] - starage[TPstart[k]])
				aver_interpulse_time.append(np.mean(aver_interpulse_time1))	
				print 'average interpulse time',aver_interpulse_time[-1]

				#### total mass lost during TPAGB stage		
				totmasslost.append( historydata.get('star_mass')[TPmods[0]] - h1_bdy[TPmods[-1]])	

				#### final envelope mass of whole evolution
				final_env_mass.append( historydata.get('star_mass')[-1]-h1_bdy[-1])

				#### final core mass of whole evolution
				final_core_mass.append(h1_bdy[-1])


		#print 'historydatas'
                #print self.run_historydata

		print '##### TP parameter ######'

		out='M_ini[Msun] &  m_c[Msun] & Log[Lsun] & R[Rsun]  & N_TP   &N_DT& t_TP[1e6yrs]&M_DUPmax[1e-2Msun]&M_tot_mDUP[1e-2Msun]&  t_ip[yr]   & Mtotlost[Msun]  & M_core_final[Msun] & M_env_final[Msun] & T_maxCEB[K] & T_maxPDCZ[K] & M_max_TP[1e-2Msun]  &  logLHe_max   & logL_max'
		print out
		f1.write(out+'\n')
		for hh in range(len(self.run_historydata)):
			out='{:.2f}'.format(minis[hh])+'   & '+	'{:.3f}'.format(h_freecore[hh])+'  & '+ '{:.2f}'.format(lum[hh])+' & '+'{:.0f}'.format(radius[hh])+'&  '+str(int(num_TP[hh]))+'    &  '+str(int(num_TP_wDUP[hh]))+' &  '+'{:.3E}'.format(timefirstTP[hh])+' &  '+'{:.3f}'.format(pulse_mDUP[hh])+'       &  '+'{:.3f}'.format(tot_mDUP[hh])+'         &  '+'{:.0f}'.format(aver_interpulse_time[hh])+'  & '+ '{:.2f}'.format(totmasslost[hh])+'       & '+ '{:.2f}'.format(final_core_mass[hh])+'          & '+ '{:.2f}'.format(final_env_mass[hh])+'         &  '+ '{:.3f}'.format(maxtemp_CEB[hh])+'  &  '+'{:.3f}'.format(T_maxPDCZ[hh])+'  &  '+'{:.3f}'.format(pulse_size_max[hh])+'  &  '+'{:.2f}'.format( LHe_max[hh])+'  &  '+'{:.2f}'.format(L_max[hh])
			print out
			f1.write(out+'\n')
		f1.close()	

	def table2_TPAGB_properties(self,filename):
		'''AGB properties table 1'''

		f1=open(filename,'w')
		t_TP_all=[]
		pulse_peak_T_all=[]
		mtot_pulse_all=[]
		temp_he_interpulse_all=[]
		temp_h_interpulse_all=[]
		temp_env_interpulse_all=[]
		m_h_pulse_all=[]
		m_env_DUP_all=[]
		minis=[]
		min_m_TP_all=[]
		DUP_parameter=[]
		for hh in range(len(self.run_historydata)):
			historydata=self.run_historydata[hh]
			sefiles_out=se(self.runs_H5_out[hh])
			sefiles=sefiles_out
			mini=float(sefiles_out.se.get('mini'))
			zini=float(sefiles_out.se.get('zini'))
			label=str(mini)+'$M_{\odot}$, Z='+str(zini)
			minis.append(mini)
			TP_start,TPmods,max_env,TP_end,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,c13_pocket_min,c13_pocket_max = self.TPAGB_properties(historydata,sefiles_out)
			
			#return
			print 'Tp_end'
			print TP_start
			print 'TP start'
			print TP_end
			#create tables for each run
			massenv=np.array(historydata.get('conv_mx1_bot'))*np.array(historydata.get('star_mass')) 
			for k in range(len(TPmods)-1):
				print 'model',TPmods[k],k+1
				print TP_end[k],TP_start[k+1]
				print min(massenv[TP_end[k]:TP_start[k+1]])

			h1_bdy=historydata.get('h1_boundary_mass')
			starage=historydata.get('star_age')

			#### get times since first TP
			t0=starage[TPmods[0]]
			t_TP=[]
			for k in range(0,len(TPmods)):
				t_TP.append(starage[TPmods[k]] -t0)
			t_TP_all.append(t_TP)

			#### Get T at bottom of PDCZ, mass coord. of he-free core during TP; stellar mass at TP 
			pulse_peak_T=[]
			mtot_pulse=[]
			m_h_pulse=[]
			#T_max
			starmass=historydata.get('star_mass')	
			for k in range(len(TPmods)):
				model=int(float(TPmods[k]))
				T_profile=sefiles_out.get(model,'temperature')*sefiles_out.get('temperature_unit')
				mass=sefiles_out.get(model,'mass')
				idx_mass=min(range(len(mass)), key=lambda i: abs(mass[i]-min_m_TP[k]))
				pulse_peak_T.append(np.log10(max(T_profile[idx_mass:])))
				mtot_pulse.append(starmass[TPmods[k]])	
				m_h_pulse.append(h1_bdy[TPmods[k]])
                        m_h_pulse_all.append(m_h_pulse)
			pulse_peak_T_all.append(pulse_peak_T)
			mtot_pulse_all.append(mtot_pulse)
	
			##### Interpulse time-averaged quantities:
			##### Get T of he_shel, take bottom of PDCZ mass coordinate (temp_he_interpulse)
			##### Get T at h-free core coordinate (temp_h_interpulse)
			##### Get T and the bottom of the conv. envelope (temp_env_interpulse)
			##### Get deepest extend of TDUP (m_env_DUP)

			temp_he_interpulse=[]
			temp_h_interpulse=[]
			temp_env_interpulse=[]
			m_env_DUP=[]
			print 'calculating the interpulse mean values'
			sparsity_interpulse_phase=100
			#for the interpulse phase parameter, take the tempearture at the deepest extend of the TDUP
			for k in range(len(TPmods)-1):
				#between end of last TP and begin of new TP
				models_interpulse=range(TP_end[k],TP_start[k+1],sparsity_interpulse_phase)
				#this is for the set1.2 M5 star only
				if ((mini-5.)<0.1) and ((zini-0.02)<0.001):
					models_interpulse11=[]
					for tt in range(len(models_interpulse)):
                                 		if models_interpulse[tt]>104200 and models_interpulse[tt]<105001:
							continue
						models_interpulse11.append(models_interpulse[tt])
					models_interpulse=models_interpulse11
	
				timesteps=10**historydata.get('log_dt')
				temp_he_interpulse1=[]
				temp_h_interpulse1=[]
				temp_env_interpulse1=[]
				sum_timesteps=0
				for model_mid in models_interpulse:
					#print 'try model',model_mid
					mass=sefiles.get(model_mid,'mass')	
					temp=np.log10(sefiles.get('temperature_unit')*np.array(sefiles.get(model_mid,'temperature')))
					#min_m_TP[k]
					idx_mass=min(range(len(mass)), key=lambda i: abs(mass[i]-min_m_TP[k]))
					while True:
						#for some weird error
						if (temp[idx_mass]*timesteps[model_mid]) >0:
							break
						model_mid+=20
		                                mass=sefiles.get(model_mid,'mass')
						idx_mass=min(range(len(mass)), key=lambda i: abs(mass[i]-min_m_TP[k]))
	                                        temp=np.log10(sefiles.get('temperature_unit')*np.array(sefiles.get(model_mid,'temperature')))
					sum_timesteps+=timesteps[model_mid]						
					temp_he_interpulse1.append( temp[idx_mass]*timesteps[model_mid]  )
					idx_mass=min(range(len(mass)), key=lambda i: abs(mass[i]-h1_bdy[model_mid]))	
					temp_h_interpulse1.append(temp[idx_mass]*timesteps[model_mid])
					idx_mass=min(range(len(mass)), key=lambda i: abs(mass[i]-(massenv[model_mid]) ))
					temp_env_interpulse1.append(temp[idx_mass]*timesteps[model_mid])
				temp_he_interpulse.append( sum(temp_he_interpulse1)/sum_timesteps)
				temp_h_interpulse.append( sum(temp_h_interpulse1)/sum_timesteps)	
				temp_env_interpulse.append(sum(temp_env_interpulse1)/sum_timesteps)

				#calc DUP depth 
				m_env_DUP.append(min(massenv[TP_end[k]:TP_start[k+1]])) #*historydata.get('star_mass')[DUPmods[k]])
				if temp_h_interpulse[-1]<0:
					print 'temp_h_interpulse########################3'
					print temp_h_interpulse
					print 'timesteps sum',sum_timesteps
					print 'model mid'
					print model_mid
					print 'temp_env_interpulse1'
					print temp_env_interpulse1

			temp_he_interpulse.append(0)
			temp_h_interpulse.append(0)
			temp_env_interpulse.append(0)
			m_env_DUP.append(0)	
			temp_he_interpulse_all.append(temp_he_interpulse)
			temp_h_interpulse_all.append(temp_h_interpulse)
			temp_env_interpulse_all.append(temp_env_interpulse)
			m_env_DUP_all.append(m_env_DUP)
	

			#### Get DUP parameter lambda
			DUP_parameter1=[0]
			for k in range(1,len(TPmods)):
				if k>(len(DUPmods)-1):
                                        DUP_parameter1.append(0)
                                        continue					
				if (DUPmods[k]==-1):
					DUP_parameter1.append(0)
					continue
				#first  pulse where last pulse has no TDUP
				if (DUPmods[k-1]==-1):
					m_h=h1_bdy[TPmods[k]] - h1_bdy[TPmods[k-1]]
				else:
					m_h=h1_bdy[TPmods[k]] - h1_bdy[DUPmods[k-1]]
				m_dup=h1_bdy[TPmods[k]] - h1_bdy[DUPmods[k]]
				DUP_parameter1.append(m_dup/m_h)
			DUP_parameter.append(DUP_parameter1)

			####get min mass coord of TP
			min_m_TP_all.append(min_m_TP)


		out= '##### TP parameter ######'
		print out
		f1.write(out)
		for hh in range(len(self.run_historydata)):
			out='\n'+'Initial mass '+str(minis[hh])+'\n'
			print out
			f1.write(out)
			out= 'TP&DUP_lambda &t_TP[year]  &  T_FBOT[K]&  T_HES[k] &  T_HS[K]      &   T_CEB[K] &m_FBOT[Msun]& m_HTP[Msun]      & m_Dmax[Msun]       & Mstar'
			print out
			f1.write(out+'\n')
			'''
			print len(DUP_parameter[hh])
			print len(t_TP_all[hh])
			print len(pulse_peak_T_all[hh])
			print len(temp_he_interpulse_all[hh])
			print len(temp_h_interpulse_all[hh])
			print len(temp_env_interpulse_all[hh])
			print len(min_m_TP_all[hh])
			print len(m_h_pulse_all[hh])
			print len(m_env_DUP_all[hh])
			print len(mtot_pulse_all[hh])
			'''
			for k in range(len(pulse_peak_T_all[hh])):
				out= str((k+1))+abs(2-len(str(k+1)))*' '+'&'+'{:.2f}'.format(DUP_parameter[hh][k])+'  & '+'{:.2E}'.format(t_TP_all[hh][k])+'  &  '+	'{:.2f}'.format(pulse_peak_T_all[hh][k])+'&  '+ '{:.2f}'.format(temp_he_interpulse_all[hh][k])+'&  '+'{:.2f}'.format(temp_h_interpulse_all[hh][k])+'    &  '+'{:.2f}'.format(temp_env_interpulse_all[hh][k])+' &  '+'{:.4f}'.format(min_m_TP_all[hh][k])+' &  '+'{:.4f}'.format(m_h_pulse_all[hh][k])+'       &  '+'{:.4f}'.format(m_env_DUP_all[hh][k])+'         &  '+'{:.3f}'.format(mtot_pulse_all[hh][k])
				f1.write(out+'\n')
				print out
		f1.close()

	

	def TPAGB_properties(self,historydata,sefiles):
		
		peak_lum_model,h1_mass_min_DUP_model=historydata.find_TP_attributes(t0_model=historydata.find_first_TP(),fig=3,color='r', marker_type='o')

		print 'first tp'
		print historydata.find_first_TP()
		print 'peak lum mmmodel'
		print  peak_lum_model
		print h1_mass_min_DUP_model

		TPmods=peak_lum_model
	
		DUPmods=h1_mass_min_DUP_model	
		DUPmods1=[]
		for k in range(len(DUPmods)):

                        testdupmod=int(float(DUPmods[k]))+100
                        if not testdupmod > historydata.get('model_number')[-1]:
                                #DUPmods1.append(self.get('model_number')[-1])
                                #else:
                                DUPmods1.append(int(float(DUPmods[k]))+100) #to exclude HBB? effects
		DUPmods=DUPmods1
		


		TPstart=[]
		#find beginning of TP, goes from TP peak backwards
		# find end of PDCZ by seeking from TP peak and checking mx2_bot:
		models=historydata.get('model_number')
		mx2b_array=historydata.get('conv_mx2_bot')
		mx2t_array=historydata.get('conv_mx2_top')
		massbot=mx2b_array#*historydata.header_attr['initial_mass']
		masstop=mx2t_array#*historydata.header_attr['initial_mass']
		massenv=np.array(historydata.get('conv_mx1_bot'))*np.array(historydata.get('star_mass'))   #*historydata.header_attr['initial_mass']
		
		h1_bdy=historydata.get('h1_boundary_mass')
		for k in range(len(TPmods)):
			idx=list(models).index(TPmods[k])
			mx2b=mx2b_array[:idx]
			for i in range(len(mx2b)-1,0,-1):
				if mx2b[i]==0.:
				    startTP=models[i]
				    TPstart.append(int(float(startTP)))
				    break
		#Find end of TP, goes from TP forwards:
		TPend=[]
		max_m_TP=[]
		min_m_TP=[]
		DUP_m=[]
		TP_max_env=[]
		DUPm_min_h=[]
		c13_pocket_min=[]
		c13_pocket_max=[]
		flagdecline=False
		for k in range(len(TPmods)):
			idx=list(models).index(TPmods[k])
			mx2b=mx2b_array[idx:]
			mx2t=mx2t_array[idx:]
			refsize=mx2t[0]-mx2b[0]
			for i in range(len(mx2b)):
				if i==0:
					continue
				if ((mx2t[i]-mx2b[i])<(0.5*refsize)) and (flagdecline==False):
					flagdecline=True
					refmasscoord=mx2t[i]
					print 'flagdecline to true'
					continue
				if flagdecline==True:
					if (mx2t[i]-mx2b[i])<(0.1*refsize):
						#for the massive and HDUP AGB's where PDCZ conv zone becomes the Hdup CONV ZONE
						if refmasscoord<mx2t[i]:
							endTP=models[idx+i-1]
							TPend.append(int(float(endTP)))
							#print 'HDUp, TP end',endTP
							break
						if (mx2t[i]-mx2b[i])<1e-5:
							endTP=models[idx+i-1]
							TPend.append(int(float(endTP)))
							#print 'normal TPend',endTP
							break
				'''
				if max(mx2t[0:(i-1)])>mx2t[i]:
					(max(mx2t[0:(i-1)]) - min(mx2b[0:(i-1)]))
					flag=True
					continue
				if flag==True:
                                    endidx=idx+i
                                    endTP=models[endidx]
                                    TPend.append(int(float(endTP)))								
		
				if (mx2t[i]-mx2b[i])<1e-5:			#mx2b[i])==0.:
				    endidx=idx+i
				    endTP=models[endidx]
				    TPend.append(int(float(endTP)))
				    break
				'''
			#print 'found TP boundaries',TPstart[-1],TPend[-1]
		#find max and minimum mass coord of TP at max Lum
			mtot=historydata.get('star_mass')
			masstop_tot=np.array(masstop)*np.array(mtot)
			idx_tpext=list(masstop_tot).index(max(masstop_tot[TPstart[k]:(TPend[k]-10)]))
			#print 'TP',k+1,TPmods[k]
			##print TPstart[k],TPend[k]
			#print 'INDEX',idx_tpext,models[idx_tpext]
			#print max(masstop_tot[TPstart[k]:(TPend[k]-10)])
			mtot=historydata.get('star_mass')[idx_tpext]
			max_m_TP.append(masstop[idx_tpext]*mtot)
			min_m_TP.append(massbot[idx_tpext]*mtot)
			
			TP_max_env.append(massenv[idx_tpext])#*mtot)
			if k> (len(DUPmods)-1):
				continue		
			idx=list(models).index(DUPmods[k])
			mtot=historydata.get('star_mass')[idx]
			#DUP_m.append(h1_bdy[idx])#*mtot)
		#######identify if it is really a TDUP, Def.
			#if not (no c13 pocket)
			if h1_bdy[idx]>=max_m_TP[-1]:
				print 'Pulse',k+1,'model',TPmods[k],'skip'
				print h1_bdy[idx],max_m_TP[-1]
				DUPmods[k] = -1
				DUPm_min_h.append( -1)  
				c13_pocket_min.append(-1)
				c13_pocket_max.append(-1)
				continue

			DUPm_min_h.append(h1_bdy[idx])
			###find c13 pocket area from se files (precise values)		
			h1=sefiles.get(int(models[idx]),'H-1')
			mass=sefiles.get(int(models[idx]),'mass')
			for tt in range(len(h1)):
				if h1[tt]>1e-9:
					c13_pocket_min.append(mass[tt] -0.00015)
					break
			for tt in range(len(h1)):
				if h1[tt]>0.999*h1[-1]:
					c13_pocket_max.append(mass[tt])
					break
		for k in range(len(TPmods)):
			print '#############'
			print 'TP ',k+1
			print 'Start: ',TPstart[k]
			print 'Peak' , TPmods[k],TP_max_env[k]
			print '(conv) PDCZ size: ',min_m_TP[k],' till ',max_m_TP[k]
			print 'End',TPend[k]
			if k <=(len(DUPmods)-1):
				print len(DUPmods),k
				print 'DUP max',DUPmods[k]
				print DUPm_min_h[k]
				print 'C13 mass range',c13_pocket_min[k],c13_pocket_max[k]
			else:
				print 'no DUP'

		return TPstart,TPmods,TP_max_env,TPend,min_m_TP,max_m_TP,DUPmods,DUPm_min_h,c13_pocket_min,c13_pocket_max







#######################Global functions ###########################

def get_ini_elem_abu(iniabupath='/home/christian/PPN/forum.astro.keele.ac.uk/frames/mppnp/USEEPP/iniab2.0E-02GN93.ppn',element='H'):
	'''

		Return abundance of element from .ppn file iniabupath
	'''

	abu=0
	a=u.iniabu(iniabupath)
	for iso in a.names:
		iso_namescheme=iso.capitalize()[:2].split()[0]+"-"+iso[3:].split()[-1]
		if element == iso_namescheme.split('-')[0]:
			abu+=a.habu[iso]
	print element,abu
	return abu




	isotopes=sefiles.se.isotopes
	if type(cycles)==int:
		len_cycles=1
	else:
		len_cycles=len(cycles)
	type_hdf5= type(sefiles.get(sefiles.se.cycles[0],"H-1"))
	if type_hdf5 == np.float64:
		abu=np.zeros(len_cycles)
		for i in range(len(isotopes)):
			if element == isotopes[i]:
				abu+=sefiles.get(cycles,isotopes[i])
	## hdf5 out
	#else:
		#abu

	return abu


def get_stable(specie_list,get_elements=True):

	'''
		Input isotope list or element list, get stable list back
		if get_elements=True return elements of all stables in isotope_list

	'''

	stable_list=[]
	for specie in specie_list:
		if is_stable(specie) == 't':
			if get_elements==True and  len(specie.split("-"))==2:
				if specie.split("-")[0] not in stable_list:
					stable_list.append(specie.split("-")[0])	
			else:
				stable_list.append(specie)
	
	return stable_list


def is_stable(isotope):

	'''
		Input either isotope e.g. "C-12" or element "C".
		if stable returns "t"
		Recent correction in this list: Na-23 was set f
	'''                
	isotopes={"N-1":"t",
	"H-1":"t","H-2":"t","He-3":"t","He-4":"t","Be-7":"f",
	"B-8":"f","Li-7":"t","C-11":"f","B-11":"t","C-12":"t",
	"C-13":"t","N-13":"f","N-14":"t","C-14":"f","N-15":"t",
	"O-16":"t","O-17":"t","O-18":"t","F-17":"f","F-18":"f",
	"F-19":"t","Ne-20":"t","Ne-21":"t","Ne-22":"t","Na-22":"f",
	"Na-23":"t","Mg-23":"f","Mg-24":"t","Mg-25":"t","Mg-26":"t",
	"Al-26":"f","Al-27":"t","Si-27":"f","Si-28":"t","Si-29":"t",
	"Si-30":"t","P-31":"t","S-31":"f","Be-8":"f","O-14":"f",
	"O-15":"f","Na-21":"f","Al-25":"f","P-29":"f","P-30":"f",
	"Pb-206":"t","Pb-207":"t","Bi-211":"f","Po-210":"f","F-20":"f",
	"Ne-19":"f","Na-24":"f","Mg-27":"f","Mg-28":"f","Al-28":"f",
	"Al-29":"f","Si-31":"f","Si-32":"f","P-32":"f","P-33":"f",
	"P-34":"f","P-35":"f","S-32":"t","S-33":"t","S-34":"t",
	"S-35":"f","S-36":"t","S-37":"f","S-38":"f","Cl-34":"f",
	"Cl-35":"t","Cl-36":"l","Cl-37":"t","Cl-38":"f","Cl-39":"f",
	"Cl-40":"f","Ar-35":"f","Ar-36":"t","Ar-37":"f","Ar-38":"t",
	"Ar-39":"f","Ar-40":"t","Ar-41":"f","Ar-42":"f","Ar-43":"f",
	"Ar-44":"f","K-38":"f","K-39":"t","K-40":"t","K-41":"t",
	"K-42":"f","K-43":"f","K-44":"f","K-45":"f","K-46":"f",
	"Ca-39":"f","Ca-40":"t","Ca-41":"f","Ca-42":"t","Ca-43":"t",
	"Ca-44":"t","Ca-45":"f","Ca-46":"t","Ca-47":"f","Ca-48":"t",
	"Ca-49":"f","Sc-43":"f","Sc-44":"f","Sc-45":"t","Sc-46":"f",
	"Sc-47":"f","Sc-48":"f","Sc-49":"f","Sc-50":"f","Ti-44":"f",
	"Ti-45":"f","Ti-46":"t","Ti-47":"t","Ti-48":"t","Ti-49":"t",
	"Ti-50":"t","Ti-51":"f","Ti-52":"f","V-47":"f","V-48":"f",
	"V-49":"f","V-50":"t","V-51":"t","V-52":"f","V-53":"f",
	"Cr-48":"f","Cr-49":"f","Cr-50":"t","Cr-51":"f","Cr-52":"t",
	"Cr-53":"t","Cr-54":"t","Cr-55":"f","Cr-56":"f","Mn-51":"f",
	"Mn-52":"f","Mn-53":"f","Mn-54":"f","Mn-55":"t","Mn-56":"f",
	"Mn-57":"f","Fe-52":"f","Fe-53":"f","Fe-54":"t","Fe-55":"f",
	"Fe-56":"t","Fe-57":"t","Fe-58":"t","Fe-59":"f","Fe-60":"f",
	"Fe-61":"f","Co-55":"f","Co-56":"f","Co-57":"f","Co-58":"f",
	"Co-59":"t","Co-60":"f","Co-61":"f","Co-62":"f","Co-63":"f",
	"Ni-56":"f","Ni-57":"f","Ni-58":"t","Ni-59":"l","Ni-60":"t",
	"Ni-61":"t","Ni-62":"t","Ni-63":"l","Ni-64":"t","Ni-65":"f",
	"Ni-66":"f","Ni-67":"f","Ni-68":"f","Cu-60":"f","Cu-61":"f",
	"Cu-62":"f","Cu-63":"t","Cu-64":"f","Cu-65":"t","Cu-66":"f",
	"Cu-67":"f","Cu-68":"f","Cu-69":"f","Cu-70":"f","Cu-71":"f",
	"Zn-62":"f","Zn-63":"f","Zn-64":"t","Zn-65":"f","Zn-66":"t",
	"Zn-67":"t","Zn-68":"t","Zn-69":"f","Zn-70":"t","Zn-71":"f",
	"Zn-72":"f","Zn-73":"f","Zn-74":"f","Ga-65":"f","Ga-66":"f",
	"Ga-67":"f","Ga-68":"f","Ga-69":"t","Ga-70":"f","Ga-71":"t",
	"Ga-72":"f","Ga-73":"f","Ga-74":"f","Ga-75":"f","Ge-66":"f",
	"Ge-67":"f","Ge-68":"f","Ge-69":"f","Ge-70":"t","Ge-71":"f",
	"Ge-72":"t","Ge-73":"t","Ge-74":"t","Ge-75":"f","Ge-76":"t",
	"Ge-77":"f","Ge-78":"f","As-69":"f","As-70":"f","As-71":"f",
	"As-72":"f","As-73":"f","As-74":"f","As-75":"t","As-76":"f",
	"As-77":"f","As-78":"f","As-79":"f","As-80":"f","As-81":"f",
	"Se-72":"f","Se-73":"f","Se-74":"t","Se-75":"f","Se-76":"t",
	"Se-77":"t","Se-78":"t","Se-79":"l","Se-80":"t","Se-81":"f",
	"Se-82":"t","Se-83":"f","Se-84":"f","Br-74":"f","Br-75":"f",
	"Br-76":"f","Br-77":"f","Br-78":"f","Br-79":"t","Br-80":"f",
	"Br-81":"t","Br-82":"f","Br-83":"f","Br-84":"f","Br-85":"f",
	"Br-86":"f","Br-87":"f","Kr-76":"f","Kr-77":"f","Kr-78":"t",
	"Kr-79":"f","Kr-80":"t","Kr-81":"f","Kr-82":"t","Kr-83":"t",
	"Kr-84":"t","Kr-85":"l","Kr-86":"t","Kr-87":"f","Kr-88":"f",
	"Kr-89":"f","Kr-90":"f","Rb-79":"f","Rb-80":"f","Rb-81":"f",
	"Rb-82":"f","Rb-83":"f","Rb-84":"f","Rb-85":"t","Rb-86":"f",
	"Rb-87":"t","Rb-88":"f","Rb-89":"f","Rb-90":"f","Rb-91":"f",
	"Sr-80":"f","Sr-81":"f","Sr-82":"f","Sr-83":"f","Sr-84":"t",
	"Sr-85":"f","Sr-86":"t","Sr-87":"t","Sr-88":"t","Sr-89":"f",
	"Sr-90":"f","Sr-91":"f","Sr-92":"f","Sr-93":"f","Sr-94":"f",
	"Y-85":"f","Y-86":"f","Y-87":"f","Y-88":"f","Y-89":"t",
	"Y-90":"f","Y-91":"f","Y-92":"f","Y-93":"f","Y-94":"f",
	"Y-95":"f","Y-96":"f","Zr-86":"f","Zr-87":"f","Zr-88":"f",
	"Zr-89":"f","Zr-90":"t","Zr-91":"t","Zr-92":"t","Zr-93":"l",
	"Zr-94":"t","Zr-95":"f","Zr-96":"t","Zr-97":"f","Zr-98":"f",
	"Nb-89":"f","Nb-90":"f","Nb-91":"f","Nb-92":"f","Nb-93":"t",
	"Nb-94":"f","Nb-95":"f","Nb-96":"f","Nb-97":"f","Nb-98":"f",
	"Nb-99":"f","Mo-90":"f","Mo-91":"f","Mo-92":"t","Mo-93":"f",
	"Mo-94":"t","Mo-95":"t","Mo-96":"t","Mo-97":"t","Mo-98":"t",
	"Mo-99":"f","Mo-100":"t","Mo-101":"f","Mo-102":"f","Tc-93":"f",
	"Tc-94":"f","Tc-95":"f","Tc-96":"f","Tc-97":"f","Tc-98":"l",
	"Tc-99":"f","Tc-100":"f","Tc-101":"f","Tc-102":"f","Tc-103":"f",
	"Tc-104":"f","Tc-105":"f","Ru-94":"f","Ru-95":"f","Ru-96":"t",
	"Ru-97":"f","Ru-98":"t","Ru-99":"t","Ru-100":"t","Ru-101":"t",
	"Ru-102":"t","Ru-103":"f","Ru-104":"t","Ru-105":"f","Ru-106":"f",
	"Rh-98":"f","Rh-99":"f","Rh-100":"f","Rh-101":"f","Rh-102":"f",
	"Rh-103":"t","Rh-104":"f","Rh-105":"f","Rh-106":"f","Rh-107":"f",
	"Rh-108":"f","Pd-99":"f","Pd-100":"f","Pd-101":"f","Pd-102":"t",
	"Pd-103":"f","Pd-104":"t","Pd-105":"t","Pd-106":"t","Pd-107":"l",
	"Pd-108":"t","Pd-109":"f","Pd-110":"t","Pd-111":"f","Pd-112":"f",
	"Ag-101":"f","Ag-102":"f","Ag-103":"f","Ag-104":"f","Ag-105":"f",
	"Ag-106":"f","Ag-107":"t","Ag-108":"f","Ag-109":"t","Ag-110":"f",
	"Ag-111":"f","Ag-112":"f","Ag-113":"f","Cd-102":"f","Cd-103":"f",
	"Cd-104":"f","Cd-105":"f","Cd-106":"t","Cd-107":"f","Cd-108":"t",
	"Cd-109":"f","Cd-110":"t","Cd-111":"t","Cd-112":"t","Cd-113":"t",
	"Cd-114":"t","Cd-115":"f","Cd-116":"t","Cd-117":"f","Cd-118":"f",
	"In-106":"f","In-107":"f","In-108":"f","In-109":"f","In-110":"f",
	"In-111":"f","In-112":"f","In-113":"t","In-114":"f","In-115":"t",
	"In-116":"f","In-117":"f","In-118":"f","In-119":"f","Sn-108":"f",
	"Sn-109":"f","Sn-110":"f","Sn-111":"f","Sn-112":"t","Sn-113":"f",
	"Sn-114":"t","Sn-115":"t","Sn-116":"t","Sn-117":"t","Sn-118":"t",
	"Sn-119":"t","Sn-120":"t","Sn-121":"f","Sn-122":"t","Sn-123":"f",
	"Sn-124":"t","Sn-125":"f","Sn-126":"f","Sn-127":"f","Sn-128":"f",
	"Sn-129":"f","Sn-130":"f","Sb-112":"f","Sb-113":"f","Sb-114":"f",
	"Sb-115":"f","Sb-116":"f","Sb-117":"f","Sb-118":"f","Sb-119":"f",
	"Sb-120":"f","Sb-121":"t","Sb-122":"f","Sb-123":"t","Sb-124":"f",
	"Sb-125":"f","Sb-126":"f","Sb-127":"f","Sb-128":"f","Sb-129":"f",
	"Sb-130":"f","Sb-131":"f","Sb-132":"f","Sb-133":"f","Te-114":"f",
	"Te-115":"f","Te-116":"f","Te-117":"f","Te-118":"f","Te-119":"f",
	"Te-120":"t","Te-121":"f","Te-122":"t","Te-123":"t","Te-124":"t",
	"Te-125":"t","Te-126":"t","Te-127":"f","Te-128":"t","Te-129":"f",
	"Te-130":"t","Te-131":"f","Te-132":"f","Te-133":"f","Te-134":"f",
	"I-117":"f","I-118":"f","I-119":"f","I-120":"f","I-121":"f",
	"I-122":"f","I-123":"f","I-124":"f","I-125":"f","I-126":"f",
	"I-127":"t","I-128":"f","I-129":"f","I-130":"f","I-131":"f",
	"I-132":"f","I-133":"f","I-134":"f","I-135":"f","Xe-118":"f",
	"Xe-119":"f","Xe-120":"f","Xe-121":"f","Xe-122":"f","Xe-123":"f",
	"Xe-124":"t","Xe-125":"f","Xe-126":"t","Xe-127":"f","Xe-128":"t",
	"Xe-129":"t","Xe-130":"t","Xe-131":"t","Xe-132":"t","Xe-133":"f",
	"Xe-134":"t","Xe-135":"f","Xe-136":"t","Xe-137":"f","Xe-138":"f",
	"Cs-123":"f","Cs-124":"f","Cs-125":"f","Cs-126":"f","Cs-127":"f",
	"Cs-128":"f","Cs-129":"f","Cs-130":"f","Cs-131":"f","Cs-132":"f",
	"Cs-133":"t","Cs-134":"f","Cs-135":"f","Cs-136":"f","Cs-137":"f",
	"Cs-138":"f","Cs-139":"f","Ba-124":"f","Ba-125":"f","Ba-126":"f",
	"Ba-127":"f","Ba-128":"f","Ba-129":"f","Ba-130":"t","Ba-131":"f",
	"Ba-132":"t","Ba-133":"f","Ba-134":"t","Ba-135":"t","Ba-136":"t",
	"Ba-137":"t","Ba-138":"t","Ba-139":"f","Ba-140":"f","Ba-141":"f",
	"Ba-142":"f","La-127":"f","La-128":"f","La-129":"f","La-130":"f",
	"La-131":"f","La-132":"f","La-133":"f","La-134":"f","La-135":"f",
	"La-136":"f","La-137":"f","La-138":"t","La-139":"t","La-140":"f",
	"La-141":"f","La-142":"f","La-143":"f","Ce-130":"f","Ce-131":"f",
	"Ce-132":"f","Ce-133":"f","Ce-134":"f","Ce-135":"f","Ce-136":"t",
	"Ce-137":"f","Ce-138":"t","Ce-139":"f","Ce-140":"t","Ce-141":"f",
	"Ce-142":"t","Ce-143":"f","Ce-144":"f","Ce-145":"f","Ce-146":"f",
	"Pr-133":"f","Pr-134":"f","Pr-135":"f","Pr-136":"f","Pr-137":"f",
	"Pr-138":"f","Pr-139":"f","Pr-140":"f","Pr-141":"t","Pr-142":"f",
	"Pr-143":"f","Pr-144":"f","Pr-145":"f","Pr-146":"f","Pr-147":"f",
	"Pr-148":"f","Pr-149":"f","Nd-134":"f","Nd-135":"f","Nd-136":"f",
	"Nd-137":"f","Nd-138":"f","Nd-139":"f","Nd-140":"f","Nd-141":"f",
	"Nd-142":"t","Nd-143":"t","Nd-144":"t","Nd-145":"t","Nd-146":"t",
	"Nd-147":"f","Nd-148":"t","Nd-149":"f","Nd-150":"t","Nd-151":"f",
	"Nd-152":"f","Pm-137":"f","Pm-138":"f","Pm-139":"f","Pm-140":"f",
	"Pm-141":"f","Pm-142":"f","Pm-143":"f","Pm-144":"f","Pm-145":"l",
	"Pm-146":"l","Pm-147":"l","Pm-148":"f","Pm-149":"f","Pm-150":"f",
	"Pm-151":"f","Pm-152":"f","Pm-153":"f","Pm-154":"f","Sm-140":"f",
	"Sm-141":"f","Sm-142":"f","Sm-143":"f","Sm-144":"t","Sm-145":"f",
	"Sm-146":"l","Sm-147":"t","Sm-148":"t","Sm-149":"t","Sm-150":"t",
	"Sm-151":"l","Sm-152":"t","Sm-153":"f","Sm-154":"t","Sm-155":"f",
	"Sm-156":"f","Sm-157":"f","Sm-158":"f","Eu-143":"f","Eu-144":"f",
	"Eu-145":"f","Eu-146":"f","Eu-147":"f","Eu-148":"f","Eu-149":"f",
	"Eu-150":"l","Eu-151":"t","Eu-152":"f","Eu-153":"t","Eu-154":"f",
	"Eu-155":"f","Eu-156":"f","Eu-157":"f","Eu-158":"f","Eu-159":"f",
	"Gd-144":"f","Gd-145":"f","Gd-146":"f","Gd-147":"f","Gd-148":"l",
	"Gd-149":"f","Gd-150":"l","Gd-151":"f","Gd-152":"t","Gd-153":"f",
	"Gd-154":"t","Gd-155":"t","Gd-156":"t","Gd-157":"t","Gd-158":"t",
	"Gd-159":"f","Gd-160":"t","Gd-161":"f","Gd-162":"f","Tb-147":"f",
	"Tb-148":"f","Tb-149":"f","Tb-150":"f","Tb-151":"f","Tb-152":"f",
	"Tb-153":"f","Tb-154":"f","Tb-155":"f","Tb-156":"f","Tb-157":"l",
	"Tb-158":"f","Tb-159":"t","Tb-160":"f","Tb-161":"f","Tb-162":"f",
	"Tb-163":"f","Tb-164":"f","Tb-165":"f","Dy-148":"f","Dy-149":"f",
	"Dy-150":"f","Dy-151":"f","Dy-152":"f","Dy-153":"f","Dy-154":"l",
	"Dy-155":"f","Dy-156":"t","Dy-157":"f","Dy-158":"t","Dy-159":"f",
	"Dy-160":"t","Dy-161":"t","Dy-162":"t","Dy-163":"t","Dy-164":"t",
	"Dy-165":"f","Dy-166":"f","Dy-167":"f","Dy-168":"f","Ho-153":"f",
	"Ho-154":"f","Ho-155":"f","Ho-156":"f","Ho-157":"f","Ho-158":"f",
	"Ho-159":"f","Ho-160":"f","Ho-161":"f","Ho-162":"f","Ho-163":"f",
	"Ho-164":"f","Ho-165":"t","Ho-166":"f","Ho-167":"f","Ho-168":"f",
	"Ho-169":"f","Er-154":"f","Er-155":"f","Er-156":"f","Er-157":"f",
	"Er-158":"f","Er-159":"f","Er-160":"f","Er-161":"f","Er-162":"t",
	"Er-163":"f","Er-164":"t","Er-165":"f","Er-166":"t","Er-167":"t",
	"Er-168":"t","Er-169":"f","Er-170":"t","Er-171":"f","Er-172":"f",
	"Er-173":"f","Er-174":"f","Er-175":"f","Tm-159":"f","Tm-160":"f",
	"Tm-161":"f","Tm-162":"f","Tm-163":"f","Tm-164":"f","Tm-165":"f",
	"Tm-166":"f","Tm-167":"f","Tm-168":"f","Tm-169":"t","Tm-170":"f",
	"Tm-171":"f","Tm-172":"f","Tm-173":"f","Tm-174":"f","Tm-175":"f",
	"Tm-176":"f","Yb-160":"f","Yb-161":"f","Yb-162":"f","Yb-163":"f",
	"Yb-164":"f","Yb-165":"f","Yb-166":"f","Yb-167":"f","Yb-168":"t",
	"Yb-169":"f","Yb-170":"t","Yb-171":"t","Yb-172":"t","Yb-173":"t",
	"Yb-174":"t","Yb-175":"f","Yb-176":"t","Yb-177":"f","Yb-178":"f",
	"Yb-179":"f","Yb-180":"f","Lu-165":"f","Lu-166":"f","Lu-167":"f",
	"Lu-168":"f","Lu-169":"f","Lu-170":"f","Lu-171":"f","Lu-172":"f",
	"Lu-173":"f","Lu-174":"f","Lu-175":"t","Lu-176":"t","Lu-177":"f",
	"Lu-178":"f","Lu-179":"f","Lu-180":"f","Lu-181":"f","Lu-182":"f",
	"Hf-166":"f","Hf-167":"f","Hf-168":"f","Hf-169":"f","Hf-170":"f",
	"Hf-171":"f","Hf-172":"f","Hf-173":"f","Hf-174":"t","Hf-175":"f",
	"Hf-176":"t","Hf-177":"t","Hf-178":"t","Hf-179":"t","Hf-180":"t",
	"Hf-181":"f","Hf-182":"f","Hf-183":"f","Hf-184":"f","Hf-185":"f",
	"Ta-169":"f","Ta-170":"f","Ta-171":"f","Ta-172":"f","Ta-173":"f",
	"Ta-174":"f","Ta-175":"f","Ta-176":"f","Ta-177":"f","Ta-178":"f",
	"Ta-179":"f","Ta-180":"t","Ta-181":"t","Ta-182":"f","Ta-183":"f",
	"Ta-184":"f","Ta-185":"f","Ta-186":"f","W-172":"f","W-173":"f",
	"W-174":"f","W-175":"f","W-176":"f","W-177":"f","W-178":"f",
	"W-179":"f","W-180":"t","W-181":"f","W-182":"t","W-183":"t",
	"W-184":"t","W-185":"f","W-186":"t","W-187":"f","W-188":"f",
	"W-189":"f","W-190":"f","Re-175":"f","Re-176":"f","Re-177":"f",
	"Re-178":"f","Re-179":"f","Re-180":"f","Re-181":"f","Re-182":"f",
	"Re-183":"f","Re-184":"f","Re-185":"t","Re-186":"l","Re-187":"t",
	"Re-188":"f","Re-189":"f","Re-190":"f","Re-191":"f","Os-179":"f",
	"Os-180":"f","Os-181":"f","Os-182":"f","Os-183":"f","Os-184":"t",
	"Os-185":"f","Os-186":"t","Os-187":"t","Os-188":"t","Os-189":"t",
	"Os-190":"t","Os-191":"f","Os-192":"t","Os-193":"f","Os-194":"f",
	"Os-195":"f","Os-196":"f","Ir-181":"f","Ir-182":"f","Ir-183":"f",
	"Ir-184":"f","Ir-185":"f","Ir-186":"f","Ir-187":"f","Ir-188":"f",
	"Ir-189":"f","Ir-190":"f","Ir-191":"t","Ir-192":"f","Ir-193":"t",
	"Ir-194":"f","Ir-195":"f","Ir-196":"f","Ir-197":"f","Pt-184":"f",
	"Pt-185":"f","Pt-186":"f","Pt-187":"f","Pt-188":"f","Pt-189":"f",
	"Pt-190":"t","Pt-191":"f","Pt-192":"t","Pt-193":"f","Pt-194":"t",
	"Pt-195":"t","Pt-196":"t","Pt-197":"f","Pt-198":"t","Pt-199":"f",
	"Pt-200":"f","Pt-201":"f","Pt-202":"f","Au-185":"f","Au-186":"f",
	"Au-187":"f","Au-188":"f","Au-189":"f","Au-190":"f","Au-191":"f",
	"Au-192":"f","Au-193":"f","Au-194":"f","Au-195":"f","Au-196":"f",
	"Au-197":"t","Au-198":"f","Au-199":"f","Au-200":"f","Au-201":"f",
	"Au-202":"f","Au-203":"f","Hg-189":"f","Hg-190":"f","Hg-191":"f",
	"Hg-192":"f","Hg-193":"f","Hg-194":"f","Hg-195":"f","Hg-196":"t",
	"Hg-197":"f","Hg-198":"t","Hg-199":"t","Hg-200":"t","Hg-201":"t",
	"Hg-202":"t","Hg-203":"f","Hg-204":"t","Hg-205":"f","Hg-206":"f",
	"Hg-207":"f","Hg-208":"f","Tl-192":"f","Tl-193":"f","Tl-194":"f",
	"Tl-195":"f","Tl-196":"f","Tl-197":"f","Tl-198":"f","Tl-199":"f",
	"Tl-200":"f","Tl-201":"f","Tl-202":"f","Tl-203":"t","Tl-204":"l",
	"Tl-205":"t","Tl-206":"f","Tl-207":"f","Tl-208":"f","Tl-209":"f",
	"Tl-210":"f","Pb-193":"f","Pb-194":"f","Pb-195":"f","Pb-196":"f",
	"Pb-197":"f","Pb-198":"f","Pb-199":"f","Pb-200":"f","Pb-201":"f",
	"Pb-202":"f","Pb-203":"f","Pb-204":"t","Pb-205":"l","Pb-208":"t",
	"Pb-209":"f","Pb-210":"f","Pb-211":"f","Bi-202":"f","Bi-203":"f",
	"Bi-204":"f","Bi-205":"f","Bi-206":"f","Bi-207":"f","Bi-208":"f",
	"Bi-209":"t","Bi-210":"f","Al-*26":"f","Kr-*85":"f","Cd-*15":"f"}        

	##input is isotope
	if len(isotope.split("-"))==2:
		isotope="".join(isotope.split())
		for i in range(len(isotopes.keys())):
                	#print "|"+isotope+stable_list[i][0]+"|"
			if isotope == isotopes.keys()[i]:
				return isotopes[isotope]
	#if element
	if len(isotope.split("-"))==1:
		for i in range(len(isotopes)):
			if isotope == isotopes.keys()[i].split("-")[0]  and isotopes[isotopes.keys[i]] == "t":
				return "t"
