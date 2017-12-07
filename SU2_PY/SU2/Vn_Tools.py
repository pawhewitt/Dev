# Tools for the Vn Process used for CAD optimisation
import sys,os
sys.path.append(os.environ['SU2_RUN'])
import SU2
import numpy as np
import pickle 



class Vn_Tools(object):

	def __init__(self,config):
		self.config=config
		self.marker=config['DEFINITION_DV']['MARKER'][0][0]
		self.mesh=config['MESH_FILENAME']
		self.home=os.getenv("HOME")
		
	# def Mesh_Out(self):
	# 	mesh=self.mesh
	# 	marker=self.marker
	# 	home=self.home
	# 	# read entire mesh 
	# 	meshdata=SU2.mesh.tools.read(mesh)
	# 	# get marker points
	# 	marker_points,marker_nodes=SU2.mesh.tools.get_markerPoints(meshdata,marker)
	# 	# Points don't need sorted as Vn_program uses kd-tree
		
	
	# 	np.savetxt(open(home+"/Dropbox/Opt_Sync/Initial_Design.txt",'w'),marker_points)
		
	# 	return

	def Make_Pickle(self,Vn_data):
		home=self.home
		print "Exporting Design Vn_Data"
		pickle.dump(Vn_data,open(home+'/Dropbox/Opt_Sync/Vn_Data.pkl','w'))

		return

	def Remove_Pickle(self):
		home=self.home
		print "Removing Vn_Data.pkl--delete"
		os.remove(home+'/Dropbox/Opt_Sync/Vn_Data.pkl--delete')
		
		return

	def Remove_Results(self,module):
		home=self.home
		print "Removing Results file"
		if module =="CFD":
			os.remove(home+'/Dropbox/Opt_Sync/Sens.txt')
		if module =="DEF":
			os.remove(home+'/Dropbox/Opt_Sync/Disp.txt')

		return