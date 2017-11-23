# Tools for the Vn Process used for CAD optimisation
import sys,os
sys.path.append(os.environ['SU2_RUN'])
import SU2
import numpy as np



class Vn_Tools(object):

	def __init__(self,config):
		self.config=config
		self.marker=config['DEFINITION_DV']['MARKER'][0][0]
		self.mesh=config['MESH_FILENAME']
		
	def Mesh_Out(self):
		mesh=self.mesh
		marker=self.marker
		# read entire mesh 
		meshdata=SU2.mesh.tools.read(mesh)
		# get marker points
		marker_points,marker_nodes=SU2.mesh.tools.get_markerPoints(meshdata,marker)
		# Points don't need sorted as Vn_program uses kd-tree
		np.savetxt(open("/home/phewitt/Dropbox/Opt_Sync/Initial_Coords.txt",'w'),marker_points)
		
		return

	def Make_Pickle(self,Vn_data):
		print Vn_data['wantdisp']





		return