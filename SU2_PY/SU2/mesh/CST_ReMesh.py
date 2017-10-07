# CST generator
import sys,os
sys.path.append(os.environ['SU2_RUN'])
import SU2 # import all the python scripts in /usr/local/SU2_RUN/SU2_PY/SU2
import numpy as np
from math import factorial as fac
import matplotlib.pyplot as plt
from scipy.optimize import fmin_slsqp

class CST_ReMesh(object):

	def __init__(self,Config):

		# Order determined by the number of design variables supplied
		# Note that it's assumed that the order is identical for both surfaces
		self.Config=Config
		self.Order=int(0.5*len(Config['DEFINITION_DV']['PARAM'])-1)
		self.Marker=Config['DEFINITION_DV']['MARKER'][0][0]
		self.Mesh=Config['MESH_FILENAME']
		self.Mesh_out=self.Mesh[:-4]+'_CST'+self.Mesh[-4:]

		# TODO
		# Why are the x components different in the CST and 
		# the original foil for the RAE?


	def Update_Params(self,Au,Al):
		dvs=[0.0]*((len(Au)+len(Al)))
		#dvs=[0.2,0.0,0.2,0.0,0.0,0.0,0.0,0.0]
		self.Config.unpack_dvs(dvs)

		j=0
		k=0
		for i in range(len(Au)*2):
			if self.Config['DEFINITION_DV']['PARAM'][i][0]==1:
				self.Config['DEFINITION_DV']['PARAM'][i][1]=Au[j]
				self.Config['DEFINITION_DV']['KIND'][i]="CST"
				j+=1
			else:
				self.Config['DEFINITION_DV']['PARAM'][i][1]=Al[k]
				self.Config['DEFINITION_DV']['KIND'][i]="CST"
				k+=1

		# Change Mesh out filename
		self.Config['MESH_OUT_FILENAME']=self.Mesh_out

		SU2.io.config.dump_config('temp_config.cfg',self.Config)
		return

	def Update_Config(self):#
		# Use the new mesh
		self.Config['MESH_FILENAME']=self.Mesh_out
		# Reset the Meshout
		self.Config['MESH_OUT_FILENAME']='mesh_out.su2'
		return self.Config


	def Read_Mesh(self):

		Mesh=self.Mesh
		Marker=self.Marker
	
		# Using the Su2 python scripts for reading mesh
		Meshdata=SU2.mesh.tools.read(Mesh) # read the mesh

		# sort airfoil coords to be arrange clockwise from trailing edge
		Points,Loop=SU2.mesh.tools.sort_airfoil(Meshdata,Marker)
		
		# get the points for the surface marker
		Foil_Points,Foil_Nodes=SU2.mesh.tools.get_markerPoints(Meshdata,Marker)

		#Get the sorted points 
		Coords=np.zeros([len(Points),2])
		for i in range(len(Points)):
			Coords[i][0]=Foil_Points[Loop[i]][0]
			Coords[i][1]=Foil_Points[Loop[i]][1]
		# Divide coords for surfaces
		U_Coords,L_Coords=self.Split(Coords)

		return U_Coords,L_Coords

	def Fit(self,U_Coords,L_Coords):
	# initial coefficents set for upper (u) and lower (l) surfaces
		Au=np.ones(self.Order+1)# one more than the order
		Al=np.ones(self.Order+1)*-1 
		# Upper -  note that a seperator is needed after args var for an unknown reason
		Au=fmin_slsqp(self.Get_L2,Au,args=(U_Coords,),iprint=0)
		# Lower
		Al=fmin_slsqp(self.Get_L2,Al,args=(L_Coords,),iprint=0)

		return Au,Al 

	def Get_L2(self,A,Coords): 

		CST_Coords=self.CST(Coords,A)
		# Calculate the current L2 norm 
		L2=np.linalg.norm(CST_Coords - Coords,ord=2)

		return L2



	def Bi_Coeff(self,Order): 
		#compute the binomial coefficient
		K=np.zeros(Order+1)
		for i in range(len(K)):
			K[i]=fac(Order)/(fac(i)*(fac(Order-i)))
		return K


	def C_n1n2(self,Coords): 
		# class function
		n1=0.5
		n2=1.0
		C=np.zeros(len(Coords))
		for i in range(len(C)):
			C[i]=(Coords[i][0]**n1)*(1-Coords[i][0]**n2)
		return C

	def Total_Shape(self,Coords,A): 
		# Total shape function
		S=np.zeros(len(Coords))
		# Component Shape Function
		S_c=self.Comp_Shape(Coords)

		S_c=np.transpose(S_c)
		for  i in range(len(Coords)):
			S[i]+=np.dot(A,S_c[i])

		return S

	def Comp_Shape(self,Coords):
		Order=self.Order
		# Component Shape function
		K=self.Bi_Coeff(Order)
		# compute the Binomial Coefficient
		S_c=np.zeros([Order+1,len(Coords)])

		for i in range(Order+1): # order loop
			for j in range(len(Coords)): # point loop
				S_c[i][j]=(K[i]*(Coords[j][0]**i))*((1-Coords[j][0])**(Order-i))
		
		return S_c

	def CST(self,Coords,A): 
		CST_Coords=np.zeros([len(Coords),2])
		# Compute Class Function
		C=self.C_n1n2(Coords)

		# Compute the Shape Function
		S=self.Total_Shape(Coords,A)
		# evaluate the CST function
		for i in range(len(Coords)):
			CST_Coords[i][1]=C[i]*S[i]
			CST_Coords[i][0]=Coords[i][0]

		return CST_Coords


	def Split(self,Coords):
			# Spilt the surfaces according to the z component of the normal

		U_Coords=[]
		L_Coords=[]
		Normals=self.Get_Normal(Coords)

		for i in range(len(Coords)):
			if Normals[i][1]<0:
				L_Coords.append(Coords[i])
			else:
				U_Coords.append(Coords[i])

		# Convert to numpy array
		L_Coords=np.array(L_Coords)
		U_Coords=np.array(U_Coords)

		return U_Coords,L_Coords

	def Get_Normal(self,Coords):
		# Compute the normals

		# TODO 
		# Clean this up

		Normals=np.zeros([len(Coords),2])
		for i in range(len(Coords)):
			if i==0:
				dx_1=Coords[i][0]-Coords[len(Coords)-1][0]
				dy_1=Coords[i][1]-Coords[len(Coords)-1][1]
				dx_2=Coords[i+1][0]-Coords[i][0]
				dy_2=Coords[i+1][1]-Coords[i][1]
			
			elif i==len(Coords)-1:
				dx_1=Coords[i][0]-Coords[i-1][0]
				dy_1=Coords[i][1]-Coords[i-1][1]
				dx_2=Coords[0][0]-Coords[i][0]
				dy_2=Coords[0][1]-Coords[i][1]

			else:
				dx_1=Coords[i][0]-Coords[i-1][0]
				dy_1=Coords[i][1]-Coords[i-1][1]
				dx_2=Coords[i+1][0]-Coords[i][0]
				dy_2=Coords[i+1][1]-Coords[i][1]
		
			norm_1=-dy_1,dx_1
			norm_2=-dy_2,dx_2

			Normals[i][0]=0.5*(norm_1[0]+norm_2[0])
			Normals[i][1]=0.5*(norm_1[1]+norm_2[1])

		return Normals

	def Re_Mesh(self,Al,Au):
		# Update the Config and dump a raw temp copy for SU2_DEF
		self.Update_Params(Al,Au) 
		# Call Mesh Deformation Code
		os.system("SU2_DEF "+"temp_config.cfg")
		# remove the temp file
		os.system("rm temp_config.cfg")

		return 



	def Plot(self,U_Coords,L_Coords,Au,Al):

		u_coords=np.transpose(U_Coords)
		l_coords=np.transpose(L_Coords)
		
		# Temp Code
		Au[0]=Au[0]+0.3

		CST_Upper=CST(U_Coords,Au)
		CST_Lower=CST(L_Coords,Al)

		cst_upper=np.transpose(CST_Upper)
		cst_lower=np.transpose(CST_Lower)


		fig=plt.figure()
		ax=fig.add_subplot(1,1,1)

		
		plt.plot(cst_upper[0],cst_upper[1],'o',label='CST',color='blue',markersize=5)
		plt.plot(cst_lower[0],cst_lower[1],'o',color='blue',markersize=5)
		

		plt.plot(u_coords[0],u_coords[1],label='Baseline',color='green')
		plt.plot(l_coords[0],l_coords[1],color='green')

		plt.grid()
		ax.legend(loc='best')
		ax.set_xlabel('x/c')
		ax.set_ylabel('z/c')
		plt.title('Curve Fitting')
		
		filename='./Foils_Plot.png'
		plt.savefig(filename,dpi=150)
		plt.close()

		# temp - print out coords to compare with DEF mesh 
		Coords=np.append(CST_Upper,CST_Lower,axis=0)
		np.savetxt("CST_Coords.dat",Coords)

		return



