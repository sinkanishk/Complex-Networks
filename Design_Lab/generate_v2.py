import os, sys, time
import subprocess

#new_dir = sys.argv[1]
#dir_lfr = sys.argv[2]

#new_dir = "./netsForcomparingBaseline/"
new_dir = "./newnetsForcomparingBaseline/"
#dir_lfr = "./netsForcomparingBaseline/dir_lfr/"
dir_lfr = "./newnetsForcomparingBaseline/dir_lfr/"


if not os.path.exists(new_dir):
	os.makedirs(new_dir)
if not os.path.exists(new_dir + "_networks"):
	os.makedirs(new_dir + "_networks")


#n=1000
#k=11
#maxk=30
n=100		#number of nodes in each layer
k=6     	#Average degree	
maxk=10		#Max degree
mu =0.05

# old experiments (lfr_multilayer_v1)
	# list_alpha = [0.4, 0.6, 0.8, 1.0]
	# list_p = [0.1,0.25,0.4,0.6,0.8]
	# list_density = [0.004, 0.01, 0.025, 0.04, 0.055, 0.07]
	# list_mu=[0.05, 0.20,0.75,0.55, 0.4, 0.6]

	# list_alpha = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
	# list_mu = [0.05, 0.20, 0.4, 0.55, 0.6, 0.75]
	# list_p = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
	# list_p1 = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]
	# list_p2 = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]

	# list_alpha = [0.8]
	# list_mu = [0.05]#, 0.20, 0.4, 0.55, 0.6, 0.75]
	# list_p = [0.7]
	# #list_density = [0.9]#, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
	# list_p1 = [0.1]#, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]
	# list_p2 = [0.1]# 0.3, 0.5, 0.7, 0.8, 0.9, 1.0]

# config1
# list_alpha = [0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9]
# #list_alpha = [0.1, 0.2]
# list_mu = [0.05]
# #list_p = [1.0,0.95,0.9,0.85,0.8,0.75,0.7]
# list_p = [0.6,0.7,0.8,0.9]
# list_p1 = [0.7]
# list_p2 = [0.0]

#commented code
	# config2
	# list_alpha = [0.6]
	# list_mu = [0.05]
	# list_p = [0.6]
	# list_p1 = [0.1, 0.3, 0.6, 0.8]
	# list_p2 = [0.3]

	# config3
	# list_alpha = [0.6]
	# list_mu = [0.05]
	# list_p = [0.6]
	# list_p1 = [0.7]
	# list_p2 = [0.1, 0.3, 0.6, 0.8]

	 
	# config4
	# list_alpha = [0.6]
	# list_mu = [0.05]
	# list_p = [0.2,0.4,0.6,0.8]
	# list_p1 = [0.7]
	#list_p2 = [0.3]

	# config5
	# list_alpha = [0.6]
	# list_mu = [0.05, 0.2, 0.4, 0.55]
	# list_p = [0.6]
	# list_p1 = [0.7]
	# list_p2 = [0.3]

#for new lets just take one value for each parameter 
# list_alpha = [0.2,0.4,0.6,0.8]
# list_mu = [0.05,0.1,0.2,0.3]
# list_p = [0.2,0.4,0.6,0.8]
# list_p1 = [0.2,0.4,0.6,0.8]
# list_p2 = [0.0,0.1,0.2,0.3]

list_alpha = [0.4,0.6,0.8]
list_mu = [0.05]
list_p = [0.2,0.6,0.8]
list_p1 = [0.2,0.8]
list_p2 = [0.0]

def calculateDm():
	dm=0
	filepath="./netsForDtDmDb/_networks/network_0.6_0.7_0.05_1.0_0.0"
	f=open(filepath,'r')
	f.readline()
	f.readline()
	for j in range(0,3):
		line = f.readline()
		line = line.strip()
		itr = int(float(line))
		if(j==2):
			return itr/100.0
		for i in range(0,itr):
			line=f.readline()
			line=line.rstrip()
			line=line.split()
			n1=int(float(line[0]))
			n2=int(float(line[1]) )
			#print(n1,n2)
			if ((n1<=100 and n2>100) or (n1>100 and n2<=100)):
				dm+=1
				break
		if(j==0):
			f.readline()
		elif(j==1):
			f.readline()
			f.readline()
	dm =dm/100.0
	return dm  

def _plot(avgdegree):
	#maxk = avgdegree +5
	i=0
	nb_iteration = len(list_alpha) * len(list_mu) * len(list_p) * len(list_p1) * len(list_p2)

	for alpha in list_alpha:
		new_dir_alpha = new_dir + "/alpha-" + str(alpha)
		if not os.path.exists(new_dir_alpha):
			os.makedirs(new_dir_alpha)

		for p in list_p:
			new_dir_p = new_dir_alpha + "/p-" + str(p)
			if not os.path.exists(new_dir_p):
				os.makedirs(new_dir_p)
				
			for mu in list_mu:
				new_dir_mu = new_dir_p + "/mu-" + str(mu)
				if not os.path.exists(new_dir_mu):
					os.makedirs(new_dir_mu)	
			
				for p1 in list_p1:
					new_dir_p1 = new_dir_mu + "/p1-" + str(p1)
					if not os.path.exists(new_dir_p1):
						os.makedirs(new_dir_p1)
					
					for p2 in list_p2:
						new_dir_p2 = new_dir_p1 + "/p2-" + str(p2)
						if not os.path.exists(new_dir_p2):
							os.makedirs(new_dir_p2)

						i += 1
						print "%i/%i" % (i, nb_iteration)

						#new_dir_p2 = new_dir +"testingnetgeneration"
						#dir_lfr = new_dir +"testingnetgeneration/"
						#cmd = "python2.7 lfr_multilayer_v3.py %s %s -p %f -a %f -p1 %f -p2 %f" % (new_dir_p2, dir_lfr, p, alpha, p1, p2)
						#cmd = "python2.7 lfr_multilayer_v3.py %s %s -p %f -a %f -p1 %f -p2 %f -n %f -k %f -maxk %f -mu %f" % (new_dir_p2, dir_lfr, p, alpha, p1, p2, n, avgdegree, maxk, mu)
						cmd = "python2.7 temp_lfr_multilayer_v3.py %s %s -p %f -a %f -p1 %f -p2 %f -n %f -k %f -maxk %f -mu %f" % (new_dir_p2, dir_lfr, p, alpha, p1, p2, n, avgdegree, maxk, mu)
						print cmd
						#sys.exit(1)
						process = os.popen(cmd)
						preprocessed = process.read()
						process.close()
						
						if os.path.isfile(new_dir_p2 + "/new_format"):
							cmd = "cp " + new_dir_p2 + "/new_format " + new_dir + "_networks/network_" + str(alpha)  + "_" + str(p) + "_" + str(mu) + "_" + str(p1) + "_" + str(p2)+ "_" + str(avgdegree)
			 				os.popen(cmd)
							print "File copied ----------"
						time.sleep(0.15)
						
	


def plot(dm):
	jumps=1
	for i in range(3,100,jumps):
		if(i==10):
			jumps = i
		
		dt = (i)*dm
		avgdegree=round(dt)
		_plot(avgdegree)

_plot(k)
#dm = calculateDm()
#print("dm: ",round(dm))
#plot(round(dm))


