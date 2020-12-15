from __future__ import print_function
import csv
import numpy as np
from copy import deepcopy as copy
#functions including
#1. simple classes to store all the information about each species and reaction.
#2. Functions to read in the species and reaction file and check for sanity
#3. Functions to write out files necessary for UCLCHEM


##########################################################################################
#1. simple classes to store all the information about each species and reaction.
#largely just to make the other functions more readable.
##########################################################################################
reaction_types=['PHOTON','CRP','CRPHOT','FREEZE','THERM','DESOH2','DESCR','DEUVCR',"H2FORM","ER","ERDES","LH","LHDES"]
#these reaction types removed as UCLCHEM does not handle them. 'CRH','PHOTD','XRAY','XRSEC','XRLYA','XRPHOT'
elementList=['H','D','HE','C','N','O','F','P','S','CL','LI','NA','MG','SI','PAH','15N','13C','18O']
elementMass=[1,2,4,12,14,16,19,31,32,35,3,23,24,28,420,15,13,18]
symbols=['#','+','-','(',')']

class Species:
	def __init__(self,inputRow):
		self.name=inputRow[0]
		self.mass=inputRow[1]
		self.bindener=float(inputRow[2])
		self.solidFraction=float(inputRow[3])
		self.monoFraction=float(inputRow[4])
		self.volcFraction=float(inputRow[5])
		self.enthalpy=float(inputRow[6])

	def is_grain_species(self):
		return self.name[0]=='#'

	def is_ion(self):
		return (self.name[-1]=="+" or self.name[-1]=="-")

	def find_constituents(self):
		speciesName=self.name[:]
		i=0
		atoms=[]
		bracket=False
		bracketContent=[]
		#loop over characters in species name to work out what it is made of
		while i<len(speciesName):
			#if character isn't a #,+ or - then check it otherwise move on
			if speciesName[i] not in symbols:
				if i+1<len(speciesName):
					#if next two characters are (eg) 'MG' then atom is Mg not M and G
					if speciesName[i:i+3] in elementList:
						j=i+3
					elif speciesName[i:i+2] in elementList:
						j=i+2
					#otherwise work out which element it is
					elif speciesName[i] in elementList:
						j=i+1

				#if there aren't two characters left just try next one
				elif speciesName[i] in elementList:
					j=i+1
				#if we've found a new element check for numbers otherwise print error
				if j>i:
					if bracket:
						bracketContent.append(speciesName[i:j])
					else:
						atoms.append(speciesName[i:j])#add element to list
					if j<len(speciesName):
						if is_number(speciesName[j]):
							if int(speciesName[j])>1:
								for k in range(1,int(speciesName[j])):
									if bracket:
										bracketContent.append(speciesName[i:j])
									else:
										atoms.append(speciesName[i:j])
								i=j+1
							else:
								i=j
						else:
							i=j
					else:
						i=j
				else:
					print(speciesName[i])
					print("\t{0} contains elements not in element list:".format(speciesName))
					print(elementList)
			else:
				#if symbol is start of a bracketed part of molecule, keep track
				if (speciesName[i]=="("):
					bracket=True
					bracketContent=[]
					i+=1
				#if it's the end then add bracket contents to list
				elif speciesName[i]==")":
					if is_number(speciesName[i+1]):
						for k in range(0,int(speciesName[i+1])):
							atoms.extend(bracketContent)
						i+=2
					else:
						atoms.extend(bracketContent)
						i+=1
				#otherwise move on
				else:
					i+=1

		self.n_atoms=len(atoms)
		mass=0
		for atom in atoms:
			mass+=elementMass[elementList.index(atom)]
		if mass!=float(self.mass):
			print(f"Input mass of {self.name} does not match calculated mass of constituents")
			print("using calculated mass")
			self.mass=str(mass)


class Reaction:
	def __init__(self,inputRow):
		self.reactants=[inputRow[0],inputRow[1],self.NANCheck(inputRow[2])]
		self.products=[inputRow[3],self.NANCheck(inputRow[4]),self.NANCheck(inputRow[5]),self.NANCheck(inputRow[6])]
		self.alpha=float(inputRow[7])
		self.beta=float(inputRow[8])
		self.gamma=float(inputRow[9])
		self.templow=float(inputRow[10])
		self.temphigh=float(inputRow[11])
		self.reac_type=self.get_reaction_type()
		self.duplicate=False

		#body_count is the number of factors of density to include in ODE
		#we drop a factor of density from both the LHS and RHS of ODES
		#So reactions with 1 body have no factors of density which we manage by counting from -1
		self.body_count=-1
		for reactant in self.reactants:
			if (reactant not in reaction_types) and reactant!="NAN":
				self.body_count+=1
			if reactant in ["DESOH2","FREEZE"]:
				self.body_count+=1

	def NANCheck(self,a):
		aa  = a if a else 'NAN'
		return aa

	def get_reaction_type(self):
		if (self.reactants[2] in reaction_types):
			return self.reactants[2]
		else:
			if self.reactants[1] in reaction_types:
				return self.reactants[1]
			else:
				return "TWOBODY"

	def same_reaction(self,other):
		if set(self.reactants)==set(other.reactants):
			if set(self.products)==set(other.products):
				return True
		return False

	def print(self):
		print(" + ".join(self.reactants), "->","+".join(self.products))


##########################################################################################
#2. Functions to read in the species and reaction file and check for sanity
##########################################################################################

# Read the entries in the specified species file
def read_species_file(fileName):
	speciesList=[]
	f = open(fileName, 'r')
	reader = csv.reader(f, delimiter=',', quotechar='|')
	for row in reader:
		if row[0]!="NAME" and "!" not in row[0] :
			speciesList.append(Species(row))
	nSpecies = len(speciesList)
	return nSpecies,speciesList

# Read the entries in the specified reaction file and keep the reactions that involve the species in our species list
def read_reaction_file(fileName, speciesList, ftype):
	reactions=[]
	dropped_reactions=[]
	keepList=['','NAN','#','E-','e-','ELECTR']
	keepList.extend(reaction_types)

	for species in speciesList:
		keepList.append(species.name)			                                  
	if ftype == 'UMIST': # if it is a umist database file
		f = open(fileName, 'r')
		reader = csv.reader(f, delimiter=':', quotechar='|')
		for row in reader:
			if all(x in keepList for x in [row[2],row[3],row[4],row[5],row[6],row[7]]): #if all the reaction elements belong to the keeplist
				#umist file doesn't have third reactant so add space and has a note for how reactions there are so remove that
				reactions.append(Reaction(row[2:4]+['']+row[4:8]+row[9:]))
	if ftype == 'UCL':	# if it is a ucl made (grain?) reaction file
		f = open(fileName, 'r')
		reader = csv.reader(f, delimiter=',', quotechar='|')
		for row in reader:
			if all(x in keepList for x in row[0:7]):	#if all the reaction elements belong to the keeplist
				if row[10]=="":
					row[10]=0.0
					row[11]=10000.0
				reactions.append(Reaction(row))	
			else:
				dropped_reactions.append(row)

	nReactions = len(reactions)
	return nReactions, reactions, dropped_reactions

def remove_duplicate_species(speciesList):
	#check for duplicate species
	duplicates=0
	duplicate_list=[]
	for i in range(0,len(speciesList)):
		for j in range(0,len(speciesList)):
			if speciesList[i].name==speciesList[j].name:
				if (j!=i) and speciesList[i].name not in duplicate_list:
					print("\t {0} appears twice in input species list".format(speciesList[i].name))
					duplicate_list.append(speciesList[i].name)

	for duplicate in duplicate_list:
		removed=False
		i=0
		while not removed:
			if speciesList[i].name==duplicate:
				del speciesList[i]
				print("\tOne entry of {0} removed from list".format(duplicate))
				removed=True
			else:
				i+=1
	return speciesList

#Look for possibly incorrect parts of species list
def check_and_filter_species(speciesList,reactionList):
	#check for species not involved in any reactions
	lostSpecies=[]
	for species in speciesList:
		keepFlag=False
		for reaction in reactionList:
			if species.name in reaction.reactants or species.name in reaction.products:
				keepFlag=True
		if not keepFlag:
			lostSpecies.append(species.name)
			speciesList.remove(species)

	print('\tSpecies in input list that do not appear in final list:')
	print('\t',lostSpecies)
	print('\n')
	for species in speciesList:
		species.find_constituents()
	return speciesList

#All species should freeze out at least as themselves and all grain species should desorb according to their binding energy
#This function adds those reactions automatically to slim down the grain file
def add_desorb_reactions(speciesList,reactionList):
	desorb_reacs=['DESOH2',"DESCR","DEUVCR","THERM"]

	for species in speciesList:
		if species.is_grain_species():
			for reacType in desorb_reacs:
				newReaction=Reaction([species.name,reacType,'NAN',species.name[1:],'NAN','NAN','NAN',1,0,species.bindener,0.0,10000.0])
				reactionList.append(newReaction)
	return reactionList


def add_chemdes_reactions(speciesList,reactionList):
	new_reacs=[]
	for reaction in reactionList:
		if reaction.reac_type in ["LH","ER"]:
			new_reac=copy(reaction)
			new_reac.reac_type=new_reac.reac_type+"DES"
			new_reac.reactants[2]=new_reac.reactants[2]+"DES"
			for i,product in enumerate(new_reac.products):
				if ("#" in product):
					new_reac.products[i]=new_reac.products[i][1:]
				else:
					if product!="NAN":
						print("All Langmuir-Hinshelwood and Eley-Rideal reactions should be input with products on grains only.")
						print("The fraction of products that enter the gas is dealt with by Makerates and UCLCHEM.")
						print("the following reaction caused this warning")
						reaction.print()
			new_reacs.append(new_reac)

	reactionList=reactionList+new_reacs
	return reactionList

#check reactions to alert user of potential issues including repeat reactions
#and multiple freeze out routes
def reaction_check(speciesList,reactionList,freeze_check=True):


	#first check for multiple freeze outs so user knows to do alphas
	print("\tSpecies with multiple freeze outs, check alphas:")
	for spec in speciesList:
		freezes=0
		for reaction in reactionList:
			if (spec.name in reaction.reactants and 'FREEZE' in reaction.reactants):
				freezes+=1
		if (freezes>1):
			print("\t{0} freezes out through {1} routes".format(spec.name,freezes))
		if freezes<1 and not spec.is_grain_species() and freeze_check:
			print("\t{0} does not freeze out".format(spec.name,freezes))

	#now check for duplicate reactions
	duplicate_list=[]
	print("\n\tPossible duplicate reactions for manual removal:")
	duplicates=False
	for i, reaction1 in enumerate(reactionList):
		#if i not in duplicate_list:
			for j, reaction2 in enumerate(reactionList):
				if i!=j:
					if reaction1.same_reaction(reaction2):
						print("\tReactions {0} and {1} are possible duplicates".format(i+1,j+1))
						reaction1.print()
						reaction2.print()
						duplicates=True
						#adjust temperatures so temperature ranges are adjacent
						if reaction1.temphigh > reaction2.temphigh:
							if reaction1.templow<reaction2.temphigh:
								print(f"\tReactions {i+1} and {j+1} have non-adjacent temperature ranges")
						reaction1.duplicate=True
						reaction2.duplicate=True
	
	if (not duplicates):
		print("\tNone")

#capitalize files
def make_capitals(fileName):
	a=open(fileName).read()
	output = open(fileName, mode='w')
	output.write(a.upper())
	output.close()


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
##########################################################################################
#3. Functions to write out files necessary for UCLCHEM
##########################################################################################

# Write the species file in the desired format
def write_species(fileName, speciesList):
	f= open(fileName,'w')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')		
	nSpecies = len(speciesList)
	for species in speciesList:
		writer.writerow([species.name,species.mass,species.n_atoms])

# Write the reaction file in the desired format
def write_reactions(fileName, reactionList):
	f = open(fileName, 'w')
	writer = csv.writer(f,delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	nReactions = len(reactionList)
	for reaction in reactionList:
		#if statement changes beta for ion freeze out to 1. This is how ucl_chem recognises ions when calculating freeze out rate
		if ('FREEZE' in reaction.reactants and reaction.reactants[0][-1]=='+'):
			reaction.beta=1
		writer.writerow(reaction.reactants+reaction.products+[reaction.alpha,reaction.beta,reaction.gamma,reaction.templow,reaction.temphigh])

def write_odes_f90(fileName, speciesList, reactionList):
	output = open(fileName, mode='w')
	#go through every species and build two strings, one with eq for all destruction routes and one for all formation
	ydotString=build_ode_string(speciesList,reactionList)
	output.write(ydotString)
	output.close()    

def build_ode_string(speciesList, reactionList):
	species_names=[]
	for i, species in enumerate(speciesList):
		species_names.append(species.name)
		species.losses=""
		species.gains=""
	for i,reaction in enumerate(reactionList):
		ODE_BIT=f"+RATE({i+1})"
		
		#every body after the first requires a factor of density
		for body in range(reaction.body_count):
			ODE_BIT=ODE_BIT+"*D"
		
		#then bring in factors of abundances
		for species in reaction.reactants:
			if species in species_names:
				ODE_BIT=ODE_BIT+f"*Y({species_names.index(species)+1})"
			if "H2FORM" in reaction.reactants:
				#only 1 factor of H abundance in Cazaux & Tielens 2004 H2 formation so stop looping after first iteration
				break

		#now add to strings
		for species in reaction.reactants:
			if species in species_names:
				#Eley-Rideal reactions take a share of total freeze out rate which is already accounted for
				#so we add as a loss term to the frozen version of the species rather than the gas version
				if ("ELEYRIDEAL" in reaction.reactants) and (not speciesList[species_names.index(species)].is_grain_species()):
					speciesList[species_names.index("#"+species)].losses+=ODE_BIT
				else:
					speciesList[species_names.index(species)].losses+=ODE_BIT
		for species in reaction.products:
			if species in species_names:
				speciesList[species_names.index(species)].gains+=ODE_BIT

	ode_string=""
	for n,species in enumerate(speciesList):


		if species.losses != '':
			loss_string = '    LOSS = '+species.losses[1:]+'\n'
			loss_string = truncate_line(loss_string)
			ode_string+=loss_string
		if species.gains != '':
			prod_string = '    PROD = '+species.gains[1:]+'\n'
			prod_string = truncate_line(prod_string)
			ode_string+=prod_string
		
		#start with empty string and add production and loss terms if they exists
		ydot_string=''
		if prod_string != '':
			ydot_string += 'PROD'
			if loss_string != '':
				ydot_string += '-LOSS'

		#if we have prod and/or loss add ydotstring to odes
		if ydot_string=='':
			ydot_string="0.0"
		ydot_string=f"    YDOT({n+1}) = {ydot_string}\n"
		ydot_string = truncate_line(ydot_string)
		ode_string+=ydot_string
	return ode_string

#create a file containing length of each list of moleculetypes and then the two lists (gas and grain) of species in each type
#as  well as fraction that evaporated in each type of event
def write_evap_lists(openFile,speciesList):
	grainlist=[];mgrainlist=[];solidList=[];monoList=[];volcList=[]
	bindEnergyList=[];enthalpyList=[]

	for i,species in enumerate(speciesList):
		if species.name[0]=='#':
			#find gas phase version of grain species. For #CO it looks for first species in list with just CO and then finds the index of that
			try:
				j=speciesList.index(next((x for x in speciesList if x.name==species.name[1:]))) 
			except:
				print("\n**************************************\nWARNING\n**************************************")
				print("{0} has no gas phase equivalent in network. Every species should at least freeze out and desorb.".format(species.name))
				print("ensure {0} is in the species list, and at least one reaction involving it exists and try again".format(species.name[1:]))
				print("Alternatively, provide the name of the gas phase species you would like {0} to evaporate as".format(species.name))
				input=raw_input("type x to quit Makerates or any species name to continue\n")
				if input.lower()=="x":
					exit()
				else:
					j=speciesList.index(next((x for x in speciesList if x.name==input.upper())))					

			#plus ones as fortran and python label arrays differently
			mgrainlist.append(i+1)
			grainlist.append(j+1)
			solidList.append(species.solidFraction)
			monoList.append(species.monoFraction)
			volcList.append(species.volcFraction)
			bindEnergyList.append(species.bindener)
			enthalpyList.append(species.enthalpy)

	openFile.write(array_to_string("gasGrainList",grainlist,type="int"))
	openFile.write(array_to_string("grainList",mgrainlist,type="int"))
	openFile.write(array_to_string("solidFractions",solidList,type="float"))
	openFile.write(array_to_string("monoFractions",monoList,type="float"))
	openFile.write(array_to_string("volcanicFractions",volcList,type="float"))
	openFile.write(array_to_string("bindingEnergy",bindEnergyList,type="float",parameter=False))
	openFile.write(array_to_string("formationEnthalpy",enthalpyList,type="float"))

	return len(grainlist)

# Create the appropriate multiplication string for a given number
def multiple(number):
    if number == 1: return ''
    else: return str(number)+'*'

def truncate_line(input, codeFormat='F90', continuationCode=None):
	lineLength = 72
	maxlines=300
	lines=0
	result = ''
	i=0
	j=0
	while i+j<len(input):
		j+=1
		if j>lineLength:
			#important not to break entries so split lines at ,s
			try:
				k=input[i+j-16:i+j].index(",")
			except:
				try:
					k=input[i+j-16:i+j].index("*")
				except:
					k=input[i+j-16:i+j].index(")")
			j=j-16+k
			result+=input[i:i+j]+"&\n    &"
			i=i+j
			j=0
	result+=input[i:i+j]
	return result    


def write_network_file(fileName,speciesList,reactionList):
	openFile=open(fileName,"w")
	openFile.write("MODULE network\nUSE constants\nIMPLICIT NONE\n")
	openFile.write("    INTEGER, PARAMETER :: nSpec={0}, nReac={1}\n".format(len(speciesList),len(reactionList)))

	#write arrays of all species stuff
	names=[]
	atoms=[]
	masses=[]
	for species in speciesList:
		names.append(species.name)
		masses.append(float(species.mass))
		atoms.append(species.n_atoms)

	speciesIndices=""
	for element in ["E-","C+","H+","H2","SI+","S+","CL+","CO","HE+","#H","#H2","#N","#O",'#OH']+elementList:
		try:
			species_index=names.index(element)+1
	
		except:
			print(element," not in network, adding dummy index")
			species_index=len(speciesList)
		name=element.lower().replace("+","x").replace("e-","elec").replace("#","g")
		speciesIndices+="n{0}={1},".format(name,species_index)
	if len(speciesIndices)>72:
		speciesIndices=truncate_line(speciesIndices)
	speciesIndices=speciesIndices[:-1]+"\n"
	openFile.write("    INTEGER, PARAMETER ::"+speciesIndices)
	openFile.write(array_to_string("    specname",names,type="string"))
	openFile.write(array_to_string("    mass",masses,type="float"))
	openFile.write(array_to_string("    atomCounts",atoms,type="int"))


	#then write evaporation stuff
	nGrain=write_evap_lists(openFile,speciesList)

	#finally all reactions
	reactant1=[]
	reactant2=[]
	reactant3=[]
	prod1=[]
	prod2=[]
	prod3=[]
	prod4=[]
	alpha=[]
	beta=[]
	gama=[]
	reacTypes=[]
	duplicates=[]
	tmins=[]
	tmaxs=[]
	#store important reactions
	reactionIndices=""

	for i,reaction in enumerate(reactionList):
		if ("CO" in reaction.reactants) and ("PHOTON" in reaction.reactants):
			if "O" in reaction.products and "C" in reaction.products:
				reactionIndices+="nR_CO_hv={0},".format(i+1)
		if ("C" in reaction.reactants) and ("PHOTON" in reaction.reactants):
			reactionIndices+="nR_C_hv={0},".format(i+1)
		if ("H2FORM" in reaction.reactants):
			reactionIndices+=f"nR_H2Form_CT={i+1},"
		if (("H" in reaction.reactants) and ("#H" in reaction.reactants)):
			if "H2" in reaction.products:
				reactionIndices+=f"nR_H2Form_ERDes={i+1},"
			elif "#H2" in reaction.products:
				reactionIndices+=f"nR_H2Form_ER={i+1},"
		if ((reaction.reactants.count("#H")==2) and ("LH" in reaction.reactants)):
			reactionIndices+=f"nR_H2Form_LH={i+1},"
		if ((reaction.reactants.count("#H")==2) and ("LHDES" in reaction.reactants)):
			reactionIndices+=f"nR_H2Form_LHDes={i+1},"
		if (("H" in reaction.reactants) and ("FREEZE" in reaction.reactants)):
			reactionIndices+=f"nR_HFreeze={i+1},"
		if (("E-" in reaction.reactants) and ("FREEZE" in reaction.reactants)):
			reactionIndices+=f"nR_EFreeze={i+1},"
	reactionIndices=reactionIndices[:-1]


	if len(reactionIndices)>60:
		reactionIndices=reactionIndices[:60]+"&\n&"+reactionIndices[60:]
	reactionIndices=reactionIndices[:-1]+"\n"
	openFile.write("    INTEGER, PARAMETER ::"+reactionIndices)

	for i,reaction in enumerate(reactionList):
		if "CO" in reaction.reactants and "PHOTON" in reaction.reactants:
			if "O" in reaction.products and "C" in reaction.products:
				openFile.write(f"INTEGER, PARAMETER :: nrco={i+1}\n")
		reactant1.append(find_reactant(names,reaction.reactants[0]))
		reactant2.append(find_reactant(names,reaction.reactants[1]))
		reactant3.append(find_reactant(names,reaction.reactants[2]))
		prod1.append(find_reactant(names,reaction.products[0]))
		prod2.append(find_reactant(names,reaction.products[1]))
		prod3.append(find_reactant(names,reaction.products[2]))
		prod4.append(find_reactant(names,reaction.products[3]))
		alpha.append(reaction.alpha)
		beta.append(reaction.beta)
		gama.append(reaction.gamma)
		if reaction.duplicate:
			duplicates.append(i+1)
			tmaxs.append(reaction.temphigh)
			tmins.append(reaction.templow)
		reacTypes.append(reaction.reac_type)
	if len(duplicates)==0:
		duplicates=[9999]
		tmaxs=[0]
		tmins=[0]

	openFile.write(array_to_string("\tre1",reactant1,type="int"))
	openFile.write(array_to_string("\tre2",reactant2,type="int"))
	openFile.write(array_to_string("\tre3",reactant3,type="int"))
	openFile.write(array_to_string("\tp1",prod1,type="int"))
	openFile.write(array_to_string("\tp2",prod2,type="int"))
	openFile.write(array_to_string("\tp3",prod3,type="int"))
	openFile.write(array_to_string("\tp4",prod4,type="int"))
	openFile.write(array_to_string("\talpha",alpha,type="float",parameter=False))
	openFile.write(array_to_string("\tbeta",beta,type="float",parameter=False))
	openFile.write(array_to_string("\tgama",gama,type="float",parameter=False))
	openFile.write(array_to_string("\tduplicates",duplicates,type="int",parameter=True))
	openFile.write(array_to_string("\tminTemps",tmins,type="float",parameter=True))
	openFile.write(array_to_string("\tmaxTemps",tmaxs,type="float",parameter=True))
	
	reacTypes=np.asarray(reacTypes)
	for reaction_type in reaction_types+["TWOBODY"]:
		list_name=reaction_type.lower()+"Reacs"
		indices=np.where(reacTypes==reaction_type)[0]
		if len(indices>1):
			indices=[indices[0]+1,indices[-1]+1]
			openFile.write(array_to_string("\t"+list_name,indices,type="int",parameter=True))
	openFile.write("END MODULE network")
	openFile.close()

def find_reactant(species_list,reactant):
	try:
		return species_list.index(reactant)+1
	except:
		return 9999



def array_to_string(name,array,type="int",parameter=True):
	if parameter:
		outString=", PARAMETER :: "+name+" ({0})=(/".format(len(array))
	else:
		outString=" :: "+name+" ({0})=(/".format(len(array))
	if type=="int":
		outString="INTEGER"+outString
		for value in array:
			outString+="{0},".format(value)
	elif type=="float":
		outString="REAL(dp)"+outString
		for value in array:
			outString+="{0:.4e},".format(value)
	elif type=="string":
		strLength=len(max(array, key=len))
		outString="CHARACTER(Len={0:.0f})".format(strLength)+outString
		for value in array:
			outString+="\""+value.ljust(strLength)+"\","
	else:
		print("Not a valid type for array to string")
	outString=outString[:-1]+"/)\n"
	outString=truncate_line(outString)
	return outString
