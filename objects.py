#!/usr/bin/env python

import sys

class Interaction(object):
    """docstring for Interaction"""
    def __init__(self, chain1=None, res1=None, atom1=None, chain2=None, \
                 res2=None, atom2=None, distance=None):
        super(Interaction, self).__init__()

        self.chain1 = chain1
        self.chain2 = chain2
        self.res1 = res1
        self.res2 = res2
        self.atom1 = atom1
        self.atom2 = atom2
        self.distance = distance

        self.id = '{0}:{1}-{2}:{3}'.format(chain1,res1,chain2,res2)


    def load_from_PISA_line(self,PISA_line):
        # PISA line is as shown in the PISA server results page
        # 
        # example line:- H:GLN  40[ NE2]     2.83    L:GLN  40[ OE1]
        # 
        # The above example line lists a hydrogen bond between NE2 atom of 
        # 40_GLN of chain H and OE1 atom of 40_GLN of chain L at 2.83 distance

        # split the line first
        items = PISA_line.split('\t')

        # ideally if the line follows the example format, this can be done
        chain1, rest = items[0].split(':')
        res1, atom1 = rest.split('[')
        # there will be another ] , remove it
        atom1 = atom1.strip()[:-1]
        res1 = res1.strip()
        res1, resnum1 = res1.split()
        res1 = res1+'_'+resnum1
        
        ## repeat for chain2
        chain2, rest = items[-1].split(':')
        res2, atom2 = rest.split('[')
        # there will be another ] , remove it
        atom2 = atom2.strip()[:-1]
        res2 = res2.strip()
        res2, resnum2 = res2.split()
        res2 = res2+'_'+resnum2
        
        # distance 
        distance = float(items[1])
        
        self.__init__(chain1, res1, atom1, chain2, res2, atom2, distance)

        return True
        
    def load_from_PISA_xml_bond(self, xml_bond):
        ''' 
        xml_bond is an Element object which is a child of ElementTree object
        of the xml library
        type:- xml.etree.ElementTree.Element
        '''

        # first check if the type is matching
        if not (str(type(xml_bond))=="<class 'xml.etree.ElementTree.Element'>"\
                and xml_bond.tag=='bond' ):
            return False

        else:
            chain1 = xml_bond[0].text
            res1 = xml_bond[1].text+'_'+xml_bond[2].text
            atom1 = xml_bond[4].text

            chain2 = xml_bond[5].text
            res2 = xml_bond[6].text+'_'+xml_bond[7].text
            atom2 = xml_bond[9].text

            distance = float(xml_bond[-1].text)
            
            self.__init__(chain1, res1, atom1, chain2, res2, atom2, distance)

            return True

class Interaction_Set(object):
    """docstring for Interaction_Set"""
    def __init__(self, interaction_type, n_interactions=0):
        super(Interaction_Set, self).__init__()
        
        if interaction_type in ['h-bonds', 'salt-bridges']:
            self.interaction_type = interaction_type

        self.n_interactions = int(n_interactions)
        self.interactions = []

    def load_from_PISA_xml_bonds(self, xml_bonds):
        # xml_bonds is a set of xml.etree.ElementTree.Element type objects
        # see Interaction.load_from_PISA_xml_bond

        for bond in xml_bonds:
            interaction = Interaction()
            if interaction.load_from_PISA_xml_bond(bond):
                self.interactions.append(interaction)

    def calc_frequency_by_residue(self):
        ''' 
        there may be several residues making more than one interactions 
        frequency by residue will give us these numbers
        '''

        frequency = {}
        for interaction in self.interactions:
            side1, side2 = interaction.id.split('-')
            if side1 in frequency:
                frequency[side1]+=1
            else:
                frequency[side1]=1
            if side2 in frequency:
                frequency[side2]+=1
            else:
                frequency[side2]=1

        return frequency
    
    def __repr__(self):
        return 'Interaction set {0} with {1} interactions'\
            .format(self.interaction_type,self.n_interactions)
        

class Interface(object):
    """docstring for Interface"""
    def __init__(self, pdb_id, chain1=None, chain2=None, int_area=0, \
                 interaction_sets = []):
        
        super(Interface, self).__init__()
        self.pdb_id = pdb_id
        self.chain1 = chain1
        self.chain2 = chain2
        self.id = '{0}_{1}'.format(chain1,chain2)
        self.int_area = float(int_area)
        self.interaction_sets = interaction_sets
    
    def load_from_PISA_xml_interface(self, xml_interface):

        # first check if the type is matching
        if not (str(type(xml_interface))=="<class 'xml.etree.ElementTree.Element'>" \
                and xml_interface.tag=='interface' ):
            return False

        else:
            self.int_area = float(xml_interface[3].text)
            self.css = float(xml_interface[7].text)
            hbonds = xml_interface[11]
            salts = xml_interface[12]

            mol1 = xml_interface[15]
            mol2 = xml_interface[16]
            
            # See if it is a protein interface or a ligand is involved
            if mol1[2].text=='Ligand' or mol2[2].text=='Ligand':
                self.type = 'Ligand'
            else:
                self.type = 'Protein'

            self.chain1 = mol1[1].text
            self.chain2 = mol2[1].text
            
            # Number of interacting residues from each side
            self.n_int_res = {}
            self.n_int_res[self.chain1] = int(mol1[21].text)
            self.n_int_res[self.chain2] = int(mol2[21].text)
            
            # Number of interacting atoms from each side
            self.n_int_atoms = {}
            self.n_int_atoms[self.chain1] = int(mol1[20].text)
            self.n_int_atoms[self.chain2] = int(mol2[20].text)
                
            self.id = '{0}_{1}'.format(self.chain1,self.chain2)
            
            hbonds_set = Interaction_Set(interaction_type='h-bonds',\
                                         n_interactions=hbonds[0].text)
                
            hbonds_set.load_from_PISA_xml_bonds(hbonds[1:])
            hbonds_freq = hbonds_set.calc_frequency_by_residue()


            salts_set = Interaction_Set(interaction_type='salt-bridges',\
                                        n_interactions=salts[0].text)
                
            salts_set.load_from_PISA_xml_bonds(salts[1:])
            salts_freq = salts_set.calc_frequency_by_residue()

            self.interaction_sets = [hbonds_set, salts_set]

            return hbonds_freq, salts_freq

    def __repr__(self):
        return 'Interface {0} with interface area {1}'\
            .format(self.id,self.int_area)

try:
    import wget
except:
    try:
        import requests
    except:
        print('Neither wget nor requests library could be imported. Exiting.. Bye!')
        sys.exit(0)
    
try:
    import xml.etree.ElementTree as ET
except:
    print('xml library could not be imported. Exiting.. Bye!')
    sys.exit(0)

class PISA_Data(object):
    def __init__(self,pdb_id):
        LINK = 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?{0}'
        if wget:
            wget.download(LINK.format(pdb_id),out='{0}.xml'.format(pdb_id))
        elif requests:
            url = LINK.format(pdb_id)
            outf = '{0}.xml'.format(pdb_id))
            r = requests.get(url)
            open(outf,'wb').write(r.content)
            
        self.pdb_id = pdb_id
        tree = ET.parse('{0}.xml'.format(pdb_id))

        root = tree.getroot()

        # root.tag ==> pisa_interfaces (i.e. the start of the XML tree should be <pisa_interfaces>
        # there must only be one such tag in the xmlcache file (a tree can only have one root, hence
        # the need for 50 caches per interfaces$i-$j.xml file)
        interfacestatus = root[0].text

        self.protein_interfaces = {}
        self.ligand_interfaces = {}


        if interfacestatus == 'Ok':
            
            # define dictionary here just to avoid doing it more than once
            # trans = {'CYS':'C', 'ASP':'D', 'SER':'S', 'GLN':'Q', 'LYS':'K', 'ILE':'I', 'PRO':'P', 'THR':'T', 'PHE':'F', 'ASN':'N', 'GLY':'G', 'HIS':'H', 'LEU':'L', 'ARG':'R', 'TRP':'W', 'ALA':'A', 'VAL':'V', 'GLU':'E', 'TYR':'Y', 'MET':'M'}
            
            for pdbentry in root:
                # the first tag in root is actually <status> but having dealt with it in the if statement, ignore it
                if pdbentry.tag == 'status':
                    # the first child of root isn't actually a pdbentry, it's the status - treat this separately
                    if pdbentry.text != 'Ok':
                        # handle N/A value output
                        # no PDB code to report
                        errormsg = "Something wrong | Exiting.. "
                        print(errormsg)
                        break
                else:
                    # parse the <pdb_entry>
                    pdb_id = pdbentry[0].text.upper()
                    entrystatus = pdbentry[1].text
                    if entrystatus == 'Ok':
                        n_interfaces = pdbentry[2].text
                        print(n_interfaces, 'interfaces found!')
                        for i in range(0,int(n_interfaces)):
                            interface = pdbentry[3+i]
        
                            interface_obj = Interface(pdb_id)
                            interface_obj.load_from_PISA_xml_interface(interface)
                            
                            if interface_obj.type=='Ligand':
                                self.ligand_interfaces[interface_obj.id] =  interface_obj
                            
                            elif interface_obj.type=='Protein':
                                self.protein_interfaces[interface_obj.id] =  interface_obj
        
                    else:
                        # no interface data on this PDB structure :-(
                        # log the pdb_id and entrystatus
                        pdbcount = 0
                        # count the number of PDB entries for the log message
                        for tags in root:
                            if tags.tag == "pdb_entry":
                                pdbcount+=1
                        errormsg = "Something wrong | Exiting.. "
                        print(errormsg)
                        # e.g. Interface 1 in interfaces1-3.xml, entry 12 (3KCH) : Entry not found
        else:
            # log the interfacestatus
            errormsg = "Bad interface status | Exiting.. "
            print(errormsg)
    
    def __repr__(self):
        return '{0} protein interfaces\n{1} protein interfaces'\
            .format(len(self.protein_interfaces),len(self.ligand_interfaces))
    
    
