import configparser

### config file looks like that:
# [SECTION_NAME]
# key1 = value1
# key2 = value2

config = configparser.RawConfigParser()
config.optionxform = str ### preserve case sensitivity

def read(process):
   fname = '../setup/config.'+process+'.txt'
   print("Reading configuration from: ",fname)
   config.read(fname)

def get(section,var):
   return float(dict(config.items(section))[var])

def set(process,sides,doprint=False):
   map = {}
   ### read
   read(process)
   ### set
   map.update( {'me':get('Electron mass GeV','me')} )
   map.update( {'me2':map['me']*map['me']} )
   
   map.update( {'cm2m':get('Units','cm2m')} )
   map.update( {'cm2um':get('Units','cm2um')} )
   map.update( {'um2cm':get('Units','um2cm')} )
   
   map.update( {'B':get('Magnetic field Tesla, Meter','B')} )
   map.update( {'LB':get('Magnetic field Tesla, Meter','LB')} )
   
   map.update( {'Emax':get('Possible energies GeV','Emax')} )
   map.update( {'Emin':get('Possible energies GeV','Emin')} )

   map.update( {'Rbeampipe':get('General geometry cm','Rbeampipe')} )
   map.update( {'zDipoleExit':get('General geometry cm','zDipoleExit')} )
   map.update( {'xDipoleExitMinAbs':get('General geometry cm','xDipoleExitMinAbs')} )
   map.update( {'xDipoleExitMaxAbs':get('General geometry cm','xDipoleExitMaxAbs')} )
   map.update( {'yDipoleExitMin':get('General geometry cm','yDipoleExitMin')} )
   map.update( {'yDipoleExitMax':get('General geometry cm','yDipoleExitMax')} )
   map.update( {'xAbsMargins':get('General geometry cm','xAbsMargins')} )
   map.update( {'yAbsMargins':get('General geometry cm','yAbsMargins')} )
   
   map.update( {'Hstave':get('Stave geometry cm','Hstave')} )
   map.update( {'Lstave':get('Stave geometry cm','Lstave')} )
   map.update( {'RoffsetBfield':get('Stave geometry cm','RoffsetBfield')} )
   
   map.update( {'xPsideL':-map['RoffsetBfield']-map['Lstave']} )
   map.update( {'xPsideR':-map['RoffsetBfield']} )
   map.update( {'xEsideL':+map['RoffsetBfield']} )
   map.update( {'xEsideR':+map['RoffsetBfield']+map['Lstave']} )
   map.update( {'yUp':+map['Hstave']/2.} )
   map.update( {'yDn':-map['Hstave']/2.} )
   
   map.update( {'detXmin':map['xPsideL']} )
   map.update( {'detXmax':map['xEsideR']} )
   if(process=="trident"): map['detXmax'] = map['xPsideR']
   if(process=="bppp" and sides=="e+"): map['detXmax'] = map['xPsideR']
   if(process=="bppp" and sides=="e-"): map['detXmin'] = map['xEsideL']
   
   if(doprint):
      print("Configuration map:",map)
      print("")
   return map

# set("bppp","e+e-")