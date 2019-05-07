#!/usr/bin/env python3

#Example usage from GK_HEAD_DIR
#> scripts/config_generation.py $(find ./src ./externals -name '*.[fF]??')
#Will produce one file for each detected config type in the source
#with the boiler plate routines and canonical form of the type.

typeIndent="   "
indent="  "
declString = "{typ} :: {name}"
assignmentString = " = {value}"

readParamsSubString = """
type({module}_config_type), intent(in), optional :: config_in

if (present(config_in)) config_internal = config_in

call config_internal%init(name = '{namelist}', requires_index = {req})

! Copy out internal values into module level parameters
{COPYOUT_LIST}

exist = config_internal%exist"""

readSubString = """        
!> Reads in the {namelist} namelist and populates the member variables
subroutine read_{module}_config(self)
  use file_utils, only: input_unit_exist, get_indexed_namelist_unit
  use mp, only: proc0
  implicit none
  class({module}_config_type), intent(in out) :: self
  logical :: exist
  integer :: in_file, iostat

  ! Note: When this routine is in the module where these variables live
  ! we shadow the module level variables here. This is intentional to provide
  ! isolation and ensure we can move this routine to another module easily.    
{DECL_LIST}

  namelist /{namelist}/ {NAMELIST_LIST}

  ! Only proc0 reads from file
  if (.not. proc0) return

  ! First set local variables from current values
{COPYOUT_LIST}

  ! Now read in the main namelist
  if (self%requiresIndex) then
    call get_indexed_namelist_unit(in_file, trim(self%name), self%index)
    read(in_file, nml = {namelist}, iostat = iostat)
    exist = .true.
    if (iostat .ne.0) exist = .false.
    close(unit = in_file, iostat = iostat)
  else
    in_file = input_unit_exist(self%name, exist)
    if (exist) read(in_file, nml = {namelist})
  endif

  ! Now copy from local variables into type members
{COPYIN_LIST}

  self%exist = exist
end subroutine read_{module}_config"""

writeSubString = """
!> Writes out a namelist representing the current state of the config object
subroutine write_{module}_config(self, unit)
  implicit none
  class({module}_config_type), intent(in) :: self
  integer, intent(in) , optional:: unit
  integer :: unit_internal

  unit_internal = 6 ! @todo -- get stdout from file_utils
  if (present(unit)) then
     unit_internal = unit
  endif

  call self%write_namelist_header(unit_internal)
{KEYVAL_LIST}
  call self%write_namelist_footer(unit_internal)
end subroutine write_{module}_config"""

resetSubString = """
!> Resets the config object to the initial empty state
subroutine reset_{module}_config(self)
  class({module}_config_type), intent(in out) :: self
  type({module}_config_type) :: empty
  select type (self)
  type is ({module}_config_type)
     self = empty
  end select
end subroutine reset_{module}_config"""

broadcastSubString = """
!> Broadcasts all config parameters so object is populated identically on
!! all processors
subroutine broadcast_{module}_config(self)
  use mp, only: broadcast
  implicit none
  class({module}_config_type), intent(in out) :: self
{BROADCAST_LIST}

  call broadcast(self%exist)
end subroutine broadcast_{module}_config"""

class var():
    assignInString = "{instanceName}%{name} = {name}"
    assignOutString = "{name} = {instanceName}%{name}"
    broadcastString = "call broadcast({instanceName}%{name})"
    writeKeyValString = "call {instanceName}%write_key_val("+'"{name}"'+", {instanceName}%{name}, unit_internal)"    
    
    def __init__(self, name = None, typ = None, default = None, comment = "FIXME: Add documentation"):
        self.name = name.lower()
        self.typ  = typ.lower()
        if isinstance(comment, list):
            comment = "\n".join(comment)
        
        self.comment = comment.split("\n")
        self.default = default

        #Bit of massaging for string inputs
#        if isinstance(self.default, str):
        if self.typ.startswith("character"):
            #Already quoted so do nothing
            if self.default[0] == '"' or self.default[0] == "'":
                pass
            else:
                #Not quoted so add quotes
                self.default = '"' + self.default + '"'
                
    def __str__(self):
        st = typeIndent+"!> "+("\n"+typeIndent+"!> ").join(self.comment)+"\n"
        st += typeIndent+declString.format(typ=self.typ, name=self.name)
        if self.default is not None:
            st+= assignmentString.format(value = self.default)
        return st

    def assignIn(self, instanceName = 'self'):
        return self.assignInString.format(instanceName = instanceName, name = self.name)

    def assignOut(self, instanceName = 'self'):
        return self.assignOutString.format(instanceName = instanceName, name = self.name)

    def broadcast(self, instanceName = 'self'):
        return self.broadcastString.format(instanceName = instanceName, name = self.name)

    def keyval(self, instanceName = 'self'):
        return self.writeKeyValString.format(instanceName = instanceName, name = self.name)
    
class config():
    startTypeString = "type, extends(abstract_config_type) :: {module}_config_type"
    endTypeString   = "end type {module}_config_type"
    boundString     = """
 contains
{indent}procedure, public :: read => read_{module}_config
{indent}procedure, public :: write => write_{module}_config
{indent}procedure, public :: reset => reset_{module}_config
{indent}procedure, public :: broadcast => broadcast_{module}_config"""

    instanceString  = "type({module}_config_type) :: {instance_name}" 

    subComment = """
!---------------------------------------
! Following is for the config_type
!---------------------------------------
"""
    
    def __init__(self, varList = [], moduleName = None, namelistName = None, requiresIndex = False, comment = "FIXME: Add documentation"):
        self.varList = set(varList)
        self.moduleName = moduleName.lower()
        self.namelistName = namelistName.lower()
        self.requiresIndex = requiresIndex
        if isinstance(comment, list):
            comment = "\n".join(comment)
        
        self.comment = comment.split("\n")

        #Now sort vars
        self.sortVars()
        
    def __str__(self):
        st = "!> "+"\n!> ".join(self.comment)+"\n"
        st += self.startTypeString.format(module = self.moduleName) + "\n"
        st += typeIndent+"! namelist : "+self.namelistName + "\n"
        st += "\n".join([str(x) for x in self.varList])
        st += self.boundString.format(module = self.moduleName, indent = typeIndent)
        st += "\n"+self.endTypeString.format(module = self.moduleName)
        return st

    def addVar(self, varIn):
        self.varList.add(varIn)
        self.sortVars()

    def sortVars(self):
        self.varList = sorted(self.varList, key = lambda x: x.name)
        
    def namelistList(self):
        return ", ".join([x.name for x in self.varList])

    def copyInList(self, **kwargs):
        st = [v.assignIn(**kwargs) for v in self.varList]
        return indent+("\n"+indent).join(st)

    def copyOutList(self, **kwargs):
        st = [v.assignOut(**kwargs) for v in self.varList]
        return indent+("\n"+indent).join(st)

    def keyValList(self, **kwargs):
        st = [v.keyval(**kwargs) for v in self.varList]
        return indent+("\n"+indent).join(st)

    def broadcastList(self, **kwargs):
        st = [v.broadcast(**kwargs) for v in self.varList]
        return indent+("\n"+indent).join(st)
    
    def declarationList(self):
        perTypList = {}
        #First get an empty list per type
        for var in self.varList:
            perTypList[var.typ] = []
            
        #Now add each var to the appropriate list
        for var in self.varList:
            perTypList[var.typ].append(var.name)

        #Now convert to comma separated list
        st = []
        for k, v in perTypList.items():
            st.append(declString.format(typ = k, name = ", ".join(v)))

        return indent+("\n"+indent).join(st)
        
        
    def decl(self):
        return (str(self))

    def instance(self, instance_name = "config_internal"):
        return self.instanceString.format(module = self.moduleName, instance_name = instance_name)

    def readParams(self, **kwargs):
        req = ".false."
        if self.requiresIndex:
            req = ".true."
            
        return readParamsSubString.format(module = self.moduleName, namelist = self.namelistName, COPYOUT_LIST = self.copyOutList(instanceName='config_internal', **kwargs), req=req)
    
    def readConfig(self, **kwargs):
        return readSubString.format(namelist = self.namelistName, module = self.moduleName,
                                    COPYIN_LIST = self.copyInList(**kwargs), COPYOUT_LIST = self.copyOutList(**kwargs),
                                    NAMELIST_LIST = self.namelistList(**kwargs), DECL_LIST=self.declarationList(**kwargs))

    def writeConfig(self, **kwargs):
        return writeSubString.format(namelist = self.namelistName, module = self.moduleName,KEYVAL_LIST = self.keyValList(**kwargs))

    def resetConfig(self, **kwargs):
        return resetSubString.format(module = self.moduleName)

    def broadcastConfig(self, **kwargs):

        return broadcastSubString.format(module = self.moduleName, BROADCAST_LIST = self.broadcastList(**kwargs))
    

    
    def printAll(self,**kwargs):
        st = []
        st.append(self.decl(**kwargs))
        st.append("\n"+self.instance(**kwargs))
        st.append(self.readParams())        
        st.append(self.subComment)
        st.append(self.readConfig())
        st.append(self.writeConfig())
        st.append(self.resetConfig())
        st.append(self.broadcastConfig())
        return "\n".join(st)

    def dumpToFile(self,fileOut=None):
        if fileOut is None:
            fileOut = self.moduleName+"_placeholders.f90"

        with open(fileOut, 'w') as f:
            f.writelines(self.printAll())
        
def constructFromFile(fileIn=None, verbose=False):
    # Currently not able to get value of requires_index automatically.

    if fileIn is None:
        raise ValueError("Need to specify which file to read")
    
    import re

    commentPattern = re.compile(r"""\s*![!>]""", re.VERBOSE)            
    def getComments(theList, commentPattern = commentPattern , stripCom = False):

        #Search backwards for all connected comments
        comments = []
        commentMatchFound = True
        i = 1
        while commentMatchFound:
            commentMatch = commentPattern.match(theList[-i])
            if commentMatch is not None:
                txt = theList[-i]
                if stripCom:
                    txt = txt.strip()
                    txt = txt.split("!>")
                    if len(txt) > 1:
                        txt = ''.join(txt[1:])
                    else:
                        txt = txt[0]
                    txt = txt.split("!!")
                    if len(txt) > 1:
                        txt = ''.join(txt[1:])
                    else:
                        txt = txt[0]

                txt = txt.lstrip()
                comments.append(txt)
                i += 1
            else:
                commentMatchFound = False

        if len(comments) == 0:
            return "FIXME: Add documentation"
        return comments[::-1]
    
    startPattern = re.compile(r"""\s*type,\s*extends\s*\(abstract_config_type\)\s*::\s*(?P<moduleName>.*?)_config_type""", re.VERBOSE)
    endPattern = re.compile(r"""\s*end\s*type\s*(?P<moduleName>.*?)_config_type""", re.VERBOSE)

    declPattern = re.compile(r"""\s*(?P<dtype>(double precision|real|logical|character|integer|type)+?)(?P<dtypeExtra>.*?)::\s*(?P<vname>[^\s=]*)\s*(=\s*(?P<default>.*))?""",re.VERBOSE)
    nmlPattern = re.compile(r"""\s*![^!>]\s*namelist\s*:\s*(?P<namelist>[^>\s]+)""", re.VERBOSE)
    
    result = None
    theLines = {}
    theVars = {}
    theConfigs = {}
    
    with open(fileIn, 'r') as f:
        found = False
        currName = 'ERROR'
        otherLines = []
        namelistName = 'ERROR'
        theComments = None
        
        for line in f.readlines():
            if found:
                theLines[currName].append(line)
            else:
                otherLines.append(line)
                
            startMatch = startPattern.match(line)
            endMatch   = endPattern.match(line)
            
            if startMatch:
                name = startMatch.group("moduleName")                
                currName = name

                theComments = getComments(otherLines[:-1], stripCom = True)
                ncom = len(theComments)
                theVars[currName]  = []
                theLines[currName] = otherLines[-1-ncom:-1]
                theLines[currName].append(line)
                
                found = True
                
            elif endMatch:
                if not found: continue
                #Note assumption one line can't match start and end
                name = endMatch.group("moduleName")
                if name != currName:
                    raise ValueError("Mismatched type names {a} and {b}".format(a=name, b=currName))

                theConfigs[currName] = config(varList = theVars[currName], moduleName = currName, namelistName = namelistName, comment = theComments)
                
                found = False
                currName = 'ERROR'
                namelistName = 'ERROR'
                theComments = None
                
            elif found:
                # Currently inside type so work out what sort of line we have
                declMatch = declPattern.match(line)
                nmlMatch = nmlPattern.match(line)
                
                if declMatch is not None:
                    varComments = getComments(theLines[currName][:-1], stripCom = True)
                    vv = var(name = declMatch.group('vname'),
                             typ = declMatch.group('dtype')+declMatch.group('dtypeExtra').strip(),
                             default = declMatch.group('default'),
                             comment = varComments)

                    theVars[currName].append(vv)
                elif nmlMatch is not None:
                    namelistName = nmlMatch.group('namelist')
                    
    if(len(theLines.keys())>0):
        result = theConfigs
        if verbose:        
            print("In file {f} found the following :".format(f=fileIn))
            
            for k,v in theLines.items():
                print(k+" : ")
                print(''.join(v))
                print('')
                print('\n'.join(str(x) for x in theVars[k]))
    else:
        if verbose: print("No config types found in {f}".format(f=fileIn))


    return result

if __name__ == "__main__":
    import sys

    probablyFiles = sys.argv[1:]

    for fil in probablyFiles:
        tmp = constructFromFile(fileIn=fil)
        if tmp is None: continue
        for k, v in tmp.items():
            print("Config '{k}' in file '{f}'".format(k=k,f=fil))
            v.dumpToFile()
            
