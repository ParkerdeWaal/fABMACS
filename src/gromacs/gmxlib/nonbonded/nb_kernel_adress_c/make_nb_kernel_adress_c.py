#!/usr/bin/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

import sys
import os
sys.path.append ( "../preprocessor" )
from gmxpreprocess import gmxpreprocess

# "The happiest programs are programs that write other programs."
#
#
# This script controls the generation of Gromacs nonbonded kernels.
#
# We no longer generate kernels on-the-fly, so this file is not run
# during a Gromacs compile - only when we need to update the kernels (=rarely).
#
# To maximize performance, each combination of interactions in Gromacs
# has a separate nonbonded kernel without conditionals in the code.
# To avoid writing hundreds of different routines for each architecture,
# we instead use a custom preprocessor so we can encode the conditionals
# and expand for-loops (e.g, for water-water interactions)
# from a general kernel template. While that file will contain quite a
# few preprocessor directives, it is still an order of magnitude easier
# to maintain than ~200 different kernels (not to mention it avoids bugs).
#
# To actually generate the kernels, this program iteratively calls the
# preprocessor with different define settings corresponding to all
# combinations of coulomb/van-der-Waals/geometry options.
#
# A main goal in the design was to make this new generator _general_. For
# this reason we have used a lot of different fields to identify a particular
# kernel and interaction. Basically, each kernel will have a name like
#
# nbkernel_ElecXX_VdwYY_GeomZZ_VF_QQ()
#
# Where XX/YY/ZZ/VF are strings to identify what the kernel computes.
#
# Elec/Vdw describe the type of interaction for electrostatics and van der Waals.
# The geometry settings correspond e.g. to water-water or water-particle kernels,
# and finally the VF setting is V,F,or VF depending on whether we calculate
# only the potential, only the force, or both of them. The final string (QQ)
# is the architecture/language/optimization of the kernel.
#
Arch       = 'c'

# Explanation of the 'properties':
#
# It is cheap to compute r^2, and the kernels require various other functions of r for
# different kinds of interaction. Depending on the needs of the kernel and the available
# processor instructions, this will be done in different ways.
#
# 'rinv' means we need 1/r, which is calculated as 1/sqrt(r^2).
# 'rinvsq' means we need 1/(r*r). This is calculated as rinv*rinv if we already did rinv, otherwise 1/r^2.
# 'r' is similarly calculated as r^2*rinv when needed
# 'table' means the interaction is tabulated, in which case we will calculate a table index before the interaction
# 'shift' means the interaction will be modified by a constant to make it zero at the cutoff.
# 'cutoff' means the interaction is set to 0.0 outside the cutoff
#

FileHeader = \
'/*\n' \
' * This file is part of the GROMACS molecular simulation package.\n' \
' *\n' \
' * Copyright (c) 2012, by the GROMACS development team, led by\n' \
' * David van der Spoel, Berk Hess, Erik Lindahl, and including many\n' \
' * others, as listed in the AUTHORS file in the top-level source\n' \
' * directory and at http://www.gromacs.org.\n' \
' *\n' \
' * GROMACS is free software; you can redistribute it and/or\n' \
' * modify it under the terms of the GNU Lesser General Public License\n' \
' * as published by the Free Software Foundation; either version 2.1\n' \
' * of the License, or (at your option) any later version.\n' \
' *\n' \
' * GROMACS is distributed in the hope that it will be useful,\n' \
' * but WITHOUT ANY WARRANTY; without even the implied warranty of\n' \
' * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n' \
' * Lesser General Public License for more details.\n' \
' *\n' \
' * You should have received a copy of the GNU Lesser General Public\n' \
' * License along with GROMACS; if not, see\n' \
' * http://www.gnu.org/licenses, or write to the Free Software Foundation,\n' \
' * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.\n' \
' *\n' \
' * If you want to redistribute modifications to GROMACS, please\n' \
' * consider that scientific software is very special. Version\n' \
' * control is crucial - bugs must be traceable. We will be happy to\n' \
' * consider code for inclusion in the official distribution, but\n' \
' * derived work must not be called official GROMACS. Details are found\n' \
' * in the README & COPYING files - if they are missing, get the\n' \
' * official version at http://www.gromacs.org.\n' \
' *\n' \
' * To help us fund GROMACS development, we humbly ask that you cite\n' \
' * the research papers on the package. Check out http://www.gromacs.org.\n' \
' */\n' \
'/*\n' \
' * Note: this file was generated by the GROMACS '+Arch+' kernel generator.\n' \
' */\n'

###############################################
# ELECTROSTATICS
# Interactions and flags for them
###############################################
ElectrostaticsList = {
    'None'                    : [],
    'Coulomb'                 : ['rinv','rinvsq'],
    'ReactionField'           : ['rinv','rinvsq'],
    'GeneralizedBorn'         : ['rinv','r'],
    'CubicSplineTable'        : ['rinv','r','table'],
    'Ewald'                   : ['rinv','rinvsq','r'],
}


###############################################
# VAN DER WAALS
# Interactions and flags for them
###############################################
VdwList = {
    'None'                    : [],
    'LennardJones'            : ['rinvsq'],
    'Buckingham'              : ['rinv','rinvsq','r'],
    'CubicSplineTable'        : ['rinv','r','table'],
}


###############################################
# MODIFIERS
# Different ways to adjust/modify interactions to conserve energy
###############################################
ModifierList = {
    'None'                    : [],
    'ExactCutoff'             : ['exactcutoff'],        # Zero the interaction outside the cutoff, used for reaction-field-zero
    'PotentialShift'          : ['shift','exactcutoff'],
    'PotentialSwitch'         : ['rinv','r','switch','exactcutoff']
}


###############################################
# GEOMETRY COMBINATIONS
###############################################
GeometryNameList = [
    [ 'Particle' , 'Particle' ],
    [ 'Water3'   , 'Particle' ],
    [ 'Water3'   , 'Water3'   ],     
    [ 'Water4'   , 'Particle' ],
    [ 'Water4'   , 'Water4'   ]
]


###############################################
# POTENTIAL / FORCE
###############################################
VFList = [
    'PotentialAndForce',
# 'Potential',   # Not used yet
    'Force'
]


###############################################
# GEOMETRY PROPERTIES
###############################################
# Dictionaries with lists telling which interactions are present
# 1,2,3 means particles 1,2,3 (but not 0) have electrostatics!
GeometryElectrostatics = {
    'Particle'  : [ 0 ],
    'Particle2' : [ 0 , 1 ],
    'Particle3' : [ 0 , 1 , 2 ],
    'Particle4' : [ 0 , 1 , 2 , 3 ],
    'Water3'    : [ 0 , 1 , 2 ],
    'Water4'    : [ 1 , 2 , 3 ]
}

GeometryVdw = {
    'Particle'  : [ 0 ],
    'Particle2' : [ 0 , 1 ],
    'Particle3' : [ 0 , 1 , 2 ],
    'Particle4' : [ 0 , 1 , 2 , 3 ],
    'Water3'    : [ 0 ],
    'Water4'    : [ 0 ]
}


###############################################
# ADRESS PROPERTIES
###############################################
AdressList = [
    'CG', 'EX'
]

# Dictionary to abbreviate all strings (mixed from all the lists)
Abbreviation = {
    'None'                    : 'None',
    'Coulomb'                 : 'Coul',
    'Ewald'                   : 'Ew',
    'ReactionField'           : 'RF',
    'GeneralizedBorn'         : 'GB',
    'CubicSplineTable'        : 'CSTab',
    'LennardJones'            : 'LJ',
    'Buckingham'              : 'Bham',
    'PotentialShift'          : 'Sh',
    'PotentialSwitch'         : 'Sw',
    'ExactCutoff'             : 'Cut',
    'PotentialAndForce'       : 'VF',
    'Potential'               : 'V',
    'Force'                   : 'F',
    'Water3'                  : 'W3',
    'Water4'                  : 'W4',
    'Particle'                : 'P1',
    'Particle2'               : 'P2',
    'Particle3'               : 'P3',
    'Particle4'               : 'P4'
}


###############################################
# Functions
###############################################

# Return a string with the kernel name from current settings
def MakeKernelFileName(KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelRes):
    ElecStr = 'Elec' + Abbreviation[KernelElec]
    if(KernelElecMod!='None'):
        ElecStr = ElecStr + Abbreviation[KernelElecMod]
    VdwStr  = 'Vdw'  + Abbreviation[KernelVdw]
    if(KernelVdwMod!='None'):
        VdwStr = VdwStr + Abbreviation[KernelVdwMod]
    GeomStr = 'Geom' + Abbreviation[KernelGeom[0]] + Abbreviation[KernelGeom[1]]
    return 'nb_kernel_' + ElecStr + '_' + VdwStr + '_' + GeomStr + '_' + KernelRes + '_' + Arch

def MakeKernelName(KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelVF,KernelRes):
    ElecStr = 'Elec' + Abbreviation[KernelElec]
    if(KernelElecMod!='None'):
        ElecStr = ElecStr + Abbreviation[KernelElecMod]
    VdwStr  = 'Vdw'  + Abbreviation[KernelVdw]
    if(KernelVdwMod!='None'):
        VdwStr = VdwStr + Abbreviation[KernelVdwMod]
    GeomStr = 'Geom' + Abbreviation[KernelGeom[0]] + Abbreviation[KernelGeom[1]]
    VFStr   = Abbreviation[KernelVF]
    return 'nb_kernel_' + ElecStr + '_' + VdwStr + '_' + GeomStr + '_' + VFStr + '_' + KernelRes + '_' + Arch

# Return a string with a declaration to use for the kernel;
# this will be a sequence of string combinations as well as the actual function name
# Dont worry about field widths - that is just pretty-printing for the header!
def MakeKernelDecl(KernelName,KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelOther,KernelVF):
    KernelStr   = '\"'+KernelName+'\"'
    ArchStr     = '\"'+Arch+'\"'
    ElecStr     = '\"'+KernelElec+'\"'
    ElecModStr  = '\"'+KernelElecMod+'\"'
    VdwStr      = '\"'+KernelVdw+'\"'
    VdwModStr   = '\"'+KernelVdwMod+'\"'
    GeomStr     = '\"'+KernelGeom[0]+KernelGeom[1]+'\"'
    OtherStr    = '\"'+KernelOther+'\"'
    VFStr       = '\"'+KernelVF+'\"'

    ThisSpec = ArchStr+', '+ElecStr+', '+ElecModStr+', '+VdwStr+', '+VdwModStr+', '+GeomStr+', '+OtherStr+', '+VFStr
    ThisDecl = '    { '+KernelName+', '+KernelStr+', '+ThisSpec+' }'
    return ThisDecl


# Returns 1 if this kernel should be created, 0 if we should skip it
# This routine is not critical - it is not the end of the world if we create more kernels,
# but since the number is pretty large we save both space and compile-time by reducing it a bit.
def KeepKernel(KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelVF):

    # No need for kernels without interactions
    if(KernelElec=='None' and KernelVdw=='None'):
        return 0

    # No need for modifiers without interactions
    if((KernelElec=='None' and KernelElecMod!='None') or (KernelVdw=='None' and KernelVdwMod!='None')):
        return 0

    # No need for LJ-only water optimization, or water optimization with implicit solvent.
    if('Water' in KernelGeom[0] and (KernelElec=='None' or 'GeneralizedBorn' in KernelElec)):
        return 0

    # Non-matching table settings are pointless
    if( ('Table' in KernelElec) and ('Table' in KernelVdw) and KernelElec!=KernelVdw ):
        return 0

    # Try to reduce the number of different switch/shift options to get a reasonable number of kernels
    # For electrostatics, reaction-field can use 'exactcutoff', and ewald can use switch or shift.
    if(KernelElecMod=='ExactCutoff' and KernelElec!='ReactionField'):
        return 0
    if(KernelElecMod in ['PotentialShift','PotentialSwitch'] and KernelElec!='Ewald'):
        return 0
    # For Vdw, we support switch and shift for Lennard-Jones/Buckingham
    if((KernelVdwMod=='ExactCutoff') or
       (KernelVdwMod in ['PotentialShift','PotentialSwitch'] and KernelVdw not in ['LennardJones','Buckingham'])):
        return 0

    # Choose either switch or shift and don't mix them...
    if((KernelElecMod=='PotentialShift' and KernelVdwMod=='PotentialSwitch') or
       (KernelElecMod=='PotentialSwitch' and KernelVdwMod=='PotentialShift')):
        return 0

    # Don't use a Vdw kernel with a modifier if the electrostatics one does not have one
    if(KernelElec!='None' and KernelElecMod=='None' and KernelVdwMod!='None'):
        return 0

    # Don't use an electrostatics kernel with a modifier if the vdw one does not have one,
    # unless the electrostatics one is reaction-field with exact cutoff.
    if(KernelVdw!='None' and KernelVdwMod=='None' and KernelElecMod!='None'):
        if(KernelElec=='ReactionField' and KernelVdw!='CubicSplineTable'):
            return 0
        elif(KernelElec!='ReactionField'):
            return 0

    return 1



#
# The preprocessor will automatically expand the interactions for water and other
# geometries inside the kernel, but to get this right we need to setup a couple
# of defines - we do them in a separate routine to keep the main loop clean.
#
# While this routine might look a bit complex it is actually quite straightforward,
# and the best news is that you wont have to modify _anything_ for a new geometry
# as long as you correctly define its Electrostatics/Vdw geometry in the lists above!
#
def SetDefines(KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelVF,defines):
    # What is the _name_ for the i/j group geometry?
    igeometry            = KernelGeom[0]
    jgeometry            = KernelGeom[1]
    # define so we can access it in the source when the preprocessor runs
    defines['GEOMETRY_I'] = igeometry
    defines['GEOMETRY_J'] = jgeometry

    # For the i/j groups, extract a python list of which sites have electrostatics
    # For SPC/TIP3p this will be [1,1,1], while TIP4p (no elec on first site) will be [0,1,1,1]
    ielec                = GeometryElectrostatics[igeometry]
    jelec                = GeometryElectrostatics[jgeometry]
    # Zero out the corresponding lists in case we dont do Elec
    if(KernelElec=='None'):
        ielec = []
        jelec = []

    # Extract similar interaction lists for Vdw interactions (example for SPC: [1,0,0])
    iVdw                 = GeometryVdw[igeometry]
    jVdw                 = GeometryVdw[jgeometry]

    # Zero out the corresponding lists in case we dont do Vdw
    if(KernelVdw=='None'):
        iVdw = []
        jVdw = []

    # iany[] and jany[] contains lists of the particles actually used (for interactions) in this kernel
    iany = list(set(ielec+iVdw))  # convert to+from set to make elements unique
    jany = list(set(jelec+jVdw))

    defines['PARTICLES_ELEC_I'] = ielec
    defines['PARTICLES_ELEC_J'] = jelec
    defines['PARTICLES_VDW_I']  = iVdw
    defines['PARTICLES_VDW_J']  = jVdw
    defines['PARTICLES_I']      = iany
    defines['PARTICLES_J']      = jany

    # elecij,Vdwij are sets with pairs of particles for which the corresponding interaction is done
    # (and anyij again corresponds to either electrostatics or Vdw)
    elecij = []
    Vdwij  = []
    anyij  = []

    for i in ielec:
        for j in jelec:
            elecij.append([i,j])

    for i in iVdw:
        for j in jVdw:
            Vdwij.append([i,j])

    for i in iany:
        for j in jany:
            if [i,j] in elecij or [i,j] in Vdwij:
                anyij.append([i,j])

    defines['PAIRS_IJ']     = anyij

    # Make an 2d list-of-distance-properties-to-calculate for i,j
    ni = max(iany)+1
    nj = max(jany)+1
    # Each element properties[i][j] is an empty list
    properties = [ [ [] for j in range(0,nj) ] for i in range (0,ni) ]
    # Add properties to each set
    for i in range(0,ni):
        for j in range(0,nj):
            if [i,j] in elecij:
                properties[i][j] = properties[i][j] + ['electrostatics'] + ElectrostaticsList[KernelElec] + ModifierList[KernelElecMod]
            if [i,j] in Vdwij:
                properties[i][j] = properties[i][j] + ['vdw'] + VdwList[KernelVdw] + ModifierList[KernelVdwMod]
            # Add rinv if we need r
            if 'r' in properties[i][j]:
                properties[i][j] = properties[i][j] + ['rinv']
            # Add rsq if we need rinv or rinsq
            if 'rinv' in properties[i][j] or 'rinvsq' in properties[i][j]:
                properties[i][j] = properties[i][j] + ['rsq']

    defines['INTERACTION_FLAGS']    = properties



def PrintStatistics(ratio):
    ratio = 100.0*ratio
    print '\rGenerating %s nonbonded kernels... %5.1f%%' % (Arch,ratio),
    sys.stdout.flush()



defines = {}
kerneldecl = []

cnt     = 0.0
nelec   = len(ElectrostaticsList)
nVdw    = len(VdwList)
nmod    = len(ModifierList)
ngeom   = len(GeometryNameList)

ntot    = nelec*nmod*nVdw*nmod*ngeom

numKernels = 0

fpdecl = open('nb_kernel_' + Arch + '.c','w')
fpdecl.write( FileHeader )
fpdecl.write( '#ifndef nb_kernel_' + Arch + '_h\n' )
fpdecl.write( '#define nb_kernel_' + Arch + '_h\n\n' )
fpdecl.write( '#include "../nb_kernel.h"\n\n' )

for KernelElec in ElectrostaticsList:
    defines['KERNEL_ELEC'] = KernelElec

    for KernelElecMod in ModifierList:
        defines['KERNEL_MOD_ELEC'] = KernelElecMod

        for KernelVdw in VdwList:
            defines['KERNEL_VDW'] = KernelVdw

            for KernelVdwMod in ModifierList:
                defines['KERNEL_MOD_VDW'] = KernelVdwMod

                for KernelGeom in GeometryNameList:

                    cnt += 1
                    KernelFilename = MakeKernelFileName(KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom) + '.c'
                    fpkernel = open(KernelFilename,'w')
                    defines['INCLUDE_HEADER'] = 1  # Include header first time in new file
                    DoHeader = 1

                    for KernelVF in VFList:

                        KernelName = MakeKernelName(KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelVF)

                        defines['KERNEL_NAME'] = KernelName
                        defines['KERNEL_VF']   = KernelVF

                        # Check if this is a valid/sane/usable combination
                        if not KeepKernel(KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelVF):
                            continue;

                        # The overall kernel settings determine what the _kernel_ calculates, but for the water
                        # kernels this does not mean that every pairwise interaction has e.g. Vdw interactions.
                        # This routine sets defines of what to calculate for each pair of particles in those cases.
                        SetDefines(KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelVF,defines)

                        if(DoHeader==1):
                            fpkernel.write( FileHeader )

                        gmxpreprocess('nb_kernel_template_' + Arch + '.pre', KernelName+'.tmp' , defines, force=1,contentType='C')
                        numKernels = numKernels + 1

                        defines['INCLUDE_HEADER'] = 0   # Header has been included once now
                        DoHeader=0

                        # Append temp file contents to the common kernelfile
                        fptmp = open(KernelName+'.tmp','r')
                        fpkernel.writelines(fptmp.readlines())
                        fptmp.close()
                        os.remove(KernelName+'.tmp')

                        # Add a declaration for this kernel
                        fpdecl.write('nb_kernel_t ' + KernelName + ';\n');

                        # Add declaration to the buffer
                        KernelOther=''
                        kerneldecl.append(MakeKernelDecl(KernelName,KernelElec,KernelElecMod,KernelVdw,KernelVdwMod,KernelGeom,KernelOther,KernelVF))

                    filesize = fpkernel.tell()
                    fpkernel.close()
                    if(filesize==0):
                        os.remove(KernelFilename)

                    PrintStatistics(cnt/ntot)
                pass
            pass
        pass
    pass
pass

# Write out the list of settings and corresponding kernels to the declaration file
fpdecl.write( '\n\n' )
fpdecl.write( 'nb_kernel_info_t\n' )
fpdecl.write( '    kernellist_'+Arch+'[] =\n' )
fpdecl.write( '{\n' )
for decl in kerneldecl[0:-1]:
    fpdecl.write( decl + ',\n' )
fpdecl.write( kerneldecl[-1] + '\n' )
fpdecl.write( '};\n\n' )
fpdecl.write( 'int\n' )
fpdecl.write( '    kernellist_'+Arch+'_size = sizeof(kernellist_'+Arch+')/sizeof(kernellist_'+Arch+'[0]);\n\n')
fpdecl.write( '#endif\n')
fpdecl.close()
