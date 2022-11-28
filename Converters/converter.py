#!/usr/bin/env python
# -----------------------------------------------------------------------
# General converter script
# -----------------------------------------------------------------------

import os,re,sys,optparse
from math import *

# -----------------------------------------------------------------------
parser = optparse.OptionParser()
(options, args) = parser.parse_args()

def finder(place,begin,end):
    motif = re.compile(begin + '.*' + end)
    obj = re.search(motif, place)
    texte=''
    if obj is None:
        motif = re.compile(begin + '.*' + end,re.DOTALL)
        obj = re.search(motif, place)
        if obj is not None:
            texte = obj.group()
            texte = re.sub(begin,'',texte)
            texte = re.sub(end,'',texte)
            texte = re.sub(r'\\','',texte)
            texte = re.sub('\\n','',texte)
            texte = re.sub(' +',' ',texte)
            texte = texte.strip()
            return texte
        else:
            return texte
    else:
        texte = obj.group()
        texte = re.sub(begin,'',texte)
        texte = re.sub(end,'',texte)
        texte = re.sub(r'\\','',texte)
        texte = re.sub('\\n','',texte)
        texte = re.sub(' +',' ',texte)
        texte = texte.strip()
        return texte


def abort(message):
    print message
    print "Exiting Converter"
    sys.exit()
# -----------------------------------------------------------------------



# -----------------------------------------------------------------------
# DMol3 or VASP ?
# -----------------------------------------------------------------------

if len(args)==0:                 # need at least a flag
   abort("Argument missing!\nUsage: %prog # [name]\n#: 1 for DMol3 (then name of structure needed), 2 for VASP")
elif len(args)==1:              # if only one, then it has to be VASP
    if not args[0]==str(2):
       abort("Usage: %prog # [name]\n#: 1 for DMol3 (then name of structure needed), 2 for VASP")
    flag=int(args[0])
elif len(args)==2:              # if two, then FLAG and NAME
    flag=int(args[0])
    name=args[1]
else:
    abort("Too many argument")

#########################################################################################
#########################################################################################
#########################################################################################
    

if flag==1:
    
# -----------------------------------------------------------------------
# DMol3 routine
# -----------------------------------------------------------------------
         
    if not os.path.exists(name+'_potential.grd'):               #check if name_potential.grd exists
       abort("No "+name+"_potential.grd could be found.")
    if not os.path.exists(name+'.car'):                         #check for name.car
       abort('No '+name+'.car was found')

    filin1=open(name+'.car','r')             #open files
    filin2=open(name+'_potential.grd','r')
    filout=open(name+'.cube','w')
    fulpath_out = name+".cube"

    # -----------------------------------------------------------------------
    # READING CAR DATAS
    # -----------------------------------------------------------------------

    data=filin1.read().split('\n')
    coord={}
    kind={}
    n=1
    for line in data:
        if re.search('PBC',line):
            boxcar=finder(line,'PBC *','$').split(' ')[0:3]             #info about the box
            anglecar=finder(line,'PBC *','$').split(' ')[3:]
        elif re.search('XXXX',line):
            coord[n]=finder(line,'','XXXX').split(' ')                  # coord and nature of the atoms
            kind[n]=finder(line,'xx *',' ')
            del coord[n][0]
            n=n+1

    # -----------------------------------------------------------------------
    # READING GRID DATAS
    # -----------------------------------------------------------------------

    data2=filin2.read().split('\n')
    boxgrid=finder(data2[2],'^','$').split(' ')[0:3]            #info about the box (to compare later with .car datas)
    anglegrid=finder(data2[2],'^','$').split(' ')[3:]
    grid=finder(data2[3],'','$').split(' ')                     #number of grid spaces
    info=finder(data2[4],'','$').split(' ')                     #info : fastest changing index and origin

    v={}
    for x in range(1,int(grid[0])+2):
        v[x]={}
        for y in range(1,int(grid[1])+2):
            v[x][y]={}

    count=5
    if int(info[0])==1:                     #extract the ESP, if x is the fastest moving...
        for z in range(1,int(grid[2])+2):
            for y in range(1,int(grid[1])+2):
                for x in range(1,int(grid[0])+2):
                    v[x][y][z]=data2[count]
                    count=count+1
    elif int(info[0])==3:                   #...if z is...(y impossible)
        for x in range(1,int(grid[0])+2):
            for y in range(1,int(grid[1])+2):
                for z in range(1,int(grid[2])+2):
                    v[x][y][z]=data2[count]
                    count=count+1

    # -----------------------------------------------------------------------
    # POSSIBLE ERRORS
    # -----------------------------------------------------------------------

    if [float(boxgrid[0]),float(boxgrid[1]),float(boxgrid[2])]<>[float(boxcar[0]),float(boxcar[1]),float(boxcar[2])]:
        abort('WARNING : grd and car don''t have the same box')
    if [float(anglegrid[0]),float(anglegrid[1]),float(anglegrid[2])]<>[float(anglecar[0]),float(anglecar[1]),float(anglecar[2])]:
        abort('WARNING : grd and car don''t have the same box')

    if len(coord)<>len(kind):
        abort('WARNING : problem when reading the .car file')

    if [float(anglegrid[0]),float(anglegrid[1]),float(anglegrid[2])]<>[90,90,90]:
        abort('WARNING : not orthogonal box: script not valid')

    if count-5<>len(data2)-6:
        abort('WARNING : problem with the number of grid points')

    # -----------------------------------------------------------------------
    # NEEDED CUBE DATAS
    # -----------------------------------------------------------------------

    nbr=len(coord)
    origin=[0,0,0]
    norm=[0,0,0]
    for i in range(0,3):                #step-size (conversion to Bohr)
        norm[i]=float(boxgrid[i])*1.889725985/(float(grid[i])+1)
        origin[i]=norm[i]*float(info[2*i+1])*(-1)

    
    atomic={'H':'1','He':'2','Li':'3','Be':'4','B':'5','C':'6','N':'7','O':'8','F':'9','Ne':'10','Na':'11','Mg':'12','Al':'13','Si':'14','P':'15','S':'16','Cl':'17','Ar':'18','K':'19','Ca':'20','Sc':'21','Ti':'22','V':'23','Cr':'24','Mn':'25','Fe':'26','Co':'27','Ni':'28','Cu':'29','Zn':'30','Ga':'31','Ge':'32','As':'33','Se':'34','Br':'35','Kr':'36','Rb':'37','Sr':'38','Y':'39','Zr':'40','Nb':'41','Mo':'42','Tc':'43','Ru':'44','Rh':'45','Pd':'46','Ag':'47','Cd':'48','In':'49','Sn':'50','Sb':'51','Te':'52','I':'53','Xe':'54','Cs':'55','Ba':'56','La':'57','Ce':'58','Pr':'59','Nd':'60','Pm':'61','Sm':'62','Eu':'63','Gd':'64','Th':'65','Dy':'66','Ho':'67','Er':'68','Tm':'69','Yb':'70','Lu':'71','Hf':'72','Ta':'73','W':'74','Re':'75','Os':'76','Ir':'77','Pt':'78','Au':'79','Hg':'80','Tl':'81','Pb':'82','Bi':'83','Po':'84','At':'85','Rn':'86','Fr':'87','Ra':'88','Ac':'89','Th':'90','Pa':'91','U':'92','Np':'93','Pu':'94','Am':'95','Cm':'96','Bk':'97','Cf':'98','Es':'99','Fm':'100','Md':'101','No':'102','Lw':'103'}  #convert nature of the atom into atomic number

    # -----------------------------------------------------------------------
    # WRITE FIRST LINES + COORD (conversion to Bohr)
    # -----------------------------------------------------------------------

    line0='--------------------REPEAT charges--------------------'
    line1='------cube file created from DMol3 (grd and car)------'
    line2='% 5s' %nbr+'% 12.6f' %origin[0]+'% 12.6f' %origin[1]+'% 12.6f' %origin[2]
    vect1='% 5s' %(int(grid[0])+1)+'% 12.6f' %float(norm[0])+'% 12.6f' %0+'% 12.6f' %0
    vect2='% 5s' %(int(grid[1])+1)+'% 12.6f' %0+'% 12.6f' %float(norm[1])+'% 12.6f' %0
    vect3='% 5s' %(int(grid[2])+1)+'% 12.6f' %0+'% 12.6f' %0+'% 12.6f' %float(norm[2])
    filout.write(line0+'\n'+line1+'\n'+line2+'\n'+vect1+'\n'+vect2+'\n'+vect3+'\n')
    for i in range(1,len(kind)+1):
        begin='%5.0f' %float(atomic[kind[i]])+'% 12.6f' %float(atomic[kind[i]])
        end='% 12.6f' %(float(coord[i][0])*1.889725985)+'% 12.6f' %(float(coord[i][1])*1.889725985)+'% 12.6f' %(float(coord[i][2])*1.889725985)
        filout.write(begin+end+'\n')
        
    # -----------------------------------------------------------------------
    # WRITE ESP
    # -----------------------------------------------------------------------

    n_loop=floor((int(grid[2])+1)/6)
    i_extra=(int(grid[2])+1)%6

    for x in range(1,int(grid[0])+2):
        for y in range(1,int(grid[1])+2):
            if i_extra==0:
                line=''
                for i in range(1,int(n_loop+1)):
                    for k in range(1,7):
                       line=line+'%13.5E' %float(v[x][y][(i-1)*6+k])
                    filout.write(line+'\n')
                    line=''
            else:
                line=''
                for i in range(1,int(n_loop+1)):
                    for k in range(1,7):
                        line=line+'%13.5E' %float(v[x][y][(i-1)*6+k])
                    filout.write(line+'\n')
                    line=''
                for k in range(1,i_extra+1):
                    line=line+'%13.5E' %float(v[x][y][(i-1)*6+k])
                filout.write(line+'\n')
                line=''


    # -----------------------------------------------------------------------
    # CONCLUSION
    # -----------------------------------------------------------------------

    filin1.close()
    filin2.close()
    filout.close()

    print "Returned file",fulpath_out


#########################################################################################
#########################################################################################
#########################################################################################

    
if flag==2:
    
# -----------------------------------------------------------------------
# VASP routine
# -----------------------------------------------------------------------

    if not os.path.exists('LOCPOT'):               #check if LOCPOT exists
       abort("No LOCPOT could be found.")
    if not os.path.exists('OUTCAR'):               #check for OUTCAR
       abort('No OUTCAR was found')

    filin1=open('LOCPOT','r')             #open files
    filin2=open('OUTCAR','r')

    # -----------------------------------------------------------------------
    # READING LOCPOT DATAS
    # -----------------------------------------------------------------------

    data2=filin1.read().split('\n')

    name=finder(data2[0],'^','$')    #cube file
    filout=open(name+'.cube','w')
    fulpath_out=name+'.cube'

    if float(data2[1])<>1:
        print('WARNING: Scaling factor is not 1: script might not be valid')

    lattice={}
    for i in range(0,3):
        lattice[i]=finder(data2[i+2],'^','$').split(' ')   #lattice vectors

    info=finder(data2[5],'^','$').split(' ')
    n_types=len(info)     #number of different type of atom
    n_atoms=0
    for i in info:
        n_atoms=n_atoms+int(i)  #total number of atoms

    coord_type=finder(data2[6],'^','$')

    coord={}
    coord_resc={}
    for i in range(1,n_atoms+1):
        coord[i]=finder(data2[i+6],'^','$').split(' ')
        coord_resc[i]=[0,0,0]
        if coord_type=='Direct':
            for i_dim in range(0,3):
                coord_resc[i][i_dim]=0
                for j_dim in range(0,3):
                    coord_resc[i][i_dim]=coord_resc[i][i_dim]+float(coord[i][j_dim])*float(lattice[j_dim][i_dim])
        else:
            coord_resc[i]=coord[i]
            print('WARNING: coord type is not Direct: script might not be valid')

    grid=finder(data2[n_atoms+8],'','$').split(' ')     #number of grid points

    v={}
    for x in range(1,int(grid[0])+1):
        v[x]={}
        for y in range(1,int(grid[1])+1):
            v[x][y]={}
#            for z in range(1, int(grid[2])+1):
#                v[x][y][z]=0

    count=1
    for z in range(1,int(grid[2])+1):
        for y in range(1,int(grid[1])+1):
            for x in range(1,int(grid[0])+1):
                line=int(floor(count/5))
                nbr=count%5
                if nbr==0:
                    v[x][y][z]=finder(data2[line+9+n_atoms-1],'^','$').split(' ')[4]
                else:
                    v[x][y][z]=finder(data2[line+9+n_atoms],'^','$').split(' ')[nbr-1]
                count=count+1

    # -----------------------------------------------------------------------
    # READING OUTCAR DATAS
    # -----------------------------------------------------------------------

    data1=filin2.read().split('\n')
    kind={}
    for line in data1:                #read the nature of the unique atoms
        if re.search('INCAR',line):
            n=data1.index(line)
            for i in range(1,n_types+1):
                kind[i]=finder(data1[n+i],'POTCAR:','$').split(' ')[1]

    for i in kind.keys():
        if kind[i]=='Al_h':
            kind[i]='Al'


    # -----------------------------------------------------------------------
    # NEEDED CUBE DATAS
    # -----------------------------------------------------------------------

    origin=[0,0,0]
    cube_lattice={}
    for i in range(0,3):
        cube_lattice[i]=[]
        for j in range(0,3):
            cube_lattice[i].append(float(lattice[i][j])*1.889725985/float(grid[i]))   #step-size (conversion to Bohr)

    factor_e=-1/27.212
    atomic={'H':'1','He':'2','Li':'3','Be':'4','B':'5','C':'6','N':'7','O':'8','F':'9','Ne':'10','Na':'11','Mg':'12','Al':'13','Si':'14','P':'15','S':'16','Cl':'17','Ar':'18','K':'19','Ca':'20','Sc':'21','Ti':'22','V':'23','Cr':'24','Mn':'25','Fe':'26','Co':'27','Ni':'28','Cu':'29','Zn':'30','Ga':'31','Ge':'32','As':'33','Se':'34','Br':'35','Kr':'36','Rb':'37','Sr':'38','Y':'39','Zr':'40','Nb':'41','Mo':'42','Tc':'43','Ru':'44','Rh':'45','Pd':'46','Ag':'47','Cd':'48','In':'49','Sn':'50','Sb':'51','Te':'52','I':'53','Xe':'54','Cs':'55','Ba':'56','La':'57','Ce':'58','Pr':'59','Nd':'60','Pm':'61','Sm':'62','Eu':'63','Gd':'64','Th':'65','Dy':'66','Ho':'67','Er':'68','Tm':'69','Yb':'70','Lu':'71','Hf':'72','Ta':'73','W':'74','Re':'75','Os':'76','Ir':'77','Pt':'78','Au':'79','Hg':'80','Tl':'81','Pb':'82','Bi':'83','Po':'84','At':'85','Rn':'86','Fr':'87','Ra':'88','Ac':'89','Th':'90','Pa':'91','U':'92','Np':'93','Pu':'94','Am':'95','Cm':'96','Bk':'97','Cf':'98','Es':'99','Fm':'100','Md':'101','No':'102','Lw':'103'}  #convert nature of the atom into atomic number

    # -----------------------------------------------------------------------
    # WRITE FIRST LINES + COORD (conversion to Bohr)
    # -----------------------------------------------------------------------

    line0='--------------------REPEAT charges--------------------'
    line1='---cube file created from VASP (LOCPOT and OUTCAR)----'
    line2='% 5s' %n_atoms+'% 12.6f' %origin[0]+'% 12.6f' %origin[1]+'% 12.6f' %origin[2]
    vect1='% 5s' %int(grid[0])+'% 12.6f' %float(cube_lattice[0][0])+'% 12.6f' %float(cube_lattice[0][1])+'% 12.6f' %float(cube_lattice[0][2])
    vect2='% 5s' %int(grid[1])+'% 12.6f' %float(cube_lattice[1][0])+'% 12.6f' %float(cube_lattice[1][1])+'% 12.6f' %float(cube_lattice[1][2])
    vect3='% 5s' %int(grid[2])+'% 12.6f' %float(cube_lattice[2][0])+'% 12.6f' %float(cube_lattice[2][1])+'% 12.6f' %float(cube_lattice[2][2])
    filout.write(line0+'\n'+line1+'\n'+line2+'\n'+vect1+'\n'+vect2+'\n'+vect3+'\n')

    count=1
    for i in range(0, n_types):
        for j in range(0, int(info[i])):
            begin='%5.0f' %float(atomic[kind[i+1]])+'% 12.6f' %float(atomic[kind[i+1]])
            end='% 12.6f' %(float(coord_resc[count][0])*1.889725985)+'% 12.6f' %(float(coord_resc[count][1])*1.889725985)+'% 12.6f' %(float(coord_resc[count][2])*1.889725985)
            count=count+1
            filout.write(begin+end+'\n')


    # -----------------------------------------------------------------------
    # WRITE ESP
    # -----------------------------------------------------------------------

    n_loop=floor(int(grid[2])/6)
    i_extra=int(grid[2])%6
#
#    for x in range(1,int(grid[0])+1):
#        for y in range(1,int(grid[1])+1):
#            if i_extra==0:
#                line=''
#                for i in range(1,int(n_loop+1)):
#                    for k in range(1,7):
#                       line=line+'%13.5E' %(float(v[x][y][(i-1)*6+k])*factor_e)
#                    filout.write(line+'\n')
#                    line=''
#            else:
#                line=''
#                for i in range(1,int(n_loop+1)):
#                    for k in range(1,7):
#                        line=line+'%13.5E' %(float(v[x][y][(i-1)*6+k])*factor_e)
#                    filout.write(line+'\n')
#                    line=''
#                for k in range(1,i_extra+1):
#                    line=line+'%13.5E' %(float(v[x][y][(i-1)*6+k])*factor_e)
#                filout.write(line+'\n')
#                line=''


    for x in range(1,int(grid[0])+1):
        for y in range(1,int(grid[1])+1):
            for z in range(1,int(grid[2])+1):
                filout.write('%13.5E' %(float(v[x][y][z])*factor_e))
                if z%6==0:
                    filout.write('\n')
            filout.write('\n')
               

    # -----------------------------------------------------------------------
    # CONCLUSION
    # -----------------------------------------------------------------------

    filin1.close()
    filin2.close()
    filout.close()

    print "Returned file",fulpath_out


