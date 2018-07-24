

# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 17:52:33 2017

@author: Multielio
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 18:00:02 2017

@author: Multielio
"""
import math
import time
import functools
import numpy as np
from random import *
from visual import *
from visual.controls import *
from visual.graph import *
import wx
import wx.propgrid

__version__ = 6.1
nscene = None
gvel = 1
gw = None
gel = 0.4
gely = 0.4
gr = 0.01225
clocksp = 0.0001
distance = 0.067
colorb = [color.red,color.yellow,color.green,color.orange,color.blue,color.cyan]


#######################################################################
# Manip coord
#############################################################################


def fx2(vo,t,p):
    xo,yo,zo= p
    vxo,vyo,vzo = vo
    return (vxo*t)+xo
def fy2(vo,t,p):
    xo,yo,zo= p
    vxo,vyo,vzo = vo
    return (-1.0/2.0)*9.81*(t**2)+(vyo*t)+yo

def fvx2(vo,t,p):
    xo,yo,zo= p
    vxo,vyo,vzo = vo
    return vxo
def fvy2(vo,t,p):
    xo,yo,zo= p
    vxo,vyo,vzo = vo
    return (-9.81*t)+vyo

def traj(radius,xplan,vo,wo,p,dt,elas,time_gen_afimpact,debug=False):
    global distance
    xo,yo,zo = p
    vxo,vyo,vzo = vo
    b,y = fictif_wcollide(xplan,radius,vo,p)
    print("Ypos when collide = {}".format(y))
    pos = []
    vel = []
    fictif_timpact = (xplan-radius-xo)/vxo # moment à la balle touche une ligne fictive en x=xplan-radius (Formule obtenue via PFD)
    timeimpact = (xplan-radius-xo)/vxo 
    pimpact = None #Position du centre de la balle lors de l'impact
    impactloc = None #Point de contact entre la balle et le bord lors de l'impact
    velimpact = None
    if not(b): # Si la balle tape dans les y<=0
        #On doit calculer la traj de t=0 t=impact (timpact obtenu par PFD)
        pos.extend([(fx2(vo,t,p),fy2(vo,t,p),0) for t in np.linspace(0,fictif_timpact,fictif_timpact/dt)]) #Génération des trajectoires
        vel.extend([(fvx2(vo,t,p),fvy2(vo,t,p),0) for t in np.linspace(0,fictif_timpact,fictif_timpact/dt)])
        pimpact = (xplan-radius,y,0)
        impactloc = (xplan,y,0)
        velimpact = (vxo,(-9.81*fictif_timpact)+vyo,0)
    else:
        timeimpact,boo = dicosolver(fictif_timpact,((xplan-xo)/vxo),vxo,vyo,xplan,radius,0.0001*distance/5,xo,yo)
        if boo:
            print("Temps impact (cercle): {} s".format(timeimpact))
            pos.extend([(fx2(vo,t,p),fy2(vo,t,p),0) for t in np.linspace(0,timeimpact,timeimpact/dt)]) #Génération des trajectoires
            vel.extend([(fvx2(vo,t,p),fvy2(vo,t,p),0) for t in np.linspace(0,timeimpact,timeimpact/dt)])
            pimpact = pos[-1]
            impactloc = (xplan,0,0) 
            velimpact = (vxo,(-9.81*timeimpact)+vyo,0)
        else:
            print("Impossible d'effectuer la dichotomie !")
            pos.extend([(fx2(vo,t,p),fy2(vo,t,p),0) for t in np.linspace(0,time_gen_afimpact,(time_gen_afimpact)/dt)])
            vel.extend([(fvx2(vo,t,p),fvy2(vo,t,p),0) for t in np.linspace(0,time_gen_afimpact,(time_gen_afimpact)/dt)])
            return pos,vel,None,0,0
    iimpact = len(pos) #On stocke le nombre de positions que l'on a calculé
    if pimpact != None and impactloc != None and velimpact != None:
        px,py,pz = pimpact
        x,y,z= impactloc
        vix,viy,viz = velimpact
        ex,ey = elas 
        ### CALCUL DU TETA PLAN D'IMPACT 
        tetaplan = np.pi/2 
        if xplan-radius-px != 0:
            tetaplan = np.pi/2 -np.arccos((xplan-px)/radius) # Calcul de l'angle entre le plan d'impact et l'horizontale
        print("Teta plan d'incidence: "+str((tetaplan*180)/np.pi)+" degres" ) # On affiche cela dans la console
        #On calcule le plan d'impact à afficher
        plan = None
        if tetaplan != np.pi/2:
            plan = planhorizontal_teta(distance,distance*5.7/5,(-1)*distance/5,distance/5,30,tetaplan)   
        ###########################################
        ### PROJECTION VERS LE PLAN D'IMPACT
        vixp, viyp= transfo1(vix,viy,tetaplan)
        #####
        S = (radius * wo)/vixp # On calcule le paramètre de spin 
        alpha = 2.0/5.0
        vx2 = ( ((1-alpha*ex)/(1+alpha))  + ((alpha*(1+ex)*S)/(1+alpha)) )*vixp #On utilise l'équation donné par l'article
        neww = ((alpha-ex)/(1+alpha) + (1+ex)/((1+alpha)*S))*wo # On utilise l'équation donnant la nouvelle vitesse de rotation donnée par l'article
        vy2 = -viyp*ey 
        #### PROJECTION VERS LE PLAN CARTE
        vx2t,vy2t = transfo1(vx2,vy2,-tetaplan)
        ####### Nouvelles vitesses dans le repère initial
        newvel = (vx2t,vy2t,0)
        #{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}#
        
        # On calcule maintenant la traj après impact.
        pos.extend([(fx2(newvel,t,pimpact),fy2(newvel,t,pimpact),0) for t in np.linspace(0,time_gen_afimpact,time_gen_afimpact/dt)])
        vel.extend([(fvx2(newvel,t,pimpact),fvy2(newvel,t,pimpact),0) for t in np.linspace(0,time_gen_afimpact,time_gen_afimpact/dt)])
        if debug:
            file =open("debug.txt", "w") 
            file.write(str(pos))
            file.close()
            print("pimpact",pimpact)
        return pos,vel,plan,timeimpact,iimpact

def fictif_wcollide(xplan,radius,vo,p):
    #Attention la fonction marche uniquement quand x0=y0=0
    # Renvoie true si la balle va collide sur l'arc de cercle. + le y si y<0
    xo,yo,zo = p
    vxo,vyo,vzo =vo
    y = fy2(vo,(xplan-radius-xo)/vxo,p)
    return (y >0,y)
def u(t,vxo,vyo,xplan,xo,yo):
    return ((-9.81/2)*(t**2)+vyo*t+yo)**2 + (vxo*t+xo-xplan)**2  # Rayon du cercle de centre (xplan,0) passant par le point de la trajectoire en t
def dicosolver(t1,t2,vxo,vyo,xplan,radius,pres,xo,yo):
    umid = u((t2+t1)/2,vxo,vyo,xplan,xo,yo)
    u1 = u(t1,vxo,vyo,xplan,xo,yo)
    u2 = u(t2,vxo,vyo,xplan,xo,yo)
    print("t1: {} / t2: {} / umid: {} / u1: {} / u2: {} / radiuscarre: {} ".format(t1,t2,umid,u1,u2,radius**2))
    if radius**2< u2 or radius**2 > u1:
        return (t1,false)
    if u1<u2:
        u1,u2 = u2,u1
    
    if abs(umid-radius**2) < pres:
        return ((t2+t1)/2,true)
    if radius**2 <umid and u2 <radius**2 :
        return dicosolver((t2+t1)/2,t2,vxo,vyo,xplan,radius,pres,xo,yo)
    if radius**2>=umid and radius**2<=u1:
        return dicosolver(t1,(t2+t1)/2,vxo,vyo,xplan,radius,pres,xo,yo)
        
        
def planvertical(ys,ye,zs,ze,x,pres):
    return [(x,y,z) for y in np.linspace(ye,ys,pres) for z in np.linspace(ze,zs,pres)]
    
def planhorizontal(xs,xe,zs,ze,y,pres):
    return [(x,y,z) for x in np.linspace(xe,xs,pres-10) for z in np.linspace(ze,zs,pres)]

def planhorizontal_teta(xs,xe,zs,ze,pres,teta):
    a= [(x,(x-xs)*np.tan(teta),z) for x in np.linspace(xe,xs,pres-10) for z in np.linspace(ze,zs,pres)]
    return a
def transfo(compo, teta_plan_dinci):
    return (compo*np.cos(teta_plan_dinci),compo*np.sin(teta_plan_dinci))
def transfo1(compovx,compovy, teta_plan_dinci):
    return ((compovx*np.cos(teta_plan_dinci)+compovy*np.sin(teta_plan_dinci)),(-compovx*np.sin(teta_plan_dinci)+compovy*np.cos(teta_plan_dinci)))


#####################################################################
#EXTRAIT D'UNE API
#############################################################################
def axes( frame, colour, sz, posn ): 
    directions = [vector(sz,0,0), vector(0,sz,0), vector(0,0,sz)]
    texts = ["X","Y","Z"]
    posn = vector(posn)
    for i in range (3):
       curve( frame = frame, color = colour, pos= [ posn, posn+directions[i]])
       label( frame = frame,color = colour,  text = texts[i], pos = posn+ directions[i],opacity = 0, box = False )
#####################################################################
#OUTILS
######################################################################

def random_color():
    global colorb
    if colorb ==[]:
        colorb = [color.red,color.yellow,color.green,color.orange,color.blue,color.cyan]
    g = colorb[random.randint(0,len(colorb))]  
    colorb.remove(g)
    return g



#######################################################################
# FONCTIONS APPELLE PAR LES BOUTONS
#############################################################################
def start(evt):
    global gvel
    global gel
    global gw
    global gr
    global gely
    global nscene
    global clocksp
    global distance
    nscene.time = 0
    if gw == None:
        gw = float(gvel)/float(gr)
    Sphe = Sphere(colorss= random_color(), pos=(0,float(gr),0),radius = float(gr), vel = (float(gvel),0,0),elasticity=(float(gel),float(gely)),mass =0.028, w = float(gw))
    nscene.add_obj(Sphe)
    Sphe.VpyCreate()
    Sphe.Gen_traj(distance,clocksp,0.4)
def call(evt):
    global nscene
    nscene.time = 0
    nscene.func = gcurve(color=random_color())
    nscene.time = 0
    for i in nscene.objects:
        if i.objtype == 0:
            i.Removeobj(true)
            nscene.remove_obj(i)
        if i.objtype == 1:
            if not(i.solid):
                i.Removeobj()
                nscene.remove_obj(i)
def stopp(evt):
    global nscene
    nscene.time = 0
    nscene.func = gcurve(color=random_color())
    nscene.time = 0
    for i in nscene.objects:
        if i.objtype == 0:
            i.Removeobj(false)
            nscene.remove_obj(i)
        
def delplan(evt):
    global nscene
    for i in nscene.objects:
        if i.objtype == 1:
            if not(i.solid):
                i.Removeobj()
                nscene.remove_obj(i)
def setvel(evt):
    global gvel
    global nscene
    gvel = nscene.bvel.GetValue()
    print("Nouvelle vitesse: "+gvel+" m/s")
  
def setel(evt):
    global gel
    global nscene
    gel = nscene.bel.GetValue()
    print("Nouvelle elasticite x: "+str(gel))
 
def setely(evt):
    global gely
    global nscene
    gely = nscene.bely.GetValue()
    print("Nouvelle elasticite y: "+str(gely))
    
def setw(evt):
    global gw
    global nscene
    gw = nscene.bw.GetValue()
    print("Nouvelle rotation: "+str(gw)+" rad/s")
def setr(evt):
    global gr
    global nscene
    gr = nscene.rw.GetValue()
    print("Nouveau rayon: "+str(gr))
#############################################################################
# OBJETS 
#############################################################################
class Scene(object):
    """
    L'object scène permet de gérer l'ensemble des objets ainsi que l'écoulement du temps + affichage des boutons
    """
    def __init__(self, time=0,clockspeed=0.0001, objects = []):
        self.time = float(time) #Temps 
        self.objects = objects # Objets attachés à la scène
        self.clockspeed = clockspeed # C'est le dt ajouté à chaque uptade de la scene
        self.win = None # Les variables ci dessous s'occupent de rendre les objets facilement accessibles une fois créé
        self.bvel = None
        self.bel = None
        self.bely = None
        self.bw = None
        self.rw = None
        self.bline = None
        self.func = None
        self.text_clock = None #
       
    def add_obj(self,ob)   :
        self.objects.append(ob) # On ajoute un objet à la liste qui stoke les objets associés à la scene
    def exist_obj(self,ob):
        # On vérifie grace à l'unicité des id pour chaque objet l'existance d'un objets dans la scene
        for i in self.objects:
            if i == ob.ids :
                return true
        return false
    def remove_obj(self,ob):
        if self.exist_obj(ob):
            for i in self.objects:
                if i.ids == ob.ids:
                    self.objects.remove(i)
    def update_scene(self):
        # A chaque fois que cette fonction est appellée tout les objets associés à la scene sont update et le temps
        # est avancé de dt
        self.time += self.clockspeed #On avance le temps de dt
        text1 = " {} s".format(self.time) #On met à jour le temps sur l'interface utilisateur
        self.text_clock.SetLabel(text1) #
        flag = False
        for i in self.objects : #On parcours tout les objets attachés à la scène
            if not(flag) :
                if i.objtype ==0: # L'id 0 correspond aux sphères, on s'occupe donc d'actualiser uniquement les sphères
                    if i.exist == True: # On vérifie par mesure de sécurité si la sphère existe encore, elle aurait pû etre détruite.
                        flag = True #Permet de pas faire bug la courbes des vitesse en temps réel en affichant uniquement la vitesse de la dernière sphère ajoutée
                        self.func.plot( pos=(self.time, i.getNormeVel()) ) #On affiche la courbe des vitesse
            i.update_object(self.clockspeed,self) # Même si le flag est déclenché on update l'objet grace à la fonction update associé à l'objet

    def start(self):
        global distance
        # Fonction faite pour générer la fenetre, les boutons etc (On remplit aussi les variables initialisées sur None)
        
        #Ci dessous, les variables permettent de modifier l'allure du menu très rapidement/ Espaces entre les différents boutons etc
        taille = 850
        ecart_h = 100
        ecart_v1 = 35
        ecart_v = 70
        
        block1_hi = 25
        block2_hi = 280
        text_ecart = 20
        #
        
        self.win =window(menus=False, title="Simulateur de Collisions V6.1", x=0, y=0, width=1150, height=800) # Création de la fenêtre
        scene= display(window=self.win, x=50, y=30, width=taille, height=450) # Creation de la scène à l'interieur de la fenêtre
        scene.autoscale = scene.autocenter = False 
        scene.center = (distance+(2*distance/10),distance/5,(3*distance)/5)
        gdisplay(window=self.win,x=50, y=500, width=taille,title='Evolution de la vitesse en fonction du temps', xtitle='t(s)', ytitle='v(m/s)', height=200) # On génère la courbe des vitesses en temps réel
        self.func = gcurve(color=color.cyan)
      
        p = self.win.panel 
        axes( None, color.white, distance/3, (0,0,0)) #Création d'un repère XYZ pour que l'on voit mieux ce qu'il se passe
        
        
        #########################################
        ### BOUTON
        ########################################
        
        text_clock = "{} s".format(0)
        self.text_clock = wx.StaticText(self.win.panel, pos=(960,185L), label=text_clock )
        
        deca_h = ecart_h + taille 
        left = wx.Button(p, label='START', pos=(deca_h,block1_hi))
        left.Bind(wx.EVT_BUTTON, start)
        
        left = wx.Button(p, label='STOP', pos=(deca_h,block1_hi+ecart_v1))
        left.Bind(wx.EVT_BUTTON, stopp)
        left = wx.Button(p, label='DEL PLAN', pos=(deca_h,block1_hi+2*ecart_v1))
        left.Bind(wx.EVT_BUTTON, delplan)
        left = wx.Button(p, label='CLEAR ALL', pos=(deca_h,block1_hi+3*ecart_v1))
        left.Bind(wx.EVT_BUTTON, call)
        
        wx.StaticText(p, pos=(deca_h,block2_hi), label="Vel en m/s:" )
        self.bvel = wx.TextCtrl(p, pos=(deca_h,block2_hi+text_ecart), value='',size=(100,20),style=wx.TE_PROCESS_ENTER)
        self.bvel.Bind(wx.EVT_TEXT_ENTER, setvel)
        
        wx.StaticText(p, pos=(deca_h,block2_hi+ecart_v), label="Elasticite tang (Ex):" )
        self.bel = wx.TextCtrl(p, pos=(deca_h,block2_hi+ecart_v+text_ecart), value='',size=(100,20),style=wx.TE_PROCESS_ENTER)
        self.bel.Bind(wx.EVT_TEXT_ENTER, setel)
        
        wx.StaticText(p, pos=(deca_h,block2_hi+2*ecart_v), label="Elasticite verti (Ey):" )
        self.bely = wx.TextCtrl(p, pos=(deca_h,block2_hi+2*ecart_v+text_ecart), value='',size=(100,20),style=wx.TE_PROCESS_ENTER)
        self.bely.Bind(wx.EVT_TEXT_ENTER, setely)
        
        wx.StaticText(p, pos=(deca_h,block2_hi+3*ecart_v), label="W en rad/s:" )
        self.bw = wx.TextCtrl(p, pos=(deca_h,block2_hi+3*ecart_v+text_ecart), value='',size=(100,20),style=wx.TE_PROCESS_ENTER)
        self.bw.Bind(wx.EVT_TEXT_ENTER, setw)
        
        wx.StaticText(p, pos=(deca_h,block2_hi+4*ecart_v), label="Rayon (def: 0.2):" )
        self.rw = wx.TextCtrl(p, pos=(deca_h,block2_hi+4*ecart_v+text_ecart), value='',size=(100,20),style=wx.TE_PROCESS_ENTER)
        self.rw.Bind(wx.EVT_TEXT_ENTER, setr)
        self.bline = wx.StaticLine(p,1,style= wx.LI_VERTICAL,pos=(0,0))
        self.bline.SetSize((5,1060))
        self.bline = wx.StaticLine(p,1,style= wx.LI_VERTICAL,pos=(deca_h-30,0))
        self.bline.SetSize((4,1060))
        self.bline = wx.StaticLine(p,1,style= wx.LI_VERTICAL,pos=(deca_h+130,0))
        self.bline.SetSize((4,1060))
 
        
####################################################################################################
#################################################################################################### 
####################################################################################################
class Solid_Rectangle(object):
    """
        Pour créer des solides rectangulaires
    """
    def __init__(self,objtype =1,pos=(0,0,0),points=[]):
        u = random.randint(0,99999999) #Génération aléatoire de l'id que l'objet va avoir
        self.pos = pos 
        self.points = points #Ensemble des points associés au solide
        self.objtype = objtype # Les solid rectangle sont des objets de type 1 
        self.obj = None
        self.ids = u
        self.colorr = color.white #Couleur des points du solide
        self.solid = True # Définit la consistance du solide : les autres objets peuvent t-ils passer à travers ?
        self.exist= True # Existance de l'objet, utile pour détruire les objets
    def set_not_solid(self):
        self.solid = False
        
    def set_color(self,c):
        self.colorr = c
    def Removeobj(self):
        if self.exist:
            print("[INFO] Deleting obj / TYPE PLAN  /"+ str(self.ids))
            self.exist = False
            self.obj.visible =False
            del(self.obj)
    def add_points(self,p):
        self.points.append(p)
    def VpyCreate(self):
        self.obj =points(color= self.colorr, pos= self.points,size=0.01)
    def update_object(self,dt,scene):
        pass

####################################################################################################
####################################################################################################    
####################################################################################################      
class Sphere(object):
    """
        Objet sphère, permettant de faire bouger la sphère + gestion des collisions
    """
    def __init__(self,objtype =0,ids=0, pos=(0,0,0), colorss= color.orange, radius =1, vel=(0,0,0), elasticity=(1,1),mass = 1,w = 0,forces=[], obj = None):
        u = random.randint(0,99999999) #Génération de l'id de manière aléatoire, pas envie de créer une variable qui s'incrémente
        
        self.objtype = objtype #Les sphères sont des objets de type 0 (choisi arbitrairement)
        self.exist = True
        self.ids = u
        self.nbrcoll = 0 # Nombre de contacts avec un autre objet après le premier contact avec celui ci
        self.colorss = colorss
        
        self.traj = None
        self.traji = 0
        self.pos = pos #Pos est un triplet (x,y,z) qui donne la position de l'objet (est mis à jour par update())
        self.radius = radius # Rayon de la sphère
        self.vel = vel  # Triplet qui donne la vitesse selon les trois composantes
        self.elasticity = elasticity   #Couple d'elasticitée
        self.forces = forces # N'est plus utilisé (je dois le supprimer un jour :) 
        self.mass = mass # Masse de la sphère
        self.w = w # Rotation de la sphère
        
        
        self.obj = obj # Texture de associé à l'objet
        self.text1 = None #Boutons / TEXTES associés à l'objet
        self.text2 = None
        self.collid = False #
        print("[INFO] Created obj / TYPE SPHERE / ID : "+str(self.ids)+" / ELAS: "+str(self.elasticity)+" / VEL: "+str(self.vel)) # On affiche ce que l'on créer dans la console histoire de pouvoir debug 
    def Gen_traj(self,xplan,dt,time_gen_afimpact):
        print("Self.pos = {}".format(self.pos))
        self.traj = traj(self.radius,xplan,self.vel,self.w,self.pos,dt,self.elasticity,time_gen_afimpact,debug=True)
        
    def VpyCreate(self):
        # Création de la texture / animation visuelle associé à l'objet
        self.obj = sphere(pos=self.pos, radius=self.radius,color=self.colorss, material = materials.BlueMarble, make_trail=True ) #On génère l'objet texturé et on le stocke histoire de pouvoir le faire bouger sur l'ecran de l'utilisateur plus tard
    def TxtCreator(self,scene):
        # Quelques textes affichés à l'utilisateur pour savoir ce qu'il se passe précisement
        text = "XYZ = ({},{},{})".format(0,0,0)
        self.text1 = wx.StaticText(scene.win.panel, pos=(960,210L), label=text )
        text2 = "VXYZ = ({},{},{})".format(0,0,0)
        self.text2 = wx.StaticText(scene.win.panel, pos=(960,235L), label=text2 )
    def Removeobj(self,brut):
        # Fonction de destruction de l'objet
        if self.exist:
            print("[INFO] Deleting obj / TYPE SPHERE / "+ str(self.ids)+ " / Nbr coll: "+str(self.nbrcoll))
            self.exist = False
            self.obj.visible =False
            if brut:
                self.obj.make_trail = False
                self.text1.Hide()
                self.text2.Hide()
            del(self.obj) 
    def getNormeVel(self):
        # Juste pour récuperer la norme de la vitesse
        vx,vy,vz = self.vel
        return math.sqrt(vx**2 + vy**2 + vz**2)
    
    def display_plan(self,plan):
        if plan != None:
            NSr = Solid_Rectangle(points= plan) #Affichage purement esthétique du plan d'impact
            NSr.set_color(color.cyan)
            NSr.set_not_solid() 
            NSr.VpyCreate()
            nscene.add_obj(NSr)
     

    def update_object(self,dt,scene):
        #Mise à jour de l'objet : on recup les paramètres de l'objet et on l
        if not(self.exist):
            return None
        try:
            if(self.text1 == None):
                self.TxtCreator(scene) #
            pos,vel,plan,timeimpact,iimpact = self.traj
            if self.traji == iimpact:
                self.display_plan(plan)
                print("Position centre de masse à l'impact: {}".format(pos[self.traji]))
            self.pos = pos[self.traji]
            self.vel =vel[self.traji]
            self.traji += 1
            self.VpyUpdate()
            x,y,z = self.pos
            vx,vy,vz = self.vel
            text1 = "XYZ = ({},{},{})".format(round(x,3),round(y,3),round(z,3))
            self.text1.SetLabel(text1)
            text2 = "VXYZ = ({},{},{})".format(round(vx,2),round(vy,2),round(z,2))
            self.text2.SetLabel(text2)
        except Exception as e:
            print(str(e))
    def VpyUpdate(self):
        try:
            x,y,z = self.pos
            self.obj.pos = (x,y,z)
            self.obj.rotate(angle=self.w, origin=self.pos)
        except Exception as e:
            print(str(e))
        
####################################################################################################
# UTILISATION DES OBJETS
#################################################################################################### 
def startsimul():
    global nscene
    global distance
    nscene = Scene(clockspeed=clocksp)
    nscene.start()
    Sr = Solid_Rectangle(points= planvertical(((-2)*distance)/5,0,(-1)*distance/5,distance/5,distance,50)+planhorizontal(distance,5.2*distance/5,(-1)*distance/5,distance/5,0,50))
    Sr.VpyCreate()
    nscene.add_obj(Sr)
    print("------------------------ \n Simulateur de Collisions \n   v. %r \n ------------------------ \n" %__version__)
    while 1 :
        rate(100)
        nscene.update_scene()
startsimul()
