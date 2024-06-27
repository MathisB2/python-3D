import pygame
from pygame.locals import *
from math import*
from random import*

####################################### FONCTIONS #######################################

def equa_plan(n,p):
    #n(a,b,c) = vecteur normal du plan
    #p(xp,yp,zp) = point appartenant au plan
    # equation de la forme ax + by + cz + d = 0

    a=n[0]
    b=n[1]
    c=n[2]
    xp=p[0]
    yp=p[1]
    zp=p[2]

    d=-(a*xp+b*yp+c*zp)

##    print("L'équation est : {}x + {}y + {}z + {} = 0".format(a,b,c,d))
    return [a,b,c,d]





def equa_droite(v,j):
    #v = vecteur dirrecteur de la droite
    #j = point appartenant à la droite

    #équation de la forme :
    #   x = at + xj
    #   y = bt + yj
    #   z = ct + zj

    a=v[0]
    b=v[1]
    c=v[2]

    xj=j[0]
    yj=j[1]
    zj=j[2]


    x=[a,xj]
    y=[b,yj]
    z=[c,zj]

##    print("L'équation de la droite est :")
##    print("x = {}t + {}".format(x[0],x[1]))
##    print("y = {}t + {}".format(y[0],y[1]))
##    print("z = {}t + {}".format(z[0],z[1]))

    return [x,y,z]



def coor_V(N,O):
    #p est léquation de plan connue

    #on cherche le point V(xv,yv,xn)


    xn=N[0]
    yn=N[1]
    zn=N[2]

    xo=O[0]
    yo=O[1]
    zo=O[2]


    if xn==0 and yn>0:
       xv=-1
       yv=yn
    elif xn==0 and yn<0:
       xv=1
       yv=yn
    elif xn>0 and yn==0:
       xv=xn
       yv=1
    elif xn<0 and yn==0:
       xv=xn
       yv=-1


##    #solution 1
    else :
       xv=(-yo + xn*sqrt(xo**2 - 2*xo*yn + yo**2 - 2*yo*yn + xn**2 + yn**2) + yn)/(sqrt(xo**2 - 2*xo*yn + yo**2 - 2*yo*yn + xn**2 + yn**2))
       yv=yn+((xo-xn)/(sqrt(xo**2 - 2*xo*yn + yo**2 - 2*yo*yn + xn**2 + yn**2)))



    #solution 2
##    else:
##       xv=(yo + xn*sqrt(xo**2 - 2*xo*yn + yo**2 - 2*yo*yn + xn**2 + yn**2) - yn)/(sqrt(xo**2 - 2*xo*yn + yo**2 - 2*yo*yn + xn**2 + yn**2))
##       yv=yn-((xo-xn)/(sqrt(xo**2 - 2*xo*yn + yo**2 - 2*yo*yn + xn**2 + yn**2)))



    zv=zn


    return [xv,yv,zv] #renvoie un point

##    return [xv-xn,yv-yn,zv-zn] #renvoie un vecteur



def coor_U(O,N,V):
    #O=[xo,yo,zo] = coordonnées du point O
    #N=[xn,yn,zn] = coordonnées du point N
    #V=[xv,yv,zv] = coordonnées du point V renvoyées par la fonction coor_V

    #on cherche le point U(xu,yu,zu)

    xn=N[0]
    yn=N[1]
    zn=N[2]

    xv=V[0]
    yv=V[1]
    zv=V[2]

    xo=O[0]
    yo=O[1]
    zo=O[2]




    #solution 1
    xu=( xn * sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2)) +(yn*zv)-(zn*yv)) / (sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2)))

    yu=(-xn*zv+yn*(sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2)))+(zn*xv)) / (sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2)))

    zu=zn+(((xn*yv)-(yn*xv)) / (sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2))))



    #solution 2
##    xu=( xn * sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2)) -(yn*zv)+(zn*yv)) / (sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2)))
##
##    yu=(xn*zv+yn*(sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2)))-(zn*xv)) / (sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2)))
##
##    zu=zn-(((xn*yv)-(yn*xv)) / (sqrt((xn**2)*(zv**2)+(xn**2)*(yv**2)-(2*xn*yn*xv*yv)-(2*xn*zn*xv*zv)+(yn**2)*(xv**2)+(yn**2)*(zv**2)-(2*yn*zn*zv*yv)+(zn**2)*(xv**2)+(zn**2)*(yv**2))))



    #utiliser la sollution 1 pour z_nu>0


    return [xu,yu,zu] #renvoie un point



def coor_M(zoom,O,N):
    #zoom = distance MN
    #O=[xo,yo,zo] = coordonnées du point O
    #N=[xn,yn,zn] = coordonnées du point N

    xn=N[0]
    yn=N[1]
    zn=N[2]

    xo=O[0]
    yo=O[1]
    zo=O[2]

    NO=vecteur(N,O)



    a=zoom/(sqrt((xo-xn)**2+(yo-yn)**2+(zo-zn)**2))     #coefficient tel que a*||NO||=zoom

    xm=xn-a*NO[0]
    ym=yn-a*NO[1]
    zm=zn-a*NO[2]
    M=[xm,ym,zm]

    return M



def rotation(a,b,N,O):
    #O=[xo,yo,zo] = coordonnées du point O
    #N=[xn,yn,zn] = coordonnées du point N

    xn=N[0]
    yn=N[1]
    zn=N[2]

    xo=O[0]
    yo=O[1]
    zo=O[2]

    dist1=sqrt((xo-xn)**2+(yo-yn)**2+(zo-zn)**2)

    vect_NO=vecteur(N,O)
    plan=equa_plan(vect_NO,N)

    V=coor_V(N,O)
    U=coor_U(O,N,V)

    xu=U[0]
    yu=U[1]
    zu=U[2]

    xv=V[0]
    yv=V[1]
    zv=V[2]

    a/=-60
    b/=-60
##    c=(a*xn-a*xv+b*xn-b*xu+xo-xn)/(xo-xn)
    c=0

    xn+=a*(xv-xn)+b*(xu-xn)+c*(xo-xn)
    yn+=a*(yv-yn)+b*(yu-yn)+c*(yo-yn)
    zn+=a*(zv-zn)+b*(zu-zn)+c*(zo-zn)

    dist2=sqrt((xo-xn)**2+(yo-yn)**2+(zo-zn)**2)
    while dist2-dist1>0.00001:
        c=0.0001
        xn+=c*(xo-xn)
        yn+=c*(yo-yn)
        zn+=c*(zo-zn)
        dist2=sqrt((xo-xn)**2+(yo-yn)**2+(zo-zn)**2)


    return [xn,yn,zn]








def lim_plan(n,u,v):
    #renvoie les coordonnées des quatre coins du plan (A,B,C,D)
    #n=[xn,yn,zn] = coordonnées du point N
    #u=[xu,yu,zu] = coordonnées du point U
    #v=[xv,yv,zv] = coordonnées du point V

    xn=n[0]
    yn=n[1]
    zn=n[2]

    xu=u[0]
    yu=u[1]
    zu=u[2]

    xv=v[0]
    yv=v[1]
    zv=v[2]

    #point I
    xa=-0.8*(xv-xn)+0.45*(xu-xn)+xn
    ya=-0.8*(yv-yn)+0.45*(yu-yn)+yn
    za=-0.8*(zv-zn)+0.45*(zu-zn)+zn

    #point J
    xb=+0.8*(xv-xn)+0.45*(xu-xn)+xn
    yb=+0.8*(yv-yn)+0.45*(yu-yn)+yn
    zb=+0.8*(zv-zn)+0.45*(zu-zn)+zn

    #point J
    xc=+0.8*(xv-xn)-0.45*(xu-xn)+xn
    yc=+0.8*(yv-yn)-0.45*(yu-yn)+yn
    zc=+0.8*(zv-zn)-0.45*(zu-zn)+zn

    #point L
    xd=-0.8*(xv-xn)-0.45*(xu-xn)+xn
    yd=-0.8*(yv-yn)-0.45*(yu-yn)+yn
    zd=-0.8*(zv-zn)-0.45*(zu-zn)+zn

    return [[xa,ya,za],[xb,yb,zb],[xc,yc,zc],[xd,yd,zd]]


def coor_2D(inter,n,u,v):
    #renvoie les coordonnées en 2D du point inter dans le plan (N,V,U)
    #N(0,0), A(-0.8,0.45), B(-0.8,0.45), C(0.8,-0.45), D(-0.8,-0.45)
    #inter=[xi,yi,zi] = coordonnées du point d'intersection I
    #n=[xn,yn,zn] = coordonnées du point N
    #u=[xu,yu,zu] = coordonnées du point U
    #v=[xv,yv,zv] = coordonnées du point V

    #on cherche a et b tels que NI est une combinaison linéaire de NV et de NU


    xn=n[0] #d
    yn=n[1] #e
    zn=n[2] #f

    xu=u[0] #g
    yu=u[1] #h
    zu=u[2] #i

    xv=v[0] #j
    yv=v[1] #k
    zv=v[2] #l

    xi=inter[0] #m
    yi=inter[1] #n
    zi=inter[1] #o


    a=(-((xn*(-xn*yv +xn*yi +yn*xv -yn*xi -xv*yi +yv*xi)) / (xn*yu -xn*yv -yn*xu +yn*xv +xu*yv -yu*xv)) +((xu*(-xn*yv +xn*yi +yn*xv -yn*xi -xv*yi +yv*xi)) / (xn*yu -xn*yv -yn*xu +yn*xv +xu*yv -yu*xv)) + xn -xi) / (xn-xv)
    b=(-xn*yv +xn*yi +yn*xv -yn*xi -xv*yi +yv*xi) / (xn*yu -xn*yv -yn*xu +yn*xv +xu*yv -yu*xv)

    return [a,b]








def intersection(droite,plan):
    #p = équation cartésienne de plan, de la forme ax + by +cz + d =0
    #p = (a,b,c,d)
    #d = équation paramétriquer de droite, de la forme :
    #   x= e*t + xj
    #   y= f*t + yj
    #   z= g*t + zj
    #g = [(e,xj),(f,yj),(g,zj)]

    #on cherche le point d'intesection i(x,y,z) de la droite d avec le plan p


    a=plan[0]
    b=plan[1]
    c=plan[2]
    d=plan[3]


    e=droite[0][0]
    xj=droite[0][1]
    f=droite[1][0]
    yj=droite[1][1]
    g=droite[2][0]
    zj=droite[2][1]


    t=(-a*xj-b*yj-c*zj-d)/(a*e+b*f+c*g)

    x=e*t+xj
    y=f*t+yj
    z=g*t+zj

    i=[x,y,z]

    return i


def vecteur(depart,arrivee):
    #depart=[xd,yd,zd]=point de départ du vecteur
    #arrive=[xa,ya,za]=point d'arrivée du vecteur
    xd=depart[0]
    yd=depart[1]
    zd=depart[2]

    xa=arrivee[0]
    ya=arrivee[1]
    za=arrivee[2]

    return [xa-xd,ya-yd,za-zd]


def verif_pos(J,N,O,U,V):
    #N=[xn,yn,zn] = coordonnées du point N
    #O=[xo,yo,zo] = coordonnées du point O

    #V=[xv,yv,zv] = coordonnées du point V
    #U=[xu,yu,zu] = coordonnées du point U
    #J=[xj,yj,zj] = coordonnées du point testé


    d=xn=N[0] #d
    x=yn=N[1] #x
    f=zn=N[2] #f

    g=xu=U[0] #g
    h=yu=U[1] #h
    y=zu=U[2] #y

    z=xv=V[0] #z
    k=yv=V[1] #k
    l=zv=V[2] #l

    m=xj=J[0] #m
    n=yj=J[1] #n
    o=zj=J[1] #o

    p=xo=O[0] #p
    q=yo=O[1] #q
    r=zo=O[1] #r

    vect_JN=vecteur(J,N)


##    m-d=a*(p-d)+b*(z-d)+c*(g-d)
##    n-x=a*(q-x)+b*(k-x)+c*(h-x)
##    o-f=a*(r-f)+b*(l-f)+c*(y-f)

##    a=(xn*yu*zv-xn*yu*zj+xn*yv*zu-xn*zv*yj+xn*yj*zu+zn*xu*yv-zn*xu*yj+zn*yu*xj-zn*yu*xv-zn*yv*xj+zn*yj*xv-xu*yv*zj + xu*zv*yj-xu*zv*yn+xu*zj*yn-yu*zv*xj+yu*zj*xv+yv*xj*zu+zv*xj*yn-xj*yn*zu-yj*zu*xv-zj*yn*xv+yn*zu*xv)/(xn*yu*zv-xn*yu*zo+xn*yv*zo-xn*yv*zu-xn*zv*yo+xn*yo*zu+zn*xu*yv-zn*xu*yo+zn*yu*xo-zn*yu*xv-zn*yv*xo+zn*yo*xv-xu*yv*zo+xu*zv*yo - xu*zv*yn+xu*zo*yn-yu*zv*xo+yu*zo*xv+yv*xo*zu+zv*xo*yn-xo*yn*zu-yo*zu*xv-zo*yn*xv+yn*zu*xv)


    a2=(d*h*l-d*h*o+d*k*o-d*k*y-d*l*n+d*n*y+f*k*g-f*g*n+f*h*m-f*h*z-f*k*m+f*n*z-g*k*o+g*l*n-g*l*x+g*o*x-h*l*m+h*o*z+k*m*y+l*m*x-m*x*y-n*y*z-o*x*z+x*y*z) / (d*h*l-d*h*r+d*k*r-d*k*y-d*l*q+d*q*y+f*g*k-f*g*q+f*h*p-f*h*z-f*k*p+f*q*z-g*k*r+g*l*q-g*l*x+g*r*x-h*l*p+h*r*z+k*p*y+l*p*x-p*x*y-q*y*z-r*x*z+x*y*z)


    if a2<0:
        return False    #renvoie True si le point est devant la caméra, False sinon
    else:
        return True






def refresh(N,O):
    #N=[xn,yn,zn] = coordonnées du point N
    #O=[xo,yo,zo] = coordonnées du point O

    #V=[xv,yv,zv] = coordonnées du point V
    #U=[xu,yu,zu] = coordonnées du point U


    xn=N[0]
    yn=N[1]
    zn=N[2]

    xo=O[0]
    yo=O[1]
    zo=O[2]

    vect_NO=vecteur(N,O)
    plan=equa_plan(vect_NO,N)

    V=coor_V(N,O)
    U=coor_U(O,N,V)

    M=coor_M(zoom,O,N)


    #calculer les limites du plan
    face=lim_plan(N,U,V)

    lim_x1=coor_2D(intersection(equa_droite(vect_NO,face[0]),plan),N,U,V)[0]*echelle+screen_l/2
    lim_x2=coor_2D(intersection(equa_droite(vect_NO,face[1]),plan),N,U,V)[0]*echelle+screen_l/2
    lim_y1=coor_2D(intersection(equa_droite(vect_NO,face[2]),plan),N,U,V)[1]*echelle+screen_h/2
    lim_y2=coor_2D(intersection(equa_droite(vect_NO,face[1]),plan),N,U,V)[1]*echelle+screen_h/2


    faces=[]
    #on regroupe toutes les faces dans une seule liste
    for i in range(len(scene)):
        #i = un objet de la scene
        for j in range(len(scene[i])):
            #scene[i][j] = une face de l'objet

            #déterminer si un face doit etre affichée

            R=0
            for k in range(len(scene[i][j])):
                J=scene[i][j][k]
                if verif_pos(J,N,O,U,V)==False:
                   R+=1
            if R==0:
               faces.append(scene[i][j])     #à modifier !



    points_moy=[]
    fenetre.fill(col[bg_num])



    #on détermine l'ordre d'affichage des faces
    i=0
    for face in faces:
        if fill==True:
            x_moy=0
            y_moy=0
            z_moy=0
            n=0
            for l in face:
               x_moy+=l[0]
               y_moy+=l[1]
               z_moy+=l[2]
               n+=1
            x_moy/=n
            y_moy/=n
            z_moy/=n

            dist=sqrt((xn-x_moy)**2+(yn-y_moy)**2+(zn-z_moy)**2)
            points_moy.append(dist)
        else:
            points_moy.append(i)
            i+=1
    ordre=tri(points_moy)


    for i in range(len(faces)):
        if fill==True:
            face=faces[ordre[i]]
        else:
            face=faces[i]

        points_2D=[]

        for k in face:
            #k = un point de la face
            #on détermine les coordonnées 2D
            vect=vecteur(M,k)
            droite=equa_droite(vect,k)
            I=intersection(droite,plan)
            points_2D.append(coor_2D(I,N,U,V))

##            print(points_2D)
        inter=0
        if fill==True:
            #on affiche les faces pleines
            points_fill=[]
            for k in range(len(points_2D)):
                ax=screen_l/2+(points_2D[k][0]*echelle)
                ay=screen_h/2+(points_2D[k][1]*echelle)
                if lim_x1<=ax<=lim_x2 and lim_y1<=ay<=lim_y2:
                   inter+=1
                points_fill.append([ax,ay])
            if inter>0:
               pygame.draw.polygon(fenetre, col[fill_num], points_fill)


        for k in range(len(points_2D)):
            #on affiche les arretes

            dx=screen_l/2+(points_2D[k-1][0]*echelle)
            dy=screen_h/2+(points_2D[k-1][1]*echelle)
            ax=screen_l/2+(points_2D[k][0]*echelle)
            ay=screen_h/2+(points_2D[k][1]*echelle)

            if (lim_x1<=ax<=lim_x2 and lim_y1<=ay<=lim_y2) or (lim_x1<=dx<=lim_x2 and lim_y1<=dy<=lim_y2):
                   inter+=1

                   pygame.draw.line(fenetre, col[edge_num], (dx,dy), (ax,ay), edge_width)



    #afficher le plan
    face=lim_plan(N,U,V)
    points_2D=[]

    for k in face:
        #k = un point de la face
        #on détermine les coordonnées 2D
        droite=equa_droite(vect_NO,k)
        I=intersection(droite,plan)
        points_2D.append(coor_2D(I,N,U,V))


    for k in range(4):
        dx=screen_l/2+(points_2D[k-1][0]*echelle)
        dy=screen_h/2+(points_2D[k-1][1]*echelle)
        ax=screen_l/2+(points_2D[k][0]*echelle)
        ay=screen_h/2+(points_2D[k][1]*echelle)
        pygame.draw.line(fenetre, col[lim_num], (dx,dy), (ax,ay), 4)



    pygame.display.flip()




def rech_max(L):
    #Recherche la valeur maximale de la liste L et retourne son indice
    if len(L)>0:
        m=max(L)
        for i in range(len(L)):
            if L[i]==m:
               return i
    else:
        return 0


def tri(dist):
    #dist=liste des distances des faces
    indices=[]
    for i in range(len(dist)):
        n=rech_max(dist)
        dist[n]=min(dist)-1
        indices.append(n)

    return indices







####################################### VARIABLES #######################################


pygame.init()
screen_h,screen_l=810,1440

fenetre=pygame.display.set_mode((screen_l,screen_h))

mode=1

blanc=(255,255,255)
gris=(184,192,202)

noir_dark=(0,0,0)
noir=(0,22,54)
noir75=(59, 77, 104)
violet=(151,107,244)
bleu=(50,144,255)
vert=(84,205,161)
jaune=(248,197,67)
rouge=(220,105,105)

col=[blanc,gris,noir75,noir,noir_dark,violet,bleu,vert,jaune,rouge]





#objets de la scene
A=(-0.5,-0.5,-0.5)
B=(0.5,-0.5,-0.5)
C=(0.5,0.5,-0.5)
D=(-0.5,0.5,-0.5)
E=(-0.5,-0.5,0.5)
F=(0.5,-0.5,0.5)
G=(0.5,0.5,0.5)
H=(-0.5,0.5,0.5)

face1=[A,B,F,E]
face2=[B,C,G,F]
face3=[C,G,H,D]
face4=[D,A,E,H]
face5=[A,B,C,D]
face6=[E,F,G,H]

cube=[face1,face2,face3,face4,face5,face6]


A=(2.5,-0.5,-0.5)
B=(3.5,-0.5,-0.5)
C=(3.5,0.5,-0.5)
D=(2.5,0.5,-0.5)
E=(2.5,-0.5,0.5)
F=(3.5,-0.5,0.5)
G=(3.5,0.5,0.5)
H=(2.5,0.5,0.5)

face1=[A,B,F,E]
face2=[B,C,G,F]
face3=[C,G,H,D]
face4=[D,A,E,H]
face5=[A,B,C,D]
face6=[E,F,G,H]

cube2=[face1,face2,face3,face4,face5,face6]







scene=[cube,cube2]



#points essentiels

point_N=[3,2,1]   #attention ! au moins une coordonnée doit etre différente de 0
point_O=[0,0,0]
vecteur_NO=vecteur(point_N,point_O)

zoom=1    #distance MN
point_M=coor_M(zoom,point_O,point_N)

zoom_fac=0.1   #sensibilité du zoom

echelle=900   #échelle d'affichage


fill=True
fill_num=1  #numéro de couleur de remplissage : 1 = gris
edge_num=9 #couleur des arretes : 9 = rouge
lim_num=9  #couleur de la limite : 9 = rouge
bg_num=2  #couleur de fond : 3 = noir75

edge_width=2
lim_width=4


##plan=equa_plan(vecteur_NO,point_N)
##
##
##
##print('point N :',point_N)
##point_V=coor_V(point_N,point_O)
##
##print('point V :',point_V)
##
##
##point_U=coor_U(point_O,point_N,point_V)
##print('point U :',point_U)
##
##
##
##point_A=lim_plan(point_N,point_U,point_V)[0]
##print("point A :",point_A)
##
##point_B=lim_plan(point_N,point_U,point_V)[1]
##print("point B :",point_B)
##
##point_C=lim_plan(point_N,point_U,point_V)[2]
##print("point C :",point_C)
##
##
##
##point_D=lim_plan(point_N,point_U,point_V)[3]
##print("point D :",point_D)
##
##
##print("coor 2D :",coor_2D(point_N,point_N,point_U,point_V))

print("-------------- Contrôles --------------")
print()
print("- Clic gauche pour tourner la caméra")
##print("- Clic droit pour dépalcer la caméra")
print("- Mollette pour raprocher la caméra")
print("- Ctrl + Mollette pour augmenter la distance focale")
print("- Maj + Mollette pour augmenter l'échelle d'afficchage")
print("- Tab pour changer le mode de vue (plein ou fillaire)")
print("- [1] pour changer la couleur de remplissage")
print("- [2] pour changer la couleur des arrêtes")
print("- [3] pour changer la couleur de fond")
print("- [4] pour changer la couleur des limites de la caméra")







####################################### PROGRAMME #######################################

refresh(point_N,point_O)

while mode!=-1:
    for event in pygame.event.get():



##      croix de la fenetre
        if event.type == QUIT:
            mode = -1

##      rotation
        if event.type ==MOUSEBUTTONDOWN and event.button == 1:
            press=True
            ox=event.pos[0]
            oy=event.pos[1]
            while press:
                for event in pygame.event.get():
                    if event.type ==MOUSEBUTTONUP and event.button == 1:
                       press=False

                    try:
                        ax=event.pos[0]
                        ay=event.pos[1]

                        a=ax-ox
                        b=ay-oy

                        point_N=rotation(a,b,point_N,point_O)
                        vecteur_NO=vecteur(point_N,point_O)

                        refresh(point_N,point_O)

                        ox=event.pos[0]
                        oy=event.pos[1]
                    except:
                        pass

##      zoom : distance NO
        if event.type ==MOUSEBUTTONDOWN and event.button == 4:
            print("zoom +")
            dx=vecteur_NO[0]*zoom_fac
            dy=vecteur_NO[1]*zoom_fac
            dz=vecteur_NO[2]*zoom_fac



            xn=point_N[0]+dx
            yn=point_N[1]+dy
            zn=point_N[2]+dz

            point_N=[xn,yn,zn]


            vecteur_NO=vecteur(point_N,point_O)

            refresh(point_N,point_O)




        if event.type ==MOUSEBUTTONDOWN and event.button == 5:
            print("zoom -")
            dx=vecteur_NO[0]*zoom_fac
            dy=vecteur_NO[1]*zoom_fac
            dz=vecteur_NO[2]*zoom_fac

            xn=point_N[0]-dx
            yn=point_N[1]-dy
            zn=point_N[2]-dz

            point_N=[xn,yn,zn]

            vecteur_NO=vecteur(point_N,point_O)



            refresh(point_N,point_O)

##      échelle d'affichage
        if event.type==KEYDOWN and event.key==K_LSHIFT:
           press=True
           while press:
                for event in pygame.event.get():
                    if event.type==KEYUP and event.key==K_LSHIFT:
                       press=False

                    if event.type ==MOUSEBUTTONDOWN and event.button == 4:
                        print("échelle +")
                        echelle+=20

                        refresh(point_N,point_O)




                    if event.type ==MOUSEBUTTONDOWN and event.button == 5:
                        print("échelle -")
                        echelle-=20

                        refresh(point_N,point_O)


##      FOV : distance MN
        if event.type==KEYDOWN and event.key==K_LCTRL:
           press=True
           while press:
                for event in pygame.event.get():
                    if event.type==KEYUP and event.key==K_LCTRL:
                       press=False

                    if event.type ==MOUSEBUTTONDOWN and event.button == 4:
                        print("FOV +")
                        zoom=1.1*zoom

                        refresh(point_N,point_O)




                    if event.type ==MOUSEBUTTONDOWN and event.button == 5:
                        print("FOV -")
                        zoom=0.9*zoom

                        refresh(point_N,point_O)



##      chager le mode d'affichage : fillaire ou rempli
        if event.type==KEYDOWN and event.key==K_TAB:
           fill=not fill
           refresh(point_N,point_O)





##      réinitialiser la vue
        if event.type==KEYDOWN and event.key==K_SPACE:
            point_N=[3,2,1]   #attention ! au moins une coordonnée doit etre différente de 0
            point_O=[0,0,0]
            vecteur_NO=vecteur(point_N,point_O)
            zoom=1    #distance MN
            point_M=coor_M(zoom,point_O,point_N)

            zoom_fac=0.1   #sensibilité du zoom

            echelle=900   #échelle d'affichage

            fill=True
            fill_col=noir

            refresh(point_N,point_O)


##      changer la couler de remplissage
        if event.type==KEYDOWN and event.key==K_1:
            if fill_num==len(col)-1:
               fill_num=0
            else:
                fill_num+=1
            refresh(point_N,point_O)



##      changer la couler des arretes
        if event.type==KEYDOWN and event.key==K_2:
            if edge_num==len(col)-1:
               edge_num=0
            else:
                edge_num+=1

            refresh(point_N,point_O)


##      changer la couler de fond
        if event.type==KEYDOWN and event.key==K_3:
            if bg_num==len(col)-1:
               bg_num=0
            else:
                bg_num+=1
            refresh(point_N,point_O)



##      changer la couler de la limite
        if event.type==KEYDOWN and event.key==K_4:
            if lim_num==len(col)-1:
               lim_num=0
            else:
                lim_num+=1

            refresh(point_N,point_O)




pygame.quit()