def reatancia(f,tipo,valor):
    if tipo == -1:
        x_c = 1/(2*(np.pi)*f*valor)
        return -1j*x_c
    if tipo == 1:
        x_l = 2*(np.pi)*f*valor
        return 1j*x_l

def polar(z):
    return (round(np.absolute(z),3),round(np.angle(z,deg=True),3),'polar')

def div_pol(z1,z2):
    return (np.absolute(z1)/np.absolute(z2),np.angle(z1,deg=True)-np.angle(z2,deg=True),'polar')

def imp_mod_ang(mod,fp,tipo):
    if tipo == 1:
        return np.round(op_rot(mod,(180/np.pi)*np.arccos(fp)),3)
    if tipo == -1:
        return np.round(op_rot(mod,-(180/np.pi)*np.arccos(fp)),3)    

def retangular(z):
    ang_rad_comp = (np.pi/180)*z[1]*1j
    return z[0]*np.exp(ang_rad_comp)

def cap_crrg_fp(s,v,f,fp_i,fp_n):
    ang_i = np.arccos(fp_i)
    ang_f = np.arccos(fp_n)
    p = s*fp_i
    delta_q = p*(np.tan(ang_i) - np.tan(ang_f))
    vel_ang = 2*np.pi*f
    return (delta_q)/((v**2)*vel_ang)

def zeq_p(list_imp,pol=0):
    s_p = 0
    for i in list_imp:
        s_p += 1/i
    if pol == 1:
        return polar(1/s_p)
    return np.round(1/s_p,3)

def zeq_s(list_imp,pol=0):
    s_p = 0
    for i in list_imp:
        s_p += i
    if pol == 1:
        return polar(s_p)
    return s_p

def op_rot(modulo,angulo):
    ang_rad_comp = 1j*((np.pi/180)*angulo)
    return np.round(modulo*np.exp(ang_rad_comp),3)

def delta_y(ra,rb,rc):
    r1 = (rb*rc)/(ra+rb+rc)
    r2 = (ra*rc)/(ra+rb+rc)
    r3 = (ra*rb)/(ra+rb+rc)
    print(f'r1 = {r1}')
    print(f'r2 = {r2}')
    print(f'r3 = {r3}')
    return [r1,r2,r3]

def y_delta(r1,r2,r3):
    ra = (r1*r2 + r2*r3 + r1*r3)/r1
    rb = (r1*r2 + r2*r3 + r1*r3)/r2
    rc = (r1*r2 + r2*r3 + r1*r3)/r3
    print(f'ra = {ra}')
    print(f'rb = {rb}')
    print(f'rc = {rc}')
    return [ra,rb,rc]   

def delta_3fv(vp,seq,fase_ref,prt = True):
    a = op_rot(1,120)
    vf = op_rot(vp,fase_ref)
    if seq == 1:
        vab = vf*(a**0)
        vbc = vf*(a**2)
        vca = vf*(a**1)
        if prt == False:
            return[vab,vbc,vca]
        print(f'vab = {polar(vab)}')
        print(f'vbc = {polar(vbc)}')
        print(f'vca = {polar(vca)}')
        return[vab,vbc,vca]
    if seq == -1:
        vab = vf*(a**0)
        vbc = vf*(a**1)
        vca = vf*(a**2)
        if prt == False:
            return[vab,vbc,vca]
        print(f'vab = {polar(vab)}')
        print(f'vbc = {polar(vbc)}')
        print(f'vca = {polar(vca)}')
        return[vab,vbc,vca]
    
def y_3fv(vp,seq,fase_ref,prt=True):
    a = op_rot(1,120)
    b = op_rot(np.sqrt(3),30)
    vf = op_rot(vp,fase_ref)
    if seq == 1:
        van,vbn,vcn = vf*(a**0),vf*(a**2),vf*(a**1)
        vab,vbc,vca = b*van,b*vbn,b*vcn
        if prt == False:
            return [van,vbn,vcn],[vab,vbc,vca]
        print(f'van = {polar(van)}| vbn = {polar(vbn)}| vcn = {polar(vcn)}')
        print(f'vab = {polar(vab)}| vbc = {polar(vbc)}| vca = {polar(vca)}')
        return [van,vbn,vcn],[vab,vbc,vca]
    if seq == -1:
        van,vbn,vcn = vf*(a**0),vf*(a**1),vf*(a**2)
        vab,vbc,vca = b*van,b*vbn,b*vcn
        if prt == False:
            return [van,vbn,vcn],[vab,vbc,vca]
        print(f'van = {polar(van)}| vbn = {polar(vbn)}| vcn = {polar(vcn)}')
        print(f'vab = {polar(vab)}| vbc = {polar(vbc)}| vca = {polar(vca)}')
        return [van,vbn,vcn],[vab,vbc,vca]
    
import numpy as np

def delta_3fi(vp,seq,fase_ref,z,prt=True):
    v3f = delta_3fv(vp,seq,fase_ref,False)
    iab,ibc,ica = v3f[0]/z,v3f[1]/z,v3f[2]/z
    ia,ib,ic = iab - ica, ibc - iab, ica - ibc
    if prt == False:
        return [iab,ibc,ica],[ia,ib,ic]
    print(f'iab = {polar(iab)}| ibc = {polar(ibc)}| ica = {polar(ica)}')
    print(f'ia = {polar(ia)}| ib = {polar(ib)}| ic = {polar(ic)}')
    return [iab,ibc,ica],[ia,ib,ic]

def y_3fi(vp,seq,fase_ref,z,prt=True):
    v3f = y_3fv(vp,seq,fase_ref,False)
    ia,ib,ic = v3f[0][0]/z,v3f[0][1]/z,v3f[0][2]/z
    i_n = ia + ib + ic
    if prt == False:
        return [ia,ib,ic,i_n]
    print(f'ia = {polar(ia)}| ib = {polar(ib)}| ic = {polar(ic)}| in = {polar(i_n)}')
    return [ia,ib,ic,i_n]

def ref_trafo(valor,grandeza,n1,n2,lado,pol=1):
    n = pol*n2/n1
    if lado == 'p':
        if grandeza == 'v':
            return valor/n
        if grandeza == 'i':
            return valor*n
        if grandeza == 'z':
            return valor/(n**2)
    if lado == 's':
        if grandeza == 'v':
            return valor*n
        if grandeza == 'i':
            return valor/n
        if grandeza == 'z':
            return valor*(n**2)   
        
def rend_trafo(n1,n2,lado,rp,xp,rs,xs,rc,xm,z_l,pot,fp=1):
    s_carga = pot/fp, 
    v_carga = np.sqrt(s_carga)
    pass
