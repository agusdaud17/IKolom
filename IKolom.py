#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 22:01:15 2022

metric mm, MPa, kN, kNm
imperial inch, psi, kips, kip_ft

@author: agus daud
"""
# %config InlineBackend.figure_formats = ['svg']
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from math import pi

import shapely.geometry as sg
from shapely.geometry import MultiLineString
import shapely.ops as so



bar = {'#3':0.375,
       '#4':0.500,
       '#5':0.625,
       '#6':0.750,
       '#7':0.875,
       '#8':1.000,
       '#9':1.128,
       '#10':1.270,
       '#11':1.410,
       '#14':1.693,
       '#18':2.257}
       
class Rect:
    def __init__(self,fy=420,fc=25,Es=200000,ey=0.002,small_axial = False,Rf=1):

        self.unit = 'metric'
        self.Fy = fy
        self.fc = fc
        self.Es = Es
        self.ey = ey
        self.small_axial = small_axial
        self.Rf = Rf

    def units(self,unit):
        if unit == 'imperial':
            self.unit = unit
            self.fc = self.fc #ksi
            self.Fy = self.Fy #ksi
            self.Es = self.Es #ksi
            self.Dm = bar[f'{self.long_bar_no}']
            self.Dv = bar[f'{self.stir_bar_no}']
            self.unit_force = "kips"
            self.unit_moment = "kips-ft"
            self.unit_length = "in"
            self.unit_area = "in^2"
        else:
            self.unit_force = "kips"
            self.unit_moment = "kips-ft"
            self.unit_length = "in"
            self.unit_area = "in^2"
            pass
    
    def dimension(self,b=400,h=400,ds=40,n=0,nx=4,ny=4,Dm=19,Dv=10,long_bar_no='#6',stir_bar_no='#3',digit_round=2,factor=1):
        self.b = b
        self.h = h
        self.c = [b,h] # b == arah x; dan h == arah y
        self.ds = ds
        if n > 0: self.n = [int((n/4) +1), int((n/4) +1)]
        else: self.n = [nx,ny]
        self.Dm = Dm
        self.Dv = Dv
        self.long_bar_no = long_bar_no
        self.stir_bar_no = stir_bar_no
        self.digit_round = digit_round
        self.factor = factor

        # check if units in metric or imperial
        if type(self.Dm) == str:
            self.long_bar_no = self.Dm
            self.stir_bar_no = self.Dv
            self.units('imperial')
            self.to_kip_ft = 1/12/1000
        else:
            self.units('metric')
            self.to_kip_ft = 1/1000000

        # Area of single bar
        self.Atul = round(0.25*pi*self.Dm**2,2)

        # spacing (s), and effectif height(d1) for each side
        self.s = []
        self.d1 = []

        for i in range(2):
            self.s.append((self.c[i] - (2*(self.ds+self.Dv) + self.Dm)) / (self.n[i]-1))
            self.d1.append(self.c[i] - (self.ds + self.Dv + self.Dm/2))

        self.nt = 2*((self.n[0]-2) + (self.n[1]))
        self.Ast = round(self.nt*self.Atul,self.digit_round)
        self.Ag = round(self.b*self.h,self.digit_round)
        self.rho = round(self.Ast/self.Ag,5)

    def reduction_factor(self,es1):
        if abs(es1)<=self.ey: phi = 0.65
        elif abs(es1) < self.ey+0.003: phi = 0.65 + (0.25*(abs(es1)-self.ey)/(0.005-self.ey))
        else: phi = 0.9
        return(phi)

    def beta1(self):
        if self.unit == 'imperial':
            if self.fc <= 4000: beta1 = 0.85
            elif self.fc < 8000: beta1 = 0.85 - (0.05*((self.fc-4000)/1000))
            else: beta1 = 0.65
        else:
            if self.fc <= 28: beta1 = 0.85
            elif self.fc < 55: beta1 = 0.85 - (0.05*((self.fc-28)/7))
            else: beta1 = 0.65
        return(round(beta1,2))


    def table_result(self,direction):
        # side 1 = arah_x; side 2 = arah_y
        if direction == 'x':
            self.direction = direction 
            x = 1
            y = 0
        elif direction == 'y':
            self.direction = direction
            x = 0
            y = 1

        self.Po = round(((0.85*self.fc*(self.Ag-self.Ast)) + (self.Fy*self.Ast)) / 1000,self.digit_round)
        self.phi_Po = round(0.65*self.Po, self.digit_round)
        self.phi_Pn_max = round(0.8*self.phi_Po,self.digit_round)

        list_Z = [0.5]
        list_phi = [0.65]
        list_Pn = []
        list_Mn = [0]
        list_phi_Mn = [0]
        list_phi_Pn_unf = [self.phi_Po]
        list_phi_Pn = [self.phi_Pn_max] # true value; had been cut off!!
        Pmin = round((0.10*self.fc*self.Ag)/1000,self.digit_round) # minimum value ro increase value of phi
        Z = 0.5 # Nilai Z awal = 0.5

        beta1 = self.beta1()

        while True:
            # untuk row pertama
            es1 = Z*self.ey
            c = (0.003 / (0.003 - es1)) * self.d1[y]
            
            if abs(es1) >= self.ey:
                fs1 = (es1/abs(es1)) * self.Fy
            else:
                fs1 = es1*self.Es

            a = beta1*c

            if a > self.c[y]:
                a == self.c[y]
            
            if a >= self.d1[y]:
                Fs1 = (fs1 - (0.85*self.fc))*self.Atul*self.n[x]
            else:
                Fs1 = fs1*self.Atul*self.n[x]

            Cc = 0.85*self.fc*a*self.c[x]
            # Cc = (0.85*self.fc*(a*self.b - self.Ast)) # yang pertama sebelum diedit
            Mn = (Cc*(self.c[y]/2 - a/2)) + (Fs1*(self.c[y]/2 - self.d1[y]))
            Fs = Cc + Fs1
            
            for row in range(2,self.n[y]+1,1):
                d = self.d1[y] - ((row-1)*self.s[y]) # jarak tulangan dari sisi tinjauan
                es = 0.003 * (c-d) / c
                if abs(es) >= self.ey:
                    fs = (es/abs(es)) * self.Fy
                else:
                    fs = es*self.Es
                
                # jika row terakhir
                if row == self.n[y]:
                    As = self.Atul*self.n[x]
                else:
                    As = self.Atul*2
                
                if a >= d:
                    Fs_add = (fs - (0.85*self.fc))*As
                else:
                    Fs_add = fs*As
                    
                Fs += Fs_add
                
                Mn += Fs_add*(self.c[y]/2 - d)
            
            Pn = Fs

            #considering small axial force ACI 318-99 Point 9.3.2.2
            if self.small_axial ==True:
                if Z == -1:
                    self.Pb = round(Pn/1000,self.digit_round)
                else:
                    self.Pb = 10000000000

                self.Pmin = min(Pmin,self.Pb)

                # faktor reduksi
                if Pn/1000 < self.Pmin:
                    phi = 0.9 - ((Pn/1000)*(0.9-0.65)/self.Pmin)
                else:
                    phi = 0.65
            else:
                phi = self.reduction_factor(es1)

            if Mn <= -100000 or Pn <= -100000:
                pass
            else:
                Mn = round(Mn*self.to_kip_ft,self.digit_round)
                Pn = round(Pn/1000,self.digit_round)
                phi_Mn = phi*Mn
                phi_Pn = phi*Pn

                list_Z.append(Z)
                list_phi.append(phi)
                list_Pn.append(Pn)
                list_Mn.append(Mn)
                list_phi_Mn.append(phi_Mn)
                list_phi_Pn_unf.append(phi_Pn)

                # cut for phi_Pn_max
                if phi_Pn >= self.phi_Pn_max:
                    list_phi_Pn.append(self.phi_Pn_max)
                else:
                    list_phi_Pn.append(phi_Pn)

            if Z <= -20:
                break
                
            Z -= 0.02

        # retun table format Z, phi, Pn, Mn
        result = zip(list_Z,list_phi,list_Pn,list_Mn,list_phi_Pn_unf,list_phi_Mn,list_phi_Pn)

        return(result)

    def plot(self,direction='y',grid=False,title="Interaction Diagram",figsize = (5,5),**kwargs):
        table = list(self.table_result(direction))
        self.__df = pd.DataFrame(table, columns = ['Z','phi','Pn','Mn','phi Pn uncut','phi Mn','phi Pn'])

        if self.unit == 'imperial':
            fac = 1/25.4
            xlabel = "$\phi M_n$ (kip-ft)"
            ylabel = "$\phi P_n$ (kips)"
        else:
            fac = 1
            xlabel = "$\phi M_n$ (kNm)"
            ylabel = "$\phi P_n$ (kN)"

        plt.figure(figsize=figsize)
        ax = plt.gca()
        self.__df.plot(x='Mn',y='Pn',ax=ax,style='--', label=r'$P_n - M_n$', color='black')
        self.__df.plot(x='phi Mn',y='phi Pn uncut',ax=ax, style='--',color='grey')
        self.__df.plot(x='phi Mn',y='phi Pn',ax=ax,label=r"$\phi P_n - \phi M_n$",color='black')

        if 'data' in kwargs:
            plot_value = kwargs['data']
            self.DCRatio = []
            for i, val in enumerate(plot_value):
                Mn = round(self.interp_Mn(Pu=val[0]),self.digit_round)
                R = round(val[1]/(self.Rf*Mn),self.digit_round)
                self.DCRatio.append(R)
                if R>1: color_plot = "red"
                else: color_plot = "blue"

                x = val[1]
                y = val[0]

                plt.plot(x,y,'x',color=color_plot)
                label = f'Pu={y}\nMu={x}'
                # plt.annotate(label,(x,y),textcoords="offset points",xytext=(-5,0),ha='right',fontsize=8,color=color_plot)
                plt.annotate(f"({i+1})",(x,y),textcoords="offset points",xytext=(0,5),ha='center',fontsize=8,color=color_plot)

                x, y = figsize 
                plt.text(0.18*(x+0.2), 0.15*(y+(-i)), f"Plot-{i+1}\nPu = {val[0]} {self.unit_force}\nMu = {val[1]} {self.unit_moment}\nphi_Mn = {Mn} {self.unit_moment}\nD/C Ratio = {R}", fontsize=8, color=color_plot, transform=plt.gcf().transFigure)
        
        # plot section control
        plot_loc_z = [0,-0.5,-1]
        for z in plot_loc_z:  
            val = self.interp_z(z=z,x_val='Z',df=self.__df,y1_val='Mn',y2_val='Pn')
            x,y = (val[1],val[0])
            plt.plot([0,x],[0,y],'--',color='grey',marker='.',linewidth=0.5,markersize=10)
            
            sym = 'fy'
            if z == 0: sym = ''
            elif z == -1: z = '-'

            label = f'fs={z}{sym}'
            plt.annotate(label,(x,y),textcoords="offset points",xytext=(-10,0),ha='right',fontsize=8)

        self.list_Mpr = []
        if 'to_Mpr' in kwargs:
            list_data = kwargs['to_Mpr']
            for Pu in list_data:
                Mpr = self.interp_Mpr(Pu=Pu)[1]
                self.list_Mpr.append([Pu,Mpr])
                x,y = (Mpr,Pu)
                plt.plot([Mpr],[Pu],marker='+',color='red',markersize=10)
                label = f'Pu={Pu}\nMpr={Mpr}'
                plt.annotate(label,(x,y),textcoords="offset points",xytext=(10,-5),ha='left',fontsize=8)

        plt.xlabel(f"{xlabel}")
        plt.ylabel(f"{ylabel}")
        # plt.legend(bbox_to_anchor = (1.05, 0.6))
        plt.grid(grid)
        plt.xlim(xmin=0.0)
        plt.ylim(ymin=0.0)
        plt.title(title,fontweight = 'bold')
        plt.show()


    def vertical_dimension(self,start,end,dim,scale_dim, axs,fac):
        # Arrow start
        x_arrow_start = [start[0]-(scale_dim*10*fac), start[0], start[0]+(scale_dim*10*fac)]
        y_arrow_start = [start[1]+(scale_dim*10*fac), start[1], start[1]+(scale_dim*10*fac)]
        # Arrow end
        x_arrow_end = [end[0]-(scale_dim*10*fac), end[0], end[0]+(scale_dim*10*fac)]
        y_arrow_end = [end[1]-(scale_dim*10*fac), end[1], end[1]-(scale_dim*10*fac)]
        # Dimension line
        x1 = [start[0],end[0]]
        y1 = [start[1],end[1]]
        
        axs.fill(x1,y1,color='black',linewidth=1)
        axs.fill(x_arrow_start,y_arrow_start,color='black',linewidth=1)
        axs.fill(x_arrow_end,y_arrow_end,color='black',linewidth=1)
        axs.text(x1[1]-50*fac, y1[1]/2, dim, rotation=90,va='center',fontsize = 10)

    def horizontal_dimension(self,start,end,dim,scale_dim, axs,fac):
        # Arrow start
        x_arrow_start = [start[0]+(scale_dim*10*fac), start[0], start[0]+(scale_dim*10*fac)]
        y_arrow_start = [start[1]+(scale_dim*10*fac), start[1], start[1]-(scale_dim*10*fac)]
        # Arrow end
        x_arrow_end = [end[0]-(scale_dim*10*fac), end[0], end[0]-(scale_dim*10*fac)]
        y_arrow_end = [end[1]-(scale_dim*10*fac), end[1], end[1]+(scale_dim*10*fac)]
        # Dimension line
        x1 = [start[0],end[0]]
        y1 = [start[1],end[1]]
        
        axs.fill(x1,y1,color='black',linewidth=1)
        axs.fill(x_arrow_start,y_arrow_start,color='black',linewidth=1)
        axs.fill(x_arrow_end,y_arrow_end,color='black',linewidth=1)
        axs.text(x1[1]/2, y1[1]-40*fac, dim, ha='center',fontsize = 10)
            
    def interp_Mn(self,Pu=None):
    	func = interp1d(self.__df['phi Pn'],self.__df['phi Mn'])
    	Mn = func(Pu)
    	return(round(Mn.tolist(),self.digit_round))

    def interp_Pn(self,Mu=None):
        func = interp1d(self.__df['phi Mn'],self.__df['phi Pn'])
        Pn = func(Mu)
        return(round(Pn.tolist(),self.digit_round))

    def interp_Mpr(self,Pu=None):
        func = interp1d(self.__df['Pn'],self.__df['Mn'])
        Mpr = func(Pu)
        return(Pu,round(Mpr.tolist(),self.digit_round))
        
    def interp_z(self,z,df,x_val,y1_val,y2_val):
        func_Mn = interp1d(df[x_val],df[y1_val])
        func_Pn = interp1d(df[x_val],df[y2_val])
        Mn = func_Mn(z)
        Pn = func_Pn(z)
        # return(round(Pn.tolist(),2),round(Mn.tolist(),2))
        return(round(Pn.tolist(),self.digit_round),round(Mn.tolist(),self.digit_round))

    def control_value(self):

        print("Control value:")
        print(f'phi Pn_max = {self.phi_Pn_max} {self.unit_force}')
        zero_fy = self.interp_z(z=0,x_val='Z',df=self.__df,y1_val='phi Mn',y2_val='phi Pn')
        half_fy = self.interp_z(z=-0.5,x_val='Z',df=self.__df,y1_val='phi Mn',y2_val='phi Pn')
        min_fy = self.interp_z(z=-1,x_val='Z',df=self.__df,y1_val='phi Mn',y2_val='phi Pn')
        print("----fs = 0:----")
        print(f'phi Pn = {zero_fy[0]} {self.unit_force}')
        print(f'phi Mn = {zero_fy[1]} {self.unit_moment}')

        print("----fs=-0.5fy:----")
        print(f'phi Pn = {half_fy[0]} {self.unit_force}')
        print(f'phi Mn = {half_fy[1]} {self.unit_moment}')

        print("----fs=-fy:----")
        print(f'phi Pn = {min_fy[0]} {self.unit_force}')
        print(f'phi Mn = {min_fy[1]} {self.unit_moment}')

        pure_bending = self.interp_Mn(Pu=0)
        print(f'Pure bending; Mn = {pure_bending} {self.unit_moment}')

        control_data = [self.phi_Pn_max,zero_fy[0],zero_fy[1],half_fy[0],half_fy[1],min_fy[0],min_fy[1],pure_bending]

    #---------------------------------------------------------------------------------------------
    def section(self,scale = 0.8,scale_dim = 0.8):
        # side 1 = arah_x; side 2 = arah_y
        if self.direction == 'x':
            xx = 1
            yy = 0
        elif self.direction == 'y':
            xx = 0
            yy = 1

        # self.analyze()
        fig, axs = plt.subplots()
        axs.set_aspect('equal', 'datalim')

        # Outer
        outer = sg.box(self.b*scale,self.h*scale,0,0)
        x, y = outer.exterior.xy
        axs.fill(x, y, alpha=0.7, fc='white', ec='black')

        d_ef = (self.ds + self.Dv + self.Dm/2)

        # Bottom side
        for i in range(0,self.n[xx]):
            x_coord = (d_ef+(self.s[xx]*i))*scale
            y_coord = d_ef*scale
            bar = sg.Point(x_coord, y_coord).buffer(self.Dm*scale/2)
            x, y = bar.exterior.xy
            axs.fill(x, y, alpha=1, fc='lightgrey', ec='black')
        
        # Left side
        for i in range(0,self.n[yy]):
            x_coord = d_ef*scale
            y_coord = (d_ef+(self.s[yy]*i))*scale
            bar = sg.Point(x_coord, y_coord).buffer(self.Dm*scale/2)
            x, y = bar.exterior.xy
            axs.fill(x, y, alpha=1, fc='lightgrey', ec='black')

        # Right side
        for i in range(0,self.n[yy]):
            x_coord = (self.b-d_ef)*scale
            y_coord = (d_ef+(self.s[yy]*i))*scale
            bar = sg.Point(x_coord, y_coord).buffer(self.Dm*scale/2)
            x, y = bar.exterior.xy
            axs.fill(x, y, alpha=1, fc='lightgrey', ec='black')
        
        # Bottom side
        for i in range(0,self.n[xx]):
            x_coord = (d_ef+(self.s[xx]*i))*scale
            y_coord = (self.h-d_ef)*scale
            bar = sg.Point(x_coord, y_coord).buffer(self.Dm*scale/2)
            x, y = bar.exterior.xy
            axs.fill(x, y, alpha=1, fc='lightgrey', ec='black')

        if self.unit == 'imperial':
            dim_unit = " inch"
            fac = 1/25.4
            bar = self.long_bar_no
        else:
            dim_unit = " mm"
            fac = 1
            bar = f'D{self.Dm}'

        # Dimensi tinggi
        gap = 40*fac
        start = [(0-gap)*scale, 0]
        end = [(0-gap)*scale, self.h*scale]
        # plt.plot(start, end,color='black')
        self.vertical_dimension(start,end,f'{self.h}{dim_unit}',scale_dim,axs,fac)

        # Dimensi lebar
        start = [0, (0-gap)*scale]
        end = [self.b*scale, (0-gap)*scale]
        self.horizontal_dimension(start,end,f'{self.b}{dim_unit}',scale_dim,axs,fac)

        axs.text(self.b/2*scale, (self.h/2 + gap)*scale, f'{self.nt}{bar}', ha='center',fontsize = 11)
        rasio = round(self.rho*100,3)
        axs.text(self.b/2*scale, (self.h/2 - gap)*scale, f'R={rasio}%', ha='center',fontsize = 11)

        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        ax.axis("off")
        plt.show()

    def properties(self):
        if self.unit == 'imperial':
            length = 'inch'
            stress = 'psi'
            axial = 'kips'
            moment = 'kips-ft'
        else:
            length = 'mm'
            stress = 'MPa'
            axial = 'kN'
            moment = 'kNm'

        print("Material:")
        print(f"fc = {self.fc} {stress}")
        print(f"fy = {self.Fy} {stress}")
        print(f"Es = {self.Es} {stress}")
        print("\n")
        print("Dimension:")
        print(f"b = {self.b} {length}")
        print(f"h = {self.h} {length}")
        print(f"ds = {self.ds} {length}")
        print("\n")
        print("Reinforcement bar:")
        print(f"Dm = {self.Dm} {length}")
        print(f"Dv = {self.Dv} {length}")
        print(f"rho = {round(self.rho*100,3)} %")
