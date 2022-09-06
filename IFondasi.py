from math import*
import pandas as pd
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px

a_factors = {
    'clay soft':3,
    'clay firm':5,
    'silt soft':7,
    'silt loose':12,
    'silt medium':15,
    'silt dense':20,
    'silty sand loose':20,
    'sand loose':22,
    'sand medium':28,
    'sand dense':35,
    'gravel loose':35,
    'gravel dense':45}
# table for calculate Nc*, for Qb on clay; Reese and O'Neil 1999
list_cu = [24,48,96,192]
list_Nc = [6.55,8.01,8.69,8.94]

xR = np.array([-2,-1.5,-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.5,2])
zR = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.5,2,2.5])
Iz = [0.5,0.48403,0.46782,0.45132,0.43452,0.41747,0.40028,0.36713,0.36589,0.34892,0.33223,0.30016,0.25624,0.196,0.15104]

# unit conversion
tsf = 107.25177
ft = 0.0328084 # cm
m = 0.3048 # ft

class Stresses:
    """The distribution of stresses within a soil from applied surface loads or stresses is determined
    by assuming that the soil is a semi-infinite, homogeneous, linear, isotropic, elastic material."""

    def __init__(self):
        pass

    def rect(z,B,L,q):
        """Asumsi pondasi dibagi empat"""
        Lx = B/2
        Ly = L/2
        m = Lx/z
        n = Ly/z

        A = 2*m*n
        B = m**2 + n**2 + 1
        C = m**2 + n**2 + 2
        D = (m**2) * (n**2)
        I = (1/(4*pi))*(((A*sqrt(B) / B+D)*(C/B)) + atan((A*sqrt(B)) / (B-D)))

        delta_p_z = 4*I*q
        return round(delta_p_z,3)

    def circular(z,D,q):
        r = D/2
        Rz = (1+((r/z)**2))**(3/2)
        I = 1-(1/Rz)
        pz = q*I
        return round(pz,2)

    def edge_circular(z,D,q):
        r = D/2
        zRi = z/r
        I = np.interp(zRi,zR,Iz)
        qz = q*I
        return round(qz,2)

class Profil:
    """class untuk menghandle masalah parameter untuk setiap lapisan tanah"""

    def __init__(self, id_="Profil 1"):
        self.id_ = id_

    def read_excel(self,name):
        self.data = pd.read_excel(name)

    def De(self,Cr,Cc,pc,po,Delta_p):
        p1 = po+Delta_p

        if po == pc:
            De = Cc*log10(p1/po)
        elif pc > po:
            if p1 < pc:
                De = Cr*log10(p1/po)
            else:
                De = (Cr*log10(pc/po)) + (Cc*log10(p1/pc))
        return De

    def rect_fondation(self,B,L,Df,q):
        self.q = q
        self.Df = Df
        self.B = B
        self.L = L

        df1 = pd.concat([self.data[:1]]*2, ignore_index=True)
        df2 = pd.concat([self.data[1:]], ignore_index=True)

        # Menggabungkan dua data frame
        df3 = pd.concat([df1,df2])
        df3.reset_index(inplace=True,drop=True)

        # modifikasi data frame: memasukkan 1 lapisan fondasi
        df3.at[0,'z_end (m)'] = Df
        df3.at[1,'z_start (m)'] = Df
        df3.insert(3,'h (m)',df3['z_end (m)'] - df3['z_start (m)'])
        df3.insert(2,'z_mid (m)',df3['z_start (m)'] + df3['h (m)']/2)

        #reset index
        df3.reset_index(inplace=True,drop=True)

        # menghitung tekanan overburden
        po_start = []
        po_mid = []
        po_end = []

        for i, gamma in enumerate(df3['g_ef (kN/m3)']):
            if i == 0:
                po_start.append(gamma*df3.loc[i,'z_start (m)'])
                po_end.append(po_start[i] + gamma*df3.loc[i,'h (m)'])
                po_mid.append(po_end[i] - gamma*df3.loc[i,'h (m)']/2)
            else:
                po_start.append(po_end[i-1])
                po_end.append(po_start[i] + gamma*df3.loc[i,'h (m)'])
                po_mid.append(po_end[i] - gamma*df3.loc[i,'h (m)']/2)

        df3['po_start'] = po_start
        df3['po_mid'] = po_mid
        df3['po_end'] = po_end
        self.qn = self.q - df3.loc[0,'po_end']
        self.data_modified = df3

    def primary_consolidation(self, method = 'classic'):
        delta_p = []
        delta_e = []
        Sc_z = []

        for i, z in enumerate(self.data_modified['z_mid (m)']):
            if z < self.Df:
                dp = None
                Delta_e = None
                Sc = None
            else:
                z = z-self.Df
                dp = Stresses.rect(z,self.B,self.L,self.qn)
                Cr = self.data_modified.loc[i,'Cr']
                Cc = self.data_modified.loc[i,'Cc']
                po = self.data_modified.loc[i,'po_mid']
                e0 = self.data_modified.loc[i,'e0']
                H = self.data_modified.loc[i,'h (m)']
                pc = self.data_modified.loc[i,'pc']

                if Cc == 0 and Cr == 0:
                    Sc = 0
                else:
                    if method == 'classic':
                        Delta_e = self.De(Cr,Cc,pc,po,dp)
                        Sc = Delta_e/(1+e0) * H
                    else:
                        Sc = self.janbu_cohesive(H,po,pc,dp,Cr,Cc,e0)

            delta_p.append(dp)
            delta_e.append(Delta_e)
            Sc_z.append(Sc)

        self.data_modified['Delta_p_z'] = delta_p
        self.data_modified['Delta_e'] = delta_e
        self.data_modified['Sc_z'] = Sc_z

        # Penurunan konsolidasi total
        self.Sc = round(self.data_modified['Sc_z'].sum(),5)

    def janbu_elastic(self):
        list_Se = []

        for i, dp in enumerate(self.data_modified['Delta_p_z']):
            if self.data_modified.loc[i,'z_mid (m)'] <= self.Df or dp == None:
                list_Se.append(0)
            else:
                j = 0.5 # untuk pasir
                H = self.data_modified.loc[i,'h (m)']
                po = self.data_modified.loc[i,'po_mid']
                D = self.data_modified.loc[i,'E (kN/m)']
                Se = self.janbu_sand(H,po,dp,j,D)
                list_Se.append(Se)

        self.data_modified['Se'] = list_Se
        self.Se = round(self.data_modified['Se'].sum(),5)

    def elastic_settlement(self):
        # average modulus elastic
        E_h = 0
        mu_h = 0
        z_bar = 0
        for i, E in enumerate(self.data_modified['E (kN/m)']):
            mu = self.data_modified.loc[i,'mu']
            h = self.data_modified.loc[i,'h (m)']
            E_h += E*h
            mu_h += mu*h
            z_bar += h

        self.E_avg = round(E_h/z_bar,2)
        self.mu_avg = round(mu_h/z_bar,2)

        # Budhu p450
        self.Ab = self.B*self.L
        self.Pn = self.Ab*self.qn
        AbL = self.Ab/(4*self.L**2)
        mu_s = 0.45*AbL**(-0.38)
        mu_emb = 1 - ((0.04*self.Df/self.B)*(1+(4/3 * AbL)))
        Eu = self.E_avg
        mu = self.mu_avg
        self.Se = round(self.Pn/(Eu*self.L) * (1-mu**2) *mu_s*mu_emb,5)

    def janbu_sand(self,h,po,dp,j,D):
        p1 = po+dp
        pr = 100 #kPa
        m = D/(j*dp) * ((p1/pr)**j - (po/pr)**j)
        strain = 1/m*j * ((p1/pr)**j - (po/pr)**j)
        Se = round(strain*h,5)
        return Se

    def janbu_cohesive(self,H,po,pc,dp,Cr,Cc,e0):
        p1 = po+dp
        m = log(10)*((1+e0)/Cc)
        mr = log(10)*((1+e0)/Cr)

        if po == pc:
            strain = (1/m)*log(p1/po)
        elif pc > po:
            if p1 < pc:
                strain = (1/mr)*log(p1/po)
            else:
                strain = ((1/mr)*log(pc/po)) + ((1/m)*log(p1/pc))
        Sc = strain*H
        return round(Sc,5)

    def m_qc(self,po,qt,Ko,soiltype):
        # empirical correlation Massarsch 1994
        a = a_factors[soiltype]
        pr = 100
        pv = po
        Qtn = qt/sqrt(pr*pv)
        K = sqrt(3/(1+(2*Ko)))
        m = a*sqrt(Qtn*K)

class Shallow:
    """Bearing capacity shallow foundation"""

    def __init__(self,name="P1",B=1,L=0,Df=0,shape="rectangular"):
        self.name = name
        self.B = B
        self.shape = shape
        self.Df = Df

        if shape == "circular":
            self.D = B
            self.L = B
        else:
            if L == 0:
                self.L = B
            else:
                self.L = L

        self.stress_influence_df = pd.read_excel("ahlvin_ulery.xlsx")

    def checkLocalShear(self):
        if self.phi < 36:
            self.phi = round(np.arctan((2/3)*np.tan(np.radians(self.phi)))*180/np.pi,3)

    def soil_ESA(self,GWL=99,phi=30,g=18,phio=0,go=0, CPT=None, method="Davis and Booker",local_shear=False):
        self.analysis_type = "ESA"
        self.local_shear = local_shear
        if CPT is None:
            self.GWL = GWL
            self.g = g
            self.phi = phi
            if go == 0:
                self.go = g
            else:
                self.go = go

            if phio == 0:
                self.phio = phi
            else:
                self.phio = phio

            if self.local_shear == True:
                self.checkLocalShear()
        else:
            self.df = CPT.df
            self.GWL = CPT.GWL
            i_start = int(self.Df/0.2)
            i_end = int((self.B+self.Df)/0.2)
            self.phi = round(np.mean(self.df.loc[i_start:i_end,'phi']),3)
            self.phio = round(np.mean(self.df.loc[:i_start,'phi']),3)
            self.g = round((np.mean(self.df.loc[i_start:i_end,'g'])),2)
            self.go = round((np.mean(self.df.loc[:i_start,'g'])),2)
            self.Su = round((np.mean(self.df.loc[i_start:i_end,'Su'])),2)

            if self.local_shear == True:
                self.checkLocalShear()

        if local_shear == True:
            if self.phi <= 36.0:
                phi_local = round(np.arctan((2/3)*np.tan(np.radians(self.phi)))*180/np.pi,2)
                self.phi = phi_local

        self.bearingCapacityFactors(method=method)

    def soil_TSA(self,GWL=99,Su=10,g=19,go=18,CPT=None):
        self.analysis_type = "TSA"

        if CPT is None:
            self.Su = Su
            self.g = g
            self.go = go
            self.GWL = GWL
            if go == 0:
                self.go = g
            else:
                self.go = go
        else:
            self.df = CPT.df
            self.GWL = CPT.GWL
            i_start = int(self.Df/0.2)
            i_end = int((self.B+self.Df)/0.2)
            self.g = round((np.mean(self.df.loc[i_start:i_end,'g'])),2)
            self.go = round((np.mean(self.df.loc[:i_start,'g'])),2)
            self.Su = round((np.mean(self.df.loc[i_start:i_end,'Su'])),2)

    def bearingCapacityFactors(self,method="Davis and Booker"):
        phi = np.radians(self.phi)
        self.Nq = round(np.e**(np.pi*np.tan(phi)) * np.tan(np.radians(45) + phi/2)**2,2)

        if method == "Davis and Booker":
            # Ng = 0.0663*np.exp(9.3*self.phi) # for smooth footing
            self.Ng = round(0.1054*np.exp(9.6*np.radians(self.phi)),2) # for rough footing
        elif method == "Vesic":
            self.Ng = round(2*(self.Nq+1)*np.tan(phi),2)
        elif method == "Meyerhof":
            self.Ng = round((self.Nq-1)*np.tan(1.4*phi),2)
        else:
            raise("method error")

    def geometricFactorsTSA(self,H):
        sc = 1 + 0.2*(self.Be/self.Le)
        dc = 1 + 0.33*np.arctan(self.Df/self.Be)

        nb = (2 + self.Be/self.Le) / (1 + self.Be/self.Le)
        ic = 1 - (nb*H)/(5.14*self.Su*self.Be*self.Le)

        eta = 0 # base slope angle. beta < phi; eta + beta < 90 deg
        beta = 0 # slope angle. beta < phi; eta + beta < 90 deg
        bc = 1 - (eta/147)
        gc = 1 - (beta/147)

        factor = [sc,dc,ic,bc,gc]

        return factor

    def geometricFactorsESA(self,Vn,H):
        phi = np.radians(self.phi)
        if H == self.Hb:
            n = (2 + self.Be/self.Le) / (1 + self.Be/self.Le)
        elif H == self.Hl:
            n = (2 + self.Le/self.Be) / (1 + self.Le/self.Be)

        eta = 0 # base slope angle. beta < phi; eta + beta < 90 deg
        beta = 0 # slope angle. beta < phi; eta + beta < 90 deg

        sq = round(1 + (self.Be/self.Le)*np.tan(phi),2)
        dq = round(1 + (2*np.tan(phi) * (1-np.sin(phi))**2 * (self.Df/self.Be)),2)
        iq = round((1 - H/Vn)**n,2)
        bq = round((1 - eta*np.tan(phi))**2,2)
        gq = round((1 - np.tan(beta))**2,2)

        sg = round(1 - (0.4*(self.Be/self.Le)),2)
        dg = 1
        ig = round((1 - H/Vn)**(n+1),2)
        bg = bq
        gg = gq

        factor_q = [sq,dq,iq,bq,gq]
        factor_g = [sg,dg,ig,bg,gg]

        return factor_q, factor_g

    def circular(self,Vn,H,M):
        self.Vn = Vn
        self.H = H
        self.M = M
        self.e = M/Vn

        if self.e > self.D/6:
            raise("Eksentrisitas besar! revisi dimensi")

        self.Be = self.D - (2*self.e)
        self.Le = self.Be
        self.A = round(0.25*np.pi*self.D**2,3)
        self.Ae = round((self.D**2 / 2) * (np.arccos((2*self.e)/self.D) - ((self.e/self.D) * (1-(2*self.e/self.D)**2)**0.5)),3)

    def rectangular(self,Vn,Hb,Hl,Mb,Ml):
        self.eb = Mb/Vn
        self.el = Ml/Vn

        if self.eb > self.B/6:
            raise("e > B/6; revisi lebar fondasi (B)")
        elif self.el > self.L/6:
            raise("e > L/6; revisi panjang fondasi (L)")

        self.Vn = Vn
        self.Hb = Hb
        self.Hl = Hl
        self.Mb = Mb
        self.Ml = Ml

        self.Be = self.B - (2*self.eb)
        self.Le = self.L - (2*self.el)
        self.Ae = self.Be*self.Le
        self.A = self.B*self.L

        self.maxStressX = round((self.Vn/self.A) * (1 + (6*self.eb / self.B)),3)
        self.minStressX = round((self.Vn/self.A) * (1 - (6*self.eb / self.B)),3)

        self.maxStressY = round((self.Vn/self.A) * (1 + (6*self.el / self.L)),3)
        self.minStressY = round((self.Vn/self.A) * (1 - (6*self.el / self.L)),3)

    def groundWaterEffect(self):
        if self.B+self.Df < self.GWL:
            self.qo = round(self.go*self.Df,3)
        elif self.GWL == 0:
            self.g = self.g-9.81
            self.go = self.go-9.81
            self.qo = round(self.go*self.Df,3)
        elif self.Df > self.GWL:
            dw = self.Df-self.GWL
            self.g = self.g-9.81
            self.qo = ((self.go-9.81)*(self.Df-dw)) + (self.go*dw)
        else:
            z = self.GWL - self.Df
            if z < self.B:
                self.g = round((self.g-9.81) + ((z/self.B)*(self.g-(self.g-9.81))),3)
                self.qo = round(self.go*self.Df,3)
            else:
                raise("Error GWL")

    def loads(self,Vn=1,H=0,M=0,Hb=0,Hl=0,Mb=0,Ml=0):
        self.Vn = Vn
        self.H = H
        if self.shape == "circular":
            self.H = H
            self.Hb = Hb
            self.M = M
            self.circular(Vn,Hb,M)
            self.geometricFactorsESA(Vn,Hb)
        else:
            self.Hb = Hb
            self.Hl = Hl
            self.Mb = Mb
            self.Ml = Ml
            self.rectangular(Vn,Hb,Hl,Mb,Ml)

    def ESA(self,Be,factor):
        factor_q, factor_g = factor
        sq,dq,iq,bq,gq = factor_q
        sg,dg,ig,bg,gg = factor_g
        qu = (self.qo*(self.Nq-1)*sq*dq*iq*bq*gq) + (0.5*self.g*Be*self.Ng*sg*dg*ig*bg*gg)
        return qu

    def TSA(self,Be,factor):
        sc,dc,ic,bc,gc = factor
        qu = 5.14*self.Su*sc*dc*ic*bc*gc # p.432 Budhu
        return qu

    def checkSuitability(self):
        if self.Pa > self.Vn:
            self.suitability = f'Vn = {self.Vn} kN < Pa = {self.Pa} kN; (ecceptable)'
        else:
            self.suitability = f'Vn = {self.Vn} kN > Pa = {self.Pa} kN; (Unecceptable)'

    def analyze(self,FS=3):
        self.groundWaterEffect()

        if self.analysis_type == "ESA":
            if self.shape == "rectangular":
                self.factorx = self.geometricFactorsESA(self.Vn,self.Hl)
                self.qu_x = round(self.ESA(self.Be,self.factorx),0)
                self.qult_x = round(self.qu_x+self.qo,0)
                self.qa_x = round(self.qu_x/FS + self.qo,0)
                self.Pa_x = round(self.qa_x*self.Ae,0)

                self.factory = self.geometricFactorsESA(self.Vn,self.Hb)
                self.qu_y = round(self.ESA(self.Le,self.factory),0)
                self.qult_y = round(self.qu_y+self.qo,0)
                self.qa_y = round(self.qu_y/FS + self.qo,0)
                self.Pa_y = round(self.qa_y*self.Ae,0)

                self.Pa = min(self.Pa_x,self.Pa_y)
                self.qa = min(self.qa_x,self.qa_y)
                self.qu = min(self.qu_x,self.qu_y)
                self.qult = min(self.qult_x,self.qult_y)
                self.checkSuitability()
            else:
                self.factor = self.geometricFactorsESA(self.Vn,self.Hb)
                self.qu = round(self.ESA(self.Be,self.factor),0)
                self.qult = round(self.qu+self.qo,0)
                self.qa = round(self.qu/FS + self.qo,0)
                self.Pa = round(self.qa*self.Ae,0)
                self.checkSuitability()

        elif self.analysis_type == "TSA":
            self.factor = self.geometricFactorsTSA(self.H)
            self.qu = round(self.TSA(self.Be,self.factor),3)
            self.qult = round(self.qu+self.qo,3)
            self.qa = round(self.qu/FS + self.qo,3)
            self.Pa = round(self.qa*self.Ae,3)

    def settlement_DeBeer(self,h,qc,svo_ef,p1):
        if svo_ef == 0:
            Si = 0
        else:
            C = (1.5*qc)/svo_ef
            Si = (h/C)*np.log10(p1/svo_ef)
        return round(Si*1000,2)

    def settlement_Schmertmann(self,z,svo_ef,qt,dz,q_net):
        cd = max(1 - (0.5*(svo_ef/q_net)),0.5) # depth factor
        time = 10 # year
        ct = 1 + 0.2*np.log10(time/0.1) # creep factor
        I_cp = 0.5 + (0.1*np.sqrt(q_net/svo_ef))
        if z/self.B < 0.5:
            I_co = 0.1 + (2*(I_cp-0.1)*(z/self.B))
        elif 2 > z/self.B:
            I_co = I_cp*(1-(2/3 * ((z/self.B)-0.5)))
        else:
            I_co = 0
        beta = 2.5 # axisymetric
        S = (cd*ct/beta) * q_net * (I_co/qt) * dz

        return round(S*1000,2)

    def settlement_Janbu(self,p1,dz,m,j,svo_ef):
        pr = 100 # kPa
        M = m*j
        P0 = (svo_ef/pr)**j
        P1 = (p1/pr)**j
        strain_e = (1/M)*(P1-P0)
        Sz = strain_e*dz
        return round(Sz*1000,2)

    def settlement(self,method="Janbu",q=0):
        self.q_net = (self.Vn/self.Ae)
        if q == 0:
            q_net = self.q_net
        else:
            q_net = q

        self.settl_method = method
        list_settl_mid = []
        list_settl_edge = []

        for i, z in enumerate(self.df['z']):
            if i == 0:
                Sz = 0
                Sz_edge = 0
                S_db = 0
                S_sm = 0
            else:
                if z > self.B or z < self.Df:
                    Sz = 0
                    Sz_edge = 0
                    S_db = 0
                    S_sm = 0
                else:
                    if self.shape == 'circular':
                        pz = Stresses.circular(z,self.B,q_net)
                        pz_edge = Stresses.edge_circular(z,self.B,q_net)
                    else:
                        pz = Stresses.rect(z,self.B,self.L,q_net)
                        pz_edge = pz
                    qt = self.df.loc[i,'qt']
                    svo_ef = self.df.loc[i,'svo_ef']
                    m = self.df.loc[i,'m']
                    j = self.df.loc[i,'j']
                    dz = 0.2 # m; only for sondir
                    p1 = pz + svo_ef
                    p1_edge = pz_edge + svo_ef

                    if self.settl_method == "Janbu":
                        Sz = self.settlement_Janbu(p1,dz,m,j,svo_ef)
                        Sz_edge = self.settlement_Janbu(p1_edge,dz,m,j,svo_ef)
                    elif self.settl_method == "Schmertmann":
                        Sz = self.settlement_Schmertmann(z,svo_ef,qt,dz,q_net)
                        Sz_edge = 0
                    elif self.settl_method == "De Beer":
                        Sz = self.settlement_DeBeer(dz,qt,svo_ef,p1)
                        Sz_edge = 0

            list_settl_mid.append(Sz)
            list_settl_edge.append(Sz_edge)
        return list_settl_mid, list_settl_edge

    def solve_settl(self,method="Janbu"):
        list_settl_mid, list_settl_edge = self.settlement(method=method)

        self.df['Sz'] = list_settl_mid
        self.df['Sz edge'] = list_settl_edge

        cumulative_Sz = []
        cumulative_Sz_edge = []
        for i, z in enumerate(self.df['z']):
            if z < self.Df or z > self.B:
                Sz_cum = 0
                Sz_cum_edge = 0
            else:
                Sz_cum = sum(self.df.loc[i:,'Sz'])
                Sz_cum_edge = sum(self.df.loc[i:,'Sz edge'])
            cumulative_Sz.append(Sz_cum)
            cumulative_Sz_edge.append(Sz_cum_edge)
        self.df['Sz_tot'] = cumulative_Sz
        self.df['Sz_tot edge'] = cumulative_Sz_edge

        self.S = round(sum(self.df['Sz']),0)
        self.S_edge = round(sum(self.df['Sz edge']),0)

    def qa_settlement(self,method="Janbu",settl_target=25):
        q = 1
        while True:
            list_settl_mid, list_settl_edge = self.settlement(method=method,q=q)
            S = sum(list_settl_mid)
            if abs(S-settl_target) < 1:
                break
            else:
                q += (q*(abs(settl_target-S)/settl_target))/2
        return round(q,2)

    def qa_varwidth(self,width=[1,18],interval=0.5,settl_target=25,FS=2.5):
        self.settl_target = settl_target
        list_B = np.arange(width[0],width[1],interval)

        list_qa_shear = []
        list_qa_settl = []

        Df = self.Df
        q = 1

        for B in list_B:
            if self.analysis_type == "ESA":
                self.B = B
                self.Be = B
                qa_settl = self.qa_settlement(settl_target=settl_target)
                list_qa_settl.append(qa_settl)

                qu = round(self.ESA(B,self.factor),3)
                qa = round(qu/FS + self.qo,3)
                list_qa_shear.append(qa)
            elif self.analysis_type == "TSA":
                pass

        self.qa_var = pd.DataFrame()
        self.qa_var['B (m)'] = list_B
        self.qa_var['qa shear (kPa)'] = list_qa_shear
        self.qa_var['qa settl (kPa)'] = list_qa_settl

    def plot_qa_var(self):
        f, ax = plt.subplots()
        var = self.qa_var
        ax.plot(var["B (m)"],var["qa shear (kPa)"], label = r"$q_a$ shear")
        ax.plot(var["B (m)"],var["qa settl (kPa)"], label = f"$q_a$ settl. {self.settl_target} mm")
        ax.legend()
        ax.set_xlabel(r'Width, B (m)')
        ax.set_ylabel(r'$q_a$ (kPa)')
        # ax.grid(True)

    def settl_contour(self):
        list_x = xR*self.B/2
        list_y = zR*self.B/2

        z_new = np.arange(0,max(self.df['z'])+0.2,0.2)
        df_Influence = pd.DataFrame(np.zeros([len(z_new),len(list_x)]))
        #------------------------------------------------------
        for i in range(17):
            zp = list_y
            Ip = self.stress_influence_df.iloc[:,i]
            Ip_new = np.interp(z_new, zp, Ip)
            df_Influence[i] = Ip_new
        #------------------------------------------------------
        q_n = self.q_net
        df_pz = df_Influence * q_n
        df_settl = df_pz
        for i in range(len(z_new)):
            list_Sz = []
            svo_ef = self.df.loc[i,'svo_ef']
            m = self.df.loc[i,'m']
            j = self.df.loc[i,'j']
            qt = self.df.loc[i,'qt']
            z = self.df.loc[i,'z']
            dz = 0.2

            for xi in range(len(list_x)):
                if abs(list_x[xi]) < self.B/2 and z_new[i] < self.Df:
                    Sz = None
                else:
                    pz = df_pz.iloc[i,xi]
                    p1 = pz + svo_ef
                    if self.settl_method == "Janbu":
                        Sz = self.settlement_Janbu(p1,dz,m,j,svo_ef)
                    elif self.settl_method == "De Beer":
                        Sz = self.settlement_DeBeer(dz,qt,svo_ef,p1)
                list_Sz.append(Sz)
            df_settl.loc[i] = list_Sz

        settl_cumulative = df_settl
        for i in range(df_settl.shape[0]):
            i = df_settl.shape[0]- 1 - i
            if i == df_settl.shape[0]-1:
                pass
            else:
                settl_cumulative.loc[i] =settl_cumulative.loc[i] + settl_cumulative.loc[i+1]

        colorscale = [[0, 'lightsalmon'], [0.5, 'mediumturquoise'], [1, 'gold']]
        fig = go.Figure(data = go.Contour(z=settl_cumulative.values.tolist(),x=list_x,y=z_new, colorscale=colorscale,contours=dict(
                            showlabels = True, # show labels on contours
                            labelfont = dict(size = 12,color = 'white'))))
        fig.update_yaxes(autorange="reversed")
        fig.update_layout(title="Settlement (mm)",xaxis_title = "Distance (m)", yaxis_title = "Depth (m)",autosize=False,width=800,height=500)
        fig.show()

    def info(self):
        print(f"qu = {self.qu} kPa")
        print(f"qult = {self.qult} kPa")
        print(f"qa = {self.qa} kPa")
        print(f"Pa = {self.Pa} kN")
