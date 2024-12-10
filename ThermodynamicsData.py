import pandas as pd
import matplotlib.pyplot as plt
import prettytable as pt
import seaborn as sns
import numpy as np
from numpy import *
import math

class ThermodynamicsData():
    def __init__(self):        
        pd.set_option('display.unicode.east_asian_width',True)
        pd.set_option('display.width', 1000)  # 设置字符显示宽度
        pd.set_option('display.max_rows', 1000)  # 设置显示最大行
        pd.set_option('display.max_columns', 1000)  #最多显示五列

        self.R=8.314

        ##All T below is K


    def COgibbs(self,T):

        if T<=1700:
            self.CO_gibbs=-203680.536 + 609.354681*T + 3.117571299E-03*T**2\
                           +770627.667*T**(-1) + 32063.3157*math.log(T) - 10357.0197*T**0.5\
                           - 90.7535842*T*math.log(T)

        elif T<=6000:
            self.CO_gibbs=-224220.539 + 142.214896*T + 5816543.22*T**(-1)\
                           + 20893.9080*math.log(T) - 2749.10737*T**0.5 - 44.0863203*T*math.log(T)

        return self.CO_gibbs

    def CO2gibbs(self,T):

        if T<=1900:
            self.CO2_gibbs=-415578.806 + 642.658490*T + 2.371303142E-03*T**2\
                            +20124.5222*T**(-1) - 6993.14893*T**0.5 + 11004.7406*math.log(T)\
                            -103.344603*T*math.log(T)
        elif T<=6000:
            self.CO2_gibbs=402917.740 - 923.173242*T - 21647150.0*T**(-1)\
                            -0.824488778*T**1.5 + 30319.4765*T**0.5 - 196911.682*math.log(T)\
                            + 49.1972413*T*math.log(T)

        return self.CO2_gibbs

    def Cgibbs(self,T):

        if T<=1000:
            self.C_gibbs=20210.9773 + 515.671145*T + 4.008845366E-03*T**2\
                          -412784.987*T**(-1) - 4096.05580*T**0.5 + 1100795.94*T**(-2)\
                          -61.1790976*T*math.log(T) 


        elif T<=6000:
            self.C_gibbs=-13835.8469 + 176.395502*T - 4.793618847E-04*T**2\
                          + 1441704.16*T**(-1) - 158.769410*T**0.5 - 24.7780746*T*math.log(T)

        return self.C_gibbs

    def Fegibbs(self,T):  ##

        if T<=1184.81:  #S1
            self.Fe_gibbs=1225.70000 + 124.134000*T - 4.397520000E-03*T**2\
                           + 77358.5000*T**(-1) - 5.892690000E-08*T**3 - 23.5143000*T*math.log(T)\
                           +self.GmagFe(T)

        elif T<=1667.47: #S2
            self.Fe_gibbs=- 236.700000 + 132.416000*T - 3.757520000E-03*T**2\
                           + 77358.5000*T**(-1) - 5.892690000E-08*T**3 - 24.6643000*T*math.log(T)\
                           +self.GmagFe(T)

        elif T<=1810.95: #S1
            self.Fe_gibbs=1225.70000 + 124.134000*T - 4.397520000E-03*T**2\
                           + 77358.5000*T**(-1) - 5.892690000E-08*T**3 - 23.5143000*T*math.log(T)\
                           +self.GmagFe(T)

        elif T<=1812:  #S1-L1
            self.Fe_gibbs=- 24287.8308 + 298.768006*T - 46.0000000*T*math.log(T)\
                           +self.GmagFe(T)
            
            
        elif T<=3131.93: #L1
            self.Fe_gibbs=-10838.8389 + 291.301997*T - 46.0000000*T*math.log(T)

        return self.Fe_gibbs


    def FeOgibbs(self,T):

        if T<=1644.15:
            self.FeO_gibbs=- 322147.542 - 330.687435*T - 1.530402997E-02*T**2\
                            + 1266650.00*T**(-1) + 6003.60001*T**0.5 + 18.0244741*T*math.log(T)
            
        elif T<=2000:
            self.FeO_gibbs=- 268094.665 + 398.288735*T - 68.1992000*T*math.log(T)

        else:
            print('no data for FeO above 2000 K')

        '''
        elif T<=3900:
            self.FeO_gibbs=- 835030.307 + 1185.81676*T + 1.368558597E-03*T**2\
                            + 31173911.5*T**(-1) + 253414.297*math.log(T) - 37141.2387*T**0.5\
                            - 140.037464*T*math.log(T)
        '''
        return self.FeO_gibbs


    def Fe3O4gibbs(self,T):

        if T<=848.12:
            self.Fe3O4_gibbs=- 1200277.81 + 1282.95850*T - 1.901728207E-02*T**2\
                              + 3621668.52*T**(-1) + 3.967903931E-05*T**3 - 3.104596095E-08*T**4\
                              - 110726059*T**(-2) - 207.930830*T*math.log(T)

        elif T<=1870.00:
            self.Fe3O4_gibbs=- 1185253.60 + 1258.71640*T + 3621668.52*T**(-1)\
                              - 110726059*T**(-2) - 207.930830*T*math.log(T)     
            
        elif T<=2500:
            self.Fe3O4_gibbs=- 1053600.58 + 1230.41497*T - 213.384000*T*math.log(T)

        return self.Fe3O4_gibbs

    def Fe2O3gibbs(self,T): ##

        if T<=2500:
            self.Fe2O3_gibbs=- 861183.055 + 828.050052*T + 1453820.00*T**(-1)\
                              - 137.008930*T*math.log(T)\
                              +self.GmagFe2O3(T)

        return self.Fe2O3_gibbs

    def GmagFe(self,T):

        if T<=1184.81:  #S1
            Tc=1043.00
            beta=2.22
            P=0.4
            S=1

        elif T<=1667.47: #S2     
            Tc=67.00
            beta=0.7
            P=0.28
            S=0.3333

        elif T<=1810.95: #S1
            Tc=1043.00
            beta=2.22
            P=0.4
            S=1

        elif T<=1812:  #S1-L1
            Tc=1043.00
            beta=2.22
            P=0.4
            S=1

        tao=T/Tc
        
        D=518/1125+11692/15975*(P**(-1)-1)

        if tao<=1:
            g=1-(79*tao**(-1)/(140*P)+474/497*(P**(-1)-1)*(tao**(3)/6+tao**(9)/135+tao**(15)/600))/D
        
        elif tao>1:
            g=-(tao**(-5)/10+tao**(-15)/315+tao**(-25)/1500)/D
            
        self.G_mag=self.R*T*math.log(beta+1)*g

        return self.G_mag

    def GmagFe2O3(self,T):

        if T<=2500:  #S1
            Tc=955.67
            beta=8.3667
            P=0.28
            S=0.3333

        tao=T/Tc
        
        D=518/1125+11692/15975*(P**(-1)-1)

        if tao<=1:
            g=1-(79*tao**(-1)/(140*P)+474/497*(P**(-1)-1)*(tao**(3)/6+tao**(9)/135+tao**(15)/600))/D
        
        elif tao>1:
            g=-(tao**(-5)/10+tao**(-15)/315+tao**(-25)/1500)/D
            
        self.G_mag=self.R*T*math.log(beta+1)*g

        return self.G_mag


    def H2gibbs(self,T):

        if T<=1200:
            self.H2_gibbs=- 13779.8231 - 17.8247812*T - 1.538863732E-03*T**2\
                           + 147589.949*T**(-1) - 2.383074898E-07*T**3 + 779.444975*T**(0.5)\
                           - 19.8256305*T*math.log(T) 

        elif T<=4100:
            self.H2_gibbs=330318.973 - 121.445297*T - 1.032053794E-03*T**2\
                           - 10416577.5*T**(-1) - 74353.4040*math.log(T) + 8539.33826*T**(0.5)\
                           - 14.3874604*T*math.log(T) 
        return self.H2_gibbs


    def H2Ogibbs(self,T):

        if T<=373.13:
            self.H2O_gibbs=- 256638.942 - 1118.62474*T - 0.760349801*T**2\
                            - 1924378.83*T**(-1) + 5.318874400E-04*T**3 - 2.059132015E-07*T**4\
                            + 203.118982*T*math.log(T)

        elif T<=1100:
            self.H2O_gibbs=- 255475.808 - 15.2825865*T - 7.474858139E-03*T**2\
                            + 13999.6597*T**(-1) + 9.205931575E-08*T**3 + 1107.27182*math.log(T)\
                            - 25.7816397*T*math.log(T) 
            

        elif T<=4000:
            self.H2O_gibbs=152152.281 + 164.707573*T - 8.054003823E-05*T**2\
                            - 12075583.2*T**(-1) - 83128.2757*math.log(T) + 5947.37003*T**0.5\
                            - 53.1457895*T*math.log(T)

        return self.H2O_gibbs


    def O2gibbs(self,T):

        if T<=1000:
            self.O2_gibbs=-5219.33235 - 12.1798565*T - 8.489340766E-03*T**2\
                          - 114664.628*T**(-1) + 1.127694209E-06*T**3 - 316.646634*T**0.5\
                          - 26.9240574*T*math.log(T)
        elif T<=4000:
            self.O2_gibbs=- 389938.784  + 638.676784*T + 7.237224420E-04*T**2\
                          + 9341342.96*T**(-1) - 16506.1488*T**0.5 + 95803.9595*math.log(T)\
                          - 89.6813271*T*math.log(T)
        return self.O2_gibbs


if __name__ == '__main__':
    A=ThermodynamicsData()
    print(A.Fegibbs(800))
