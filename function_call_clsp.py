#CLSP-Rmin with back-ordering example
from gurobipy import *
import time
import math
import numpy as np


def compute_values(n,d,f,h,hR,tp,tpR,ts,tsR,b,k,M,s0,sR0,T,m,BO0,Pro,Pe,O,Rp,CostBo,realizationO):
    #Decision variables default settings are 0.0 lb for all continous variables (positive values)
    x = n.addVars(Pro, Pe, vtype=GRB.BINARY, name="x")
    s = n.addVars(Pro, Pe, vtype=GRB.CONTINUOUS, name="s")
    y = n.addVars(Pro, Pe, vtype=GRB.INTEGER, name="y")
    yS = n.addVars(Pro, Pe, vtype=GRB.INTEGER, name="yS")
    yR = n.addVars(Pro, Pe, vtype=GRB.INTEGER, name="yR")
    sR = n.addVars(Pro, Pe, vtype=GRB.CONTINUOUS, name="sR")
    xR = n.addVars(Pro, Pe, vtype=GRB.BINARY, name="xR")
    R = n.addVars(Pro, Pe, vtype=GRB.INTEGER, name="R")
    BO = n.addVars(Pro, Pe, vtype=GRB.CONTINUOUS, name="BO")
    #n.update()

    #Objective
    obj= quicksum(((f[j]*x[j,t]+h[j]*s[j,t])+(f[j]*xR[j,t]+hR[j]*sR[j,t])+((h[j]*m)*BO[j,t])) for j in Pro for t in Pe)
    n.setObjective(obj, GRB.MINIMIZE) 
    #n.update()

    #initialized back order variable NO NEED FOR CODING THIS CUS PYTHON INDEX 0 REFERS TO PERIOD 1
    #for j in Pro:
    #    if t==0:
    #        n.addConstr(BO[j,0]==BO0[j])
    #    else:
    #        n.addConstr(BO[j,t]==BO[j,t])
    #n.addConstrs((BO[j,0]==BO0[j]) for j in Pro) BO DOESNT NEED TO BE 0 AT THE END OF 1ST PERIOD

    #n.update()
    #s[j,T] the last period of product j equals 0
    n.addConstrs((s[j,T-1] == 0) for j in Pro)    #T is defined as T=len(b) to count how many periods we have =10                 
    n.addConstrs((sR[j,T-1] == 0) for j in Pro)
    if CostBo==True:
        n.addConstrs((BO[j,T-1] == 0) for j in Pro)
    else:
        pass
    #Inventory balance constraints
    for j in Pro:
        for t in Pe:
            if t==0:
                n.addConstr(s[j,t]-BO[j,t]==s0[j]+yS[j,t]+yR[j,t]-d[j][t]-BO0[j])
            else:
                n.addConstr(s[j,t]-BO[j,t]==s[j,t-1]+yS[j,t]+yR[j,t]-d[j][t]-BO[j,t-1])
    #n.update()
    for j in Pro:
        for t in Pe:
            if t==0:
                n.addConstr(sR[j,t]==sR0[j]+R[j,t]-yR[j,t])
            else:
                n.addConstr(sR[j,t]==sR[j,t-1]+R[j,t]-yR[j,t])


    #Capacity constraints
    n.addConstrs(
        (quicksum(((tp[j]*y[j,t]+ts[j]*x[j,t])+(tpR[j]*yR[j,t]+tsR[j]*xR[j,t])) for j in Pro) <= b[t]
            for t in Pe), "Capacity-Ct")

    #Big M+ Big M for rework #change from d[j][t:] to d[j][:]
    n.addConstrs((y[j,t]<= sum(d[j][:])*x[j,t] for j in Pro for t in Pe), "BigM-Ct")
    n.addConstrs((yR[j,t]<= sum(d[j][:])*xR[j,t] for j in Pro for t in Pe), "BigM-Ct for Rework")

    #!production amount at least as size of lot size
    for j in Pro:
        for t in Pe:
            n.addConstr(y[j,t]>=k[j]*x[j,t])
            #floor brackets  
            #m.addConstr(yS[j,t] ==math.floor((y[j,t]*(1-O[t])))) # this doesn't work with Gurobi, so use two constraints below
            n.addConstr(yS[j,t] >= (y[j,t]*(1-O[t]) - 0.9999))
            n.addConstr(yS[j,t] <= (y[j,t]*(1-O[t])))
            #ceil brackets
            #m.addConstr(R[j,t]==math.ceil((y[j,t]*O[t]))) # same here
            n.addConstr(R[j,t]>=(y[j,t]*O[t]))
            n.addConstr(R[j,t]<=(y[j,t]*O[t] + 0.9999))

    #!production amount at least as size of lot size
    for j in Pro:
        for t in Pe:
            n.addConstr(y[j,t]>=k[j]*x[j,t])
            #floor brackets  
            #m.addConstr(yS[j,t] ==math.floor((y[j,t]*(1-O[t])))) # this doesn't work with Gurobi, so use two constraints below
            n.addConstr(yS[j,t] >= (y[j,t]*(1-O[t]) - 0.9999))
            n.addConstr(yS[j,t] <= (y[j,t]*(1-O[t])))
            #ceil brackets
            #m.addConstr(R[j,t]==math.ceil((y[j,t]*O[t]))) # same here
            n.addConstr(R[j,t]>=(y[j,t]*O[t]))
            n.addConstr(R[j,t]<=(y[j,t]*O[t] + 0.9999))

            
    #Set time limit
    #n.Params.timeLimit =60.0
    #logging.basicConfig(filename='info.log', level=logging.INFO)

    n.optimize()

    #n.printAttr('x','y*')
    #n.printAttr('x','sR*')
    #n.update()

    #PRINT ALL RESULTS BELOW
       #holding cost
    pp=0
    for j in Pro:
        for t in Pe:
            pp= pp + (h[j]*s[j,t])
    #holding cost for rework
    ppp=0
    for j in Pro:
        for t in Pe:
            ppp= ppp + (hR[j]*sR[j,t])
    #setup cost combined
    pppp=0
    for j in Pro:
        for t in Pe:
            pppp=pppp+(f[j]*x[j,t]+f[j]*xR[j,t])
    #setup cost for production
    fcost=0
    for j in Pro:
        for t in Pe:
            fcost=fcost+(f[j]*x[j,t])
    #Nr of setups for production
    xprod=0
    for j in Pro:
        for t in Pe:
            xprod += x[j,t].x
    #Setup costs for rework
    fcostr=0
    for j in Pro:
        for t in Pe:
            fcostr=fcostr+(f[j]*xR[j,t])
    xxx=0
    #Capacity Utilization
    for j in Pro:
        for t in Pe:
            xxx=xxx+(((tp[j]*y[j,t]+ts[j]*x[j,t])+(tpR[j]*yR[j,t]+tsR[j]*xR[j,t]))/b[0])*0.1
            
    #Nr of setups for rework
    xrework=0
    for j in Pro:
        for t in Pe:
            xrework += xR[j,t].x
    #Back order costs
    cbcost=0
    for j in Pro:
        for t in Pe:
            cbcost=cbcost+(h[j]*m)*BO[j,t].x
            
 #realized calulations.................................................................................  

    yS = np.zeros((Pro[-1]+1, Pe[-1]+1))
    R = np.zeros((Pro[-1]+1, Pe[-1]+1))
    
    s = np.zeros((Pro[-1]+1, Pe[-1]+1))
    BO = np.zeros((Pro[-1]+1, Pe[-1]+1))
    
    sR = np.zeros((Pro[-1]+1, Pe[-1]+1))
    
    
    
    for j in Pro:
        for t in Pe:
            #floor brackets  
            yS[j,t] = math.floor((y[j,t].x*(1-realizationO[t])))  #math.floor works, no gurobi used here
            #ceil brackets
            R[j,t] = math.ceil((y[j,t].x*realizationO[t]))
            print('O', realizationO[t])
    
#     print(yS)
#     print(R)    
        
    
    for j in Pro:
        for t in Pe:
            if t==0:
                inv = s0[j]+yS[j,t]+yR[j,t]-d[j][t]-BO0[j]
            else:
                inv = s[j,t-1]+yS[j,t]+yR[j,t]-d[j][t]-BO[j,t-1]

 #inv calculated from gurobi values/expressions & no matrix created so i use getValue()

            if inv.getValue() > 0:
                s[j,t] = inv.getValue()
            else:
                BO[j,t] = -inv.getValue()
                
    for j in Pro:
        for t in Pe:
            if t==0:
                sR[j,t] = max(sR0[j] + R[j,t] - yR[j,t].x, 0)
            else:
                sR[j,t] = max(sR[j,t-1] + R[j,t] - yR[j,t].x, 0)
    
    #Back order realization costs 
    cbcost2=0
    for j in Pro:
        for t in Pe:
            cbcost2=cbcost2+(h[j]*m)*BO[j,t]
            
    cos=[] # empty BO0 by wyciagnac okrojone wartosci BO

    for j in Pro:
        cos.append(BO[j,Rp-1])
        
    mag=[] #s0
    for j in Pro:
        mag.append(s[j,Rp-1])
    
    magr=[] #sR0
    for j in Pro:
        magr.append(sR[j,Rp-1])
    
    
    #Objective fuction
    ZZ=n.ObjVal
    #Cost reduction
    #cred=0
    #for j in Pro:
    #    for t in Pe:
    #        cred=((ZZ*100)/(ZZ+cbcost))    
    


    #print(pp.getValue(), "Table 4 h Holding cost")
#     print(ppp.getValue(), "Table 4 hR Holding cost for rework")
#     print(pppp.getValue(), "Setup cost combined")
#     print(fcost.getValue(), "Table 4 f Setup cost for production")
#     print(xprod, "Table 4 x Nr of set-ups for production")
#     print(fcostr.getValue(), "Table 4 f Setup cost for rework")
#     print(xrework, "Table 4 xR Setup cost for rework")
#     print(xxx.getValue(), "Table 4 Capacity Utilization")
#     print(ZZ, "Table 4 Z objective function value")
#     print(cbcost,"Table 4 Back order costs cb")
#     print(cbcost2,"Table 4 Back order realization costs cb")
    #print(cred,"Table 4 Cost reduction")
    #print("Total time:", time.time()-start)
    
    return cos, mag, magr, pp.getValue(), ppp.getValue(), pppp.getValue(), fcost.getValue(), xprod, fcostr.getValue(), xrework, xxx.getValue(), ZZ, cbcost, cbcost2, yS, yR, s, sR, BO, x, xR

