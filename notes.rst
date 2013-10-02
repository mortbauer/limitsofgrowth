Notes
#####
Frage 2
*******
Frage 2
* in welchem Zeitbereich soll sich das Gleichgewicht einstellen?
* Sensitivitätsanalyse, was verstehen sie hier genau?

Sensitivitätsanalyse
********************
x1
    Bevölkerung
x2
    Umweltbelastung
x3
    Wirtschaft
p1 
    Geburtenrate
p2 
    Sterberate
p3
    Regenerationsrate
p4
    Belastungsrate
p5
    Wirtschaftsziel
p6 
    Zuwachsrate


Ableitungen von dx, bezeichnet als f, nach x:
=============================================

f1,1 = Geburtenrate * Umweltqualität * Wirtschaft - Sterberate * Umweltbelastung
f1,1 = p1 * 1/x2 * x3 - p2 * x2
f1,2 = Geburtenrate * Bevölkerung * (-1) * Umweltqualität**2 * Wirtschaft - Bevölkerung * Sterberate
f1,2 = p1 * x1 * (-1) * x2**2 * x3 - x1 * p2
f1,3 = Geburtenrate * Bevölkerung * Umweltqualität 
f1,3 = p1 * x1 * 1/x2

f2,1 = Belastungsrate * Wirtschaft 
f2,1 = p4 * x3
f2,2 = Regenerationsrate if Umweltqualität > 1 else 0
f2,2 = p3 if 1/x2 > 1 else 0
f2,3 = Belastungsrate * Bevölkerung 
f2,3 = p4 * x1

f3,1 = 0
f3,2 = Zuwachsrate * Wirtschaft - 2 * Zuwachsrate * Wirtschaft**2 * Umweltbelastung / Wirtschaftsziel
f3,2 = p6 * x3 - 2 * p6 * x3**2 * x2 / p5
f3,3 = Zuwachsrate * Umweltbelastung  - 2 * Zuwachsrate * Umweltbelastung**2  * Wirtschaft / Wirtschaftsziel
f3,3 = p6 * x2 - 2 * p6 * x2**2 * x3 / p5

Ableitungen von dx, bezeichnet als f, nach p:
=============================================

f1,p1 = Bevölkerung * Umweltqualität * Wirtschaft
f1,p1 = x1 * 1/x2 * x3
f1,p2 = - Bevölkerung * Umweltbelastung
f1,p2 = - x1 * x2
f1,p3 = 0
f1,p4 = 0
f1,p5 = 0
f1,p6 = 0

f2,p1 = 0
f2,p2 = 0
f2,p3 = - Umweltbelastung if Umweltqualität > 1 else -1
f2,p3 = - x2 if 1/x2 > 1 else -1
f2,p4 = Wirtschaft * Bevölkerung 
f2,p4 = x3 * x1
f2,p5 = 0
f2,p6 = 0

f3,p1 = 0
f3,p2 = 0
f3,p3 = 0
f3,p4 = 0
f3,p5 = Zuwachsrate * Wirtschaft**2 * Umweltbelastung**2 / Wirtschaftsziel**2
f3,p5 = p6 * x3**2 * x2**2 / p5**2
f3,p6 = Wirtschaft * Umweltbelastung * (1-(Wirtschaft * Umweltbelastung)/Wirtschaftsziel)
f3,p6 = x3 * x2 * (1-(x3*x2/p5))


Ableitungen von x0 nach p:
==========================
x,p = 0 (da immer x0 = 1)
