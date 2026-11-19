from ngsolve import *
from ngsolve.webgui import Draw
from netgen.occ import *
from netgen.geom2d import SplineGeometry
import numpy as np
from math import sqrt
import inspect
"""
der Code simuliert die Kirchhoff-love Gl mit HHJ-methode, 
1. definiert Auflagerpunkte
2. schneidet diese als kleine Kreise aus Geometrie aus und setzt diese als Dirichlet gleich 0
3. simuliert
4. speichert und retourniert xyz-File der GravitationsWerte
gibt außerdem noch maximale Durchbiegung aus (größter Wert - kleinster Wert)

"""
       
def gravitationsEinfluss_Auflager(l,b,t,step,maxh,order):
    #region GEOMETRIE
    shape = MoveTo(0,0).Rectangle(l,b).Face()
    shape.edges.name="dirichlet"
    geo = OCCGeometry(shape,dim=2)

    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    mesh.Curve(3)
    print("geometry - done")
    #Draw(mesh)
    #endregion

    #region MATERIALKONSTANTEN
    E  = 70e6      #Glas ~ N'/mm² Elastizitätsmodul - ACHTUNG: N' nicht gleich N!!!      
    """
    UMRECHNUNG (m -> mm)

    Einheit Pascal: 1 Pa = 1 N/m² 
    Einheit Newton: 1 N = 1 kg*m/s²     - dh müssen auch Newton umrechnen

    Definiere:
    N' := kg*mm/s²      - in Millimeter!!!!
        -> 1 N = 1000 N'

    1 Pa = 1 N/m² = 10e-6 N/mm² = 10e-6 * 10e3 N'/mm² = 10e-3 N'/mm²

    Glas hat E-Modul von 70e9 Pa (70 GPa) = 70e9 N/m² = 70e9 * 10^(-3) N'/mm² = 70e6 N'/mm²         #neuer Wert
    """
    nu = 0.23       #dimensionslos, also bei m und mm gleich
    rho = 2.5e-6   #kg/mm³ Dichte 
    """
    Dichte von Glas 2500 kg/m³ = 2.5 * 10e3 kg/m³
    1kg/m³ = 1000g/10e9mm³ = 10e(-6)g/mm³
    -> 2.5 * 10e3 * 10e(-6) g/mm³ = 2.5 * 10e(-3) g/mm³ = 0.0025 g/mm³

    """

    g = 9810     # Erdbeschleunigung mm/s²
    """
    Erdbeschleunigung g = 9.81 m/s² = 9810 mm/s²
    """   

    q = rho * t * g     #Eigengewicht der Platte punktweise!!! - rechte Seite der PDE, also unser f_PDE
    print("materialkonstanten - done")
    #endregion

    #region FUNKTIONEN für PDE
    Db = E*t**3/(12*(1-nu**2))

    def D(A):
        return Db *((1-nu)*A+ nu*Trace(A)*Id(2))

    def Dinv(A):
        return 1/Db * (1/(1-nu)*A - nu/(1-nu**2)*Trace(A)*Id(2))
    print("funktionen deklarieren - done")
    #endregion

    #region SIMULATION: FUNKTIONENRAUM,SCHWACHE FORM,VISUALISIERUNG
    Q = H1(mesh, order=order+1, dirichlet="dirichlet")
    V = HDivDiv(mesh, order=order,dirichlet="dirichlet")
    X = V * Q

    (sigma,w),(tau,v)= X.TnT()

    n = specialcf.normal(2)

    def tang(u):
        return u - (u*n)*n

    a = BilinearForm(X,symmetric=True)
    a += InnerProduct(Dinv(sigma),tau) * dx
    a += div(sigma)*Grad(v) * dx
    a += div(tau) * Grad(w) * dx
    a += - (sigma[n,:] * tang(Grad(v)) + tau[n,:] * tang(Grad(w))  ) * dx(element_boundary=True)
    a.Assemble()

    L = LinearForm(X)
    L += q * v *dx
    L.Assemble()
    print("a und L assemble - done")

    gf_solution = GridFunction(X)
    gf_solution.vec.data = a.mat.Inverse(X.FreeDofs(),inverse="") * L.vec

    gf_sigma, gf_w = gf_solution.components
    #Draw(gf_sigma, mesh,name="sigma")
    #Draw(100*gf_w, mesh, name="disp",deformation=True, euler_angles=[-60,5,30])
    print("simulation - done")
    #endregion
    
    #region xyzFiles
    #erzeugt gitter für gewünschte werte
    x_points_ungerundet = np.linspace(0,l, step)
    x_points = np.round(x_points_ungerundet, 2)
    y_points_ungerundet = np.linspace(0,b, step)
    y_points = np.round(y_points_ungerundet, 2)

    X_mesh, Y_mesh = np.meshgrid(x_points, y_points)

    #U
    #speichert (x,y,z) in gravitation_matrix für alle (x,y) in X_mesh,Y_mesh
    u_matrix = np.zeros((step, step),dtype=object)
    for i in range(step):
        for j in range(step):
            try:
                zEintrag_u = gf_w(X_mesh[i, j], Y_mesh[i, j])
            except:
                zEintrag_u = 0
            u_matrix[i, j] = (float(X_mesh[i,j]), float(Y_mesh[i,j]), zEintrag_u)

    #*.xyz file erstellen
    filename_u = f"xyzFiles/h{maxh}_s{step}_o{order}_u.xyz"
    with open(filename_u, "w") as f:
        for i in range(step):
            for j in range(step):
                x, y, w = u_matrix[i, j]
                f.write(f"{x:.10f}\t{y:.10f}\t{w:.10f}\n")


    #SIGMA
    sigma_matrix = np.zeros((step, step),dtype=object)
    for i in range(step):
        for j in range(step):
            try:
                zEintrag_sigma = gf_sigma(X_mesh[i, j], Y_mesh[i, j])
            except:
                zEintrag_sigma = 0
            sigma_matrix[i, j] = (float(X_mesh[i,j]), float(Y_mesh[i,j]), zEintrag_sigma[0],zEintrag_sigma[1],zEintrag_sigma[3])
    # print(sigma_matrix[30,30])
    # print(sigma_matrix[round(step/2),round(step/2)])

    #*.xyz file erstellen
    filename_sigma = f"xyzFiles/h{maxh}_s{step}_o{order}_sigma.xyz"
    with open(filename_sigma, "w") as f:
        for i in range(step):
            for j in range(step):
                x, y, sigma_a,sigma_b,sigma_c = sigma_matrix[i, j]

                f.write(
                    f"{x:.10f}\t{y:.10f}\t"
                    f"{sigma_a:.10f}\t{sigma_b:.10f}\t{sigma_c:.10f}\n"
                )
                # x, y, w = sigma_matrix[i, j]
                # f.write(f"{x:.10f}\t{y:.10f}\t{w:.10f}\n")
                
    print("xyzFiles - done")
    #endregion

    # returnWert ist egl nicht wichtig, weil wir uns den dateinamen eh aus den anderen variablen bestimmen können
    return filename_u, filename_sigma

#hier egl unnoetig
def calc_L2norm_abstand(filename_ref,maxh, step, order, mass_Omega, function):
    """
    wir berechnen die L2norm über ein funktionenraster, dabei erhalten wir folgende formel
    ||uh-uref||_2 = sqrt( sum_{i=1}^{anzahl_punkte} (u_h(i)-u_ref(i))^2 * maß(omega)/anzahl_pkt)
    = sqrt( summe[ (uh(i)-uref(i))^2 ] * maß_Omega / anzahl_pkte )
    """
    filename_func = f"xyzFiles/h{maxh}_s{step}_o{order}_{function}.xyz"
    data_func = np.loadtxt(filename_func)
    #x_uh = data_func[:,0]    
    #y_uh = data_func[:,1]
    z_uh = data_func[:,2]

    data_u_ref = np.loadtxt(filename_ref)
    #x_u_ref = data_u_ref[:,0]    
    #y_u_ref = data_u_ref[:,1]
    z_u_ref = data_u_ref[:,2]

    normierungsZahl = step*step
    sum = 0
    maxDiff = 0

    for i in range(normierungsZahl):
        diff = abs(z_uh[i]-z_u_ref[i])
        #berechnen Maximalabstand
        if diff > maxDiff:
            maxDiff = diff

        #quadrierter summand für integral-summe
        sum += diff*diff
    integral = sum/normierungsZahl*mass_Omega
    L2_norm = sqrt(integral)
    # print("L2norm - done")
    return L2_norm, maxDiff

if __name__ == "__main__":
    l = 200
    b = 100
    t = 5
    step = 40
    maxh = 100 # in {10, 50, 100, 200}
    order = 4
    mass_Omega = l*b
    print("gestartet!")
    # maxh_array = [1,5,10,50,100,200] 
    maxh_array = [500,200,100,50,10]
    # maxh_array = [500]
    #zuerst ersten teil ausführen, dann den zweiten code für alle fälle ausführen
    #entweder: referenzdaten berechnen

    for maxh_ref in maxh_array:
        print(maxh_ref)
        filename_u_ref, filename_sigma_ref = gravitationsEinfluss_Auflager(l,b,t,step,maxh_ref,order)
        print(f"calculation of maxh={maxh_ref} - done")
    #oder 
    #andere funktion berechnen und fehlerauswertung
    # for maxh_h in maxh_array:
    #     print("maxh = ",maxh_h)
    #     # filename_u_ref = f"xyzFiles/h{maxh_ref}_s{step}_o{order}_u.xyz"
    #     # filename_sigma_ref = f"xyzFiles/h{maxh_ref}:_s{step}_o{order}_sigma.xyz"
    #     filename_u_ref = f"xyzFiles/h1_s40_o3_u.xyz"
    #     filename_sigma_ref = f"xyzFiles/h1_s40_o3_sigma.xyz"
    #     #filename_u,filename_sigma = gravitationsEinfluss_Auflager(l,b,t,step,maxh,order)

    #     L2fehler_u, maxDiff_u = calc_L2norm_abstand(filename_u_ref, maxh_h, step, order, mass_Omega, "u")
    #     L2fehler_sigma, maxDiff_sigma = calc_L2norm_abstand(filename_sigma_ref, maxh_h, step, order, mass_Omega, "sigma")

    #     print(" - L2-fehler u: ",L2fehler_u,"\n - max Abstand u: ",maxDiff_u,"\n - L2-fehler sigma: ",L2fehler_sigma,"\n - max Abstand sigma: ",maxDiff_sigma)

