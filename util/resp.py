#!/usr/bin/python
"""
 Fits atomic charges according to the RESP procedure
"""

import numpy, re

def read_potential(filename):
    """Read the electrostatic potential data"""
    # RE querries
    querry_np = re.compile("(\d*) points will be used for fitting atomic charges")
    querry_na = re.compile("NAtoms= *(\d*)")
    querry_pot= re.compile("Electrostatic Properties \(Atomic Units\).* Leave", re.DOTALL)
    querry_pts= re.compile("Electrostatic Properties Using The SCF Density.* Entering OneElI...", re.DOTALL)


    # read the file
    d = open(filename)
    data = d.read(); d.close()

    # find the number of points
    r = querry_np.findall(data)
    np = int(r[0])
    V = numpy.zeros((np, 4), numpy.float64)

    # find the number of atoms
    r = querry_na.findall(data)
    n_all_atoms = int(r[0])
    R = numpy.zeros((n_all_atoms, 3), numpy.float64)
   
    # find the potential points coordinates
    r = querry_pts.findall(data)
    r = r[0].split('\n')
    r_pot = r[4+n_all_atoms:-1]
    r_at  = r[4:4+n_all_atoms]
    assert len(r_pot) == np, "Error: number of potential points is not matching!"
    assert len(r_at)  == n_all_atoms, "Error: number of atoms do not match!"
 
    # find the potential values
    r = querry_pot.findall(data)
    r = r[0].split('\n')
    r = r[6+n_all_atoms:-2]
    assert len(r) == np, "Error: number of potential points is not matching!"
 
    # Fill the data in the R and V matrixes
    for i in range(n_all_atoms):
        R[i] = numpy.array(r_at[i].split()[-3:], numpy.float64)
    for i in range(np):
        V[i,:3] = numpy.array(r_pot[i].split()[-3:], numpy.float64)
        V[i, 3] = numpy.float64(r[i].split()[-1])
    return R, V


class RESP(object):
    """
 RESP Method.

 Usage:

 solver = RESP(R, V, Q)
 solver.run(a = 0.0005, b = 0.1, c = 'all')
 print solver

 where:

 g - degeneracy vector
 d - capping vector
 R - array of atomic positions
 V - array of position vectors of all electrostatic potential evaluation centres and their values
 G - grouping information

 a - hyperbolic restraint strength
 b - hyperbolic restraint tightness
 c - atoms which are restrained during fitting
"""

    def __init__(self, R, V, Q):
        """Assign the molecule properties: atoms, electrostatic potential and total charge"""
        self.__R = R
        self.__V = V
        self.__Q = Q
        self.__log = ""
        self.__nat= len(R)
        self.__np = len(V)

        self.__vec_rai, self.__vec_frai, self.__mat_rarb = self._compute_rarb(self.__V)

    def _compute_rarb(self, V):
        """Compute the atom-point and atom-atom potential-summed distance matrixes"""
        rai = numpy.array(self.__R)[:,numpy.newaxis] - numpy.array(self.__V[:,:3])[numpy.newaxis:]
        rai = rai*rai
        rai = numpy.sqrt(rai.sum(axis=2))
        rai = 1./rai
        frai = numpy.dot(rai, self.__V[:,3])
        rarb = numpy.dot(rai, rai.T)
        return rai, frai, rarb

    def __repr__(self):
        return str(self.__log)

    def fit(self, G, cap = False, resp = False, a = 0.0005, b = 0.10, conv = 1e-8):
        self._make_A_matrix(G, cap, resp, a, b, conv)
        return

    def _make_A_matrix(self, G, cap, resp, resp_a, resp_b, resp_conv):
        n_all_unique_atoms     = len(G["All Unique Atom Groups"])
        n_frozen_unique_atoms  = sum(G["Frozen Unique Atom Groups"])
        n_active_unique_atoms  = n_all_unique_atoms - n_frozen_unique_atoms

        # degeneracy vector for all unique groups
        vec_g_all_groups = numpy.array([ len(x) for x in G["Grouping"] ], dtype = numpy.float64)

        # charges of unique atoms
        vec_q_all_groups = numpy.array(G["Unique Charges"], dtype=numpy.float64)

        # indices
        idx_act = numpy.where(numpy.array(G["Frozen Unique Atom Groups"])== 0)[0]
        idx_fro = numpy.where(numpy.array(G["Frozen Unique Atom Groups"]) > 0)[0]

        # degeneracy vectors for active and frozen spaces
        vec_g_act = vec_g_all_groups[idx_act]
        vec_g_fro = vec_g_all_groups[idx_fro]

        # charges for active and frozen spaces
        vec_q_act = vec_q_all_groups[idx_act]
        vec_q_fro = vec_q_all_groups[idx_fro]

        idx_atom_groups = numpy.array(G["Grouping"])
        idx_atom_groups_act = idx_atom_groups[idx_act]
        idx_atom_groups_fro = idx_atom_groups[idx_fro]

        # capping vector for all atoms
        if cap:
           vec_c_all_groups = numpy.array( G["Capping"], dtype = numpy.float64 )
           vec_c_act        = vec_c_all_groups[idx_act]

        # allocate memory for matrices and vectors
        NW    = n_active_unique_atoms + 1 if not cap else n_active_unique_atoms + 2
        # H matrix
        mat_H = numpy.zeros( (NW, NW) , dtype = numpy.float64)
        # active charges and lagrange multipliers
        vec_q = numpy.zeros(  NW      , dtype = numpy.float64)
        # b (reference) vector
        vec_b = numpy.zeros(  NW      , dtype = numpy.float64)
        # frozen constraint vector
        vec_f = numpy.zeros(  NW      , dtype = numpy.float64)
        # restraint vector (for RESP)
        vec_r = numpy.zeros(  NW      , dtype = numpy.float64)

        # compute A matrix within H matrix
        for i in range(n_active_unique_atoms):
            g_i = vec_g_act[i]
            for j in range(n_active_unique_atoms):
                g_j = vec_g_act[j]
                mat_H[i,j] = g_i * g_j * self.__mat_rarb[idx_atom_groups_act[i]] [:,idx_atom_groups_act[j]] . sum()
        mat_H *= 2.0

        # finish with H matrix
        mat_H[:n_active_unique_atoms,n_active_unique_atoms] = vec_g_act
        mat_H[n_active_unique_atoms,:n_active_unique_atoms] = vec_g_act
        if cap:
           mat_H[:n_active_unique_atoms,n_active_unique_atoms+1] = vec_g_act * vec_c_act
           mat_H[n_active_unique_atoms+1,:n_active_unique_atoms] = vec_g_act * vec_c_act

        # inverse H matrix
        mat_Hinv = numpy.linalg.inv(mat_H)

        # build f (freezing constraint) vector
        for i in range(n_active_unique_atoms):
            g_i = vec_g_act[i]
            fi  = 0.0
            for u in range(n_frozen_unique_atoms):
                g_u = vec_g_fro[u]
                q_u = vec_q_fro[u]
                fi += g_i * g_u * q_u * self.__mat_rarb[idx_atom_groups_act[i]] [:,idx_atom_groups_fro[u]] . sum()
            vec_f[i] = fi
        vec_f *= 2.0

        # finish with f vector
        vec_f[n_active_unique_atoms] = numpy.dot(vec_q_fro, vec_g_fro) . sum()   # total frozen charge

        # build b (reference) vector
        for i in range(n_active_unique_atoms):
            g_i = vec_g_act[i]
            vec_b[i] = g_i * self.__vec_frai[idx_atom_groups_act[i]] . sum()
        vec_b *= 2.0

        # finish with b vector
        vec_b[n_active_unique_atoms] = self.__Q


        # Perform final fitting of the charges
        vec_x = numpy.dot(mat_Hinv, vec_b - vec_f)
        vec_q_act = vec_x[:n_active_unique_atoms].copy()
        print numpy.dot(vec_q_act, vec_g_act)
        if cap: print numpy.dot(vec_q_act*vec_c_act, vec_g_act)
        print vec_q_act

        # compute the modelled potential
        v_model =  numpy.zeros(self.__np)
        for a in range(n_active_unique_atoms):  # contribution from active charges
            v_model += vec_q_act[a] * vec_g_act[a] * self.__vec_rai[idx_atom_groups_act[a]] . sum(axis=0)
        for u in range(n_frozen_unique_atoms):  # contribution from frozen charges
            v_model += vec_q_fro[u] * vec_g_fro[u] * self.__vec_rai[idx_atom_groups_fro[u]] . sum(axis=0)
        rms_esp = self._rms(v_model, self.__V[:,3])
        print " ESP Done. ESP RMS: %f" % rms_esp


        # constraint vector and RESP 
        if resp:
           for i in range(n_active_unique_atoms):
               q_i      = vec_q_act[i] * vec_g_act[i]
               vec_r[i] = q_i / numpy.sqrt(q_i*q_i + resp_b*resp_b)
           vec_r *= resp_a
       
           vec_q_act_old = vec_q_act.copy()                                                          
           vec_q_act_new = numpy.dot(mat_Hinv, vec_b - vec_f + vec_r) [:n_active_unique_atoms]
           resp_rms = self._rms(vec_q_act_new, vec_q_act_old)
           while resp_rms > resp_conv:
                 for i in range(n_active_unique_atoms):
                     q_i      = vec_q_act_new[i] * vec_g_act[i]
                     vec_r[i] = q_i / numpy.sqrt(q_i*q_i + resp_b*resp_b)
                 vec_r[i] *= resp_a
                 vec_q_act_old = vec_q_act_new.copy()
                 vec_q_act_new = numpy.dot(mat_Hinv, vec_b - vec_f + vec_r) [:n_active_unique_atoms]
                 resp_rms = self._rms(vec_q_act_new, vec_q_act_old)
                 print resp_rms
                 print vec_q_act_new
                 

        # start iterative RESP procedure
        

        return 0

    def _rms(self, vec_new, vec_old):
        return numpy.sqrt( sum( (vec_new - vec_old)**2) )

    def run(self, a = 0.0005, b = 0.10, c = 'all'):
        """Perform single stage of RESP fit"""
        H = linalg
        return



