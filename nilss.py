import os
import struct
from subprocess import *

import f90nml
from numpy import *
from scipy import sparse
import scipy.sparse.linalg as splinalg

# ====================== time horizon parameters ============================= #
## time horizon start
m0 = 10000
## time segment size
dm = 1000
## number of time segments
K = 10
# number of homogeneous adjoints
p = 2
# ====================== time horizon parameters ============================= #

def primal():
    os.chdir('flow')
    call(["./flow"])
    with open('objective.bin', 'rb') as f:
        J = frombuffer(f.read(), dtype='>d')
    os.chdir('..')
    return J

def adjoint(t_strt,t_end,w_end,rhs=0):
    os.chdir('adj')
    with open('input.bin', 'wb') as f:
        f.write(asarray(w_end, dtype='>d').tobytes())
    call(["./adj",str(t_end),str(t_strt),str(rhs)])
    with open('output.bin', 'rb') as f:
        w = frombuffer(f.read(), dtype='>d')
    with open('dxdt.bin', 'rb') as f:
        dxdt = frombuffer(f.read(), dtype='>d')
    g = loadtxt("grad.txt")
    os.chdir('..')
    return w,g,dxdt

# system size (# of DoF)
def get_fun3d_dof():
    nml = f90nml.read('fun3d.nml')
    mesh_file = nml['project']['project_rootname'] + '.b8.ugrid'
    with open(mesh_file, 'rb') as f:
        n = struct.unpack('>I', f.read(4))[0]
    return n * 5

n = get_fun3d_dof()
n = int(open('flow/n').read())

# compute check points
t_chkpts = m0 + dm * arange(K+1)

# run primal
J = primal()
J_mean = J[m0:m0+dm*K].mean()

# Set Adjoint Terminal Conditions
#random.seed(12)
W = random.rand(n,p)
W = vstack([W, W])
[QT,RT] = linalg.qr(W)
W = QT

# Loop over all time segments, solve p homogeneous and 1 inhomogeneous
# adjoint on each

RTs = zeros([K,p,p])
bs = zeros([K+1,p])
gs = zeros([K,p])
h = 0.0

wh = random.rand(n)
wh = hstack([wh, wh])

dXdt = zeros([2*n,p])

_, _, dxdt = adjoint(t_chkpts[K],t_chkpts[K],wh)
for i in range(K-1,-1,-1):
    t_strt = t_chkpts[i]
    t_end = t_chkpts[i+1]
    print "segment {0}: [{1},{2}]".format(str(i), t_strt, t_end)

    dxdt_normalized = dxdt / linalg.norm(dxdt)
    P = eye(2*n) - outer(dxdt_normalized, dxdt_normalized)
    W = dot(P, W)
    wh = dot(P, wh) - (J[t_strt] - J_mean) * dxdt / dot(dxdt, dxdt)

    # solve homogeneous adjoints
    for j in range(p):
        W[:,j], gs[i,j], _ = adjoint(t_strt,t_end,W[:,j])

    # solve inhomogeneous adjoint
    wh, hi, dxdt = adjoint(t_strt,t_end,wh,1)
    h = h + hi

    dxdt_normalized = dxdt / linalg.norm(dxdt)
    P = eye(2*n) - outer(dxdt_normalized, dxdt_normalized)
    W = dot(P, W)
    wh = dot(P, wh) + (J[t_strt] - J_mean) * dxdt / dot(dxdt, dxdt)

    # QR decomposition
    [Q,R] = linalg.qr(W)
    RTs[i,:,:] = R.T
    bs[i,:] = dot(Q.T, wh)

    # Set terminal conditions for next segment
    W = Q
    wh = wh - dot(Q, bs[i,:])

    # print R, bs[i,:], wh

# form KKT system and solve
d = ones(K)
eyes = eye(p,p) / d[:,newaxis,newaxis]
B = -sparse.bsr_matrix((eyes, r_[1:K+1], r_[:K+1]))\
          + sparse.bsr_matrix((RTs, r_[:K],r_[:K+1]),\
          shape=(K*p, (K+1)*p))
B = B.tocsr()
A = B * B.T
rhs = -B * ravel(bs)

# print B.todense()

'''
import matplotlib.pyplot as plt

plt.figure()
plt.spy(A)
plt.show()
'''

alpha = splinalg.spsolve(A, rhs)

# compute sensitivities
T = 10
grad = (dot(alpha,ravel(gs)) + h) / T
print(grad)

v = loadtxt('flow/fd.txt')
plot(v[:,0], v[:,1], 'o')
plot([v[0,0], v[0,0] + 1], [v[0,1], v[0,1] + grad], '-')
