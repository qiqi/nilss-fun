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

def adjoint(t_strt,t_end,w_end,rhs=0):
    with open('input.bin', 'wb') as f:
        f.write(asarray(w_end, dtype='>d').tobytes())
    call(["./adj",str(t_end),str(t_strt),str(rhs)])
    with open('output.bin', 'rb') as f:
        w = frombuffer(f.read(), dtype='>d')
    g = loadtxt("grad.txt")
    return w,g

# system size (# of DoF)
def get_fun3d_dof():
    nml = f90nml.read('fun3d.nml')
    mesh_file = nml['project']['project_rootname'] + '.b8.ugrid'
    with open(mesh_file, 'rb') as f:
        n = struct.unpack('>I', f.read(4))[0]
    return n * 5

n = get_fun3d_dof()
print n
n = 3

# compute check points
t_chkpts = m0 + dm * arange(K+1)

# run primal
os.chdir("flow/")
call("./flow")

# Set Adjoint Terminal Conditions
#random.seed(12)
W = random.rand(n,p)
W = vstack([W, W])
[QT,RT] = linalg.qr(W)
W = QT

# Loop over all time segments, solve p homogeneous and 1 inhomogeneous
# adjoint on each

os.chdir("../adj/")

RTs = zeros([K,p,p])
bs = zeros([K+1,p])
gs = zeros([K,p])
h = 0.0

wh = random.rand(n)
wh = hstack([wh, wh])

for i in range(K-1,-1,-1):
    t_strt = t_chkpts[i]
    t_end = t_chkpts[i+1]
    print "segment {0}: [{1},{2}]".format(str(i), t_strt, t_end)

    # solve homogeneous adjoints
    for j in range(p):
        W[:,j], gs[i,j] = adjoint(t_strt,t_end,W[:,j])

    # solve inhomogeneous adjoint
    wh, hi = adjoint(t_strt,t_end,wh,1)
    h = h + hi

    # QR decomposition
    [Q,R] = linalg.qr(W)
    RTs[i,:,:] = R.T
    bs[i,:] = dot(Q.T, wh)

    # Set terminal conditions for next segment
    W = Q
    wh = wh - dot(Q, bs[i,:])

    # print R, bs[i,:], wh

os.chdir("..")

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
grad = dot(alpha,ravel(gs)) + h

print grad

