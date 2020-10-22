import numpy as np

def rotation_matrix(theta, r):
	# returns matrix R, such that Ru is a rotation of angle theta about r.
	R = np.zeros((3, 3))
	# numpy [2, 0] is down column; eg [:, 1] is the 2nd column.
	S = np.sin(theta)
	C = np.cos(theta)
	t = 1 - C
	ux = r[0]
	uy = r[1]
	uz = r[2]

	R[0,0] = t * ux * ux + C
	R[1,0] = t * ux * uy + S * uz
	R[2,0] = t * ux * uz - S * uy
	R[0,1] = t * ux * uy - S * uz
	R[1,1] = t * uy * uy + C
	R[2,1] = t * uy * uz + S * ux
	R[0,2] = t * ux * uz + S * uy
	R[1,2] = t * uy * uz - S * ux
	R[2,2] = t * uz * uz + C
	return R

def vrot(A, B, C, theta, r):
	# performs a rotation of [A B C] about r, where all A, B, C are column vectors
	R = rotation_matrix(theta, r)
	u = np.zeros((3, 3))
	u[:,0] = A
	u[:,1] = B
	u[:,2] = C
	v = np.matmul(R, u)
	return v[:,0], v[:,1], v[:,2]


def get_fixed_axes(t):
    global coords
    global view
    A = coords[t, 0] - ((coords[t, 2] + coords[t, 3]) / 2.) # eg vector between CG and the midpoint of the epsilons
    A = A / np.linalg.norm(A) # normalise
    D = coords[t, 1] - coords[t, 4] # eg vector across
    D = D / np.linalg.norm(D)
    B = np.cross(A, D) # cross product gives out of plane vector
    B = B / np.linalg.norm(B)
    C = np.cross(A, B) # cross gives in plane vector perpendicular to CG-epsilon
    C = C / np.linalg.norm(C)
    mid_all =  coords[t, 0] + coords[t,1] + coords[t,2]+ coords[t,3]+ coords[t,4]
    mid_all = mid_all / 5.
    
    return [A, C, B, mid_all] # A C B so X is parallel to ring axis, Z is perpendicular to plane

def get_rotated_vector(t, A, B, C, theta, r):
    # returns X, Y, Z, corresponding to A, B, C rotated about r an amount theta
    x, y, z = calc.vrot(A, B, C, theta, r)
    return [x, y, z]

def draw_fixed_axes(t):
    global coords
    global view
    [A, B, C, mid_all ] = get_fixed_axes(t)
    m_A = mid_all + A * 3
    m_B = mid_all + B * 3
    m_C = mid_all + C * 3
    print("Adding arrows")
    print(mid_all)
    view.shape.add_arrow(mid_all.tolist(), m_A.tolist(), [1, 0, 0], 0.2)
    view.shape.add_arrow(mid_all.tolist(), m_B.tolist(), [0, 1, 0], 0.2)
    view.shape.add_arrow(mid_all.tolist(), m_C.tolist(), [0, 0, 1], 0.2)
    return [A, B, C, mid_all]

def draw_rotated_vector(D, mid):
    global view
    m_D = mid + D * 3.
    view.shape.add_arrow(mid.tolist(), m_D.tolist(), [0, 1, 1], 0.2)
    
def draw_cart_vector(X, Y, Z, mid):
    global view
    m_X = mid + X * 2.
    m_Y = mid + Y * 2.
    m_Z = mid + Z * 2.
    view.shape.add_arrow(mid.tolist(), m_X.tolist(), [1, 0, 0], 0.2)
    view.shape.add_arrow(mid.tolist(), m_Y.tolist(), [0, 1, 0], 0.2)
    view.shape.add_arrow(mid.tolist(), m_Z.tolist(), [0, 0, 1], 0.2)
    
def get_efg_tensors(t, theta):
    # returns the X, Y, Z components for time t, rotated about the z axis an amount theta
    [Xm, Ym, Zm, mid] = get_fixed_axes(t)
    [X, Y, Z] = get_rotated_vector(t, Xm, Ym, Zm, theta, Zm)
    return [X, Y, Z]
    
def get_order_parameters(k, theta):
    # k bounded by 1 to n_iter/2.
    global coords
    sf = 10
    q = np.shape(coords)
    n_iter = int(q[0] / float(sf))
    mul = (1. / float((n_iter - k)))
    su = np.array([0, 0, 0])
    for i in range(0, n_iter - k):
        [Xik, Yik, Zik] = get_efg_tensors((i + k) * sf, theta)
        [Xi, Yi, Zi] = get_efg_tensors(i * sf, theta)
        dp = np.array([0.,0.,0.])
        dp[0] = np.dot(Xik, Xi)
        dp[1] = np.dot(Yik, Yi)
        dp[2] = np.dot(Zik, Zi)
        #print(Xik)
        #print(Xi)
        #print(np.dot(Xik, Xi))
        #print(dp[0])
        #print("\n\n")
        comp = (3 * dp * dp - 1) / 2.
        su = su + comp
    order_params = mul * su
    return order_params
    

def apply_wigner(X, Y, Z, alpha, beta, gamma):
	X, Y, Z = vrot(X, Y, Z, alpha, Z)
	X, Y, Z = vrot(X, Y, Z, beta, Y)
	X, Y, Z = vrot(X, Y, Z, gamma, Z)
	return X, Y, Z
