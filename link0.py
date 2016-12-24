from sympy import *
import numpy
import time
start = time.time()
theta = symbols('theta')
a = 32 #input
b = 90 #intermediate
c = 52 #output
d = 90 #fixed)
gamma = atan(3/4)
pre_ans = numpy.radians(80)
csvfile = open('link.csv','a')
csvfile.write("a={0},b={1},c={2},d={3}\r\n".format(a,b,c,d))
csvfile.write("phi,theta,alpha,torque\r\n")
for phi_deg in range(78,103):
	phi = numpy.radians(phi_deg)
	s_gamma = sin(gamma)
	c_gamma = cos(gamma)
	cos_angle_abc = (a**2+b**2-(d*s_gamma/sin(phi)-c)**2-(d*s_gamma/tan(phi)+d*c_gamma)**2+2*(d*s_gamma/sin(phi)-c)*(d*s_gamma/tan(phi)+d*c_gamma))/(2*a*b)
	angle_abc = acos(cos_angle_abc)
	t_trad = tan(theta)
	s_angle_abc = sin(angle_abc)
	bc_prime = a - sqrt(((d*s_gamma-c*sin(phi))/(t_trad))**2+(d*s_gamma-c*sin(phi))**2)
	sin_alpha = (bc_prime*s_angle_abc)/(c*cos(phi)+d*c_gamma-(d*s_gamma-c*sin(phi)/(t_trad)))
	alpha = asin(sin_alpha)
	c_trad = cos(theta)
	l1 = (d*s_gamma/tan(phi)+d*c_gamma)
	l2 = (d*s_gamma/sin(phi)-c)
	L = sqrt(a**2+l1**2-2*a*l1*c_trad)
	x1 = (l2**2+L**2-b**2)/(2*l2*L)
	y1 = (l1**2+L**2-b**2)/(2*l1*L)
	x2 = sqrt(1-x1**2)
	y2 = sqrt(1-y1**2)
	#eq = x1*y1-x2*y2-cos(phi)
	#macluarin
	m_cp = 1 - phi**2/2 + phi**4/24 
	m_sp = phi - phi**3/6 + phi**5/120 
	m_tp = phi + phi**3/3 + 2*phi**5/15 
	m_ct = 1 - theta**2/2
	# + theta**4/24 
	m_st = theta - theta**3/6 + theta**5/120 
	m_tt = theta + theta**3/3 + 2*theta**5/15 
	eq2 = -sqrt(1 - (a**2 - 2*a*(d*sin(gamma)/m_tp + d*cos(gamma))*m_ct - b**2 + (-c + d*sin(gamma)/m_sp)**2 + (d*sin(gamma)/m_tp + d*cos(gamma))**2)**2/((-2*c + 2*d*sin(gamma)/m_sp)**2*(a**2 - 2*a*(d*sin(gamma)/m_tp + d*cos(gamma))*m_ct + (d*sin(gamma)/m_tp + d*cos(gamma))**2)))*sqrt(1 - (a**2 - 2*a*(d*sin(gamma)/m_tp + d*cos(gamma))*m_ct - b**2 + 2*(d*sin(gamma)/m_tp + d*cos(gamma))**2)**2/((2*d*sin(gamma)/m_tp + 2*d*cos(gamma))**2*(a**2 - 2*a*(d*sin(gamma)/m_tp + d*cos(gamma))*m_ct + (d*sin(gamma)/m_tp + d*cos(gamma))**2))) - m_cp + (a**2 - 2*a*(d*sin(gamma)/m_tp + d*cos(gamma))*m_ct - b**2 + 2*(d*sin(gamma)/m_tp + d*cos(gamma))**2)*(a**2 - 2*a*(d*sin(gamma)/m_tp + d*cos(gamma))*m_ct - b**2 + (-c + d*sin(gamma)/m_sp)**2 + (d*sin(gamma)/m_tp + d*cos(gamma))**2)/((-2*c + 2*d*sin(gamma)/m_sp)*(2*d*sin(gamma)/m_tp + 2*d*cos(gamma))*(a**2 - 2*a*(d*sin(gamma)/m_tp + d*cos(gamma))*m_ct + (d*sin(gamma)/m_tp + d*cos(gamma))**2))
	
	#newton's method loop
	diff_eq2 = diff(eq2,theta) 
	abs = 1
	f = eq2/diff_eq2
	theta_val = pre_ans
	epsilon = 10**(-1) #limit of error
	first_start = time.time()
	while abs > epsilon:
		appro = theta_val - f.subs([(theta,theta_val)])
		theta_ans = theta_val
		abs = numpy.abs(appro - theta_val)
		theta_val = abs
	theta_val = theta_ans
	first_end = time.time()
	first_stage = first_end - first_start
	print("phi:{0},1st:{1}".format(phi_deg,first_stage))
	epsilon = 10**(-2) #limit of error
	second_start = time.time()
	while abs > epsilon:
		appro = theta_val - f.subs([(theta,theta_val)]) 
		theta_ans = theta_val
		abs = numpy.abs(appro - theta_val)
		theta_val = abs
	theta_val = theta_ans
	second_end = time.time()
	second_stage = second_end - second_start
	print("phi:{0},2nd:{1}".format(phi_deg,second_stage))
	epsilon = 10**(-3) #limit of error
	third_start = time.time()
	while abs > epsilon:
		appro = theta_val - f.subs([(theta,theta_val)])
		theta_ans = theta_val
		abs = numpy.abs(appro - theta_val)
		theta_val = abs
	theta_val = theta_ans
	third_end = time.time()
	third_stage = third_end -third_start
	print("phi:{0},3rd:{1}".format(phi_deg,third_stage))
	"""
	epsilon = 10**(-4) #limit of error
	while abs > epsilon:
		appro = theta_val- f.subs([(theta,theta_val)])
		theta_ans = theta_val
		abs = numpy.abs(appro - theta_val)
		theta_val = abs
	theta_val = theta_ans
	"""
	pre_ans = theta_ans
	theta_output =180 * theta_ans/numpy.pi
	alpha_output = 180*(alpha.subs([(theta,theta_ans)]))/numpy.pi
	torque_output = (c*sin(phi-alpha.subs([(theta,theta_ans)])))/(a*sin(theta_ans+alpha.subs([(theta,theta_ans)])))
	#output data
	csvfile.write("{0}, {1},{2},{3}\r\n".format(phi_deg,theta_output,alpha_output,torque_output))	
csvfile.close()
elapsed_time = time.time() -start
print ("elapsed_time:{0}".format(elapsed_time))
