function compute_error_norms(u_num, u_exact, weights, semi)

	mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
	@unpack inverse_jacobian = cache.elements

	# Set up data structures
	l2_error = 0.0
	linf_error = copy(l2_error)
	total_volume = zero(real(mesh))

	# Iterate over all elements for error calculations
	for element in eachelement(solver, cache)

		for j in eachnode(solver), i in eachnode(solver)
			u_exact_loc = real(u_exact[i, j, element])
			u_num_loc = u_num[i, j, element]
			diff = u_exact_loc - u_num_loc

			volume_jacobian = abs(inv(cache.elements.inverse_jacobian[i, j, element]))

			l2_error += diff .^ 2 * (weights[i] * weights[j] * volume_jacobian)
			linf_error = @. max(linf_error, abs(diff))
			total_volume += weights[i] * weights[j] * volume_jacobian
		end
	end
	# For L2 error, divide by total volume
	l2_error = @. sqrt(l2_error / total_volume)

	return l2_error, linf_error
end

struct parameters_problem{RealT <: Real}
	H::RealT
	L::RealT
	R::RealT
	T0::RealT
	g::RealT
	delta::RealT
	ps::RealT
	xc::RealT
	d::RealT
	deltaT::RealT
	Nxmax::Int64
end

function parametercollect(Nxmax)
	H = 10000.0
	L = 300000.0
	R = 287.0
	T0 = 250.0
	g = 9.81
	delta = g / (T0 * R)
	ps = 100_000.0
	xc = 100_000.0
	d = 5000.0
	deltaT = 0.001
	return parameters_problem(H, L, R, T0, g, delta, ps, xc, d, deltaT, Nxmax)
end

function computeLF(alpha, beta, t)

	alpha_sqr = alpha * alpha
	beta_sqr  = beta * beta

	b_p_a = beta + alpha
	b_m_a = beta - alpha

	sin_at = sin(alpha * t)
	sin_bt = sin(beta * t)
	cos_at = cos(alpha * t)
	cos_bt = cos(beta * t)

	b_p_a = beta + alpha
	b_m_a = beta - alpha

	epsilon = 1.0e-20

	if (abs(b_p_a) > epsilon) && (abs(b_m_a) > epsilon)
		denom_ab = 1.0 / (b_p_a * b_m_a)

		is_alpha_null = abs(alpha) < 1.0e-30
		is_beta_null = abs(beta) < 1.0e-30

		if !is_alpha_null && !is_beta_null

			LT_m1 = (-cos_at / alpha_sqr + cos_bt / beta_sqr) * denom_ab + 1.0 / (alpha_sqr * beta_sqr)
			LT_0 = (sin_at / alpha - sin_bt / beta) * denom_ab

		elseif is_alpha_null && !is_beta_null

			b_t = beta * t
			LT_m1 = (0.5 * b_t^2 + cos(b_t) - 1.0) / beta^4
			LT_0 = (b_t - sin(b_t)) / beta^3

		elseif !is_alpha_null && is_beta_null

			a_t = alpha * t
			LT_m1 = (0.5 * a_t^2 + cos(a_t) - 1.0) / alpha^4
			LT_0 = (a_t - sin(a_t)) / alpha^3
		end

		LT_1 = (cos_at - cos_bt) * denom_ab
		LT_2 = (-sin_at * alpha + sin_bt * beta) * denom_ab
		LT_3 = (-cos_at * alpha^2 + cos_bt * beta^2) * denom_ab

	else
		if abs(alpha) > 1.0e-30

			alpha_sqr = alpha^2
			a_t = alpha * t
			cos_at = cos(a_t)
			sin_at = sin(a_t)

			LT_m1 = 1.0 / (alpha_sqr^2) * (-0.5 * a_t * sin_at - cos_at + 1.0)
			LT_0 = 0.5 / alpha_sqr * (sin_at / alpha - t * cos_at)
			LT_1 = 0.5 / alpha * t * sin_at
			LT_2 = 0.5 / alpha * (sin_at + a_t * cos_at)
			LT_3 = cos_at - 0.5 * a_t * sin_at

		else

			LT_3 = 1.0
			LT_2 = t
			LT_1 = LT_2 * t * 0.5
			LT_0 = LT_1 * t / 3.0
			LT_m1 = LT_0 * t * 0.25
		end
	end

	return LT_m1, LT_0, LT_1, LT_2, LT_3
end

function FourierSolution(kx, kz, rho_b0, t, pt, fac_kz)

	@unpack R, T0, g, delta, ps = pt
	f = 0
	rhos = ps / (T0 * R)
	cp = 1004.0
	cv = cp - R
	cs = sqrt(cp / cv * R * T0)

	p1 = cs^2 * (kx^2 + kz^2 + delta^2 / 4) + f^2
	q1 = g * kx^2 * (cs^2 * delta - g) + cs^2 * f^2 * (kz^2 + delta^2 / 4)

	alfa2 = p1 / 2 - sqrt(p1^2 / 4 - q1)
	beta2 = p1 / 2 + sqrt(p1^2 / 4 - q1)
	alpha = sqrt(alfa2)
	beta = sqrt(beta2)

	Lm1, L0, L1, L2, L3 = computeLF(alpha, beta, t)

	kd_p_m = 1im * kz - 0.5 * delta
	kd_p_p = 1im * kz + 0.5 * delta
	kd_m_p = -kd_p_m
    kd_m_m = -kd_p_p

	a2p = g - cs^2 * kd_p_p
	a2m = g - cs^2 * kd_m_p
	v_fac = g * rho_b0 / rhos

	hilf_pm = 1im * kx * L0 * v_fac
	u_FT_p = hilf_pm * a2p

	u = u_FT_p

	w_FT_p = -(L2 + (kx^2 * cs^2 + f^2) * L0) * v_fac

	hilf_pm = 1im * kx * L0 * v_fac
	# Calculate solutions
	u_bp = hilf_pm * a2p
	u_bm = -hilf_pm * a2m
	u_b = (u_bp * fac_kz + u_bm/fac_kz)
	w_b = -(L2 + (f^2 + cs^2 * kx^2) * L0) * g * rho_b0 / rhos

	rho_bp = (L3 + (p1 + g * kd_p_m) * L1 + (cs^2 * (kz^2 + delta^2 / 4) + g * kd_p_m) * f^2 * Lm1) * rho_b0
	rho_bm = - (L3 + (p1 + g * kd_m_m) * L1 + (cs^2 * (kz^2 + delta^2 / 4) + g * kd_p_m) * f^2 * Lm1) * rho_b0

	rho_b = (rho_bp * fac_kz + rho_bm / fac_kz)
	hilf_pm = (L1 + f^2 * Lm1)*g*rho_b0
	p_bp = -a2p * hilf_pm
	p_bm = a2m * hilf_pm
	p_b = (p_bp * fac_kz + p_bm/fac_kz)
	w_b = (w_FT_p * fac_kz -w_FT_p/fac_kz)
	return u_b, w_b, rho_b, p_b
end



function ComputeAnalyticalSolution( Nxmax,  semi, cells_per_dimension, polydeg, t)

	# dx e dz definisco la risoluzione dell'integrale
	# da automatizzare per ottenere la risoluzione richiesta

	# Nxmax e Nzmax definiscono il numero massimo di lunghezze d'onda nel dominio spettrale

	pt = parametercollect(Nxmax)

	@unpack L, H, delta, xc, ps, T0, R, deltaT, d = pt
	mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)

	ub = zeros(Complex{Float64}, polydeg + 1, polydeg + 1, prod(cells_per_dimension))

	wb = zeros(Complex{Float64}, polydeg + 1, polydeg + 1, prod(cells_per_dimension))

	pb = zeros(Complex{Float64}, polydeg + 1, polydeg + 1, prod(cells_per_dimension))

	rhob = zeros(Complex{Float64}, polydeg + 1, polydeg + 1, prod(cells_per_dimension))

	Tb = zeros(Complex{Float64}, polydeg + 1, polydeg + 1, prod(cells_per_dimension))

	thetab = zeros(Complex{Float64}, polydeg + 1, polydeg + 1, prod(cells_per_dimension))

	rhos = ps / (T0 * R)
	
	dkx = 2 * pi / L
	dkz = pi / H
	kz = dkz
	delta_rho = -rhos / T0 * deltaT
	@unpack node_coordinates = cache.elements
	kmax = Integer(Nxmax / 2)
	for element in eachelement(solver, cache)
		for j in eachnode(solver)
			for i in eachnode(solver)
				x_z_local = Trixi.get_node_coords(node_coordinates, equations, solver, i, j, element)
				xloc = x_z_local[1] - t * 20.0 - xc
				fac_kz = exp(1im * kz * x_z_local[2])
				breth_transform = exp(-0.5 * delta * x_z_local[2])
				p0 = ps * breth_transform*breth_transform
				rho0 = p0 / (R *T0)
				theta0 = T0*(ps/p0)^(R/1004)
				for nx in -kmax:kmax-1
					kx = dkx * nx
					fac_kx = exp(1im * kx * xloc)
					FT_gauss_x = d / (2.0 * sqrt(pi)) * exp(-0.25 * kx^2 * d^2)
					FT_sin_z = -0.5 * 1im
					rho_FT_0 = delta_rho * FT_gauss_x * FT_sin_z
					u_b, w_b, rho_b, p_b = FourierSolution(kx, kz, rho_FT_0, t, pt, fac_kz)
					phase = fac_kx * dkx
					ub[i, j, element] = ub[i, j, element] + u_b * phase / breth_transform
					wb[i, j, element] = wb[i, j, element] + w_b * phase / breth_transform
					pb[i, j, element] = pb[i, j, element] + p_b  * phase * breth_transform
					rhob[i, j, element] = rhob[i, j, element] + rho_b * phase * breth_transform			

				end
				Tb[i, j, element] = Tb[i, j, element] + T0 * (pb[i,j,element]/p0 - rhob[i,j,element]/rho0)
				thetab[i,j,element] = thetab[i,j,element] + theta0*(Tb[i,j,element] /T0 - R/1004 * pb[i,j,element]/p0)
			end
		end
	end
	
	return ub, wb, pb, rhob, Tb, thetab
end
