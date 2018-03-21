

module funcapp

	using ApproxFun
	using ApproXD
	using Base.Test
	using FastGaussQuadrature
	using Plots

	# use chebyshev to interpolate this:
	function q1(n=15)
		#define function
		f(x) = x + 2x^2 - exp(-x)
		#degree of approx
		deg_approx = n-1
		#interval bounds
		a, b = -3, 3 #endpoints

		# Compute chebyshev nodes
		zc = gausschebyshev(n) #compute nodes on [-1,1]
		xc = (a+b)/2 + (b-a)/2*zc[1] #map nodes to [-3,3] interval
		yc = [f(x) for x in xc] #function value at nodes
		# alternative = 'by hand'
		#xc = [(a+b)/2 + (b-a)/2*cos(pi*(2k-1)/2n) for k in 1:15]

		#Construct interpolation matrix Phi at z
		T(z,j) = cos.(acos.(z)*j)
		phi = [T(z,j) for z in zc[1], j in 0:n-1]

		#Coefficient vector
		c = phi \ yc

		#predict other values
		n_new = 100
		x_ab = linspace(a, b, n_new)
		z_ab = 2*(x_ab - a) / (b-a) - 1 #map to [-1:1] interval
		y_ab = [f(x) for x in x_ab]
		#basis function for new points
		phi_pred = [T(x,j) for x in z_ab, j in 0:n-1]
		y_hat = phi_pred * c # new predicted values
		dev = y_ab - y_hat # compute deviation

		#plot results
		p1 = plot(x_ab, y_ab, linewidth=2, linecolor=:red, label="f(x)", title = "Approximation of f")
		scatter!(x_ab, y_hat, markershape=:circle, markersize=1, markercolor = :black, label="approx f'(x)")
		p2 = plot(x_ab, dev, label="deviation", title = "Approximation error")
		plot(p1,p2,layout=(2,1))

		# Test approximation accuracy
		@test findmax(dev)[1] < 10.0^-9

		return Dict(:err => 1.0) ## WTH ??

	end

	function q2(n)
		f(x) = x + 2x^2 - exp(-x)
		n = 100
		x = linspace(-3, 3, n)

		#approximation using package
		y_hat = Fun(f, Chebyshev(-3..3))

		# true function
	    y = [f(x) for x in x]
		# deviation
		dev = y - y_hat

		#plot results
	    p3 = plot(x, y, label = "truth" )
		scatter!(x, y_hat, markershape=:circle, markersize=.5, markercolor = :black, label="approx.")
	    #p4 = plot(x, y - y_hat, yformatter = :scientific, label = "deviation")
		#plot(p3,p4,layout=(2,1))
	end


	# plot the first 9 basis Chebyshev Polynomial Basis Fnctions
	function q3()
		a, b = -3.0, 3.0#bounds
		x = linspace(a, b, 500)
		z(x) = 2*(x - a) / (b-a) - 1 #map to [-1,1] interval
		T(z,j) = cos.(acos.(z)*j)

		poly = [T(z(x),j) for x in x, j in 1:9]
		p5 = plot(x,poly,layout=(3,3),grid=false,legend=false,linewidth=2,linecolor=:black)
		return p5
	end

	ChebyT(x,deg) = cos(acos(x)*deg)
	unitmap(x,lb,ub) = 2.*(x.-lb)/(ub.-lb) - 1	#[a,b] -> [-1,1]

	type ChebyType
		f::Function # fuction to approximate
		nodes::Union{Vector,LinSpace} # evaluation points
		basis::Matrix # basis evaluated at nodes
		coefs::Vector # estimated coefficients

		deg::Int 	# degree of chebypolynomial
		lb::Float64 # bounds
		ub::Float64

		# constructor
		function ChebyType(_nodes::Union{Vector,LinSpace},_deg,_lb,_ub,_f::Function)
			n = length(_nodes)
			y = _f(_nodes)
			_basis = Float64[ChebyT(unitmap(_nodes[i],_lb,_ub),j) for i=1:n,j=0:_deg]
			_coefs = _basis \ y  # type `?\` to find out more about the backslash operator. depending the args given, it performs a different operation
			# create a ChebyType with those values
			new(_f,_nodes,_basis,_coefs,_deg,_lb,_ub)
		end
	end

	# function to predict points using info stored in ChebyType
	function predict(Ch::ChebyType,x_new)

		true_new = Ch.f(x_new)
		basis_new = Float64[ChebyT(unitmap(x_new[i],Ch.lb,Ch.ub),j) for i=1:length(x_new),j=0:Ch.deg]
		basis_nodes = Float64[ChebyT(unitmap(Ch.nodes[i],Ch.lb,Ch.ub),j) for i=1:length(Ch.nodes),j=0:Ch.deg]
		preds = basis_new * Ch.coefs
		preds_nodes = basis_nodes * Ch.coefs

		return Dict("x"=> x_new,"truth"=>true_new, "preds"=>preds, "preds_nodes" => preds_nodes)
	end

	function q4a(deg=(5,9,15),lb=-1.0,ub=1.0)
		deg = (5,9,15)
		lb = -1.0
		ub = 1.0
		f(x) = 1/(1+25x^2)

		u_nodes = [linspace(lb, ub, n+1) for n in deg]
		cheb_nodes = [gausschebyshev(n+1) for n in deg]

		# linspace approx
		#lin_approx = Any[]
		#for n in deg
		  #y_hat = ChebyType(collect(linspace(-1.0,1.0, 6)), 5, -1.0, 1.0, f)
		#end

		println("Could not solve the question")
	end

	function q4b()
		println("Could not solve this question")
	end

	function q5()
		f(x) = sqrt(sqrt(x^2))
		lb = -1.0
		ub = 1.0
		k = 13
		n = 65

		#plot true function
		x = linspace(lb,ub, n)
		plot(x, f, label ="true function")

		#spline with uniform knot vector
		bs_unif = BSpline(k, 3, lb, ub)
		#y_unif = [f(x) for ]
		B_unif = full(getBasis(collect(x),bs_unif)) #basis function
		c_unif= B_unif \ g(collect(y)) #coeff

		println("Could not solve the rest")
	end


		# function to run all questions
	function runall()
		println("running all questions of HW-funcapprox:")
		q1(15)
		q2(15)
		q3()
		q4a()
		q4b()
		q5()
	end


end
