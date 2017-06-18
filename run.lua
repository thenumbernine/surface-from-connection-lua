#! /usr/bin/env luajit
require 'ext'
local bit = require 'bit'
local gl = require 'gl'
local glCall = require 'gl.call'
local matrix = require 'matrix'
local complex = require 'symmath.complex'
local symmath = require 'symmath'
local gnuplot = require 'gnuplot'
local ImGuiApp = require 'imguiapp'
local View = require 'glapp.view'
local Orbit = require 'glapp.orbit'


local function int_fe(r, q, dq_dr, dr)
	return q + dq_dr(r, q) * dr, r + dr
end

local function int_rk4(r, q, dq_dr, dr)
	local k1 = dq_dr(r, q)
	local k2 = dq_dr(r + .5 * dr, q + .5 * dr * k1)
	local k3 = dq_dr(r + .5 * dr, q + .5 * dr * k2)
	local k4 = dq_dr(r + dr, q + dr * k3)
	return q + (k1 + 2 * k2 + 2 * k3 + k4) * dr / 6, r + dr
end

local function permutations(args)
	local parity = args.parity or (#args.elements%2==0 and 1 or -1)
	local p = args.elements
	local callback = args.callback
	local index = args.index
	if args.size then
		if index and #index == args.size then
			return callback(index, parity)
		end
	else
		if #p == 0 then
			return callback(index, parity)
		end
	end
	for i=1,#p do
		local subset = table(p)
		local subindex = table(index)
		subindex:insert(subset:remove(i))
		parity = parity * -1		-- TODO times -1 or times the distance of the swap?
		if permutations{
			elements = subset, 
			callback = callback, 
			index = subindex,
			size = args.size,
			parity = parity,
		} then 
			return true 
		end
	end
end

function matrix.det(m)
	local dim = m:size()
	local volume = table.combine(dim, function(a,b) return a * b end) or 0
	if volume == 0 then return 1 end
	if #dim ~= 2 then error("determinant only works on degree-2 matrices") end
	if dim[1] ~= dim[2] then error("determinant only works on square matrices") end

	local n = dim[1]
	if n == 1 then return m[1][1] end

	-- any further -- use recursive case
	local result = 0
	permutations{
		elements = range(n),
		callback = function(index, parity)
			local product
			for i=1,n do
				local entry = m[i][index[i]]
				if not product then
					product = entry
				else
					product = product * entry
				end
			end
			result = result + parity * product
		end,
	}
	return result
end
--matrix.__tostring = tolua

function matrix:ident()
	assert(#self == 2 and self:degree() == 1)
	return self:lambda(function(i,j)
		return i == j and 1 or 0
	end)
end

function matrix:inv()
	local size = self:size()
	assert(#size == 2, "must be a matrix, not a vector or higher dimension array")
	assert(size[1] == size[2], "must be a square matrix, not a rectangular matrix")
	local n = size[1]
	if n == 2 then
		local a,b = self[1]:unpack()
		local c,d = self[2]:unpack()
		local det = a * d - b * c 
		return matrix{
			{d, -b},
			{-c, a}
		} / det
	elseif n == 3 then
		-- transpose, +-+- sign stagger, for each element remove that row and column and 
		return matrix{
			{self[2][2]*self[3][3]-self[2][3]*self[3][2], self[1][3]*self[3][2]-self[1][2]*self[3][3], self[1][2]*self[2][3]-self[1][3]*self[2][2]},
			{self[2][3]*self[3][1]-self[2][1]*self[3][3], self[1][1]*self[3][3]-self[1][3]*self[3][1], self[1][3]*self[2][1]-self[1][1]*self[2][3]},
			{self[2][1]*self[3][2]-self[2][2]*self[3][1], self[1][2]*self[3][1]-self[1][1]*self[3][2], self[1][1]*self[2][2]-self[1][2]*self[2][1]},
		} / self:det()
	else
		error"idk how to invert this"
	end
end

-- used for providing initial values for the metric
-- and verifying accuracy of numerical calculations
local Geometry = class()

function Geometry:init(app)
	self.app = app
	self.xmin = matrix(self.xmin)
	self.xmax = matrix(self.xmax)

	self.coordVars = table.map(self.coords, function(name)
		return symmath.var(name)
	end)

	local Tensor = symmath.Tensor
	do --if self.coordVars then
		Tensor.coords{{variables=self.coordVars}}
		local n = #self.coords
		local g = self.createMetric and self:createMetric() or Tensor('_ab', table.unpack(symmath.Matrix.identity(n,n)))
		local gU = Tensor('^ab', table.unpack( (symmath.Matrix.inverse(g)) ))
		local dg = Tensor('_abc', table.unpack(g'_ab,c'()))
		local ConnL = Tensor('_abc', table.unpack( ((dg'_abc' + dg'_acb' - dg'_bca')/2)() ))
		local Conn = Tensor('^a_bc', table.unpack( (gU'^ad' * ConnL'_dbc')() ))	
		self.calc = {
			g = self:compileTensor(g),
			Conn = self:compileTensor(Conn),
		}
	end
end

function Geometry:testExact()
	local app = self.app
	local exactConns = self:calc_conns()
	local connNumConnDiff = app.size:lambda(function(i,j)
		return (app.conns[i][j] - exactConns[i][j]):norm()
	end)

	local x = matrix{app.size[1]-2}:lambda(function(i) return app.xs[i+1][1][1] end)
	local y = matrix{app.size[2]-2}:lambda(function(j) return app.xs[1][j+1][2] end)
	local z = (app.size-2):lambda(function(i,j) return connNumConnDiff[i+1][j+1] end)
	gnuplot{
		output = 'conn numeric vs analytic.png',
		style = 'data lines',
		xlabel = self.coords and self.coords[1] or 'x^1',
		ylabel = self.coords and self.coords[2] or 'x^2',
		--log = 'z',
		griddata = {x = x, y = y, z},
		{splot=true, using='1:2:3', title = 'Γ^a_b_c |analytic - numeric|'},
	}
	print('max Γ^a_bc |analytic - numeric|', z:normLInf())
end

function Geometry:compileTensor(expr)
	local size = table.map(expr:dim(), function(i) return i.value end)
	local fs = matrix(size):lambda(function(...)
		return expr[{...}]:compile(self.coordVars)
	end)
	return function(x)
		return matrix(size):lambda(function(...)
			return fs[{...}]( table.unpack(x) )
		end)
	end
end

function Geometry:calc_gs()
	if not self.calc.g then return end
	local app = self.app
	return app.size:lambda(function(...)
		return self.calc.g(self.app.xs[{...}])
	end)
end

-- this is only used for comparing the error between the numerical and the exact connections
function Geometry:calc_conns()
	if not self.calc.Conn then return end
	return self.app.size:lambda(function(...)
		return self.calc.Conn(self.app.xs[{...}])
	end)
end


--[[
what each geometry subclass needs:
* coords, to determine the dimension of the surface
* xmin, xmax
* startCoord (where to start the surface generation at.  this is very important.
* one of the two:
	* createMetric, if the geometry is to use an analytical metric 
	* calc_conns.  otherwise these return identity g_ab's and zero Conn^a_bc's.
--]]


-- 2D geometries


local PolarHolGeom = class(Geometry)
PolarHolGeom.coords = {'r', 'θ'}
-- [[ initialize our basis of e=I at r=1 ...
PolarHolGeom.xmin = matrix{1, 0}
PolarHolGeom.xmax = matrix{10, 2 * math.pi}
PolarHolGeom.startCoord = {1,0}
--]]
--[[ applying a rotation, e=I at r=1 and theta!=0 still produces the same shape.
PolarHolGeom.xmin = matrix{1, 0}
PolarHolGeom.xmax = matrix{10, 2 * math.pi}
PolarHolGeom.startCoord = {1, math.pi/2}
--]]
--[[ ...otherwise the shape gets messed up -- for r=2
PolarHolGeom.xmin = matrix{1, 0}
PolarHolGeom.xmax = matrix{10, 2 * math.pi}
PolarHolGeom.startCoord = {2,0}
--]]
--[[ ...otherwise the shape gets messed up -- for r=1/2 (needs the range readjusted so rmin isn't 1)
PolarHolGeom.xmin = matrix{.1, 0}
PolarHolGeom.xmax = matrix{2, 2 * math.pi}
PolarHolGeom.startCoord = {.5,0}
--]] 
function PolarHolGeom:createMetric()
	local r, theta = self.coordVars:unpack()
	return symmath.Tensor('_ab', {1, 0}, {0, r^2})
end


-- The thing about non-holonomic geometry is
-- it needs commutation coefficients as well.
-- Otherwise how does it know how far to integrate the geodesics
-- to get to the next coordinate location?
-- This information is typically stored in the metric of the holonomic coordinate map.
local PolarNonHolGeom = class(Geometry)
PolarNonHolGeom.coords = {'r', 'θ'}
PolarNonHolGeom.xmin = matrix{1, 0}
PolarNonHolGeom.xmax = matrix{10, 2 * math.pi}
function PolarNonHolGeom:calc_conns()
	return self.app.size:lambda(function(i,j)
		local r = self.app.xs[i][j][1]
		-- Γ^θ_rθ = -Γ^r_θθ = 1/r
		return matrix{ {{0,0},{0,-1/r}}, {{0,1/r},{0,0}} }
	end)
end


-- cyl surface can't be reconstructed
-- because it needs extrinsic curvature information
-- it has a connection of zero

-- sphere surface likewise is a 2 dimensional system inside 3 dimensions
-- like cyl surface, it needs extrinsic curvature information to be properly rebuilt 
local SphereSurfaceHolGeom = class(Geometry)
local eps = .01
SphereSurfaceHolGeom.coords = {'θ', 'phi'}
SphereSurfaceHolGeom.xmin = matrix{eps, eps}
SphereSurfaceHolGeom.xmax = matrix{math.pi-eps, 2*math.pi-eps}
SphereSurfaceHolGeom.startCoord = {math.pi/2, math.pi}
function SphereSurfaceHolGeom:createMetric()
	local theta, phi = self.coordVars:unpack()
	local r = 1
	return symmath.Tensor('_ab', {r^2, 0}, {0, r^2 * symmath.sin(theta)^2})
end


local PoincareDisk = class(Geometry)
PoincareDisk.coords = {'u', 'v'}
PoincareDisk.xmin = {-1, -1}
PoincareDisk.xmax = {1, 1}
PoincareDisk.startCoord = {0, 0}
function PoincareDisk:createMetric()
	local u, v = self.coordVars:unpack()
	return symmath.Tensor('_ab', {4 / (1 - u^2 - v^2), 0}, {0, 4 / (1 - u^2 - v^2)})
end

-- 3D geometries


local CylHolGeom = class(Geometry)
CylHolGeom.coords = {'r', 'θ', 'z'}
CylHolGeom.xmin = {1, 0, -5}
CylHolGeom.xmax = {10, 2*math.pi, 5}
CylHolGeom.startCoord = {1,math.pi,0}
function CylHolGeom:createMetric()
	local r, theta, z = self.coordVars:unpack()
	return symmath.Tensor('_ab', {1, 0, 0}, {0, r^2, 0}, {0, 0, 1})
end


local SphereHolGeom = class(Geometry)
local eps = .01	-- the holonomic connection gets singularities (cot theta = inf) at the boundaries
				-- this could be avoided if the metric was evaluated at grid centers instead of vertices.
SphereHolGeom.coords = {'r', 'θ', 'φ'}
SphereHolGeom.xmin = {1, eps, -math.pi + eps}
SphereHolGeom.xmax = {10, math.pi - eps, math.pi - eps}
SphereHolGeom.startCoord = {1, math.pi/2, 0}
function SphereHolGeom:createMetric()
	local r, theta, phi = self.coordVars:unpack()
	return symmath.Tensor('_ab', {1,0,0}, {0, r^2, 0}, {0, 0, r^2 * symmath.sin(theta)^2})
end


local function I(x)
	return function()
		return x
	end
end


local App = class(Orbit(View.apply(ImGuiApp)))

App.title = 'reconstruct surface from geodesics' 
App.viewDist = 10

function App:initGL()
	gl.glEnable(gl.GL_DEPTH_TEST)

	-- 2D
	self.geom = PolarHolGeom(self)
	--self.geom = PolarNonHolGeom(self)
	--self.geom = SphereSurfaceHolGeom(self)	-- in absence of extrinsic curvature, this creates a polyconic projection
	--self.geom = ParabolaHolGeom(self)
	--self.geom = PoincareDisk(self)
	-- 3D
	--self.geom = CylHolGeom(self)
	--self.geom = SphereHolGeom(self)

	local n = #self.geom.coords
	self.size = matrix{n}:lambda(I( ({[2]=64, [3]=16})[n] ))

	self.xmin = self.geom and self.geom.xmin or matrix{n}:lambda(I(-1))
	self.xmax = self.geom and self.geom.xmax or matrix{n}:lambda(I(1))
	
	local n = #self.size
	
	self.view.ortho = n == 2
	
	--[[ cell centered, including borders
	self.dx = matrix{n}:ones():emul(self.xmax - self.xmin):ediv(self.size)
	self.xs = self.size:lambda(function(...) return (matrix{...} - .5):emul(self.dx) + self.xmin end)
	--]]
	-- [[ vertex centered, excluding borders, so position 2,2 is at xmin (useful for centering the corner vertex)
	self.dx = matrix{n}:ones():emul(self.xmax - self.xmin):ediv(self.size-3)
	self.xs = self.size:lambda(function(...) return (matrix{...} - 2):emul(self.dx) + self.xmin end)
	--]]

	-- if we are calculating the connection from discrete derivatives of the metric ...
	if self.geom.calc_gs then
		local gs = self.geom:calc_gs() 
		--or self.size:lambda(function(...) return matrix{n,n}:ident() end)
		local gUs = self.size:lambda(function(...)
			return gs[{...}]:inv()
		end)
		local dgs = self.size:lambda(function(...)
			local i = matrix{...}
			-- dg_abc = dg_bc/dx^a
			local dg = matrix{n}:lambda(function(a)
				-- using a constant border
				local ip = matrix(i) ip[a] = math.min(ip[a]+1, self.size[a])
				local im = matrix(i) im[a] = math.max(im[a]-1, 1)
				-- using a first-order derivative
				return (gs[ip] - gs[im]) / (self.xs[ip][a] - self.xs[im][a])
			end)
			-- reorientate dg_abc = g_ab,c
			return matrix{n,n,n}:lambda(function(a,b,c) return dg[c][a][b] end)
		end)

		self.conns = self.size:lambda(function(...)
			local i = matrix{...}
			local dg = dgs[i]
			local numConnLower = matrix{n,n,n}:lambda(function(a,b,c)
				return .5 * (dg[a][b][c] + dg[a][c][b] - dg[b][c][a])
			end)
			local gU = gUs[i]
			local check1 = gU * numConnLower
			local check2 = matrix{n,n,n}:lambda(function(a,b,c)
				local s = 0
				for d=1,n do
					s = s + gU[a][d] * numConnLower[d][b][c]
				end
				return s
			end)
			local err = (check1 - check2):norm()
			if err ~= 0 
			and err == err	-- skip the nans
			then
				print('check1')
				print(check1)
				print('check2')
				print(check2)
				error('norms differ by '..err)
			end
			return check2
		end)
		if self.geom and self.geom.calc_conns then
			self.geom:testExact()
		end
	else
		-- if calc_gs isn't there, then rely on calc_conns for providing the initial connections
		self.conns = self.geom:calc_conns()
	end

	-- embedded space position
	self.Xs = self.size:lambda(function(...)
		return matrix{n}:zeros()
	end)
	self.es = self.size:lambda(function(...)
		return matrix{n,n}:zeros()
	end)

	local xInit = matrix(assert(self.geom.startCoord))
	local i = ((xInit - self.geom.xmin):ediv(self.geom.xmax - self.geom.xmin):emul(self.size-2) + 2.5):map(math.floor)
	self.es[i] = matrix{n,n}:ident()
	local function withinEdge(index,size)
		for i=1,n do
			if index[i] <= 1 or index[i] >= size[i] then return false end
		end
		return true
	end

print'building surface...'
	-- now to reconstruct the es based on the conns ...
	-- [=[ flood fill 
	local todo = table{i}
	local sofar = {[tostring(index)] = true}
	while #todo > 0 do
		local index = todo:remove(1)
		local conn = self.conns[index]
		-- for each direction ...
		local _ = matrix.index
		for k=1,n do
			local connk = conn(_,_,k)
	
			for dir=-1,1,2 do 
				local nextIndex = matrix(index)
				nextIndex[k] = nextIndex[k] + dir
				-- skip the edges
				if withinEdge(nextIndex, self.size)
				and not sofar[tostring(nextIndex)]
				then
					local nextConnK = self.conns[nextIndex](_,_,k)
				
					-- derivative along coordinate
					local ds = dir * self.dx[k]
				
					local eOrig = self.es[index]
					local XOrig = self.Xs[index]
					
					local e = matrix(eOrig)
					local X = matrix(XOrig)
					
					--[[ forward-euler
					e = e + (e * connk) * ds
					X = X + e(_,k) * ds
					--]]
					-- [[ rk4 ...
					e = int_rk4(0, e, function(s, e)
						-- treating connections as constant
						--return e * connk
						-- interpolating connections between cells
						return e * (nextConnK * s + connk * (ds - s)) / ds
						-- using analytical solution
						--local x = self.xs[index]
						--local xNext = self.xs[nextIndex]
						--return e * self.geom.calc.Conn( (x * s + xNext * (ds - s)) / ds  )(_,k)
					end, ds)
					X = int_rk4(0, X, function(s, X)
						local f = s / ds
						-- constant:
						--return eOrig(_,k)
						-- linear interpolation
						return e(_,k) * f + eOrig(_,k) * (1 - f)
					end, ds)
					--]]
					--[[
					do
						-- if d/ds e = conn e then ...
						-- then we can solve this as a linear dynamic system!
						-- even though this assumes conn is constant (which it's not)
						-- for constant connections (like integrating polar coordinates around theta) this is no problem
						-- de/e = conn ds
						-- e = C exp(conn s) = C Q exp(λ s) Q^-1
						--
						-- 2D eigenvalues
						-- (a-λ)(d-λ) - bc = 0 <=> λ^2 - (a+d) λ + (ad-bc) = 0 <=> λ = 1/2 ( (a+d) +- sqrt( a^2 + 2ad + d^2 - 4 (ad - bc) ) )
						-- λ = (a+d)/2 +- sqrt( ((a-d)/2)^2 + bc ) )
						--print('conn '..self.geom.coords[k]..':\n'..connk)
						local a,b = connk[1]:unpack()
						local c,d = connk[2]:unpack()
						local asym = (a - d) / 2
						local discr = asym^2 + b * c
						
						local de

						-- [a b]
						-- [c a]
						if a == d and b ~= 0 and c ~= 0 then
							if b * c >= 0 then
								local l1, l2 = a + sd, a - sd
								local evR = matrix{{1, 1}, {math.sqrt(c/b), -math.sqrt(c/b)}}
								local evL = matrix{{1, math.sqrt(b/c)}, {1, -math.sqrt(b/c)}} / 2
								de = function(ds)
									return evR * matrix{
										{math.exp(ds * l1), 0},
										{0, math.exp(ds * l2)}
									} * evL
								end
							else
								-- b c < = means either b or c < 0 but not both
								local theta = math.sqrt(-b * c)
								local ratio = math.sqrt(-b / c)
								de = function(ds)
									local costh = math.cos(theta * ds)
									local sinth = math.sin(theta * ds)
									return matrix{
										{costh, -sinth / ratio},
										{sinth / ratio, costh}
									} * math.exp(a * ds)
								end
							end
						-- [a b]
						-- [0 a] for a real and b nonzero 
						elseif a == d and b ~= 0 then
							-- TODO solve this without eigen-decomposition
							error("defective matrix "..connk)
						-- [a 0]
						-- [0 d] for a,b any real
						elseif b == 0 and c == 0 then
							local l1, l2 = a, d
							local evR = matrix{{1,0},{0,1}}
							local evL = matrix(evR)
							de = function(ds)
								return evR * matrix{
									{math.exp(ds * l1), 0},
									{0, math.exp(ds * l2)}
								} * evL
							end
						-- [a b]
						-- [0 d]
						elseif c == 0 then
							local l1, l2 = a, d
							local evR = matrix{{1, 0}, {b, d-a}}:T()
							local evL = matrix{{1, b/(a-d)}, {0, -1/(a-d)}}
							de = function(ds)
								return evR * matrix{
									{math.exp(ds * l1), 0},
									{0, math.exp(ds * l2)}
								} * evL						
							end
						-- [a 0]
						-- [c d]
						elseif b == 0 then
							local l1, l2 = a, d
							local evR = matrix{{a-d, c}, {0, 1}}:T()
							local evL = matrix{{1/(a-d), 0}, {-c/(a-d), 1}}
							de = function(ds)
								return evR * matrix{
									{math.exp(ds * l1), 0},
									{0, math.exp(ds * l2)}
								} * evL						
							end
						elseif discr == 0 then	-- means (a-d)^2 = 4*b*c, then we have multiplicity 2
							error"here"
						elseif discr > 0 then
							-- 2D eigenvectors using the smaller eigenvalue
							-- [a-λ, b]    [ a - (a+d)/2 + sqrt( ((a-d)/2)^2 + bc ),    b ]    [ (a-d)/2 + sqrt( ((a-d)/2)^2 + bc ), b ]
							-- [c, d-λ] => [ c,    d - (a+d)/2 + sqrt( ((a-d)/2)^2 + bc ) ] => [ c, (d-a)/2 + sqrt( ((a-d)/2)^2 + bc ) ]
							--
							-- [ (sqrt( ((a-d)/2)^2 + bc ) + (a-d)/2)(sqrt( ((a-d)/2)^2 + bc ) - (a-d)/2), b((d-a)/2 + sqrt( ((a-d)/2)^2 + bc )) ]
							-- [ b c, b((d-a)/2 + sqrt( ((d-a)/2)^2 + bc )) ]
							--
							-- c x + ((d-a)/2 + sqrt( ((a-d)/2)^2 + bc )) y = 0
							-- y = t
							-- x = 1/c ((a-d)/2 - sqrt( ((a-d)/2)^2 + bc )) t
							local avg = (a + d) / 2
							local sd = complex.sqrt(discr) 
							local l1, l2 = avg + sd, avg - sd
							local evR = matrix{{asym + sd, asym - sd}, {c, c}}
							local evL = matrix{
								{-1, (asym + sd)/c},
								{1, (-asym + sd)/c},
							} / (2 * sd)
							de = function(ds)
								return evR * matrix{
									{math.exp(ds * l1), 0},
									{0, math.exp(ds * l2)}
								} * evL
							end
						else -- discr < 0	-- complex eigenvectors
							error"here"
						end
						
						X = int_rk4(0, X, function(s) return (e * de(s))(_,k) end, ds)
						e = e * de(ds)
					end
					--]]
					
					--[[ normalize columns
					e = e:T()
					for k=1,n do
						e[k] = e[k] / e[k]:norm()
					end
					e = e:T()
					--]]

--[[
print()
print('index='..index)
--print('self.dx='..self.dx)
print('nextIndex='..nextIndex)
print('ds='..ds)
print('x='..self.xs[index])
print('conn[k]='..connk)
print('XOrig='..XOrig)
print('X='..X)
print('eOrig='..eOrig)
print('e='..e)
--]]

					self.es[nextIndex] = e
					self.Xs[nextIndex] = X

--print('e2 from '..eOrig(_,2)..' to '..e(_,2)..' changing by '..(eOrig(_,2) - e(_,2)):norm())
--print('|e2| from '..eOrig(_,2):norm()..' to '..e(_,2):norm()..' changing by '..(e(_,2):norm() - eOrig(_,2):norm()))
--print('X from '..XOrig..' to '..X..' changing by '..(X - XOrig):norm())
					
					todo:insert(nextIndex)
					--todo:insert(math.random(#todo+1), nextIndex)
					--todo:insert(math.floor((#todo+1)/2), nextIndex)
					-- the linear dynamic system method only works right for polar coordinates when we use this order
					--todo:insert(1, nextIndex)
					sofar[tostring(nextIndex)] = true
				end
			end	
		end
	end
	--]]
end

local function glColor(m)
	if #m == 2 then
		gl.glColor3d(m[1], m[2], .5)
	elseif #m == 3 then
		gl.glColor3d(m:unpack())
	else
		error"can't color this many dimensions"
	end
end

local function glVertex(m)
	assert(({
		[2] = gl.glVertex2d,
		[3] = gl.glVertex3d,
	})[#m])(m:unpack())
end

function App:drawGrid()
	local xmin, xmax, ymin, ymax = self.view:getBounds(self.width / self.height)
	xmin = xmin + self.view.pos[1]
	ymin = ymin + self.view.pos[2]
	xmax = xmax + self.view.pos[1]
	ymax = ymax + self.view.pos[2]
	
	gl.glColor3f(.1, .1, .1)
	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)
	gl.glBegin(gl.GL_LINES)
	for x=xticmin,xticmax do
		gl.glVertex2f(x*xstep,ymin)
		gl.glVertex2f(x*xstep,ymax)
	end
	gl.glEnd()
	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	gl.glBegin(gl.GL_LINES)
	for y=yticmin,yticmax do
		gl.glVertex2f(xmin,y*ystep)
		gl.glVertex2f(xmax,y*ystep)
	end
	gl.glEnd()

	gl.glColor3f(.5, .5, .5)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex2f(xmin, 0)
	gl.glVertex2f(xmax, 0)
	gl.glVertex2f(0, ymin)
	gl.glVertex2f(0, ymax)
	gl.glEnd()
end

function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
		
	App.super.update(self)

	local n = #self.size

	self.list = self.list or {}
	glCall(self.list, function()
print'building GL call list...'

		--gl.glColor3f(0,1,1)
		gl.glBegin(gl.GL_LINES)
		local sizeMinusOne = self.size-1
		for indexMinusOne in (sizeMinusOne-1):range() do
			local index = indexMinusOne+1
			for k=1,n do
				if index[k] < sizeMinusOne[k] then
					local nextIndex = matrix(index)
					nextIndex[k] = nextIndex[k] + 1
					glColor((index-1):ediv(sizeMinusOne))
					glVertex(self.Xs[index])
					glColor((nextIndex-1):ediv(sizeMinusOne))
					glVertex(self.Xs[nextIndex])
				end
			end
		end
		gl.glEnd()

--[[
		local colors = matrix{
			{1,0,0},
			{0,1,0},
			{0,0,1},
		}

		local scale = .05
		gl.glMatrixMode(gl.GL_MODELVIEW_MATRIX)
		gl.glPushMatrix()
		gl.glTranslatef(0,0,.1 * scale)
		gl.glBegin(gl.GL_LINES)
		for indexMinusOne in sizeMinusOne:range() do
			local index = indexMinusOne+1
			local u = self.Xs[index]
			local e = self.es[index]:T()
			for k=1,n do
				glColor(colors[k])
				glVertex(u)
				glVertex(u + scale * e[k])
			end
		end
		gl.glEnd()
		gl.glPopMatrix()
--]]	
	end)

	self:drawGrid()
end

App():run()
