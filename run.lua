#! /usr/bin/env luajit
require 'ext'
local bit = require 'bit'
local gl = require 'gl'
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


local App = class(Orbit(View.apply(ImGuiApp)))

App.title = 'reconstruct surface from geodesics' 
App.viewDist = 1

--matrix.__tostring = tolua

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
	else
		error"idk how to invert this"
	end
end


-- used for providing initial values for the metric
-- and verifying accuracy of numerical calculations
local Geometry = class()

function Geometry:init(app)
	self.app = app

	local Tensor = symmath.Tensor
	if self.coordVars then
		Tensor.coords{{variables=self.coordVars}}
		local g = self:createMetric()
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
		xlabel = self.coords[1],
		ylabel = self.coords[2],
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


local PolarHolGeom = class(Geometry)
	
PolarHolGeom.coords = {'r', 'θ'}
PolarHolGeom.coordVars = {symmath.vars('r', 'theta')}

PolarHolGeom.umin = matrix{1, 0}
PolarHolGeom.umax = matrix{3, 2 * math.pi}

function PolarHolGeom:createMetric()
	local r, theta = table.unpack(self.coordVars)
	return symmath.Tensor('_ab', {1, 0}, {0, r^2})
end

function PolarHolGeom:calc_gs()
	local app = self.app
	return app.size:lambda(function(i,j)
		-- manually
		--local r = self.app.xs[i][j][1]
		--return matrix{{1,0},{0,r^2}}
		-- automatically 
		return self.calc.g(self.app.xs[i][j])
	end)
end

function PolarHolGeom:calc_conns()
	return self.app.size:lambda(function(i,j)
		-- manually
		--local r = self.app.xs[i][j][1]
		-- 1/Γ^θ_θr = 1/Γ^θ_rθ = -Γ^r_θθ = r
		--return matrix{ {{0,0},{0,-r}}, {{0,1/r},{1/r,0}} }
		-- automatically
		return self.calc.Conn(self.app.xs[i][j])
	end)
end

-- The thing about non-holonomic geometry is
-- it needs commutation coefficients as well.
-- Otherwise how does it know how far to integrate the geodesics
-- to get to the next coordinate location?
-- This information is typically stored in the metric of the holonomic coordinate map.
local PolarNonHolGeom = class(Geometry)

function PolarNonHolGeom:calc_conns()
	return self.app.size:lambda(function(i,j)
		local r = app.xs[i][j][1]
		-- Γ^θ_rθ = -Γ^r_θθ = 1/r
		return matrix{ {{0,0},{0,-1/r}}, {{0,1/r},{0,0}} }
	end)
end

function App:initGL()
	self.size = matrix{5,16}
		
	self.geom = PolarHolGeom(self)
	--self.geom = PolarNonHolGeom(self)
	
	self.umin = self.geom and self.geom.umin or matrix{-1, -1}
	self.umax = self.geom and self.geom.umax or matrix{1, 1}
	
	local n = #self.size
	
	--[[ cell centered, including borders
	self.dx = matrix{n}:ones():emul(self.umax - self.umin):ediv(self.size)
	self.xs = self.size:lambda(function(i,j) return matrix{i-.5, j-.5}:emul(self.dx) + self.umin end)
	--]]
	-- [[ vertex centered, excluding borders, so position 2,2 is at umin (useful for centering the corner vertex)
	self.dx = matrix{n}:ones():emul(self.umax - self.umin):ediv(self.size-3)
	self.xs = self.size:lambda(function(i,j) return matrix{i-2, j-2}:emul(self.dx) + self.umin end)
	--]]
	
	self.gs = self.geom and self.geom.calc_gs and self.geom:calc_gs() 
		or self.size:lambda(function(i,j) return matrix{n,n}:ident() end)
	self.gUs = self.size:lambda(function(i,j)
		return self.gs[i][j]:inv()
	end)
	self.dgs = self.size:lambda(function(i,j)
		-- dg_abc = dg_bc/dx^a
		local dg = matrix{n}:lambda(function(a)
			-- using a constant border
			local ip = matrix{i,j} ip[a] = math.min(ip[a]+1, self.size[a])
			local im = matrix{i,j} im[a] = math.max(im[a]-1, 1)
			-- using a first-order derivative
			return (self.gs[ip] - self.gs[im]) / (self.xs[ip][a] - self.xs[im][a])
		end)
		-- reorientate dg_abc = g_ab,c
		return matrix{n,n,n}:lambda(function(a,b,c) return dg[c][a][b] end)
	end)
	self.conns = self.size:lambda(function(i,j)
		local dg = self.dgs[i][j]
		local numConnLower = matrix{n,n,n}:lambda(function(a,b,c)
			return .5 * (dg[a][b][c] + dg[a][c][b] - dg[b][c][a])
		end)
		local gU = self.gUs[i][j]
		return matrix{n,n,n}:lambda(function(a,b,c)
			local s = 0
			for d=1,n do
				s = s + gU[a][d] * numConnLower[d][b][c]
			end
			return s
		end)
	end)

	print(self.geom.coords[1]..': '..matrix{self.size[1]}:lambda(function(i) return self.xs[i][1][1] end))
	print(self.geom.coords[2]..': '..matrix{self.size[2]}:lambda(function(i) return self.xs[1][i][2] end))

	if self.geom then
		self.geom:testExact()
	end

	-- embedded space position
	self.Xs = self.size:lambda(function(i,j)
		return matrix{n}:zeros()
	end)
	self.es = self.size:lambda(function(i,j)
		return matrix{n,n}:zeros()
	end)
	
	--local i,j = (self.size/2):map(math.floor):unpack()
	local i,j = 2,2 
	self.es[i][j] = matrix{n,n}:lambda(function(i,j) return i == j and 1 or 0 end)

	-- now to reconstruct the es based on the conns ...
	-- [=[ flood fill 
	local todo = table{ matrix{i,j} }
	local sofar = table()
	while #todo > 0 do
		local index = todo:remove(1)
		sofar:insert(index)
		local i,j = index:unpack()
		local conn = self.conns[i][j]
		-- for each direction ...
		local x = self.xs[i][j]
		local g = self.gs[i][j]	
		local _ = matrix.index
		for k=1,n do
			local connk = conn(_,_,k)
	
			for dir=-1,1,2 do 
				local nextIndex = matrix(index)
				nextIndex[k] = nextIndex[k] + dir
				local ni, nj = nextIndex:unpack()
				-- skip the edges
				if 2 <= ni and ni <= self.size[1]-1
				and 2 <= nj and nj <= self.size[2]-1
				and not table.find(sofar, nextIndex) 
				and not table.find(todo, nextIndex)
				then
					local nextConnK = self.conns[ni][nj](_,_,k)
					
					-- ds = sqrt(g_ab dx^a dx^b) 
					local dx = matrix{n}:zeros()
					dx[k] = dir * self.dx[k]
					local ds = math.sqrt(dx * g * dx)
				
					local eOrig = self.es[i][j]
					local XOrig = self.Xs[i][j]
					
					local e = matrix(eOrig)
					local X = matrix(XOrig)

					--[[ forward-euler
					e = e + (e * connk) * ds
					X = X + e(_,k) * ds
					--]]
					-- [[ rk4 ...
					e = int_rk4(0, e, function(s, e)
						local f = s / ds
						-- treating connections as constant
						--return connk * e
						-- interpolating connections between cells
						-- (notice, when integrating polar coordinates around theta, r remains constant, and so do the connection coefficients)
						local conn = nextConnK * f + connk * (1 - f)
						return e * connk
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

					if ni == 3 then
						print()
						print('index='..index)
						--print('self.dx='..self.dx)
						print('nextIndex='..nextIndex)
						print('dx='..dx)
						print('ds='..ds)
						print('x='..self.xs[i][j])
						print('conn[k]='..connk)
						print('XOrig='..XOrig)
						print('X='..X)
						print('eOrig='..eOrig)
						print('e='..e)
					end

					self.es[ni][nj] = e
					self.Xs[ni][nj] = X

--print('e2 from '..eOrig(_,2)..' to '..e(_,2)..' changing by '..(eOrig(_,2) - e(_,2)):norm())
--print('|e2| from '..eOrig(_,2):norm()..' to '..e(_,2):norm()..' changing by '..(e(_,2):norm() - eOrig(_,2):norm()))
--print('X from '..XOrig..' to '..X..' changing by '..(X - XOrig):norm())
					
					--todo:insert(nextIndex)
					--todo:insert(math.random(#todo+1), nextIndex)
					--todo:insert(math.floor((#todo+1)/2), nextIndex)
					-- the linear dynamic system method only works right for polar coordinates when we use this order
					todo:insert(1, nextIndex)
				end
			end	
		end
	end
	--]]
end

function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	App.super.update(self)

	local n = #self.size

	--gl.glColor3f(0,1,1)
	gl.glBegin(gl.GL_LINES)
	local size = self.size-1
	for i=2,size[1] do
		for j=2,size[2] do
			-- show the parameter space 
			if i < size[1] then
				gl.glColor3f((i-1)/size[1], (j-1)/size[2], .5)
				gl.glVertex2d(self.Xs[i][j]:unpack())
				gl.glColor3f(i/size[1], (j-1)/size[2], .5)
				gl.glVertex2d(self.Xs[i+1][j]:unpack())
			end
			if j < size[2] then
				gl.glColor3f((i-1)/size[1], (j-1)/size[2], .5)
				gl.glVertex2d(self.Xs[i][j]:unpack())
				gl.glColor3f((i-1)/size[1], j/size[2], .5)
				gl.glVertex2d(self.Xs[i][j+1]:unpack())
			end
		end
	end
	gl.glEnd()

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
	for i=2,self.size[1]-1 do
		for j=2,self.size[2]-1 do
			local u = self.Xs[i][j]
			local e = self.es[i][j]:T()
			for k=1,n do
				gl.glColor3f(colors[k]:unpack()) -- color by axis	
				gl.glVertex2d(u:unpack())
				local ek = e[k]
				--[[ normalize or not
				ek = ek / ek:norm()
				--]]
				gl.glVertex2d((u + scale * ek):unpack())
			end
		end
	end
	gl.glEnd()
	gl.glPopMatrix()
end

App():run()
