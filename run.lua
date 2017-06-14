#! /usr/bin/env luajit
require 'ext'
local bit = require 'bit'
local gl = require 'gl'
local matrix = require 'matrix'
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

function App:testExact()
	local exactConns = self.size:lambda(function(i,j)
		-- flat space:
		--return matrix{n,n,n}:zeros()
		--[[ anholonomic polar in coordinates r, theta:
		local r = self.xs[i][j][1]
		-- Gamma^theta_r_theta = -Gamma^r_theta_theta = 1/r
		return matrix{ {{0,0},{0,-1/r}}, {{0,1/r},{0,0}} }
		--]]
		-- [[ holonomic polar
		local r = self.xs[i][j][1]
		-- 1/Gamma^theta_theta_r = 1/Gamma^theta_r_theta = -Gamma^r_theta_theta = r
		return matrix{
			{
				{0,0},
				{0,-r},
			},
			{
				{0,1/r},
				{1/r,0},
			}
		}
		--]]
	end)
	local connNumConnDiff = self.size:lambda(function(i,j)
		return (self.conns[i][j] - exactConns[i][j]):norm()
	end)

	local gnuplot = require 'gnuplot'
	local x = matrix{self.size[1]-2}:lambda(function(i) return self.xs[i+1][1][1] end)
	local y = matrix{self.size[2]-2}:lambda(function(j) return self.xs[1][j+1][2] end)
	local z = (self.size-2):lambda(function(i,j) return connNumConnDiff[i+1][j+1] end)
	gnuplot{
		output = 'conn numeric vs analytic.png',
		style = 'data lines',
		xlabel = self.coordNames[1],
		ylabel = self.coordNames[2],
		--log = 'z',
		griddata = {x = x, y = y, z},
		{splot=true, using='1:2:3', title = ' Γ^a_b_c |analytic - numeric|'},
	}
	print('max error between analytical and numerical', z:normLInf())
end

function App:initGL()
	self.coordNames = {'r', 'θ'}
	self.size = matrix{16,16}
	-- [[ polar 
	self.umin = matrix{.1, 0}
	self.umax = matrix{1, 2 * math.pi}
	--]]
	--[[
	self.umin = matrix{-1, -1}
	self.umax = matrix{1, 1}
	--]]
	local n = #self.size
	self.dx = matrix{n}:ones():emul(self.umax - self.umin):ediv(self.size)
	self.xs = self.size:lambda(function(i,j)
		return matrix{i-.5, j-.5}:emul(self.dx) + self.umin
	end)
	self.gs = self.size:lambda(function(i,j)
		-- holonomic polar
		local r = self.xs[i][j][1]
		return matrix{{1,0},{0,r^2}}
	end)
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
	
	self:testExact()
	
	-- embedded space position
	self.Xs = self.size:lambda(function(i,j)
		return matrix{n}:zeros()
	end)
	self.es = self.size:lambda(function(i,j)
		return matrix{n,n}:zeros()
	end)
	
	local i,j = (self.size/2):map(math.floor):unpack()
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
		--for k=1,n do
		for k=2,2 do
			local connk = conn(_,_,k)
	
			for dir=-1,1,2 do 
				local nextIndex = matrix(index)
				nextIndex[k] = nextIndex[k] + dir
				local ni, nj = nextIndex:unpack()
				if 1 <= ni and ni <= self.size[1]
				and 1 <= nj and nj <= self.size[2]
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
					e = e + (connk * e) * ds
					X = X + e(_,k) * ds
					--]]
					-- [[ rk4 ...
					e = int_rk4(0, e, function(s, e)
						local f = s / ds
						-- treating connections as constant
						return connk * e
						-- interpolating connections between cells
						--local conn = nextConnK * f + connk * (1 - f)
						--return conn * e
					end, ds)
					
					X = int_rk4(0, X, function(s, X)
						local f = s / ds
						-- constant:
						--return eOrig(_,k)
						-- linear interpolation
						return e(_,k) * f + eOrig(_,k) * (1 - f)
					end, ds)
					--]]
					-- if d/ds e = conn e then ...
					-- then we can solve this as a linear dynamic system!
					-- even though this assumes conn is constant (which it's not)
					--[[
					e = int_rk4(s, e, function(s, e)
						return connk * e
					end, ds)
					--]]
					
					--[[ normalize columns
					e = e:T()
					for k=1,n do
						e[k] = e[k] / e[k]:norm()
					end
					e = e:T()
					--]]
					
					self.es[ni][nj] = e
					self.Xs[ni][nj] = X

print('e2 from '..eOrig(_,2)..' to '..e(_,2)..' changing by '..(eOrig(_,2) - e(_,2)):norm())
print('|e2| from '..eOrig(_,2):norm()..' to '..e(_,2):norm()..' changing by '..(e(_,2):norm() - eOrig(_,2):norm()))
print('X from '..XOrig..' to '..X..' changing by '..(X - XOrig):norm())
					
					todo:insert(nextIndex)
					--todo:insert(math.random(#todo+1), nextIndex)
					--todo:insert(math.floor((#todo+1)/2), nextIndex)
					--todo:insert(1, nextIndex)
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

	gl.glColor3f(0,1,1)
	gl.glBegin(gl.GL_LINES)
	for i=1,self.size[1] do
		for j=1,self.size[2] do
			-- show the parameter space 
			if i < self.size[1] then
				gl.glVertex2d(self.Xs[i][j]:unpack())
				gl.glVertex2d(self.Xs[i+1][j]:unpack())
			end
			if j < self.size[2] then
				gl.glVertex2d(self.Xs[i][j]:unpack())
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
	for i=1,self.size[1] do
		for j=1,self.size[2] do
			local u = self.Xs[i][j]
			local e = self.es[i][j]:T()
			for k=1,n do
				gl.glColor3f(colors[k]:unpack())
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
