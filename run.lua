#! /usr/bin/env luajit
require 'ext'
local bit = require 'bit'
local gl = require 'gl'
local matrix = require 'matrix'
local ImGuiApp = require 'imguiapp'

local View = require 'glapp.view'
local Orbit = require 'glapp.orbit'

local App = class(Orbit(View.apply(ImGuiApp)))

App.title = 'reconstruct surface from geodesics' 
App.viewDist = 10

function App:initGL()
	self.size = matrix{16,16}
	self.umin = matrix{-10, -10}
	self.umax = matrix{10, 10}
	local n = #self.size
	self.du = matrix{n}:ones():emul(self.umax - self.umin):ediv(self.size)
	self.us = self.size:lambda(function(i,j)
		return matrix{i-.5, j-.5}:emul(self.du) + self.umin
	end)
	self.conns = self.size:lambda(function(i,j)
		-- flat space:
		--return matrix{n,n,n}:zeros()
		-- anholonomic cylinder in coordinates r, theta:
		local r = self.us[i][j][1]
		-- Gamma^theta_r_theta = -Gamma^r_theta_theta = 1/r
		return matrix{ {{0,0},{0,-1/r}}, {{0,1/r},{0,0}} }
	end)
	self.es = self.size:lambda(function(i,j)
		return matrix{n,n}:zeros()
	end)
	local i,j = (self.size/2):map(math.floor):unpack()
	self.es[i][j] = matrix{n,n}:lambda(function(i,j) return i == j and 1 or 0 end)
	-- now to reconstruct the es based on the conns ...
	local todo = table{ matrix{i,j} }
	local sofar = table()
	while #todo > 0 do
		local index = todo:remove(1)
		local i,j = index:unpack()
		local conn = self.conns[i][j]
		-- for each direction ...
		local _ = matrix.index
		local r = self.us[i][j][1]
		print('r = '..r)
		print('1/r = '..(1/r))
		print()
		print'Γ:'
		print(conn)
		print('Γ^θ_rθ = '..conn[2][1][2])
		print('Γ^r_θθ = '..conn[1][2][2])
		local coordNames = {'r', 'θ'}
		for k=1,n do
			local connk = conn(_,_,k)
			print()
			print('Γ^a_b'..coordNames[k]..':')
			print(connk)
	
			for dir=-1,1,2 do 
				local nextIndex = matrix(index)
				nextIndex[k] = nextIndex[k] + dir
				local ni, nj = nextIndex:unpack()
				if 1 <= ni and ni <= self.size[1]
				and 1 <= nj and nj <= self.size[2]
				and not table.find(sofar, nextIndex) 
				then
					self.es[ni][nj] = self.es[i][j]
						+ (connk * self.es[i][j]:T()):T() * self.du[k]
print()
print(nextIndex..':')
print(self.es[ni][nj])
					todo:insert(nextIndex)
				end
			end	
		end
	end
end

function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	App.super.update(self)

-- [[ connect cell centers
	gl.glColor3f(0,1,1)
	gl.glBegin(gl.GL_LINES)
	for i=1,self.size[1] do
		for j=1,self.size[2] do
			-- show the parameter space 
			if i < self.size[1] then
				gl.glVertex2d(self.us[i][j]:unpack())
				gl.glVertex2d(self.us[i+1][j]:unpack())
			end
			if j < self.size[2] then
				gl.glVertex2d(self.us[i][j]:unpack())
				gl.glVertex2d(self.us[i][j+1]:unpack())
			end
		end
	end
	gl.glEnd()
--]]
--[[ connect cell vertices
	local corners = matrix{
		{-.5, -.5},
		{-.5, .5},
		{.5, .5},
		{.5, -.5},
	}

	gl.glColor3f(0,1,1)
	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
	gl.glBegin(gl.GL_QUADS)
	for i=1,self.size[1]-1 do
		for j=1,self.size[2]-1 do
			for _,corner in ipairs(corners) do
				gl.glVertex2d((self.us[i][j] + self.du:emul(corner)):unpack())
			end
		end
	end
	gl.glEnd()
	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
--]]

	gl.glMatrixMode(gl.GL_MODELVIEW_MATRIX)
	gl.glPushMatrix()
	gl.glTranslatef(0,0,.1)
	gl.glColor3f(1,1,0)
	gl.glBegin(gl.GL_LINES)
	for i=1,self.size[1] do
		for j=1,self.size[2] do
			local u = self.us[i][j]
			local e = self.es[i][j]
			gl.glVertex2d(u:unpack())
			gl.glVertex2d((u + e[1]):unpack())
			gl.glVertex2d(u:unpack())
			gl.glVertex2d((u + e[2]):unpack())
		end
	end
	gl.glEnd()
	gl.glPopMatrix()
end

App():run()
