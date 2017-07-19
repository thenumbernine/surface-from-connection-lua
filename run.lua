#! /usr/bin/env luajit
require 'ext'
local bit = require 'bit'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local gl = require 'gl'
local glCall = require 'gl.call'
local GLProgram = require 'gl.program'
local matrix = require 'matrix'
local complex = require 'symmath.complex'
local gnuplot = require 'gnuplot'
local ImGuiApp = require 'imguiapp'
local View = require 'glapp.view'
local Orbit = require 'glapp.orbit'
local symmath = require 'symmath'
local template = require 'template'
local clnumber = require 'cl.obj.number'
local CLEnv = require 'cl.obj.env'

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

function matrix:isfinite()
	for i in self:iter() do
		if not math.isfinite(self[i]) then return false end
	end
	return true
end

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


--local surfaceBuildOrder = 'first'
local surfaceBuildOrder = 'last'
--local surfaceBuildOrder = 'random'
--local surfaceBuildOrder = 'middle'


-- used for providing initial values for the metric
-- and verifying accuracy of numerical calculations
local Geometry = class()

function Geometry:init(app)
	self.app = app
	self.xmin = matrix(self.xmin)
	self.xmax = matrix(self.xmax)
	self.startCoord = matrix(self.startCoord)

	self.coordVars = table.map(self.coords, function(name)
		return symmath.var(name)
	end)

	if self.createMetric then
		symmath.Tensor.coords{{variables=self.coordVars}}
		local n = #self.coords
		local g = self:createMetric()
		local gU = symmath.Tensor('^ab', table.unpack( (symmath.Matrix.inverse(g)) ))
		local dg = symmath.Tensor('_abc', table.unpack(g'_ab,c'()))
		local ConnL = symmath.Tensor('_abc', table.unpack( ((dg'_abc' + dg'_acb' - dg'_bca')/2)() ))
		local Conn = symmath.Tensor('^a_bc', table.unpack( (gU'^ad' * ConnL'_dbc')() ))	
		self.calc = {
			g = self:compileTensor(g),
			Conn = self:compileTensor(Conn),
		}
	end
end

function Geometry:testExact()
	local app = self.app
	local exactConns = self:calcFromEqns_conns()
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

function Geometry:calcFromEqns_gs()
	if not self.calc.g then return end
	local app = self.app
	return app.size:lambda(function(...)
		return self.calc.g(self.app.xs[{...}])
	end)
end

-- this is only used for comparing the error between the numerical and the exact connections
function Geometry:calcFromEqns_conns()
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
	* create_conns.  otherwise these return identity g_ab's and zero Conn^a_bc's.
--]]


-- 2D geometries


local Polar = class(Geometry)
Polar.coords = {'r', 'θ'}
-- [[ initialize our basis of e=I at r=1 ...
Polar.xmin = matrix{1, 0}
Polar.xmax = matrix{10, 2 * math.pi}
Polar.startCoord = {1,0}
--]]
--[[ applying a rotation, e=I at r=1 and theta!=0 still produces the same shape.
Polar.xmin = matrix{1, 0}
Polar.xmax = matrix{10, 2 * math.pi}
Polar.startCoord = {1, math.pi/2}
--]]
--[[ ...otherwise the shape gets messed up -- for r=2
Polar.xmin = matrix{1, 0}
Polar.xmax = matrix{10, 2 * math.pi}
Polar.startCoord = {2,0}
--]]
--[[ ...otherwise the shape gets messed up -- for r=1/2 (needs the range readjusted so rmin isn't 1)
Polar.xmin = matrix{.1, 0}
Polar.xmax = matrix{2, 2 * math.pi}
Polar.startCoord = {.5,0}
--]] 
-- [[ analytically
function Polar:createMetric()
	local r, theta = self.coordVars:unpack()
	return symmath.Tensor('_ab', {1, 0}, {0, r^2})
end
--]]
--[[ purely numerically
function Polar:create_conns()
	return self.app.size:lambda(function(...)
		local x = self.app.xs[{...}]
		local r, theta = x:unpack()
		-- 1/Γ^θ_θr = 1/Γ^θ_rθ = -Γ^r_θθ = r
		return matrix{
			{
				{0, 0},
				{0, -r},
			},
			{
				{0, 1/r},
				{1/r, 0},
			},
		}
	end)
end
--]]


-- The thing about non-holonomic geometry is
-- it needs commutation coefficients as well.
-- Otherwise how does it know how far to integrate the geodesics
-- to get to the next coordinate location?
-- This information is typically stored in the metric of the holonomic coordinate map.
local PolarAnholonomic = class(Geometry)
PolarAnholonomic.coords = {'r', 'θ'}
PolarAnholonomic.xmin = matrix{1, 0}
PolarAnholonomic.xmax = matrix{10, 2 * math.pi}
PolarAnholonomic.startCoord = {1,0}
function PolarAnholonomic:create_conns()
	return self.app.size:lambda(function(i,j)
		local r = self.app.xs[i][j][1]
		-- Γ^θ_rθ = -Γ^r_θθ = 1/r
		return matrix{ {{0,0},{0,-1/r}}, {{0,1/r},{0,0}} }
	end)
end
-- if you can get the lengths holonomic basis' vectors then reconstructing the original shape is much easier 
function PolarAnholonomic:get_basis_lengths(r, theta)
	return matrix{1, r}
end


-- cyl surface can't be reconstructed
-- because it needs extrinsic curvature information
-- it has a connection of zero

-- sphere surface likewise is a 2 dimensional system inside 3 dimensions
-- like cyl surface, it needs extrinsic curvature information to be properly rebuilt 
local SphereSurface = class(Geometry)
local eps = .01
SphereSurface.coords = {'θ', 'φ'}
SphereSurface.xmin = matrix{eps, eps}
SphereSurface.xmax = matrix{math.pi-eps, 2*math.pi-eps}

SphereSurface.startCoord = {math.pi/2, math.pi}	-- poles along x axis
--SphereSurface.startCoord = {2*eps, math.pi}			-- stretched to infinite becuase of infinite connections
--SphereSurface.startCoord = {math.pi/4, math.pi}
--SphereSurface.startCoord = {math.pi/2, 2*eps}		-- mostly y>0
--SphereSurface.startCoord = {math.pi/2, math.pi*2-2*eps}	-- mostly y<0
function SphereSurface:createMetric()
	local theta, phi = self.coordVars:unpack()
	local r = 1
	return symmath.Tensor('_ab', {r^2, 0}, {0, r^2 * symmath.sin(theta)^2})
end


local TorusSurface = class(Geometry)
TorusSurface.coords = {'θ', 'φ'}
TorusSurface.xmin = {-math.pi, -math.pi}
TorusSurface.xmax = {math.pi, math.pi}
TorusSurface.startCoord = {0, 0}
function TorusSurface:createMetric()
	local theta, phi = self.coordVars:unpack()
	local r = 2
	local R = 5
	return symmath.Tensor('_ab', {r^2, 0}, {0, (R + r * symmath.sin(theta))^2})
end

local PoincareDisk2D = class(Geometry)
PoincareDisk2D.coords = {'u', 'v'}
PoincareDisk2D.xmin = {-1, -1}
PoincareDisk2D.xmax = {1, 1}
PoincareDisk2D.startCoord = {0, 0}
function PoincareDisk2D:createMetric()
	local u, v = self.coordVars:unpack()
	return symmath.Tensor('_ab', {4 / (1 - u^2 - v^2), 0}, {0, 4 / (1 - u^2 - v^2)})
end


-- hmm, how to incorporate signature into the metric ...
local Minkowski2D = class(Geometry)
Minkowski2D.coords = {'t', 'x'}
Minkowski2D.xmin = {-1, -1}
Minkowski2D.xmax = {1, 1}
Minkowski2D.startCoord = {0,0}
function Minkowski2D:createMetric()
	return symmath.Tensor('_ab', {-1, 0}, {0, 1})
end


-- here's Schwarzschild in time and in radial 
-- it is treating Rs as constant, which means this metric is true for the spacetime *outside* of the massive body
local Schwarzschild1Plus1 = class(Geometry)
Schwarzschild1Plus1.coords = {'r', 't'}
Schwarzschild1Plus1.xmin = {-2, -2}
Schwarzschild1Plus1.xmax = {2, 2}
Schwarzschild1Plus1.startCoord = {2, 0}
function Schwarzschild1Plus1:createMetric()
	local r, t = self.coordVars:unpack()
	local mass = .1
	return symmath.Tensor('_ab', {1/(1 - 2 * mass / r), 0}, {0, 1 - 2 * mass / r})
end

-- Schwarzschild in time and radial
-- backwards from traditional order: r, t (so that time can point upwards)
-- except for within and outside of the matter source
-- this is my first metric that is based on numerically specifying g_ab, then computing Conn^a_bc
-- since this doesn't have an extra dimension to anchor it, as the spacetime grows from r+ to r- it gets really twisted
local Schwarzschild1Plus1EOS = class(Geometry)
Schwarzschild1Plus1EOS.coords = {'r', 't'}
Schwarzschild1Plus1EOS.xmin = {-2, -2}
Schwarzschild1Plus1EOS.xmax = {2, 2}
Schwarzschild1Plus1EOS.startCoord = {2, 0}
function Schwarzschild1Plus1EOS:create_gs()
	return self.app.size:lambda(function(...)
		local r, t = self.app.xs[{...}]:unpack()
		local mass = .1
		local radius = 1
		local radialFrac = math.max(math.abs(r), radius) / radius
		local massWithinRadius = radialFrac * mass
		return matrix{
			{1/(1 - 2 * massWithinRadius / r), 0},
			{0, 1 - 2 * massWithinRadius / r},
		}
	end)
end

-- here's one from section 2.5 of the "Covariant Loop Quantum Gravity" book
local LagrangianTotalEnergy = class(Geometry)
LagrangianTotalEnergy.coords = {'a', 'b'}
LagrangianTotalEnergy.xmin = {-2, -2}
LagrangianTotalEnergy.xmax = {2, 2}
LagrangianTotalEnergy.startCoord = {0, 0}
function LagrangianTotalEnergy:create_gs()
	return self.app.size:lambda(function(...)
		local a, b = self.app.xs[{...}]:unpack()
		local E = 10
		local r = 2 * E - a^2 - b^2
		return matrix{{r, 0}, {0, r}}
	end)
end


-- 3D geometries


local Cylinder = class(Geometry)
Cylinder.coords = {'r', 'θ', 'z'}
Cylinder.xmin = {1, 0, -5}
Cylinder.xmax = {10, 2*math.pi, 5}
Cylinder.startCoord = {1,math.pi,0}
function Cylinder:createMetric()
	local r, theta, z = self.coordVars:unpack()
	return symmath.Tensor('_ab', {1, 0, 0}, {0, r^2, 0}, {0, 0, 1})
end


local Sphere = class(Geometry)
local eps = .05	-- the holonomic connection gets singularities (cot theta = inf) at the boundaries
				-- this could be avoided if the metric was evaluated at grid centers instead of vertices.
Sphere.coords = {'r', 'θ', 'φ'}
Sphere.xmin = {1, eps, -math.pi + eps}
Sphere.xmax = {10, math.pi - eps, math.pi - eps}
Sphere.startCoord = {1, math.pi/2, 0}
--Sphere.startCoord = {2, math.pi/2, 0}	-- squashed to an ellipsoid, just like the polar case
--Sphere.startCoord = {1, math.pi/2, math.pi-2*eps}	-- changing phi_0 doesn't affect it at all though
--Sphere.startCoord = {2, math.pi/2, math.pi-2*eps}	-- ... though it does a tiny bit (makes some waves in the coordinate system) if r_0 is not 1 
function Sphere:createMetric()
	local r, theta, phi = self.coordVars:unpack()
	return symmath.Tensor('_ab', {1,0,0}, {0, r^2, 0}, {0, 0, r^2 * symmath.sin(theta)^2})
end


local Torus = class(Geometry)
Torus.coords = {'r', 'θ', 'φ'}
Torus.xmin = {1, -math.pi, -math.pi}
Torus.xmax = {2, math.pi, math.pi}
Torus.startCoord = {1, -math.pi, -math.pi}
function Torus:createMetric()
	local r, theta, phi = self.coordVars:unpack()
	local R = 5
	return symmath.Tensor('_ab', {1, 0, 0}, {0, r^2, 0}, {0, 0, (R + r * symmath.sin(theta))^2})	-- does sin(theta) work as well?
end


local PoincareDisk3D = class(Geometry)
PoincareDisk3D.coords = {'u', 'v', 'w'}
PoincareDisk3D.xmin = {-1, -1, -1}
PoincareDisk3D.xmax = {1, 1, 1}
PoincareDisk3D.startCoord = {0, 0, 0}
function PoincareDisk3D:createMetric()
	local u, v, w = self.coordVars:unpack()
	return symmath.Tensor('_ab', {4 / (1 - u^2 - v^2 - w^2), 0, 0}, {0, 4 / (1 - u^2 - v^2 - w^2), 0}, {0, 0, 4 / (1 - u^2 - v^2 - w^2)})
end


-- Schwarzschild in time and radial
-- backwards from traditional order: r, t (so that time can point upwards)
-- except for within and outside of the matter source
-- this is my first metric that is based on numerically specifying g_ab, then computing Conn^a_bc
-- since this doesn't have an extra dimension to anchor it, as the spacetime grows from r+ to r- it gets really twisted
local Schwarzschild2Plus1EOS = class(Geometry)
Schwarzschild2Plus1EOS.coords = {'t', 'x', 'y'}
Schwarzschild2Plus1EOS.xmin = {-2, -2, -2}
Schwarzschild2Plus1EOS.xmax = {2, 2, 2}
--Schwarzschild2Plus1EOS.startCoord = {0, 2, 2}
Schwarzschild2Plus1EOS.startCoord = {0, 0, 0}
function Schwarzschild2Plus1EOS:create_gs()
	return self.app.size:lambda(function(...)
		local t, x, y = self.app.xs[{...}]:unpack()
		local mass = .05
		local radius = 1
		local r = math.sqrt(x*x + y*y)
		local radialFrac = math.max(math.abs(r), radius) / radius
		local massWithinRadius = radialFrac * mass
		return matrix{
			{1 - 2 * massWithinRadius / r, 0, 0},
			{0, 1/(1 - 2 * massWithinRadius / r), 0},
			{0, 0, 1/(1 - 2 * massWithinRadius / r)},
		}
	end)
end


local SchwarzschildSphere2Plus1EOS = class(Geometry)
SchwarzschildSphere2Plus1EOS.coords = {'t', 'r', 'φ'}
SchwarzschildSphere2Plus1EOS.xmin = {-2, .1, -math.pi}
SchwarzschildSphere2Plus1EOS.xmax = {2, 2, math.pi}
SchwarzschildSphere2Plus1EOS.startCoord = {0, 1, 0}
function SchwarzschildSphere2Plus1EOS:create_gs()
	return self.app.size:lambda(function(...)
		local t, r, phi = self.app.xs[{...}]:unpack()
		local mass = .05
		local radius = 1
		local radialFrac = math.max(r, radius) / radius
		local massWithinRadius = radialFrac * mass
		return matrix{
			{1 - 2 * massWithinRadius / r, 0, 0},
			{0, 1/(1 - 2 * massWithinRadius / r), 0},
			{0, 0, r*r},
		}
	end)
end


-- here's a connection coefficient that gives rise to the stress-energy of a uniform electric field 
-- it's based on an analytical connection, but I haven't made a metric for it just yet
local UniformElectricFieldNumeric = class(Geometry)
UniformElectricFieldNumeric.coords = {'t', 'x', 'y', 'z'}		-- ut oh, now we introduce metric signatures ... 
UniformElectricFieldNumeric.xmin = {-1, -1, -1}
UniformElectricFieldNumeric.xmax = {1, 1, 1}
UniformElectricFieldNumeric.startCoord = {0, 0, 0}
function UniformElectricFieldNumeric:create_conns()
	local E = 1
	return self.app.size:lambda(function(...)
		return matrix{
			{
				{0,E,0,0},
				{E,0,0,0},
				{0,0,0,0},
				{0,0,0,0},
			},
			{
				{-E,0,0,0},
				{0,0,0,0},
				{0,0,E,0},
				{0,0,0,E},
			},
			{
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0},
			},
		}
	end)
end

-- Here's connections that give rise to the stress-energy for the magnetic field
--  around an infinite wire with no net charge.
-- Still no analytical metric for this either.
-- I got it from my file at "symmath/tests/electrovacuum/infinite wire no charge.lua"
-- notice it is 4D, so I need to start thinking of how to handle that.
local InfiniteWireMagneticFieldNumeric = class(Geometry)
InfiniteWireMagneticFieldNumeric.coords = {'t', 'r', 'φ', 'z'}
InfiniteWireMagneticFieldNumeric.xmin = {-1, -1, -math.pi, -1}
InfiniteWireMagneticFieldNumeric.xmax = {1, 1, math.pi, 1}
InfiniteWireMagneticFieldNumeric.startCoord = {0, 0, 0, 0}
function InfiniteWireMagneticFieldNumeric:create_conns()
	local I = 1
	return self.app.size:lambda(function(...)
		local t, r, phi, z = self.app.xs[{...}]:unpack()
		return matrix{
			{
				{0,0,2*I,0},
				{0,0,0,0},
				{2*I,0,0,0},
				{0,0,0,0},
			},
			{
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0},
			},
			{
				{2*I/r^2,0,0,0},
				{0,2*I/r^2,0,0},
				{0,0,0,0},
				{0,0,0,2*I/r^2},
			},
			{
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0},
			},
		}
	end)
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
	App.super.initGL(self)
	gl.glEnable(gl.GL_DEPTH_TEST)

	self.controlsOpened = ffi.new('bool[1]', true)
	self.geomID = ffi.new('int[1]', 0)

	self:buildSurface'Polar'
end

local geomClassesForName = table{
	-- 2D
	{Polar = Polar},
	{PolarAnholonomic = PolarAnholonomic},
	{SphereSurface = SphereSurface},
	{TorusSurface = TorusSurface},
	{PoincareDisk2D = PoincareDisk2D},
	{Minkowski2D = Minkowski2D},
	{Schwarzschild1Plus1 = Schwarzschild1Plus1},
	{Schwarzschild1Plus1EOS = Schwarzschild1Plus1EOS},
	{LagrangianTotalEnergy = LagrangianTotalEnergy},
	-- 3D
	{Cylinder = Cylinder},
	{Sphere = Sphere},
	{Torus = Torus},
	{PoincareDisk3D = PoincareDisk3D},
	{UniformElectricFieldNumeric = UniformElectricFieldNumeric},
	{Schwarzschild2Plus1EOS = Schwarzschild2Plus1EOS},
	{SchwarzschildSphere2Plus1EOS = SchwarzschildSphere2Plus1EOS},
	-- 4D
	{InfiniteWireMagneticFieldNumeric = InfiniteWireMagneticFieldNumeric},
}
local geomClassNames = geomClassesForName:map(function(kv) return (next(kv)) end)

function App:buildSurface(geomName)
	assert(geomName)	
	local loc, geomClass = geomClassesForName:find(nil, function(kv)
		return next(kv) == geomName
	end)
	assert(geomClass, "couldn't find geometry named "..geomName)
	geomClass = assert(select(2, next(geomClass)))
	self.geomID[0] = loc - 1

	self.env = CLEnv()

	self.animShader = GLProgram{
		vertexCode = [[
varying vec4 color;
uniform float t;
void main() {
	color = gl_Color;
	
	vec3 pos1 = gl_Vertex.xyz; 
	vec3 pos2 = gl_MultiTexCoord0.xyz;
	vec3 pos = mix(pos1, pos2, t);
	gl_Position = gl_ModelViewProjectionMatrix * vec4(pos, 1.);
}
]],
		fragmentCode = [[
varying vec4 color;
void main() {
	gl_FragColor = color;
}
]],
		uniforms = {t = 0},
	}

	self.geom = geomClass(self)



	local n = #self.geom.coords
	self.size = matrix{n}:lambda(I( ({[2]=64, [3]=16, [4]=8})[n] ))

	self.domain = self.env:domain{size=self.size}

	self.xmin = self.geom and self.geom.xmin or matrix{n}:lambda(I(-1))
	self.xmax = self.geom and self.geom.xmax or matrix{n}:lambda(I(1))
	
	self.view.ortho = n == 2

--[=[ opencl code
	-- regenerate these to prevent ffi cdef redefs	
	self.rank1Type = 'rank1_t'
	self.rank2Type = 'rank2_t'
	local rank1TypeCode = 'typedef real '..self.rank1Type..'['..n..'];'
	local rank2TypeCode = 'typedef real '..self.rank2Type..'['..n..'];'
	ffi.cdef(rank1TypeCode)
--]=]

	--[[ cell centered, including borders
	self.dx = matrix{n}:ones():emul(self.xmax - self.xmin):ediv(self.size)
	self.xs = self.size:lambda(function(...) return (matrix{...} - .5):emul(self.dx) + self.xmin end)
	--]]
	-- [[ vertex centered, excluding borders, so position 2,2 is at xmin (useful for centering the corner vertex)
	self.dx = matrix{n}:ones():emul(self.xmax - self.xmin):ediv(self.size-3)
	self.xs = self.size:lambda(function(...) return (matrix{...} - 2):emul(self.dx) + self.xmin end)

--[=[ opencl code
	self.headerCode = table{
		rank1TypeCode,
		template([[
<? for j=0,n-1 do ?>
#define size_<?=j?> <?=size[j+1]?>
#define xmin_<?=j?> <?=clnumber(xmin[j+1])?>
#define xmax_<?=j?> <?=clnumber(xmax[j+1])?>
#define dx_<?=j?> ((xmax_<?=j?> - xmin_<?=j?>) / (real)(size_<?=j?> - 3))
<? end ?>
]], 	{
			n = n,
			clnumber = clnumber,
			size = self.size,
			xmin = self.geom.xmin,
			xmax = self.geom.xmax,
		}),
	}:concat'\n'


	self.xsBuf = self.domain:buffer{name='xs', type=self.rank1Type}
	self.domain:kernel{
		header = self.headerCode,
		argsOut = {self.xsBuf},
		body = template([[
<? for j=0,n-1 do 
?>	xs[index][<?=j?>] = (i.s<?=j?> - 2) * dx_<?=j?> + xmin_<?=j?>;
<? end 
?>]], {n = n}),
	}()
--]=]	
	--]]

--[=[ opencl
	if self.geom.createMetric 
	or self.geom.create_gs
	then
		self.gsBuf = self.domain:buffer{name='gs', type=self.rank2Type}
		self.gUsBuf = self.domain:buffer{name='gUs', type=self.rank2Type}
	end
--]=]

	self:rebuildSurface()
end

function App:rebuildSurface()

	local n = #self.size

	if self.list and self.list.id then
		gl.glDeleteLists(self.list.id, 1)
	end
	self.list = {}

	-- if we are calculating the connection from discrete derivatives of the metric ...
	if self.geom.createMetric 
	or self.geom.create_gs
	then
		local gs = self.geom.create_gs and self.geom:create_gs() or self.geom:calcFromEqns_gs() 
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
		if self.geom and self.geom.create_conns then
			self.geom:testExact()
		end
	
	else
		-- if calcFromEqns_gs isn't there, then rely on create_conns for providing the initial connections
		self.conns = self.geom:create_conns()
	end

	-- embedded space position
	self.Xs = self.size:lambda(function(...)
		return matrix{n}:zeros()
	end)
	self.es = self.size:lambda(function(...)
		return matrix{n,n}:zeros()
	end)

	local xInit = matrix(assert(self.geom.startCoord))
	local i = ((xInit - self.geom.xmin):ediv(self.geom.xmax - self.geom.xmin):emul(self.size-3) + 2.5):map(math.floor)
	for j=1,n do
		i[j] = math.clamp(i[j], 2, self.size[j]-1)
		-- TODO reset the start coord in the GUI if it has to be clamped
	end
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
			
					-- cheating to get around anholonomic coordinate systems
					-- technically I should be using commutation coefficients or something
					local len = 1
					if self.geom.get_basis_lengths then
						local lens = self.geom:get_basis_lengths(self.xs[index]:unpack())
						len = len * lens[k]
					end
		
					-- derivative along coordinate
					local ds = dir * self.dx[k] * len
				
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

					if surfaceBuildOrder == 'last' then
						todo:insert(nextIndex)
					elseif surfaceBuildOrder == 'random' then
						todo:insert(math.random(#todo+1), nextIndex)
					elseif surfaceBuildOrder == 'middle' then
						todo:insert(math.floor((#todo+1)/2), nextIndex)
					elseif surfaceBuildOrder == 'first' then
						-- the linear dynamic system method only works right for polar coordinates when we use this order
						todo:insert(1, nextIndex)
					end
					sofar[tostring(nextIndex)] = true
				end
			end	
		end
	end
	--]]

	-- recenter ...
	local com = matrix{n}:lambda(I(0))
	for i in self.size:range() do
		if self.Xs[i]:isfinite() then
			com = com + self.Xs[i]
		end
	end
	com = com / self.size:prod()
	for i in self.size:range() do
		self.Xs[i] = self.Xs[i] - com
	end
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


local function glColor(m)
	if #m == 2 then
		gl.glColor3d(m[1], m[2], .5)
	elseif #m == 3 then
		gl.glColor3d(m:unpack())
	elseif #m == 4 then
		gl.glColor3d(table.unpack(m,1,3))
	else
		error"can't color this many dimensions"
	end
end

local function glTexCoord(m)
	assert(({
		[2] = gl.glTexCoord2d,
		[3] = gl.glTexCoord3d,
		[4] = gl.glTexCoord4d,
	})[#m])(m:unpack())
end

local function glVertex(m)
	assert(({
		[2] = gl.glVertex2d,
		[3] = gl.glVertex3d,
		[4] = gl.glVertex4d,
	})[#m])(m:unpack())
end

local animating = false
local animTime = ffi.new('float[1]', 0)
local lastTime = 0
function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

	local n = #self.size

	local thisTime = os.clock()	
	if animating then
		local deltaTime = (thisTime - lastTime) / math.pi
		animTime[0] = (((animTime[0] + deltaTime) + 1) % 2) - 1
	end
	lastTime = thisTime
	
	self.animShader:use()

	gl.glUniform1f(self.animShader.uniforms.t.loc, .5 - .5 * math.cos(math.pi * animTime[0]))

	glCall(self.list, function()
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
					glTexCoord(self.xs[index])
					glVertex(self.Xs[index])
					glColor((nextIndex-1):ediv(sizeMinusOne))
					glTexCoord(self.xs[nextIndex])
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

	GLProgram:useNone()

	if n == 2 then
		self:drawGrid()
	end
		
	App.super.update(self)
end

local function hoverTooltip(name)
	if ig.igIsItemHovered() then
		ig.igBeginTooltip()
		ig.igText(name)
		ig.igEndTooltip()
	end
end

local function wrapTooltip(fn)
	return function(name, ...)
		ig.igPushIdStr(name)
		local result = ig[fn]('', ...)
		hoverTooltip(name)
		ig.igPopId()
		return result
	end
end

local sliderTooltip = wrapTooltip'igSliderFloat'
local comboTooltip = wrapTooltip'igCombo'
local inputFloatTooltip = wrapTooltip'igInputFloat'

local float = ffi.new('float[1]', 0)
function App:updateGUI()
	if ig.igBegin('Controls', self.controlsOpened) then
		if ig.igButton(animating and 'Stop Animation' or 'Start Animation') then
			animating = not animating
		end
		
		sliderTooltip('animation coefficient', animTime, -1, 1)

		if comboTooltip('coordinate system', self.geomID, geomClassNames) then
			self:buildSurface(geomClassNames[self.geomID[0]+1])
		end

		local n = #self.size
		for _,field in ipairs{'xmin', 'xmax', 'startCoord'} do
			for j=1,n do
				float[0] = self.geom[field][j]
				if inputFloatTooltip(field..' '..j, float, 0, 0, -1, ig.ImGuiInputTextFlags_EnterReturnsTrue) then
					self.geom[field][j] = float[0]
					self:rebuildSurface()
				end
				if j < n then ig.igSameLine() end
			end
		end
	end
	ig.igEnd()
end

App():run()
