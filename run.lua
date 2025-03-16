#! /usr/bin/env luajit
-- out of date with the current run.rua
local math = require 'ext.math'
local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'imgui'
local gl = require 'gl'
local matrix = require 'matrix'
local symmath = require 'symmath'
local ffi = require 'ffi'
local vec3f = require 'vec-ffi.vec3f'
local GLSceneObject = require 'gl.sceneobject'

--[[
local CLEnv = require 'cl.obj.env'
--]]

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

function matrix:isfinite()
	for i in self:iter() do
		if not math.isfinite(self[i]) then return false end
	end
	return true
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
		local var = symmath.var(name)
		for _,p in ipairs(require 'symmath.tensor.symbols'.greekSymbolsAndNames) do
			local k,v = next(p)
			if var.name == v then
				var:nameForExporter('Lua', k)
			end
		end
		return var
	end)

	if self.createMetric then
		local chart = symmath.Tensor.Chart{coords=self.coordVars}
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
	local gnuplot = require 'gnuplot'
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
	local size = expr:dim()
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
	* create_gs - to calculate the metric numerically
	* createMetric - to calculate the metric analytically based on the coordinate chart
	* create_conns - to use an identity metric and numerical connections
	* ... otherwise these return identity g_ab's and zero Γ^a_bc's.
--]]


-- 2D geometries


local Polar = Geometry:subclass()
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

-- [[ purely numerically
-- ... if this and create_gs/createMetric is enabled then this is only used to calculate numerical error of connections
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
local PolarAnholonomic = Geometry:subclass()
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
local SphereSurface = Geometry:subclass()
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


local TorusSurface = Geometry:subclass()
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

local PoincareDisk2D = Geometry:subclass()
PoincareDisk2D.coords = {'u', 'v'}
PoincareDisk2D.xmin = {-1, -1}
PoincareDisk2D.xmax = {1, 1}
PoincareDisk2D.startCoord = {0, 0}
function PoincareDisk2D:createMetric()
	local u, v = self.coordVars:unpack()
	return symmath.Tensor('_ab', {4 / (1 - u^2 - v^2), 0}, {0, 4 / (1 - u^2 - v^2)})
end


-- hmm, how to incorporate signature into the metric ...
local Minkowski2D = Geometry:subclass()
Minkowski2D.coords = {'t', 'x'}
Minkowski2D.xmin = {-1, -1}
Minkowski2D.xmax = {1, 1}
Minkowski2D.startCoord = {0,0}
function Minkowski2D:createMetric()
	return symmath.Tensor('_ab', {-1, 0}, {0, 1})
end


-- here's Schwarzschild in time and in radial
-- it is treating Rs as constant, which means this metric is true for the spacetime *outside* of the massive body
local Schwarzschild1Plus1 = Geometry:subclass()
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
-- this is my first metric that is based on numerically specifying g_ab, then computing Γ^a_bc
-- since this doesn't have an extra dimension to anchor it, as the spacetime grows from r+ to r- it gets really twisted
local Schwarzschild1Plus1EOS = Geometry:subclass()
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
local LagrangianTotalEnergy = Geometry:subclass()
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


local Cylinder = Geometry:subclass()
Cylinder.coords = {'r', 'θ', 'z'}
Cylinder.xmin = {1, 0, -5}
Cylinder.xmax = {10, 2*math.pi, 5}
Cylinder.startCoord = {1,math.pi,0}
function Cylinder:createMetric()
	local r, theta, z = self.coordVars:unpack()
	return symmath.Tensor('_ab', {1, 0, 0}, {0, r^2, 0}, {0, 0, 1})
end


local Sphere = Geometry:subclass()
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


local Torus = Geometry:subclass()
Torus.coords = {'r', 'θ', 'φ'}
Torus.xmin = {1, -math.pi, -math.pi}
Torus.xmax = {2, math.pi, math.pi}
Torus.startCoord = {1, -math.pi, -math.pi}
function Torus:createMetric()
	local r, theta, phi = self.coordVars:unpack()
	local R = 5
	return symmath.Tensor('_ab', {1, 0, 0}, {0, r^2, 0}, {0, 0, (R + r * symmath.sin(theta))^2})	-- does sin(theta) work as well?
end


local PoincareDisk3D = Geometry:subclass()
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
-- this is my first metric that is based on numerically specifying g_ab, then computing Γ^a_bc
-- since this doesn't have an extra dimension to anchor it, as the spacetime grows from r+ to r- it gets really twisted
local Schwarzschild2Plus1EOS = Geometry:subclass()
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


local SchwarzschildSphere2Plus1EOS = Geometry:subclass()
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
local UniformElectricFieldNumericIn2Plus1D = Geometry:subclass()
UniformElectricFieldNumericIn2Plus1D.coords = {'t', 'x', 'y'}		-- ut oh, now we introduce metric signatures ...
UniformElectricFieldNumericIn2Plus1D.xmin = {-1, -1, -1}
UniformElectricFieldNumericIn2Plus1D.xmax = {1, 1, 1}
UniformElectricFieldNumericIn2Plus1D.startCoord = {0, 0, 0}
function UniformElectricFieldNumericIn2Plus1D:create_conns()
	local E = 1
	return self.app.size:lambda(function(...)
		return matrix{
			{
				{0,E,0},
				{E,0,0},
				{0,0,0},
			},
			{
				{-E,0,0},
				{0,0,0},
				{0,0,E},
			},
			{
				{0,0,0},
				{0,0,0},
				{0,0,0},
			},
		}
	end)
end


-- here's a connection coefficient that gives rise to the stress-energy of a uniform electric field
-- it's based on an analytical connection, but I haven't made a metric for it just yet
local UniformElectricFieldNumeric = Geometry:subclass()
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
local InfiniteWireMagneticFieldNumeric = Geometry:subclass()
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


local App = require 'imguiapp.withorbit'()
App.title = 'reconstruct surface from geodesics'
App.viewDist = 10

function App:initGL()
	App.super.initGL(self)
	gl.glEnable(gl.GL_DEPTH_TEST)

	self.controlsOpened = true
	self.geomID = 1

	self.lineVtxs = ffi.new'vec3f_t[2]'
	self.lineObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			vertexCode = [[
in vec3 vertex;
uniform mat4 mvProjMat;
void main() {
	gl_Position = mvProjMat * vec4(vertex, 1.);
}
]],
			fragmentCode = [[
out vec4 fragColor;
uniform vec3 color;
void main() {
	fragColor = vec4(color, 1.);
}
]],
		},
		vertexes = {
			data = ffi.cast('float*', self.lineVtxs),
			size = ffi.sizeof'vec3f_t' * 2,
			dim = 3,
			count = 2,
		},
		geometry = {
			mode = gl.GL_LINES,
		},
	}

	self:buildSurface'Polar'
end

local geomClassesForName = table{
	-- 2D
	{Polar = Polar},
	{PolarAnholonomic = PolarAnholonomic},
	{SphereSurface = SphereSurface},
	{TorusSurface = TorusSurface},
--	{PoincareDisk2D = PoincareDisk2D},	-- crashing
	{Minkowski2D = Minkowski2D},
	{Schwarzschild1Plus1 = Schwarzschild1Plus1},
	{Schwarzschild1Plus1EOS = Schwarzschild1Plus1EOS},
	{LagrangianTotalEnergy = LagrangianTotalEnergy},
	-- 3D
	{Cylinder = Cylinder},
	{Sphere = Sphere},
	{Torus = Torus},
--	{PoincareDisk3D = PoincareDisk3D},	-- crashing
	{UniformElectricFieldNumericIn2Plus1D = UniformElectricFieldNumericIn2Plus1D },
	{Schwarzschild2Plus1EOS = Schwarzschild2Plus1EOS},
	{SchwarzschildSphere2Plus1EOS = SchwarzschildSphere2Plus1EOS},
	-- 4D
--	{UniformElectricFieldNumeric = UniformElectricFieldNumeric},	-- crashing
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
	self.geomID = loc

--[[
	self.env = CLEnv()
--]]

	self.geom = geomClass(self)

	local n = #self.geom.coords
	self.size = matrix{n}:lambda(I( ({[2]=64, [3]=16, [4]=8})[n] ))

--[[
	self.domain = self.env:domain{size=self.size}
--]]

	self.xmin = self.geom and self.geom.xmin or matrix{n}:lambda(I(-1))
	self.xmax = self.geom and self.geom.xmax or matrix{n}:lambda(I(1))

	self.view.ortho = n == 2
	-- reset view every time you change charts
	self.view.pos:set(0,0,self.viewDist)
	self.view.orbit:set(0,0,0)
	self.view.angle:set(0,0,0,1)

--[=[ opencl code
	-- regenerate these to prevent ffi cdef redefs
	self.rank1Type = 'rank1_t'
	self.rank2Type = 'rank2_t'
	local rank1TypeCode = 'typedef real '..self.rank1Type..'['..n..'];'
	local rank2TypeCode = 'typedef real '..self.rank2Type..'['..n..'];'
	local ffi = require 'ffi'
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
	local template = require 'template'
	local clnumber = require 'cl.obj.number'

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
print()
print('App:rebuildSurface...')
	local n = #self.size

	-- if we are calculating the connection from discrete derivatives of the metric ...
	if self.geom.createMetric
	or self.geom.create_gs
	then
print('creating metric...')
		local gs
		if self.geom.create_gs then
print('creating metric numerically...')
			gs = self.geom:create_gs()
		else
print('creating metric analytically...')
			gs = self.geom:calcFromEqns_gs()
		end
		--or self.size:lambda(function(...) return matrix{n,n}:eye() end)
print('creating metric inverse numerically...')
		local gUs = self.size:lambda(function(...)
			return gs[{...}]:inv()
		end)
print('creating metric derivatives numerically...')
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
			-- reshape dg_abc = g_ab,c
			return matrix{n,n,n}:lambda(function(a,b,c) return dg[c][a][b] end)
		end)

print('creating connections numerically...')
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
print('*** comparing numeric to exact ***')
			self.geom:testExact()
		end
	else
print('creating connections numerically...')
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
	self.es[i] = matrix{n,n}:eye()
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
							local complex = require 'complex'
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

	local vertexes = table()
	local vertex2s = table()
	local colors = table()
	local sizeMinusOne = self.size-1
	for indexMinusOne in (sizeMinusOne-1):range() do
		local index = indexMinusOne+1
		for k=1,n do
			if index[k] < sizeMinusOne[k] then
				local nextIndex = matrix(index)
				nextIndex[k] = nextIndex[k] + 1

				local color = (index-1):ediv(sizeMinusOne)
				colors:append{color[1] or 0, color[2] or 0, color[3] or .5}
				local vertex2 = self.xs[index]
				vertex2s:append{vertex2[1] or 0, vertex2[2] or 0, vertex2[3] or 0}
				local vertex = self.Xs[index]
				vertexes:append{vertex[1] or 0, vertex[2] or 0, vertex[3] or 0}

				local color = (nextIndex-1):ediv(sizeMinusOne)
				colors:append{color[1] or 0, color[2] or 0, color[3] or .5}
				local vertex2 = self.xs[nextIndex]
				vertex2s:append{vertex2[1] or 0, vertex2[2] or 0, vertex2[3] or 0}
				local vertex = self.Xs[nextIndex]
				vertexes:append{vertex[1] or 0, vertex[2] or 0, vertex[3] or 0}
			end
		end
	end

	self.surfaceObj = GLSceneObject{
		program = {
			version = 'latest',
			precision = 'best',
			vertexCode = [[
in vec3 vertex;
in vec3 vertex2;
in vec3 color;
out vec3 colorv;
uniform float t;
uniform mat4 mvProjMat;
void main() {
	colorv = color;
	vec3 pos1 = vertex.xyz;
	vec3 pos2 = vertex2.xyz;
	vec3 pos = mix(pos1, pos2, t);
	gl_Position = mvProjMat * vec4(pos, 1.);
}
]],
			fragmentCode = [[
in vec3 colorv;
out vec4 fragColor;
void main() {
	fragColor = vec4(colorv, 1.);
}
]],
			uniforms = {t = 0},
		},
		geometry = {
			mode = gl.GL_LINES,
		},
		vertexes = {
			data = vertexes,
			count = #vertexes / 3,
			dim = 3,
		},
		attrs = {
			vertex2 = {
				buffer = {
					data = vertex2s,
					count = #vertex2s / 3,
					dim = 3,
				},
			},
			color = {
				buffer = {
					data = colors,
					count = #colors / 3,
					dim = 3,
				},
			},
		},
	}
end

function App:drawGrid()
	local xmin, xmax, ymin, ymax = self.view:getBounds(self.width / self.height)
	xmin = xmin + self.view.pos.x
	ymin = ymin + self.view.pos.y
	xmax = xmax + self.view.pos.x
	ymax = ymax + self.view.pos.y

	self.lineObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
	self.lineObj.uniforms.color = {.1, .1, .1}

	local xrange = xmax - xmin
	local xstep = 10^math.floor(math.log(xrange, 10) - .5)
	local xticmin = math.floor(xmin/xstep)
	local xticmax = math.ceil(xmax/xstep)

	for x=xticmin,xticmax do
		self.lineVtxs[0]:set(x*xstep, ymin, 0)
		self.lineVtxs[1]:set(x*xstep, ymax, 0)
		self.lineObj.vertexes:bind():updateData():unbind()
		self.lineObj:draw()
	end

	local yrange = ymax - ymin
	local ystep = 10^math.floor(math.log(yrange, 10) - .5)
	local yticmin = math.floor(ymin/ystep)
	local yticmax = math.ceil(ymax/ystep)
	for y=yticmin,yticmax do
		self.lineVtxs[0]:set(xmin, y*ystep, 0)
		self.lineVtxs[1]:set(xmax, y*ystep, 0)
		self.lineObj.vertexes:bind():updateData():unbind()
		self.lineObj:draw()
	end

	self.lineObj.uniforms.color = {.5, .5, .5}

	self.lineVtxs[0]:set(xmin, 0, 0)
	self.lineVtxs[1]:set(xmax, 0, 0)
	self.lineObj.vertexes:bind():updateData():unbind()
	self.lineObj:draw()

	self.lineVtxs[0]:set(0, ymin, 0)
	self.lineVtxs[1]:set(0, ymax, 0)
	self.lineObj.vertexes:bind():updateData():unbind()
	self.lineObj:draw()
end


local animating = false
local lastTime = 0
local gui = {
	animTime = 0,
}
function App:update()
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

	local n = #self.size

	local thisTime = os.clock()
	if animating then
		local deltaTime = (thisTime - lastTime) / math.pi
		gui.animTime = (((gui.animTime + deltaTime) + 1) % 2) - 1
	end
	lastTime = thisTime

	self.surfaceObj.uniforms.t = .5 - .5 * math.cos(math.pi * gui.animTime)
	self.surfaceObj.uniforms.mvProjMat = self.view.mvProjMat.ptr
	self.surfaceObj:draw()

--[[ draw the basis at every point?
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
				glColor3f(colors[k])
				glVertex(u)
				glVertex(u + scale * e[k])
			end
		end
		gl.glEnd()
		gl.glPopMatrix()
--]]

	if n == 2 then
		self:drawGrid()
	elseif n == 3 then
		-- draw the xyz basis
		-- TODO viewport instead?
		local hudView = require 'glapp.view'{
			angle = self.view.angle,
		}
		local vec3d = require 'vec-ffi.vec3d'
		local aspectRatio = self.width / self.height
		hudView.pos = hudView.angle:rotate(vec3d(4 * aspectRatio, 4, 5))
		hudView:setup(aspectRatio)

		-- draw basis
		self.lineObj.uniforms.mvProjMat = hudView.mvProjMat.ptr

		for i=1,3 do
			self.lineObj.uniforms.color = {0, 0, 0}
			self.lineObj.uniforms.color[i] = 1

			self.lineVtxs[0]:set(0,0,0)
			self.lineVtxs[1]:set(0,0,0)
			self.lineVtxs[1].s[i-1] = 1
			self.lineObj.vertexes:bind():updateData():unbind()
			self.lineObj:draw()
		end
	end

	App.super.update(self)
end

function App:updateGUI()
	if ig.luatableBegin('Controls', self, 'controlsOpened') then
		if ig.igButton(animating and 'Stop Animation' or 'Start Animation') then
			animating = not animating
		end

		ig.luatableTooltipSliderFloat('animation coefficient', gui, 'animTime', -1, 1)

		if ig.luatableTooltipCombo('coordinate system', self, 'geomID', geomClassNames) then
			self:buildSurface(geomClassNames[self.geomID])
		end

		local n = #self.size
		for _,field in ipairs{'xmin', 'xmax', 'startCoord'} do
			for j=1,n do
				if ig.luatableTooltipInputFloat(field..' '..j, self.geom[field], j, .01, .1, '%f', ig.ImGuiInputTextFlags_EnterReturnsTrue) then
					self:rebuildSurface()
				end
				if j < n then ig.igSameLine() end
			end
		end
	end
	ig.igEnd()
end

App():run()
