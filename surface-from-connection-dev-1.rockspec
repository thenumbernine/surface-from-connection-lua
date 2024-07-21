package = "surface-from-connection"
version = "dev-1"
source = {
	url = "git+https://github.com/thenumbernine/surface-from-connection-lua"
}
description = {
	detailed = "reconstructing surfaces from their coordinate connection coefficients.",
	homepage = "https://github.com/thenumbernine/surface-from-connection-lua",
	license = "MIT"
}
dependencies = {
	"lua >= 5.1"
}
build = {
	type = "builtin",
	modules = {
		["surface-from-connection.run"] = "run.lua"
	},
}
