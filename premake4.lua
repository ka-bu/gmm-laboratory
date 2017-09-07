solution "GMMLab"
   configurations { "Debug" }

   location "build"

project "gmmlab"
kind "ConsoleApp"
language "C++"
location "build"
targetdir "bin"

includedirs { "include"}
libdirs { "lib" }



links { "dl", "X11", "Xrandr", "Xi", "Xxf86vm","Xinerama", "Xcursor", "GL", "GLU", "rt",
	"pthread", "AntTweakBar", "glfw3",
	"boost_program_options",
	"boost_system",
	"boost_filesystem"
}

configuration { "Debug" }
    flags { "Symbols" }
    configuration {}

buildoptions { "-std=c++0x -O2" }

files
{
    "src/**.h",
    "src/**.cpp"
}

excludes
{
    "src/main_testlab.cpp",
    

    
}

project "testlab"
kind "ConsoleApp"
language "C++"
location "build"
targetdir "bin"

includedirs { "include" }
libdirs { "lib" }
links {
	--"boost_regex",
	"boost_program_options",
	"boost_system",
	"boost_filesystem"
}

configuration { "Debug" }
    flags { "Symbols" }
    configuration {}

buildoptions { "-std=c++0x -O2" }

files
{
    "src/**.h",
    "src/**.cpp"
}

excludes
{
    "src/main_gmmlab.cpp",
    "src/display/gmmdisplay.*",
    "src/display/displaybar.*",
    "src/display/observationbar.*",
        
    
}
